//! Compute molecular consequence of variants.

use std::{collections::HashMap, rc::Rc};

use hgvs::{
    data::interface::{Provider, TxForRegionRecord},
    mapper::assembly::{Config as AssemblyConfig, Mapper as AssemblyMapper},
    parser::{Accession, GenomeInterval, GenomeLocEdit, HgvsVariant, Mu, NaEdit},
    static_data::Assembly,
};

use crate::db::create::txs::data::{Strand, TranscriptBiotype};

use super::{
    ann::{Allele, AnnField, Consequence, FeatureBiotype, FeatureType, Pos, Rank, SoFeature},
    provider::MehariProvider,
};

/// A variant description how VCF would do it.
#[derive(Debug, PartialEq, Eq, Clone, Default)]
pub struct VcfVariant {
    /// Chromosome name.
    pub chromosome: String,
    /// 1-based position on the chromosome of first base of `reference`.
    pub position: i32,
    /// Reference bases.
    pub reference: String,
    /// Alternative bases.
    pub alternative: String,
}

pub struct ConsequencePredictor {
    /// The internal transcript provider for locating transcripts.
    provider: Rc<MehariProvider>,
    /// Assembly mapper for variant consequence prediction.
    mapper: AssemblyMapper,
    /// Mapping from chromosome name to accession.
    chrom_to_acc: HashMap<String, String>,
}

/// Padding to look for genes upstream/downstream.
pub const PADDING: i32 = 5_000;
/// Generally used alternative alignment method.
pub const ALT_ALN_METHOD: &str = "splign";

impl ConsequencePredictor {
    pub fn new(provider: Rc<MehariProvider>, assembly: Assembly) -> Self {
        let acc_to_chrom = provider.get_assembly_map(assembly);
        let mut chrom_to_acc = HashMap::new();
        for (acc, chrom) in &acc_to_chrom {
            let chrom = if chrom.starts_with("chr") {
                chrom.strip_prefix("chr").unwrap()
            } else {
                chrom
            };
            chrom_to_acc.insert(chrom.to_string(), acc.clone());
            chrom_to_acc.insert(format!("chr{}", chrom), acc.clone());
        }

        let config = AssemblyConfig {
            replace_reference: false,
            strict_bounds: false,
            renormalize_g: false,
            ..Default::default()
        };
        let mapper = AssemblyMapper::new(config, provider.clone());

        ConsequencePredictor {
            provider,
            mapper,
            chrom_to_acc,
        }
    }

    pub fn predict(&self, var: &VcfVariant) -> Result<Option<Vec<AnnField>>, anyhow::Error> {
        // Obtain accession from chromosome name.
        let chrom_acc = self.chrom_to_acc.get(&var.chromosome);
        let chrom_acc = if let Some(chrom_acc) = chrom_acc {
            chrom_acc
        } else {
            tracing::debug!(
                "Could not determine chromosome accession for {:?}; giving up on annotation",
                &var
            );
            return Ok(None);
        };

        // Get all affected transcripts.
        let var_start = var.position;
        let var_end = var.position + var.reference.len() as i32;
        let qry_start = var_start - PADDING;
        let qry_end = var_end + PADDING;
        let mut txs =
            self.provider
                .get_tx_for_region(chrom_acc, ALT_ALN_METHOD, qry_start, qry_end)?;
        txs.sort_by(|a, b| a.tx_ac.cmp(&b.tx_ac));

        // Generate `AnnField` records for each transcript.
        Ok(Some(
            txs.into_iter()
                .map(|tx| self.build_ann_field(var, tx, chrom_acc.clone(), var_start - 1, var_end))
                .collect::<Result<Vec<_>, _>>()?,
        ))
    }

    fn build_ann_field(
        &self,
        var: &VcfVariant,
        tx_record: TxForRegionRecord,
        chrom_acc: String,
        var_start: i32,
        var_end: i32,
    ) -> Result<AnnField, anyhow::Error> {
        // NB: The coordinates of var_start, var_end, and the transcript are all 0-based.

        let tx = self.provider.get_tx(&tx_record.tx_ac).unwrap();
        let mut consequences: Vec<Consequence> = Vec::new();

        let alignment = tx.genome_alignments.first().unwrap();

        let mut min_start = None;
        let mut max_end = None;

        // Find first exon that overlaps with variant or intron that contains the variant.
        //
        // Note that exons are stored in genome position order.
        let mut prev_end = None;
        let mut rank = Rank::default();
        let mut is_exonic = false;
        let mut is_intronic = false;
        let mut distance = None;
        let mut tx_len = 0;
        for exon_alignment in &alignment.exons {
            tx_len += exon_alignment.alt_end_i - exon_alignment.alt_start_i;

            let exon_start = exon_alignment.alt_start_i;
            let exon_end = exon_alignment.alt_end_i;
            let intron_start = prev_end;
            let intron_end = exon_start;

            // Check the cases where the variant overlaps with the exon or is contained within an
            // intron.
            if var_start < exon_end && exon_start < var_end {
                // overlaps with exon
                rank = Rank {
                    ord: exon_alignment.ord + 1,
                    total: alignment.exons.len() as i32,
                };
                is_exonic = true;

                if var_start <= exon_start || var_end >= exon_end {
                    // overlaps with exon/intron boundary
                    distance = Some(0);
                } else {
                    let dist_start = -(var_start - exon_start + 1);
                    let dist_end = exon_end - var_end + 1;
                    if dist_end >= dist_start.abs() {
                        distance = Some(dist_end);
                    } else {
                        distance = Some(dist_start);
                    }
                }
            } else if let Some(intron_start) = intron_start {
                if var_start > intron_start && var_end < intron_end {
                    // contained within intron: cannot be in next exon
                    if !is_exonic {
                        rank = Rank {
                            ord: exon_alignment.ord + 1,
                            total: alignment.exons.len() as i32 - 1,
                        };
                        is_intronic = true;

                        let dist_start = -(var_start - intron_start);
                        let dist_end = intron_end - var_end;
                        if dist_end >= dist_start.abs() {
                            distance = Some(dist_end);
                        } else {
                            distance = Some(dist_start);
                        }
                    }
                }
            }

            // Check the cases where the variant overlaps with whole exon.
            if var_start <= exon_start && var_end >= exon_end {
                consequences.push(Consequence::ExonLossVariant);
            }
            // Check the cases where the variant overlaps with the splice acceptor/donor site.
            if let Some(intron_start) = intron_start {
                if var_start < intron_start + 2 && var_end > intron_start {
                    // Left side, is acceptor/donor depending on transcript's strand.
                    match alignment.strand {
                        Strand::Plus => consequences.push(Consequence::SpliceDonorVariant),
                        Strand::Minus => consequences.push(Consequence::SpliceAcceptorVariant),
                    }
                }
            }
            // Check the case where the variant overlaps with the splice donor site.
            if var_start < intron_end && var_end > intron_end - 2 {
                // Left side, is acceptor/donor depending on transcript's strand.
                match alignment.strand {
                    Strand::Plus => consequences.push(Consequence::SpliceAcceptorVariant),
                    Strand::Minus => consequences.push(Consequence::SpliceDonorVariant),
                }
            }
            // Check the case where the variant overlaps with the splice region (1-3 bases in exon
            // or 3-8 bases in intron).
            if let Some(intron_start) = intron_start {
                if (var_start < exon_end && var_end > exon_end - 3)
                    || (var_start < intron_start + 8 && var_end > intron_start + 2)
                    || (var_start < intron_end - 8 && var_end > intron_end - 2)
                    || (var_start < intron_end + 3 && var_end > intron_end)
                {
                    consequences.push(Consequence::SpliceRegionVariant);
                }
            }
            // Check the case where the variant fully contains the stop codon, based on `var_c` coordinates.
            todo!();
            // Check the case where the variant fully contains the start codon, based on `var_c` coordinates.
            todo!();

            min_start = Some(std::cmp::min(min_start.unwrap_or(exon_start), exon_start));
            max_end = Some(std::cmp::min(max_end.unwrap_or(exon_end), exon_end));
            prev_end = Some(exon_end);
        }

        let min_start = min_start.expect("must have seen exon");
        let max_end = max_end.expect("must have seen exon");

        let is_upstream = var_end < min_start;
        let is_downstream = var_start > max_end;
        if is_exonic {
            consequences.push(Consequence::ExonVariant);
        } else if is_intronic {
            consequences.push(Consequence::IntronVariant);
        } else if is_upstream {
            let val = -(min_start - var_end + 1);
            if val.abs() < 5_000 {
                match alignment.strand {
                    Strand::Plus => consequences.push(Consequence::UpstreamGeneVariant),
                    Strand::Minus => consequences.push(Consequence::DownstreamGeneVariant),
                }
            }
            distance = Some(val);
        } else if is_downstream {
            let val = var_start - max_end + 1;
            if val.abs() < 5_000 {
                match alignment.strand {
                    Strand::Plus => consequences.push(Consequence::DownstreamGeneVariant),
                    Strand::Minus => consequences.push(Consequence::UpstreamGeneVariant),
                }
            }
            distance = Some(val);
        }

        let feature_biotype = match tx.biotype {
            TranscriptBiotype::Coding => FeatureBiotype::Coding,
            TranscriptBiotype::NonCoding => FeatureBiotype::Noncoding,
        };

        let var_g = HgvsVariant::GenomeVariant {
            accession: Accession { value: chrom_acc },
            gene_symbol: None,
            loc_edit: GenomeLocEdit {
                loc: Mu::Certain(GenomeInterval {
                    start: Some(var.position),
                    end: Some(var.position + var.reference.len() as i32 - 1),
                }),
                edit: Mu::Certain(NaEdit::RefAlt {
                    reference: var.reference.clone(),
                    alternative: var.alternative.clone(),
                }),
            },
        };

        let var_n = self.mapper.g_to_n(&var_g, &tx.id)?;
        let tx_pos = match &var_n {
            HgvsVariant::TxVariant { loc_edit, .. } => Some(Pos {
                ord: loc_edit.loc.inner().start.base,
                total: Some(tx_len),
            }),
            _ => panic!("Invalid tx position: {:?}", &var_n),
        };

        let (var_t, _var_p, hgvs_p, cds_pos, protein_pos) = match feature_biotype {
            FeatureBiotype::Coding => {
                let cds_len = tx.stop_codon.unwrap() - tx.start_codon.unwrap();
                let prot_len = cds_len / 3;

                let var_c = self.mapper.n_to_c(&var_n)?;
                let var_p = self.mapper.c_to_p(&var_c)?;
                let hgvs_p = Some(format!("{}", &var_p));
                let cds_pos = match &var_c {
                    HgvsVariant::CdsVariant { loc_edit, .. } => Some(Pos {
                        ord: loc_edit.loc.inner().start.base,
                        total: Some(cds_len),
                    }),
                    _ => panic!("Invalid CDS position: {:?}", &var_n),
                };
                let protein_pos = match &var_p {
                    HgvsVariant::ProtVariant { loc_edit, .. } => match &loc_edit {
                        hgvs::parser::ProtLocEdit::Ordinary { loc, .. } => Some(Pos {
                            ord: loc.inner().start.number,
                            total: Some(prot_len),
                        }),
                        _ => None,
                    },
                    _ => panic!("Not a protein position: {:?}", &var_n),
                };
                (var_c, Some(var_p), hgvs_p, cds_pos, protein_pos)
            }
            FeatureBiotype::Noncoding => (var_n, None, None, None, None),
        };
        let hgvs_t = Some(format!("{}", &var_t));

        // Take a highest-ranking consequence and derive putative impact from it.
        consequences.sort();
        if consequences.is_empty() {
            tracing::warn!("No consequences for {:?} on {}", var, &tx_record.tx_ac);
        }
        let putative_impact = (*consequences.first().unwrap()).into();

        // Build and return ANN field from the information derived above.
        Ok(AnnField {
            allele: Allele::Alt {
                alternative: var.alternative.clone(),
            },
            consequences,
            putative_impact,
            gene_symbol: tx.gene_name.clone(),
            gene_id: format!("HGNC:{}", &tx.gene_id),
            feature_type: FeatureType::SoTerm {
                term: SoFeature::Transcript,
            },
            feature_id: tx.id.clone(),
            feature_biotype,
            rank,
            hgvs_t,
            hgvs_p,
            tx_pos,
            cds_pos,
            protein_pos,
            distance,
            messages: None,
        })
    }
}

#[cfg(test)]
mod test {
    use std::{fs::File, io::BufReader};

    use csv::ReaderBuilder;
    use pretty_assertions::assert_eq;
    use serde::Deserialize;

    use crate::annotate::seqvars::{ann::PutativeImpact, load_tx_db};

    use super::*;

    #[test]
    fn annotate_snv_brca1_one_variant() -> Result<(), anyhow::Error> {
        let tx_path = "tests/data/annotate/db/seqvars/grch37/txs.bin";
        let tx_db = load_tx_db(tx_path, 5_000_000)?;
        let provider = Rc::new(MehariProvider::new(tx_db, Assembly::Grch37p10));

        let predictor = ConsequencePredictor::new(provider, Assembly::Grch37p10);

        let res = predictor
            .predict(&VcfVariant {
                chromosome: String::from("17"),
                position: 41_197_701,
                reference: String::from("G"),
                alternative: String::from("C"),
            })?
            .unwrap();

        assert_eq!(res.len(), 4);
        assert_eq!(
            res[0],
            AnnField {
                allele: Allele::Alt {
                    alternative: String::from("C")
                },
                consequences: vec![Consequence::ExonVariant],
                putative_impact: PutativeImpact::Modifier,
                gene_symbol: String::from("BRCA1"),
                gene_id: String::from("HGNC:1100"),
                feature_type: FeatureType::SoTerm {
                    term: SoFeature::Transcript
                },
                feature_id: String::from("NM_007294.4"),
                feature_biotype: FeatureBiotype::Coding,
                rank: Rank { ord: 23, total: 23 },
                hgvs_t: Some(String::from("NM_007294.4:c.5586C>G")),
                hgvs_p: Some(String::from("NP_009225.1:p.His1862Gln")),
                tx_pos: Some(Pos {
                    ord: 5699,
                    total: Some(7088)
                }),
                cds_pos: Some(Pos {
                    ord: 5586,
                    total: Some(5592)
                }),
                protein_pos: Some(Pos {
                    ord: 1862,
                    total: Some(1864)
                }),
                distance: Some(-1390),
                messages: None,
            }
        );

        Ok(())
    }

    #[derive(Debug, Deserialize)]
    struct Record {
        pub var: String,
        pub tx: String,
        pub csq: String,
    }

    #[test]
    fn annotate_brca1_clinvar_vars() -> Result<(), anyhow::Error> {
        let tx_path = "tests/data/annotate/db/seqvars/grch37/txs.bin";
        let tx_db = load_tx_db(tx_path, 5_000_000)?;
        let provider = Rc::new(MehariProvider::new(tx_db, Assembly::Grch37p10));
        let predictor = ConsequencePredictor::new(provider, Assembly::Grch37p10);

        let txs = vec![
            String::from("NM_007294.4"),
            String::from("NM_007297.4"),
            String::from("NM_007298.3"),
            String::from("NM_007299.4"),
            String::from("NM_007300.4"),
        ];
        let path_tsv = "tests/data/annotate/vars/clinvar.excerpt.snpeff.tsv";

        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(File::open(path_tsv).map(BufReader::new)?);

        for record in reader.deserialize() {
            let record: Record = record?;
            if txs.contains(&record.tx) {
                let arr = record.var.split('-').collect::<Vec<_>>();
                let anns = predictor
                    .predict(&VcfVariant {
                        chromosome: arr[0].to_string(),
                        position: arr[1].parse::<i32>()?,
                        reference: arr[2].to_string(),
                        alternative: arr[3].to_string(),
                    })?
                    .unwrap();
                for ann in &anns {
                    if ann.feature_id == record.tx {
                        assert_eq!(
                            ann.consequences
                                .iter()
                                .map(|c| c.to_string())
                                .collect::<Vec<_>>()
                                .join("&"),
                            record.csq,
                            "variant: {}, tx: {}, hgvs_c: {:?}, hgvs_p: {:?}",
                            record.var,
                            record.tx,
                            ann.hgvs_t.as_ref(),
                            ann.hgvs_p.as_ref()
                        );
                    }
                }
            }
        }

        Ok(())
    }
}
