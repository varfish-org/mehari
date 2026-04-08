use crate::annotate::genotype_string;
use crate::annotate::seqvars::ann::{AnnField, FeatureBiotype};
use crate::annotate::seqvars::{AnnotatedVariant, VariantAnnotation};
use crate::common::contig::ContigManager;
use crate::common::noodles::{NoodlesVariantReader, open_variant_reader};
use crate::common::{TsvContigStyle, guess_assembly_from_vcf};
use crate::ped::{PedigreeByName, Sex};
use anyhow::Error;
use biocommons_bioutils::assemblies::Assembly;
use clap::Parser;
use flate2::Compression;
use flate2::write::GzEncoder;
use futures::TryStreamExt;
use noodles::vcf::Header as VcfHeader;
use noodles::vcf::variant::record::samples::keys::key::{
    CONDITIONAL_GENOTYPE_QUALITY, GENOTYPE, READ_DEPTH,
};
use noodles::vcf::variant::record_buf::info::field;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::str::FromStr;

#[derive(Parser, Debug)]
#[command(about = "Export annotated VCF to VarFish TSV format", long_about = None)]
pub struct Args {
    /// Path to the annotated input VCF file.
    #[arg(long, required = true)]
    pub input: String,

    /// Path to the input PED file.
    #[arg(long, required = true)]
    pub pedigree: String,

    /// Path to the output TSV file.
    #[arg(short = 'o', long, required = true)]
    pub output: String,

    /// Path to HGNC TSV file.
    #[arg(long, required = true)]
    pub hgnc: String,

    /// Style for contig names in TSV output.
    #[arg(long, value_enum, default_value_t = TsvContigStyle::Auto)]
    pub tsv_contig_style: TsvContigStyle,
}

pub async fn run(_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!("Starting VarFish TSV export...");

    // 1. Load context files
    let hgnc_map = load_hgnc_map(&args.hgnc)?;
    tracing::info!("Loading pedigree from file {}", args.pedigree);
    let pedigree = PedigreeByName::from_path(&args.pedigree)?;

    // 2. Open VCF Reader and extract assembly
    let mut reader = open_variant_reader(&args.input).await?;
    let header = reader.read_header().await?;
    let assembly = guess_assembly_from_vcf(&header, true, None)?;

    // 3. Initialize TSV Writer
    let mut writer = VarFishSeqvarTsvWriter::with_path(&args.output, args.tsv_contig_style)?;
    writer.set_hgnc_map(hgnc_map);
    writer.set_pedigree(&pedigree);
    writer.set_assembly(assembly);
    writer.write_header()?;

    // 4. Process Records Stream
    let mut records = reader.records(&header).await;
    while let Some(record) = records.try_next().await? {
        // Parse ANN fields into the Sidecar
        let consequences =
            if let Some(Some(field::Value::Array(field::value::Array::String(arr)))) =
                record.info().get("ANN")
            {
                arr.iter()
                    .filter_map(|s| s.as_ref().map(|s| AnnField::from_str(s).unwrap()))
                    .collect::<Vec<_>>()
            } else {
                vec![]
            };

        let annotated_variant = AnnotatedVariant {
            vcf: record,
            annotation: VariantAnnotation {
                consequences,
                frequencies: None, // TSV writer pulls freqs directly from INFO below
                clinvar: None,
            },
        };

        writer.write_annotated_record(&header, &annotated_variant)?;
    }

    writer.flush()?;
    tracing::info!("... done writing VarFish TSV");
    Ok(())
}

fn load_hgnc_map(path: &str) -> anyhow::Result<FxHashMap<String, HgncRecord>> {
    tracing::info!("Loading HGNC map...");
    let mut result = FxHashMap::default();
    let tsv_file = File::open(path)?;
    let mut tsv_reader = csv::ReaderBuilder::new()
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .from_reader(tsv_file);
    for record in tsv_reader.deserialize() {
        let record: HgncRecord = record?;
        result.insert(record.hgnc_id.clone(), record);
    }
    tracing::info!("... done loading HGNC map");
    Ok(result)
}

// =========================================================================================
// RECOVERED STRUCTS FROM DIFF
// =========================================================================================

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct HgncRecord {
    pub hgnc_id: String,
    pub ensembl_gene_id: String,
    pub entrez_id: String,
    #[serde(alias = "symbol")]
    pub gene_symbol: String,
}

#[derive(Debug, Default)]
struct GenotypeInfo {
    pub name: String,
    pub gt: Option<String>,
    pub dp: Option<i32>,
    pub ad: Option<i32>,
    pub gq: Option<i32>,
    pub sq: Option<f32>,
}

#[derive(Debug, Default)]
struct GenotypeCalls {
    pub entries: Vec<GenotypeInfo>,
}

impl GenotypeCalls {
    pub fn for_tsv(&self) -> String {
        let mut result = String::new();
        result.push('{');
        let mut first = true;
        for entry in &self.entries {
            if first {
                first = false;
            } else {
                result.push(',');
            }
            result.push_str(&format!("\"\"\"{}\"\"\":{{", entry.name));

            let mut prev = false;
            if let Some(gt) = &entry.gt {
                prev = true;
                result.push_str(&format!("\"\"\"gt\"\"\":\"\"\"{}\"\"\"", gt));
            }
            if let Some(ad) = &entry.ad {
                if prev {
                    result.push(',');
                }
                prev = true;
                result.push_str(&format!("\"\"\"ad\"\"\":{}", ad));
            }
            if let Some(dp) = &entry.dp {
                if prev {
                    result.push(',');
                }
                prev = true;
                result.push_str(&format!("\"\"\"dp\"\"\":{}", dp));
            }
            if let Some(gq) = &entry.gq {
                if prev {
                    result.push(',');
                }
                result.push_str(&format!("\"\"\"gq\"\"\":{}", gq));
            } else if let Some(sq) = &entry.sq {
                if prev {
                    result.push(',');
                }
                result.push_str(&format!("\"\"\"gq\"\"\":{}", sq.round() as i32));
            }
            result.push('}');
        }
        result.push('}');
        result
    }
}

pub struct VarFishSeqvarTsvWriter {
    inner: Box<dyn Write>,
    assembly: Option<Assembly>,
    pedigree: Option<PedigreeByName>,
    header: Option<VcfHeader>,
    hgnc_map: Option<FxHashMap<String, HgncRecord>>,
    contig_manager: Option<ContigManager>,
    tsv_contig_style: TsvContigStyle,
}

impl VarFishSeqvarTsvWriter {
    pub fn with_path<P: AsRef<Path>>(
        p: P,
        tsv_contig_style: TsvContigStyle,
    ) -> anyhow::Result<Self> {
        let path = p.as_ref().to_path_buf();
        let inner: Box<dyn Write> = if path.extension().unwrap_or_default() == "gz" {
            Box::new(GzEncoder::new(File::create(&path)?, Compression::default()))
        } else {
            Box::new(File::create(&path)?)
        };
        Ok(Self {
            inner,
            assembly: None,
            pedigree: None,
            header: None,
            hgnc_map: None,
            contig_manager: None,
            tsv_contig_style,
        })
    }

    pub fn set_hgnc_map(&mut self, hgnc_map: FxHashMap<String, HgncRecord>) {
        self.hgnc_map = Some(hgnc_map)
    }

    pub fn set_assembly(&mut self, assembly: Assembly) {
        self.assembly = Some(assembly);
        self.contig_manager = Some(ContigManager::new(assembly));
    }

    pub fn set_pedigree(&mut self, pedigree: &PedigreeByName) {
        self.pedigree = Some(pedigree.clone())
    }

    pub fn write_header(&mut self) -> Result<(), anyhow::Error> {
        let header = &[
            "release",
            "chromosome",
            "chromosome_no",
            "start",
            "end",
            "bin",
            "reference",
            "alternative",
            "var_type",
            "case_id",
            "set_id",
            "info",
            "genotype",
            "num_hom_alt",
            "num_hom_ref",
            "num_het",
            "num_hemi_alt",
            "num_hemi_ref",
            "in_clinvar",
            "exac_frequency",
            "exac_homozygous",
            "exac_heterozygous",
            "exac_hemizygous",
            "thousand_genomes_frequency",
            "thousand_genomes_homozygous",
            "thousand_genomes_heterozygous",
            "thousand_genomes_hemizygous",
            "gnomad_exomes_frequency",
            "gnomad_exomes_homozygous",
            "gnomad_exomes_heterozygous",
            "gnomad_exomes_hemizygous",
            "gnomad_genomes_frequency",
            "gnomad_genomes_homozygous",
            "gnomad_genomes_heterozygous",
            "gnomad_genomes_hemizygous",
            "refseq_gene_id",
            "refseq_transcript_id",
            "refseq_transcript_coding",
            "refseq_hgvs_c",
            "refseq_hgvs_p",
            "refseq_effect",
            "refseq_exon_dist",
            "ensembl_gene_id",
            "ensembl_transcript_id",
            "ensembl_transcript_coding",
            "ensembl_hgvs_c",
            "ensembl_hgvs_p",
            "ensembl_effect",
            "ensembl_exon_dist",
        ];
        writeln!(self.inner, "{}", header.join("\t"))
            .map_err(|e| anyhow::anyhow!("Error writing VarFish TSV header: {}", e))
    }

    pub fn flush(&mut self) -> Result<(), anyhow::Error> {
        self.inner.flush()?;
        Ok(())
    }

    pub fn write_annotated_record(
        &mut self,
        header: &VcfHeader,
        record: &AnnotatedVariant,
    ) -> Result<(), anyhow::Error> {
        self.header = Some(header.clone());
        let mut tsv_record = VarFishSeqvarTsvRecord::default();
        if !self.fill_coords(&record.vcf, &mut tsv_record)? {
            return Ok(());
        }
        self.fill_genotype_and_freqs(&record.vcf, &mut tsv_record)?;
        self.fill_bg_freqs(&record.vcf, &mut tsv_record)?;
        self.fill_clinvar(&record.vcf, &mut tsv_record)?;
        self.expand_refseq_ensembl_and_write(record, &mut tsv_record)
    }

    fn fill_coords(
        &self,
        record: &noodles::vcf::variant::RecordBuf,
        tsv_record: &mut VarFishSeqvarTsvRecord,
    ) -> Result<bool, Error> {
        let assembly = self.assembly.expect("assembly must have been set");
        let contig_manager = self
            .contig_manager
            .as_ref()
            .expect("contig manager must be set");

        tsv_record.release = match assembly {
            Assembly::Grch37 | Assembly::Grch37p10 => String::from("GRCh37"),
            Assembly::Grch38 => String::from("GRCh38"),
        };
        let name = record.reference_sequence_name();

        if let Some(contig_info) = contig_manager.get_contig_info(name) {
            tsv_record.chromosome_no = contig_info.chrom_no;
            tsv_record.chromosome = match self.tsv_contig_style {
                TsvContigStyle::Passthrough => name.to_string(),
                TsvContigStyle::WithChr => contig_info.name_with_chr,
                TsvContigStyle::WithoutChr => contig_info.name_without_chr,
                TsvContigStyle::Auto => {
                    if assembly == Assembly::Grch38 {
                        if ContigManager::is_mitochondrial(contig_info.chrom_no) {
                            "chrM".into()
                        } else {
                            contig_info.name_with_chr
                        }
                    } else {
                        contig_info.name_without_chr
                    }
                }
            };
        } else {
            return Ok(false);
        }

        tsv_record.reference = record.reference_bases().to_string();
        tsv_record.alternative = record.alternate_bases().as_ref()[0].to_string();
        tsv_record.start = record
            .variant_start()
            .expect("Telomeres unsupported")
            .into();
        tsv_record.end = tsv_record.start + tsv_record.reference.len() - 1;
        tsv_record.bin = crate::annotate::seqvars::binning::bin_from_range(
            tsv_record.start as i32 - 1,
            tsv_record.end as i32,
        )? as u32;

        tsv_record.var_type = if tsv_record.reference.len() == tsv_record.alternative.len() {
            if tsv_record.reference.len() == 1 {
                String::from("snv")
            } else {
                String::from("mnv")
            }
        } else {
            String::from("indel")
        };

        Ok(true)
    }

    fn fill_genotype_and_freqs(
        &self,
        record: &noodles::vcf::variant::RecordBuf,
        tsv_record: &mut VarFishSeqvarTsvRecord,
    ) -> Result<(), anyhow::Error> {
        use noodles::vcf::variant::record_buf::samples::sample::Value;
        use noodles::vcf::variant::record_buf::samples::sample::value::Array;

        let hdr = self.header.as_ref().expect("VCF header must be set");
        let file_format_gt = noodles::vcf::header::FileFormat::new(4, 3);
        let mut gt_calls = GenotypeCalls::default();
        let samples = record.samples();
        let sample_names = hdr.sample_names().iter();

        let genotypes = samples.select(GENOTYPE);
        let get_genotype = |sample_idx| {
            genotypes.as_ref().and_then(|gt| {
                gt.get(sample_idx).map(|value| match value {
                    Some(Value::String(s)) => s.to_owned(),
                    Some(Value::Genotype(gt)) => genotype_string(gt, file_format_gt),
                    _ => ".".into(),
                })
            })
        };

        let read_depths = samples.select(READ_DEPTH);
        let get_read_depth = |sample_idx| {
            read_depths.as_ref().and_then(|rd| {
                rd.get(sample_idx)
                    .map(|value| match value {
                        Some(Value::Integer(i)) => Ok(*i),
                        None => Ok(-1),
                        _ => anyhow::bail!(format!("invalid DP value {:?}", value)),
                    })
                    .transpose()
                    .unwrap()
            })
        };

        let allele_depths = samples.select("AD");
        let get_allele_depths = |sample_idx| {
            allele_depths.as_ref().and_then(|ad| {
                ad.get(sample_idx)
                    .map(|value| match value {
                        Some(Value::Array(Array::Integer(arr))) => Ok(arr[1].unwrap_or(0)),
                        None => Ok(-1),
                        _ => anyhow::bail!(format!("invalid AD value {:?}", value)),
                    })
                    .transpose()
                    .unwrap()
            })
        };

        let conditional_gt_quality = samples.select(CONDITIONAL_GENOTYPE_QUALITY);
        let get_conditional_gt_quality = |sample_idx| {
            conditional_gt_quality.as_ref().and_then(|cgq| {
                cgq.get(sample_idx)
                    .map(|value| match value {
                        Some(Value::Integer(i)) => Ok(*i),
                        None => Ok(-1),
                        _ => anyhow::bail!("invalid GQ value"),
                    })
                    .transpose()
                    .unwrap()
            })
        };

        let sq = samples.select("SQ");
        let get_sq = |sample_idx| {
            sq.as_ref().and_then(|sq| {
                sq.get(sample_idx)
                    .map(|value| match value {
                        Some(Value::Float(f)) => Ok(*f),
                        Some(Value::Array(Array::Float(f))) => {
                            Ok(f[0].expect("SQ should be a single float"))
                        }
                        None => Ok(-1.0),
                        _ => anyhow::bail!("invalid SQ value"),
                    })
                    .transpose()
                    .unwrap()
            })
        };

        for (i, name) in sample_names.enumerate() {
            let mut gt_info = GenotypeInfo {
                name: name.clone(),
                ..Default::default()
            };
            if let Some(gt) = get_genotype(i) {
                let pedigree = self
                    .pedigree
                    .as_ref()
                    .ok_or_else(|| anyhow::anyhow!("Pedigree has not been set on the writer"))?;
                let individual = pedigree.individuals.get(name).ok_or_else(|| {
                    anyhow::anyhow!("Sample '{}' found in VCF but not in PED file", name)
                })?;
                if ContigManager::is_chr_x(tsv_record.chromosome_no) {
                    match individual.sex {
                        Sex::Male => {
                            if gt.contains('1') {
                                tsv_record.num_hemi_alt += 1;
                            } else {
                                tsv_record.num_hemi_ref += 1;
                            }
                        }
                        _ => {
                            let matches_1 = gt.matches('1').count();
                            let matches_0 = gt.matches('0').count();
                            if matches_0 == 2 {
                                tsv_record.num_hom_ref += 1;
                            } else if matches_1 == 1 {
                                tsv_record.num_het += 1;
                            } else if matches_1 == 2 {
                                tsv_record.num_hom_alt += 1;
                            }
                        }
                    }
                } else if ContigManager::is_chr_y(tsv_record.chromosome_no) {
                    if individual.sex == Sex::Male {
                        if gt.contains('1') {
                            tsv_record.num_hemi_alt += 1;
                        } else if gt.contains('0') {
                            tsv_record.num_hemi_ref += 1;
                        }
                    }
                } else {
                    let matches_1 = gt.matches('1').count();
                    let matches_0 = gt.matches('0').count();
                    if matches_0 == 2 {
                        tsv_record.num_hom_ref += 1;
                    } else if matches_1 == 1 {
                        tsv_record.num_het += 1;
                    } else if matches_1 == 2 {
                        tsv_record.num_hom_alt += 1;
                    }
                }
                gt_info.gt = Some(gt);
            }
            if let Some(dp) = get_read_depth(i) {
                gt_info.dp = Some(dp);
            }
            if let Some(ad) = get_allele_depths(i) {
                gt_info.ad = Some(ad);
            }
            if let Some(gq) = get_conditional_gt_quality(i) {
                gt_info.gq = Some(gq);
            }
            if let Some(sq) = get_sq(i) {
                gt_info.sq = Some(sq);
            }
            gt_calls.entries.push(gt_info);
        }
        tsv_record.genotype = gt_calls.for_tsv();
        Ok(())
    }

    fn fill_bg_freqs(
        &self,
        record: &noodles::vcf::variant::RecordBuf,
        tsv_record: &mut VarFishSeqvarTsvRecord,
    ) -> Result<(), anyhow::Error> {
        let gnomad_exomes_an = record
            .info()
            .get("gnomad_exomes_an")
            .unwrap_or_default()
            .map(|v| match v {
                field::Value::Integer(value) => *value,
                _ => 0,
            })
            .unwrap_or_default();
        if gnomad_exomes_an > 0 {
            tsv_record.gnomad_exomes_homozygous = record
                .info()
                .get("gnomad_exomes_hom")
                .unwrap_or(Some(&field::Value::Integer(0)))
                .map(|v| match v {
                    field::Value::Integer(value) => *value,
                    _ => 0,
                })
                .unwrap_or_default();
            tsv_record.gnomad_exomes_heterozygous = record
                .info()
                .get("gnomad_exomes_het")
                .unwrap_or(Some(&field::Value::Integer(0)))
                .map(|v| match v {
                    field::Value::Integer(value) => *value,
                    _ => 0,
                })
                .unwrap_or_default();
            tsv_record.gnomad_exomes_hemizygous = record
                .info()
                .get("gnomad_exomes_hemi")
                .unwrap_or(Some(&field::Value::Integer(0)))
                .map(|v| match v {
                    field::Value::Integer(value) => *value,
                    _ => 0,
                })
                .unwrap_or_default();
            tsv_record.gnomad_exomes_frequency = (tsv_record.gnomad_exomes_hemizygous
                + tsv_record.gnomad_exomes_heterozygous
                + tsv_record.gnomad_exomes_homozygous * 2)
                as f32
                / gnomad_exomes_an as f32;
        }
        let gnomad_genomes_an = record
            .info()
            .get("gnomad_genomes_an")
            .unwrap_or_default()
            .map(|v| match v {
                field::Value::Integer(value) => *value,
                _ => 0,
            })
            .unwrap_or_default();
        if gnomad_genomes_an > 0 {
            tsv_record.gnomad_genomes_homozygous = record
                .info()
                .get("gnomad_genomes_hom")
                .unwrap_or(Some(&field::Value::Integer(0)))
                .map(|v| match v {
                    field::Value::Integer(value) => *value,
                    _ => 0,
                })
                .unwrap_or_default();
            tsv_record.gnomad_genomes_heterozygous = record
                .info()
                .get("gnomad_genomes_het")
                .unwrap_or(Some(&field::Value::Integer(0)))
                .map(|v| match v {
                    field::Value::Integer(value) => *value,
                    _ => 0,
                })
                .unwrap_or_default();
            tsv_record.gnomad_genomes_hemizygous = record
                .info()
                .get("gnomad_genomes_hemi")
                .unwrap_or(Some(&field::Value::Integer(0)))
                .map(|v| match v {
                    field::Value::Integer(value) => *value,
                    _ => 0,
                })
                .unwrap_or_default();
            tsv_record.gnomad_genomes_frequency = (tsv_record.gnomad_genomes_hemizygous
                + tsv_record.gnomad_genomes_heterozygous
                + tsv_record.gnomad_genomes_homozygous * 2)
                as f32
                / gnomad_genomes_an as f32;
        }
        Ok(())
    }

    fn fill_clinvar(
        &self,
        record: &noodles::vcf::variant::RecordBuf,
        tsv_record: &mut VarFishSeqvarTsvRecord,
    ) -> Result<(), anyhow::Error> {
        tsv_record.in_clinvar = record
            .info()
            .get("clinvar_germline_classification")
            .is_some();
        Ok(())
    }

    fn expand_refseq_ensembl_and_write(
        &mut self,
        record: &AnnotatedVariant,
        tsv_record: &mut VarFishSeqvarTsvRecord,
    ) -> Result<(), anyhow::Error> {
        let empty_hgnc_record = HgncRecord {
            hgnc_id: "".to_string(),
            ensembl_gene_id: "".to_string(),
            entrez_id: "".to_string(),
            gene_symbol: "".to_string(),
        };

        if !record.annotation.consequences.is_empty() {
            let anns = record.annotation.consequences.clone();
            let mut anns_by_gene: FxHashMap<String, Vec<AnnField>> = FxHashMap::default();
            for ann in anns {
                let gene_id = ann.gene_id.clone();
                anns_by_gene.entry(gene_id).or_default().push(ann);
            }

            for anns in anns_by_gene.values_mut() {
                anns.sort_by_key(|ann| ann.consequences[0]);
            }

            let mut gene_ids: Vec<_> = anns_by_gene.keys().cloned().collect();
            gene_ids.sort_by_key(|id| {
                self.hgnc_map
                    .as_ref()
                    .unwrap()
                    .get(id)
                    .map(|r| r.entrez_id.clone())
                    .unwrap_or_else(|| id.clone())
            });

            for hgnc_id in gene_ids {
                let hgnc_record = self
                    .hgnc_map
                    .as_ref()
                    .unwrap()
                    .get(&hgnc_id)
                    .unwrap_or(&empty_hgnc_record);
                tsv_record.clear_refseq_ensembl();

                let anns = &anns_by_gene[&hgnc_id];

                let worst_refseq_ann = anns
                    .iter()
                    .find(|ann| ann.feature_id.starts_with("N") || ann.feature_id.starts_with("X"));
                let worst_ensembl_ann = anns.iter().find(|ann| ann.feature_id.starts_with("ENST"));

                let (refseq_source_ann, ensembl_source_ann) =
                    match (worst_refseq_ann, worst_ensembl_ann) {
                        (Some(refseq), Some(ensembl)) => (Some(refseq), Some(ensembl)),
                        (Some(refseq), None) => (Some(refseq), Some(refseq)),
                        (None, Some(ensembl)) => (Some(ensembl), Some(ensembl)),
                        (None, None) => (None, None),
                    };

                if let Some(ann) = ensembl_source_ann {
                    tsv_record.ensembl_gene_id = Some(hgnc_record.ensembl_gene_id.clone());
                    tsv_record.ensembl_transcript_id = Some(ann.feature_id.clone());
                    tsv_record.ensembl_transcript_coding =
                        Some(ann.feature_biotype.contains(&FeatureBiotype::Coding));
                    tsv_record.ensembl_hgvs_c.clone_from(&ann.hgvs_c);
                    tsv_record.ensembl_hgvs_p.clone_from(&ann.hgvs_p);
                    if !ann.consequences.is_empty() {
                        tsv_record.ensembl_effect = Some(
                            ann.consequences
                                .iter()
                                .map(|c| format!("\"{}\"", &c))
                                .collect(),
                        );
                    }
                    tsv_record.ensembl_exon_dist = ann.distance;
                }

                if let Some(ann) = refseq_source_ann {
                    tsv_record.refseq_gene_id = Some(hgnc_record.entrez_id.clone());
                    tsv_record.refseq_transcript_id = Some(ann.feature_id.clone());
                    tsv_record.refseq_transcript_coding =
                        Some(ann.feature_biotype.contains(&FeatureBiotype::Coding));
                    tsv_record.refseq_hgvs_c.clone_from(&ann.hgvs_c);
                    tsv_record.refseq_hgvs_p.clone_from(&ann.hgvs_p);
                    if !ann.consequences.is_empty() {
                        tsv_record.refseq_effect = Some(
                            ann.consequences
                                .iter()
                                .map(|c| format!("\"{}\"", &c))
                                .collect(),
                        );
                    }
                    tsv_record.refseq_exon_dist = ann.distance;
                }
                writeln!(self.inner, "{}", tsv_record.to_tsv().join("\t"))?;
            }
            Ok(())
        } else {
            writeln!(self.inner, "{}", tsv_record.to_tsv().join("\t"))?;
            Ok(())
        }
    }
}

#[derive(Debug, Default)]
pub struct VarFishSeqvarTsvRecord {
    pub release: String,
    pub chromosome: String,
    pub chromosome_no: u32,
    pub start: usize,
    pub end: usize,
    pub bin: u32,
    pub reference: String,
    pub alternative: String,
    pub var_type: String,
    pub genotype: String,
    pub num_hom_alt: u32,
    pub num_hom_ref: u32,
    pub num_het: u32,
    pub num_hemi_alt: u32,
    pub num_hemi_ref: u32,
    pub in_clinvar: bool,
    pub gnomad_exomes_frequency: f32,
    pub gnomad_exomes_homozygous: i32,
    pub gnomad_exomes_heterozygous: i32,
    pub gnomad_exomes_hemizygous: i32,
    pub gnomad_genomes_frequency: f32,
    pub gnomad_genomes_homozygous: i32,
    pub gnomad_genomes_heterozygous: i32,
    pub gnomad_genomes_hemizygous: i32,
    pub refseq_gene_id: Option<String>,
    pub refseq_transcript_id: Option<String>,
    pub refseq_transcript_coding: Option<bool>,
    pub refseq_hgvs_c: Option<String>,
    pub refseq_hgvs_p: Option<String>,
    pub refseq_effect: Option<Vec<String>>,
    pub refseq_exon_dist: Option<i32>,
    pub ensembl_gene_id: Option<String>,
    pub ensembl_transcript_id: Option<String>,
    pub ensembl_transcript_coding: Option<bool>,
    pub ensembl_hgvs_c: Option<String>,
    pub ensembl_hgvs_p: Option<String>,
    pub ensembl_effect: Option<Vec<String>>,
    pub ensembl_exon_dist: Option<i32>,
}

impl VarFishSeqvarTsvRecord {
    pub fn clear_refseq_ensembl(&mut self) {
        self.refseq_gene_id = None;
        self.refseq_transcript_id = None;
        self.refseq_transcript_coding = None;
        self.refseq_hgvs_c = None;
        self.refseq_hgvs_p = None;
        self.refseq_effect = None;
        self.refseq_exon_dist = None;
        self.ensembl_gene_id = None;
        self.ensembl_transcript_id = None;
        self.ensembl_transcript_coding = None;
        self.ensembl_hgvs_c = None;
        self.ensembl_hgvs_p = None;
        self.ensembl_effect = None;
        self.ensembl_exon_dist = None;
    }

    pub fn to_tsv(&self) -> Vec<String> {
        vec![
            self.release.clone(),
            self.chromosome.clone(),
            format!("{}", self.chromosome_no),
            format!("{}", self.start),
            format!("{}", self.end),
            format!("{}", self.bin),
            self.reference.clone(),
            self.alternative.clone(),
            self.var_type.clone(),
            String::from("."),
            String::from("."),
            String::from("{}"),
            self.genotype.clone(),
            format!("{}", self.num_hom_alt),
            format!("{}", self.num_hom_ref),
            format!("{}", self.num_het),
            format!("{}", self.num_hemi_alt),
            format!("{}", self.num_hemi_ref),
            if self.in_clinvar {
                "TRUE".into()
            } else {
                "FALSE".into()
            },
            String::from("0"),
            String::from("0"),
            String::from("0"),
            String::from("0"),
            String::from("0"),
            String::from("0"),
            String::from("0"),
            String::from("0"),
            format!("{}", self.gnomad_exomes_frequency),
            format!("{}", self.gnomad_exomes_homozygous),
            format!("{}", self.gnomad_exomes_heterozygous),
            format!("{}", self.gnomad_exomes_hemizygous),
            format!("{}", self.gnomad_genomes_frequency),
            format!("{}", self.gnomad_genomes_homozygous),
            format!("{}", self.gnomad_genomes_heterozygous),
            format!("{}", self.gnomad_genomes_hemizygous),
            self.refseq_gene_id.clone().unwrap_or(String::from(".")),
            self.refseq_transcript_id
                .clone()
                .unwrap_or(String::from(".")),
            self.refseq_transcript_coding
                .map(|c| if c { "TRUE".into() } else { "FALSE".into() })
                .unwrap_or(String::from(".")),
            self.refseq_hgvs_c.clone().unwrap_or(String::from(".")),
            self.refseq_hgvs_p.clone().unwrap_or(String::from(".")),
            format!(
                "{{{}}}",
                self.refseq_effect
                    .as_ref()
                    .map(|e| e.join(","))
                    .unwrap_or_default()
            ),
            self.refseq_exon_dist
                .map(|d| format!("{}", d))
                .unwrap_or(String::from(".")),
            self.ensembl_gene_id.clone().unwrap_or(String::from(".")),
            self.ensembl_transcript_id
                .clone()
                .unwrap_or(String::from(".")),
            self.ensembl_transcript_coding
                .map(|c| if c { "TRUE".into() } else { "FALSE".into() })
                .unwrap_or(String::from(".")),
            self.ensembl_hgvs_c.clone().unwrap_or(String::from(".")),
            self.ensembl_hgvs_p.clone().unwrap_or(String::from(".")),
            format!(
                "{{{}}}",
                self.ensembl_effect
                    .as_ref()
                    .map(|e| e.join(","))
                    .unwrap_or_default()
            ),
            self.ensembl_exon_dist
                .map(|d| format!("{}", d))
                .unwrap_or(String::from(".")),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use temp_testdir::TempDir;

    #[tokio::test]
    async fn smoke_test_export_tsv() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("output.tsv");

        let args_common = crate::common::Args {
            verbose: Default::default(),
        };
        let args = Args {
            input: String::from("tests/data/annotate/seqvars/brca2_zar1l/brca2_zar1l.mehari.vcf"),
            pedigree: String::from("tests/data/annotate/seqvars/brca2_zar1l/brca2_zar1l.ped"),
            output: path_out.display().to_string(),
            hgnc: String::from("tests/data/annotate/db/hgnc.tsv"),
            tsv_contig_style: TsvContigStyle::Auto,
        };

        run(&args_common, &args).await?;

        // Standard string comparison against the expected TSV snapshot
        let actual = std::fs::read_to_string(&args.output)?;
        let expected =
            std::fs::read_to_string("tests/data/annotate/seqvars/brca2_zar1l/brca2_zar1l.tsv")?;
        assert_eq!(actual, expected);

        Ok(())
    }
}
