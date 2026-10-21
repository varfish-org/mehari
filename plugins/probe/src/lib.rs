//! Test plugin: reports what every function of `mehari/wit/plugin.wit` returns, as `key=value`
//! pairs separated by `;` in its single column `probe`.
//!
//! Options:
//! - `fail-init=<message>`: `init` fails with the message.
//! - `init-read=<file>`: `init` reads the file in the data directory, or fails.
//! - `fail=<message>`: `annotate` fails with the message.
//! - `bad-index`: asks for an annotation past the last one.
//! - `field=<column>`: reports what an earlier plugin wrote into the column.
//! - `read=<file>`: reports the content of the file in the data directory.
//! - `write=<file>`: reports whether the file in the data directory could be written.

use std::cell::RefCell;
use std::fmt::Display;

wit_bindgen::generate!({
    path: "../../mehari/wit",
    world: "plugin",
});

use mehari::plugin::host;
use mehari::plugin::types::{Position, Rank};

thread_local! {
    static OPTIONS: RefCell<Vec<(String, String)>> = const { RefCell::new(Vec::new()) };
}

struct Probe;

impl Guest for Probe {
    fn init(options: Vec<(String, String)>) -> Result<Vec<String>, String> {
        if let Some((_, message)) = options.iter().find(|(key, _)| key == "fail-init") {
            return Err(message.clone());
        }
        if let Some((_, file)) = options.iter().find(|(key, _)| key == "init-read") {
            std::fs::read(format!("/data/{file}")).map_err(|e| e.to_string())?;
        }
        OPTIONS.set(options);
        Ok(vec!["probe".into()])
    }

    fn annotate(
        variants: Vec<VcfVariant>,
        annotations: &Annotations,
    ) -> Result<Vec<Vec<Option<String>>>, String> {
        let options = OPTIONS.with_borrow(Clone::clone);
        let option = |key: &str| {
            options
                .iter()
                .find(|(k, _)| k == key)
                .map(|(_, value)| value.clone())
        };
        if let Some(message) = option("fail") {
            return Err(message);
        }
        if option("bad-index").is_some() {
            annotations.feature_id(annotations.len());
        }

        let rows = (0..annotations.len()).map(|index| {
            let a = annotations;
            let mut pairs = vec![
                ("assembly", host::assembly()),
                ("version", host::version()),
                ("vep", host::vep_compatible().to_string()),
                ("variants", variants.len().to_string()),
                ("feature-id", a.feature_id(index)),
                ("feature-type", a.feature_type(index)),
                ("consequences", a.consequences(index).join("&")),
                ("impact", a.impact(index)),
                ("gene-symbol", a.gene_symbol(index)),
                ("gene-id", a.gene_id(index)),
                ("biotypes", a.biotypes(index).join("&")),
                ("tags", a.tags(index).join("&")),
                ("rank", show(a.rank(index).map(rank))),
                ("hgvs-g", show(a.hgvs_g(index))),
                ("hgvs-n", show(a.hgvs_n(index))),
                ("hgvs-c", show(a.hgvs_c(index))),
                ("hgvs-p", show(a.hgvs_p(index))),
                ("cdna-pos", show(a.cdna_pos(index).map(position))),
                ("cds-pos", show(a.cds_pos(index).map(position))),
                ("protein-pos", show(a.protein_pos(index).map(position))),
                ("distance", show(a.distance(index))),
                ("warnings", a.warnings(index).join("&")),
                ("cds-end", show(a.cds_end(index))),
                ("frameshift", show(a.frameshift(index))),
                ("first-changed-codon", show(a.first_changed_codon(index))),
            ];
            if let Some(column) = option("field") {
                pairs.push(("field", show(a.field(index, &column))));
            }
            if let Some(tx) = a.transcript(index) {
                let exon = tx.exons.first().map(|exon| {
                    format!(
                        "{}-{}:{}-{}:{}",
                        exon.tx_start, exon.tx_end, exon.genome_start, exon.genome_end, exon.cigar
                    )
                });
                pairs.extend([
                    ("tx-id", tx.id),
                    ("tx-protein", show(tx.protein)),
                    ("tx-contig", tx.contig),
                    ("tx-build", tx.genome_build),
                    ("tx-strand", format!("{:?}", tx.strand)),
                    (
                        "tx-genome-cds",
                        format!("{}-{}", show(tx.genome_cds_start), show(tx.genome_cds_end)),
                    ),
                    ("tx-exons", tx.exons.len().to_string()),
                    ("tx-first-exon", show(exon)),
                ]);
            }
            for (key, sequences) in [("ref", a.reference(index)), ("alt", a.alternative(index))] {
                if let Some(seq) = sequences {
                    let value = format!(
                        "{}/{}/{}",
                        seq.transcript.len(),
                        seq.protein.len(),
                        seq.cds_start
                    );
                    pairs.push((key, value));
                }
            }
            if let Some(translation) = a.translation(index) {
                let exceptions = translation
                    .exceptions
                    .iter()
                    .map(|exception| format!("{}{}", exception.position, exception.amino_acid))
                    .collect::<Vec<_>>();
                let value = format!("{}/{}", translation.table, exceptions.join("&"));
                pairs.push(("translation", value));
            }
            if let Some(file) = option("read") {
                let content = std::fs::read_to_string(format!("/data/{file}"));
                pairs.push(("read", content.unwrap_or_else(|e| format!("error {e}"))));
            }
            if let Some(file) = option("write") {
                let written = std::fs::write(format!("/data/{file}"), "x").is_ok();
                pairs.push(("write", written.to_string()));
            }
            let probe = pairs
                .iter()
                .map(|(key, value)| format!("{key}={value}"))
                .collect::<Vec<_>>()
                .join(";");
            vec![Some(probe)]
        });
        Ok(rows.collect())
    }
}

fn show(value: Option<impl Display>) -> String {
    value.map(|value| value.to_string()).unwrap_or_default()
}

fn rank(rank: Rank) -> String {
    format!("{}/{}", rank.ord, rank.total)
}

fn position(pos: Position) -> String {
    format!("{}/{}", pos.ord, show(pos.total))
}

export!(Probe);
