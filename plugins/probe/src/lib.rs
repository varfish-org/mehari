//! Test plugin: reports what every function of `mehari/wit/plugin.wit` returns, as `key^value`
//! pairs separated by `~`, with lists joined by `+`, in its column `probe`. Also fills a column
//! of every other type, e.g. `probe_int` and `probe_ints`.
//!
//! Options:
//! - `fail-init=<message>`: `init` fails with the message.
//! - `init-read=<file>`: `init` reads the file in the data directory, or fails.
//! - `column=<name>`: `init` also declares a string column with that name.
//! - `fail=<message>`: `annotate` fails with the message.
//! - `stale`: keeps the first annotation and reads it again in the next call.
//! - `set=<case>`: sets a value that does not fit: `separator`, `ampersand`, `list-element`,
//!   `undeclared`, `wrong-type`, `list-as-single`, `nan` or `nan-list`.
//! - `empty-list`: sets `probe_strings` to an empty list.
//! - `field=<column>`: reports the value of the column, e.g. from an earlier plugin.
//! - `read=<file>`: reports the content of the file in the data directory.
//! - `write=<file>`: reports whether the file in the data directory could be written.

use std::cell::RefCell;
use std::fmt::Display;

wit_bindgen::generate!({
    path: "../../mehari/wit",
    world: "plugin",
});

use exports::mehari::plugin::transcript_plugin::Guest;
use mehari::plugin::host;
use mehari::plugin::types::{Annotation, Position, Prediction, Rank, Setup};

thread_local! {
    static OPTIONS: RefCell<Vec<(String, String)>> = const { RefCell::new(Vec::new()) };
    static KEPT: RefCell<Option<Annotation>> = const { RefCell::new(None) };
}

fn option(key: &str) -> Option<String> {
    OPTIONS.with_borrow(|options| {
        options
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, value)| value.clone())
    })
}

struct Probe;

impl Guest for Probe {
    fn init(setup: &Setup) -> Result<(), String> {
        OPTIONS.set(setup.options());
        if let Some(message) = option("fail-init") {
            return Err(message);
        }
        if let Some(file) = option("init-read") {
            std::fs::read(format!("/data/{file}")).map_err(|e| e.to_string())?;
        }
        setup.add_string_column("probe");
        setup.add_integer_column("probe_int");
        setup.add_float_column("probe_float");
        setup.add_boolean_column("probe_bool");
        setup.add_string_list_column("probe_strings");
        setup.add_integer_list_column("probe_ints");
        setup.add_float_list_column("probe_floats");
        setup.add_boolean_list_column("probe_bools");
        if let Some(name) = option("column") {
            setup.add_string_column(&name);
        }
        Ok(())
    }

    fn annotate(prediction: &Prediction) -> Result<(), String> {
        if let Some(message) = option("fail") {
            return Err(message);
        }
        let annotations = prediction.annotations();
        if option("stale").is_some() {
            // Read the annotation kept from the previous call, then keep this call's first.
            KEPT.with_borrow(|kept| kept.as_ref().map(Annotation::feature_id));
            KEPT.set(prediction.annotations().into_iter().next());
        }
        for (index, a) in (0..).zip(&annotations) {
            match option("set").as_deref() {
                Some("separator") => a.set_string("probe", "a|b"),
                Some("ampersand") => a.set_string("probe", "a&b"),
                Some("list-element") => a.set_string_list("probe_strings", &["a|b".into()]),
                Some("undeclared") => a.set_string("undeclared", "x"),
                Some("wrong-type") => a.set_integer("probe", 1),
                Some("list-as-single") => a.set_string("probe_strings", "x"),
                Some("nan") => a.set_float("probe_float", f64::NAN),
                Some("nan-list") => a.set_float_list("probe_floats", &[0.5, f64::NAN]),
                _ => {}
            }
            a.set_string("probe", &report(prediction, a));
            a.set_integer("probe_int", -index - 1);
            a.set_float("probe_float", 0.25);
            a.set_boolean("probe_bool", index == 0);
            if option("empty-list").is_some() {
                a.set_string_list("probe_strings", &[]);
            } else {
                a.set_string_list("probe_strings", &["a".into(), "b".into()]);
            }
            a.set_integer_list("probe_ints", &[1, -2, 3]);
            a.set_float_list("probe_floats", &[0.5, 1.5]);
            a.set_boolean_list("probe_bools", &[true, false]);
        }
        Ok(())
    }
}

fn report(prediction: &Prediction, a: &Annotation) -> String {
    let exon = a.exons().first().map(|exon| {
        format!(
            "{}-{}:{}-{}:{}",
            exon.tx_start, exon.tx_end, exon.genome_start, exon.genome_end, exon.cigar
        )
    });
    let variants = prediction.variants();
    let len = |sequence: Option<String>| sequence.map(|s| s.len());
    let exceptions = a
        .translation_exceptions()
        .iter()
        .map(|exception| format!("{}{}", exception.position, exception.amino_acid))
        .collect::<Vec<_>>();
    let mut pairs = vec![
        ("assembly", host::assembly()),
        ("version", host::version()),
        ("vep", host::vep_compatible().to_string()),
        ("variants", variants.len().to_string()),
        (
            "variant-position",
            show(variants.first().map(|var| var.position)),
        ),
        ("feature-id", a.feature_id()),
        ("feature-type", a.feature_type()),
        ("consequences", a.consequences().join("+")),
        ("impact", a.impact()),
        ("gene-symbol", a.gene_symbol()),
        ("gene-id", a.gene_id()),
        ("biotypes", a.biotypes().join("+")),
        ("tags", a.tags().join("+")),
        ("rank", show(a.rank().map(rank))),
        ("hgvs-g", show(a.hgvs_g())),
        ("hgvs-n", show(a.hgvs_n())),
        ("hgvs-c", show(a.hgvs_c())),
        ("hgvs-p", show(a.hgvs_p())),
        ("cdna-pos", show(a.cdna_pos().map(position))),
        ("cds-pos", show(a.cds_pos().map(position))),
        ("protein-pos", show(a.protein_pos().map(position))),
        ("distance", show(a.distance())),
        ("warnings", a.warnings().join("+")),
        ("protein-id", show(a.protein_id())),
        ("contig", show(a.contig())),
        ("genome-build", show(a.genome_build())),
        (
            "strand",
            show(a.strand().map(|strand| format!("{strand:?}"))),
        ),
        ("exons", a.exons().len().to_string()),
        ("first-exon", show(exon)),
        (
            "genome-cds",
            format!(
                "{}-{}",
                show(a.genome_cds_start()),
                show(a.genome_cds_end())
            ),
        ),
        ("reference-transcript", show(len(a.reference_transcript()))),
        ("reference-protein", show(len(a.reference_protein()))),
        (
            "reference-cds",
            format!(
                "{}-{}",
                show(a.reference_cds_start()),
                show(a.reference_cds_end())
            ),
        ),
        (
            "alternative-transcript",
            show(len(a.alternative_transcript())),
        ),
        ("alternative-protein", show(len(a.alternative_protein()))),
        ("alternative-cds-start", show(a.alternative_cds_start())),
        ("translation-table", show(a.translation_table())),
        ("translation-exceptions", exceptions.join("+")),
        ("frameshift", show(a.frameshift())),
        ("first-changed-codon", show(a.first_changed_codon())),
    ];
    if let Some(column) = option("field") {
        pairs.push(("field", show(a.field(&column))));
    }
    if let Some(file) = option("read") {
        let content = std::fs::read_to_string(format!("/data/{file}"));
        pairs.push(("read", content.unwrap_or_else(|e| format!("error {e}"))));
    }
    if let Some(file) = option("write") {
        let written = std::fs::write(format!("/data/{file}"), "x").is_ok();
        pairs.push(("write", written.to_string()));
    }
    pairs
        .iter()
        .map(|(key, value)| format!("{key}^{}", clean(value)))
        .collect::<Vec<_>>()
        .join("~")
}

/// Replace what VCF and `ANN` use as separators.
fn clean(value: &str) -> String {
    value.replace(
        |c: char| c.is_whitespace() || matches!(c, '|' | ',' | ';' | '=' | '&'),
        "_",
    )
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
