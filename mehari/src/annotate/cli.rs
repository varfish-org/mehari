use crate::annotate::seqvars::consequence::plugin::PluginSpec;
use clap::Args as ClapArgs;
use strum::Display;

#[derive(Debug, Clone, Default, ClapArgs)]
#[group(required = true, multiple = true)]
pub struct Sources {
    /// Transcript database containing the transcript information.
    ///
    /// Pre-built databases are available at https://github.com/varfish-org/mehari-data-tx/releases
    #[arg(long)]
    pub transcripts: Option<Vec<String>>,

    /// Frequency database.
    ///
    /// The frequency database contains gnomAD frequencies for the variants.
    /// There is no pre-built release because of its size, see https://github.com/varfish-org/mehari/blob/main/docs/db_build.md#building-the-frequency-database
    #[arg(long)]
    pub frequencies: Option<Vec<String>>,

    /// ClinVar database.
    ///
    /// The ClinVar database contains clinical significance information for the variants.
    /// Pre-built databases are available at https://github.com/varfish-org/annonars-data-clinvar/releases
    #[arg(long)]
    pub clinvar: Option<Vec<String>>,

    /// CADD database path (built with `annonars cadd import`).
    #[arg(long)]
    pub cadd: Option<Vec<String>>,

    /// SpliceAI database path (built with `annonars spliceai import`).
    #[arg(long)]
    pub spliceai: Option<Vec<String>>,

    /// dbSNP database path (built with `annonars dbsnp import`).
    #[arg(long)]
    pub dbsnp: Option<Vec<String>>,

    /// Custom database paths (in name=path format, e.g. --custom-db my_db=/path/to/db).
    #[arg(long)]
    pub custom_db: Option<Vec<String>>,
}

#[derive(Debug, ClapArgs, Default, Clone)]
pub struct PredictorSettings {
    /// Enable vep compatibility mode.
    /// This enables specific normalization flags, less fine-grained vep consequence terms and
    /// VEP's rules where they differ from HGVS or SO.
    #[arg(long, default_value_t = false, hide = true)]
    pub vep_compatibility_mode: bool,

    #[clap(flatten)]
    pub transcript_settings: TranscriptSettings,

    #[clap(flatten)]
    pub reporting_settings: ReportingSettings,

    #[clap(flatten)]
    pub normalization_settings: NormalizationSettings,

    #[clap(flatten)]
    pub compound_settings: CompoundSettings,
}

impl PredictorSettings {
    pub fn do_not_normalize_variants(&self) -> bool {
        self.normalization_settings.do_not_normalize_variants
    }

    pub fn do_not_renormalize_g(&self) -> bool {
        self.normalization_settings.do_not_renormalize_g
    }

    pub fn vep_consequence_terms(&self) -> bool {
        self.reporting_settings.use_vep_consequence_terms || self.vep_compatibility_mode
    }

    /// `--hgvs-parentheses`, by default `never` in VEP compatibility mode, as VEP writes no
    /// parentheses, and `predicted` otherwise.
    pub fn hgvs_parentheses(&self) -> HgvsParentheses {
        self.reporting_settings
            .hgvs_parentheses
            .unwrap_or(if self.vep_compatibility_mode {
                HgvsParentheses::Never
            } else {
                HgvsParentheses::Predicted
            })
    }
}

#[derive(Debug, ClapArgs, Clone)]
pub struct TranscriptSettings {
    /// Report only transcripts that have one of these tags, comma-separated.
    ///
    /// A gene without such a transcript gets no annotation, unless `--restrict-fallback`
    /// is given.
    #[arg(long, value_delimiter = ',')]
    pub restrict_transcripts: Vec<TagCriterion>,

    /// Report the top-ranked transcript (see `--pick-order`) of a gene that has none of the
    /// `--restrict-transcripts` tags, marked with the message `INFO_RESTRICT_FALLBACK`.
    #[arg(long, requires = "restrict_transcripts")]
    pub restrict_fallback: bool,

    /// Report only one annotation per gene, transcript or allele, chosen by `--pick-order`.
    #[arg(long, value_enum)]
    pub pick: Option<PickUnit>,

    /// How to rank transcripts for `--pick` and `--restrict-fallback`, comma-separated.
    ///
    /// Each criterion only breaks the ties left by the ones before it, and the lowest
    /// accession breaks any remaining tie. A tag also matches its backport, a `-backport`
    /// tag only the backport.
    #[arg(long, value_delimiter = ',', default_values_t = DEFAULT_PICK_ORDER.to_vec())]
    pub pick_order: Vec<PickCriterion>,

    /// Report all annotations and mark the picked ones in the `pick` column, instead of
    /// dropping the others.
    #[arg(long, requires = "pick")]
    pub flag_pick: bool,
}

impl Default for TranscriptSettings {
    fn default() -> Self {
        Self {
            restrict_transcripts: Vec::new(),
            restrict_fallback: false,
            pick: None,
            pick_order: DEFAULT_PICK_ORDER.to_vec(),
            flag_pick: false,
        }
    }
}

#[cfg(test)]
impl TranscriptSettings {
    /// Most severe consequence per `unit`, MANE breaking ties.
    pub(crate) fn most_severe(unit: PickUnit) -> Self {
        Self {
            pick: Some(unit),
            pick_order: vec![
                PickCriterion::Severity,
                PickCriterion::ManeSelect,
                PickCriterion::ManePlusClinical,
            ],
            ..Default::default()
        }
    }
}

#[derive(Debug, ClapArgs, Default, Clone)]
pub struct ReportingSettings {
    /// Whether to keep intergenic variants.
    #[arg(long, default_value_t = false)]
    pub keep_intergenic: bool,

    /// Whether to report splice variants in UTRs.
    #[arg(long, default_value_t = false)]
    pub discard_utr_splice_variants: bool,

    /// Whether to use less fine-grained VEP consequence terms.
    #[arg(long, default_value_t = false, hide = true)]
    use_vep_consequence_terms: bool,

    /// Parentheses around predicted changes, e.g. `p.(Gly2Arg)`, as HGVS recommends. Default
    /// `predicted`, `never` with `--vep-compatibility-mode`.
    #[arg(long, value_enum)]
    pub hgvs_parentheses: Option<HgvsParentheses>,

    /// Plugin to run, adding ANN columns: tx-sequence-ref, tx-sequence-alt, aa-sequence-ref,
    /// aa-sequence-alt, nmd, or the path of a WebAssembly component (`.wasm`) followed by its
    /// options as `,key=value`; `data=<dir>` lets the component read that directory.
    #[arg(long = "plugin")]
    pub plugins: Vec<PluginSpec>,
}

#[derive(Debug, ClapArgs, Default, Clone)]
pub struct NormalizationSettings {
    /// Whether to do hgvs shifting for hgvs.g like vep does
    #[arg(long, default_value_t = false, hide = true)]
    vep_hgvs_shift: bool,

    /// Whether to skip HGVS normalization.
    #[arg(long, default_value_t = false, hide = true)]
    do_not_normalize_variants: bool,

    /// Whether to skip re-normalizing genomic variants.
    #[arg(long, default_value_t = false, hide = true)]
    do_not_renormalize_g: bool,
}

/// Value of `--hgvs-parentheses`.
#[derive(Debug, Copy, Clone, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum HgvsParentheses {
    /// Around predicted changes, i.e. every `p.` except `p.?` and `p.Met1?`.
    #[default]
    Predicted,
    /// Never.
    Never,
}

/// Unit for `--pick`.
#[derive(Debug, Copy, Clone, PartialEq, Eq, clap::ValueEnum)]
pub enum PickUnit {
    /// One annotation per gene.
    Gene,
    /// One annotation per transcript, for transcripts with several annotations.
    Transcript,
    /// One annotation per allele.
    Allele,
}

/// Transcript tag for `--restrict-transcripts`.
#[derive(Debug, Copy, Clone, PartialEq, Eq, clap::ValueEnum)]
pub enum TagCriterion {
    ManeSelect,
    ManeSelectBackport,
    ManePlusClinical,
    ManePlusClinicalBackport,
    RefseqSelect,
    RefseqSelectBackport,
    EnsemblCanonical,
    EnsemblCanonicalBackport,
    GencodePrimary,
    GencodePrimaryBackport,
    Basic,
    BasicBackport,
}

/// Criterion for `--pick-order`.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Display, clap::ValueEnum)]
#[strum(serialize_all = "kebab-case")]
pub enum PickCriterion {
    ManeSelect,
    ManeSelectBackport,
    ManePlusClinical,
    ManePlusClinicalBackport,
    RefseqSelect,
    RefseqSelectBackport,
    EnsemblCanonical,
    EnsemblCanonicalBackport,
    GencodePrimary,
    GencodePrimaryBackport,
    Basic,
    BasicBackport,
    /// Protein-coding transcript.
    Coding,
    /// No soft-filter reason from the database build, e.g. no missing stop codon.
    Unflagged,
    /// Most severe consequence.
    Severity,
    /// Longest transcript.
    Length,
}

/// Default for `--pick-order`.
pub const DEFAULT_PICK_ORDER: &[PickCriterion] = &[
    PickCriterion::ManeSelect,
    PickCriterion::ManePlusClinical,
    PickCriterion::RefseqSelect,
    PickCriterion::EnsemblCanonical,
    PickCriterion::GencodePrimary,
    PickCriterion::Basic,
    PickCriterion::Coding,
    PickCriterion::Unflagged,
    PickCriterion::Severity,
    PickCriterion::Length,
];

#[derive(Debug, clap::Args, Default, Clone)]
pub struct CompoundSettings {
    /// Experimental: Enable variant grouping to evaluate the compound effect of multiple variants on the same transcript.
    /// When disabled, Mehari evaluates each variant independently.
    #[arg(long, default_value_t = false)]
    pub enable_compound_variants: bool,

    /// Experimental: The strategy used to evaluate grouped variants for compound effects.
    #[arg(long, value_enum, default_value_t = PhasingStrategy::Strict)]
    pub phasing_strategy: PhasingStrategy,

    /// Experimental: Maximum distance in bp between neighbouring variants of a compound group.
    /// Unlimited if not set; a good value for peptide analyses is 33.
    #[arg(long)]
    pub compound_window_size: Option<u32>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, clap::ValueEnum)]
pub enum PhasingStrategy {
    /// Variants are only grouped if explicitly phased ('|') and sharing a Phase Set (PS).
    /// Unphased variants are evaluated independently.
    #[default]
    Strict,

    /// Respects explicit phasing, but treats homozygous variants as universally phased
    /// Unphased heterozygous variants remain independent.
    Relaxed,

    /// Completely ignores phasing metadata and _assumes_ all variants on the transcript are on the same haplotype.
    Ignore,
}

#[cfg(test)]
mod test {
    use super::{HgvsParentheses, PredictorSettings};

    /// `--hgvs-parentheses` defaults to `never` in VEP compatibility mode, as VEP writes none.
    #[rstest::rstest]
    #[case(None, false, HgvsParentheses::Predicted)]
    #[case(None, true, HgvsParentheses::Never)]
    #[case(Some(HgvsParentheses::Predicted), true, HgvsParentheses::Predicted)]
    #[case(Some(HgvsParentheses::Never), false, HgvsParentheses::Never)]
    fn hgvs_parentheses_default(
        #[case] option: Option<HgvsParentheses>,
        #[case] vep_compatibility_mode: bool,
        #[case] expected: HgvsParentheses,
    ) {
        let mut settings = PredictorSettings {
            vep_compatibility_mode,
            ..Default::default()
        };
        settings.reporting_settings.hgvs_parentheses = option;
        assert_eq!(settings.hgvs_parentheses(), expected);
    }
}
