Annotating sequence variants (SNV/indel/MNV).

## Sequence Variants And Their VCF Representation

Sequence variants are small variants that are typically represented by the substitution of one or more bases in the reference genome.
In a VCF file, this will typically look as follows:

```text
#CHROM  POS       ID  REF  ALT  QUAL  FILTER  INFO
17      41197700  .   A    C    .     .       .
17      41197708  .   TG   T    .     .       .
```

In the literature, there generally is made a distinction between sequence (or: small) variants and structural variants.
There, small variants are those affecting up to 50bp while structural variants are those that are larger.
Mehari does not enforce such a restriction but attempt to annotate all variants based on their VCF representation (chromosome, 1-based position, reference base and alternative bases).
You can find the VCF variant specification [here on Github](https://samtools.github.io/hts-specs/).

## Sequence Variant Annotation

Currently, Mehari will annotate variants using:

- The predicted impact on gene transcripts and the corresponding protein sequence (in the case of coding genes).
- Their frequency in the gnomAD exomes and genomes databases as well as the HelixMtDb database in the case of mitochondrial databases.
- Variant information from ClinVar, if any

## Command Line Invocation

You can invoke Mehari to annotate a VCF file `IN.vcf` creating an output file `OUT.vcf` using the built (or downloaded) databases – for example the transcript database – as follows:

```text
$ mehari annotate seqvars \
    --reference path/to/reference.fa \
    --transcripts path/to/transcripts-db \
    --input IN.vcf \
    --output OUT.vcf
```

Note that the input and output files can optionally be gzip/bgzip compressed VCF files with suffixes (`.gz` or `.bgz`) or BCF files with suffix `.bcf`.
The database assembly should match the one in the input VCF file (e.g., both should either be GRCh37/hg19 or GRCh38/hg38).
You can optionally define the assembly via `--assembly`.
Output formats can be toggled using `--output-format` (e.g., `vcf` or `jsonl`), and you can control the number of threads used for parallel processing by using the `--threads` flag (we recommend a value of around 16).

### Plugins

Plugins add columns to the `ANN` field; enable each with `--plugin`, in the order the columns should appear.
Built-in plugins are `tx-sequence-ref`, `tx-sequence-alt`, `aa-sequence-ref` and `aa-sequence-alt` (transcript and protein sequences without and with the variant), and `nmd` (whether stop-gained and frameshift variants trigger nonsense-mediated decay).
`--plugin` also takes the path of a WebAssembly component, for which you can write your own plugins in Python, Rust and other languages; see [Plugins](`self::plugins`).

### Transcript Selection

By default, Mehari reports one annotation per transcript near the variant, including transcripts up to 5 kbp upstream or downstream.
Four options narrow this down:

- `--restrict-transcripts mane-select,mane-plus-clinical` reports only transcripts with one of these tags.
  A gene without such a transcript gets no annotation, and Mehari logs at the end how many genes this affected.
  With `--restrict-fallback`, such a gene keeps its top-ranked transcript instead, marked with the message `INFO_RESTRICT_FALLBACK`.
- `--pick gene`, `--pick transcript` or `--pick allele` reports one annotation per gene, per transcript or per allele.
  `transcript` only matters when a transcript appears more than once near the variant.
- `--pick-order` ranks the candidates for `--pick` and `--restrict-fallback`.
  Each criterion only breaks the ties left by the ones before it, and the lowest accession breaks any remaining tie.
  The default is `mane-select,mane-plus-clinical,refseq-select,ensembl-canonical,gencode-primary,basic,coding,unflagged,severity,length`.
  A tag also matches its backport to GRCh37, while `mane-select-backport` and the other `-backport` values match only the backport.
  Criteria before `severity` need no annotation, so `--pick gene --pick-order mane-select` annotates only the picked transcripts.
- `--flag-pick` keeps all annotations and marks the picked ones with `1` in an extra ANN column, `pick`, after the columns of the plugins.

These options replace `--pick-transcript`, `--pick-transcript-mode` and `--report-most-severe-consequence-by`:

| Before | Now |
|---|---|
| `--pick-transcript mane-select --pick-transcript-mode all` | `--restrict-transcripts mane-select --restrict-fallback` |
| `--pick-transcript mane-select --pick-transcript length --pick-transcript-mode first` | `--pick gene --pick-order mane-select,length` |
| `--report-most-severe-consequence-by gene` | `--pick gene --pick-order severity,mane-select,mane-plus-clinical` |
| `--report-most-severe-consequence-by transcript` | `--pick transcript --pick-order severity` |
| `--report-most-severe-consequence-by allele` | `--pick allele --pick-order severity,mane-select,mane-plus-clinical` |

### Experimental: Compound Variants

Mehari includes an experimental feature that allows you to evaluate the compound effect of multiple variants on the same transcript.
You can enable this behavior by passing the `--enable-compound-variants` flag.
The evaluation strategy can be managed via `--phasing-strategy` (options include `strict` [default], `relaxed`, or `ignore`).
`--compound-window-size` sets the maximum distance in bp between neighbouring variants of a group (unlimited if not set; a good value for peptide analyses is 33).

## Interpreting Annotation Output

The variant effect/consequence will be formatted similar to the one in, following the `ANN` field standard [documented here](https://pcingola.github.io/SnpEff/se_inputoutput/#ann-field-vcf-output-files).

`Gene_ID` is the gene ID of the transcript's gene record in the annotation, with its version if it has one, e.g. `672` for BRCA1 in RefSeq or `ENSG00000182378.15` for PLCXD1 in GENCODE.
Mehari refuses a transcript database that predates annotation gene IDs, i.e., one in which no transcript has a gene ID.
Rebuild such a database with the current mehari version.

The population frequency will be written to the VCF INFO fields as follows:

- gnomAD exomes
    - `gnomad_exomes_an` -- number of observed alleles in gnomAD exomes
    - `gnomad_exomes_hom` -- number of homozygous carriers in gnomAD exomes
    - `gnomad_exomes_het` -- number of heterozygous carriers in gnomAD exomes
    - `gnomad_exomes_hemi` -- number of hemizygous carriers in gnomAD exomes
- gnomAD genomes
    - `gnomad_genomes_an` -- number of observed alleles in gnomAD genomes
    - `gnomad_genomes_hom` -- number of homozygous carriers in gnomAD genomes
    - `gnomad_genomes_het` -- number of heterozygous carriers in gnomAD genomes
    - `gnomad_genomes_hemi` -- number of hemizygous carriers in gnomAD genomes
- gnomAD mtDNA
    - `gnomad_mtdna_an` -- number of individuals with sufficient coverage in gnomAD mtDNA database
    - `gnomad_mtdna_hom` -- number of homoplasmic carriers in gnomAD mtDNA database
    - `gnomad_mtdna_het` -- number of heteroplasmic carriers in gnomAD mtDNA database
- HelixMtDb
    - `helix_an` -- number of individuals in HelixMtDb
    - `helix_hom` -- number of homoplasmic carriers in HelixMtDb
    - `helix_het` -- number of heteroplasmic carriers in HelixMtDb

From these integer numbers, the relative allele frequencies can be derived as needed.
