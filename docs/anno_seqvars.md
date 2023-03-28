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

## Command Line Invocation

You can invoke Mehari like this to annotate a VCF file `IN.vcf` to an output file `OUT.vcf` using the built (or downloaded) database as `path/to/db`.

```text
$ mehari annotate seqvars \
    --path-db path/to/db \
    --input-vcf IN.vcf \
    --output-vcf OUT.vcf
```

Note that the input and output files can optionally be gzip/bgzip compressed VCF files with suffixes (`.gz` or `.bgz`) or BCF files with suffix `.bcf`.
The database genome build should match the one in the input VCF file (e.g., both should either be GRCh37/hg19 or GRCh38/hg38).

## Interpreting Annotation Output

The variant effect/consequence will be formatted similar to the one in, following the `ANN` field standard [documented here](https://pcingola.github.io/SnpEff/se_inputoutput/#ann-field-vcf-output-files).
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
