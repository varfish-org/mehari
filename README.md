[![CI](https://github.com/bihealth/mehari/actions/workflows/rust.yml/badge.svg)](https://github.com/bihealth/mehari/actions/workflows/rust.yml)
[![codecov](https://codecov.io/gh/bihealth/mehari/branch/main/graph/badge.svg?token=B1dfb7N2n8)](https://codecov.io/gh/bihealth/mehari)

# Mehari

<img align="right" width="200" height="200" src="misc/camel.jpeg">

Mehari is a software package for annotating VCF files with variant effect.
The program uses [hgvs-rs](https://crates.io/crates/hgvs) for projecting genomic variants to transcripts and proteins and thus has high prediction quality.

## Supported Sequence Variant Frequency Databases

Mehari can import public sequence variant frequency databases.
The supported set slightly differs between import for GRCh37 and GRCh38.

**GRCh37**

- gnomAD r2.1.1 Exomes [`gnomad.exomes.r2.1.1.sites.vcf.bgz`](https://gnomad.broadinstitute.org/downloads#v2)
- gnomAD r2.1.1 Genomes [`gnomad.genomes.r2.1.1.sites.vcf.bgz`](https://gnomad.broadinstitute.org/downloads#v2)
- gnomAD v3.1 mtDNA [`gnomad.genomes.v3.1.sites.chrM.vcf.bgz`](https://gnomad.broadinstitute.org/downloads#v3-mitochondrial-dna)
- HelixMTdb `HelixMTdb_20200327.tsv`

**GRCh38**

- gnomAD r2.1.1 lift-over Exomes [`gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz`](https://gnomad.broadinstitute.org/downloads#v2)
- gnomAD v3.1 Genomes [`gnomad.genomes.v3.1.2.sites.$CHROM.vcf.bgz`](https://gnomad.broadinstitute.org/downloads#v3)
- gnomAD v3.1 mtDNA [`gnomad.genomes.v3.1.sites.chrM.vcf.bgz`](https://gnomad.broadinstitute.org/downloads#v3-mitochondrial-dna)
- HelixMTdb `HelixMTdb_20200327.tsv`

## Internal Notes

```
rm -rf /tmp/out ; cargo run -- db create seqvar-freqs --path-output-db /tmp/out --genome-release grch38 --path-helix-mtdb ~/Downloads/HelixMTdb_20200327.vcf.gz --path-gnomad-mtdna ~/Downloads/gnomad.genomes.v3.1.sites.chrM.vcf.bgz --path-gnomad-exomes-xy tests/data/db/create/seqvar_freqs/xy-38/gnomad.exomes.r2.1.1.sites.chrX.vcf --path-gnomad-exomes-xy tests/data/db/create/seqvar_freqs/xy-38/gnomad.exomes.r2.1.1.sites.chrY.vcf --path-gnomad-genomes-xy tests/data/db/create/seqvar_freqs/xy-38/gnomad.genomes.r3.1.1.sites.chrX.vcf --path-gnomad-genomes-xy tests/data/db/create/seqvar_freqs/xy-38/gnomad.genomes.r3.1.1.sites.chrY.vcf

rm -rf /tmp/out ; cargo run -- db create seqvar-freqs --path-output-db /tmp/out --genome-release grch37 --path-gnomad-mtdna ~/Downloads/gnomad.genomes.v3.1.sites.chrM.vcf.bgz --path-gnomad-exomes-xy tests/data/db/create/seqvar_freqs/xy-37/gnomad.exomes.r2.1.1.sites.chrX.vcf --path-gnomad-exomes-xy tests/data/db/create/seqvar_freqs/xy-37/gnomad.exomes.r2.1.1.sites.chrY.vcf --path-gnomad-genomes-xy tests/data/db/create/seqvar_freqs/xy-37/gnomad.genomes.r2.1.1.sites.chrX.vcf
```
