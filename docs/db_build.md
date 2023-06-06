Building databases.

This chapter describes how to build the databases needed to run Mehari.

Note that we aim to provide prebuilt databases with Zenodo later on so you don't have to.

# Frequency Database Files

- HelixMtDb can be downloaded from [helix.com](https://www.helix.com/pages/mitochondrial-variant-database)
- gnomAD has a [downloage page](https://gnomad.broadinstitute.org/downloads) where you can find download URLs, e.g., from Microsoft Azure

## Prepare HelixMtDb

You can convert the HelixMtDb TSV file to VCF with the `misc/helix-to-vcf.py` script from the Mehari repository.
You will first have to install a few dependencies.
For this, we create a new conda environment using `mamba`.

```text
$ mamba create -y -n helix-to-vcf python=3.10 vcfpy htslib
$ cat HelixMTdb_20200327.tsv \
  | python misc/helix-to-vcf.py \
  | bgzip -c \
  > HelixMTdb_20200327.tsv.gz
$ tabix -f HelixMTdb_20200327.tsv.gz
```

## Building the Frequency Database

This is done with [annona-rs](https://github.com/bihealth/annona-rs).
The `annonars` crate is a Rust create that ships with a binary for building genome annotation databases as RocksDB databases.
You can install it using `cargo install annonars` or Bioconda (`conda install -c bioconda annonars`).
The `mehari` crate links to the `annonars` library for later accessing the data.
The main advantage is centralized maintanence of the RocksDB related code and the ability for fast import.

```text
# annonars freqs import \
    --genome-release grch38 \
    --path-out-rocksdb annonars-freqs-grch38 \
    --gnomad-genomes-version 3.1.2 \
    --gnomad-exomes-version 2.1.1 \
    --gnomad-mtdna-version 3.1 \
    --helixmtdb-version 20200327 \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr2.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr3.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr4.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr5.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr6.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr7.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr8.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr9.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr10.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr11.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr12.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr13.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr14.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr15.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr16.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr17.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr18.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr19.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr20.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr21.vcf.bgz \
    --path-gnomad-genomes-auto annos/grch38/gnomad_genomes/download/gnomad.genomes.v3.1.2.sites.chr22.vcf.bgz \
    --gnomad-exomes-xy annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.X.vcf.bgz \
    --gnomad-exomes-xy annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.Y.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.1.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.10.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.11.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.12.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.13.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.14.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.15.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.16.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.17.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.18.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.19.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.2.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.20.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.21.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.22.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.3.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.4.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.5.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.6.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.7.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.8.liftover_grch38.vcf.bgz \
    --gnomad-exomes-auto annos/grch38/gnomad_exomes/download/gnomad.exomes.r2.1.1.sites.9.liftover_grch38.vcf.bgz \
    --path-gnomad-genomes-xy annos/grch38/gnomad_genomes/download/gnomad.genomes.r2.1.1.sites.X.vcf.bgz \
    --gnomad-mtdna annos/grch38/gnomad_mtdna/gnomad_mtdna.vcf.gz \
    --path-helixmtdb annos/grch38/helixmtdb/helixmtdb.vcf.gz
```

### Optional: Strip gnomAD VCF Files Before Import

You can strip greatly reduce the nuclear variant files using the following [bcftools](https://samtools.github.io/bcftools/bcftools.html) command line:

```text
$ bcftools \
    annotate \
    --threads 4 \
    -x QUAL,FILTER,ID,^INFO/nhomalt,INFO/nhomalt_female,INFO/nhomalt_male,INFO/nhomalt_XX,INFO/nhomalt_XY,INFO/nonpar,INFO/AN,INFO/AC,INFO/AC_het,INFO/AC_hom,INFO/AC_female,INFO/AC_male,INFO/AC_XX,INFO/AC_XY \
    -O z \
    IN.vcf \
    OUT.vcf.gz
```

The reduction of this is quite big.
For GRCh37, the gnomAD raw download files go down from 449GB to 2.8GB and for GRCh38 from 2.2TB to 4.2GB.

# Transcript Database Files

You can download the gzip-ed JSON files from [SACGF/cdot](https://github.com/SACGF/cdot/releases/tag/v0.2.14).

## Creating SeqRepo

You will have to create a [seqrepo](https://github.com/biocommons/biocommons.seqrepo) repository that fits to the cdot files.
Since NCBI does not store old releases, your best bet is to get the latest cdot JSON and a recent copy of NCBI and ENSEMBL transcripts.
Note that cdot also provides a way to [load transcript FASTA](https://github.com/SACGF/cdot/wiki/FastaSeqFetcher) from the genome but in general, the NCBI transcript sequence differs from the genome.

First, create a conda environment for importing seqrepo:

```text
$ mamba create -y -n seqrepo python=3.8
$ conda activate seqrepo
```

Then, download the transcript FASTA files into the current directory:

```text
$ wget https://ftp.ensembl.org/pub/grch37/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz
$ wget https://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.{1..12}.rna.fna.gz \
    https://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.files.installed
```

Now, initialize a new seqrepo from the biocommons server:

```text
$ mkdir -p seqrepo-data
$ export SEQREPO_ROOT_DIR=$PWD/seqrepo-data
$ seqrepo init -i master
```

And then load the FASTA files we downloaded above:

```
$ seqrepo load --instance-name master --namespace ENSEMBL \
    Homo_sapiens.GRCh37.cdna.all.fa.gz
$ seqrepo load --instance-name master --namespace RefSeq \
            human.{1..12}.rna.fna.gz
```

Now you can create the transcript database (see below).
For some transcripts, no sequence will be available and they will be skipped.
This will be OK as there will be a more recent version available.

## Building Transcript Database

You can build the transcript database protocolbuffers binary using the following command:

```text
$ mehari db create txs \
    --path-out output/db/grch37/txs.bin.zst \
    \
    --path-seqrepo-instance path/to/seqrepo-data/master \
    \
    --path-cdot-json cdot-0.2.12.refseq.grch37_grch38.json \
    --path-cdot-json cdot-0.2.12.ensembl.grch37_grch38.json \
    \
    --path-seqrepo-instance path/to/seqrepo-data/master \
    \
    --genome-release grch37
```

You will have to build the transcript database for each genome release that you want and manually specify the release to `--genome-release`.
For GRCh38, simply use `--genome-release grch38`.

You can enable compression by using the suffix `.gz` for gzip compression and `.zst` for zstandard compression.

# Building ClinVar Database

This assumes that you have converted a recent ClinVar XML file to TSV using [clinvar-tsv](https://github.com/bihealth/clinvar-tsv).

```
$ mehari db create seqvar-clinvar \
    --path-output-db ~/Data/mehari/db/grch37/seqsvars/clinvar \
    --path-clinvar-tsv path/to/clinvar_seqvars.b37.tsv.gz
```

You can specify an optional `--genome-release grch37` argument that will be used to check the ClinVar database to be compatible with your data.

# Getting HGNC Cross-Link TSV File

You need to provide a TSV file that maps HGNC IDs to gene symbols, RefSeq, and Ensembl gene IDs.
You will need the [jq](https://stedolan.github.io/jq/) command line tool to generate this file.
You can use the provided `misc/genes-xlink-hgnc.jq` script to generate the file.

```
$ wget --no-check-certificate \
    -O /tmp/hgnc_complete_set.json \
    http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/json/hgnc_complete_set.json
$ jq \
    --raw-output \
    --from-file misc/genes-xlink-hgnc.jq \
  > ~/Data/mehari/db/hgnc.tsv
```
