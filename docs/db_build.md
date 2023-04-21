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

## Post-Processing Database Files

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

## Building Frequency Databases

You can use the following command for importing the frequencies from gnomAD and HelixMtDb into a RocksDB file to be used for `mehari annotate seqvars` (GRCh37):

```text
$ mehari db create seqvar-freqs \
    --path-output-db output/b37/seqvar/freqs \
    \
    $(for i in {1..22}; do \
        echo --path-gnomad-exomes-auto gnomad.exomes.r2.1.1.sites.chr${i}.vcf.bgz; \
        echo --path-gnomad-genomes-auto gnomad.genomes.r2.1.1.sites.chr${i}.vcf.bgz; \
    done) \
    \
    --path-gnomad-exomes-xy $BASEDIR/input/b37/gnomad.exomes.r2.1.1.sites.chrX.stripped.vcf.bgz \
    --path-gnomad-exomes-xy $BASEDIR/input/b37/gnomad.exomes.r2.1.1.sites.chrY.stripped.vcf.bgz \
    \
    --path-gnomad-genomes-xy $BASEDIR/input/b37/gnomad.genomes.r2.1.1.sites.chrX.stripped.vcf.bgz \
    \
    --path-gnomad-mtdna $BASEDIR/input/gnomad.genomes.v3.1.sites.chrM.vcf.bgz \
    --path-helix-mtdb $BASEDIR/input/HelixMTdb_20200327.vcf.gz
```

Note that you have to build the frequency for the genome build(s) that you want to use (b37 or b38).
Mehari will derive the genome build automatically from the contig lengths.
Note well that for GRCh37, there is no chrY data for gnomAD genomes (r2.1.1).
For GRCh38/r3.1.1, you will also get counts for chrY.

We recommend that you import the r2.1.1 lift-over of gnomAD exomes for GRCh38.

The corresponding command for GRCh38 is:

```text
$ mehari db create seqvar-freqs \
    --path-output-db output/b38/seqvar/freqs \
    \
    $(for i in {1..22}; do \
        echo --path-gnomad-exomes-auto gnomad.exomes.r2.1.1.sites.chr${i}.vcf.bgz; \
        echo --path-gnomad-genomes-auto gnomad.genomes.r3.1.1.sites.chr${i}.vcf.bgz; \
    done) \
    \
    --path-gnomad-exomes-xy $BASEDIR/input/b37/gnomad.exomes.r2.1.1.sites.chrX.stripped.vcf.bgz \
    --path-gnomad-exomes-xy $BASEDIR/input/b37/gnomad.exomes.r2.1.1.sites.chrY.stripped.vcf.bgz \
    \
    --path-gnomad-genomes-xy $BASEDIR/input/b37/gnomad.genomes.r3.1.1.sites.chrX.stripped.vcf.bgz \
    --path-gnomad-genomes-xy $BASEDIR/input/b37/gnomad.genomes.r3.1.1.sites.chrY.stripped.vcf.bgz \
    \
    --path-gnomad-mtdna $BASEDIR/input/gnomad.genomes.v3.1.sites.chrM.vcf.bgz \
    --path-helix-mtdb $BASEDIR/input/HelixMTdb_20200327.vcf.gz
```

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
    --path-out output/b37/txs.bin \
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

You can enable compression by using the suffix `.gz` for gzip compression and `.zstd` for zstandard compression.

# Building ClinVar Database

This assumes that you have converted a recent ClinVar XML file to TSV using [clinvar-tsv](https://github.com/bihealth/clinvar-tsv).

```
$ mehari db create seqvar-clinvar \
    --path-output-db ~/Data/mehari/db/seqvars/grch37/clinvar \
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
