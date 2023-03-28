Building databases.

This chapter describes how to build the databases needed to run Mehari.

Note that we aim to provide prebuilt databases with Zenodo later on so you don't have to.

# Frequency Database Files

- HelixMtDb can be downloaded from [helix.com](https://www.helix.com/pages/mitochondrial-variant-database)
- gnomAD has a [downloage page](https://gnomad.broadinstitute.org/downloads) where you can find download URLs, e.g., from Microsoft Azure

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

For example, you can reduce the gnomAD genomes r2.1.1 file for chr1 from 31GB with full annotations to 384MB which greatly speeds up the import into a Mehari database.

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
$ wget https://ftp.ensembl.org/pub/grch37/release-108/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz
$ wget https://hgvs-rs-data/mirror/ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.files.installed \
    https://hgvs-rs-data/mirror/ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.{1..12}.rna.fna.gz
```

Now, pull the current seqrepo from the biocommons server:

```text
$ mkdir -p seqrepo-data
$ export SEQREPO_ROOT_DIR=$PWD/seqrepo-data
$ seqrepo pull -i 2021-01-29
```

And then load the FASTA files we downloaded above:

```
$ seqrepo --instance-name master --namespace ENSEMBL \
    Homo_sapiens.GRCh37.cdna.all.fa.gz
$ seqrepo --instance-name master --namespace RefSeq \
            human.{1..12}.rna.fna.gz
```

Now you can create the transcript database (see below).
For some transcripts, no sequence will be available and they will be skipped.
This will be OK as there will be a more recent version available.

## Building Transcript Database

You can build the transcript database flatbuffers binary using the following command:

```text
$ mehari db create txs \
    --path-out output/b37/txs.bin \
    \
    --cdot-json cdot-0.2.12.refseq.grch37_grch38.json \
    --cdot-json cdot-0.2.12.ensembl.grch37_grch38.json \
    \
    --path-seqrepo-instance path/to/seqrepo-data/master \
    \
    --genome-release grch37
```

You will have to build the transcript database for each genome release that you want and manually specify the release to `--genome-release`.
