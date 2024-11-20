Getting Started.

# Installation

## via bioconda
As a prerequisite, [follow the bioconda getting started guide](http://bioconda.github.io/#usage).

Then, create a new environment;

```sh
conda create -n mehari -y mehari
conda activate mehari
```

The `mehari` executable is now available from within the activated `mehari` conda environment:

```sh
mehari --help
```

## via docker
Docker images of mehari are available from ghcr.io, see [ghcr.io/varfish-org/mehari](https://github.com/varfish-org/mehari/pkgs/container/mehari).


# Downloading Prebuilt Databases

- transcript database releases: https://github.com/varfish-org/mehari-data-tx/releases
- gnomAD frequency database releases: TODO
- clinVar database releases: https://github.com/varfish-org/annonars-data-clinvar/releases

# Annotating Example VCF Files

You can obtain an example file like this:

```sh
wget https://raw.githubusercontent.com/varfish-org/mehari/main/tests/data/db/create/seqvar_freqs/db-rs1263393206/input.vcf -O example.vcf
```

Now, annotate it using Mehari:

```sh
mehari annotate seqvars \
    --transcripts path/to/mehari-transcript-db \
    --frequencies path/to/mehari-frequency-db \
    --clinvar path/to/mehari-clinvar-db \
    --path-input-vcf example.vcf \
    --path-output-vcf example.out.vcf
```
