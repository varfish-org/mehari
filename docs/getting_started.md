Getting Started.

# Installation

You most likely want to install via bioconda.
As a prerequisite, [follow the bioconda getting started guide](http://bioconda.github.io/#usage).

Then, create a new environment (use the `mamba` if you are as impatient as us).

```text
$ mamba create -y mehari mehari
$ conda activate mehari
```

The `mehari` executable is now available:

```
$ mehari --help
```

# Downloading Prebuilt Databases

TODO: not yet available

# Annotating Example VCF Files

You can obtain an example file like this:

```text
$ wget https://raw.githubusercontent.com/bihealth/mehari/main/tests/data/db/create/seqvar_freqs/db-rs1263393206/input.vcf \
    -O example.vcf
```

Now, annotate it using Mehari:

```text
$ mehari annotate seqvars \
    --path-db path/to/mehari-db/b37 \
    --path-input-vcf example.vcf \
    --path-output-vcf example.out.vcf
$ grep -v ^# example.out.vcf
TODO: output line
```
