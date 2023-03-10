#!/usr/bin/env bash

base=/data/sshfs/data/gpfs-1/groups/cubi/work/projects/2021-07-20_varfish-db-downloader-holtgrewe/varfish-db-downloader/

tabix --print-header \
    $base/GRCh37/gnomAD_genomes/r2.1.1/download/gnomad.genomes.r2.1.1.sites.chr1.vcf.bgz 1:13656-13656 \
>> tests/data/db/create/seqvar_freqs/12-37/gnomad.genomes.r2.1.1.sites.chr1-rs1263393206.vcf

tabix --print-header \
    $base/GRCh37/gnomAD_exomes/r2.1.1/download/gnomad.exomes.r2.1.1.sites.chr1.vcf.bgz 1:13656-13656 \
>> tests/data/db/create/seqvar_freqs/12-37/gnomad.exomes.r2.1.1.sites.chr1-rs1263393206.vcf

mkdir -p tests/data/db/create/seqvar_freqs/db-rs1263393206/seqvars/freqs

cargo run -- \
    db create seqvar-freqs \
        -v \
        --path-output-db tests/data/db/create/seqvar_freqs/db-rs1263393206/seqvars/freqs \
        --genome-release grch37 \
        --path-gnomad-exomes-auto tests/data/db/create/seqvar_freqs/12-37/gnomad.exomes.r2.1.1.sites.chr1-rs1263393206.vcf \
        --path-gnomad-genomes-auto tests/data/db/create/seqvar_freqs/12-37/gnomad.genomes.r2.1.1.sites.chr1-rs1263393206.vcf \
