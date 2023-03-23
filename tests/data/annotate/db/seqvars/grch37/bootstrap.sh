# Extract transcripts for BRCA1 from the GRCh37 reference genome

cargo \
    run \
    --release \
    -- \
    db \
    create \
    txs \
    --path-out \
    tests/data/annotate/db/seqvars/grch37/txs.bin \
    --path-cdot-json \
    ../cdot-0.2.12.refseq.grch37_grch38.json \
    --path-seqrepo-instance \
    ../hgvs-rs-data/seqrepo-data/master/master \
    --genome-release \
    grch37 \
    --gene-symbols BRCA1
