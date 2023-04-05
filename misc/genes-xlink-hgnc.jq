# Create header
["hgnc_id", "ensembl_gene_id", "entrez_id", "gene_symbol"],
# Then, process the file
(
    # For each entry in .response.docs
    .response.docs[] |
        # Generate array of ensembl_gene_id, entrez_id, and symbol, using
        # the empty string for the default
        [
            .hgnc_id // "",
            .ensembl_gene_id // "",
            .entrez_id // "",
            .symbol // ""
        ]
)
# Convert everything to TSV
| @tsv
