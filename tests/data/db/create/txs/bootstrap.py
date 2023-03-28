#!/usr/bin/env python

import gzip
import json

select_symbols = ("BRCA1", "OPA1")

with gzip.open("../../../../../../cdot-0.2.12.refseq.grch37_grch38.json.gz", "rt") as inputf:
    data = json.load(inputf)

data["transcripts"] = {
    k: v
    for k, v in data["transcripts"].items()
    if v["gene_name"] in select_symbols
}

data["genes"] = {
    k: v
    for k, v in data["genes"].items()
    if v["gene_symbol"] in select_symbols
}

with open("cdot-0.2.12.refseq.grch37_grch38.brca1_opa1.json", "wt") as outputf:
    json.dump(data, outputf, indent=2)
