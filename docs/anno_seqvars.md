Annotating sequence variants (SNV/indel/MNV).

## Sequence Variants And Their VCF Representation

Sequence variants are small variants that are typically represented by the substitution of one or more bases in the reference genome.
In a VCF file, this will typically look as follows:

```text
#CHROM  POS       ID  REF  ALT  QUAL  FILTER  INFO
17      41197700  .   A    C    .     .       .
17      41197708  .   TG   T    .     .       .
```

In the literature, there generally is made a distinction between sequence (or: small) variants and structural variants.
There, small variants are those affecting up to 50bp while structural variants are those that are larger.
Mehari does not enforce such a restriction but attempt to annotate all variants based on their VCF representation (chromosome, 1-based position, reference base and alternative bases).
Symbolic alleles are annotated as the change they stand for, from the base behind `POS` up to `SVLEN` or `END` (the `SVLEN` of the allele counts first, as in VCF 4.5, unless it is 0, the length difference of an `<INV>` in older files):

- `<DEL>` like the deletion of these bases, `<INV>` like their inversion. With a reference, a precise inversion leaves out the outer bases that equal the complement of the base at the other end, as the explicit allele does.
- `<DUP>` like the tandem duplication of these bases if the copy follows them: with an adjacency claim (`SVCLAIM` `J` or `DJ`), or as `<DUP:TANDEM>` without `SVCLAIM`. Otherwise the copy lies elsewhere, so the transcript stays as it is and its copy number rises (`copy_number_increase`, or `transcript_amplification` for a dup of all exons); HGVS describes the copy as an insertion at an unknown place.
- `<INS>` like the insertion of `SVLEN` unknown bases, `N[len]`, with the frame terms of that length, or of an unknown number (`N[?]`, also for `SVLEN=0`), which elongates the transcript from within (`internal_feature_elongation`) and changes the CDS in an unknown way.
- `<CNV>` as a copy number change without direction (`copy_number_change`, no HGVS). With `CN`, it is a copy number call of a deletion or duplication, as is `<DEL>` or `<DUP>` with `SVCLAIM=D`. From VCF 4.4 on, `CN` is the copy number of the allele, rounded: 0 is a deletion, 1 no change, 2 or more a duplication. In older files, it is the total copy number: 0 and 1 are a deletion, 2 no change, 3 or more a duplication. A `<CNV>` without change has no annotations, as VCF 4.5 treats it like `GT=0`. `<CN0>` is a deletion, `<CN1>` no change and `<CN2>` or more a duplication, as older files write them. `<CNV:TR>`, the copy number of a tandem repeat, is not supported.
- Uncertain breakpoints keep the terms of the positions as called, which are not shifted. They have no `c.`, `p.?` and no positions in the transcript, and write `g.` with the intervals, e.g. `g.(1111_1131)_(1240_1280)del`, with `?` for an unknown one, and with the position alone for an interval of one base. Overlapping intervals are cut so that the positions ascend; if the start interval lies behind the end interval, there is no `g.`. An `<INS>` without `CIPOS` has no `g.`. Their annotations carry the message `INFO_IMPRECISE_BREAKPOINTS`. As in VCF 4.5, breakpoints are uncertain with `IMPRECISE`, unless `CIPOS` and `CIEND` (which defaults to `CIPOS`) are both 0; without `CIPOS`, the interval around `POS` is unknown. Without `IMPRECISE`, `CIPOS` and `CIEND` describe the micro-homology of a precise call, which is annotated as such. Breakpoints are also uncertain with `SVCLAIM=D`, and for every `<CNV>` and `<CNn>`.
- A change that reaches beyond the contig ends at its last base, with a warning.
- The change of a symbolic allele has no bases, so Mehari takes those of the transcript for `p.`, which is the same as for the explicit allele.
- A junction that joins two genes gives both their rows `gene_fusion` with its direction: `unidirectional_gene_fusion` where the 5' part of one gene reads into the 3' part of the other on the joined chromosome, `bidirectional_gene_fusion` where two 5' or two 3' parts meet (head to head or tail to tail), and `transcript_fusion` for the junction of a deletion. The junctions are the flanks of a deletion, the end and start of a tandem duplication's bases, both ends of an inversion, and the adjacency of a breakend; a copy number call (`SVCLAIM=D`, `<CNV>`, `<CN0>`) claims none. Genes that both span the junction in the reference, e.g. overlapping genes and readthrough genes, are joined already, so they are no fusion. The record's INFO field `FUSION`, the key `fusions` in JSONL output, lists the gene pairs as `GENE1::GENE2` (symbols, or IDs without a symbol): the 5' gene first, or, for a bidirectional fusion, which has no 5' gene, the gene before the junction on the joined chromosome. `--pick` and plugins see the rows with the fusion terms; a picked transcript can still be one of the gene that the junction does not fuse, while `FUSION` names the gene. A unidirectional fusion also gets its reading frame where both junction bases lie in the CDS of their transcript, from the start codon to the stop codon, in an exon or an intron: `inframe_unidirectional_gene_fusion` if the CDS bases that the 5' transcript keeps, plus the bases inserted at the junction, less those that the 3' transcript skips in front of its first joined exon base, are a multiple of three, else `frameshift_unidirectional_gene_fusion`. With a junction in an intron, the frame assumes normal splicing; the message `INFO_ASSUMES_NORMAL_SPLICING` says so. A junction outside the CDS (e.g. in a UTR), a 5' transcript that keeps less than its start codon or its whole stop codon, a 3' transcript joined from its stop codon on, a non-coding transcript and an imprecise call give no frame term. A transcript that fuses with several transcripts of the other gene can get both frame terms, one per partner.
- A deletion, inversion or tandem duplication of whole exons with both ends in introns (the essential splice sites intact) cannot be translated by HGVS (`p.?`). Its protein terms come from the protein under normal splicing of the remaining exons, with the duplicated exons repeated; the message `INFO_ASSUMES_NORMAL_SPLICING` says so. They are those of the same change with bases: all bases of the exons, also those of a UTR, decide the frame (`frameshift_variant` with `frameshift_truncation` or `frameshift_elongation`, or a conservative or disruptive in-frame term), a deletion through the stop codon gets the stop terms of the protein that reads the 3' UTR, and the loss of all CDS exons gets `start_lost` and `p.0?`. This also applies to such a deletion given with its bases, except where HGVS translates it (`p.0?` for the whole CDS with UTR bases on both sides), which keeps its terms without the message.

A symbolic allele other than these, e.g. `<MYDEL>`, counts as its `SVTYPE` if that is one of the above.
An `SVCLAIM` other than `D`, `J` and `DJ` counts as none, with a warning.

A breakend, e.g. `A]13:123456]`, joins the base at `POS` to a piece of the mate contig, as in VCF 4.5. Each record is annotated on its own:

- A transcript with `POS` inside loses the side that the breakend does not keep. It gets the terms of the deletion of that side over its end, with `feature_truncation`. Transcripts nearby stay whole, up- or downstream of `POS`.
- `HGVS.g` describes the derivative chromosome as in the translocation examples of HGVS: the lost side becomes the inserted bases and the joined piece, e.g. `g.1001_qterdelins[NC_000013.11:g.pter_123456inv]` for `A]13:123456]` at position 1000. There is no `c.`, `n.` or `p.`.
- A single breakend (`G.`, `.A`) gets the same terms without HGVS, and a breakend at a telomere (`POS` 0) gets none.
- A breakend with `IMPRECISE` and a `CIPOS` other than `0,0` or none is imprecise. It keeps the terms at `POS` and carries the message `INFO_IMPRECISE_BREAKPOINTS`. Its `HGVS.g` has the first or last lost base within `CIPOS`, e.g. `g.(501_1501)_qterdelins[...]` for `CIPOS=-500,500`, and is empty without `CIPOS`. Without `IMPRECISE`, `CIPOS` marks the micro-homology of a precise call, as in VCF 4.5, so the breakend is precise.
- The transcripts at the mate get their annotations from the record of the mate. A caller that writes one record per adjacency with `CHR2` and `END` instead leaves them without annotations.

The breakend records of one event are also annotated together, as the symbolic allele of the event:

- Records with one `EVENT` form a group, and so do records that `MATEID` links. A single breakend is in no group.
- Two records of a group pair as mates of each other: by `MATEID` both ways if the mate position of each lies within the `CIPOS` of the other, as VCF 4.5 writes the mate position of an imprecise breakend where the mate would lie if the own position were exact, else by their mate positions alone.
- With `POS` a < b on one contig, the pair `t[b[` at a and `]a]t` at b is the deletion of a+1 to b-1 (`<DEL>`), the pair `]b]t` at a and `t[a[` at b the tandem duplication of a to b (`<DUP:TANDEM>`), and the pairs `t]b]` at a with `t]a]` at b and `[b+1[t` at a+1 with `[a+1[t` at b+1 the inversion of a+1 to b (`<INV>`), whose two pairs must share a group, e.g. by `EVENT`. Bases inserted at a junction do not change the event, but count in the reading frame of a gene fusion that it makes.
- The records of an event get its annotations after their own, with their own `ALT` as the allele, and its gene fusions in `FUSION`. So a transcript between the two junctions of a deletion gets `transcript_ablation`, which neither breakend reports. A group with several pairs, e.g. a complex `EVENT`, gives each pair its own event. With `--pick` and without `--flag-pick`, a transcript keeps one row per record: the row of higher impact, that of the breakend on a tie. In compound mode, the records of an event reach its transcripts, as the symbolic allele does.
- A group waits for mates until the input passes the last mate position on the contig that its records name, plus their largest `CIPOS` maximum. Until then, its records and all records between them are held back, so that the output keeps the order of the input. Then, or at the next contig or the end of the input, the group closes. A record whose mate has not arrived by then keeps its own annotations; Mehari warns about the number of such records. A later record that names the `EVENT` or an ID of a closed group starts a new group.
- A record with a mate on another contig, e.g. of a translocation, gets no event and does not wait for the mate: it keeps its own annotations, which already describe both sides. So do pairs whose junctions describe none of these events, e.g. one junction of an inversion alone.

Records with another symbolic allele, with a symbolic allele without a length, with a malformed breakend, or with another allele at position 0 have no such representation: they stay in the output without annotations, and Mehari warns about their number.
Records without an alternate allele (`ALT` `.`) stay in the output, too.
Mehari stops at a record with more than one alternate allele. Split such records beforehand with `bcftools norm -m-any`, without `-f`, which rejects `N` as `REF`.

The REST endpoint `/api/v1/seqvars/csq` of `mehari server` and the Python bindings annotate symbolic alleles and breakends the same way, and list the gene fusions in `fusions`.
They take one record at a time, so a breakend is annotated on its own, without its mates: the groups of breakends and their events are a feature of `annotate seqvars`.
A symbolic allele takes the values of its INFO fields `END`, `SVLEN`, `SVTYPE`, `SVCLAIM`, `CN`, `CIPOS`, `CIEND` and `IMPRECISE` from the parameters `end`, `sv_length`, `sv_type`, `sv_claim`, `cn`, `cipos`, `ciend` and `imprecise`, e.g. `alternative=<DEL>&end=41197760&cipos=-500,500`.
`sv_type` counts for a custom allele such as `<FOO>` only, as `SVTYPE` does in `annotate seqvars`.
`cn` is always the copy number of the allele as in VCF 4.4, rounded, since no file header gives the VCF version; for an older file, convert its total copy number first.
`cipos` and `ciend` are pairs `min,max` in the query string, and pairs `(min, max)` in Python.
A breakend needs only its `alternative`, and takes `imprecise` and `cipos` for the uncertainty of its position.
A `<CNV>` with `cn` 1 is no variant and gets no annotations.
`*`, `<*>`, `<NON_REF>` and `.` are no variants of their own and give no annotations.
A symbolic allele without an end past `POS`, an unsupported symbolic allele, a malformed breakend, an `sv_claim` other than `D`, `J` and `DJ`, an `sv_type` other than `DEL`, `INS`, `DUP`, `INV`, `CNV` and `BND` for a custom allele, a negative `cn`, and a `cipos` or `ciend` that is not two integers or has its minimum above its maximum give status 400 with the reason, and a `ValueError` in Python.
In a Python DataFrame, such a row gets no annotations and a warning instead, as `annotate seqvars` carries on past such records.
`annotate seqvars` reads an `SVCLAIM` other than `D`, `J` and `DJ` as none, with a warning.
You can find the VCF variant specification [here on Github](https://samtools.github.io/hts-specs/).

## Sequence Variant Annotation

Currently, Mehari will annotate variants using:

- The predicted impact on gene transcripts and the corresponding protein sequence (in the case of coding genes).
- Their frequency in the gnomAD exomes and genomes databases as well as the HelixMtDb database in the case of mitochondrial databases.
- Variant information from ClinVar, if any

## Command Line Invocation

You can invoke Mehari to annotate a VCF file `IN.vcf` creating an output file `OUT.vcf` using the built (or downloaded) databases – for example the transcript database – as follows:

```text
$ mehari annotate seqvars \
    --reference path/to/reference.fa \
    --transcripts path/to/transcripts-db \
    --input IN.vcf \
    --output OUT.vcf
```

Note that the input and output files can optionally be gzip/bgzip compressed VCF files with suffixes (`.gz` or `.bgz`) or BCF files with suffix `.bcf`.
The database assembly should match the one in the input VCF file (e.g., both should either be GRCh37/hg19 or GRCh38/hg38).
You can optionally define the assembly via `--assembly`.
Output formats can be toggled using `--output-format` (e.g., `vcf` or `jsonl`), and you can control the number of threads used for parallel processing by using the `--threads` flag (we recommend a value of around 16).

### Plugins

Plugins add columns to the `ANN` field; enable each with `--plugin`, in the order the columns should appear.
Built-in plugins are `tx-sequence-ref`, `tx-sequence-alt`, `aa-sequence-ref` and `aa-sequence-alt` (transcript and protein sequences without and with the variant), and `nmd` (whether stop-gained and frameshift variants trigger nonsense-mediated decay).
`--plugin` also takes the path of a WebAssembly component, for which you can write your own plugins in Python, Rust and other languages; see [Plugins](`self::plugins`).

### Transcript Selection

By default, Mehari reports one annotation per transcript near the variant, including transcripts up to 5 kbp upstream or downstream.
A change that contains transcripts whole, from their first to their last exon, gives one annotation per gene for them, as they all get the same terms: a deletion, duplication, inversion or copy number change, explicit or symbolic, also the event of grouped breakends, but no insertion and no breakend on its own.
It is the annotation of the transcript that `--pick-order` ranks first, marked with the message `INFO_ONE_ROW_PER_GENE`.
The ranking leaves out `severity`, which cannot tell these transcripts apart, so only the kept transcript is annotated.
Transcripts that the change does not contain whole keep their annotations, also those of the same gene, and VEP compatibility mode reports every transcript.
Plugins and the options below see these annotations.
Four options narrow this down:

- `--restrict-transcripts mane-select,mane-plus-clinical` reports only transcripts with one of these tags.
  A gene without such a transcript gets no annotation, and Mehari logs at the end how many genes this affected.
  With `--restrict-fallback`, such a gene keeps its top-ranked transcript instead, marked with the message `INFO_RESTRICT_FALLBACK`.
- `--pick gene`, `--pick transcript` or `--pick allele` reports one annotation per gene, per transcript or per allele.
  `transcript` only matters when a transcript appears more than once near the variant.
- `--pick-order` ranks the candidates for `--pick`, `--restrict-fallback` and the annotation of a gene that a change contains.
  Each criterion only breaks the ties left by the ones before it, and the lowest accession breaks any remaining tie.
  The default is `mane-select,mane-plus-clinical,refseq-select,ensembl-canonical,gencode-primary,basic,coding,unflagged,severity,length`.
  A tag also matches its backport to GRCh37, while `mane-select-backport` and the other `-backport` values match only the backport.
  Criteria before `severity` need no annotation, so `--pick gene --pick-order mane-select` annotates only the picked transcripts.
- `--flag-pick` keeps all annotations and marks the picked ones with `1` in an extra ANN column, `pick`, after the columns of the plugins.

In these options, a gene is a gene entry of the transcript database: the transcripts of one gene ID without version, e.g. `ENSG00000182378` for all versions of PLCXD1.
So `--pick gene` reports a RefSeq gene and an Ensembl gene separately, also if they have the same HGNC ID.
This also applies where `--pick gene` replaces `--report-most-severe-consequence-by gene` in the table below.

These options replace `--pick-transcript`, `--pick-transcript-mode` and `--report-most-severe-consequence-by`:

| Before | Now |
|---|---|
| `--pick-transcript mane-select --pick-transcript-mode all` | `--restrict-transcripts mane-select --restrict-fallback` |
| `--pick-transcript mane-select --pick-transcript length --pick-transcript-mode first` | `--pick gene --pick-order mane-select,length` |
| `--report-most-severe-consequence-by gene` | `--pick gene --pick-order severity,mane-select,mane-plus-clinical` |
| `--report-most-severe-consequence-by transcript` | `--pick transcript --pick-order severity` |
| `--report-most-severe-consequence-by allele` | `--pick allele --pick-order severity,mane-select,mane-plus-clinical` |

### Experimental: Compound Variants

Mehari includes an experimental feature that allows you to evaluate the compound effect of multiple variants on the same transcript.
You can enable this behavior by passing the `--enable-compound-variants` flag.
The evaluation strategy can be managed via `--phasing-strategy` (options include `strict` [default], `relaxed`, or `ignore`).
`--compound-window-size` sets the maximum distance in bp between neighbouring variants of a group (unlimited if not set; a good value for peptide analyses is 33).

A variant without bases, e.g., a structural variant or a breakend, is not evaluated together with others. If the phasing strategy places it on a haplotype, all its own annotations get the message `WARNING_COMPOUND_INCOMPLETE`, whether another variant shares that haplotype or not. If it shares a haplotype with a variant on a transcript, all annotations of that variant on that transcript get the message too. So do the annotations of variants in a group that fails validation.

## Interpreting Annotation Output

The variant effect/consequence will be formatted similar to the one in, following the `ANN` field standard [documented here](https://pcingola.github.io/SnpEff/se_inputoutput/#ann-field-vcf-output-files).

The `##INFO=<ID=ANN,...>` header line carries the VCF specification's `Source="mehari"` and `Version` (the mehari version) attributes, so that tools can tell it from SnpEff's `ANN` field.
An `ANN` field in the input is therefore removed when annotating with transcripts, also for records that get no consequences, and mehari warns about this once per run.

`--ann-key` writes the annotations to another INFO key, e.g. to keep the `ANN` field of another tool next to mehari's.
Mehari then only removes and warns about an input field under that key.
The key must match the pattern of the VCF specification, `^([A-Za-z_][0-9A-Za-z_.]*|1000G)$`, and must differ from the other INFO fields that the run writes.
`mehari postprocess proteome` and `mehari postprocess varfish-seqvars` read `ANN` by default; pass them the same `--ann-key`.

`Gene_ID` is the gene ID of the transcript's gene record in the annotation, with its version if it has one, e.g. `672` for BRCA1 in RefSeq or `ENSG00000182378.15` for PLCXD1 in GENCODE.
`HGNC_ID`, which mehari adds after `Gene_ID`, is the HGNC ID of the gene, e.g. `HGNC:1100` for BRCA1, and empty for a gene without one.
Mehari refuses a transcript database that predates annotation gene IDs, i.e., one in which no transcript has a gene ID.
Rebuild such a database with the current mehari version.

The population frequency will be written to the VCF INFO fields as follows:

- gnomAD exomes
    - `gnomad_exomes_an` -- number of observed alleles in gnomAD exomes
    - `gnomad_exomes_hom` -- number of homozygous carriers in gnomAD exomes
    - `gnomad_exomes_het` -- number of heterozygous carriers in gnomAD exomes
    - `gnomad_exomes_hemi` -- number of hemizygous carriers in gnomAD exomes
- gnomAD genomes
    - `gnomad_genomes_an` -- number of observed alleles in gnomAD genomes
    - `gnomad_genomes_hom` -- number of homozygous carriers in gnomAD genomes
    - `gnomad_genomes_het` -- number of heterozygous carriers in gnomAD genomes
    - `gnomad_genomes_hemi` -- number of hemizygous carriers in gnomAD genomes
- gnomAD mtDNA
    - `gnomad_mtdna_an` -- number of individuals with sufficient coverage in gnomAD mtDNA database
    - `gnomad_mtdna_hom` -- number of homoplasmic carriers in gnomAD mtDNA database
    - `gnomad_mtdna_het` -- number of heteroplasmic carriers in gnomAD mtDNA database
- HelixMtDb
    - `helix_an` -- number of individuals in HelixMtDb
    - `helix_hom` -- number of homoplasmic carriers in HelixMtDb
    - `helix_het` -- number of heteroplasmic carriers in HelixMtDb

From these integer numbers, the relative allele frequencies can be derived as needed.
