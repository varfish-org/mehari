Annotating structural variants.

TODO: implement functionality and documentation.

# Information Extracted from Caller VCF

This section describes the information that is extracted from the VCF files of the SV callers.

Independent of the used caller, the following information is extracted:

- Sample names are extracted from the VCF header.
- Chromosome and start position are extracted from the VCF `CHROM` and `POS` columns.
- End position for inversion and deletions are extracted from the `INFO/END` field.
- Only one alternate allele is supported in the `ALT` column.
- The SV type is extracted from the `ALT` allele and the `INFO/SVTYPE` field.
- In the case of breakends, the alternate allele is parsed, e.g., `[17:198982[A` to derive the target chromosome and position as well as the orientation (e.g., 3' to 3' in this case).
- Unless specified otherwise, confidence intervals around start and end positions are extracted from `INFO/CIPOS` and `INFO/CIEND`.
- The genotype for each sample is extracted from the standard field `FORMAT/GT`.
- The genotype quality is extracted from the standard field `FORMAT/GQ`.
- Any genotype-specific filter values are extracted from the standard field `FORMAT/FT`.

The following section describes the SV caller specific information that is extracted.

## Delly

The Delly caller and version is identified by the `INFO/SVMETHOD` field of the first record.

- Paired-end variant read count is extracted from `FORMAT/DV`.
- Paired-end reference read count is extracted from `FORMAT/DR`.
- Split-read variant read count is extracted from `FORMAT/RV`.
- Split-read reference read count is extracted from `FORMAT/RR`.

## Manta

The Manta caller and version is identified by the `##source=` VCF header.
The string `GenerateSVCandidates <VERSION>` is used to identify Manta and the used version.

- Paired-end reference and variant read count is extracted from `FORMAT/PR[0]` and `FORMAT/PR[1]`.
- Split-end reference and variant read count is extracted from `FORMAT/SR[0]` and `FORMAT/SR[1]`.

## GATK gCNV

GATK gCNV is identified by looking at the `##source=` header line.
The tool is identified by `PostprocessGermlineCNVCalls` but the version cannot be identified automatically.

- Copy number is extracted from `FORMAT/CN`
- Number of points is extracted from `FORMAT/NP`

## PopDel

The PopDel caller and version is identified by the `INFO/SVMETHOD` of the first record.

- Paired-end reference support is extracted from `FORMAT/DAD[0]`.
- Paired-end variant support is extracted from `FORMAT/DAD[3]`.

## Dragen SV

The Dragen SV caller is identical to Manta.
However, here the caller and version are identified by looking at the `##source=DRAGEN <VERSION>` line.

## Dragen CNV

The Dragen CNV caller is identified by the header line `##DRAGENVersion=<ID=dragen,Version="SW: <VERSION>, HW: <IGNORED>">`.

- Copy number is extracted from `FORMAT/CN`.
- Number of points is extracted from `FORMAT/PE`.
