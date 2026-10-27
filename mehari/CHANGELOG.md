# Changelog

## [0.48.1](https://github.com/varfish-org/mehari/compare/mehari-v0.48.0...mehari-v0.48.1) (2026-10-27)


### Bug Fixes

* **annotate:** consequence terms per SO and HGVS ([#1147](https://github.com/varfish-org/mehari/issues/1147)) ([c9f4470](https://github.com/varfish-org/mehari/commit/c9f4470c0dd4966a653e8573aaea3dc08bbac49c))

## [0.48.0](https://github.com/varfish-org/mehari/compare/mehari-v0.47.0...mehari-v0.48.0) (2026-10-26)


### ⚠ BREAKING CHANGES

* **deps:** bump hgvs to 4.0 ([#1125](https://github.com/varfish-org/mehari/issues/1125))

### Features

* **annotate:** warn about the chrM of hg19 ([#1141](https://github.com/varfish-org/mehari/issues/1141)) ([386e517](https://github.com/varfish-org/mehari/commit/386e517854d71b782da37c98288404390c28905f))
* **annotate:** warn when REF differs from genome ([#1139](https://github.com/varfish-org/mehari/issues/1139)) ([815e059](https://github.com/varfish-org/mehari/commit/815e059a234d2d12e8bc39383b1b3f6c8a0805d1))
* limit repeated warnings per call site ([#1138](https://github.com/varfish-org/mehari/issues/1138)) ([2bf07d6](https://github.com/varfish-org/mehari/commit/2bf07d64312e06a458e59d96a5c1e80328956ce7))


### Bug Fixes

* **annotate:** conservative in-frame insertions ([#1132](https://github.com/varfish-org/mehari/issues/1132)) ([ecd1753](https://github.com/varfish-org/mehari/commit/ecd1753ab4b13572c6df7c05c86aebd7830994df))
* **annotate:** drop transcripts out of reach ([#1126](https://github.com/varfish-org/mehari/issues/1126)) ([cd058d2](https://github.com/varfish-org/mehari/commit/cd058d2daaaf259be264906d5abf8c65b579fbdb))
* **annotate:** give ambiguous bases VEP's terms ([#1128](https://github.com/varfish-org/mehari/issues/1128)) ([d297790](https://github.com/varfish-org/mehari/commit/d297790524b6963cd4e07cf0a6661dec50e42cee))
* **annotate:** match VEP's polypyrimidine tract ([#1127](https://github.com/varfish-org/mehari/issues/1127)) ([1e5e68c](https://github.com/varfish-org/mehari/commit/1e5e68cf63540c7e17fe7573e146068dd4cc897e))
* **annotate:** read N in ALT as one unknown base ([#1130](https://github.com/varfish-org/mehari/issues/1130)) ([832d5a5](https://github.com/varfish-org/mehari/commit/832d5a5c13b62c88e9d7cdd865a7b44511b7b3b7))
* **annotate:** read Xaa as unknown, not as a stop ([#1131](https://github.com/varfish-org/mehari/issues/1131)) ([fa70f6c](https://github.com/varfish-org/mehari/commit/fa70f6c631b55c858a96a6cedbc9e8ae80d59d57))
* **annotate:** skip filtered transcripts in lookup ([#1136](https://github.com/varfish-org/mehari/issues/1136)) ([5ae4c7b](https://github.com/varfish-org/mehari/commit/5ae4c7bec7fe4f715443dd5e971194ba67c66bd3))
* **annotate:** take VEP's in-frame insertion terms ([#1129](https://github.com/varfish-org/mehari/issues/1129)) ([513fb99](https://github.com/varfish-org/mehari/commit/513fb9936a10e39173ec25882318a6ff9ec93ec5))
* **db:** merge transcripts without duplicates ([#1140](https://github.com/varfish-org/mehari/issues/1140)) ([c4d4797](https://github.com/varfish-org/mehari/commit/c4d47976290820607d38586e601d1fae5d0d61af))
* **db:** shift cdot codons only before data 0.2.34 ([#1135](https://github.com/varfish-org/mehari/issues/1135)) ([16fa64d](https://github.com/varfish-org/mehari/commit/16fa64db8c4a0bc34e4d92f5134143806fa782cc))
* turn logging off with -qqq ([#1137](https://github.com/varfish-org/mehari/issues/1137)) ([478e990](https://github.com/varfish-org/mehari/commit/478e990aef5b47cd509bf0a030c4c5b126fbc011))


### Performance Improvements

* **annotate:** format hgvs_g once per strand ([#1133](https://github.com/varfish-org/mehari/issues/1133)) ([a6e4d5a](https://github.com/varfish-org/mehari/commit/a6e4d5a235aa2989d85cfe419410ff3de5418721))


### Build System

* **deps:** bump hgvs to 4.0 ([#1125](https://github.com/varfish-org/mehari/issues/1125)) ([d83a51d](https://github.com/varfish-org/mehari/commit/d83a51d3fedee62e8f1e627bf8889a485b9eef39))

## [0.47.0](https://github.com/varfish-org/mehari/compare/mehari-v0.46.0...mehari-v0.47.0) (2026-10-22)


### ⚠ BREAKING CHANGES

* **db:** take gene_id from the annotation, add hgnc_id ([#1085](https://github.com/varfish-org/mehari/issues/1085))
* **db:** import GFF3 transcripts like cdot ([#1065](https://github.com/varfish-org/mehari/issues/1065))
* **annotate:** redesign transcript picking options ([#1093](https://github.com/varfish-org/mehari/issues/1093)) (#1096)
* WebAssembly plugins ([#426](https://github.com/varfish-org/mehari/issues/426)) (#1091)
* add transcript plugins and an NMD plugin ([#246](https://github.com/varfish-org/mehari/issues/246)) (#1107)
* normalize MT variants across the circular origin ([#1087](https://github.com/varfish-org/mehari/issues/1087))
* **deps:** bump hgvs to 2.0 and the MSRV to 1.92 ([#1080](https://github.com/varfish-org/mehari/issues/1080))
* restrict mehari db to transcripts ([#1033](https://github.com/varfish-org/mehari/issues/1033)) (#1074)

### Features

* add transcript plugins and an NMD plugin ([#246](https://github.com/varfish-org/mehari/issues/246)) ([#1107](https://github.com/varfish-org/mehari/issues/1107)) ([4c79261](https://github.com/varfish-org/mehari/commit/4c792612f0c66974dc2b61158fcf70888a93af7f))
* **annotate:** add --compound-window-size ([#1109](https://github.com/varfish-org/mehari/issues/1109)) ([c478355](https://github.com/varfish-org/mehari/commit/c47835512d2c46a0093f24e4aa3f8929034b8643))
* **annotate:** redesign transcript picking options ([#1093](https://github.com/varfish-org/mehari/issues/1093)) ([#1096](https://github.com/varfish-org/mehari/issues/1096)) ([5da296b](https://github.com/varfish-org/mehari/commit/5da296b33cbe99fddbc81def5b3a5c837cceb1e4))
* **annotate:** warn when no reference genome is given ([#1112](https://github.com/varfish-org/mehari/issues/1112)) ([b9fb18c](https://github.com/varfish-org/mehari/commit/b9fb18c47e0fe2acc9369d656139ba945db1b80c))
* **annotate:** write predicted protein changes in parentheses ([#1115](https://github.com/varfish-org/mehari/issues/1115)) ([2037881](https://github.com/varfish-org/mehari/commit/2037881525dfb8ceb34dc6c8711eb0409953b1b6))
* apply translation exceptions and complete stop codons only with evidence ([#1068](https://github.com/varfish-org/mehari/issues/1068)) ([b87f3af](https://github.com/varfish-org/mehari/commit/b87f3af6eadbaf5c5b89a758b11ce9ead5ce3c69))
* **db:** import GFF3 transcripts like cdot ([#1065](https://github.com/varfish-org/mehari/issues/1065)) ([c0849bf](https://github.com/varfish-org/mehari/commit/c0849bf6c179ee1032f37611789b478d7a39f5db))
* **db:** take gene_id from the annotation, add hgnc_id ([#1085](https://github.com/varfish-org/mehari/issues/1085)) ([776f2e0](https://github.com/varfish-org/mehari/commit/776f2e035147c68d54e39e2db46af2f13d4692ea))
* **server:** report missing transcripts as message ([#386](https://github.com/varfish-org/mehari/issues/386)) ([#1095](https://github.com/varfish-org/mehari/issues/1095)) ([2cee16b](https://github.com/varfish-org/mehari/commit/2cee16b827baac51cc34c29a321015548ce1e0c1))
* WebAssembly plugins ([#426](https://github.com/varfish-org/mehari/issues/426)) ([#1091](https://github.com/varfish-org/mehari/issues/1091)) ([6bab41a](https://github.com/varfish-org/mehari/commit/6bab41a013ca4d1b1a166a8f9723bf24f471561e))


### Bug Fixes

* **annotate:** apply VEP rules only in compatibility mode ([#1111](https://github.com/varfish-org/mehari/issues/1111)) ([936c7a2](https://github.com/varfish-org/mehari/commit/936c7a2f61f50c27adaaebc9cfe096aee31db9e9))
* **annotate:** correct indel, UTR and start codon annotations ([#1067](https://github.com/varfish-org/mehari/issues/1067)) ([3953d9d](https://github.com/varfish-org/mehari/commit/3953d9dc23d6761a7951aaa40e00d49111dd6448))
* **annotate:** decide insertions at exon edges by sequence ([#1114](https://github.com/varfish-org/mehari/issues/1114)) ([62373fa](https://github.com/varfish-org/mehari/commit/62373fa065bb0b503a90166916e8ea782ad0a62b))
* **annotate:** derive SV type from ALT if no SVTYPE ([#557](https://github.com/varfish-org/mehari/issues/557)) ([#1090](https://github.com/varfish-org/mehari/issues/1090)) ([e0a174b](https://github.com/varfish-org/mehari/commit/e0a174b85d64c20133bd232518288ba20ffad871))
* **annotate:** discard UTR splice terms by region ([#1082](https://github.com/varfish-org/mehari/issues/1082)) ([1829566](https://github.com/varfish-org/mehari/commit/1829566d010514af6633461b51a288935c4e65c1))
* **annotate:** no CDS.pos for UTR variants ([#1099](https://github.com/varfish-org/mehari/issues/1099)) ([#1101](https://github.com/varfish-org/mehari/issues/1101)) ([834f658](https://github.com/varfish-org/mehari/commit/834f6582026c5afcd0cc100840f80d8cb7725e51))
* **annotate:** no frameshift for indels over the CDS end ([#1081](https://github.com/varfish-org/mehari/issues/1081)) ([b1c2741](https://github.com/varfish-org/mehari/commit/b1c2741d3be9bb3ddde038d68ded412457cacfb3))
* **annotate:** normalize phased variants as one change ([#1083](https://github.com/varfish-org/mehari/issues/1083)) ([f046f88](https://github.com/varfish-org/mehari/commit/f046f88f7b57b24254f76ec26bcf903a4ba41e9e))
* **annotate:** number introns in transcript order ([#1098](https://github.com/varfish-org/mehari/issues/1098)) ([#1100](https://github.com/varfish-org/mehari/issues/1100)) ([471fc81](https://github.com/varfish-org/mehari/commit/471fc81ed36d89c38999c2939007eb6c125f6141))
* **annotate:** order compound variant groups deterministically ([#1097](https://github.com/varfish-org/mehari/issues/1097)) ([66fb709](https://github.com/varfish-org/mehari/commit/66fb7099246307f0a17e5657ca34e365aea3e25e))
* **annotate:** reconcile consequence terms with the protein sequence ([#1116](https://github.com/varfish-org/mehari/issues/1116)) ([604d350](https://github.com/varfish-org/mehari/commit/604d350364175367824646aafa94148b20aeca48))
* **annotate:** take indel terms and p. from the 3'-shifted HGVS ([#1069](https://github.com/varfish-org/mehari/issues/1069)) ([034515c](https://github.com/varfish-org/mehari/commit/034515c6418e3cc6ec784ad6f1c076676741a44b))
* **annotate:** take the terms of stop codon deletions from the protein ([#1113](https://github.com/varfish-org/mehari/issues/1113)) ([d0ee375](https://github.com/varfish-org/mehari/commit/d0ee375c9fc32053e77ea3db8385c529aa92aefb))
* **annotate:** transcript length, strand and totals ([#1066](https://github.com/varfish-org/mehari/issues/1066)) ([9f63f8d](https://github.com/varfish-org/mehari/commit/9f63f8dfe407dbbf3950f6b6bc91fe67731198e3))
* **annotate:** VEP stop terms in partial codons ([#1084](https://github.com/varfish-org/mehari/issues/1084)) ([8d4c805](https://github.com/varfish-org/mehari/commit/8d4c805faec5d86cf52359086c5068d553d1b419))
* correct edge cases found in the final review of the stack ([#1117](https://github.com/varfish-org/mehari/issues/1117)) ([b861941](https://github.com/varfish-org/mehari/commit/b861941ba466d8795100904dab035d6bbbf6da0d))
* normalize MT variants across the circular origin ([#1087](https://github.com/varfish-org/mehari/issues/1087)) ([3550666](https://github.com/varfish-org/mehari/commit/355066678a070fb666b1eced572a998912a3d3ae))


### Performance Improvements

* **annotate:** overlap writing and prepare VCF records on workers ([#1102](https://github.com/varfish-org/mehari/issues/1102)) ([1f918b6](https://github.com/varfish-org/mehari/commit/1f918b6acf7ecf03d2d41800bf833f1a0e4ff539))
* **annotate:** take the region from c_to_p and stop copying for plugins ([#1118](https://github.com/varfish-org/mehari/issues/1118)) ([02e55a4](https://github.com/varfish-org/mehari/commit/02e55a4471709815c19a4f1992d93a79a985d94c))
* pre-size transcript DB buffer from zstd header ([#1103](https://github.com/varfish-org/mehari/issues/1103)) ([0e42732](https://github.com/varfish-org/mehari/commit/0e42732ccc677f3ee13e79e14ff755f188accc69))


### Code Refactoring

* restrict mehari db to transcripts ([#1033](https://github.com/varfish-org/mehari/issues/1033)) ([#1074](https://github.com/varfish-org/mehari/issues/1074)) ([ffb923c](https://github.com/varfish-org/mehari/commit/ffb923cca27d831faedad13178c6a6f6c6692af3))


### Build System

* **deps:** bump hgvs to 2.0 and the MSRV to 1.92 ([#1080](https://github.com/varfish-org/mehari/issues/1080)) ([2228a2a](https://github.com/varfish-org/mehari/commit/2228a2a7e888234390aefec5375ec9962fbae513))

## [0.46.0](https://github.com/varfish-org/mehari/compare/mehari-v0.45.1...mehari-v0.46.0) (2026-09-22)


### Features

* `server run` cli allow user to control number of web worker ([#1058](https://github.com/varfish-org/mehari/issues/1058)) ([88a6de3](https://github.com/varfish-org/mehari/commit/88a6de344a99086139dda89731456cb9417d3704))
* **mehari-python:** add tqdm-compatible build progress ([#1055](https://github.com/varfish-org/mehari/issues/1055)) ([d54313d](https://github.com/varfish-org/mehari/commit/d54313d5e2730c8ef22d4184a0d2e4cf016c121f))


### Bug Fixes

* **db:** honor CDS phase when importing GFF3 transcripts ([#1048](https://github.com/varfish-org/mehari/issues/1048)) ([90e3d89](https://github.com/varfish-org/mehari/commit/90e3d89b0496a0d3743c7615934bb24a3e5936eb))
* **db:** merge repeated --mane-transcripts rows per transcript ([#1052](https://github.com/varfish-org/mehari/issues/1052)) ([db94ae5](https://github.com/varfish-org/mehari/commit/db94ae504e009f51399c4972cf471d058a8c0a93))
* **db:** pad incomplete CDS with CIGAR D instead of I ([#1062](https://github.com/varfish-org/mehari/issues/1062)) ([1d80aeb](https://github.com/varfish-org/mehari/commit/1d80aeb166f067c16ac6b1573994c552ef345313))
* read multi-member gzip (BGZF) annotation files completely ([#1050](https://github.com/varfish-org/mehari/issues/1050)) ([ff5766c](https://github.com/varfish-org/mehari/commit/ff5766cf59cdb3b1cb2544855a6640c6091c1b41))
* store GFF3 exons in genomic order for minus-strand transcripts ([#1046](https://github.com/varfish-org/mehari/issues/1046)) ([c61fa7d](https://github.com/varfish-org/mehari/commit/c61fa7dd86e8f367b966b6cfb1adb31ca9fe34db))

## [0.45.1](https://github.com/varfish-org/mehari/compare/mehari-v0.45.0...mehari-v0.45.1) (2026-07-03)


### Bug Fixes

* spliceAI and CADD VCF header types ([#1027](https://github.com/varfish-org/mehari/issues/1027)) ([c9b3132](https://github.com/varfish-org/mehari/commit/c9b313227c155373f1330a489964995e637316ff))

## [0.45.0](https://github.com/varfish-org/mehari/compare/mehari-v0.44.1...mehari-v0.45.0) (2026-07-02)


### ⚠ BREAKING CHANGES

* support CADD, spliceAI, dbSNP and building custom/generic lookup annotation DBs ([#1025](https://github.com/varfish-org/mehari/issues/1025))

### Features

* support CADD, spliceAI, dbSNP and building custom/generic lookup annotation DBs ([#1025](https://github.com/varfish-org/mehari/issues/1025)) ([7abb129](https://github.com/varfish-org/mehari/commit/7abb12903fde91bd4341c10726e2ca46914769a6))

## [0.44.1](https://github.com/varfish-org/mehari/compare/mehari-v0.44.0...mehari-v0.44.1) (2026-06-17)


### Bug Fixes

* exclude filtered out transcripts from gene-to-tx mapping ([#1016](https://github.com/varfish-org/mehari/issues/1016)) ([17d3bf7](https://github.com/varfish-org/mehari/commit/17d3bf77bcb50249d6d7999c3912741755bdb4f7))

## [0.44.0](https://github.com/varfish-org/mehari/compare/mehari-v0.43.3...mehari-v0.44.0) (2026-06-15)


### Features

* add regex db subset options ([#1012](https://github.com/varfish-org/mehari/issues/1012)) ([120f304](https://github.com/varfish-org/mehari/commit/120f304bd667b2a95e91578fb670b5bd04e1ab64))

## [0.43.3](https://github.com/varfish-org/mehari/compare/mehari-v0.43.2...mehari-v0.43.3) (2026-06-10)


### Bug Fixes

* more robust freqs/clinvar assembly checks ([#1008](https://github.com/varfish-org/mehari/issues/1008)) ([d0ead01](https://github.com/varfish-org/mehari/commit/d0ead012210ae2daa8947f97e6130d0b81922961))
* varfish postprocess TSV assembly capitalization ([#1007](https://github.com/varfish-org/mehari/issues/1007)) ([d3c997d](https://github.com/varfish-org/mehari/commit/d3c997d01ff5da20a493ad3c621e777e7806c988))

## [0.43.2](https://github.com/varfish-org/mehari/compare/mehari-v0.43.0...mehari-v0.43.2) (2026-05-07)


### ⚠ BREAKING CHANGES

* add experimental option to predict compound effect of multiple variants on the same transcript, allow building txdb from gff3+fasta ([#972](https://github.com/varfish-org/mehari/issues/972))

### Features

* add experimental option to predict compound effect of multiple variants on the same transcript, allow building txdb from gff3+fasta ([#972](https://github.com/varfish-org/mehari/issues/972)) ([e377c4e](https://github.com/varfish-org/mehari/commit/e377c4eeeb7da1644feba7d37be7a71dc86ada69))
* add options to include reference/alternative cDNA and protein sequences in output ([#966](https://github.com/varfish-org/mehari/issues/966)) ([ec87cb0](https://github.com/varfish-org/mehari/commit/ec87cb0e5c8bf6c70672a5d2b8fd5f76da82c0b5))
* Add python bindings via pyO3 ([#959](https://github.com/varfish-org/mehari/issues/959)) ([b2450be](https://github.com/varfish-org/mehari/commit/b2450be38bdf65aa0c21c5b77ea6770910bd46a5))
* allow db subsetting by region ([#988](https://github.com/varfish-org/mehari/issues/988)) ([09aa020](https://github.com/varfish-org/mehari/commit/09aa0208187cda4928295e7f8e54d2008059de2b))


### Bug Fixes

* hide server functionality behind feature gate (enabled by default) ([#964](https://github.com/varfish-org/mehari/issues/964)) ([0031153](https://github.com/varfish-org/mehari/commit/00311532ac5d296277f85220d415d5a1cc96dbd3))
* pick mode first and length picking ([#993](https://github.com/varfish-org/mehari/issues/993)) ([e63b7b0](https://github.com/varfish-org/mehari/commit/e63b7b0c334b092e3a2f63a37aebe861a949697a))
* use with_added_extension ([#989](https://github.com/varfish-org/mehari/issues/989)) ([0dd76f9](https://github.com/varfish-org/mehari/commit/0dd76f90046badb479b0784c25d79bf891bf6705))


### Miscellaneous Chores

* bump version for annonars update ([1e625f9](https://github.com/varfish-org/mehari/commit/1e625f9a3a90854427669609296b93e8b7c964ed))
* bump version for jemalloc default feature ([517e1d3](https://github.com/varfish-org/mehari/commit/517e1d3236d39d8670e57d6870bb1f48b6dcf4e0))
* Release-As: 0.43.2 ([e170371](https://github.com/varfish-org/mehari/commit/e170371ae6fcd5cc4cab847d20433087fa474513))

## [0.35.0](https://github.com/varfish-org/mehari/compare/mehari-v0.43.0...mehari-v0.35.0) (2026-05-07)


### ⚠ BREAKING CHANGES

* add experimental option to predict compound effect of multiple variants on the same transcript, allow building txdb from gff3+fasta ([#972](https://github.com/varfish-org/mehari/issues/972))

### Features

* add experimental option to predict compound effect of multiple variants on the same transcript, allow building txdb from gff3+fasta ([#972](https://github.com/varfish-org/mehari/issues/972)) ([e377c4e](https://github.com/varfish-org/mehari/commit/e377c4eeeb7da1644feba7d37be7a71dc86ada69))
* add options to include reference/alternative cDNA and protein sequences in output ([#966](https://github.com/varfish-org/mehari/issues/966)) ([ec87cb0](https://github.com/varfish-org/mehari/commit/ec87cb0e5c8bf6c70672a5d2b8fd5f76da82c0b5))
* Add python bindings via pyO3 ([#959](https://github.com/varfish-org/mehari/issues/959)) ([b2450be](https://github.com/varfish-org/mehari/commit/b2450be38bdf65aa0c21c5b77ea6770910bd46a5))
* allow db subsetting by region ([#988](https://github.com/varfish-org/mehari/issues/988)) ([09aa020](https://github.com/varfish-org/mehari/commit/09aa0208187cda4928295e7f8e54d2008059de2b))


### Bug Fixes

* hide server functionality behind feature gate (enabled by default) ([#964](https://github.com/varfish-org/mehari/issues/964)) ([0031153](https://github.com/varfish-org/mehari/commit/00311532ac5d296277f85220d415d5a1cc96dbd3))
* pick mode first and length picking ([#993](https://github.com/varfish-org/mehari/issues/993)) ([e63b7b0](https://github.com/varfish-org/mehari/commit/e63b7b0c334b092e3a2f63a37aebe861a949697a))
* use with_added_extension ([#989](https://github.com/varfish-org/mehari/issues/989)) ([0dd76f9](https://github.com/varfish-org/mehari/commit/0dd76f90046badb479b0784c25d79bf891bf6705))


### Miscellaneous Chores

* bump version for annonars update ([1e625f9](https://github.com/varfish-org/mehari/commit/1e625f9a3a90854427669609296b93e8b7c964ed))
* bump version for jemalloc default feature ([517e1d3](https://github.com/varfish-org/mehari/commit/517e1d3236d39d8670e57d6870bb1f48b6dcf4e0))

## [0.43.0](https://github.com/varfish-org/mehari/compare/v0.42.0...v0.43.0) (2026-05-06)


### Features

* allow db subsetting by region ([#988](https://github.com/varfish-org/mehari/issues/988)) ([09aa020](https://github.com/varfish-org/mehari/commit/09aa0208187cda4928295e7f8e54d2008059de2b))


### Bug Fixes

* pick mode first and length picking ([#993](https://github.com/varfish-org/mehari/issues/993)) ([e63b7b0](https://github.com/varfish-org/mehari/commit/e63b7b0c334b092e3a2f63a37aebe861a949697a))
* use with_added_extension ([#989](https://github.com/varfish-org/mehari/issues/989)) ([0dd76f9](https://github.com/varfish-org/mehari/commit/0dd76f90046badb479b0784c25d79bf891bf6705))

## [0.42.0](https://github.com/varfish-org/mehari/compare/v0.41.2...v0.42.0) (2026-04-14)


### ⚠ BREAKING CHANGES

* add experimental option to predict compound effect of multiple variants on the same transcript, allow building txdb from gff3+fasta ([#972](https://github.com/varfish-org/mehari/issues/972))

### Features

* add experimental option to predict compound effect of multiple variants on the same transcript, allow building txdb from gff3+fasta ([#972](https://github.com/varfish-org/mehari/issues/972)) ([e377c4e](https://github.com/varfish-org/mehari/commit/e377c4eeeb7da1644feba7d37be7a71dc86ada69))
* add options to include reference/alternative cDNA and protein sequences in output ([#966](https://github.com/varfish-org/mehari/issues/966)) ([ec87cb0](https://github.com/varfish-org/mehari/commit/ec87cb0e5c8bf6c70672a5d2b8fd5f76da82c0b5))

## [0.41.2](https://github.com/varfish-org/mehari/compare/v0.41.1...v0.41.2) (2026-03-22)


### Bug Fixes

* hide server functionality behind feature gate (enabled by default) ([#964](https://github.com/varfish-org/mehari/issues/964)) ([0031153](https://github.com/varfish-org/mehari/commit/00311532ac5d296277f85220d415d5a1cc96dbd3))

## [0.41.0](https://github.com/varfish-org/mehari/compare/v0.40.0...v0.41.0) (2026-03-21)


### Features

* Add python bindings via pyO3 ([#959](https://github.com/varfish-org/mehari/issues/959)) ([b2450be](https://github.com/varfish-org/mehari/commit/b2450be38bdf65aa0c21c5b77ea6770910bd46a5))
