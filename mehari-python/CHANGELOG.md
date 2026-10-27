# Changelog

## [0.49.0](https://github.com/varfish-org/mehari/compare/mehari-python-v0.48.1...mehari-python-v0.49.0) (2026-10-27)


### Features

* **annotate:** report the HGNC ID in ANN ([#1149](https://github.com/varfish-org/mehari/issues/1149)) ([7d13cde](https://github.com/varfish-org/mehari/commit/7d13cde0f17b151dcf45978ce9ca9ac7aba1ca59))

## [0.48.1](https://github.com/varfish-org/mehari/compare/mehari-python-v0.48.0...mehari-python-v0.48.1) (2026-10-27)


### Miscellaneous Chores

* **mehari-python:** Synchronize mehari-sync versions

## [0.48.0](https://github.com/varfish-org/mehari/compare/mehari-python-v0.47.0...mehari-python-v0.48.0) (2026-10-26)


### Miscellaneous Chores

* **mehari-python:** Synchronize mehari-sync versions

## [0.47.0](https://github.com/varfish-org/mehari/compare/mehari-python-v0.46.0...mehari-python-v0.47.0) (2026-10-22)


### ⚠ BREAKING CHANGES

* **db:** take gene_id from the annotation, add hgnc_id ([#1085](https://github.com/varfish-org/mehari/issues/1085))
* **annotate:** redesign transcript picking options ([#1093](https://github.com/varfish-org/mehari/issues/1093)) (#1096)
* add transcript plugins and an NMD plugin ([#246](https://github.com/varfish-org/mehari/issues/246)) (#1107)
* **deps:** bump hgvs to 2.0 and the MSRV to 1.92 ([#1080](https://github.com/varfish-org/mehari/issues/1080))
* restrict mehari db to transcripts ([#1033](https://github.com/varfish-org/mehari/issues/1033)) (#1074)

### Features

* add transcript plugins and an NMD plugin ([#246](https://github.com/varfish-org/mehari/issues/246)) ([#1107](https://github.com/varfish-org/mehari/issues/1107)) ([4c79261](https://github.com/varfish-org/mehari/commit/4c792612f0c66974dc2b61158fcf70888a93af7f))
* **annotate:** redesign transcript picking options ([#1093](https://github.com/varfish-org/mehari/issues/1093)) ([#1096](https://github.com/varfish-org/mehari/issues/1096)) ([5da296b](https://github.com/varfish-org/mehari/commit/5da296b33cbe99fddbc81def5b3a5c837cceb1e4))
* **annotate:** write predicted protein changes in parentheses ([#1115](https://github.com/varfish-org/mehari/issues/1115)) ([2037881](https://github.com/varfish-org/mehari/commit/2037881525dfb8ceb34dc6c8711eb0409953b1b6))
* **db:** take gene_id from the annotation, add hgnc_id ([#1085](https://github.com/varfish-org/mehari/issues/1085)) ([776f2e0](https://github.com/varfish-org/mehari/commit/776f2e035147c68d54e39e2db46af2f13d4692ea))


### Bug Fixes

* **annotate:** correct indel, UTR and start codon annotations ([#1067](https://github.com/varfish-org/mehari/issues/1067)) ([3953d9d](https://github.com/varfish-org/mehari/commit/3953d9dc23d6761a7951aaa40e00d49111dd6448))
* **annotate:** decide insertions at exon edges by sequence ([#1114](https://github.com/varfish-org/mehari/issues/1114)) ([62373fa](https://github.com/varfish-org/mehari/commit/62373fa065bb0b503a90166916e8ea782ad0a62b))
* **annotate:** take the terms of stop codon deletions from the protein ([#1113](https://github.com/varfish-org/mehari/issues/1113)) ([d0ee375](https://github.com/varfish-org/mehari/commit/d0ee375c9fc32053e77ea3db8385c529aa92aefb))
* correct edge cases found in the final review of the stack ([#1117](https://github.com/varfish-org/mehari/issues/1117)) ([b861941](https://github.com/varfish-org/mehari/commit/b861941ba466d8795100904dab035d6bbbf6da0d))
* **python:** bump pyo3 to 0.29 for advisories ([#1077](https://github.com/varfish-org/mehari/issues/1077)) ([c45846a](https://github.com/varfish-org/mehari/commit/c45846a50ebc0d701b75c8f620da15a384755b1d))


### Performance Improvements

* **annotate:** take the region from c_to_p and stop copying for plugins ([#1118](https://github.com/varfish-org/mehari/issues/1118)) ([02e55a4](https://github.com/varfish-org/mehari/commit/02e55a4471709815c19a4f1992d93a79a985d94c))


### Code Refactoring

* restrict mehari db to transcripts ([#1033](https://github.com/varfish-org/mehari/issues/1033)) ([#1074](https://github.com/varfish-org/mehari/issues/1074)) ([ffb923c](https://github.com/varfish-org/mehari/commit/ffb923cca27d831faedad13178c6a6f6c6692af3))


### Build System

* **deps:** bump hgvs to 2.0 and the MSRV to 1.92 ([#1080](https://github.com/varfish-org/mehari/issues/1080)) ([2228a2a](https://github.com/varfish-org/mehari/commit/2228a2a7e888234390aefec5375ec9962fbae513))

## [0.46.0](https://github.com/varfish-org/mehari/compare/mehari-python-v0.45.1...mehari-python-v0.46.0) (2026-09-22)


### Features

* **mehari-python:** add tqdm-compatible build progress ([#1055](https://github.com/varfish-org/mehari/issues/1055)) ([d54313d](https://github.com/varfish-org/mehari/commit/d54313d5e2730c8ef22d4184a0d2e4cf016c121f))

## [0.45.1](https://github.com/varfish-org/mehari/compare/mehari-python-v0.45.0...mehari-python-v0.45.1) (2026-07-03)


### Miscellaneous Chores

* **mehari-python:** Synchronize mehari-sync versions

## [0.45.0](https://github.com/varfish-org/mehari/compare/mehari-python-v0.44.1...mehari-python-v0.45.0) (2026-07-02)


### ⚠ BREAKING CHANGES

* support CADD, spliceAI, dbSNP and building custom/generic lookup annotation DBs ([#1025](https://github.com/varfish-org/mehari/issues/1025))

### Features

* support CADD, spliceAI, dbSNP and building custom/generic lookup annotation DBs ([#1025](https://github.com/varfish-org/mehari/issues/1025)) ([7abb129](https://github.com/varfish-org/mehari/commit/7abb12903fde91bd4341c10726e2ca46914769a6))

## [0.44.1](https://github.com/varfish-org/mehari/compare/mehari-python-v0.44.0...mehari-python-v0.44.1) (2026-06-17)


### Miscellaneous Chores

* **mehari-python:** Synchronize mehari-sync versions

## [0.44.0](https://github.com/varfish-org/mehari/compare/mehari-python-v0.43.3...mehari-python-v0.44.0) (2026-06-15)


### Miscellaneous Chores

* **mehari-python:** Synchronize mehari-sync versions

## [0.43.3](https://github.com/varfish-org/mehari/compare/mehari-python-v0.43.2...mehari-python-v0.43.3) (2026-06-10)


### ⚠ BREAKING CHANGES

* add experimental option to predict compound effect of multiple variants on the same transcript, allow building txdb from gff3+fasta ([#972](https://github.com/varfish-org/mehari/issues/972))

### Features

* add experimental option to predict compound effect of multiple variants on the same transcript, allow building txdb from gff3+fasta ([#972](https://github.com/varfish-org/mehari/issues/972)) ([e377c4e](https://github.com/varfish-org/mehari/commit/e377c4eeeb7da1644feba7d37be7a71dc86ada69))
* add options to include reference/alternative cDNA and protein sequences in output ([#966](https://github.com/varfish-org/mehari/issues/966)) ([ec87cb0](https://github.com/varfish-org/mehari/commit/ec87cb0e5c8bf6c70672a5d2b8fd5f76da82c0b5))
* Add python bindings via pyO3 ([#959](https://github.com/varfish-org/mehari/issues/959)) ([b2450be](https://github.com/varfish-org/mehari/commit/b2450be38bdf65aa0c21c5b77ea6770910bd46a5))


### Bug Fixes

* hide server functionality behind feature gate (enabled by default) ([#964](https://github.com/varfish-org/mehari/issues/964)) ([0031153](https://github.com/varfish-org/mehari/commit/00311532ac5d296277f85220d415d5a1cc96dbd3))
* include polars and pyarrow in mehari-python dependencies ([#962](https://github.com/varfish-org/mehari/issues/962)) ([23b75c4](https://github.com/varfish-org/mehari/commit/23b75c4331db943d53bcc6b6f832a8e7cd09ac8f))
* more robust freqs/clinvar assembly checks ([#1008](https://github.com/varfish-org/mehari/issues/1008)) ([d0ead01](https://github.com/varfish-org/mehari/commit/d0ead012210ae2daa8947f97e6130d0b81922961))
* pick mode first and length picking ([#993](https://github.com/varfish-org/mehari/issues/993)) ([e63b7b0](https://github.com/varfish-org/mehari/commit/e63b7b0c334b092e3a2f63a37aebe861a949697a))
* varfish postprocess TSV assembly capitalization ([#1007](https://github.com/varfish-org/mehari/issues/1007)) ([d3c997d](https://github.com/varfish-org/mehari/commit/d3c997d01ff5da20a493ad3c621e777e7806c988))


### Miscellaneous Chores

* bump version for annonars update ([1e625f9](https://github.com/varfish-org/mehari/commit/1e625f9a3a90854427669609296b93e8b7c964ed))
* bump version for jemalloc default feature ([517e1d3](https://github.com/varfish-org/mehari/commit/517e1d3236d39d8670e57d6870bb1f48b6dcf4e0))
* Release-As: 0.43.2 ([e170371](https://github.com/varfish-org/mehari/commit/e170371ae6fcd5cc4cab847d20433087fa474513))

## [0.43.2](https://github.com/varfish-org/mehari/compare/v0.43.0...v0.43.2) (2026-05-07)


### ⚠ BREAKING CHANGES

* add experimental option to predict compound effect of multiple variants on the same transcript, allow building txdb from gff3+fasta ([#972](https://github.com/varfish-org/mehari/issues/972))

### Features

* add experimental option to predict compound effect of multiple variants on the same transcript, allow building txdb from gff3+fasta ([#972](https://github.com/varfish-org/mehari/issues/972)) ([e377c4e](https://github.com/varfish-org/mehari/commit/e377c4eeeb7da1644feba7d37be7a71dc86ada69))
* add options to include reference/alternative cDNA and protein sequences in output ([#966](https://github.com/varfish-org/mehari/issues/966)) ([ec87cb0](https://github.com/varfish-org/mehari/commit/ec87cb0e5c8bf6c70672a5d2b8fd5f76da82c0b5))
* Add python bindings via pyO3 ([#959](https://github.com/varfish-org/mehari/issues/959)) ([b2450be](https://github.com/varfish-org/mehari/commit/b2450be38bdf65aa0c21c5b77ea6770910bd46a5))


### Bug Fixes

* hide server functionality behind feature gate (enabled by default) ([#964](https://github.com/varfish-org/mehari/issues/964)) ([0031153](https://github.com/varfish-org/mehari/commit/00311532ac5d296277f85220d415d5a1cc96dbd3))
* include polars and pyarrow in mehari-python dependencies ([#962](https://github.com/varfish-org/mehari/issues/962)) ([23b75c4](https://github.com/varfish-org/mehari/commit/23b75c4331db943d53bcc6b6f832a8e7cd09ac8f))
* pick mode first and length picking ([#993](https://github.com/varfish-org/mehari/issues/993)) ([e63b7b0](https://github.com/varfish-org/mehari/commit/e63b7b0c334b092e3a2f63a37aebe861a949697a))


### Miscellaneous Chores

* bump version for annonars update ([1e625f9](https://github.com/varfish-org/mehari/commit/1e625f9a3a90854427669609296b93e8b7c964ed))
* bump version for jemalloc default feature ([517e1d3](https://github.com/varfish-org/mehari/commit/517e1d3236d39d8670e57d6870bb1f48b6dcf4e0))
* Release-As: 0.43.2 ([e170371](https://github.com/varfish-org/mehari/commit/e170371ae6fcd5cc4cab847d20433087fa474513))

## [0.35.0](https://github.com/varfish-org/mehari/compare/v0.43.0...v0.35.0) (2026-05-07)


### ⚠ BREAKING CHANGES

* add experimental option to predict compound effect of multiple variants on the same transcript, allow building txdb from gff3+fasta ([#972](https://github.com/varfish-org/mehari/issues/972))

### Features

* add experimental option to predict compound effect of multiple variants on the same transcript, allow building txdb from gff3+fasta ([#972](https://github.com/varfish-org/mehari/issues/972)) ([e377c4e](https://github.com/varfish-org/mehari/commit/e377c4eeeb7da1644feba7d37be7a71dc86ada69))
* add options to include reference/alternative cDNA and protein sequences in output ([#966](https://github.com/varfish-org/mehari/issues/966)) ([ec87cb0](https://github.com/varfish-org/mehari/commit/ec87cb0e5c8bf6c70672a5d2b8fd5f76da82c0b5))
* Add python bindings via pyO3 ([#959](https://github.com/varfish-org/mehari/issues/959)) ([b2450be](https://github.com/varfish-org/mehari/commit/b2450be38bdf65aa0c21c5b77ea6770910bd46a5))


### Bug Fixes

* hide server functionality behind feature gate (enabled by default) ([#964](https://github.com/varfish-org/mehari/issues/964)) ([0031153](https://github.com/varfish-org/mehari/commit/00311532ac5d296277f85220d415d5a1cc96dbd3))
* include polars and pyarrow in mehari-python dependencies ([#962](https://github.com/varfish-org/mehari/issues/962)) ([23b75c4](https://github.com/varfish-org/mehari/commit/23b75c4331db943d53bcc6b6f832a8e7cd09ac8f))
* pick mode first and length picking ([#993](https://github.com/varfish-org/mehari/issues/993)) ([e63b7b0](https://github.com/varfish-org/mehari/commit/e63b7b0c334b092e3a2f63a37aebe861a949697a))


### Miscellaneous Chores

* bump version for annonars update ([1e625f9](https://github.com/varfish-org/mehari/commit/1e625f9a3a90854427669609296b93e8b7c964ed))
* bump version for jemalloc default feature ([517e1d3](https://github.com/varfish-org/mehari/commit/517e1d3236d39d8670e57d6870bb1f48b6dcf4e0))

## [0.42.1](https://github.com/varfish-org/mehari/compare/v0.42.0...v0.42.1) (2026-05-06)


### Bug Fixes

* pick mode first and length picking ([#993](https://github.com/varfish-org/mehari/issues/993)) ([e63b7b0](https://github.com/varfish-org/mehari/commit/e63b7b0c334b092e3a2f63a37aebe861a949697a))

## [0.42.0](https://github.com/varfish-org/mehari/compare/v0.41.2...v0.42.0) (2026-04-14)


### ⚠ BREAKING CHANGES

* add experimental option to predict compound effect of multiple variants on the same transcript, allow building txdb from gff3+fasta ([#972](https://github.com/varfish-org/mehari/issues/972))

### Features

* add experimental option to predict compound effect of multiple variants on the same transcript, allow building txdb from gff3+fasta ([#972](https://github.com/varfish-org/mehari/issues/972)) ([e377c4e](https://github.com/varfish-org/mehari/commit/e377c4eeeb7da1644feba7d37be7a71dc86ada69))
* add options to include reference/alternative cDNA and protein sequences in output ([#966](https://github.com/varfish-org/mehari/issues/966)) ([ec87cb0](https://github.com/varfish-org/mehari/commit/ec87cb0e5c8bf6c70672a5d2b8fd5f76da82c0b5))

## [0.41.2](https://github.com/varfish-org/mehari/compare/v0.41.1...v0.41.2) (2026-03-22)


### Bug Fixes

* hide server functionality behind feature gate (enabled by default) ([#964](https://github.com/varfish-org/mehari/issues/964)) ([0031153](https://github.com/varfish-org/mehari/commit/00311532ac5d296277f85220d415d5a1cc96dbd3))

## [0.41.1](https://github.com/varfish-org/mehari/compare/v0.41.0...v0.41.1) (2026-03-22)


### Bug Fixes

* include polars and pyarrow in mehari-python dependencies ([#962](https://github.com/varfish-org/mehari/issues/962)) ([23b75c4](https://github.com/varfish-org/mehari/commit/23b75c4331db943d53bcc6b6f832a8e7cd09ac8f))

## [0.41.0](https://github.com/varfish-org/mehari/compare/v0.40.0...v0.41.0) (2026-03-21)


### Features

* Add python bindings via pyO3 ([#959](https://github.com/varfish-org/mehari/issues/959)) ([b2450be](https://github.com/varfish-org/mehari/commit/b2450be38bdf65aa0c21c5b77ea6770910bd46a5))
