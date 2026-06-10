# Changelog

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
