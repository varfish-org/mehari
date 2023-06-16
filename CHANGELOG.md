# Changelog

## [0.5.3](https://github.com/bihealth/mehari/compare/v0.5.2...v0.5.3) (2023-06-16)


### Miscellaneous Chores

* re-release as v0.5.3 ([12aed0a](https://github.com/bihealth/mehari/commit/12aed0aa98e8714c02425cf977faf33eb7d6c312))

## [0.5.2](https://github.com/bihealth/mehari/compare/v0.5.1...v0.5.2) (2023-06-14)


### Bug Fixes

* docker build in CI ([#96](https://github.com/bihealth/mehari/issues/96)) ([ba34c03](https://github.com/bihealth/mehari/commit/ba34c034c17849b729dbc5faa9fd57bbaba913fa))

## [0.5.1](https://github.com/bihealth/mehari/compare/v0.5.0...v0.5.1) (2023-06-13)


### Features

* build Docker image in CI ([#92](https://github.com/bihealth/mehari/issues/92)) ([c72d804](https://github.com/bihealth/mehari/commit/c72d80494043e46cf80113ccd984984d730e9368))


### Bug Fixes

* add Cargo.lock to git ([#94](https://github.com/bihealth/mehari/issues/94)) ([f7c4f03](https://github.com/bihealth/mehari/commit/f7c4f03628e8657a510f471245dfbc6d3c0382c8))
* remove unused dependencies ([#95](https://github.com/bihealth/mehari/issues/95)) ([a1e1054](https://github.com/bihealth/mehari/commit/a1e10542403339c4e7e421e74eef95e74ecf5468))

## [0.5.0](https://github.com/bihealth/mehari/compare/v0.4.1...v0.5.0) (2023-06-12)


### Features

* REST server for variant consequence prediction ([#18](https://github.com/bihealth/mehari/issues/18)) ([#90](https://github.com/bihealth/mehari/issues/90)) ([b8c5a7b](https://github.com/bihealth/mehari/commit/b8c5a7b6ec60e03fe7d973f96cab6ac2e628e7ea))

## [0.4.1](https://github.com/bihealth/mehari/compare/v0.4.0...v0.4.1) (2023-06-12)


### Bug Fixes

* properly flush stream in "db create txs" ([#88](https://github.com/bihealth/mehari/issues/88)) ([95458f3](https://github.com/bihealth/mehari/commit/95458f3e0e8d598bfdb1ad24e8ae6bb930a725a9))

## [0.4.0](https://github.com/bihealth/mehari/compare/v0.3.1...v0.4.0) (2023-06-08)


### Features

* move clinvar database code (import etc.) to annonars ([#86](https://github.com/bihealth/mehari/issues/86)) ([457be46](https://github.com/bihealth/mehari/commit/457be468591ee56d7b0c6afabde4e6ef18611d31))
* switch to annonars for frequency db creation ([#81](https://github.com/bihealth/mehari/issues/81),[#77](https://github.com/bihealth/mehari/issues/77),[#20](https://github.com/bihealth/mehari/issues/20)) ([#83](https://github.com/bihealth/mehari/issues/83)) ([bb27686](https://github.com/bihealth/mehari/commit/bb2768686cbd96d20d064c58eaadbc18e09ce35c))


### Code Refactoring

* move from linked-hash-map to indexmap ([#82](https://github.com/bihealth/mehari/issues/82)) ([#87](https://github.com/bihealth/mehari/issues/87)) ([f9e354d](https://github.com/bihealth/mehari/commit/f9e354dbcbe32e25c2a41d90ad1d0b2874332a06))

## [0.3.1](https://github.com/bihealth/mehari/compare/v0.3.0...v0.3.1) (2023-05-04)


### Bug Fixes

* use copy number to fill in genotypes for CNVs ([#72](https://github.com/bihealth/mehari/issues/72)) ([5585a08](https://github.com/bihealth/mehari/commit/5585a084ea10de0612334530ff9a4a55b60518d3))

## [0.3.0](https://github.com/bihealth/mehari/compare/v0.2.1...v0.3.0) (2023-05-03)


### Features

* add genotype counting ([#70](https://github.com/bihealth/mehari/issues/70)) ([c40f83c](https://github.com/bihealth/mehari/commit/c40f83ccfa5bd669e3acec2280084b072738695b))


### Bug Fixes

* only copy samples from input VCF header ([#69](https://github.com/bihealth/mehari/issues/69)) ([316b123](https://github.com/bihealth/mehari/commit/316b1237648435121e3d2c877426dc36ca31bcdd))

## [0.2.1](https://github.com/bihealth/mehari/compare/v0.2.0...v0.2.1) (2023-05-02)


### Bug Fixes

* adjust non-adjacent exon recognition ([#67](https://github.com/bihealth/mehari/issues/67)) ([96437e4](https://github.com/bihealth/mehari/commit/96437e4a69349cabc7a165a130fd974a133fe307))
* allow multiple bases in BND ALT ([#63](https://github.com/bihealth/mehari/issues/63)) ([e5b93d5](https://github.com/bihealth/mehari/commit/e5b93d5727c16cd0881cfcc4e59777a6b9ab5923))
* properly recognize .gz/.zst extension ([#66](https://github.com/bihealth/mehari/issues/66)) ([6960521](https://github.com/bihealth/mehari/commit/696052175c1dbd5d9cc48026a4abddc33f81fd23))
* skipping "*" ALT alleles ([#65](https://github.com/bihealth/mehari/issues/65)) ([95b2b61](https://github.com/bihealth/mehari/commit/95b2b6170f1180545a7e3d327ee61d32ebf6ba59))

## [0.2.0](https://github.com/bihealth/mehari/compare/v0.1.1...v0.2.0) (2023-04-24)


### Features

* allow annotation from maelstorm cov/mq BED files ([#51](https://github.com/bihealth/mehari/issues/51)) ([#55](https://github.com/bihealth/mehari/issues/55)) ([f3966ac](https://github.com/bihealth/mehari/commit/f3966acb9a7e0d81be2f6ece51422766bafcca09))
* annotation of structural variants ([#46](https://github.com/bihealth/mehari/issues/46)) ([#52](https://github.com/bihealth/mehari/issues/52)) ([8da1d38](https://github.com/bihealth/mehari/commit/8da1d388e21d6214e0cbe96ac54f703d9a736638))
* change database path setup ([#48](https://github.com/bihealth/mehari/issues/48)) ([#60](https://github.com/bihealth/mehari/issues/60)) ([b26671f](https://github.com/bihealth/mehari/commit/b26671fb5c1b7ede8f0ede58c97d778586b761a0))
* implement command to copy rocksdb database ([#45](https://github.com/bihealth/mehari/issues/45)) ([#47](https://github.com/bihealth/mehari/issues/47)) ([3c8543f](https://github.com/bihealth/mehari/commit/3c8543fbae281e8f9ea41bbad1718608fc9f00a6))
* switch from flatbuffers to protobuf ([#15](https://github.com/bihealth/mehari/issues/15)) ([#57](https://github.com/bihealth/mehari/issues/57)) ([3fe3322](https://github.com/bihealth/mehari/commit/3fe332246fc5730a31d8a355bc92ba106e956db8))


### Bug Fixes

* issue on negative CDS position ([#49](https://github.com/bihealth/mehari/issues/49)) ([#58](https://github.com/bihealth/mehari/issues/58)) ([9a63e30](https://github.com/bihealth/mehari/commit/9a63e3004e044c65c6552e3698b3a3daeeab3a6f))

## [0.1.1](https://github.com/bihealth/mehari/compare/v0.1.0...v0.1.1) (2023-04-06)


### Bug Fixes

* release 0.1.1 ([5565d93](https://github.com/bihealth/mehari/commit/5565d93e082de77f1f9b532bd581c827390e4b68))

## 0.1.0 (2023-04-06)


### Features

* adding command "db create seqvars-clinvar" ([#24](https://github.com/bihealth/mehari/issues/24)) ([#25](https://github.com/bihealth/mehari/issues/25)) ([fd4c523](https://github.com/bihealth/mehari/commit/fd4c5238c2ff5b56577642ed4447b0abbf84d739))
* allow mehari to write out VarFish compatible TSV ([#34](https://github.com/bihealth/mehari/issues/34)) ([#35](https://github.com/bihealth/mehari/issues/35)) ([72c6832](https://github.com/bihealth/mehari/commit/72c6832636c4a9a632439a2d0d9ebc2d21192209))
* command "verify seqvars" allows to compare to VEP ([#21](https://github.com/bihealth/mehari/issues/21)) ([#37](https://github.com/bihealth/mehari/issues/37)) ([b8761e3](https://github.com/bihealth/mehari/commit/b8761e3f255b1ef9988e44d61dea5bd099bfaf44))
* frequency annotation with database ([#12](https://github.com/bihealth/mehari/issues/12)) ([#13](https://github.com/bihealth/mehari/issues/13)) ([f63f8c2](https://github.com/bihealth/mehari/commit/f63f8c2f7e10487a126557bcbf8d7853ef418e7c))
* implement seqvar frequency db construction ([#2](https://github.com/bihealth/mehari/issues/2)) ([#5](https://github.com/bihealth/mehari/issues/5)) ([01fc64c](https://github.com/bihealth/mehari/commit/01fc64c2f8aa351768bd5c1703a88d7a7cb021a8))
* implementing "annotate seqvars" ([#3](https://github.com/bihealth/mehari/issues/3)) ([#7](https://github.com/bihealth/mehari/issues/7)) ([6400ed1](https://github.com/bihealth/mehari/commit/6400ed1d79ba736658549572ba0b16b6fe4626c7))
* transcript based annotation ([#11](https://github.com/bihealth/mehari/issues/11)) ([#14](https://github.com/bihealth/mehari/issues/14)) ([6677899](https://github.com/bihealth/mehari/commit/66778992274107e9a632870ead0b6d97161f7b7e))
* transcript database building ([#1](https://github.com/bihealth/mehari/issues/1)) ([#8](https://github.com/bihealth/mehari/issues/8)) ([bcff954](https://github.com/bihealth/mehari/commit/bcff9546b535fbd98caa35f89102cf6539be33b3))


### Bug Fixes

* allow loading JSON from .gz files ([#22](https://github.com/bihealth/mehari/issues/22)) ([cefd950](https://github.com/bihealth/mehari/commit/cefd950d823f40f94b095358a1ce7cdd7d834843))
* issues occuring when annotating real data ([#28](https://github.com/bihealth/mehari/issues/28)) ([0636d49](https://github.com/bihealth/mehari/commit/0636d490d2a17a0e8688c85d5b018ab1d0ae117f))

## Changelog
