# Changelog

## [0.36.3](https://github.com/varfish-org/mehari/compare/v0.36.2...v0.36.3) (2025-07-09)


### Bug Fixes

* keep longest transcript per TranscriptSource ([#790](https://github.com/varfish-org/mehari/issues/790)) ([9d12e6b](https://github.com/varfish-org/mehari/commit/9d12e6b29a0b700397e9981253c6b80618e0b661))
* never skip records during TSV output ([#789](https://github.com/varfish-org/mehari/issues/789)) ([97827f6](https://github.com/varfish-org/mehari/commit/97827f69917ac0e08b5a6d296f4f7452c433af1e))

## [0.36.2](https://github.com/varfish-org/mehari/compare/v0.36.1...v0.36.2) (2025-07-02)


### Bug Fixes

* include both refseq and ensembl info in TSV output ([#786](https://github.com/varfish-org/mehari/issues/786)) ([cba5a25](https://github.com/varfish-org/mehari/commit/cba5a252484f6868c4cd4d91151a0d16b9a14201))

## [0.36.1](https://github.com/varfish-org/mehari/compare/v0.36.0...v0.36.1) (2025-07-01)


### Bug Fixes

* TSV output + --keep-intergenic ([#785](https://github.com/varfish-org/mehari/issues/785)) ([f44a1e5](https://github.com/varfish-org/mehari/commit/f44a1e59db2d6406546b39e148534ba7623b338c))
* Use noodle's VCF FileFormat::default() ([#778](https://github.com/varfish-org/mehari/issues/778)) ([3be028d](https://github.com/varfish-org/mehari/commit/3be028d1e25e365045198a3043f64297f087b43b))

## [0.36.0](https://github.com/varfish-org/mehari/compare/v0.35.1...v0.36.0) (2025-06-17)


### Features

* Harmonize ensembl and refseq `annotate seqvars` TSV output ([#773](https://github.com/varfish-org/mehari/issues/773)) ([b608c8e](https://github.com/varfish-org/mehari/commit/b608c8e5ca656dec9a3be823312a0e4d235238b7))


### Bug Fixes

* allow multiple references in mehari server run ([#749](https://github.com/varfish-org/mehari/issues/749)) ([ddf7d3f](https://github.com/varfish-org/mehari/commit/ddf7d3f84e2501c740d9dcfc83977f190cdc084b))

## [0.35.1](https://github.com/varfish-org/mehari/compare/v0.35.0...v0.35.1) (2025-04-24)


### Bug Fixes

* add 'symbol' as alias for 'gene_symbol' for hgnc complete set TSV ([#748](https://github.com/varfish-org/mehari/issues/748)) ([5d50862](https://github.com/varfish-org/mehari/commit/5d508620aaf96e342ee796b7e41336216d777108))
* missing /rocksdb dir in docker entrypoint script ([#745](https://github.com/varfish-org/mehari/issues/745)) ([bdc9a3c](https://github.com/varfish-org/mehari/commit/bdc9a3c5fa35bf90dbdcd7b1134ccdfcb1439728))
* point entrypoint at 'frequencies', not 'freqs' ([#747](https://github.com/varfish-org/mehari/issues/747)) ([23fbeeb](https://github.com/varfish-org/mehari/commit/23fbeebb0aa992b62b72e5dd15d67494caa53859))

## [0.35.0](https://github.com/varfish-org/mehari/compare/v0.34.0...v0.35.0) (2025-04-16)


### Miscellaneous Chores

* bump version for jemalloc default feature ([517e1d3](https://github.com/varfish-org/mehari/commit/517e1d3236d39d8670e57d6870bb1f48b6dcf4e0))

## [0.34.0](https://github.com/varfish-org/mehari/compare/v0.33.1...v0.34.0) (2025-04-16)


### Miscellaneous Chores

* bump version for annonars update ([1e625f9](https://github.com/varfish-org/mehari/commit/1e625f9a3a90854427669609296b93e8b7c964ed))

## [0.33.1](https://github.com/varfish-org/mehari/compare/v0.33.0...v0.33.1) (2025-04-14)


### Bug Fixes

* update noodles to 0.97 and annonars to 0.42.5 ([#735](https://github.com/varfish-org/mehari/issues/735)) ([f106408](https://github.com/varfish-org/mehari/commit/f1064087b2296f534d56cce7ecf7bf2527d203ad))

## [0.33.0](https://github.com/varfish-org/mehari/compare/v0.32.0...v0.33.0) (2025-04-01)


### Features

* Miscellaneous consequence prediction fixes ([#704](https://github.com/varfish-org/mehari/issues/704)) ([20b16ee](https://github.com/varfish-org/mehari/commit/20b16eeb2992b1dafb3a03f9fd6110a4366dc555))
* More flexible INFO/ANN header entry ([#727](https://github.com/varfish-org/mehari/issues/727)) ([2062c3c](https://github.com/varfish-org/mehari/commit/2062c3cbf8db9b369ac105bb78ea570b6338ad80))


### Bug Fixes

* report `transcript_ablation` for deletions starting pre start and ending post stop ([#702](https://github.com/varfish-org/mehari/issues/702)) ([e71a153](https://github.com/varfish-org/mehari/commit/e71a153953140e2783a22220045d6d9761d63929))

## [0.32.0](https://github.com/varfish-org/mehari/compare/v0.31.0...v0.32.0) (2025-02-25)


### Features

* allow specifying reference FASTA ([#684](https://github.com/varfish-org/mehari/issues/684)) ([a9dc3da](https://github.com/varfish-org/mehari/commit/a9dc3da2b22f0a1e899d291fc6b14c8c44f640f7))


### Bug Fixes

* improve consequence predictions ([#690](https://github.com/varfish-org/mehari/issues/690)) ([73d0fd3](https://github.com/varfish-org/mehari/commit/73d0fd3553cb98dd436cb225c8c392a919dbd7e2))

## [0.31.0](https://github.com/varfish-org/mehari/compare/v0.30.1...v0.31.0) (2025-02-06)


### Features

* Add rare amino acid SO term ([#658](https://github.com/varfish-org/mehari/issues/658)) ([aee7693](https://github.com/varfish-org/mehari/commit/aee76930c633ab1be6fb023da042ca3fa1255ae4))
* refactor db subset CLI, add subset by VCF and subset by TxId options ([#641](https://github.com/varfish-org/mehari/issues/641)) ([9286286](https://github.com/varfish-org/mehari/commit/92862866b41307b56c2823b93b3156a876f0bf74))


### Bug Fixes

* Consolidate CLI ([#650](https://github.com/varfish-org/mehari/issues/650)) ([08c8fb6](https://github.com/varfish-org/mehari/commit/08c8fb61de2453ba411266554dab052f7f15e9a7))
* ensure FORMAT sample order is in-sync with header sample order ([#671](https://github.com/varfish-org/mehari/issues/671)) ([1d3c230](https://github.com/varfish-org/mehari/commit/1d3c230a48e68694a9bb840f72f924b9dfed0503))
* include mehari cmd and version info in VCF header ([#673](https://github.com/varfish-org/mehari/issues/673)) ([8b9fa2e](https://github.com/varfish-org/mehari/commit/8b9fa2e91485cd0b0cede083ae679d37d65f62bd))
* inframe vs missense ([#659](https://github.com/varfish-org/mehari/issues/659)) ([36bd717](https://github.com/varfish-org/mehari/commit/36bd7177fc23c55b6c912e14d44dc1355520ac5f))

## [0.30.1](https://github.com/varfish-org/mehari/compare/v0.30.0...v0.30.1) (2024-11-21)


### Bug Fixes

* OpenAPI for /api/v1/genes/transcripts ([#628](https://github.com/varfish-org/mehari/issues/628)) ([#629](https://github.com/varfish-org/mehari/issues/629)) ([d155620](https://github.com/varfish-org/mehari/commit/d1556200b712f5751680d45c947d68e0ebf29e41))

## [0.30.0](https://github.com/varfish-org/mehari/compare/v0.29.6...v0.30.0) (2024-11-20)


### ⚠ BREAKING CHANGES

* include source information in transcript database ([#615](https://github.com/varfish-org/mehari/issues/615))
* provide sub command "server run" ([#600](https://github.com/varfish-org/mehari/issues/600)) (#602)

### Features

* /api/v1/seqvars/csq endpoint with OpenAPI ([#606](https://github.com/varfish-org/mehari/issues/606)) ([#614](https://github.com/varfish-org/mehari/issues/614)) ([2fdf76e](https://github.com/varfish-org/mehari/commit/2fdf76e6fe7a928f14ba96fe8b5a8c85cfc29fe3))
* /api/v1/strucvars/csq endpoint with OpenAPI ([#607](https://github.com/varfish-org/mehari/issues/607)) ([#612](https://github.com/varfish-org/mehari/issues/612)) ([cce8a6d](https://github.com/varfish-org/mehari/commit/cce8a6dcf65b2d5d0f3debe9351698efd0c8af86))
* include source information in transcript database ([#615](https://github.com/varfish-org/mehari/issues/615)) ([111ffcc](https://github.com/varfish-org/mehari/commit/111ffcca83c5847647da36ea11b91dc736f7cf2e))
* provide genes/transcripts endpoint with openapi ([#605](https://github.com/varfish-org/mehari/issues/605)) ([#610](https://github.com/varfish-org/mehari/issues/610)) ([9132312](https://github.com/varfish-org/mehari/commit/91323124dcb872719a07b4aeec309667e91b540a))
* provide sub command "server run" ([#600](https://github.com/varfish-org/mehari/issues/600)) ([#602](https://github.com/varfish-org/mehari/issues/602)) ([48382b6](https://github.com/varfish-org/mehari/commit/48382b668a8354d78763fd495403b56f1c04d114))


### Bug Fixes

* refactor csq analysis code ([#599](https://github.com/varfish-org/mehari/issues/599)) ([faffad6](https://github.com/varfish-org/mehari/commit/faffad63b57e42eba76b75d886745ef45cfc2ef0))

## [0.29.6](https://github.com/varfish-org/mehari/compare/v0.29.5...v0.29.6) (2024-10-16)


### Bug Fixes

* add missing librocksdb8.9 to Dockerfile ([#584](https://github.com/varfish-org/mehari/issues/584)) ([f9bd018](https://github.com/varfish-org/mehari/commit/f9bd018067ee0ad794115b5e90c16e909735652d))
* skip REF-only records ([#576](https://github.com/varfish-org/mehari/issues/576)) ([688b968](https://github.com/varfish-org/mehari/commit/688b968e0a947c0eec7a918b60a77bca26c164f5))

## [0.29.5](https://github.com/varfish-org/mehari/compare/v0.29.4...v0.29.5) (2024-10-08)


### Continuous Integration

* adding dependency libsqlite3-dev ([#574](https://github.com/varfish-org/mehari/issues/574)) ([08bdc10](https://github.com/varfish-org/mehari/commit/08bdc1020de37acb9dd3c02f614cbac441cea379))

## [0.29.4](https://github.com/varfish-org/mehari/compare/v0.29.3...v0.29.4) (2024-10-08)


### Continuous Integration

* bump ubuntu to 24.04 ([#572](https://github.com/varfish-org/mehari/issues/572)) ([0ef8f7d](https://github.com/varfish-org/mehari/commit/0ef8f7dd6d78eeb980c8cf66eb1ade09bb82dc41))

## [0.29.3](https://github.com/varfish-org/mehari/compare/v0.29.2...v0.29.3) (2024-10-08)


### Bug Fixes

* install libsnappy-dev in CI ([#570](https://github.com/varfish-org/mehari/issues/570)) ([1169479](https://github.com/varfish-org/mehari/commit/11694798a295096b29503eddc513ec6446f9a0b8))

## [0.29.2](https://github.com/varfish-org/mehari/compare/v0.29.1...v0.29.2) (2024-10-08)


### Continuous Integration

* install rocksdb dependencies in CI ([#568](https://github.com/varfish-org/mehari/issues/568)) ([c1dfe6f](https://github.com/varfish-org/mehari/commit/c1dfe6f9f3891a5a490af564d20102f675a0b73b))

## [0.29.1](https://github.com/varfish-org/mehari/compare/v0.29.0...v0.29.1) (2024-10-08)


### Bug Fixes

* percent encoding clinvar annotation ([#553](https://github.com/varfish-org/mehari/issues/553)) ([#566](https://github.com/varfish-org/mehari/issues/566)) ([6540520](https://github.com/varfish-org/mehari/commit/6540520f16daf8094772673a5d89ec8ca42ab0f4))

## [0.29.0](https://github.com/varfish-org/mehari/compare/v0.28.1...v0.29.0) (2024-10-07)


### Features

* Add catch all tag "Other" ([#558](https://github.com/varfish-org/mehari/issues/558)) ([c709a36](https://github.com/varfish-org/mehari/commit/c709a3693b12d93151898b1aa1335bb6cf85d839))

## [0.28.1](https://github.com/varfish-org/mehari/compare/v0.28.0...v0.28.1) (2024-09-02)


### Bug Fixes

* disregard tsv output order in brca2_zar1l test ([#529](https://github.com/varfish-org/mehari/issues/529)) ([d5e1066](https://github.com/varfish-org/mehari/commit/d5e10662d1a9ab89bdabca845d9b5976a9e7be81))
* make seqvar annotators public for varfish-server-worker ([#528](https://github.com/varfish-org/mehari/issues/528)) ([22610d1](https://github.com/varfish-org/mehari/commit/22610d1cb017620d82826f1554c58d3b928e8ad9))

## [0.28.0](https://github.com/varfish-org/mehari/compare/v0.27.0...v0.28.0) (2024-08-19)


### Features

* add "db merge" subcommand ([#517](https://github.com/varfish-org/mehari/issues/517)) ([3c7b514](https://github.com/varfish-org/mehari/commit/3c7b5148803cc95568a709e548ffd98d95256e55))

## [0.27.0](https://github.com/varfish-org/mehari/compare/v0.26.1...v0.27.0) (2024-08-16)


### Features

* refactor annotate seqvars cli ([#510](https://github.com/varfish-org/mehari/issues/510)) ([915b30a](https://github.com/varfish-org/mehari/commit/915b30a6ba202998870ebe7395572df159ab72cc))


### Bug Fixes

* codon deletion right in front of stop codon is not a stop_gained ([#515](https://github.com/varfish-org/mehari/issues/515)) ([48fd896](https://github.com/varfish-org/mehari/commit/48fd896966d9dc4769a64a123131c66f0a1052c6))

## [0.26.1](https://github.com/varfish-org/mehari/compare/v0.26.0...v0.26.1) (2024-08-07)


### Bug Fixes

* Add gencode primary tag ([#509](https://github.com/varfish-org/mehari/issues/509)) ([a6c8566](https://github.com/varfish-org/mehari/commit/a6c856680f28b6b6f9ef08d685e6aeba9a25320b))
* for fake-rna-* transcripts for refseq MT, add fake version as well ([#507](https://github.com/varfish-org/mehari/issues/507)) ([f296e5d](https://github.com/varfish-org/mehari/commit/f296e5dbbf089614411659e08e2dac363478a1b5))

## [0.26.0](https://github.com/varfish-org/mehari/compare/v0.25.7...v0.26.0) (2024-08-02)


### ⚠ BREAKING CHANGES

* refactor DB create ([#456](https://github.com/varfish-org/mehari/issues/456))

### Features

* Make path-input-ped optional ([#472](https://github.com/varfish-org/mehari/issues/472)) ([e6cf526](https://github.com/varfish-org/mehari/commit/e6cf526dc5667d20c7a4bf1f76a698659b3dcb28))
* refactor DB create ([#456](https://github.com/varfish-org/mehari/issues/456)) ([c0fd47a](https://github.com/varfish-org/mehari/commit/c0fd47a4668a3052087e89fc890572b5217206b1))

## [0.25.7](https://github.com/varfish-org/mehari/compare/v0.25.6...v0.25.7) (2024-06-18)


### Bug Fixes

* make cargo release work again ([#469](https://github.com/varfish-org/mehari/issues/469)) ([a3f85fd](https://github.com/varfish-org/mehari/commit/a3f85fd23f0a1b2a73c684681b5c2bf73e805c3c))

## [0.25.6](https://github.com/varfish-org/mehari/compare/v0.25.5...v0.25.6) (2024-06-18)


### Bug Fixes

* bump annonars dependency ([#462](https://github.com/varfish-org/mehari/issues/462)) ([83b47d3](https://github.com/varfish-org/mehari/commit/83b47d31ae18dccf6efddf5a59e228ea19ae5e09))
* copy over keys::Var::from_vcf_allele from annonars to avoid noodles version conflicts in the future ([#461](https://github.com/varfish-org/mehari/issues/461)) ([952d8a9](https://github.com/varfish-org/mehari/commit/952d8a9c27fa2b00ab17f2dc5906b7528c2efed3))
* make previously public code public again ([#465](https://github.com/varfish-org/mehari/issues/465)) ([0ef7652](https://github.com/varfish-org/mehari/commit/0ef7652e6306f631bdf8d7e06f14173e7cd81cdb))
* Performance improvements ([#444](https://github.com/varfish-org/mehari/issues/444)) ([b996795](https://github.com/varfish-org/mehari/commit/b996795258ce5616a387f466963b131fa0ec5ac8))

## [0.25.5](https://github.com/varfish-org/mehari/compare/v0.25.4...v0.25.5) (2024-03-29)


### Bug Fixes

* actually produce one tsv record per gene+annotations ([#413](https://github.com/varfish-org/mehari/issues/413)) ([1ccabbb](https://github.com/varfish-org/mehari/commit/1ccabbbd2f59f04acd9eeb0faa62816728a7f40c))

## [0.25.4](https://github.com/varfish-org/mehari/compare/v0.25.3...v0.25.4) (2024-03-08)


### Bug Fixes

* stop codon issue for Ensembl chrMT transcripts ([#401](https://github.com/varfish-org/mehari/issues/401)) ([895e6f6](https://github.com/varfish-org/mehari/commit/895e6f6439b0b38e62ddfa83a10eb3bde5b998ef))

## [0.25.3](https://github.com/varfish-org/mehari/compare/v0.25.2...v0.25.3) (2024-03-07)


### Bug Fixes

* variants on chrX and chrY were erroneously ignored ([#399](https://github.com/varfish-org/mehari/issues/399)) ([0acc542](https://github.com/varfish-org/mehari/commit/0acc5422b1c54c992e6f4bbfabcfafeb990e2763))

## [0.25.2](https://github.com/varfish-org/mehari/compare/v0.25.1...v0.25.2) (2024-03-06)


### Bug Fixes

* strucvars/csq should not double-HGNC prefix gene IDs ([#397](https://github.com/varfish-org/mehari/issues/397)) ([b09782c](https://github.com/varfish-org/mehari/commit/b09782c2c22026a51edb112476990f9b95a67a55))

## [0.25.1](https://github.com/varfish-org/mehari/compare/v0.25.0...v0.25.1) (2024-03-06)


### Bug Fixes

* move assertion on chrMT transcripts ([#395](https://github.com/varfish-org/mehari/issues/395)) ([6a03acf](https://github.com/varfish-org/mehari/commit/6a03acf3ca309ab490557ad73652cc8b2a5d3c94))

## [0.25.0](https://github.com/varfish-org/mehari/compare/v0.24.2...v0.25.0) (2024-03-05)


### Features

* enable using ENSEMBL chrMT transcript ([#381](https://github.com/varfish-org/mehari/issues/381)) ([#393](https://github.com/varfish-org/mehari/issues/393)) ([91b4610](https://github.com/varfish-org/mehari/commit/91b46103e2a7e5b4ffdf1f884257b20408fc7c99))

## [0.24.2](https://github.com/varfish-org/mehari/compare/v0.24.1...v0.24.2) (2024-03-04)


### Features

* report MANE transcript counts in "db create" ([#378](https://github.com/varfish-org/mehari/issues/378)) ([2826b1a](https://github.com/varfish-org/mehari/commit/2826b1a6ff0660385bdd0f80488178a26d4d21aa))


### Bug Fixes

* make server genes/tx return empty page on missing gene ([#385](https://github.com/varfish-org/mehari/issues/385)) ([#387](https://github.com/varfish-org/mehari/issues/387)) ([f5d4a2f](https://github.com/varfish-org/mehari/commit/f5d4a2fb6aa45b026418330b1b9d834d4ef80132))
* miniscule fix to example URL output ([#382](https://github.com/varfish-org/mehari/issues/382)) ([dfb8824](https://github.com/varfish-org/mehari/commit/dfb882440e1ffacee6d2c4f108d39c6574ec641d))
* warn on missing transcripts rather than failing ([#380](https://github.com/varfish-org/mehari/issues/380)) ([#384](https://github.com/varfish-org/mehari/issues/384)) ([eaa99ae](https://github.com/varfish-org/mehari/commit/eaa99ae830c4e302e919f7e0c2acbf635cd7e5c0))


### Miscellaneous Chores

* release as v0.24.2 ([576dbb6](https://github.com/varfish-org/mehari/commit/576dbb6e1153dd2ea65753db3805e281c08431d5))

## [0.24.1](https://github.com/varfish-org/mehari/compare/v0.24.0...v0.24.1) (2024-02-26)


### Miscellaneous Chores

* cleanup test db bootstrapping ([#370](https://github.com/varfish-org/mehari/issues/370)) ([#371](https://github.com/varfish-org/mehari/issues/371)) ([9479cd0](https://github.com/varfish-org/mehari/commit/9479cd03ef55327395f4f47c2ba9f14252bc0765))

## [0.24.0](https://github.com/varfish-org/mehari/compare/v0.23.2...v0.24.0) (2024-02-26)


### Features

* adding "mehari db subset" command ([#367](https://github.com/varfish-org/mehari/issues/367)) ([#368](https://github.com/varfish-org/mehari/issues/368)) ([14044e3](https://github.com/varfish-org/mehari/commit/14044e33a45779f5bbdff0935e574e1deb926fa3))

## [0.23.2](https://github.com/varfish-org/mehari/compare/v0.23.1...v0.23.2) (2024-02-22)


### Bug Fixes

* handling missing FORMAT/AD values (".") ([#365](https://github.com/varfish-org/mehari/issues/365)) ([170ee00](https://github.com/varfish-org/mehari/commit/170ee003a8c8795fb6bc96c415ce57a703d45368))

## [0.23.1](https://github.com/varfish-org/mehari/compare/v0.23.0...v0.23.1) (2024-02-22)


### Bug Fixes

* work around issue with glnexus [#362](https://github.com/varfish-org/mehari/issues/362) ([#363](https://github.com/varfish-org/mehari/issues/363)) ([77afb96](https://github.com/varfish-org/mehari/commit/77afb96812944cda51cd2843ed836f83e71d2a86))

## [0.23.0](https://github.com/varfish-org/mehari/compare/v0.22.0...v0.23.0) (2024-02-21)


### Features

* adding support for Sniffles2 ([#357](https://github.com/varfish-org/mehari/issues/357)) ([#358](https://github.com/varfish-org/mehari/issues/358)) ([01641d7](https://github.com/varfish-org/mehari/commit/01641d7c5bfa7b8b2f65417bff062ec7d95c9f6e))
* script fix_glnexus.py prepare GLNexus output for noodles ([#356](https://github.com/varfish-org/mehari/issues/356)) ([#361](https://github.com/varfish-org/mehari/issues/361)) ([98d8233](https://github.com/varfish-org/mehari/commit/98d8233bfa2b744e3ad5d1a46a3191b0111df144))

## [0.22.0](https://github.com/varfish-org/mehari/compare/v0.21.2...v0.22.0) (2024-02-09)


### Features

* adding /genes/txs endpoint for server ([#339](https://github.com/varfish-org/mehari/issues/339)) ([#340](https://github.com/varfish-org/mehari/issues/340)) ([59dbaaf](https://github.com/varfish-org/mehari/commit/59dbaaf37f85e22c33515cb26d7ade1188c788fc))

## [0.21.2](https://github.com/varfish-org/mehari/compare/v0.21.1...v0.21.2) (2024-02-01)


### Miscellaneous Chores

* applying org change (bihealth =&gt; varfish-org) ([#325](https://github.com/varfish-org/mehari/issues/325)) ([8db8db0](https://github.com/varfish-org/mehari/commit/8db8db030b2076a491b901b444ac608972e259f9))

## [0.21.1](https://github.com/varfish-org/mehari/compare/v0.21.0...v0.21.1) (2023-12-01)


### Bug Fixes

* bumping annonars dependency to fix protobuf issues ([#287](https://github.com/varfish-org/mehari/issues/287)) ([007dd12](https://github.com/varfish-org/mehari/commit/007dd1244300f1d886b4b76a2e2d529c011d7b41))

## [0.21.0](https://github.com/varfish-org/mehari/compare/v0.20.0...v0.21.0) (2023-11-22)


### ⚠ BREAKING CHANGES

* adjusting enum values ([#279](https://github.com/varfish-org/mehari/issues/279))

### Bug Fixes

* adjusting enum values ([#279](https://github.com/varfish-org/mehari/issues/279)) ([6921b39](https://github.com/varfish-org/mehari/commit/6921b3922d228ac3a5fcbfc3f6126af2032bf424))

## [0.20.0](https://github.com/varfish-org/mehari/compare/v0.19.1...v0.20.0) (2023-11-21)


### ⚠ BREAKING CHANGES

* moving crate::pbs::mehari => crate::pbs ([#277](https://github.com/varfish-org/mehari/issues/277))

### Code Refactoring

* moving crate::pbs::mehari =&gt; crate::pbs ([#277](https://github.com/varfish-org/mehari/issues/277)) ([e44211e](https://github.com/varfish-org/mehari/commit/e44211ed51ab2172afe3901f61be1d7e8f06feb4))

## [0.19.1](https://github.com/varfish-org/mehari/compare/v0.19.0...v0.19.1) (2023-11-21)


### Bug Fixes

* copy in protobufs into docker builds ([#275](https://github.com/varfish-org/mehari/issues/275)) ([f4d4c20](https://github.com/varfish-org/mehari/commit/f4d4c203c8e29487e1b20d7dab1330a9d61f68e2))

## [0.19.0](https://github.com/varfish-org/mehari/compare/v0.18.1...v0.19.0) (2023-11-21)


### ⚠ BREAKING CHANGES

* serializing protobufs pbjson-build ([#272](https://github.com/varfish-org/mehari/issues/272)) (#273)

### Features

* serializing protobufs pbjson-build ([#272](https://github.com/varfish-org/mehari/issues/272)) ([#273](https://github.com/varfish-org/mehari/issues/273)) ([0f948f7](https://github.com/varfish-org/mehari/commit/0f948f7e1b1fd45ca92f1dd7010ca36d951f88cd))

## [0.18.1](https://github.com/varfish-org/mehari/compare/v0.18.0...v0.18.1) (2023-11-19)


### Miscellaneous Chores

* update dependencies ([#268](https://github.com/varfish-org/mehari/issues/268)) ([a3c379e](https://github.com/varfish-org/mehari/commit/a3c379ea80d6931aa77ddc93941b1c8a042ff120))

## [0.18.0](https://github.com/varfish-org/mehari/compare/v0.17.1...v0.18.0) (2023-11-19)


### Features

* properly handle selenoprotein import from cdot ([#224](https://github.com/varfish-org/mehari/issues/224)) ([#265](https://github.com/varfish-org/mehari/issues/265)) ([20137ad](https://github.com/varfish-org/mehari/commit/20137adff0d9f575885c3ee1f6691353a3d6d5b2))

## [0.17.1](https://github.com/varfish-org/mehari/compare/v0.17.0...v0.17.1) (2023-11-15)


### Bug Fixes

* fixing implementation of Provider::get_tx_for_gene ([#262](https://github.com/varfish-org/mehari/issues/262)) ([42945dd](https://github.com/varfish-org/mehari/commit/42945dd3d08d4921e7ea5ab00525b07a192ad3c8))

## [0.17.0](https://github.com/varfish-org/mehari/compare/v0.16.0...v0.17.0) (2023-11-10)


### ⚠ BREAKING CHANGES

* make --path-input-ped for seqvars annotation required ([#194](https://github.com/varfish-org/mehari/issues/194)) (#255)

### Features

* adding helper script to fix FreeBayes VCF ([#252](https://github.com/varfish-org/mehari/issues/252)) ([#254](https://github.com/varfish-org/mehari/issues/254)) ([4bd5461](https://github.com/varfish-org/mehari/commit/4bd5461bdbabb1124e17c735bad0c1b282fa2712))
* adding support for ClinCNV ([#253](https://github.com/varfish-org/mehari/issues/253)) ([#257](https://github.com/varfish-org/mehari/issues/257)) ([aba47c9](https://github.com/varfish-org/mehari/commit/aba47c93f21aa8752d762e72a4635df19215f643))


### Bug Fixes

* make --path-input-ped for seqvars annotation required ([#194](https://github.com/varfish-org/mehari/issues/194)) ([#255](https://github.com/varfish-org/mehari/issues/255)) ([de832b2](https://github.com/varfish-org/mehari/commit/de832b22ec6a55429e63c365c3bb787ea2e94209))

## [0.16.0](https://github.com/varfish-org/mehari/compare/v0.15.2...v0.16.0) (2023-11-09)


### Features

* allow user to select transcript source on CLI ([#247](https://github.com/varfish-org/mehari/issues/247)) ([#249](https://github.com/varfish-org/mehari/issues/249)) ([caab82e](https://github.com/varfish-org/mehari/commit/caab82ec2666ee5c38478cc2455adcb1d0a75250))
* update tx support to latest cdot with MANE label transfer ([#245](https://github.com/varfish-org/mehari/issues/245)) ([#248](https://github.com/varfish-org/mehari/issues/248)) ([b69e80d](https://github.com/varfish-org/mehari/commit/b69e80d35501beda2f29cbf635f9e4768751975f))

## [0.15.2](https://github.com/varfish-org/mehari/compare/v0.15.1...v0.15.2) (2023-10-28)


### Bug Fixes

* remove second async runtime creation ([#234](https://github.com/varfish-org/mehari/issues/234)) ([a6c3e06](https://github.com/varfish-org/mehari/commit/a6c3e063242d9fc6abbd4887d3902c03ad821776))

## [0.15.1](https://github.com/varfish-org/mehari/compare/v0.15.0...v0.15.1) (2023-10-24)


### Bug Fixes

* properly flush and sync all writers ([#232](https://github.com/varfish-org/mehari/issues/232)) ([2772565](https://github.com/varfish-org/mehari/commit/2772565b68a528809baf04ebe9e749cddbd2f723))

## [0.15.0](https://github.com/varfish-org/mehari/compare/v0.14.3...v0.15.0) (2023-10-23)


### Features

* adding async I/O, espec. for SV caller guessing ([#230](https://github.com/varfish-org/mehari/issues/230)) ([33c0e8e](https://github.com/varfish-org/mehari/commit/33c0e8efec783d3f3d504cd0619ab4ac234ce4d6))

## [0.14.3](https://github.com/varfish-org/mehari/compare/v0.14.2...v0.14.3) (2023-10-22)


### Bug Fixes

* index out of bounds error in stop_retained variant ([#226](https://github.com/varfish-org/mehari/issues/226)) ([#228](https://github.com/varfish-org/mehari/issues/228)) ([cc8d440](https://github.com/varfish-org/mehari/commit/cc8d4405f1c96c095fc7fa7d2e6ca08c3386563a))

## [0.14.2](https://github.com/varfish-org/mehari/compare/v0.14.1...v0.14.2) (2023-10-21)


### Bug Fixes

* make distance to next exon correct ([#222](https://github.com/varfish-org/mehari/issues/222)) ([#223](https://github.com/varfish-org/mehari/issues/223)) ([cad307f](https://github.com/varfish-org/mehari/commit/cad307fe9c2c150d254d471f508601dda01df363))
* switching from unmaintained tempdir to tempfile ([#227](https://github.com/varfish-org/mehari/issues/227)) ([a7b0e08](https://github.com/varfish-org/mehari/commit/a7b0e08f69b5a3527f78119bc0cf220ddca64597))

## [0.14.1](https://github.com/varfish-org/mehari/compare/v0.14.0...v0.14.1) (2023-10-19)


### Bug Fixes

* dependency bumps, including annonars ([#220](https://github.com/varfish-org/mehari/issues/220)) ([ac57c66](https://github.com/varfish-org/mehari/commit/ac57c6646dfdac25cce891edeee5f4920fbff5ee))

## [0.14.0](https://github.com/varfish-org/mehari/compare/v0.13.0...v0.14.0) (2023-10-18)


### Features

* bump annonars, write out VCV+RCV+clinsig ([#217](https://github.com/varfish-org/mehari/issues/217)) ([0d7b6d8](https://github.com/varfish-org/mehari/commit/0d7b6d8f7ca62b7152dad0d2da6557bb62580e58))

## [0.13.0](https://github.com/varfish-org/mehari/compare/v0.12.0...v0.13.0) (2023-10-16)


### Features

* derive Clone and strum::EnumIter on more enums ([#213](https://github.com/varfish-org/mehari/issues/213)) ([cc61a66](https://github.com/varfish-org/mehari/commit/cc61a66c175a4447c0422b01bbb6c572e1687725))
* implementing Default for AnnField ([#215](https://github.com/varfish-org/mehari/issues/215)) ([73b1be4](https://github.com/varfish-org/mehari/commit/73b1be40dd6f3c7c54aaeb2d0142444a702af1e3))
* proper serde renames for Consequence and PutativeImpact ([#214](https://github.com/varfish-org/mehari/issues/214)) ([02286fa](https://github.com/varfish-org/mehari/commit/02286fa15be47053d2e80a8f707dbb2361772033))


### Bug Fixes

* documenting "AN" as alleles rather than samples ([#206](https://github.com/varfish-org/mehari/issues/206)) ([62ee63b](https://github.com/varfish-org/mehari/commit/62ee63b84f9509cfaffad33ec9aa0ee2181681e0))

## [0.12.0](https://github.com/varfish-org/mehari/compare/v0.11.0...v0.12.0) (2023-10-09)


### Features

* expose strucvars annotation for worker aggregation command ([#204](https://github.com/varfish-org/mehari/issues/204)) ([0f730ce](https://github.com/varfish-org/mehari/commit/0f730ceefb356de37debe4348f5dcf29f30facc0))

## [0.11.0](https://github.com/varfish-org/mehari/compare/v0.10.0...v0.11.0) (2023-10-06)


### Features

* allow BGZ files for guess_sv_caller ([#203](https://github.com/varfish-org/mehari/issues/203)) ([31dc430](https://github.com/varfish-org/mehari/commit/31dc43037bf81e36f5581f09e743d29c86882869))


### Bug Fixes

* use new Debug trait implementation in Sample ([#201](https://github.com/varfish-org/mehari/issues/201)) ([7fafcbf](https://github.com/varfish-org/mehari/commit/7fafcbf224859c9aeba6bcd19968a9ce0f0e4d0c))

## [0.10.0](https://github.com/varfish-org/mehari/compare/v0.9.0...v0.10.0) (2023-10-05)


### Features

* expose clinvar annotation ([#199](https://github.com/varfish-org/mehari/issues/199)) ([3bcbf32](https://github.com/varfish-org/mehari/commit/3bcbf32c4cdeace14b03e7373b7950a478c39115))


### Bug Fixes

* Make mitochondrial DRAGEN calls not fail mehari ([#195](https://github.com/varfish-org/mehari/issues/195)) ([e170120](https://github.com/varfish-org/mehari/commit/e170120bcb14f3c79661f207dc4b734ad5f890e9))

## [0.9.0](https://github.com/varfish-org/mehari/compare/v0.8.0...v0.9.0) (2023-10-05)


### Features

* exposing more functions for sequence variant annotation ([#198](https://github.com/varfish-org/mehari/issues/198)) ([06a8dee](https://github.com/varfish-org/mehari/commit/06a8dee0ea8115395648434085fd9caf886a61ab))
* exposing utility functions for pedigree /  vcf building ([#196](https://github.com/varfish-org/mehari/issues/196)) ([7d77f35](https://github.com/varfish-org/mehari/commit/7d77f3575d72324de6591744411b777fcfc6e3a8))

## [0.8.0](https://github.com/varfish-org/mehari/compare/v0.7.0...v0.8.0) (2023-10-02)


### Features

* adding REST endpoint for SV annotation ([#188](https://github.com/varfish-org/mehari/issues/188)) ([9e14c09](https://github.com/varfish-org/mehari/commit/9e14c092c07017a83201a1bfc566ef3b203e15a8))
* snake_case in REST APIs and rename /tx/csq to /seqvars/csq ([#150](https://github.com/varfish-org/mehari/issues/150)) ([#189](https://github.com/varfish-org/mehari/issues/189)) ([e1021ed](https://github.com/varfish-org/mehari/commit/e1021eddef0c322fd1b29d6950cfff92e5932603))

## [0.7.0](https://github.com/varfish-org/mehari/compare/v0.6.2...v0.7.0) (2023-09-18)


### Features

* bump dependencies (especially noodles) ([#179](https://github.com/varfish-org/mehari/issues/179)) ([34b540c](https://github.com/varfish-org/mehari/commit/34b540c0a83a585e67550f9309fe25a28527fa5a))
* support for MELT VCF files ([#181](https://github.com/varfish-org/mehari/issues/181)) ([831e9af](https://github.com/varfish-org/mehari/commit/831e9af076a4b6334508ff9a6f11940f1ba42a6d))


### Bug Fixes

* Allow '.' genotypes ([#178](https://github.com/varfish-org/mehari/issues/178)) ([6d0c5b8](https://github.com/varfish-org/mehari/commit/6d0c5b8c16762efa78cdcc9577924c098fb3e9ec))

## [0.6.2](https://github.com/varfish-org/mehari/compare/v0.6.1...v0.6.2) (2023-07-17)


### Bug Fixes

* adjust clinvar column family name ([#138](https://github.com/varfish-org/mehari/issues/138)) ([95b7ffb](https://github.com/varfish-org/mehari/commit/95b7ffb6b42ef74679b8bb1433cdf0e43029db02))

## [0.6.1](https://github.com/varfish-org/mehari/compare/v0.6.0...v0.6.1) (2023-07-11)


### Bug Fixes

* strucvars annotate TSV output header ([#80](https://github.com/varfish-org/mehari/issues/80)) ([#130](https://github.com/varfish-org/mehari/issues/130)) ([f4d3315](https://github.com/varfish-org/mehari/commit/f4d3315a477496792216af77b0aa982f125111df))

## [0.6.0](https://github.com/varfish-org/mehari/compare/v0.5.7...v0.6.0) (2023-07-04)


### Features

* mehari becomes available as a library ([#118](https://github.com/varfish-org/mehari/issues/118)) ([59bb41a](https://github.com/varfish-org/mehari/commit/59bb41a7c4576e616a788b8788bab941fc6424ab))


### Bug Fixes

* fix Dockerfile by removing two-stage process ([#128](https://github.com/varfish-org/mehari/issues/128)) ([0346977](https://github.com/varfish-org/mehari/commit/03469774d466010a72cf15291effb4a1bb89f27c))
* wrong default of PATH_DB in Dockerfile ([#111](https://github.com/varfish-org/mehari/issues/111)) ([51e10cc](https://github.com/varfish-org/mehari/commit/51e10cc764d7a16384ef573d830697381c21472c))

## [0.5.7](https://github.com/varfish-org/mehari/compare/v0.5.6...v0.5.7) (2023-06-20)


### Bug Fixes

* adding libsqlite3-0 to the Docker image ([#109](https://github.com/varfish-org/mehari/issues/109)) ([be832fb](https://github.com/varfish-org/mehari/commit/be832fb91fe62e4941cc456a6b2d72d4a4dbf9cb))

## [0.5.6](https://github.com/varfish-org/mehari/compare/v0.5.5...v0.5.6) (2023-06-19)


### Bug Fixes

* dependency bump to make code compile again ([#104](https://github.com/varfish-org/mehari/issues/104)) ([8199d3e](https://github.com/varfish-org/mehari/commit/8199d3e9e46dc4eaa2d9e1ef2099d8c1f2f2de9b))
* docker build version in CI ([#107](https://github.com/varfish-org/mehari/issues/107)) ([9a0799f](https://github.com/varfish-org/mehari/commit/9a0799fb68915bef4080f0132da3d22963057e4d))


### Build System

* some small fixes to CI ([#106](https://github.com/varfish-org/mehari/issues/106)) ([fb54df0](https://github.com/varfish-org/mehari/commit/fb54df0f67223256b63ec6d3463ca47be1a5aaaf))

## [0.5.5](https://github.com/varfish-org/mehari/compare/v0.5.4...v0.5.5) (2023-06-19)


### Build System

* fix Docker builds ([#102](https://github.com/varfish-org/mehari/issues/102)) ([6366e65](https://github.com/varfish-org/mehari/commit/6366e654ca8afb879ec2fe2668da9d27593b562d))

## [0.5.4](https://github.com/varfish-org/mehari/compare/v0.5.3...v0.5.4) (2023-06-17)


### Build System

* adjust Docker builds for PRs and branches ([#100](https://github.com/varfish-org/mehari/issues/100)) ([e38388b](https://github.com/varfish-org/mehari/commit/e38388b3fcd509b99b4128944d4a289f07b6b4fb))

## [0.5.3](https://github.com/varfish-org/mehari/compare/v0.5.2...v0.5.3) (2023-06-16)


### Miscellaneous Chores

* re-release as v0.5.3 ([12aed0a](https://github.com/varfish-org/mehari/commit/12aed0aa98e8714c02425cf977faf33eb7d6c312))

## [0.5.2](https://github.com/varfish-org/mehari/compare/v0.5.1...v0.5.2) (2023-06-14)


### Bug Fixes

* docker build in CI ([#96](https://github.com/varfish-org/mehari/issues/96)) ([ba34c03](https://github.com/varfish-org/mehari/commit/ba34c034c17849b729dbc5faa9fd57bbaba913fa))

## [0.5.1](https://github.com/varfish-org/mehari/compare/v0.5.0...v0.5.1) (2023-06-13)


### Features

* build Docker image in CI ([#92](https://github.com/varfish-org/mehari/issues/92)) ([c72d804](https://github.com/varfish-org/mehari/commit/c72d80494043e46cf80113ccd984984d730e9368))


### Bug Fixes

* add Cargo.lock to git ([#94](https://github.com/varfish-org/mehari/issues/94)) ([f7c4f03](https://github.com/varfish-org/mehari/commit/f7c4f03628e8657a510f471245dfbc6d3c0382c8))
* remove unused dependencies ([#95](https://github.com/varfish-org/mehari/issues/95)) ([a1e1054](https://github.com/varfish-org/mehari/commit/a1e10542403339c4e7e421e74eef95e74ecf5468))

## [0.5.0](https://github.com/varfish-org/mehari/compare/v0.4.1...v0.5.0) (2023-06-12)


### Features

* REST server for variant consequence prediction ([#18](https://github.com/varfish-org/mehari/issues/18)) ([#90](https://github.com/varfish-org/mehari/issues/90)) ([b8c5a7b](https://github.com/varfish-org/mehari/commit/b8c5a7b6ec60e03fe7d973f96cab6ac2e628e7ea))

## [0.4.1](https://github.com/varfish-org/mehari/compare/v0.4.0...v0.4.1) (2023-06-12)


### Bug Fixes

* properly flush stream in "db create txs" ([#88](https://github.com/varfish-org/mehari/issues/88)) ([95458f3](https://github.com/varfish-org/mehari/commit/95458f3e0e8d598bfdb1ad24e8ae6bb930a725a9))

## [0.4.0](https://github.com/varfish-org/mehari/compare/v0.3.1...v0.4.0) (2023-06-08)


### Features

* move clinvar database code (import etc.) to annonars ([#86](https://github.com/varfish-org/mehari/issues/86)) ([457be46](https://github.com/varfish-org/mehari/commit/457be468591ee56d7b0c6afabde4e6ef18611d31))
* switch to annonars for frequency db creation ([#81](https://github.com/varfish-org/mehari/issues/81),[#77](https://github.com/varfish-org/mehari/issues/77),[#20](https://github.com/varfish-org/mehari/issues/20)) ([#83](https://github.com/varfish-org/mehari/issues/83)) ([bb27686](https://github.com/varfish-org/mehari/commit/bb2768686cbd96d20d064c58eaadbc18e09ce35c))


### Code Refactoring

* move from linked-hash-map to indexmap ([#82](https://github.com/varfish-org/mehari/issues/82)) ([#87](https://github.com/varfish-org/mehari/issues/87)) ([f9e354d](https://github.com/varfish-org/mehari/commit/f9e354dbcbe32e25c2a41d90ad1d0b2874332a06))

## [0.3.1](https://github.com/varfish-org/mehari/compare/v0.3.0...v0.3.1) (2023-05-04)


### Bug Fixes

* use copy number to fill in genotypes for CNVs ([#72](https://github.com/varfish-org/mehari/issues/72)) ([5585a08](https://github.com/varfish-org/mehari/commit/5585a084ea10de0612334530ff9a4a55b60518d3))

## [0.3.0](https://github.com/varfish-org/mehari/compare/v0.2.1...v0.3.0) (2023-05-03)


### Features

* add genotype counting ([#70](https://github.com/varfish-org/mehari/issues/70)) ([c40f83c](https://github.com/varfish-org/mehari/commit/c40f83ccfa5bd669e3acec2280084b072738695b))


### Bug Fixes

* only copy samples from input VCF header ([#69](https://github.com/varfish-org/mehari/issues/69)) ([316b123](https://github.com/varfish-org/mehari/commit/316b1237648435121e3d2c877426dc36ca31bcdd))

## [0.2.1](https://github.com/varfish-org/mehari/compare/v0.2.0...v0.2.1) (2023-05-02)


### Bug Fixes

* adjust non-adjacent exon recognition ([#67](https://github.com/varfish-org/mehari/issues/67)) ([96437e4](https://github.com/varfish-org/mehari/commit/96437e4a69349cabc7a165a130fd974a133fe307))
* allow multiple bases in BND ALT ([#63](https://github.com/varfish-org/mehari/issues/63)) ([e5b93d5](https://github.com/varfish-org/mehari/commit/e5b93d5727c16cd0881cfcc4e59777a6b9ab5923))
* properly recognize .gz/.zst extension ([#66](https://github.com/varfish-org/mehari/issues/66)) ([6960521](https://github.com/varfish-org/mehari/commit/696052175c1dbd5d9cc48026a4abddc33f81fd23))
* skipping "*" ALT alleles ([#65](https://github.com/varfish-org/mehari/issues/65)) ([95b2b61](https://github.com/varfish-org/mehari/commit/95b2b6170f1180545a7e3d327ee61d32ebf6ba59))

## [0.2.0](https://github.com/varfish-org/mehari/compare/v0.1.1...v0.2.0) (2023-04-24)


### Features

* allow annotation from maelstorm cov/mq BED files ([#51](https://github.com/varfish-org/mehari/issues/51)) ([#55](https://github.com/varfish-org/mehari/issues/55)) ([f3966ac](https://github.com/varfish-org/mehari/commit/f3966acb9a7e0d81be2f6ece51422766bafcca09))
* annotation of structural variants ([#46](https://github.com/varfish-org/mehari/issues/46)) ([#52](https://github.com/varfish-org/mehari/issues/52)) ([8da1d38](https://github.com/varfish-org/mehari/commit/8da1d388e21d6214e0cbe96ac54f703d9a736638))
* change database path setup ([#48](https://github.com/varfish-org/mehari/issues/48)) ([#60](https://github.com/varfish-org/mehari/issues/60)) ([b26671f](https://github.com/varfish-org/mehari/commit/b26671fb5c1b7ede8f0ede58c97d778586b761a0))
* implement command to copy rocksdb database ([#45](https://github.com/varfish-org/mehari/issues/45)) ([#47](https://github.com/varfish-org/mehari/issues/47)) ([3c8543f](https://github.com/varfish-org/mehari/commit/3c8543fbae281e8f9ea41bbad1718608fc9f00a6))
* switch from flatbuffers to protobuf ([#15](https://github.com/varfish-org/mehari/issues/15)) ([#57](https://github.com/varfish-org/mehari/issues/57)) ([3fe3322](https://github.com/varfish-org/mehari/commit/3fe332246fc5730a31d8a355bc92ba106e956db8))


### Bug Fixes

* issue on negative CDS position ([#49](https://github.com/varfish-org/mehari/issues/49)) ([#58](https://github.com/varfish-org/mehari/issues/58)) ([9a63e30](https://github.com/varfish-org/mehari/commit/9a63e3004e044c65c6552e3698b3a3daeeab3a6f))

## [0.1.1](https://github.com/varfish-org/mehari/compare/v0.1.0...v0.1.1) (2023-04-06)


### Bug Fixes

* release 0.1.1 ([5565d93](https://github.com/varfish-org/mehari/commit/5565d93e082de77f1f9b532bd581c827390e4b68))

## 0.1.0 (2023-04-06)


### Features

* adding command "db create seqvars-clinvar" ([#24](https://github.com/varfish-org/mehari/issues/24)) ([#25](https://github.com/varfish-org/mehari/issues/25)) ([fd4c523](https://github.com/varfish-org/mehari/commit/fd4c5238c2ff5b56577642ed4447b0abbf84d739))
* allow mehari to write out VarFish compatible TSV ([#34](https://github.com/varfish-org/mehari/issues/34)) ([#35](https://github.com/varfish-org/mehari/issues/35)) ([72c6832](https://github.com/varfish-org/mehari/commit/72c6832636c4a9a632439a2d0d9ebc2d21192209))
* command "verify seqvars" allows to compare to VEP ([#21](https://github.com/varfish-org/mehari/issues/21)) ([#37](https://github.com/varfish-org/mehari/issues/37)) ([b8761e3](https://github.com/varfish-org/mehari/commit/b8761e3f255b1ef9988e44d61dea5bd099bfaf44))
* frequency annotation with database ([#12](https://github.com/varfish-org/mehari/issues/12)) ([#13](https://github.com/varfish-org/mehari/issues/13)) ([f63f8c2](https://github.com/varfish-org/mehari/commit/f63f8c2f7e10487a126557bcbf8d7853ef418e7c))
* implement seqvar frequency db construction ([#2](https://github.com/varfish-org/mehari/issues/2)) ([#5](https://github.com/varfish-org/mehari/issues/5)) ([01fc64c](https://github.com/varfish-org/mehari/commit/01fc64c2f8aa351768bd5c1703a88d7a7cb021a8))
* implementing "annotate seqvars" ([#3](https://github.com/varfish-org/mehari/issues/3)) ([#7](https://github.com/varfish-org/mehari/issues/7)) ([6400ed1](https://github.com/varfish-org/mehari/commit/6400ed1d79ba736658549572ba0b16b6fe4626c7))
* transcript based annotation ([#11](https://github.com/varfish-org/mehari/issues/11)) ([#14](https://github.com/varfish-org/mehari/issues/14)) ([6677899](https://github.com/varfish-org/mehari/commit/66778992274107e9a632870ead0b6d97161f7b7e))
* transcript database building ([#1](https://github.com/varfish-org/mehari/issues/1)) ([#8](https://github.com/varfish-org/mehari/issues/8)) ([bcff954](https://github.com/varfish-org/mehari/commit/bcff9546b535fbd98caa35f89102cf6539be33b3))


### Bug Fixes

* allow loading JSON from .gz files ([#22](https://github.com/varfish-org/mehari/issues/22)) ([cefd950](https://github.com/varfish-org/mehari/commit/cefd950d823f40f94b095358a1ce7cdd7d834843))
* issues occuring when annotating real data ([#28](https://github.com/varfish-org/mehari/issues/28)) ([0636d49](https://github.com/varfish-org/mehari/commit/0636d490d2a17a0e8688c85d5b018ab1d0ae117f))

## Changelog
