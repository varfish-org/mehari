Mehari end-user documentation.

## What is Mehari?

Mehari is a software for annotating genetic variants in VCF format similar to the [Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html) and [SnpEff](http://pcingola.github.io/SnpEff/).

Why another software package?

- Mehari uses the [hgvs](https://crates.io/crates/hgvs) library which produces the same transcript and protein level predictions as the [biocommons/hgvs](github.com/biocommons/hgvs)
  library.
  The latter serves as the basis for [VariantValidator.org](https://variantvalidator.org/) which is the gold standard for HGVS variant description generation and validation.
- Mehari is written in the Rust programming language which allows it to work fast, with low memory consumption (as a C++ program would) and being memory safe at the same time (as a Java/Python/Perl program would).
- As a Rust program, it can be embedded into the backend of the [VarFish](https://github.com/varfish-org/varfish-server) variant analysis platform.

## What's Next?

We recommend to read the Mehari end-user documentation in the following order:

- [Overview](`self::user_doc::getting_started`)
- [Annotating Sequence Variants](`self::user_doc::anno_seqvars`)
- [Annotating Structural Variants](`self::user_doc::anno_strucvars`)
- [Building Databases](`self::user_doc::db_build`)
- [Implementation Notes](`self::user_doc::implementation_notes`)

## Documentation Formatting / Hosting Notes

Since Mehari is written in the Rust programming language, we host the documentation on `docs.rs` written as Rust online documentation.
This has the advantage that the documentation is bundle with the program source code (and thus always up to date) and the latest documentation is always available at <https://docs.rs/mehari>.

The drawback is that the formatting of this may not be as end-user friendly as it could be but you will manage.
