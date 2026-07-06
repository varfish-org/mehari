# mehari — agent guide

Variant effect prediction (HGVS, VEP-compatible consequences) for VCF, in Rust. CLI + Rust library + REST server + Python bindings. Org `varfish-org`, license MIT.

## Layout
- Cargo **workspace** (`resolver=3`), single member `mehari/` (lib + bin). `mehari-python/` is a **separate** crate, deliberately **excluded** from the workspace (own `[workspace]` stub + `Cargo.lock`) so jemalloc isn't linked into the extension.
- `mehari/src/`: `annotate/` (seqvars, strucvars), `db/` (transcript/CADD/spliceai/dbsnp builders; keys in `db/keys.rs`), `pbs/` (generated protobuf), `server/` (actix + utoipa REST).
- Depends on the sibling crate `annonars` (ClinVar/frequency DBs).
- CLI subcommands: `db` / `annotate` / `server`.

## Build / test / lint (real commands)
- Prereqs: `protoc`; system `librocksdb-dev libsnappy-dev libsqlite3-dev`. `.cargo/config.toml` points RocksDB/snappy at `/usr/lib/` — don't break it.
- Build `cargo build` · test `cargo test` (CI) / `cargo test --all-features` · format `cargo fmt -- --check` · lint `cargo clippy --workspace --all-features` · coverage `cargo llvm-cov --all-features --workspace`.
- Run: `cargo run -p mehari --bin mehari -- <args>`.
- Python bindings: from `mehari-python/`, `uv sync` then `uv run pytest`.

## House conventions
- `anyhow` for app/CLI paths, `thiserror` for library error types. `tracing` (not `println!`). `rayon` for data parallelism. No `unwrap`/`expect`/`panic!` in library code.
- Edition 2024, MSRV 1.91. Keep clippy clean (CI here doesn't pass `-D warnings`, but treat warnings as errors anyway).
- release-please owns `CHANGELOG.md` and version fields — never hand-edit.

## Skills (`.claude/skills/`) — load the matching one before the task
- **karpathy-principle** — always, when writing/editing code.
- **rust** — Rust idioms & error handling.
- **cargo** — build/test/deps/features/manifests.
- **protobuf** — editing `.proto` or generated types.
- **openapi** — changing the REST API surface.
- **contributing** — branching, commits, `gh`, PRs.
- **issue-workflow** — creating/tracking issues, multi-step plans.
- **testing** — tests & insta snapshots.
- **python-bindings** — anything under `mehari-python/`.
