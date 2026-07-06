---
name: cargo
description: Apply for build/test/lint, dependency or feature changes, and editing Cargo.toml or .cargo config in this repo.
---

# Cargo & manifest conventions

## Layout
- Workspace (`resolver=3`), one member: `mehari/`. `mehari-python/` is **excluded** (own `[workspace]` stub + `Cargo.lock`) — never add it to root `members`. See **python-bindings**.

## Commands
- Build `cargo build` · test `cargo test` / `cargo test --all-features` · format `cargo fmt` (check `cargo fmt -- --check`) · lint `cargo clippy --workspace --all-features` · coverage `cargo llvm-cov --all-features --workspace`.
- Run the binary: `cargo run -p mehari --bin mehari -- <args>`.
- Requires `protoc` + system `librocksdb-dev libsnappy-dev libsqlite3-dev`.

## Features
- `default = ["jemalloc", "server"]`. `server` gates actix/utoipa (+ `annonars/server`); `jemalloc` the allocator; `dhat-heap` heap profiling; `documentation` doc includes.

## Manifest (TOML) conventions
- Keep `[dependencies]` alphabetically ordered (match existing).
- Use `workspace = true` inheritance where the workspace defines a dep; otherwise follow existing style.
- `.cargo/config.toml` sets `ROCKSDB_LIB_DIR`/`SNAPPY_LIB_DIR=/usr/lib/` and a thumbv7em linker — don't remove/break these.
- release-please owns the `version` field — don't bump by hand. The `sync-lockfiles` workflow keeps root and `mehari-python/Cargo.lock` aligned on release PRs.

## Trap
- `mehari/src/db/transcripts/create/build.rs` is a **regular module named build.rs**, NOT a Cargo build script. The real build script is `mehari/build.rs`.
