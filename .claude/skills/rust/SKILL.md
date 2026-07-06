---
name: rust
description: Apply when writing or editing Rust code in this crate — idioms, error handling, logging, and lint rules.
---

# Rust conventions

- **Edition 2024, MSRV 1.91.** Don't use APIs newer than the MSRV.
- **Errors:** `anyhow::Result` for application/CLI paths (add `.context(...)`); `thiserror` enums for library error types (`src/errors.rs`).
- **No panics in library code:** no `unwrap`/`expect`/`panic!`/`todo!` on reachable paths — return `Result`. Tests may `unwrap`.
- **Logging:** `tracing` spans/events, not `println!`. Verbosity via `clap-verbosity-flag`.
- **Parallelism:** `rayon` iterators for data-parallel work.
- **Newtypes:** use `nutype` where a wrapper adds real invariants/safety; don't wrap for its own sake.
- **Lints:** keep clippy clean. CI runs `cargo clippy --workspace --all-features` (no `-D warnings` here, but fix warnings anyway). Run `cargo fmt` before finishing.
- Prefer crates already in the tree (itertools, indexmap, rustc-hash, strum, enumflags2) over hand-rolling.
