---
name: python-bindings
description: Apply when working in mehari-python/ — the PyO3/maturin Python extension.
---

# Python bindings (mehari-python/)

- **Separate crate, excluded from the workspace** (own `[workspace]` stub + `Cargo.lock`) so jemalloc isn't linked into the extension. Never add it to the root workspace `members`. It depends on `mehari` with `default-features = false`.
- Stack: **PyO3 0.28** (`extension-module`), **maturin** backend (`module-name = "mehari._mehari"`), **uv** for the env. Returns annotation results as **Arrow** RecordBatches (`arrow` 58 + `serde_arrow`); Python side uses polars/pyarrow. Requires Python ≥3.12. PyPI package name is `mehari`.
- Dev loop (from `mehari-python/`): `uv sync`, then `uv run pytest`. Build the extension with `uv run maturin develop`.
- The `sync-lockfiles` workflow keeps `mehari-python/Cargo.lock` aligned with root on release PRs — don't fight it.
- release-please owns version bumps here too (`pyproject.toml` + `Cargo.toml`) — don't hand-edit.
