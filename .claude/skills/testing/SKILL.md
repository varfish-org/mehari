---
name: testing
description: Apply when writing or running tests, or updating insta snapshots, in this repo.
---

# Testing

- Run: `cargo test` (CI default) or `cargo test --all-features`. Prefer red-green (see **karpathy-principle**).
- **insta** snapshots (YAML): after intended output changes, review with `cargo insta review` (or `cargo insta accept`). Never blindly accept.
- **Selective Git LFS:** only large fixtures/specific snapshots are LFS-tracked (`tests/data/**`, `*.bin`, and named snaps in `.gitattributes`); most `.snap` files are plain git. Have `git lfs` installed so fixtures resolve.
- `rstest` for parameterized tests/fixtures; `tracing-test` to assert on logs; `temp_testdir`/`tempfile` for scratch dirs; `pretty_assertions` for diffs.
- The full UTA-backed suite needs a Postgres service (see CI); most tests don't.
