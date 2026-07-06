---
name: karpathy-principle
description: Apply when writing, editing, reviewing, or refactoring ANY code in this repo — coding-behavior rules to avoid overcomplication and make surgical, verifiable changes.
---

# Karpathy principle — coding behavior

1. **Think before coding.** State assumptions explicitly; if uncertain, ask. If multiple interpretations exist, present them — don't silently pick. If a simpler approach exists, say so. If something is unclear, stop and name it.
2. **Simplicity first.** Minimum code that solves the problem, nothing speculative: no features beyond the ask, no abstractions for single-use code, no unrequested configurability, no error handling for impossible states. If 200 lines could be 50, rewrite. Ask "would a senior engineer call this overcomplicated?".
3. **Surgical changes.** Touch only what you must. Don't "improve" adjacent code/comments/formatting, don't refactor what isn't broken, match existing style. Remove only imports/vars/functions *your* change orphaned; mention pre-existing dead code, don't delete it. Every changed line traces to the request.
4. **Goal-driven (red-green).** Turn tasks into verifiable goals. Preferred: write the failing `cargo test` (or `cargo insta` snapshot) FIRST, confirm it fails, implement, confirm it passes — especially for bug reports (a reproducing test proves the bug, then proves the fix). For multi-step work, state a brief numbered plan with a `verify:` check per step.
5. **Work with the tool.** Follow cargo/clippy/crate-ecosystem conventions and official patterns, not workarounds. Building a complex workaround? Revisit the approach.

Applies to the docs you write too: concise, high signal, no restating what the code already shows.
