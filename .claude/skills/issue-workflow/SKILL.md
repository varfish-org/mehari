---
name: issue-workflow
description: Apply when creating or tracking a GitHub issue, or planning multi-step work — the issue template and the task convergence loop.
---

# Issue workflow

Create issues with `gh issue create`. The template also lives at `.github/ISSUE_TEMPLATE/task.md`. Sections marked *(if applicable)* are omitted for simple fixes.

```markdown
## 1. Goal
[1-3 sentences: what must be achieved and why.]
**Scope**: [In scope. Explicitly deferred.]

## 2. Implementation Strategy (if applicable — complex features)
[Approach, design decisions, key changes.]

## 3. Implementation Plan / Tasks
- [ ] [Concrete, independently completable task]

## 4. Acceptance Criteria
- [ ] [Observable, verifiable condition]
- [ ] Tests pass (`cargo test --all-features`)
- [ ] `cargo fmt --check` and `cargo clippy --all-features -- -D warnings` clean
- [ ] OpenAPI schema regenerated + committed if the API surface changed

## 5. Design Rationale (if applicable)
[Why this approach over the alternatives.]
```

**Rules:** granular `[ ]` tasks; verifiable criteria that always include the cargo test/fmt/clippy baseline; omit inapplicable sections.

## Convergence loop (optional)
Keep a gitignored `issue-<n>.md` at the repo root as the working source of truth (`issue-*.md` is git-ignored). Then loop, ONE task at a time:
1. Mark the task ` ← in progress`.
2. Implement red-green (see **karpathy-principle**).
3. Run `cargo clippy` / `cargo fmt` after each change.
4. Mark `[x]`; move to the next.

When all `[x]`: confirm the full suite green, then tick §4.
