---
name: openapi
description: Apply when changing the REST API surface — actix handlers, utoipa annotations, or request/response types under server/.
---

# OpenAPI / REST schema

- Server: actix-web + utoipa, behind the `server` feature. `ApiDoc` in `mehari/src/server/run/mod.rs`; Swagger UI via `utoipa-swagger-ui`.
- Spec is checked in at `openapi.schema.yaml` (repo root). A CI **Schema** job regenerates it and `diff`s against the committed file (stripping the `version:` line), so it must stay in sync.

## Hard rule
After ANY API-surface change (routes, handlers, `#[utoipa::path]`, schema structs), regenerate and commit the schema in the same PR, or CI fails:

```
cargo run -p mehari --bin mehari -- server schema --output-file openapi.schema.yaml
```

Review the diff (the `version:` line doesn't matter — CI strips it).
