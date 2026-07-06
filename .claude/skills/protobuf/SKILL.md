---
name: protobuf
description: Apply when editing .proto files or the generated protobuf types (src/pbs/*) in this repo.
---

# Protobuf

- Schemas: `mehari/protos/mehari/*.proto` (`txs.proto`, `seqvars.proto`, `server.proto`), proto3.
- Codegen at build time via `mehari/build.rs`: `prost-build` + `pbjson-build` (proto3 JSON). Well-known types map to `pbjson_types` (`extern_path .google.protobuf`); a file descriptor set is written to `OUT_DIR`.
- Generated Rust is surfaced through `src/pbs/{txs,server,seqvars}.rs`. Don't edit generated code — edit the `.proto` and rebuild (`cargo build`, needs `protoc`).

## Safe-edit rules (wire compatibility)
- **Never renumber or reuse an existing field number.** Add new fields with new numbers.
- Add fields; don't mutate/repurpose existing ones. To retire a field, `reserved` its number/name.
- Changing a field's type or label is breaking — add a new field instead.
- `txs.proto` is the on-disk transcript DB format; changing it can invalidate already-built databases — treat as a data migration.
