// The custom build script, needed as we use prost.

fn main() {
    println!("cargo:rerun-if-changed=src/proto/data.proto3");
    prost_build::Config::new()
        .protoc_arg("-Isrc/proto")
        // Add serde serialization and deserialization to the generated code.
        .type_attribute(".", "#[derive(serde::Serialize, serde::Deserialize)]")
        // Skip serializing `None` values.
        .type_attribute(".", "#[serde_with::skip_serializing_none]")
        // Define the protobuf files to compile.
        .compile_protos(&["src/proto/data.proto3"], &["src/"])
        .unwrap();
}
