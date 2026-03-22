// The custom build script, used to (1) generate the Rust classes for the
// protobuf implementation and (2) use pbjson for proto3 JSON serialization.

use std::{env, path::PathBuf};

fn main() -> Result<(), anyhow::Error> {
    // Integration of `prost-build` and `pbjson-build`.
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("protos");
    let proto_files = ["mehari/txs.proto", "mehari/server.proto"]
        .iter()
        .map(|f| root.join(f))
        .collect::<Vec<_>>();

    // Tell cargo to recompile if any of these proto files are changed
    for proto_file in &proto_files {
        println!("cargo:rerun-if-changed={}", proto_file.display());
    }

    let descriptor_path: PathBuf =
        PathBuf::from(env::var("OUT_DIR").unwrap()).join("proto_descriptor.bin");

    prost_build::Config::new()
        // Save descriptors to file
        .file_descriptor_set_path(&descriptor_path)
        // Override prost-types with pbjson-types
        .compile_well_known_types()
        .extern_path(".google.protobuf", "::pbjson_types")
        // Define the protobuf files to compile.
        .compile_protos(&proto_files, &[root])?;

    let descriptor_set = std::fs::read(descriptor_path).unwrap();
    pbjson_build::Builder::new()
        .register_descriptors(&descriptor_set)?
        .build(&[".mehari"])?;

    // Integration of `built`. Workaround for python bindings via maturin in mehari-python
    let src = env::var("CARGO_MANIFEST_DIR").unwrap();
    let dst = PathBuf::from(env::var("OUT_DIR").unwrap()).join("built.rs");

    let mut built_opts = built::Options::default();
    built_opts.set_dependencies(true);

    if let Err(e) = built::write_built_file_with_options(&built_opts, src.as_ref(), &dst) {
        println!(
            "cargo:warning=Failed to write built file with dependencies: {}",
            e
        );
        println!(
            "cargo:warning=Falling back to building without dependencies (likely Python sdist)"
        );

        built_opts.set_dependencies(false);
        built::write_built_file_with_options(&built_opts, src.as_ref(), &dst)
            .map_err(|e| anyhow::anyhow!("Failed to write built file: {}", e))?;
    }

    Ok(())
}
