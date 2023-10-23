// The custom build script, needed as we use prost.

fn main() {
    println!("cargo:rerun-if-changed=src/db/create/txs/data.proto3");
    prost_build::compile_protos(&["src/db/create/txs/data.proto3"], &["src/"]).unwrap();
}
