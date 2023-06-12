// The custom build script, needed as we use prost.

fn main() {
    prost_build::compile_protos(&["src/db/create/txs/data.proto3"], &["src/"]).unwrap();
}
