## Building from scratch
By default, the RocksDB bundled with the `rocksdb` crate (see `librocksdb-sys` in `Cargo.lock`) is built from source.
This takes a while on the first build but guarantees that every build, including the Docker image and the bioconda package, reads and writes the same database format.

To reduce compile times locally, you can link a pre-built `rocksdb` instead by setting the environment variables `ROCKSDB_LIB_DIR` and `SNAPPY_LIB_DIR`:

```bash
export ROCKSDB_LIB_DIR=/usr/lib/
export SNAPPY_LIB_DIR=/usr/lib/
```

Note that the system package may be too old: Ubuntu 24.04 ships RocksDB 8.9, which cannot read databases written by the bundled RocksDB (table `format_version` 7).

You will need a recent version of protoc, e.g.:

```bash
bash utils/install-protoc.sh
export PATH=$PATH:$HOME/.local/share/protoc/bin
```

To build the project, run:
```bash
cargo build --release
```

To install the project locally, run:
```bash
cargo install --path .
```