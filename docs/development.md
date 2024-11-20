## Building from scratch
To reduce compile times, we recommend using a pre-built version of `rocksdb`, either from the system package manager or e.g. via `conda`:

```bash
# Ubuntu
sudo apt-get install librocksdb-dev

# Conda
conda install -c conda-forge rocksdb
```

In either case, either add
```toml
[env]
ROCKSDB_LIB_DIR = "/usr/lib/" # in case of the system package manager, adjust the path accordingly for conda
SNAPPY_LIB_DIR = "/usr/lib/"  # same as above
```
to `.cargo/config.toml` or set the environment variables `ROCKSDB_LIB_DIR` and `SNAPPY_LIB_DIR` to the appropriate paths:

```bash
export ROCKSDB_LIB_DIR=/usr/lib/
export SNAPPY_LIB_DIR=/usr/lib/
```

By default, the environment variables are defined in the `.cargo/config.toml` as described above, i.e. may need adjustments if not using the system package manager.

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