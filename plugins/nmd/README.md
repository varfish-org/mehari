# NMD plugin for mehari

The reference WebAssembly plugin for mehari's plugin interface (`mehari/wit/plugin.wit`): mehari's built-in `nmd` plugin, ported to the interface.
mehari's tests check that both predict the same.

Build it and run it:

```sh
rustup target add wasm32-wasip2
cargo build --release --target wasm32-wasip2
mehari annotate seqvars --plugin target/wasm32-wasip2/release/mehari_plugin_nmd.wasm ...
```

After changing the plugin or the interface, update the copy the tests use:

```sh
cp target/wasm32-wasip2/release/mehari_plugin_nmd.wasm ../../mehari/tests/data/plugins/nmd.wasm
```
