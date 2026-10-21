# Probe plugin for mehari's tests

A WebAssembly plugin that reports what every function of mehari's plugin interface (`mehari/wit/plugin.wit`) returns, and that can fail on request; see the options in `src/lib.rs`.
mehari's tests in `plugin/wasm.rs` use it to check the host.

After changing the plugin or the interface, rebuild it and update the copy the tests use:

```sh
cargo build --release --target wasm32-wasip2
cp target/wasm32-wasip2/release/mehari_plugin_probe.wasm ../../mehari/tests/data/plugins/probe.wasm
```
