# NMD plugin for mehari, in Python

mehari's built-in `nmd` plugin, written in Python against mehari's plugin interface (`mehari/wit/plugin.wit`).
It is an example for writing plugins in Python; see `docs/plugins.md` for the guide.

Build it with [componentize-py](https://github.com/bytecodealliance/componentize-py) and run it:

```sh
uvx componentize-py -d ../../mehari/wit -w plugin componentize app -o nmd.wasm
mehari annotate seqvars --plugin nmd.wasm ...
```

For editor support, generate the Python bindings first:

```sh
uvx componentize-py -d ../../mehari/wit -w plugin bindings .
```

The component bundles a Python interpreter (about 18 MB), so it is not committed.
After building it, check that it predicts the same as the built-in plugin:

```sh
cargo test -p mehari -- --ignored nmd_python_plugin_matches_builtin
```
