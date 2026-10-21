Plugins add your own columns to the `ANN` field of `mehari annotate seqvars`.

A plugin is a WebAssembly component: write it in Python, Rust or another language with WebAssembly component tooling, build it into a `.wasm` file, and pass that to `--plugin`.
Plugins run sandboxed and get the annotation of every transcript plus, on request, its exon structure and its transcript and protein sequences with and without the variant.

```sh
mehari annotate seqvars --plugin my_plugin.wasm ...
```

`--plugin` also takes the names of the built-in plugins (see [Annotating Sequence Variants](`self::anno_seqvars`)) and may be given several times; the columns appear in that order.

## How a plugin works

A plugin implements the `plugin` world of [`mehari/wit/plugin.wit`](https://github.com/varfish-org/mehari/blob/main/mehari/wit/plugin.wit), which defines two functions:

- `columns()` returns the names of the columns the plugin fills. mehari calls it once, when loading the plugin.
- `annotate(variants, annotations, contexts)` is called once per variant with the annotations of all its reported transcripts. It returns one list per annotation, holding one value per column; an empty value (`None` in Python) leaves the field empty.

The arguments:

- `variants` holds the variant as written in the VCF, before normalization, or the variants of a compound group (`--enable-compound-variants`).
- Each annotation has the transcript accession (`feature-id`, empty for intergenic variants), the Sequence Ontology terms (`consequences`), the HGVS descriptions (`hgvs-g`, `hgvs-n`, `hgvs-c`, `hgvs-p`) and the transcript position (`cdna-pos`). These are the values mehari writes to `ANN`: HGVS normalization shifts the variant towards the transcript's 3' end, so its position depends on the transcript's strand.
- `contexts` gives more data for the annotation at a given index, computed when the plugin first asks for it and shared with the other plugins:
  - `transcript(index)`: strand and exons (transcript and genome positions, 1-based), or none for intergenic variants;
  - `reference(index)`: transcript sequence, protein sequence and the position of the start codon;
  - `alternative(index)`: the same with the variant applied, for coding variants only.

## Writing a plugin in Python

[componentize-py](https://github.com/bytecodealliance/componentize-py) builds a Python module, together with a Python interpreter, into a component.
The commands below run it with [uv](https://docs.astral.sh/uv/); `pip install componentize-py` works as well.
They were tested with componentize-py 0.25.1.

1. Create a directory for the plugin and copy `plugin.wit` from your mehari version into its `wit` subdirectory.
2. Optionally, generate the Python bindings so that your editor knows the types:

   ```sh
   uvx componentize-py -d wit -w plugin bindings .
   ```

3. Write `app.py`. The class must be called `WitWorld`. This example adds a `lof` column for loss-of-function consequences and an `exon_count` column:

   ```python
   import wit_world
   from wit_world.imports.types import Annotation, Contexts, VcfVariant

   LOSS_OF_FUNCTION = {
       "stop_gained",
       "frameshift_variant",
       "splice_acceptor_variant",
       "splice_donor_variant",
       "start_lost",
   }


   class WitWorld(wit_world.WitWorld):
       def columns(self) -> list[str]:
           return ["lof", "exon_count"]

       def annotate(
           self,
           variants: list[VcfVariant],
           annotations: list[Annotation],
           contexts: Contexts,
       ) -> list[list[str | None]]:
           rows = []
           for index, ann in enumerate(annotations):
               lof = "yes" if LOSS_OF_FUNCTION & set(ann.consequences) else None
               # Computed on request; None for intergenic variants.
               transcript = contexts.transcript(index)
               exon_count = str(len(transcript.exons)) if transcript else None
               rows.append([lof, exon_count])
           return rows
   ```

4. Build the component and use it:

   ```sh
   uvx componentize-py -d wit -w plugin componentize app -o lof.wasm
   mehari annotate seqvars --plugin lof.wasm ...
   ```

`plugins/nmd-python` in the mehari repository is a larger example: the built-in `nmd` plugin, written in Python.

Things to know about Python plugins:

- componentize-py bundles the modules that `app.py` imports when it is loaded, so import everything at the top of the file. A module imported inside a function is missing at run time.
- Pure-Python packages can be bundled: install them into a virtual environment that is active while building, or into a directory (`pip install --target deps ...`) that you add with `-p . -p deps`. `-p` replaces the default `.`, so keep it in the list. Packages with compiled extensions need builds for WebAssembly.
- `print` output to stdout is discarded; print to `sys.stderr` instead. An exception stops mehari and shows the Python traceback.
- The component contains the Python interpreter (about 18 MB). mehari compiles it when starting, which takes a few seconds on one thread and less with more `--threads`. Each thread then holds its own instance, about 25 MB each.

## Writing a plugin in Rust

Create a library crate of type `cdylib` that depends on [wit-bindgen](https://crates.io/crates/wit-bindgen), and build it for the `wasm32-wasip2` target.

```toml
[lib]
crate-type = ["cdylib"]

[dependencies]
wit-bindgen = "0.63"
```

```rust,ignore
wit_bindgen::generate!({ path: "wit", world: "plugin" });

struct MyPlugin;

impl Guest for MyPlugin {
    fn columns() -> Vec<String> {
        vec!["exon_count".into()]
    }

    fn annotate(
        _variants: Vec<VcfVariant>,
        annotations: Vec<Annotation>,
        contexts: &Contexts,
    ) -> Vec<Vec<Option<String>>> {
        (0..)
            .zip(&annotations)
            .map(|(index, _)| vec![contexts.transcript(index).map(|tx| tx.exons.len().to_string())])
            .collect()
    }
}

export!(MyPlugin);
```

```sh
rustup target add wasm32-wasip2
cargo build --release --target wasm32-wasip2
mehari annotate seqvars --plugin target/wasm32-wasip2/release/my_plugin.wasm ...
```

`plugins/nmd` in the mehari repository is the built-in `nmd` plugin, written this way.
Rust plugins are small (under 100 KB) and load in milliseconds.

## Other languages

Any language that builds WebAssembly components for WASI 0.2 should work, for example JavaScript with [jco](https://github.com/bytecodealliance/jco), Go with [TinyGo](https://tinygo.org) or C with [wit-bindgen](https://github.com/bytecodealliance/wit-bindgen).
These have not been tested with mehari.

## Sandbox, errors and versions

- Plugins have no access to files, the network or environment variables. Their stderr is shown, their stdout is discarded.
- If a plugin fails, or returns a different number of rows or values than expected, mehari stops with an error naming the plugin.
- Two plugins cannot fill the same column.
- The interface has a version (`package mehari:plugin@0.1.0` in `plugin.wit`). A component built for another version fails to load; rebuild it with the `plugin.wit` of your mehari version.
