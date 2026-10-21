Plugins add your own columns to the `ANN` field of `mehari annotate seqvars`.

A plugin is a WebAssembly component: write it in Python, Rust or another language with WebAssembly component tooling, build it into a `.wasm` file, and pass that to `--plugin`, followed by its options:

```sh
mehari annotate seqvars --plugin my_plugin.wasm,data=/path/to/tables,threshold=0.5 ...
```

Plugins run sandboxed and see every transcript annotation of a variant, with its exon structure, its transcript and protein sequences and everything else mehari knows about it.
`--plugin` also takes the names of the built-in plugins (see [Annotating Sequence Variants](`self::anno_seqvars`)) and may be given several times; the columns appear in that order.

## How a plugin works

A plugin implements the `plugin` world of [`mehari/wit/plugin.wit`](https://github.com/varfish-org/mehari/blob/main/mehari/wit/plugin.wit), which defines two functions:

- `init(options)` gets the options from the command line as key-value pairs and returns the names of the columns the plugin fills, or an error. mehari runs a separate instance of the plugin per thread and calls `init` once in each.
- `annotate(variants, annotations)` is called once per variant and returns one list per annotation, holding one value per column, or an error. An empty value (`None` in Python) leaves the field empty.

`variants` holds the variant as written in the VCF, before normalization, or the variants of a compound group (`--enable-compound-variants`).
`annotations` has one entry per reported transcript; `annotations.len()` says how many.
For each index, it offers:

- what mehari writes to `ANN`: `feature-id` (transcript accession, empty for intergenic variants), `consequences` (Sequence Ontology terms), `impact`, `gene-symbol`, `gene-id`, `biotypes`, `tags`, `rank` (exon or intron), `hgvs-g`, `hgvs-n`, `hgvs-c`, `hgvs-p`, `cdna-pos`, `cds-pos`, `protein-pos`, `distance` and `warnings`. HGVS normalization shifts the variant towards the transcript's 3' end, so its position depends on the transcript's strand;
- `field(index, column)`: what an earlier plugin wrote into a column;
- `transcript`: strand, exons with transcript and genome positions, contig, genome build, protein accession and the genomic range of the coding sequence;
- `reference` and `alternative`: transcript sequence, protein sequence and start codon, without and with the variant; `cds-end`, `translation` (translation table and exceptions such as selenocysteine), `frameshift` and `first-changed-codon`.

mehari computes the transcript data when a plugin first asks for it and shares it with the other plugins, so ask only for what you need.
The `host` interface tells the genome assembly, the mehari version and whether `--vep-compatibility-mode` is on.
`plugin.wit` documents every function.

## Options and data files

Options follow the plugin's path as `,key=value`; the plugin decides what they mean.
One option is mehari's: `data=<directory>` lets the plugin read that directory as `/data`, e.g. to look up a table there.
Plugins cannot write it, nor read anything else.
`init` runs once per thread, so load small tables there, and look up large ones on demand.

## Writing a plugin in Python

[componentize-py](https://github.com/bytecodealliance/componentize-py) builds a Python module, together with a Python interpreter, into a component.
The commands below run it with [uv](https://docs.astral.sh/uv/); `pip install componentize-py` works as well.
They were tested with componentize-py 0.25.1.

1. Create a directory for the plugin and copy `plugin.wit` from your mehari version into its `wit` subdirectory.
2. Optionally, generate the Python bindings so that your editor knows the types:

   ```sh
   uvx componentize-py -d wit -w plugin bindings .
   ```

3. Write `app.py`. The class must be called `WitWorld`; raise `Err` to report an error. This example adds a `lof` column for loss-of-function consequences and an `in_panel` column for genes listed in a file:

   ```python
   import wit_world
   from componentize_py_types import Err
   from wit_world.imports.types import Annotations, VcfVariant

   LOSS_OF_FUNCTION = {
       "stop_gained",
       "frameshift_variant",
       "splice_acceptor_variant",
       "splice_donor_variant",
       "start_lost",
   }
   PANEL: set[str] = set()


   class WitWorld(wit_world.WitWorld):
       def init(self, options: list[tuple[str, str]]) -> list[str]:
           for key, value in options:
               if key != "panel":
                   raise Err(f"unknown option {key!r}")
               # One gene symbol per line, in the directory given with `data`.
               with open(f"/data/{value}") as panel:
                   PANEL.update(line.strip() for line in panel)
           return ["lof", "in_panel"]

       def annotate(
           self, variants: list[VcfVariant], annotations: Annotations
       ) -> list[list[str | None]]:
           rows = []
           for index in range(annotations.len()):
               consequences = set(annotations.consequences(index))
               lof = "yes" if LOSS_OF_FUNCTION & consequences else None
               in_panel = "yes" if annotations.gene_symbol(index) in PANEL else None
               rows.append([lof, in_panel])
           return rows
   ```

4. Build the component and use it:

   ```sh
   uvx componentize-py -d wit -w plugin componentize app -o lof.wasm
   mehari annotate seqvars --plugin lof.wasm,data=panels,panel=cardio.txt ...
   ```

`plugins/nmd-python` in the mehari repository is a larger example: the built-in `nmd` plugin, written in Python.

Things to know about Python plugins:

- componentize-py bundles the modules that `app.py` imports when it is loaded, so import everything at the top of the file. A module imported inside a function is missing at run time.
- Pure-Python packages can be bundled: install them into a virtual environment that is active while building, or into a directory (`pip install --target deps ...`) that you add with `-p . -p deps`. `-p` replaces the default `.`, so keep it in the list. Packages with compiled extensions need builds for WebAssembly.
- `print` output to stdout is discarded; print to `sys.stderr` instead. An exception other than `Err` stops mehari and shows the Python traceback.
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
    fn init(options: Vec<(String, String)>) -> Result<Vec<String>, String> {
        match options.first() {
            Some((key, _)) => Err(format!("unknown option {key:?}")),
            None => Ok(vec!["exon_count".into()]),
        }
    }

    fn annotate(
        _variants: Vec<VcfVariant>,
        annotations: &Annotations,
    ) -> Result<Vec<Vec<Option<String>>>, String> {
        Ok((0..annotations.len())
            .map(|index| {
                let transcript = annotations.transcript(index);
                vec![transcript.map(|tx| tx.exons.len().to_string())]
            })
            .collect())
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

- Plugins have no access to the network or environment variables, and no file access beyond `data`. Their stderr is shown, their stdout is discarded.
- If a plugin fails, asks for an annotation past the last one, or returns a different number of rows or values than expected, mehari stops with an error naming the plugin.
- Two plugins cannot fill the same column.
- The interface has a version (`package mehari:plugin@0.1.0` in `plugin.wit`). Later 0.1.x versions only add functions, so plugins built for 0.1.0 keep working. A plugin built for another version fails to load; rebuild it with the `plugin.wit` of your mehari version.
