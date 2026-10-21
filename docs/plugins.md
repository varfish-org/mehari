Plugins add your own columns to the `ANN` field of `mehari annotate seqvars`.

A plugin is a WebAssembly component: write it in Python, Rust or another language with WebAssembly component tooling, build it into a `.wasm` file, and pass that to `--plugin`, followed by its options:

```sh
mehari annotate seqvars --plugin my_plugin.wasm,data=/path/to/tables,threshold=0.5 ...
```

Plugins run sandboxed and see every transcript annotation of a variant, with its exon structure, its transcript and protein sequences and everything else mehari knows about it.
`--plugin` also takes the names of the built-in plugins (see [Annotating Sequence Variants](`self::anno_seqvars`)) and may be given several times; the columns appear in that order.

## How a plugin works

A plugin implements the `transcript-plugin` interface of [`mehari/wit/plugin.wit`](https://github.com/varfish-org/mehari/blob/main/mehari/wit/plugin.wit), which has two functions:

- `init(setup)` reads the options with `setup.options()` and declares the columns the plugin fills, each with the type of its values: `add-string-column`, `add-integer-column`, `add-float-column` or `add-boolean-column`, or `add-string-list-column` and so on for a list of such values. Columns appear in the order declared; their names may contain letters, digits, `_`, `-` and `.`, and must differ from the fields of `ANN`, such as `Allele` or `HGVS.c`. mehari runs one instance of the plugin per thread and calls `init` once in each.
- `annotate(prediction)` is called once per variant. `prediction.variants()` gives the variant as written in the VCF, before normalization, or the variants of a compound group (`--enable-compound-variants`). `prediction.annotations()` gives one annotation per reported transcript; fill its columns with `set-string`, `set-integer`, `set-float` or `set-boolean`, or `set-string-list` and so on. Fields left unset stay empty.

Both can return an error, which stops mehari with the plugin's message.

Each annotation offers:

- what mehari writes to `ANN`: `feature-id` (transcript accession, empty for intergenic variants), `consequences` (Sequence Ontology terms), `impact`, `gene-symbol`, `gene-id`, `biotypes`, `tags`, `rank` (exon or intron), `hgvs-g`, `hgvs-n`, `hgvs-c`, `hgvs-p`, `cdna-pos`, `cds-pos`, `protein-pos`, `distance` and `warnings`. HGVS normalization shifts the variant towards the transcript's 3' end, so its position depends on the transcript's strand;
- `field(column)`: the value of a column, e.g. one an earlier plugin filled;
- the transcript: `strand`, `exons` with their transcript and genome positions, `contig`, `genome-build`, `protein-id` and `genome-cds-start`/`genome-cds-end`;
- the sequences: `reference-transcript`, `reference-protein`, `reference-cds-start` and `reference-cds-end`, the same as `alternative-*` with the variant applied, and `translation-table`, `translation-exceptions` (e.g. selenocysteine), `frameshift` and `first-changed-codon`.

mehari computes the transcript data when a plugin first asks for it and shares it with the other plugins, so ask only for what you need.
The `host` interface tells the genome assembly, the mehari version and whether `--vep-compatibility-mode` is on.
`plugin.wit` documents every function.

Values go into `ANN` as text: integers and floats as numbers, booleans as `true` or `false`, and lists joined by `&`, as `ANN` does for its own lists such as the consequences.
Strings must not contain whitespace, `|`, `,`, `;`, `=` or `&`, which separate values in VCF and `ANN`, and floats must be finite.
An empty list leaves the field empty.

## Options and data files

Options follow the plugin's path as `,key=value`; the plugin decides what they mean.
`data=<directory>` is mehari's: it lets the plugin read that directory as `/data`, e.g. to look up a table there.
Plugins cannot write it, nor read anything else.
`init` runs once per thread, so load small tables there, and look up large ones on demand.
Option names starting with `mehari-` are reserved for mehari too.

## Writing a plugin in Python

[componentize-py](https://github.com/bytecodealliance/componentize-py) builds a Python module, together with a Python interpreter, into a component.
The commands below run it with [uv](https://docs.astral.sh/uv/); `pip install componentize-py` works as well.
They were tested with componentize-py 0.25.1.

1. Create a directory for the plugin and copy `plugin.wit` from your mehari version into its `wit` subdirectory.
2. Optionally, generate the Python bindings so that your editor knows the types:

   ```sh
   uvx componentize-py -d wit -w plugin bindings .
   ```

3. Write `app.py` with a class `TranscriptPlugin`; raise `Err` to report an error. This example adds a `lof` column for loss-of-function consequences and an `in_panel` column for genes listed in a file:

   ```python
   from componentize_py_types import Err
   from wit_world import exports
   from wit_world.imports.types import Prediction, Setup

   LOSS_OF_FUNCTION = {
       "stop_gained",
       "frameshift_variant",
       "splice_acceptor_variant",
       "splice_donor_variant",
       "start_lost",
   }
   PANEL: set[str] = set()


   class TranscriptPlugin(exports.TranscriptPlugin):
       def init(self, setup: Setup) -> None:
           for key, value in setup.options():
               if key != "panel":
                   raise Err(f"unknown option {key!r}")
               # One gene symbol per line, in the directory given with `data`.
               with open(f"/data/{value}") as panel:
                   PANEL.update(line.strip() for line in panel)
           setup.add_boolean_column("lof")
           setup.add_boolean_column("in_panel")

       def annotate(self, prediction: Prediction) -> None:
           for annotation in prediction.annotations():
               consequences = set(annotation.consequences())
               annotation.set_boolean("lof", bool(LOSS_OF_FUNCTION & consequences))
               annotation.set_boolean("in_panel", annotation.gene_symbol() in PANEL)
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

use exports::mehari::plugin::transcript_plugin::Guest;
use mehari::plugin::types::{Prediction, Setup};

struct MyPlugin;

impl Guest for MyPlugin {
    fn init(setup: &Setup) -> Result<(), String> {
        if let Some((key, _)) = setup.options().first() {
            return Err(format!("unknown option {key:?}"));
        }
        setup.add_integer_column("exon_count");
        Ok(())
    }

    fn annotate(prediction: &Prediction) -> Result<(), String> {
        for annotation in prediction.annotations() {
            let exons = annotation.exons();
            if !exons.is_empty() {
                annotation.set_integer("exon_count", exons.len() as i64);
            }
        }
        Ok(())
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
- mehari stops with an error naming the plugin if the plugin returns an error or fails, sets a column it did not declare or with a value of another type, or uses an annotation after the `annotate` call it came from.
- Two plugins cannot fill the same column.
- The interface has a version (`package mehari:plugin@0.1.0` in `plugin.wit`). Later 0.1.x versions only add functions, so plugins built for an earlier one keep working; a plugin that uses newer functions needs a mehari that has them. A plugin built for another version, e.g. 0.2, fails to load; rebuild it with the `plugin.wit` of your mehari version.
