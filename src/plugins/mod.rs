#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use anyhow::Result;
    use extism::convert::Json;
    use extism::{Manifest, Plugin, Wasm};
    use mehari_plugins::*;
    use serde_json;

    #[test]
    fn call_external_plugin() -> Result<()> {
        let record = Dummy {
            some_field: "before".into(),
        };
        let serialized = serde_json::to_string(&record)?;

        let mut plugin_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        plugin_path.push("target/wasm32-unknown-unknown/debug/mehari_plugin_vep_nmd.wasm");
        let url = Wasm::file(plugin_path);
        let manifest = Manifest::new([url]);
        let mut plugin = Plugin::new(&manifest, [], true)?;
        let Json(result) = plugin.call::<String, Json<Dummy>>("process", serialized)?;

        assert_ne!(result, record);
        assert_eq!(result.some_field, "after");
        Ok(())
    }
}
