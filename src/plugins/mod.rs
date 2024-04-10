#[cfg(test)]
mod tests {
    use anyhow::Result;
    use extism::convert::Json;
    use extism::{Manifest, Plugin, Wasm};
    use mehari_plugins::*;
    use std::path::PathBuf;

    #[test]
    fn nmd_plugin_header_information() -> Result<()> {
        let mut plugin_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        plugin_path.push("target/wasm32-unknown-unknown/debug/mehari_plugin_vep_nmd.wasm");
        let url = Wasm::file(plugin_path);
        let manifest = Manifest::new([url]);
        let mut plugin = Plugin::new(&manifest, [], true)?;

        let Json(header_info) = plugin.call::<(), Json<HeaderInfo>>("header_info", ())?;

        assert_eq!(
            header_info,
            HeaderInfo {
                tag: "NMD".into(),
                description: "Nonsense-mediated mRNA decay escaping variants prediction".into(),
            }
        );
        Ok(())
    }

    #[test]
    fn call_nmd_plugin() -> Result<()> {
        let tva = TranscriptVariationAlleleBuilder::default()
            .transcript(
                TranscriptBuilder::default()
                    .introns(vec![Intron::new(10, 20), Intron::new(30, 40)])
                    .exons(vec![
                        Exon::new(0, 10),
                        Exon::new(20, 30),
                        Exon::new(40, 100),
                    ])
                    .strand(-1)
                    .build()?,
            )
            .transcript_variation(TranscriptVariation::new(CodingSequence::new(0, 100)))
            .variation_feature(VariationFeature::new(SeqRegion::new(15, 25)))
            .overlap_consequences(vec![OverlapConsequence::new("stop_gained".into())])
            .build()?;

        let mut plugin_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        plugin_path.push("target/wasm32-unknown-unknown/debug/mehari_plugin_vep_nmd.wasm");
        let url = Wasm::file(plugin_path);
        let manifest = Manifest::new([url]);
        let mut plugin = Plugin::new(&manifest, [], true)?;
        if let Some(Json(result)) = plugin
            .call::<Json<TranscriptVariationAllele>, Option<Json<Annotation>>>(
                "process",
                Json(tva),
            )?
        {
            assert_eq!(result, Annotation::new("NMD".into()));
        } else {
            panic!()
        }
        Ok(())
    }
}
