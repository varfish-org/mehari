//! Dump schema of the REST API server.

use std::{
    fs::File,
    io::{self, Write},
};

use utoipa::OpenApi as _;

use crate::server::run::openapi::ApiDoc;

/// Command line arguments for `server schema` sub command.
#[derive(clap::Parser, Debug, Clone)]
#[command(author, version, about = "Dump REST API schema", long_about = None)]
pub struct Args {
    /// Path to the output file.  Use stdout if missing.
    #[arg(long)]
    pub output_file: Option<String>,
}

impl Args {
    /// Get writeable output file or stdout.
    fn get_output(&self) -> Result<Box<dyn Write>, io::Error> {
        match self.output_file {
            Some(ref path) => File::create(path).map(|f| Box::new(f) as Box<dyn Write>),
            None => Ok(Box::new(io::stdout())),
        }
    }
}

/// Main entry point for `server run` sub command.
///
/// # Errors
///
/// In the case that there is an error running the server.
pub fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!("args_common = {:?}", &args_common);
    tracing::info!("args = {:?}", &args);

    let schema_yaml = ApiDoc::openapi()
        .to_yaml()
        .map_err(|e| anyhow::anyhow!("Failed to convert OpenAPI to YAML: {}", e))?;
    let mut output = args
        .get_output()
        .map_err(|e| anyhow::anyhow!("Failed to open output file: {}", e))?;
    write!(output, "{}", &schema_yaml)
        .map_err(|e| anyhow::anyhow!("Failed to write output: {}", e))?;

    tracing::info!("All done. Have a nice day!");
    Ok(())
}
