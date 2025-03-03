//! Mehari library main entry point.

#[cfg(feature = "documentation")]
#[doc = include_str!("../docs/index.md")]
pub mod user_doc {
    #[doc = include_str!("../docs/getting_started.md")]
    pub mod getting_started {}

    #[doc = include_str!("../docs/anno_seqvars.md")]
    pub mod anno_seqvars {}

    #[doc = include_str!("../docs/anno_strucvars.md")]
    pub mod anno_strucvars {}

    #[doc = include_str!("../docs/db_build.md")]
    pub mod db_build {}

    #[doc = include_str!("../docs/implementation_notes.md")]
    pub mod implementation_notes {}
}

pub mod annotate;
pub mod common;
pub mod db;
pub mod pbs;
pub mod ped;
pub mod server;
pub mod verify;

pub mod plugins;

/// Information about the build.
pub mod built_info {
    include!(concat!(env!("OUT_DIR"), "/built.rs"));
}
