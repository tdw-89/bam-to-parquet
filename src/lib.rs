pub mod cli;
pub mod fragment;
pub mod header;
pub mod multimap;
pub mod reader;
pub mod stats;
pub mod writer;

use anyhow::Result;
use clap::Parser;

pub fn run() -> Result<()> {
    let args = cli::Args::parse();
    reader::convert(args)
}
