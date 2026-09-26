use std::{fmt, path::PathBuf, str::FromStr};

use anyhow::{Context, Result, bail};
use clap::Parser;
use clap::ValueEnum;

#[derive(Debug, Parser)]
#[command(
    name = "bam2frag",
    version,
    about = "Convert coordinate-sorted BAM files to Parquet fragment tables"
)]
pub struct Args {
    /// Coordinate-sorted input BAM.
    #[arg(short, long)]
    pub input: PathBuf,

    /// Existing output directory (defaults to <basename>/ beside the input BAM).
    #[arg(short, long)]
    pub outdir: Option<PathBuf>,

    /// Output directory and file name prefix (defaults to the input file stem).
    #[arg(short, long)]
    pub basename: Option<String>,

    /// BGZF decompression workers.
    #[arg(short, long, default_value_t = 4, value_parser = positive_usize)]
    pub threads: usize,

    /// Write one Parquet file with a chrom_id column.
    #[arg(short = '1', long)]
    pub single_file: bool,

    /// Maximum rows in each Parquet row group.
    #[arg(short, long, default_value_t = 1_000_000, value_parser = positive_usize)]
    pub row_group_size: usize,

    /// Compression: zstd[:level], snappy, or none.
    #[arg(short, long, default_value = "zstd:3")]
    pub compression: CompressionSpec,

    /// Multi-mapper detection rule; auto-detect from BAM @PG records by default.
    #[arg(short, long, value_enum, default_value = "auto")]
    pub aligner: AlignerChoice,

    /// Drop records below this MAPQ (valid scores are 1 through 254).
    #[arg(short = 'q', long, value_parser = clap::value_parser!(u8).range(1..=254))]
    pub min_mapq: Option<u8>,

    /// Drop records marked duplicate (SAM flag 0x400).
    #[arg(short = 'd', long)]
    pub drop_dups: bool,

    /// Drop secondary alignments (SAM flag 0x100).
    #[arg(short = 'D', long)]
    pub drop_secondary: bool,

    /// Stats JSON path (defaults to <outdir>/<basename>.stats.json).
    #[arg(short = 'S', long)]
    pub stats: Option<PathBuf>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, ValueEnum)]
pub enum AlignerChoice {
    #[default]
    Auto,
    Bowtie2,
    Bwa,
    Star,
    Hisat2,
}

fn positive_usize(value: &str) -> std::result::Result<usize, String> {
    let n = value.parse::<usize>().map_err(|error| error.to_string())?;
    if n == 0 {
        return Err("must be greater than zero".to_owned());
    }
    Ok(n)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CompressionSpec {
    Zstd(i32),
    Snappy,
    None,
}

impl FromStr for CompressionSpec {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self> {
        match value.to_ascii_lowercase().as_str() {
            "snappy" => Ok(Self::Snappy),
            "none" => Ok(Self::None),
            "zstd" => Ok(Self::Zstd(3)),
            other if other.starts_with("zstd:") => {
                let level = other[5..]
                    .parse::<i32>()
                    .context("Zstandard level must be an integer")?;
                if !(-131_072..=22).contains(&level) {
                    bail!("Zstandard level must be between -131072 and 22");
                }
                Ok(Self::Zstd(level))
            }
            _ => bail!("compression must be zstd[:level], snappy, or none"),
        }
    }
}

impl fmt::Display for CompressionSpec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Zstd(level) => write!(f, "zstd:{level}"),
            Self::Snappy => f.write_str("snappy"),
            Self::None => f.write_str("none"),
        }
    }
}

#[cfg(test)]
mod tests {
    use clap::Parser;

    use super::Args;

    #[test]
    fn mapq_filter_requires_a_valid_positive_score() {
        let base = ["bam2frag", "--input", "in.bam"];
        assert!(Args::try_parse_from([base.as_slice(), &["--min-mapq", "1"]].concat()).is_ok());
        assert!(Args::try_parse_from([base.as_slice(), &["--min-mapq", "254"]].concat()).is_ok());
        assert!(Args::try_parse_from([base.as_slice(), &["--min-mapq", "0"]].concat()).is_err());
        assert!(Args::try_parse_from([base.as_slice(), &["--min-mapq", "255"]].concat()).is_err());
        assert!(Args::try_parse_from([base.as_slice(), &["--min-mapq", "2.5"]].concat()).is_err());
    }

    #[test]
    fn short_forms_parse_to_the_same_arguments_as_long_forms() {
        let long = Args::try_parse_from([
            "bam2frag",
            "--input",
            "in.bam",
            "--outdir",
            "out",
            "--basename",
            "s1",
            "--threads",
            "8",
            "--single-file",
            "--row-group-size",
            "100",
            "--compression",
            "snappy",
            "--aligner",
            "star",
            "--min-mapq",
            "30",
            "--drop-dups",
            "--drop-secondary",
            "--stats",
            "s.json",
        ])
        .unwrap();
        let short = Args::try_parse_from([
            "bam2frag", "-i", "in.bam", "-o", "out", "-b", "s1", "-t", "8", "-1", "-r", "100",
            "-c", "snappy", "-a", "star", "-q", "30", "-d", "-D", "-S", "s.json",
        ])
        .unwrap();

        assert_eq!(format!("{long:?}"), format!("{short:?}"));
        assert_eq!(short.basename.as_deref(), Some("s1"));
        assert!(short.single_file && short.drop_dups && short.drop_secondary);
    }

    #[test]
    fn output_directory_is_optional() {
        let args = Args::try_parse_from(["bam2frag", "-i", "in.bam"]).unwrap();
        assert!(args.outdir.is_none());
        assert!(args.basename.is_none());
    }
}
