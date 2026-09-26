use std::{fs::File, io::Read, path::Path};

use anyhow::{Context, Result, bail};
use noodles_bam as bam;

#[derive(Clone, Debug, serde::Serialize)]
pub struct ReferenceSequence {
    pub name: String,
    pub length: u64,
    pub md5: Option<String>,
}

#[derive(Clone, Debug)]
pub struct BamHeader {
    pub program_lines: Vec<String>,
    pub references: Vec<ReferenceSequence>,
}

pub fn read_raw_header(path: &Path) -> Result<String> {
    let file = File::open(path).with_context(|| format!("opening BAM {}", path.display()))?;
    let mut reader = bam::io::Reader::new(file);
    let mut header_reader = reader.header_reader();
    header_reader
        .read_magic_number()
        .context("reading BAM magic number")?;
    let mut sam_header_reader = header_reader
        .raw_sam_header_reader()
        .context("opening raw SAM header reader")?;
    let mut raw_text = String::new();
    sam_header_reader
        .read_to_string(&mut raw_text)
        .context("reading raw SAM header")?;
    Ok(raw_text)
}

pub fn parse(raw_text: String) -> Result<BamHeader> {
    let hd_line = raw_text
        .lines()
        .find(|line| line.starts_with("@HD\t") || *line == "@HD")
        .context("BAM header is missing @HD")?;

    let sort_order = fields(hd_line)
        .find(|(key, _)| *key == "SO")
        .map(|(_, value)| value);
    if sort_order != Some("coordinate") {
        bail!("BAM header must declare @HD SO:coordinate (found {sort_order:?})");
    }

    let mut references = Vec::new();
    let mut program_lines = Vec::new();

    for line in raw_text.lines() {
        if line.starts_with("@SQ\t") {
            let mut name = None;
            let mut length = None;
            let mut md5 = None;

            for (key, value) in fields(line) {
                match key {
                    "SN" => name = Some(value.to_owned()),
                    "LN" => {
                        length = Some(
                            value
                                .parse::<u64>()
                                .with_context(|| format!("invalid @SQ length in {line:?}"))?,
                        )
                    }
                    "M5" => md5 = Some(value.to_owned()),
                    _ => {}
                }
            }

            references.push(ReferenceSequence {
                name: name.context("@SQ record is missing SN")?,
                length: length.context("@SQ record is missing LN")?,
                md5,
            });
        } else if line.starts_with("@PG\t") || line == "@PG" {
            program_lines.push(line.to_owned());
        }
    }

    Ok(BamHeader {
        program_lines,
        references,
    })
}

fn fields(line: &str) -> impl Iterator<Item = (&str, &str)> {
    line.split('\t')
        .skip(1)
        .filter_map(|field| field.split_once(':'))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_coordinate_header_and_sq_metadata() {
        let header = parse(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\tM5:abcd\n@PG\tID:bwa\tPN:bwa\tCL:bwa mem ref in.bam\n".to_owned(),
        )
        .unwrap();

        assert_eq!(header.references[0].name, "chr1");
        assert_eq!(header.references[0].length, 100);
        assert_eq!(header.references[0].md5.as_deref(), Some("abcd"));
        assert_eq!(header.program_lines.len(), 1);
    }

    #[test]
    fn rejects_missing_or_noncoordinate_sort_order() {
        assert!(parse("@HD\tVN:1.6\tSO:queryname\n".to_owned()).is_err());
        assert!(parse("@SQ\tSN:chr1\tLN:100\n".to_owned()).is_err());
    }
}
