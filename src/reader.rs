use std::{
    collections::HashMap,
    fs::{self, File},
    io::ErrorKind,
    num::NonZeroUsize,
    path::{Path, PathBuf},
    time::SystemTime,
};

use anyhow::{Context, Result, anyhow, bail};
use chrono::{DateTime, SecondsFormat, Utc};
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::alignment::record::data::field::{Tag, Value};
use noodles_sam::alignment::record::{Flags, cigar::op::Kind};
use serde::Serialize;

use crate::{
    cli::{AlignerChoice, Args},
    fragment::{self, Derivation},
    header::{self, BamHeader},
    multimap::{self, Aligner},
    stats::{Counters, ReferenceStats, SkipReason, StatsReport},
    writer::{self, MetadataContext, OutputRow, OutputWriter},
};

#[derive(Eq, Hash, PartialEq)]
enum ReadIdentity {
    Named(Vec<u8>, Vec<u8>, u8),
    Anonymous(u64),
}

#[derive(Serialize)]
struct StatsFile {
    #[serde(flatten)]
    report: StatsReport,
}

pub fn convert(args: Args) -> Result<()> {
    let basename = match args.basename.clone() {
        Some(basename) => basename,
        None => args
            .input
            .file_stem()
            .and_then(|stem| stem.to_str())
            .context("input BAM path has no UTF-8 file stem; pass --basename")?
            .to_owned(),
    };
    writer::validate_basename(&basename)?;
    let outdir = resolve_outdir(args.outdir.as_deref(), &args.input, &basename)?;

    let raw_text = header::read_raw_header(&args.input)?;
    let parsed_header = header::parse(raw_text)?;
    let file_metadata = fs::metadata(&args.input)
        .with_context(|| format!("reading metadata for {}", args.input.display()))?;
    let source_mtime = file_metadata
        .modified()
        .map(format_system_time)
        .unwrap_or_else(|_| "unknown".to_owned());

    if args.single_file && parsed_header.references.len() > 256 {
        bail!("--single-file supports at most 256 references because chrom_id is UInt8");
    }

    let aligner = selected_aligner(args.aligner, &parsed_header);
    eprintln!(
        "aligner: {}; multi-mapper rule: {}",
        aligner.name(),
        aligner.rule()
    );

    let created_utc = Utc::now().to_rfc3339_opts(SecondsFormat::Secs, true);
    let context = MetadataContext {
        basename: basename.clone(),
        source_bam: args.input.display().to_string(),
        source_bam_size_bytes: file_metadata.len(),
        source_bam_mtime: source_mtime.clone(),
        created_utc: created_utc.clone(),
        program_lines: parsed_header.program_lines.clone(),
        multi_rule: aligner.rule().to_owned(),
        git_sha: option_env!("GIT_SHA").unwrap_or("unknown").to_owned(),
    };

    let mut output = OutputWriter::new(
        outdir.clone(),
        basename.clone(),
        args.single_file,
        args.row_group_size,
        args.compression,
        context,
        parsed_header.references.clone(),
    );
    output.initialize()?;
    let mut overall = Counters::default();
    let mut reference_counts = vec![Counters::default(); parsed_header.references.len()];
    let mut frag_ids = HashMap::<ReadIdentity, i32>::new();
    let mut next_frag_id = 0i32;

    let worker_count =
        NonZeroUsize::new(args.threads).context("--threads must be greater than zero")?;
    let file =
        File::open(&args.input).with_context(|| format!("opening BAM {}", args.input.display()))?;
    let decoder = bgzf::io::MultithreadedReader::with_worker_count(worker_count, file);
    let mut reader = bam::io::Reader::from(decoder);
    let _sam_header = reader.read_header().context("reading BAM header")?;
    if _sam_header.reference_sequences().len() != parsed_header.references.len() {
        bail!("raw SAM header and decoded BAM reference dictionaries disagree");
    }

    let mut record = bam::Record::default();
    let mut last_coordinate: Option<(usize, usize)> = None;
    let mut active_reference: Option<usize> = None;
    let mut saw_unplaced_unmapped = false;

    loop {
        if reader
            .read_record(&mut record)
            .context("reading BAM record")?
            == 0
        {
            break;
        }

        overall.record_read();
        let flags = record.flags();
        let raw_flags = flags.bits();
        let record_reference_id = record
            .reference_sequence_id()
            .transpose()
            .context("decoding record reference ID")?;

        if let Some(reference_id) = record_reference_id {
            if reference_id >= parsed_header.references.len() {
                bail!("record reference ID {reference_id} is out of range for the BAM header");
            }
            reference_counts[reference_id].record_read();
        }

        let alignment_start = record
            .alignment_start()
            .transpose()
            .context("decoding alignment start")?
            .map(usize::from);
        let coordinate = match (record_reference_id, alignment_start) {
            (Some(reference_id), Some(position)) => Some((reference_id, position)),
            _ if flags.is_unmapped() => None,
            _ => bail!("mapped record has no complete reference coordinate"),
        };

        if let Some((reference_id, position)) = coordinate {
            if flags.is_unmapped() && saw_unplaced_unmapped {
                bail!("coordinate-bearing unmapped record follows unplaced unmapped records");
            }
            if !flags.is_unmapped() && saw_unplaced_unmapped {
                bail!(
                    "mapped record follows unplaced unmapped records; BAM is not coordinate-sorted"
                );
            }
            let key = (reference_id, position);
            if last_coordinate.is_some_and(|previous| key < previous) {
                bail!(
                    "BAM records are not coordinate-sorted: reference/start {key:?} follows {last_coordinate:?}"
                );
            }
            last_coordinate = Some(key);

            if active_reference.is_some_and(|previous| previous != reference_id) {
                let previous = active_reference.expect("reference change requires prior reference");
                output.finish_reference(previous, &reference_counts[previous])?;
            }
            active_reference = Some(reference_id);
        } else if flags.is_unmapped() {
            saw_unplaced_unmapped = true;
        }

        if flags.is_unmapped() {
            record_skip(
                &mut overall,
                &mut reference_counts,
                record_reference_id,
                SkipReason::Unmapped,
            );
            continue;
        }

        let (reference_id, position) = coordinate.context("mapped record has no coordinate")?;
        let start = i32::try_from(position - 1)
            .context("alignment start exceeds the Int32 Parquet schema")?;
        if flags.is_qc_fail() {
            record_skip(
                &mut overall,
                &mut reference_counts,
                Some(reference_id),
                SkipReason::QcFail,
            );
            continue;
        }
        if args.drop_secondary && flags.is_secondary() {
            record_skip(
                &mut overall,
                &mut reference_counts,
                Some(reference_id),
                SkipReason::Secondary,
            );
            continue;
        }
        if args.drop_dups && flags.is_duplicate() {
            record_skip(
                &mut overall,
                &mut reference_counts,
                Some(reference_id),
                SkipReason::Duplicate,
            );
            continue;
        }

        let mapq = record.mapping_quality().map(u8::from).unwrap_or(u8::MAX);
        if args
            .min_mapq
            .is_some_and(|threshold| mapq != u8::MAX && mapq < threshold)
        {
            record_skip(
                &mut overall,
                &mut reference_counts,
                Some(reference_id),
                SkipReason::MinMapq,
            );
            continue;
        }

        let derivation = fragment::derive(flags, record.template_length());
        let (width, pair_width) = match derivation {
            Derivation::ProperPair { width } => (width, Some(width)),
            Derivation::SkipRightmostMate => {
                record_skip(
                    &mut overall,
                    &mut reference_counts,
                    Some(reference_id),
                    SkipReason::RightmostMate,
                );
                continue;
            }
            Derivation::Singleton => {
                let ops = record
                    .cigar()
                    .iter()
                    .map(|result| result.map(|op| (op.kind(), op.len())))
                    .collect::<std::io::Result<Vec<(Kind, usize)>>>()
                    .context("decoding CIGAR operations")?;
                (fragment::cigar_reference_span(ops)?, None)
            }
        };

        let nh = multimap::nh_value(&record).context("reading NH tag")?;
        let multi =
            multimap::is_multi(&record, aligner, mapq, nh).context("reading multi-mapper tags")?;
        let frag_id = frag_id_for_record(
            &record,
            flags,
            overall.n_records_read,
            &mut frag_ids,
            &mut next_frag_id,
        )?;
        let chrom_id = if args.single_file {
            u8::try_from(reference_id).context("reference ID does not fit in UInt8 chrom_id")?
        } else {
            0
        };

        output.write_row(
            reference_id,
            OutputRow {
                start,
                width,
                mapq,
                flag: raw_flags,
                nh,
                multi,
                frag_id,
                chrom_id,
            },
        )?;
        overall.emit(pair_width, raw_flags);
        reference_counts[reference_id].emit(pair_width, raw_flags);
    }

    if let Some(reference_id) = active_reference {
        output.finish_reference(reference_id, &reference_counts[reference_id])?;
    }
    output.finish_empty_references(&reference_counts)?;
    if !overall.invariant_holds() {
        bail!(
            "record accounting invariant failed: {} emitted + {} skipped != {} records read",
            overall.n_emitted,
            overall.skips.sum(),
            overall.n_records_read
        );
    }
    output.finish(&overall)?;

    let reference_stats = parsed_header
        .references
        .into_iter()
        .zip(reference_counts)
        .map(|(reference, counts)| ReferenceStats { reference, counts })
        .collect();
    let stats = StatsFile {
        report: StatsReport {
            bam2frag_version: env!("CARGO_PKG_VERSION"),
            schema_version: "1",
            created_utc,
            basename: basename.clone(),
            source_bam: args.input.display().to_string(),
            source_bam_size_bytes: file_metadata.len(),
            source_bam_mtime: source_mtime,
            aligner: aligner.name().to_owned(),
            multi_rule: aligner.rule().to_owned(),
            min_mapq: args.min_mapq,
            drop_dups: args.drop_dups,
            drop_secondary: args.drop_secondary,
            single_file: args.single_file,
            overall,
            references: reference_stats,
        },
    };
    let stats_path = args
        .stats
        .unwrap_or_else(|| outdir.join(format!("{basename}.stats.json")));
    let stats_file = File::create(&stats_path)
        .with_context(|| format!("creating stats output {}", stats_path.display()))?;
    serde_json::to_writer_pretty(stats_file, &stats)
        .with_context(|| format!("writing stats output {}", stats_path.display()))?;

    eprintln!(
        "converted {} records into {} rows; stats: {}",
        stats.report.overall.n_records_read,
        stats.report.overall.n_emitted,
        stats_path.display()
    );
    Ok(())
}

/// An explicit `--outdir` must already exist; otherwise derive `<basename>/` beside the
/// input BAM, creating it if needed and reusing it when it is already there.
fn resolve_outdir(requested: Option<&Path>, input: &Path, basename: &str) -> Result<PathBuf> {
    if let Some(outdir) = requested {
        let metadata = fs::metadata(outdir).map_err(|error| {
            if error.kind() == ErrorKind::NotFound {
                anyhow!(
                    "output directory {} does not exist; create it or omit --outdir to derive one from --input",
                    outdir.display()
                )
            } else {
                anyhow::Error::new(error)
                    .context(format!("reading output directory {}", outdir.display()))
            }
        })?;
        if !metadata.is_dir() {
            bail!("output path {} is not a directory", outdir.display());
        }
        return Ok(outdir.to_owned());
    }

    let derived = input
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or(Path::new("."))
        .join(basename);
    fs::create_dir_all(&derived)
        .with_context(|| format!("creating output directory {}", derived.display()))?;
    eprintln!("output directory: {}", derived.display());
    Ok(derived)
}

fn selected_aligner(requested: AlignerChoice, header: &BamHeader) -> Aligner {
    match requested {
        AlignerChoice::Bowtie2 => Aligner::Bowtie2,
        AlignerChoice::Bwa => Aligner::Bwa,
        AlignerChoice::Star => Aligner::Star,
        AlignerChoice::Hisat2 => Aligner::Hisat2,
        AlignerChoice::Auto => match multimap::detect_aligner(&header.program_lines) {
            Some(aligner) => aligner,
            None => {
                eprintln!(
                    "WARNING: could not identify one unambiguous aligner from BAM @PG records; using MAPQ == 0 as the multi-mapper rule. Override with --aligner if you know the aligner."
                );
                Aligner::Unknown
            }
        },
    }
}

fn record_skip(
    overall: &mut Counters,
    references: &mut [Counters],
    reference_id: Option<usize>,
    reason: SkipReason,
) {
    overall.skip(reason);
    if let Some(counts) = reference_id.and_then(|id| references.get_mut(id)) {
        counts.skip(reason);
    }
}

fn frag_id_for_record(
    record: &bam::Record,
    flags: Flags,
    record_number: u64,
    frag_ids: &mut HashMap<ReadIdentity, i32>,
    next_frag_id: &mut i32,
) -> Result<i32> {
    let identity = match record.name() {
        Some(name) if name != b"*" => {
            let side = if flags.is_first_segment() {
                1
            } else if flags.is_last_segment() {
                2
            } else {
                0
            };
            ReadIdentity::Named(name.to_vec(), read_group_id(record)?, side)
        }
        _ => ReadIdentity::Anonymous(record_number),
    };

    if let Some(frag_id) = frag_ids.get(&identity) {
        return Ok(*frag_id);
    }

    let frag_id = *next_frag_id;
    *next_frag_id = next_frag_id
        .checked_add(1)
        .context("more unique read IDs than fit in Int32 frag_id")?;
    frag_ids.insert(identity, frag_id);
    Ok(frag_id)
}

fn read_group_id(record: &bam::Record) -> Result<Vec<u8>> {
    let Some(value) = record
        .data()
        .get(&Tag::new(b'R', b'G'))
        .transpose()
        .context("reading RG tag")?
    else {
        return Ok(Vec::new());
    };

    match value {
        Value::String(read_group) => Ok(read_group.to_vec()),
        _ => Ok(Vec::new()),
    }
}

fn format_system_time(time: SystemTime) -> String {
    DateTime::<Utc>::from(time).to_rfc3339_opts(SecondsFormat::Secs, true)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header_with(program_lines: &[&str]) -> BamHeader {
        BamHeader {
            program_lines: program_lines
                .iter()
                .map(|line| (*line).to_owned())
                .collect(),
            references: Vec::new(),
        }
    }

    #[test]
    fn unknown_header_uses_mapq_fallback() {
        assert_eq!(
            selected_aligner(AlignerChoice::Auto, &header_with(&[])),
            Aligner::Unknown
        );
        assert_eq!(
            selected_aligner(
                AlignerChoice::Auto,
                &header_with(&["@PG\tID:bwa\tCL:bwa mem ref", "@PG\tID:bowtie2\tPN:bowtie2"]),
            ),
            Aligner::Unknown
        );
    }

    #[test]
    fn explicit_aligner_overrides_header_detection() {
        let header = header_with(&["@PG\tID:bowtie2\tPN:bowtie2"]);
        assert_eq!(
            selected_aligner(AlignerChoice::Auto, &header),
            Aligner::Bowtie2
        );
        assert_eq!(
            selected_aligner(AlignerChoice::Star, &header),
            Aligner::Star
        );
    }
}
