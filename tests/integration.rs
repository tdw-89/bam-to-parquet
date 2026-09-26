use std::{fs, fs::File, path::Path};

use arrow::array::{Int32Array, UInt8Array};
use bam_to_parquet::{
    cli::{AlignerChoice, Args, CompressionSpec},
    reader::convert,
};
use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::{
        io::Write as _,
        record::{
            Flags, MappingQuality,
            cigar::{Op, op::Kind},
        },
        record_buf::Cigar,
    },
};
use parquet::{arrow::arrow_reader::ParquetRecordBatchReaderBuilder, file::reader::FileReader};
use serde_json::Value as JsonValue;
use tempfile::tempdir;

#[test]
fn converts_bam_and_round_trips_per_chromosome_and_single_file_outputs() {
    let temp = tempdir().unwrap();
    let input = temp.path().join("tiny.bam");
    write_fixture(&input, true);

    let per_chrom_dir = temp.path().join("per-chrom");
    fs::create_dir_all(&per_chrom_dir).unwrap();
    convert(test_args(&input, Some(&per_chrom_dir), false)).unwrap();

    let chr1 = per_chrom_dir.join("tiny.chr1.parquet");
    let chr2 = per_chrom_dir.join("tiny.chr2.parquet");
    let chr3 = per_chrom_dir.join("tiny.chr3.parquet");
    assert_eq!(
        read_i32_column(&chr1, "start"),
        vec![100, 300, 500, 550, 600]
    );
    assert_eq!(read_i32_column(&chr1, "width"), vec![150, 12, 10, 6, 110]);
    assert_eq!(read_i32_column(&chr1, "frag_id"), vec![0, 1, 2, 3, 0]);
    assert_eq!(read_u8_column(&chr1, "mapq"), vec![40, 40, 40, 255, 30]);
    assert_eq!(read_i32_column(&chr2, "start"), vec![4]);
    assert!(read_i32_column(&chr3, "start").is_empty());
    assert_eq!(metadata_value(&chr1, "max_width").as_deref(), Some("150"));
    assert_eq!(metadata_value(&chr1, "basename").as_deref(), Some("tiny"));
    assert_eq!(metadata_value(&chr1, "sample"), None);
    assert_eq!(metadata_value(&chr3, "n_emitted").as_deref(), Some("0"));

    let stats_path = per_chrom_dir.join("tiny.stats.json");
    let stats: JsonValue = serde_json::from_reader(File::open(stats_path).unwrap()).unwrap();
    assert_eq!(stats["basename"], "tiny");
    assert_eq!(stats["sample"], JsonValue::Null);
    assert_eq!(stats["aligner"], "bowtie2");
    assert_eq!(stats["overall"]["n_records_read"], 9);
    assert_eq!(stats["overall"]["n_emitted"], 6);
    assert_eq!(stats["overall"]["n_pairs"], 2);
    assert_eq!(stats["overall"]["n_singletons"], 4);
    assert_eq!(stats["overall"]["skips"]["unmapped"], 1);
    assert_eq!(stats["overall"]["skips"]["min_mapq"], 1);
    assert_eq!(stats["overall"]["skips"]["rightmost_mate"], 1);
    assert_eq!(
        stats["overall"]["fragment_length_histogram"]["bins_1_to_1000"][149],
        1
    );
    assert_eq!(
        stats["overall"]["fragment_length_histogram"]["bins_1_to_1000"][109],
        1
    );

    let single_dir = temp.path().join("single");
    fs::create_dir_all(&single_dir).unwrap();
    convert(test_args(&input, Some(&single_dir), true)).unwrap();
    let single_path = single_dir.join("tiny.parquet");
    assert_eq!(
        read_i32_column(&single_path, "start"),
        vec![100, 300, 500, 550, 600, 4]
    );
    assert_eq!(
        read_u8_column(&single_path, "chrom_id"),
        vec![0, 0, 0, 0, 0, 1]
    );

    let max_width = 150;
    let query_start = 200;
    let query_end = 210;
    let starts = read_i32_column(&chr1, "start");
    let widths = read_i32_column(&chr1, "width");
    let range_hits = starts
        .iter()
        .zip(widths)
        .filter(|(start, width)| {
            **start >= query_start - max_width
                && **start <= query_end
                && **start + *width > query_start
                && **start < query_end
        })
        .count();
    assert_eq!(range_hits, 1);
}

#[test]
fn rejects_records_out_of_coordinate_order() {
    let temp = tempdir().unwrap();
    let input = temp.path().join("unsorted.bam");
    let header: sam::Header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
        .parse()
        .unwrap();
    let mut writer = bam::io::Writer::new(File::create(&input).unwrap());
    writer.write_header(&header).unwrap();
    for position in [201, 101] {
        let record = make_record(
            "read",
            Flags::empty(),
            Some(0),
            Some(position),
            Some(40),
            0,
            vec![(Kind::Match, 10)],
        );
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();

    let out = temp.path().join("out");
    fs::create_dir_all(&out).unwrap();
    let mut args = test_args(&input, Some(&out), false);
    args.min_mapq = None;
    let error = convert(args).unwrap_err();
    assert!(error.to_string().contains("not coordinate-sorted"));
}

#[test]
fn derives_output_directory_from_the_input_bam_when_outdir_is_omitted() {
    let temp = tempdir().unwrap();
    let nested = temp.path().join("bams");
    fs::create_dir_all(&nested).unwrap();
    let input = nested.join("tiny.bam");
    write_fixture(&input, true);

    convert(test_args(&input, None, false)).unwrap();

    let derived = nested.join("tiny");
    assert!(derived.is_dir());
    assert_eq!(
        read_i32_column(&derived.join("tiny.chr1.parquet"), "start").len(),
        5
    );
    assert!(derived.join("tiny.stats.json").is_file());
}

#[test]
fn reuses_an_existing_derived_output_directory() {
    let temp = tempdir().unwrap();
    let input = temp.path().join("tiny.bam");
    write_fixture(&input, true);
    let derived = temp.path().join("tiny");
    fs::create_dir_all(&derived).unwrap();
    fs::write(derived.join("notes.txt"), "keep me").unwrap();

    convert(test_args(&input, None, false)).unwrap();
    convert(test_args(&input, None, false)).unwrap();

    assert_eq!(
        fs::read_to_string(derived.join("notes.txt")).unwrap(),
        "keep me"
    );
    assert_eq!(
        read_i32_column(&derived.join("tiny.chr1.parquet"), "start").len(),
        5
    );
}

#[test]
fn basename_names_both_the_derived_directory_and_the_output_files() {
    let temp = tempdir().unwrap();
    let input = temp.path().join("tiny.bam");
    write_fixture(&input, true);

    let mut args = test_args(&input, None, false);
    args.basename = Some("H3K9me3_rep1".to_owned());
    convert(args).unwrap();

    let derived = temp.path().join("H3K9me3_rep1");
    assert!(derived.is_dir());
    assert!(!temp.path().join("tiny").exists());
    assert_eq!(
        read_i32_column(&derived.join("H3K9me3_rep1.chr1.parquet"), "start").len(),
        5
    );
    assert!(derived.join("H3K9me3_rep1.stats.json").is_file());
}

#[test]
fn rejects_an_explicit_output_directory_that_does_not_exist() {
    let temp = tempdir().unwrap();
    let input = temp.path().join("tiny.bam");
    write_fixture(&input, true);

    let missing = temp.path().join("absent");
    let error = convert(test_args(&input, Some(&missing), false)).unwrap_err();
    assert!(error.to_string().contains("does not exist"));
    assert!(!missing.exists());
}

#[test]
fn rejects_an_explicit_output_directory_that_is_a_file() {
    let temp = tempdir().unwrap();
    let input = temp.path().join("tiny.bam");
    write_fixture(&input, true);
    let not_a_dir = temp.path().join("file.txt");
    fs::write(&not_a_dir, "").unwrap();

    let error = convert(test_args(&input, Some(&not_a_dir), false)).unwrap_err();
    assert!(error.to_string().contains("is not a directory"));
}

fn test_args(input: &Path, outdir: Option<&Path>, single_file: bool) -> Args {
    Args {
        input: input.to_owned(),
        outdir: outdir.map(Path::to_owned),
        basename: None,
        threads: 2,
        single_file,
        row_group_size: 2,
        compression: CompressionSpec::None,
        aligner: AlignerChoice::Auto,
        min_mapq: Some(20),
        drop_dups: false,
        drop_secondary: false,
        stats: None,
    }
}

fn write_fixture(path: &Path, with_bowtie2_program: bool) {
    let program = if with_bowtie2_program {
        "@PG\tID:bowtie2\tPN:bowtie2\tVN:2.5.1\tCL:bowtie2 --very-sensitive\n"
    } else {
        ""
    };
    let header_text = format!(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n@SQ\tSN:chr3\tLN:1000\n{program}"
    );
    let header: sam::Header = header_text.parse().unwrap();
    let mut writer = bam::io::Writer::new(File::create(path).unwrap());
    writer.write_header(&header).unwrap();

    let proper_first = Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED | Flags::FIRST_SEGMENT;
    let proper_last = Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED | Flags::LAST_SEGMENT;
    let records = [
        make_record(
            "pair-1",
            proper_first,
            Some(0),
            Some(101),
            Some(40),
            150,
            vec![(Kind::Match, 50)],
        ),
        make_record(
            "pair-1",
            proper_last,
            Some(0),
            Some(201),
            Some(40),
            -150,
            vec![(Kind::Match, 50)],
        ),
        make_record(
            "single",
            Flags::empty(),
            Some(0),
            Some(301),
            Some(40),
            0,
            vec![
                (Kind::Match, 5),
                (Kind::Insertion, 2),
                (Kind::Deletion, 3),
                (Kind::Skip, 4),
                (Kind::SoftClip, 1),
            ],
        ),
        make_record(
            "low-mapq",
            Flags::empty(),
            Some(0),
            Some(401),
            Some(10),
            0,
            vec![(Kind::Match, 10)],
        ),
        make_record(
            "duplicate",
            Flags::DUPLICATE,
            Some(0),
            Some(501),
            Some(40),
            0,
            vec![(Kind::Match, 10)],
        ),
        make_record(
            "unknown-mapq",
            Flags::empty(),
            Some(0),
            Some(551),
            None,
            0,
            vec![(Kind::Match, 6)],
        ),
        make_record(
            "pair-1",
            proper_first | Flags::SECONDARY,
            Some(0),
            Some(601),
            Some(30),
            110,
            vec![(Kind::Match, 50)],
        ),
        make_record(
            "chr2-read",
            Flags::empty(),
            Some(1),
            Some(5),
            Some(40),
            0,
            vec![(Kind::Match, 7)],
        ),
        make_record("unmapped", Flags::UNMAPPED, None, None, None, 0, Vec::new()),
    ];

    for record in records {
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
}

fn make_record(
    name: &str,
    flags: Flags,
    reference_id: Option<usize>,
    position: Option<usize>,
    mapq: Option<u8>,
    template_length: i32,
    cigar_ops: Vec<(Kind, usize)>,
) -> sam::alignment::RecordBuf {
    let cigar: Cigar = cigar_ops
        .into_iter()
        .map(|(kind, length)| Op::new(kind, length))
        .collect();
    let mut builder = sam::alignment::RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_cigar(cigar)
        .set_template_length(template_length);

    if let Some(reference_id) = reference_id {
        builder = builder.set_reference_sequence_id(reference_id);
    }
    if let Some(position) = position {
        builder = builder.set_alignment_start(Position::try_from(position).unwrap());
    }
    if let Some(mapq) = mapq {
        builder = builder.set_mapping_quality(MappingQuality::new(mapq).unwrap());
    }

    builder.build()
}

fn read_i32_column(path: &Path, column_name: &str) -> Vec<i32> {
    let builder = ParquetRecordBatchReaderBuilder::try_new(File::open(path).unwrap()).unwrap();
    let column_index = builder.schema().index_of(column_name).unwrap();
    let reader = builder.with_batch_size(2).build().unwrap();
    reader
        .flat_map(|batch| {
            let batch = batch.unwrap();
            batch
                .column(column_index)
                .as_any()
                .downcast_ref::<Int32Array>()
                .unwrap()
                .values()
                .to_vec()
        })
        .collect()
}

fn read_u8_column(path: &Path, column_name: &str) -> Vec<u8> {
    let builder = ParquetRecordBatchReaderBuilder::try_new(File::open(path).unwrap()).unwrap();
    let column_index = builder.schema().index_of(column_name).unwrap();
    let reader = builder.with_batch_size(2).build().unwrap();
    reader
        .flat_map(|batch| {
            let batch = batch.unwrap();
            batch
                .column(column_index)
                .as_any()
                .downcast_ref::<UInt8Array>()
                .unwrap()
                .values()
                .to_vec()
        })
        .collect()
}

fn metadata_value(path: &Path, key: &str) -> Option<String> {
    let reader =
        parquet::file::reader::SerializedFileReader::new(File::open(path).unwrap()).unwrap();
    reader
        .metadata()
        .file_metadata()
        .key_value_metadata()?
        .iter()
        .find(|item| item.key == key)
        .and_then(|item| item.value.clone())
}
