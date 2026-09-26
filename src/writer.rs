use std::{
    collections::HashSet,
    fs::File,
    path::{Path, PathBuf},
    sync::Arc,
};

use anyhow::{Context, Result, bail};
use arrow::{
    array::{ArrayBuilder, ArrayRef, BooleanBuilder, Int32Builder, UInt8Builder, UInt16Builder},
    datatypes::{DataType, Field, Schema},
    record_batch::RecordBatch,
};
use parquet::{
    arrow::arrow_writer::ArrowWriter,
    basic::{Compression, Encoding},
    file::{
        metadata::{KeyValue, SortingColumn},
        properties::{EnabledStatistics, WriterProperties},
    },
    schema::types::ColumnPath,
};

use crate::{cli::CompressionSpec, header::ReferenceSequence, stats::Counters};

const BATCH_SIZE: usize = 65_536;
const DATA_PAGE_SIZE: usize = 1_000_000;

#[derive(Clone, Copy, Debug)]
pub struct OutputRow {
    pub start: i32,
    pub width: i32,
    pub mapq: u8,
    pub flag: u16,
    pub nh: u8,
    pub multi: bool,
    pub frag_id: i32,
    pub chrom_id: u8,
}

#[derive(Clone, Debug)]
pub struct MetadataContext {
    pub basename: String,
    pub source_bam: String,
    pub source_bam_size_bytes: u64,
    pub source_bam_mtime: String,
    pub created_utc: String,
    pub program_lines: Vec<String>,
    pub multi_rule: String,
    pub git_sha: String,
}

pub struct OutputWriter {
    outdir: PathBuf,
    basename: String,
    single_file: bool,
    row_group_size: usize,
    compression: CompressionSpec,
    context: MetadataContext,
    references: Vec<ReferenceSequence>,
    schema: Arc<Schema>,
    current_chromosome: Option<(usize, ParquetSink)>,
    single_sink: Option<ParquetSink>,
    written_references: HashSet<usize>,
}

impl OutputWriter {
    pub fn new(
        outdir: PathBuf,
        basename: String,
        single_file: bool,
        row_group_size: usize,
        compression: CompressionSpec,
        context: MetadataContext,
        references: Vec<ReferenceSequence>,
    ) -> Self {
        let schema = Arc::new(make_schema(single_file));
        Self {
            outdir,
            basename,
            single_file,
            row_group_size,
            compression,
            context,
            references,
            schema,
            current_chromosome: None,
            single_sink: None,
            written_references: HashSet::new(),
        }
    }

    pub fn initialize(&mut self) -> Result<()> {
        if self.single_file && self.single_sink.is_none() {
            let path = self.outdir.join(format!("{}.parquet", self.basename));
            self.single_sink = Some(self.open_sink(path, None)?);
        }
        Ok(())
    }

    pub fn write_row(&mut self, reference_id: usize, row: OutputRow) -> Result<()> {
        if self.single_file {
            if self.single_sink.is_none() {
                let path = self.outdir.join(format!("{}.parquet", self.basename));
                self.single_sink = Some(self.open_sink(path, None)?);
            }
            self.single_sink
                .as_mut()
                .expect("single-file sink was just initialized")
                .write(row, Some(row.chrom_id))
        } else {
            if let Some((current_id, _)) = &self.current_chromosome {
                if *current_id != reference_id {
                    bail!(
                        "chromosome writer for reference {current_id} was not finalized before reference {reference_id}"
                    );
                }
            } else {
                let reference = self
                    .references
                    .get(reference_id)
                    .context("record reference ID is out of range")?;
                let file_stem =
                    format!("{}.{}", self.basename, safe_file_component(&reference.name));
                let path = self.outdir.join(format!("{file_stem}.parquet"));
                let sink = self.open_sink(path, Some(reference))?;
                self.current_chromosome = Some((reference_id, sink));
                self.written_references.insert(reference_id);
            }

            self.current_chromosome
                .as_mut()
                .expect("per-chromosome sink was just initialized")
                .1
                .write(row, None)
        }
    }

    pub fn finish_reference(&mut self, reference_id: usize, counts: &Counters) -> Result<()> {
        if self.single_file {
            return Ok(());
        }

        if let Some((current_id, sink)) = self.current_chromosome.take() {
            if current_id != reference_id {
                self.current_chromosome = Some((current_id, sink));
                bail!(
                    "attempted to finalize reference {reference_id} while reference {current_id} is open"
                );
            }
            sink.close(counts)?;
        }

        Ok(())
    }

    pub fn finish_empty_references(&mut self, counts: &[Counters]) -> Result<()> {
        if self.single_file {
            return Ok(());
        }
        if self.current_chromosome.is_some() {
            bail!(
                "cannot create empty chromosome outputs before the active reference is finalized"
            );
        }
        if counts.len() != self.references.len() {
            bail!("reference counters do not match the BAM header reference count");
        }

        for (reference_id, chromosome_counts) in counts.iter().enumerate() {
            if self.written_references.contains(&reference_id) {
                continue;
            }
            let reference = &self.references[reference_id];
            let file_stem = format!("{}.{}", self.basename, safe_file_component(&reference.name));
            let path = self.outdir.join(format!("{file_stem}.parquet"));
            self.open_sink(path, Some(reference))?
                .close(chromosome_counts)?;
            self.written_references.insert(reference_id);
        }

        Ok(())
    }

    pub fn finish(mut self, overall: &Counters) -> Result<()> {
        if self.current_chromosome.is_some() {
            bail!(
                "per-chromosome output must be finalized with its reference counters before closing"
            );
        }

        if let Some(sink) = self.single_sink.take() {
            sink.close_single(overall, &self.references)?;
        }

        Ok(())
    }

    fn open_sink(
        &self,
        path: PathBuf,
        reference: Option<&ReferenceSequence>,
    ) -> Result<ParquetSink> {
        let file = File::create(&path)
            .with_context(|| format!("creating Parquet output {}", path.display()))?;
        let properties =
            writer_properties(self.row_group_size, self.compression, self.single_file)?;
        let writer = ArrowWriter::try_new(file, self.schema.clone(), Some(properties))
            .with_context(|| format!("opening Parquet writer for {}", path.display()))?;
        let mut sink = ParquetSink {
            path,
            writer,
            schema: self.schema.clone(),
            single_file: self.single_file,
            builders: BatchBuilders::new(),
            last_key: None,
            max_width: 0,
        };

        for (key, value) in self.base_metadata(reference)? {
            sink.append_metadata(key, value);
        }

        Ok(sink)
    }

    fn base_metadata(
        &self,
        reference: Option<&ReferenceSequence>,
    ) -> Result<Vec<(String, String)>> {
        let mut entries = vec![
            (
                "bam2frag_version".to_owned(),
                env!("CARGO_PKG_VERSION").to_owned(),
            ),
            ("git_sha".to_owned(), self.context.git_sha.clone()),
            ("schema_version".to_owned(), "1".to_owned()),
            ("created_utc".to_owned(), self.context.created_utc.clone()),
            ("basename".to_owned(), self.context.basename.clone()),
            ("source_bam".to_owned(), self.context.source_bam.clone()),
            (
                "source_bam_size_bytes".to_owned(),
                self.context.source_bam_size_bytes.to_string(),
            ),
            (
                "source_bam_mtime".to_owned(),
                self.context.source_bam_mtime.clone(),
            ),
            (
                "bam_header_pg".to_owned(),
                serde_json::to_string(&self.context.program_lines)?,
            ),
            ("multi_rule".to_owned(), self.context.multi_rule.clone()),
        ];

        if let Some(reference) = reference {
            entries.extend([
                ("reference_name".to_owned(), reference.name.clone()),
                ("reference_length".to_owned(), reference.length.to_string()),
                (
                    "reference_md5".to_owned(),
                    reference.md5.clone().unwrap_or_default(),
                ),
            ]);
        }

        Ok(entries)
    }
}

struct ParquetSink {
    path: PathBuf,
    writer: ArrowWriter<File>,
    schema: Arc<Schema>,
    single_file: bool,
    builders: BatchBuilders,
    last_key: Option<(u8, i32)>,
    max_width: i32,
}

impl ParquetSink {
    fn write(&mut self, row: OutputRow, chrom_id: Option<u8>) -> Result<()> {
        let key = (chrom_id.unwrap_or(0), row.start);
        if self.last_key.is_some_and(|last| key < last) {
            bail!(
                "output sorting invariant failed in {}: key {key:?} follows {:?}",
                self.path.display(),
                self.last_key
            );
        }
        self.last_key = Some(key);
        self.max_width = self.max_width.max(row.width);
        self.builders.append(row, self.single_file);
        if self.builders.len() >= BATCH_SIZE {
            self.flush_batch()?;
        }
        Ok(())
    }

    fn flush_batch(&mut self) -> Result<()> {
        if self.builders.is_empty() {
            return Ok(());
        }

        let batch = self.builders.finish(&self.schema, self.single_file)?;
        self.writer
            .write(&batch)
            .with_context(|| format!("writing Parquet batch to {}", self.path.display()))
    }

    fn append_metadata(&mut self, key: String, value: String) {
        self.writer
            .append_key_value_metadata(KeyValue::new(key, Some(value)));
    }

    fn close(mut self, counts: &Counters) -> Result<()> {
        self.flush_batch()?;
        self.append_final_metadata(counts)?;
        self.writer
            .close()
            .with_context(|| format!("finalizing Parquet file {}", self.path.display()))?;
        Ok(())
    }

    fn close_single(mut self, counts: &Counters, references: &[ReferenceSequence]) -> Result<()> {
        let table = references
            .iter()
            .enumerate()
            .map(|(index, reference)| (index as u8, reference.name.as_str()))
            .collect::<Vec<_>>();
        self.append_metadata("chrom_table".to_owned(), serde_json::to_string(&table)?);
        self.flush_batch()?;
        self.append_final_metadata(counts)?;
        self.writer
            .close()
            .with_context(|| format!("finalizing Parquet file {}", self.path.display()))?;
        Ok(())
    }

    fn append_final_metadata(&mut self, counts: &Counters) -> Result<()> {
        self.append_metadata("max_width".to_owned(), self.max_width.to_string());
        append_counter_metadata(self, counts);
        Ok(())
    }
}

struct BatchBuilders {
    start: Int32Builder,
    width: Int32Builder,
    mapq: UInt8Builder,
    flag: UInt16Builder,
    nh: UInt8Builder,
    multi: BooleanBuilder,
    frag_id: Int32Builder,
    chrom_id: UInt8Builder,
}

impl BatchBuilders {
    fn new() -> Self {
        Self {
            start: Int32Builder::new(),
            width: Int32Builder::new(),
            mapq: UInt8Builder::new(),
            flag: UInt16Builder::new(),
            nh: UInt8Builder::new(),
            multi: BooleanBuilder::new(),
            frag_id: Int32Builder::new(),
            chrom_id: UInt8Builder::new(),
        }
    }

    fn len(&self) -> usize {
        self.start.len()
    }

    fn is_empty(&self) -> bool {
        self.start.is_empty()
    }

    fn append(&mut self, row: OutputRow, single_file: bool) {
        self.start.append_value(row.start);
        self.width.append_value(row.width);
        self.mapq.append_value(row.mapq);
        self.flag.append_value(row.flag);
        self.nh.append_value(row.nh);
        self.multi.append_value(row.multi);
        self.frag_id.append_value(row.frag_id);
        if single_file {
            self.chrom_id.append_value(row.chrom_id);
        }
    }

    fn finish(&mut self, schema: &Arc<Schema>, single_file: bool) -> Result<RecordBatch> {
        let mut columns: Vec<ArrayRef> = vec![
            Arc::new(self.start.finish()),
            Arc::new(self.width.finish()),
            Arc::new(self.mapq.finish()),
            Arc::new(self.flag.finish()),
            Arc::new(self.nh.finish()),
            Arc::new(self.multi.finish()),
            Arc::new(self.frag_id.finish()),
        ];
        if single_file {
            columns.push(Arc::new(self.chrom_id.finish()));
        }
        RecordBatch::try_new(schema.clone(), columns).context("building Arrow record batch")
    }
}

fn make_schema(single_file: bool) -> Schema {
    let mut fields = vec![
        Field::new("start", DataType::Int32, false),
        Field::new("width", DataType::Int32, false),
        Field::new("mapq", DataType::UInt8, false),
        Field::new("flag", DataType::UInt16, false),
        Field::new("nh", DataType::UInt8, false),
        Field::new("multi", DataType::Boolean, false),
        Field::new("frag_id", DataType::Int32, false),
    ];
    if single_file {
        fields.push(Field::new("chrom_id", DataType::UInt8, false));
    }
    Schema::new(fields)
}

fn writer_properties(
    row_group_size: usize,
    compression: CompressionSpec,
    single_file: bool,
) -> Result<WriterProperties> {
    let compression = match compression {
        CompressionSpec::Zstd(level) => {
            Compression::ZSTD(parquet::basic::ZstdLevel::try_new(level)?)
        }
        CompressionSpec::Snappy => Compression::SNAPPY,
        CompressionSpec::None => Compression::UNCOMPRESSED,
    };

    let mut sorting_columns = Vec::new();
    if single_file {
        sorting_columns.push(SortingColumn {
            column_idx: 7,
            descending: false,
            nulls_first: false,
        });
    }
    sorting_columns.push(SortingColumn {
        column_idx: 0,
        descending: false,
        nulls_first: false,
    });

    Ok(WriterProperties::builder()
        .set_compression(compression)
        .set_max_row_group_row_count(Some(row_group_size))
        .set_data_page_size_limit(DATA_PAGE_SIZE)
        .set_statistics_enabled(EnabledStatistics::Page)
        .set_column_dictionary_enabled(ColumnPath::from("start"), false)
        .set_column_dictionary_enabled(ColumnPath::from("width"), false)
        .set_column_encoding(ColumnPath::from("start"), Encoding::DELTA_BINARY_PACKED)
        .set_column_encoding(ColumnPath::from("width"), Encoding::DELTA_BINARY_PACKED)
        .set_sorting_columns(Some(sorting_columns))
        .build())
}

fn append_counter_metadata(sink: &mut ParquetSink, counts: &Counters) {
    for (key, value) in [
        ("n_records_read", counts.n_records_read),
        ("n_emitted", counts.n_emitted),
        ("n_pairs", counts.n_pairs),
        ("n_singletons", counts.n_singletons),
        ("n_secondary", counts.n_secondary),
        ("n_supplementary", counts.n_supplementary),
        ("n_duplicate", counts.n_duplicate),
        ("skip_unmapped", counts.skips.unmapped),
        ("skip_qc_fail", counts.skips.qc_fail),
        ("skip_min_mapq", counts.skips.min_mapq),
        ("skip_duplicate", counts.skips.duplicate),
        ("skip_secondary", counts.skips.secondary),
        ("skip_rightmost_mate", counts.skips.rightmost_mate),
    ] {
        sink.append_metadata(key.to_owned(), value.to_string());
    }
}

fn safe_file_component(value: &str) -> String {
    value
        .replace('%', "%25")
        .replace('/', "%2F")
        .replace('\\', "%5C")
}

pub fn validate_basename(basename: &str) -> Result<()> {
    let path = Path::new(basename);
    if basename.is_empty()
        || basename == "."
        || basename == ".."
        || path.components().count() != 1
        || path.file_name().and_then(|name| name.to_str()) != Some(basename)
    {
        bail!("basename must be a non-empty filename component");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chromosome_names_cannot_create_nested_output_paths() {
        assert_eq!(safe_file_component("chr1"), "chr1");
        assert_eq!(safe_file_component("chr/one"), "chr%2Fone");
    }

    #[test]
    fn rejects_basenames_that_escape_the_output_directory() {
        assert!(validate_basename("sample").is_ok());
        assert!(validate_basename("../outside").is_err());
        assert!(validate_basename("nested/name").is_err());
        assert!(validate_basename("").is_err());
    }
}
