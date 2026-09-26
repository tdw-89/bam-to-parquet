# bam2frag

`bam2frag` converts coordinate-sorted BAM alignments into compact Parquet fragment tables. It preserves coordinates, mapping quality, raw SAM flags, multi-mapper evidence, and fragment IDs while omitting sequences, base qualities, read names, and most auxiliary tags.

## Build and run

```sh
cargo build --release
cargo test
cargo fmt --check

cargo run --release -- --input sample.bam
```

Every option has a short form: `-i/--input`, `-o/--outdir`, `-b/--basename`, `-t/--threads`, `-1/--single-file`, `-r/--row-group-size`, `-c/--compression`, `-a/--aligner`, `-q/--min-mapq`, `-d/--drop-dups`, `-D/--drop-secondary`, `-S/--stats`.

## Output location

Only `--input` is required. Without `--outdir`, the tool creates a directory named after the input BAM, beside that BAM, and writes everything into it:

```
bam2frag -i /data/sample.bam

/data/sample.bam
/data/sample/sample.chr1.parquet
/data/sample/sample.chr2.parquet
/data/sample/sample.stats.json
```

An existing directory is reused, and files already in it are left alone unless this run overwrites them by name. Passing `--outdir` requires the directory to already exist — a missing one is an error rather than a silent `mkdir -p`, so a typo cannot scatter output into a new tree.

`--basename` names both the derived directory and the file prefix, so `-i /data/sample.bam -b H3K9me3_rep1` writes `/data/H3K9me3_rep1/H3K9me3_rep1.chr1.parquet`. It must be a single filename component, and is recorded in the Parquet footer metadata and the stats JSON under the key `basename`.

The default output is one Parquet file per reference sequence (`sample.chr1.parquet`, etc.) plus `sample.stats.json`; `--stats` moves the stats file elsewhere. `--single-file` writes `<basename>.parquet` and includes `chrom_id` values from the BAM header’s reference order; it supports up to 256 reference sequences.

Compression defaults to `zstd:3`; `--compression snappy` and `--compression none` are also supported. The default row group size is 1,000,000 rows, and the default BGZF worker count is four.

## Aligner and filtering behavior

`--aligner auto` examines BAM `@PG` records and recognizes Bowtie 2, BWA-MEM, STAR, and HISAT2. If no single supported aligner can be identified, the tool prints a warning and uses the conservative `MAPQ == 0` multi-mapper rule. Pass `--aligner bowtie2`, `bwa`, `star`, or `hisat2` to override detection. The selected rule is logged and stored in Parquet metadata.

Multi-mapper calls follow the aligner rule: Bowtie 2 uses `XS`, BWA-MEM uses `XA` or MAPQ 0, STAR and HISAT2 use `NH > 1`, and the unknown fallback uses MAPQ 0. `SA` is not treated as multi-mapper evidence. Raw `NH` is retained as `nh`, with 0 for a missing tag and saturation at 255.

Unmapped and QC-failed records are always skipped. Duplicates, secondary alignments, and supplementary alignments are retained and marked by default; `--drop-dups` and `--drop-secondary` opt into dropping the first two classes. Without `--min-mapq`, no MAPQ filtering is applied. When supplied, `--min-mapq` must be 1–254 and records below it are dropped; 255 means MAPQ is unavailable and is retained.

`frag_id` is always present. It groups emitted alignments by query name, read group when present, and first/last segment flag. Names are used only to assign IDs and are not written to Parquet. IDs for named reads are kept in memory during conversion so secondary placements can share an ID; records with no usable query name receive unique IDs.

## Fragment rows

All coordinates use 0-based, half-open intervals. `start` is the leftmost coordinate and `width` is the fragment span; `end = start + width`.

- A properly paired record with positive `TLEN` emits one fragment row at that leftmost mate, with `width = TLEN`.
- Its properly paired mate with negative `TLEN` is skipped because the fragment was already emitted.
- A proper pair with `TLEN == 0`, an unpaired record, or a non-proper pair emits one row per alignment. Width is the CIGAR reference span (`M`, `D`, `N`, `=`, and `X` consume reference).

The raw `flag` remains unchanged, so downstream code can distinguish proper-pair fragments from single-alignment rows.

The non-nullable Parquet columns are `start` (`Int32`), `width` (`Int32`), `mapq` (`UInt8`), `flag` (`UInt16`), `nh` (`UInt8`), `multi` (`Boolean`), and `frag_id` (`Int32`). Single-file mode adds `chrom_id` (`UInt8`). Each output includes provenance, reference details, the multi-mapper rule, `max_width`, and counters in footer metadata. The JSON stats file includes global and per-reference counts plus a 1 bp fragment-length histogram for proper pairs.

## Range queries

Use the file’s `max_width` when finding fragments overlapping a query window. Filtering only on `start >= qstart` misses fragments that begin before the window. For a window `[200, 210)` and `max_width = 150`, scan starts from 50 through 210, then apply the overlap check. A fragment at `start = 100`, `width = 150` ends at 250 and overlaps the window.

For a per-reference file, the DuckDB predicate is:

```sql
WHERE start BETWEEN greatest(0, qstart - max_width) AND qend
  AND start + width > qstart
  AND start < qend
```

For single-file output, also constrain `chrom_id`; rows are ordered by `chrom_id`, then `start`.

## Pre-flight before deleting a BAM

Inspect aligner metadata and whether secondary alignments and duplicate marking are present:

```sh
samtools view -c -f 256 sample.bam   # secondary alignments
samtools view -c -f 1024 sample.bam  # duplicate-marked alignments
samtools view -H sample.bam | grep '@PG'
```

If the secondary count is zero, the aligner may have reported only one location for each multi-mapping read. In that case, the Parquet output preserves the reported location and multi-mapper evidence, but cannot recover alternative locations after the BAM is deleted.
