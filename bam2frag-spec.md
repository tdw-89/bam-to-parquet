# bam2frag — BAM → Parquet fragment table

Implementation spec. Written to be handed to an agent working in a fresh repo.

---

## 1. Purpose

Convert coordinate-sorted BAM files (histone ChIP-seq, including input controls) into
compact, sorted, range-queryable Parquet "fragment tables" for downstream gene-level
signal analysis in Julia.

The motivating constraint: BAM files are pulled one at a time to a laptop with limited
disk, converted, and the BAM deleted. The output must be small enough to keep all of
them, and complete enough that the BAM never needs to be re-downloaded.

**Preserve:** alignment position, fragment span, MAPQ, SAM flags, multi-mapper evidence.
**Discard:** read sequence, base qualities, read names, CIGAR detail, most aux tags.

The marks of interest (H3K9me2/H3K9me3) sit largely in repeat-rich territory, and the
analysis compares *genes against each other* rather than conditions at one locus. That
makes multi-mapper status and MAPQ first-class data, not noise to be filtered at
conversion time. **The tool compresses; it does not make analysis decisions.** Every
filter defaults to permissive.

## 2. Non-goals

- No peak calling, normalization, binning, or bigWig output.
- No re-sorting. Coordinate-sorted input is a precondition, verified and enforced.
- No CRAM/SAM input in v1 (noodles makes it easy to add later; keep the reader behind a
  small trait boundary so it can be).
- No mate buffering or name-based pairing. See §5.

## 3. CLI

```
bam2frag --input <sample.bam> [options]

  -i, --input <path>          Coordinate-sorted BAM. Required.
  -o, --outdir <dir>          Output directory. Must already exist; a missing
                              directory is an error. Default: <basename>/ created
                              beside the input BAM, reused if already present.
  -b, --basename <name>       Names both the derived output directory and the file
                              prefix. Must be one filename component.
                              Default: input file stem.
  -t, --threads <n>           BGZF decompression workers. Default: 4.
  -1, --single-file           One Parquet file with a chrom_id column instead of
                              per-chromosome files. Default: per-chromosome.
  -r, --row-group-size <n>    Default: 1_000_000.
  -c, --compression <spec>    Default: zstd:3. Also: snappy, none.
  -a, --aligner <auto|bowtie2|bwa|star|hisat2>
                              Multi-mapper detection rule. Default: auto (from @PG).
  -q, --min-mapq <n>          Drop records below this MAPQ; valid scores are 1-254.
                              Default: unset (no MAPQ filtering).
  -d, --drop-dups             Drop flag 0x400. Default: keep, marked.
  -D, --drop-secondary        Drop flag 0x100. Default: keep, marked.
  -S, --stats <path.json>     Default: <outdir>/<basename>.stats.json.
```

`frag_id` is always emitted; there is no `--emit-frag-id` flag.

Exit non-zero with a clear message if the header is not `@HD SO:coordinate`.

## 4. Output layout and schema

Per-chromosome by default: `<outdir>/<basename>.<refname>.parquet`.

Per-chromosome files are the primary layout because downstream work iterates
annotations chromosome by chromosome. It gives free partition pruning, trivial
parallelism, and drops a column.

### Columns (all non-nullable)

| column    | Arrow type | notes |
|-----------|-----------|-------|
| `start`   | `Int32`   | 0-based, leftmost reference coordinate of the fragment |
| `width`   | `Int32`   | fragment span; `end = start + width`, half-open |
| `mapq`    | `UInt8`   | as reported; 255 = unavailable |
| `flag`    | `UInt16`  | raw SAM flag, unmodified |
| `nh`      | `UInt8`   | raw `NH` tag where present; 0 = aligner did not report; saturate at 255 |
| `multi`   | `Boolean` | derived multi-mapper call, see §6 |
| `frag_id` | `Int32`   | optional; groups alignments of one read when secondaries are present |
| `chrom_id`| `UInt8`   | only in `--single-file` mode; index into the metadata chrom table |

Note on widths: Parquet's physical type for anything ≤32 bits is `INT32` regardless, so
narrow Arrow types buy nothing on disk here — encoding does the work. `Int32` for
`start`/`width` avoids a class of overflow bugs for free. Keep `mapq`/`flag`/`nh` narrow
anyway so that materializing a chromosome to uncompressed Arrow later stays cheap.

`multi` is redundant with `nh` plus `flag` for some aligners and not for others. Store
both: `nh` is raw evidence, `multi` is the tool's judgment, and the rule behind the
judgment is recorded in metadata so it can be recomputed or overridden.

### Sorting

Rows within each file MUST be non-decreasing in `start`. This is what makes row-group
statistics useful, and it is the single invariant most likely to be broken by a
well-meaning refactor. Assert it in the writer, not just in tests.

Set Parquet's `SortingColumn` metadata for `start` ascending.

## 5. Record → row derivation

Coordinate-sorted input means mates are far apart in the file. Do **not** buffer mates.
Instead, emit each fragment once, at its leftmost mate, using `TLEN`:

**Paired, proper pair (`0x1` and `0x2` set), `TLEN > 0`**
→ emit one row. `start = alignment_start`, `width = TLEN`.
This read is the leftmost mate, so emission order follows the BAM's sort order.

**Paired, proper pair, `TLEN < 0`**
→ skip. Already emitted by the leftmost mate.

**Paired, proper pair, `TLEN == 0`**
→ degenerate (mate unmapped, or single-reference-position edge case). Fall through to
the single-read rule and count it in stats.

**Not a proper pair (discordant, mate on another chromosome), or unpaired**
→ emit one row per alignment. `width` = reference span computed from the CIGAR
(`M`, `D`, `N`, `=`, `X` consume reference; `I`, `S`, `H`, `P` do not).

Downstream can separate these cases using the `PROPER_PAIR` bit already in `flag`; no
extra column is needed. Document this in the README.

### Filtering

| condition | default |
|-----------|---------|
| unmapped (`0x4`) | always skipped |
| QC-fail (`0x200`) | skipped |
| secondary (`0x100`) | **kept**, marked in `flag` |
| supplementary (`0x800`) | **kept**, marked in `flag` |
| duplicate (`0x400`) | **kept**, marked in `flag` |
| `mapq < --min-mapq` | kept (default threshold 0) |

Keeping duplicates marked rather than removed preserves the ability to change that
decision later without the BAM.

## 6. Multi-mapper detection

Auto-detect the aligner from `@PG` lines in the header; `--aligner` overrides.

| aligner | rule for `multi = true` |
|---------|------------------------|
| bowtie2 | `XS:i` present (a second-best alignment exists) |
| bwa mem | `XA:Z` present, or `mapq == 0` |
| STAR    | `NH:i > 1` |
| hisat2  | `NH:i > 1` |
| unknown | `mapq == 0` |

Record the rule actually applied in file metadata, and log it at startup. If `--aligner
auto` cannot identify the aligner, warn loudly rather than silently falling back —
getting this wrong quietly biases exactly the comparison the data is for.

`SA:Z` indicates a supplementary/chimeric alignment, **not** multi-mapping. Do not
conflate them.

## 7. Parquet writer configuration

Use `arrow-rs` (`arrow` + `parquet` crates, pinned to the same version).

```
WriterProperties:
  compression                 ZSTD(level 3)
  max_row_group_size          1_000_000          (not the 128 MB default)
  data_page_size_limit        ~1 MB
  statistics_enabled          EnabledStatistics::Page   (emits ColumnIndex/OffsetIndex)
  dictionary_enabled          false for start/width
  column encoding: start      DELTA_BINARY_PACKED
  column encoding: width      DELTA_BINARY_PACKED or RLE — benchmark both
  sorting_columns             [start ASC]
```

Delta encoding on the sorted `start` column is where most of the size win comes from.
Page-level statistics are what let DuckDB and other readers skip below row-group
granularity.

Build batches of ~65_536 rows in Arrow builders, then hand to the writer.

## 8. Footer key-value metadata

Written into every output file. This is the provenance that makes deleting the BAM
safe, so treat it as required output, not decoration.

- `bam2frag_version`, `git_sha`, `schema_version`, `created_utc`
- `basename`, `source_bam` (path, size in bytes, mtime)
- `bam_header_pg` — the raw `@PG` lines (aligner, version, full command line)
- `reference_name`, `reference_length`, `reference_md5` (from `@SQ` `M5` if present)
- `multi_rule` — the §6 rule applied
- **`max_width`** — the largest `width` in this file
- counts: `n_emitted`, `n_pairs`, `n_singletons`, `n_secondary`, `n_supplementary`,
  `n_duplicate`, and each skip reason
- `chrom_table` — only in `--single-file` mode, mapping `chrom_id` → name

`max_width` is load-bearing. A range query for `[qstart, qend)` must scan
`start BETWEEN qstart - max_width AND qend`, then filter on `start + width > qstart`.
Without the widening, fragments that begin before the window and overlap into it are
silently lost. Call this out in the README with a worked example.

## 9. Architecture

```
src/
  main.rs         wiring, error reporting
  cli.rs          clap derive
  header.rs       sort-order check, @PG parsing, aligner detection, chrom table
  reader.rs       BGZF + BAM setup, record loop
  fragment.rs     §5 derivation — the correctness core, heavily unit-tested
  multimap.rs     §6 rules
  writer.rs       Arrow builders, Parquet writer, per-chromosome rotation
  stats.rs        counters, fragment-length histogram, JSON output
tests/
  data/tiny.bam   few hundred records, committed; PE + SE + secondary + dup cases
  integration.rs
```

### Reading

```rust
let decoder = bgzf::io::MultithreadedReader::with_worker_count(n, file);
let mut reader = bam::io::Reader::from(decoder);
let header = reader.read_header()?;
```

`with_worker_count` was deprecated around noodles-bgzf 0.48 in favor of configuring a
`rayon::ThreadPoolBuilder`. Check the API for the pinned version and use whichever form
is current. Enable the `libdeflate` feature if the pinned version exposes it.

Performance rules for the record loop, in priority order:

1. **Keep the loop sequential and ordered.** Do not use `records().par_bridge()`. It
   destroys the sort invariant of §4, which destroys the pruning of §7. Decompression
   parallelism is fine — the reader reassembles blocks in file order.
2. **Reuse one record:** `while reader.read_record(&mut rec)? != 0`, not `records()`.
   At ~30M records the per-record allocation is real.
3. **Stay lazy.** `bam::Record` is backed by the raw buffer. Read only `flags()`,
   `reference_sequence_id()`, `alignment_start()`, `mapping_quality()`,
   `template_length()`, the CIGAR when needed, and specific tags via a targeted `data()`
   lookup. Never convert to `sam::Record`, never walk all of `data()`, never touch
   sequence or quality — those are the bulk of the bytes and the whole point is to skip
   them.

Worker count plateaus early (~4). Past that the single consumer thread becomes the
bottleneck. On Apple Silicon do not use `available_parallelism()`; it counts efficiency
cores, which are poor at inflate. Default to 4.

### Writing

Optional but recommended: run the Parquet writer on its own thread behind a
`sync_channel(2)` of `RecordBatch`, so ZSTD overlaps parsing. Keep it simple — this is
not where the time goes.

Chromosome rotation: because input is coordinate-sorted, records for a reference are
contiguous. On `reference_sequence_id` change, finalize `max_width` and the counters for
the current file, write metadata, close, open the next.

## 10. Stats output

JSON alongside the Parquet. Include every counter from §8 plus a **fragment length
histogram** (proper pairs only, 1 bp bins up to 1000, then a tail bucket). That
histogram is a genuine ChIP QC signal, essentially free to compute here, and saves a
separate pass over data that will have been deleted.

Invariant to assert and report: `n_emitted + sum(skip_reasons) == n_records_read`. Any
mismatch means a case in §5 is unhandled.

## 11. Testing

- **Unit:** TLEN-based fragment derivation across all §5 branches including the `TLEN == 0`
  and discordant cases; CIGAR reference span for every operator; each §6 aligner rule.
- **Invariant:** `start` non-decreasing in every output file — assert in the writer and
  test with a shuffled-input fixture that it fails.
- **Round-trip:** convert `tiny.bam`, read back with arrow-rs, assert exact row count and
  spot-check known coordinates.
- **Cross-check:** `bedtools bamtobed -bedpe` on `tiny.bam`, compare intervals for proper
  pairs. Any disagreement is a bug in §5, not in bedtools.
- **Range-query correctness:** pick a locus whose fragments start before the window and
  overlap into it; assert the `max_width` widening rule recovers exactly them.

## 12. Milestones

Each is a commit with tests passing.

1. CLI skeleton, open BAM, read header, verify sort order, count records, exit.
2. §5 derivation + §10 stats, no Parquet yet. Validate against bedtools here — get the
   biology right before touching the output format.
3. Parquet writer, single file, full §8 metadata.
4. Per-chromosome rotation, `max_width` per file.
5. Multithreaded decode, optional writer thread, benchmark.
6. README with the Julia/DuckDB query contract from §8, plus the §13 pre-flight.

## 13. Pre-flight for the user (document in README)

Before converting and deleting a BAM, check what is actually in it:

```
samtools view -c -f 256 sample.bam    # secondary alignments present?
samtools view -c -f 1024 sample.bam   # duplicates marked, or already removed?
samtools view -H sample.bam | grep '@PG'
```

If the secondary count is 0, the aligner reported one location per multi-mapper (bowtie2's
default). The file then records only *that* a read maps elsewhere, not where — and the
reported position is an arbitrary pick among equals. This cannot be recovered after the
BAM is gone, and it rules out EM-style multi-mapper reallocation without realigning.
Worth knowing before deletion, not after.

## 14. Dependencies

```
noodles-bam, noodles-sam, noodles-bgzf   (libdeflate feature if available)
arrow, parquet                            (same pinned version)
clap (derive), anyhow, serde, serde_json
criterion (dev), optional indicatif
```

## 15. Open questions for the user

Two things the implementer should confirm rather than guess:

1. **Aligner and paired-end status** of the actual BAMs. Auto-detection covers the common
   cases, but the §6 rule is biasing if wrong.
2. **Whether secondary alignments exist** (§13). If they do, `frag_id` is required and
   grouping semantics need a decision; if they don't, drop the column entirely.
