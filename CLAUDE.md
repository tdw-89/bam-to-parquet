# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```sh
cargo build --release
cargo test                                    # 16 unit + 2 integration tests
cargo test fragment::tests                    # one module
cargo test derives_each_template_length_branch # one test by name
cargo test --test integration                 # integration tests only
cargo fmt --check                             # cargo fmt to apply
cargo clippy -- -D warnings

cargo run --release -- --input sample.bam   # or -i; writes ./sample/ beside the BAM
```

Tests need no fixture files: they synthesize BAMs in a `tempfile::tempdir()` with
`noodles_sam::alignment::RecordBuf` and read the Parquet back with arrow-rs. Do not commit BAMs.

`build.rs` injects `GIT_SHA` (via `option_env!`) into Parquet footer metadata; it is `"unknown"`
outside a git checkout.

## Documents that govern the code

- `bam2frag-spec.md` — the implementation contract: schema, record→row derivation (§5),
  multi-mapper rules (§6), writer configuration (§7), required footer metadata (§8). Changes to
  behavior belong here first.
- `README.md` — the user-facing contract, including the `max_width` range-query recipe that
  downstream Julia/DuckDB code depends on.
- `AGENTS.md` — repo conventions (commit style, test naming).

Implementation has diverged from the spec in two places, and the README is authoritative:
`frag_id` is always emitted (no `--emit-frag-id` flag), and `--min-mapq` is an `Option<u8>`
restricted to 1–254 rather than a default of 0.

`--basename` (short `-b`) sets the file prefix *and* the derived output directory name, and is
emitted under the key `basename` in both the Parquet footer metadata and the stats JSON. It was
called `sample` everywhere until the output contract was renamed to match the flag; nothing
consumed the old key, so `schema_version` stayed at `1`. Renaming an emitted key again means
bumping that.

The crate is `bam-to-parquet`; the CLI presents itself as `bam2frag`.

## Architecture

`src/reader.rs::convert` is the whole pipeline; everything else is a leaf module it calls.
`lib.rs::run` only parses args and hands off.

`resolve_outdir` decides where output goes before anything is read: an explicit `--outdir` must
already exist (missing is an error, never an implicit `mkdir -p`), while the default derives
`<basename>/` beside the input BAM and creates it, reusing an existing one.

The header is read **twice**, deliberately:

1. `header::read_raw_header` + `header::parse` pull the raw SAM text to check `@HD SO:coordinate`,
   collect `@SQ` (name/length/M5) and the raw `@PG` lines. Raw text is kept because `@PG` lines go
   verbatim into footer metadata as provenance.
2. `bam::io::Reader::read_header` during record decoding. The two reference dictionaries are
   cross-checked for length and the run aborts if they disagree.

Aligner choice (`multimap::detect_aligner` over `@PG`, or `--aligner`) picks the multi-mapper rule
before any record is read; the rule string is logged and stored in metadata. Ambiguous or
unrecognized `@PG` falls back to `Aligner::Unknown` (MAPQ == 0) with a loud warning — a silent
fallback would bias the exact comparison the data is for.

`writer.rs::OutputWriter` owns per-reference file rotation. Because input is coordinate-sorted,
records for a reference are contiguous: the reader calls `finish_reference` on a
`reference_sequence_id` change, and `OutputWriter::write_row` refuses to write to a reference other
than the open one. After the loop, `finish_empty_references` emits a zero-row Parquet file (with
full metadata) for every reference that had no emitted rows, so the per-chromosome file set always
matches the header. `--single-file` uses one `ParquetSink` instead, plus a `chrom_id` column and a
`chrom_table` metadata entry; it is capped at 256 references because `chrom_id` is `UInt8`.

`frag_id` grouping lives in `reader.rs`: a `HashMap<ReadIdentity, i32>` keyed by
(query name, RG tag, first/last-segment flag). Records with no usable name get an anonymous key
from the record ordinal so they never collide. This map grows with the number of distinct reads —
it is the one unbounded allocation in the conversion.

## Invariants asserted at runtime, not just in tests

Breaking any of these should fail the run, not just a test. Keep them that way.

- **Coordinate sort order** — the record loop rejects any `(reference_id, start)` that regresses,
  and rejects mapped or placed records appearing after unplaced unmapped records.
- **Non-decreasing `start` per output file** — `ParquetSink::write` compares against `last_key`
  (`(chrom_id, start)` in single-file mode). This is what makes row-group/page statistics prune,
  so it is checked in the writer, at the point of truth.
- **Record accounting** — `n_emitted + skips.sum() == n_records_read`, checked before outputs are
  finalized. A mismatch means a §5 branch is unhandled.
- **Path safety** — `validate_basename` and `safe_file_component` keep `--basename` and reference
  names from escaping the output directory.

## Record loop performance rules

These are load-bearing; the point of the tool is to skip the bulk of the BAM's bytes.

- Reuse a single `bam::Record` with `read_record(&mut record)`; never iterate `records()`.
- Stay lazy on the raw buffer: read only flags, reference id, start, MAPQ, TLEN, the CIGAR when a
  singleton needs a reference span, and targeted `data().get(&Tag::new(..))` lookups for NH/XS/XA/RG.
  Never convert to `sam::Record`, never walk all of `data()`, never touch sequence or qualities.
- Keep the loop sequential. BGZF decompression parallelism is fine (`MultithreadedReader`, default
  4 workers — on Apple Silicon do not use `available_parallelism()`, efficiency cores inflate
  poorly); `par_bridge()` over records would destroy the sort invariant above.

## Output-facing details worth preserving

- `max_width` in footer metadata is load-bearing for correct range queries. Any change to how
  `width` is derived must keep it accurate per file.
- `nh` is raw evidence (0 = tag absent, saturated at 255) and `multi` is the tool's judgment; both
  are stored, with the rule recorded, so `multi` can be recomputed later.
- `flag` is the unmodified SAM flag: downstream separates proper-pair fragments from single-alignment
  rows with the `PROPER_PAIR` bit, which is why no extra column exists for it.
- Filtering defaults are permissive by design (duplicates, secondary and supplementary alignments
  are kept and marked). Do not add implicit filtering.
