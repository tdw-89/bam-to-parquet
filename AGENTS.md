# Repository Guidelines

## Project Structure & Module Organization

This repository is a Rust CLI project. Keep application code in `src/` (currently `src/main.rs`); place integration tests in `tests/` and unit tests beside the code they exercise. `README.md` introduces the tool, while `bam2frag-spec.md` defines the intended BAM-to-Parquet behavior, schema, and invariants. Keep implementation and documentation aligned with that specification.

## Build, Test, and Development Commands

- `cargo build` compiles the project.
- `cargo run -- <arguments>` builds and runs the CLI; `src/main.rs` is currently a scaffold, so add supported options as the CLI is implemented.
- `cargo test` runs unit and integration tests. Add tests for record conversion, filtering, sorting, and output schema as those features land.
- `cargo fmt --check` checks standard Rust formatting; run `cargo fmt` to apply it.
- `cargo clippy -- -D warnings` checks for common Rust issues and treats warnings as errors.

## Coding Style & Naming Conventions

Use standard `rustfmt` formatting and idiomatic Rust: `snake_case` for functions, variables, and modules; `UpperCamelCase` for types and traits; and `SCREAMING_SNAKE_CASE` for constants. Keep BAM parsing, record-to-row conversion, and Parquet writing in focused modules as the implementation grows. Preserve the spec's coordinate sorting and raw alignment evidence; avoid introducing implicit filtering defaults.

## Testing Guidelines

Use Rust's built-in test harness. Name tests for the behavior they verify (for example, `skips_rightmost_proper_pair`), and cover boundary cases in fragment widths, flags, MAPQ, and multi-mapper tags. Run `cargo test` before submitting changes; include small fixtures where practical and avoid committing large BAM files.

## Commit & Pull Request Guidelines

The available Git history contains only an initial commit, so no established commit format is visible. Use short imperative subjects such as `Validate coordinate-sorted BAM input`. Pull requests should explain the behavior change, note relevant specification or schema updates, link related issues when applicable, and report the checks run. Include representative CLI output for user-visible changes.
