fn main() {
    if let Err(error) = bam_to_parquet::run() {
        eprintln!("error: {error:#}");
        std::process::exit(1);
    }
}
