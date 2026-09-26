use serde::Serialize;

#[derive(Clone, Copy, Debug)]
pub enum SkipReason {
    Unmapped,
    QcFail,
    MinMapq,
    Duplicate,
    Secondary,
    RightmostMate,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct SkipCounts {
    pub unmapped: u64,
    pub qc_fail: u64,
    pub min_mapq: u64,
    pub duplicate: u64,
    pub secondary: u64,
    pub rightmost_mate: u64,
}

impl SkipCounts {
    pub fn increment(&mut self, reason: SkipReason) {
        let count = match reason {
            SkipReason::Unmapped => &mut self.unmapped,
            SkipReason::QcFail => &mut self.qc_fail,
            SkipReason::MinMapq => &mut self.min_mapq,
            SkipReason::Duplicate => &mut self.duplicate,
            SkipReason::Secondary => &mut self.secondary,
            SkipReason::RightmostMate => &mut self.rightmost_mate,
        };
        *count += 1;
    }

    pub fn sum(&self) -> u64 {
        self.unmapped
            + self.qc_fail
            + self.min_mapq
            + self.duplicate
            + self.secondary
            + self.rightmost_mate
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct FragmentLengthHistogram {
    /// Index 0 is width 1 bp; index 999 is width 1000 bp.
    pub bins_1_to_1000: Vec<u64>,
    pub over_1000: u64,
}

impl Default for FragmentLengthHistogram {
    fn default() -> Self {
        Self {
            bins_1_to_1000: vec![0; 1000],
            over_1000: 0,
        }
    }
}

impl FragmentLengthHistogram {
    pub fn add(&mut self, width: i32) {
        if (1..=1000).contains(&width) {
            self.bins_1_to_1000[width as usize - 1] += 1;
        } else if width > 1000 {
            self.over_1000 += 1;
        }
    }
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct Counters {
    pub n_records_read: u64,
    pub n_emitted: u64,
    pub n_pairs: u64,
    pub n_singletons: u64,
    /// These flag counts describe emitted rows; dropped rows appear in `skips`.
    pub n_secondary: u64,
    pub n_supplementary: u64,
    pub n_duplicate: u64,
    pub skips: SkipCounts,
    pub fragment_length_histogram: FragmentLengthHistogram,
}

impl Counters {
    pub fn record_read(&mut self) {
        self.n_records_read += 1;
    }

    pub fn skip(&mut self, reason: SkipReason) {
        self.skips.increment(reason);
    }

    pub fn emit(&mut self, pair_width: Option<i32>, flags: u16) {
        self.n_emitted += 1;
        if let Some(width) = pair_width {
            self.n_pairs += 1;
            self.fragment_length_histogram.add(width);
        } else {
            self.n_singletons += 1;
        }
        if flags & 0x0100 != 0 {
            self.n_secondary += 1;
        }
        if flags & 0x0800 != 0 {
            self.n_supplementary += 1;
        }
        if flags & 0x0400 != 0 {
            self.n_duplicate += 1;
        }
    }

    pub fn invariant_holds(&self) -> bool {
        self.n_emitted + self.skips.sum() == self.n_records_read
    }
}

#[derive(Debug, Serialize)]
pub struct ReferenceStats {
    pub reference: crate::header::ReferenceSequence,
    #[serde(flatten)]
    pub counts: Counters,
}

#[derive(Debug, Serialize)]
pub struct StatsReport {
    pub bam2frag_version: &'static str,
    pub schema_version: &'static str,
    pub created_utc: String,
    pub basename: String,
    pub source_bam: String,
    pub source_bam_size_bytes: u64,
    pub source_bam_mtime: String,
    pub aligner: String,
    pub multi_rule: String,
    pub min_mapq: Option<u8>,
    pub drop_dups: bool,
    pub drop_secondary: bool,
    pub single_file: bool,
    pub overall: Counters,
    pub references: Vec<ReferenceStats>,
}
