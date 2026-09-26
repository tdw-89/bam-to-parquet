use noodles_bam as bam;
use noodles_sam::alignment::record::data::field::Tag;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum Aligner {
    #[default]
    Auto,
    Bowtie2,
    Bwa,
    Star,
    Hisat2,
    Unknown,
}

impl Aligner {
    pub fn name(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Bowtie2 => "bowtie2",
            Self::Bwa => "bwa",
            Self::Star => "star",
            Self::Hisat2 => "hisat2",
            Self::Unknown => "unknown",
        }
    }

    pub fn rule(self) -> &'static str {
        match self {
            Self::Bowtie2 => "bowtie2: XS tag present",
            Self::Bwa => "bwa mem: XA tag present or MAPQ == 0",
            Self::Star | Self::Hisat2 => "NH tag > 1",
            Self::Auto | Self::Unknown => "unknown aligner: MAPQ == 0",
        }
    }
}

pub fn detect_aligner(program_lines: &[String]) -> Option<Aligner> {
    let mut detected = None;
    for aligner in program_lines
        .iter()
        .filter_map(|line| identify_program_line(line))
    {
        if detected.is_some_and(|previous| previous != aligner) {
            return None;
        }
        detected = Some(aligner);
    }
    detected
}

fn identify_program_line(line: &str) -> Option<Aligner> {
    let lower = line.to_ascii_lowercase();
    if lower.contains("bowtie2") {
        Some(Aligner::Bowtie2)
    } else if lower.contains("hisat2") {
        Some(Aligner::Hisat2)
    } else if lower.contains("star") {
        Some(Aligner::Star)
    } else if lower.contains("bwa mem") || lower.contains("bwa-mem") {
        Some(Aligner::Bwa)
    } else {
        None
    }
}

pub fn is_multi(record: &bam::Record, aligner: Aligner, mapq: u8, nh: u8) -> anyhow::Result<bool> {
    let (xs, xa) = match aligner {
        Aligner::Bowtie2 => (has_tag(record, Tag::new(b'X', b'S'))?, false),
        Aligner::Bwa if mapq != 0 => (false, has_tag(record, Tag::new(b'X', b'A'))?),
        _ => (false, false),
    };
    Ok(multi_call(aligner, mapq, nh, xs, xa))
}

fn multi_call(aligner: Aligner, mapq: u8, nh: u8, xs: bool, xa: bool) -> bool {
    match aligner {
        Aligner::Bowtie2 => xs,
        Aligner::Bwa => xa || mapq == 0,
        Aligner::Star | Aligner::Hisat2 => nh > 1,
        Aligner::Auto | Aligner::Unknown => mapq == 0,
    }
}

pub fn nh_value(record: &bam::Record) -> anyhow::Result<u8> {
    let Some(value) = record.data().get(&Tag::new(b'N', b'H')).transpose()? else {
        return Ok(0);
    };

    let Some(nh) = value.as_int() else {
        return Ok(0);
    };

    Ok(saturate_nh(nh))
}

fn saturate_nh(nh: i64) -> u8 {
    nh.clamp(0, i64::from(u8::MAX)) as u8
}

fn has_tag(record: &bam::Record, tag: Tag) -> anyhow::Result<bool> {
    Ok(record.data().get(&tag).transpose()?.is_some())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detects_supported_aligners_from_program_headers() {
        let cases = [
            ("@PG\tID:bowtie2\tPN:bowtie2\tVN:2.5", Aligner::Bowtie2),
            ("@PG\tID:bwa\tPN:bwa\tCL:bwa mem ref", Aligner::Bwa),
            ("@PG\tID:STAR\tPN:STAR\tVN:2.7", Aligner::Star),
            ("@PG\tID:hisat2\tPN:hisat2\tVN:2.2", Aligner::Hisat2),
        ];

        for (line, expected) in cases {
            assert_eq!(identify_program_line(line), Some(expected));
        }
    }

    #[test]
    fn last_recognized_program_is_selected() {
        let lines = vec![
            "@PG\tID:bowtie2\tPN:bowtie2".to_owned(),
            "@PG\tID:samtools\tPN:samtools".to_owned(),
        ];
        assert_eq!(detect_aligner(&lines), Some(Aligner::Bowtie2));
        assert_eq!(detect_aligner(&[]), None);
    }

    #[test]
    fn conflicting_aligner_headers_request_fallback() {
        let lines = vec![
            "@PG\tID:bwa\tPN:bwa\tCL:bwa mem ref".to_owned(),
            "@PG\tID:bowtie2\tPN:bowtie2".to_owned(),
        ];
        assert_eq!(detect_aligner(&lines), None);
    }

    #[test]
    fn reports_rules_for_each_aligner() {
        assert_eq!(Aligner::Bowtie2.rule(), "bowtie2: XS tag present");
        assert_eq!(Aligner::Bwa.rule(), "bwa mem: XA tag present or MAPQ == 0");
        assert_eq!(Aligner::Star.rule(), "NH tag > 1");
        assert_eq!(Aligner::Hisat2.rule(), "NH tag > 1");
        assert_eq!(Aligner::Unknown.rule(), "unknown aligner: MAPQ == 0");
    }

    #[test]
    fn applies_each_multi_mapper_rule() {
        assert!(multi_call(Aligner::Bowtie2, 60, 1, true, false));
        assert!(!multi_call(Aligner::Bowtie2, 0, 4, false, true));
        assert!(multi_call(Aligner::Bwa, 60, 1, false, true));
        assert!(multi_call(Aligner::Bwa, 0, 1, false, false));
        assert!(multi_call(Aligner::Star, 60, 2, false, false));
        assert!(multi_call(Aligner::Hisat2, 60, 2, false, false));
        assert!(multi_call(Aligner::Unknown, 0, 1, false, false));
        assert!(!multi_call(Aligner::Unknown, 60, 4, false, false));
    }

    #[test]
    fn saturates_raw_nh_without_losing_multiple_mappings() {
        assert_eq!(saturate_nh(1), 1);
        assert_eq!(saturate_nh(255), 255);
        assert_eq!(saturate_nh(300), 255);
        assert_eq!(saturate_nh(-1), 0);
    }
}
