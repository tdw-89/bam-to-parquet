use noodles_sam::alignment::record::{Flags, cigar::op::Kind};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Derivation {
    ProperPair { width: i32 },
    SkipRightmostMate,
    Singleton,
}

pub fn derive(flags: Flags, template_length: i32) -> Derivation {
    if flags.is_segmented() && flags.is_properly_segmented() {
        match template_length.cmp(&0) {
            std::cmp::Ordering::Greater => Derivation::ProperPair {
                width: template_length,
            },
            std::cmp::Ordering::Less => Derivation::SkipRightmostMate,
            std::cmp::Ordering::Equal => Derivation::Singleton,
        }
    } else {
        Derivation::Singleton
    }
}

pub fn cigar_reference_span(ops: impl IntoIterator<Item = (Kind, usize)>) -> anyhow::Result<i32> {
    let span = ops
        .into_iter()
        .filter(|(kind, _)| kind.consumes_reference())
        .try_fold(0usize, |sum, (_, len)| sum.checked_add(len))
        .ok_or_else(|| anyhow::anyhow!("CIGAR reference span overflow"))?;

    i32::try_from(span).map_err(|_| anyhow::anyhow!("CIGAR reference span exceeds Int32"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derives_each_template_length_branch() {
        let proper = Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED;
        assert_eq!(derive(proper, 120), Derivation::ProperPair { width: 120 });
        assert_eq!(derive(proper, -120), Derivation::SkipRightmostMate);
        assert_eq!(derive(proper, 0), Derivation::Singleton);
        assert_eq!(derive(Flags::SEGMENTED, 120), Derivation::Singleton);
        assert_eq!(derive(Flags::empty(), 120), Derivation::Singleton);
    }

    #[test]
    fn cigar_span_counts_only_reference_consuming_operations() {
        let ops = [
            (Kind::Match, 10),
            (Kind::Insertion, 2),
            (Kind::Deletion, 3),
            (Kind::Skip, 4),
            (Kind::SoftClip, 5),
            (Kind::HardClip, 6),
            (Kind::Pad, 7),
            (Kind::SequenceMatch, 8),
            (Kind::SequenceMismatch, 9),
        ];
        assert_eq!(cigar_reference_span(ops).unwrap(), 34);
    }
}
