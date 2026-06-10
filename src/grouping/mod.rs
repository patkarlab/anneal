//! Read grouping by UMI barcode and genomic coordinates.
//!
//! Reads sharing the same [Barcode]_[Chr]_[Pos]_[MateChr]_[MatePos]_[Strand]_[ReadNum]
//! tag originated from the same strand of the same original DNA molecule.
//! These families are collapsed into Single-Strand Consensus Sequences (SSCS).

pub mod families;
