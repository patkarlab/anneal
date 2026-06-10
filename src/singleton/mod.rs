//! Singleton correction: rescue unpaired reads.
//!
//! Two strategies:
//!   1. Singleton + complementary SSCS: use opposite strand SSCS to correct
//!   2. Singleton + complementary singleton: combine both strands

pub mod correction;
