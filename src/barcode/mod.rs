//! Barcode extraction from FASTQ reads.
//!
//! Handles IDT xGen Duplex Seq adapter layout:
//!   Read structure: [UMI_bases][Spacer_bases][Genomic_insert]
//!
//! After extraction, R1_UMI + R2_UMI are concatenated and stored in the
//! FASTQ header so downstream alignment preserves the UMI tag.

pub mod extract;
