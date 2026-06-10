//! Consensus sequence generation.
//!
//! Implements both SSCS (Single-Strand Consensus Sequence) and
//! DCS (Duplex Consensus Sequence) generation with CPU and GPU paths.

pub mod config;
pub mod cpu;
pub mod dcs;
pub mod pipeline;
pub mod sscs;
