//! Consensus calling configuration.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConsensusConfig {
    /// Minimum fraction of identical bases at a position to call consensus (0.0-1.0).
    /// Default: 0.7 (per Kennedy et al.)
    pub cutoff: f64,

    /// Minimum Phred quality score to include a base in consensus.
    /// Default: 30 (Q30)
    pub min_qual: u8,

    /// Minimum number of reads to form an SSCS family.
    /// Default: 1 (singletons are handled separately)
    pub min_family_size: usize,

    /// Maximum family size (randomly downsample beyond this).
    /// Default: 1000
    pub max_family_size: usize,

    /// Whether to use GPU acceleration.
    pub use_gpu: bool,

    /// CUDA device index.
    pub gpu_device: usize,

    /// Enable singleton correction workflow.
    pub singleton_correction: bool,
}

impl Default for ConsensusConfig {
    fn default() -> Self {
        Self {
            cutoff: 0.7,
            min_qual: 30,
            min_family_size: 1,
            max_family_size: 1000,
            use_gpu: false,
            gpu_device: 0,
            singleton_correction: true,
        }
    }
}
