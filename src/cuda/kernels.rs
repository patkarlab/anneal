#![cfg(feature = "gpu")]
//! CUDA kernel loading and launch for consensus calling.

use anyhow::{Context, Result};
use cudarc::driver::safe::LaunchConfig;
use cudarc::driver::LaunchAsync;
use cudarc::driver::CudaDevice;
use cudarc::nvrtc::Ptx;
use std::sync::Arc;

/// Embedded PTX compiled from consensus_kernel.cu
const CONSENSUS_PTX: &str = include_str!("consensus_kernel.ptx");

pub const BLOCK_SIZE: u32 = 256;

/// Loaded GPU kernels ready for launch.
pub struct ConsensusKernels {
    device: Arc<CudaDevice>,
}

impl ConsensusKernels {
    pub fn new(device: Arc<CudaDevice>) -> Result<Self> {
        let ptx = Ptx::from_src(CONSENSUS_PTX);
        device
            .load_ptx(ptx, "consensus", &["consensus_kernel", "duplex_kernel"])
            .context("Failed to load consensus PTX")?;

        Ok(Self { device })
    }

    /// Run SSCS consensus on GPU for a batch of families.
    pub fn run_consensus_batch(
        &self,
        packed_seqs: &[u32],
        qualities: &[u8],
        family_meta: &[(u32, u32, u32)],  // (seq_offset_words, num_reads, read_len)
        cutoff: f64,
        min_qual: u8,
    ) -> Result<(Vec<u8>, Vec<u8>)> {
        let num_families = family_meta.len();
        if num_families == 0 {
            return Ok((Vec::new(), Vec::new()));
        }

        // Flatten family metadata
        let meta_flat: Vec<u32> = family_meta
            .iter()
            .flat_map(|&(off, nr, rl)| vec![off, nr, rl])
            .collect();

        // Compute output offsets
        let mut output_offsets = Vec::with_capacity(num_families);
        let mut total_out: u32 = 0;
        for &(_, _, rl) in family_meta {
            output_offsets.push(total_out);
            total_out += rl;
        }

        // Copy to GPU
        let d_packed = self.device.htod_sync_copy(packed_seqs)
            .context("htod packed_seqs")?;
        let d_quals = self.device.htod_sync_copy(qualities)
            .context("htod qualities")?;
        let d_meta = self.device.htod_sync_copy(&meta_flat)
            .context("htod family_meta")?;
        let d_offsets = self.device.htod_sync_copy(&output_offsets)
            .context("htod output_offsets")?;
        let mut d_cons_out = self.device.alloc_zeros::<u8>(total_out as usize)
            .context("alloc consensus_out")?;
        let mut d_qual_out = self.device.alloc_zeros::<u8>(total_out as usize)
            .context("alloc quality_out")?;

        let cutoff_scaled = (cutoff * 1000.0) as u32;
        let num_fam_u32 = num_families as u32;

        // Shared memory for base counts
        let max_rl = family_meta.iter().map(|m| m.2).max().unwrap_or(0);
        let shared_mem = max_rl * 5 * 4;

        let cfg = LaunchConfig {
            grid_dim: (num_fam_u32, 1, 1),
            block_dim: (BLOCK_SIZE, 1, 1),
            shared_mem_bytes: shared_mem,
        };

        let func = self.device
            .get_func("consensus", "consensus_kernel")
            .context("get consensus_kernel")?;

        unsafe {
            func.launch(
                cfg,
                (
                    &d_packed,
                    &d_quals,
                    &d_meta,
                    num_fam_u32,
                    cutoff_scaled,
                    min_qual,
                    &mut d_cons_out,
                    &mut d_qual_out,
                    &d_offsets,
                ),
            )
        }
        .context("launch consensus_kernel")?;

        self.device.synchronize().context("sync after consensus")?;

        let cons_out = self.device.dtoh_sync_copy(&d_cons_out)
            .context("dtoh consensus")?;
        let qual_out = self.device.dtoh_sync_copy(&d_qual_out)
            .context("dtoh quality")?;

        Ok((cons_out, qual_out))
    }

    /// Run duplex consensus on GPU for paired SSCS.
    pub fn run_duplex_batch(
        &self,
        pos_seqs: &[u8],
        pos_quals: &[u8],
        neg_seqs: &[u8],
        neg_quals: &[u8],
        pair_meta: &[(u32, u32)],  // (offset, length)
    ) -> Result<(Vec<u8>, Vec<u8>)> {
        let num_pairs = pair_meta.len();
        if num_pairs == 0 {
            return Ok((Vec::new(), Vec::new()));
        }

        let meta_flat: Vec<u32> = pair_meta
            .iter()
            .flat_map(|&(off, len)| vec![off, len])
            .collect();

        let mut output_offsets = Vec::with_capacity(num_pairs);
        let mut total_out: u32 = 0;
        for &(_, len) in pair_meta {
            output_offsets.push(total_out);
            total_out += len;
        }

        let d_pos_seq = self.device.htod_sync_copy(pos_seqs)?;
        let d_pos_qual = self.device.htod_sync_copy(pos_quals)?;
        let d_neg_seq = self.device.htod_sync_copy(neg_seqs)?;
        let d_neg_qual = self.device.htod_sync_copy(neg_quals)?;
        let d_meta = self.device.htod_sync_copy(&meta_flat)?;
        let d_offsets = self.device.htod_sync_copy(&output_offsets)?;
        let mut d_dcs_seq = self.device.alloc_zeros::<u8>(total_out as usize)?;
        let mut d_dcs_qual = self.device.alloc_zeros::<u8>(total_out as usize)?;

        let num_pairs_u32 = num_pairs as u32;

        let cfg = LaunchConfig {
            grid_dim: (num_pairs_u32, 1, 1),
            block_dim: (BLOCK_SIZE, 1, 1),
            shared_mem_bytes: 0,
        };

        let func = self.device
            .get_func("consensus", "duplex_kernel")
            .context("get duplex_kernel")?;

        unsafe {
            func.launch(
                cfg,
                (
                    &d_pos_seq,
                    &d_pos_qual,
                    &d_neg_seq,
                    &d_neg_qual,
                    &d_meta,
                    num_pairs_u32,
                    &mut d_dcs_seq,
                    &mut d_dcs_qual,
                    &d_offsets,
                ),
            )
        }
        .context("launch duplex_kernel")?;

        self.device.synchronize()?;

        let dcs_seq = self.device.dtoh_sync_copy(&d_dcs_seq)?;
        let dcs_qual = self.device.dtoh_sync_copy(&d_dcs_qual)?;

        Ok((dcs_seq, dcs_qual))
    }
}
