#![cfg(feature = "gpu")]
#![allow(dead_code)]
//! GPU context initialization and management.

use anyhow::{Context, Result};
use log::info;
use std::sync::Arc;
use super::kernels::ConsensusKernels;

pub struct GpuContext {
    device: Arc<cudarc::driver::CudaDevice>,
    device_name: String,
    pub kernels: ConsensusKernels,
}

impl GpuContext {
    pub fn new(device_id: usize) -> Result<Self> {
        let device = cudarc::driver::CudaDevice::new(device_id)
            .context("Failed to initialize CUDA device")?;

        let device_name = format!("CUDA Device {}", device_id);
        info!("GPU initialized: {}", device_name);

        let kernels = ConsensusKernels::new(device.clone())
            .context("Failed to load CUDA kernels")?;

        Ok(Self {
            device,
            device_name,
            kernels,
        })
    }

    pub fn device_name(&self) -> &str {
        &self.device_name
    }

    pub fn device(&self) -> &Arc<cudarc::driver::CudaDevice> {
        &self.device
    }
}

/// 2-bit DNA encoding for compact GPU representation.
pub mod encoding {
    #[inline]
    pub fn encode_base_2bit(base: u8) -> u8 {
        match base {
            b'A' | b'a' => 0b00,
            b'C' | b'c' => 0b01,
            b'G' | b'g' => 0b10,
            b'T' | b't' => 0b11,
            _ => 0b00,
        }
    }

    #[inline]
    pub fn decode_base_2bit(bits: u8) -> u8 {
        match bits & 0b11 {
            0b00 => b'A',
            0b01 => b'C',
            0b10 => b'G',
            0b11 => b'T',
            _ => unreachable!(),
        }
    }

    pub fn pack_sequence(seq: &[u8]) -> Vec<u32> {
        let num_words = (seq.len() + 15) / 16;
        let mut packed = vec![0u32; num_words];
        for (i, &base) in seq.iter().enumerate() {
            let word_idx = i / 16;
            let bit_offset = (i % 16) * 2;
            packed[word_idx] |= (encode_base_2bit(base) as u32) << bit_offset;
        }
        packed
    }

    pub fn unpack_sequence(packed: &[u32], len: usize) -> Vec<u8> {
        let mut seq = Vec::with_capacity(len);
        for i in 0..len {
            let word_idx = i / 16;
            let bit_offset = (i % 16) * 2;
            let bits = ((packed[word_idx] >> bit_offset) & 0b11) as u8;
            seq.push(decode_base_2bit(bits));
        }
        seq
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn test_roundtrip() {
            let seq = b"ACGTACGTNNACGTTT";
            let packed = pack_sequence(seq);
            let unpacked = unpack_sequence(&packed, seq.len());
            assert_eq!(unpacked[0], b'A');
            assert_eq!(unpacked[1], b'C');
            assert_eq!(unpacked[2], b'G');
            assert_eq!(unpacked[3], b'T');
        }

        #[test]
        fn test_pack_density() {
            let seq = vec![b'A'; 100];
            let packed = pack_sequence(&seq);
            assert_eq!(packed.len(), 7);
        }
    }
}
