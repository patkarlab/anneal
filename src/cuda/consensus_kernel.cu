extern "C" __global__ void consensus_kernel(
    const unsigned int* __restrict__ packed_seqs,
    const unsigned char* __restrict__ qualities,
    const unsigned int* __restrict__ family_meta,
    const unsigned int num_families,
    const unsigned int cutoff_scaled,
    const unsigned char min_qual,
    unsigned char* __restrict__ consensus_out,
    unsigned char* __restrict__ quality_out,
    const unsigned int* __restrict__ output_offsets
) {
    unsigned int family_idx = blockIdx.x;
    if (family_idx >= num_families) return;

    unsigned int seq_offset = family_meta[family_idx * 3];
    unsigned int num_reads  = family_meta[family_idx * 3 + 1];
    unsigned int read_len   = family_meta[family_idx * 3 + 2];
    unsigned int out_offset = output_offsets[family_idx];

    extern __shared__ unsigned int base_counts[];

    for (unsigned int i = threadIdx.x; i < read_len * 5; i += blockDim.x) {
        base_counts[i] = 0;
    }
    __syncthreads();

    unsigned int words_per_read = (read_len + 15) / 16;
    unsigned int total_work = num_reads * read_len;

    for (unsigned int work_idx = threadIdx.x; work_idx < total_work; work_idx += blockDim.x) {
        unsigned int read_idx = work_idx / read_len;
        unsigned int pos = work_idx % read_len;

        unsigned int read_start = seq_offset + read_idx * words_per_read;
        unsigned int word_idx = pos / 16;
        unsigned int bit_offset = (pos % 16) * 2;
        unsigned int base_2bit = (packed_seqs[read_start + word_idx] >> bit_offset) & 0x3;

        unsigned int qual_idx = (seq_offset * 16) + read_idx * read_len + pos;
        unsigned char qual = qualities[qual_idx];

        if (qual >= min_qual) {
            atomicAdd(&base_counts[pos * 5 + base_2bit], 1);
        }
    }
    __syncthreads();

    for (unsigned int pos = threadIdx.x; pos < read_len; pos += blockDim.x) {
        unsigned int best_base = 4;
        unsigned int best_count = 0;
        unsigned int total_valid = 0;

        for (unsigned int b = 0; b < 4; b++) {
            unsigned int count = base_counts[pos * 5 + b];
            total_valid += count;
            if (count > best_count) {
                best_count = count;
                best_base = b;
            }
        }

        unsigned char cons_base;
        unsigned char cons_qual;

        if (total_valid > 0 && best_count * 1000 >= cutoff_scaled * total_valid) {
            const unsigned char bases[] = {'A', 'C', 'G', 'T'};
            cons_base = bases[best_base];
            unsigned int q = (best_count * 30 < 60) ? best_count * 30 : 60;
            cons_qual = (unsigned char)(q + 33);
        } else {
            cons_base = 'N';
            cons_qual = 33;
        }

        consensus_out[out_offset + pos] = cons_base;
        quality_out[out_offset + pos] = cons_qual;
    }
}

extern "C" __global__ void duplex_kernel(
    const unsigned char* __restrict__ sscs_pos_seq,
    const unsigned char* __restrict__ sscs_pos_qual,
    const unsigned char* __restrict__ sscs_neg_seq,
    const unsigned char* __restrict__ sscs_neg_qual,
    const unsigned int* __restrict__ pair_meta,
    const unsigned int num_pairs,
    unsigned char* __restrict__ dcs_seq_out,
    unsigned char* __restrict__ dcs_qual_out,
    const unsigned int* __restrict__ output_offsets
) {
    unsigned int pair_idx = blockIdx.x;
    if (pair_idx >= num_pairs) return;

    unsigned int offset = pair_meta[pair_idx * 2];
    unsigned int length = pair_meta[pair_idx * 2 + 1];
    unsigned int out_off = output_offsets[pair_idx];

    for (unsigned int pos = threadIdx.x; pos < length; pos += blockDim.x) {
        unsigned char base_pos = sscs_pos_seq[offset + pos];
        unsigned char base_neg = sscs_neg_seq[offset + pos];
        unsigned char qual_pos = sscs_pos_qual[offset + pos];
        unsigned char qual_neg = sscs_neg_qual[offset + pos];

        unsigned char out_base;
        unsigned char out_qual;

        if (base_pos == 'N' && base_neg == 'N') {
            out_base = 'N'; out_qual = 33;
        } else if (base_pos == 'N') {
            out_base = base_neg; out_qual = qual_neg;
        } else if (base_neg == 'N') {
            out_base = base_pos; out_qual = qual_pos;
        } else if (base_pos == base_neg) {
            out_base = base_pos;
            unsigned int combined = (unsigned int)qual_pos + (unsigned int)qual_neg - 33;
            out_qual = (unsigned char)((combined > 93) ? 93 : combined);
        } else {
            out_base = 'N'; out_qual = 33;
        }

        dcs_seq_out[out_off + pos] = out_base;
        dcs_qual_out[out_off + pos] = out_qual;
    }
}
