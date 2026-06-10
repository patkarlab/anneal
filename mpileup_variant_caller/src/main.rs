use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::process;

// ---------------------------------------------------------------------------
// mpileup variant caller -- Rust reimplementation
//
// Faithfully replicates the logic of call_variants_from_mpileup.pl:
//   1. Parse each mpileup line (chr, pos, ref, depth, pileup, quals).
//   2. Walk the pileup string character-by-character, maintaining an
//      index_difference so that quality look-ups stay in sync.
//   3. Count SNVs and indels whose base quality exceeds the threshold;
//      decrement depth for bases that fail the filter.
//   4. Emit one VCF record per observed variant allele.
//
// Usage:
//   call_variants <mpileup> <min_base_qual_char> <output.vcf>
//
// The min_base_qual is an ASCII character (e.g. '5' = Phred 20 in Phred+33
// encoding), matching the Perl script's character-level comparison.
// ---------------------------------------------------------------------------

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        eprintln!(
            "Usage: {} <mpileup> <min_base_qual_char> <output.vcf>",
            args[0]
        );
        eprintln!();
        eprintln!("Arguments:");
        eprintln!("  mpileup             samtools mpileup output file");
        eprintln!("  min_base_qual_char  Minimum base quality as ASCII char (e.g. '5' = Phred 20)");
        eprintln!("  output.vcf          Output VCF file");
        process::exit(1);
    }

    let mpileup_path = &args[1];
    let min_qual_char = args[2].chars().next().unwrap_or('5');
    let output_path = &args[3];

    let reader: Box<dyn BufRead> = if mpileup_path == "-" {
        Box::new(BufReader::new(io::stdin()))
    } else {
        let f = File::open(mpileup_path).unwrap_or_else(|e| {
            eprintln!("ERROR: Cannot open {}: {}", mpileup_path, e);
            process::exit(1);
        });
        Box::new(BufReader::with_capacity(1024 * 1024, f))
    };

    let outfile = File::create(output_path).unwrap_or_else(|e| {
        eprintln!("ERROR: Cannot create {}: {}", output_path, e);
        process::exit(1);
    });
    let mut writer = BufWriter::with_capacity(512 * 1024, outfile);

    // VCF header
    writeln!(writer, "##fileformat=VCFv4.2").unwrap();
    writeln!(
        writer,
        "##source=mpileup_variant_caller_rust_v1.0.0"
    )
    .unwrap();
    writeln!(
        writer,
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth (quality-filtered)\">"
    )
    .unwrap();
    writeln!(
        writer,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )
    .unwrap();
    writeln!(
        writer,
        "##FORMAT=<ID=ALT,Number=1,Type=Integer,Description=\"Alternate allele observation count\">"
    )
    .unwrap();
    writeln!(
        writer,
        "##FORMAT=<ID=TOT,Number=1,Type=Integer,Description=\"Total depth (quality-filtered)\">"
    )
    .unwrap();
    writeln!(
        writer,
        "##FORMAT=<ID=FRAC,Number=1,Type=Float,Description=\"Variant allele fraction\">"
    )
    .unwrap();
    writeln!(
        writer,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    )
    .unwrap();

    let mut total_lines: u64 = 0;
    let mut total_variants: u64 = 0;

    for line_result in reader.lines() {
        let line = match line_result {
            Ok(l) => l,
            Err(_) => continue,
        };

        if line.is_empty() {
            continue;
        }

        total_lines += 1;

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            continue;
        }

        let chr = fields[0];
        let pos = fields[1];
        let ref_base = fields[2].to_uppercase();
        let mut depth: i64 = fields[3].parse().unwrap_or(0);
        let pileup = fields[4];
        let quals = fields[5];

        if depth == 0 || pileup.is_empty() {
            continue;
        }

        // Parse pileup string, counting variants
        let mut variant_counts: HashMap<String, i64> = HashMap::new();
        let mut index_diff: i64 = 0;
        let pileup_chars: Vec<char> = pileup.chars().collect();
        let qual_chars: Vec<char> = quals.chars().collect();
        let pileup_len = pileup_chars.len() as i64;

        let mut i: i64 = 0;
        while i < pileup_len {
            let c = pileup_chars[i as usize];

            match c {
                // Read start marker: skip ^ and the following mapping quality char
                '^' => {
                    i += 1; // skip the mapping quality character
                    index_diff += 2;
                }
                // Read end marker
                '$' => {
                    index_diff += 1;
                }
                // Reference match (forward/reverse)
                '.' | ',' => {
                    let qual_idx = (i - index_diff) as usize;
                    if qual_idx < qual_chars.len() {
                        if qual_chars[qual_idx] < min_qual_char {
                            depth -= 1;
                        }
                    }
                }
                // SNV (forward strand)
                'A' | 'C' | 'G' | 'T' => {
                    let qual_idx = (i - index_diff) as usize;
                    if qual_idx < qual_chars.len() {
                        if qual_chars[qual_idx] >= min_qual_char {
                            let key = c.to_string();
                            *variant_counts.entry(key).or_insert(0) += 1;
                        } else {
                            depth -= 1;
                        }
                    }
                }
                // SNV (reverse strand)
                'a' | 'c' | 'g' | 't' => {
                    let qual_idx = (i - index_diff) as usize;
                    if qual_idx < qual_chars.len() {
                        if qual_chars[qual_idx] >= min_qual_char {
                            let key = c.to_uppercase().to_string();
                            *variant_counts.entry(key).or_insert(0) += 1;
                        } else {
                            depth -= 1;
                        }
                    }
                }
                // Insertion
                '+' => {
                    let (indel_str, skip) = parse_indel(&pileup_chars, (i + 1) as usize);
                    let qual_idx = (i - index_diff) as usize;
                    if qual_idx < qual_chars.len() && qual_chars[qual_idx] >= min_qual_char {
                        let key = format!("+{}", indel_str.to_uppercase());
                        *variant_counts.entry(key).or_insert(0) += 1;
                    }
                    index_diff += 1 + skip as i64;
                    i += skip as i64;
                }
                // Deletion
                '-' => {
                    let (indel_str, skip) = parse_indel(&pileup_chars, (i + 1) as usize);
                    let qual_idx = (i - index_diff) as usize;
                    if qual_idx < qual_chars.len() && qual_chars[qual_idx] >= min_qual_char {
                        let key = format!("-{}", indel_str.to_uppercase());
                        *variant_counts.entry(key).or_insert(0) += 1;
                    }
                    index_diff += 1 + skip as i64;
                    i += skip as i64;
                }
                // Deletion placeholder in reference
                '*' | '#' => {
                    let qual_idx = (i - index_diff) as usize;
                    if qual_idx < qual_chars.len() {
                        if qual_chars[qual_idx] < min_qual_char {
                            depth -= 1;
                        }
                    }
                }
                // Skip N and other characters
                'N' | 'n' | '>' | '<' => {
                    depth -= 1;
                }
                _ => {}
            }

            i += 1;
        }

        // Emit variant records
        if depth <= 0 {
            continue;
        }

        for (allele, count) in &variant_counts {
            if *count == 0 {
                continue;
            }

            let af = *count as f64 / depth as f64;

            let (vcf_ref, vcf_alt) = if allele.starts_with('+') {
                // Insertion: REF=ref_base, ALT=ref_base+inserted_bases
                (ref_base.clone(), format!("{}{}", ref_base, &allele[1..]))
            } else if allele.starts_with('-') {
                // Deletion: REF=ref_base+deleted_bases, ALT=ref_base
                (format!("{}{}", ref_base, &allele[1..]), ref_base.clone())
            } else {
                // SNV
                (ref_base.clone(), allele.clone())
            };

            let gt = "0/1";

            writeln!(
                writer,
                "{}\t{}\t.\t{}\t{}\t.\t.\tDP={}\tGT:ALT:TOT:FRAC\t{}:{}:{}:{:.4}",
                chr, pos, vcf_ref, vcf_alt, depth, gt, count, depth, af
            )
            .unwrap();

            total_variants += 1;
        }
    }

    writer.flush().unwrap();
    eprintln!(
        "Processed {} positions, wrote {} variant records to {}",
        total_lines, total_variants, output_path
    );
}

/// Parse an indel from the pileup string starting after +/-.
/// Returns (indel_sequence, number_of_characters_consumed).
fn parse_indel(chars: &[char], start: usize) -> (String, usize) {
    // Read the numeric length
    let mut num_str = String::new();
    let mut pos = start;
    while pos < chars.len() && chars[pos].is_ascii_digit() {
        num_str.push(chars[pos]);
        pos += 1;
    }

    let length: usize = num_str.parse().unwrap_or(0);
    let mut seq = String::with_capacity(length);

    for j in 0..length {
        if pos + j < chars.len() {
            seq.push(chars[pos + j]);
        }
    }

    let total_consumed = num_str.len() + length;
    (seq, total_consumed)
}
