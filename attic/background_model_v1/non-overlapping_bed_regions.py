#!/usr/bin/env python3

import sys

def fix_overlaps(input_bed, output_bed):
	prev_end = {} 

	with open(input_bed) as bed_in, open(output_bed, "w") as bed_out:
		for line in bed_in:
			fields = line.strip().split("\t")

			chrom = fields[0]
			start = int(fields[1])
			end   = int(fields[2])

			if chrom in prev_end and start <= prev_end[chrom]:
				start = prev_end[chrom] + 1
			else:
				pass

			prev_end[chrom] = end

			if start < end:
				bed_out.write(f"{chrom}\t{start}\t{end}\n")
			else:
				sys.stderr.write(f"Invalid region skipped: {chrom}\t{start}\t{end}\n")
					

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print(f"Usage: {sys.argv[0]} <input.bed> <output.bed>")
		sys.exit(1)

	fix_overlaps(sys.argv[1], sys.argv[2])
