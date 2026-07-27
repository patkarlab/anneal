#!/usr/bin/env python3
# This script will take input files generated from fill_empty_mips.pl
# The output of the script will be a matrix file for calculate_beta_P_values_vcf.pl

import argparse
import sys
import linecache
import statistics

def count_beta_dist(arguments):
	# print ("samples are", arguments.samples, type(arguments.samples))
	# Check if all samples have the same no. of lines else exit
	line_count = 0
	No_of_samples = len(arguments.samples)
	for index, sample_file in enumerate (arguments.samples):
		print ("Sample no", index + 1, ":", sample_file )
		with open (sample_file) as infile:
			num_lines = sum(1 for _ in infile)	# Counting the no. of lines in the input file
			if ( index == 0):
				line_count = num_lines
			else:
				if (num_lines != line_count):
					print ("No. of lines in", sample_file, "not same as previous file, Stopping the execution")
					sys.exit () 	# Stop the analysis if inputs are not of same length
				else:
					pass
	
	# Opening the output file
	print ("Writing output to", arguments.output)
	outfile = open(arguments.output, "w")
	print ("chr", "pos", "alpha A", "beta A", "alpha T", "beta T", "alpha G", "beta G", "alpha C", "beta C", sep='\t', file=outfile)

	# The no. of columns are 9 , their names are : chr, pos, ref, tagcount, A, T, G, C , "null"
	mean_percent_error_A = []
	mean_percent_error_T = []
	mean_percent_error_G = []
	mean_percent_error_C = []
	total_tag_count = []
	total_A_count = []
	total_T_count = [] 
	total_G_count = []
	total_C_count = []
	total_all_mutations_list = []
	stdev_A = []
	stdev_T = []
	stdev_G = []
	stdev_C = []
	A_beta_dist_alpha = []
	A_beta_dist_beta = []
	T_beta_dist_alpha = []
	T_beta_dist_beta = []
	G_beta_dist_alpha = []
	G_beta_dist_beta = []
	C_beta_dist_alpha =[]
	C_beta_dist_beta =[]
	default_error_rate = (1 / 15000)	# Taken from supplementary of Haematologica 2017, Volume 102(9):1549-1557
	mean_default_error = default_error_rate
	stddevp_default_error = mean_default_error * 2
	# Iterate over line numbers
	for lines in range(2, line_count+1):	# Ignore the first line with column names
		# The following list are for samples
		tag_count_list = []
		a_count_list = []
		t_count_list = []
		g_count_list = []
		c_count_list = []
		a_error_list = []
		t_error_list = []
		g_error_list = []
		c_error_list = []
		for ind, samp_file in enumerate (arguments.samples):
			line_info = linecache.getline(samp_file, lines).split('\t')
			chr_id = line_info[0]
			pos = line_info[1]
			tag_count = int (line_info[3])
			a_count = int (line_info[4])
			t_count = int (line_info[5])
			g_count = int (line_info[6])
			c_count = int (line_info[7])
			a_error_count = a_count / tag_count
			t_error_count = t_count / tag_count
			g_error_count = g_count / tag_count
			c_error_count = c_count / tag_count

			tag_count_list.append(tag_count)
			a_count_list.append(a_count)
			t_count_list.append(t_count)
			g_count_list.append(g_count)
			c_count_list.append(c_count)
			a_error_list.append(a_error_count)
			t_error_list.append(t_error_count)
			g_error_list.append(g_error_count)
			c_error_list.append(c_error_count)

		# total tagcount and counts
		# total_tag_count.append(sum(tag_count_list))
		# total_A_count.append(sum(a_count_list))
		# total_T_count.append(sum(t_count_list))
		# total_G_count.append(sum(g_count_list))
		# total_C_count.append(sum(c_count_list))

		# Mean % error . if the sum of errors is < 1 in 15000, use default else use the original value
		sum_a_error_list = sum(a_error_list)
		print (sum_a_error_list)
		sum_t_error_list = sum(t_error_list)
		sum_g_error_list = sum(g_error_list)
		sum_c_error_list = sum(c_error_list)

		if ( sum_a_error_list / No_of_samples < default_error_rate):
			mean_a_error = mean_default_error
			stddevp_a_error = stddevp_default_error
		else:
			mean_a_error = (sum_a_error_list / No_of_samples)
			stddevp_a_error = statistics.pstdev(a_error_list)

		if ( sum_t_error_list / No_of_samples < default_error_rate ):
			mean_t_error = mean_default_error
			stddevp_t_error = stddevp_default_error
		else:
			mean_t_error = (sum_t_error_list / No_of_samples)
			stddevp_t_error = statistics.pstdev(t_error_list)
		
		if ( sum_g_error_list / No_of_samples < default_error_rate ):
			mean_g_error = mean_default_error
			stddevp_g_error = stddevp_default_error
		else:
			mean_g_error = (sum_g_error_list / No_of_samples)
			stddevp_g_error = statistics.pstdev(g_error_list)
		
		if ( sum_c_error_list / No_of_samples < default_error_rate):
			mean_c_error = mean_default_error
			stddevp_c_error = stddevp_default_error
		else:
			mean_c_error = ( sum_c_error_list / No_of_samples)
			stddevp_c_error = statistics.pstdev(c_error_list)

		a_beta_dist_alpha = ((1 - mean_a_error) / (stddevp_a_error ** 2)) * (mean_a_error ** 2)
		a_beta_dist_beta = a_beta_dist_alpha * ((1 / mean_a_error ) - 1)
		t_beta_dist_alpha = ((1 - mean_t_error) / (stddevp_t_error ** 2)) * (mean_t_error ** 2)
		t_beta_dist_beta = t_beta_dist_alpha * ((1 / mean_t_error) - 1)
		g_beta_dist_alpha = ((1 - mean_g_error) / (stddevp_g_error ** 2)) * (mean_g_error ** 2)
		g_beta_dist_beta = g_beta_dist_alpha * ((1 / mean_g_error) - 1)
		c_beta_dist_alpha = ((1 - mean_c_error) / (stddevp_c_error ** 2)) * (mean_c_error ** 2)
		c_beta_dist_beta = c_beta_dist_alpha * ((1 / mean_c_error) - 1)

		# mean_percent_error_A.append(mean_a_error)
		# mean_percent_error_T.append(mean_t_error)
		# mean_percent_error_G.append(mean_g_error)
		# mean_percent_error_C.append(mean_c_error)

		# Total all mutations
		# total_all_mutations_list.append( (sum(a_count_list) + sum(t_count_list) + sum(g_count_list) + sum(c_count_list)) / sum(tag_count_list))

		# Standard deviations of errors
		# stdev_A.append(stddevp_a_error)
		# stdev_T.append(stddevp_t_error)
		# stdev_G.append(stddevp_g_error)
		# stdev_C.append(stddevp_c_error)

		# Beta distribution parameters
		# A_beta_dist_alpha.append(a_beta_dist_alpha)
		# A_beta_dist_beta.append(a_beta_dist_beta)
		# T_beta_dist_alpha.append(t_beta_dist_alpha)
		# T_beta_dist_beta.append(t_beta_dist_beta)
		# G_beta_dist_alpha.append(g_beta_dist_alpha)
		# G_beta_dist_beta.append(g_beta_dist_beta)
		# C_beta_dist_alpha.append(c_beta_dist_alpha)
		# C_beta_dist_beta.append(c_beta_dist_beta)

		# chr_id and pos are obtained from the last input file
		print (chr_id, pos, a_beta_dist_alpha, a_beta_dist_beta, t_beta_dist_alpha, t_beta_dist_beta, g_beta_dist_alpha, g_beta_dist_beta, c_beta_dist_alpha, c_beta_dist_beta, sep='\t', file=outfile)

	outfile.close()

def parse_arguments():
	parser = argparse.ArgumentParser(description="Calculation of beta distribution matrix")
	parser.add_argument('--output', help='name of the output matrix file eg: beta_matrix.txt')	# positional argument
	parser.add_argument('--samples', nargs='*', help='list of files to make the beta matrix')
	return parser.parse_args()

def main ():
	args = parse_arguments()
	count_beta_dist(args)

if __name__ == "__main__":
	main()