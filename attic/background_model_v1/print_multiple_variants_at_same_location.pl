#!/usr/bin/perl

## perl print_multiple_variants_at_same_location.pl <vcf file>
## input is a sorted vcf file the program prints all unqiue locations with more than one variant

use warnings;
use strict;

my $infile =$ARGV[0];
my @line_split=();
my $previous_chrom="";
my $previous_coordinate="";
my $previous_ref="";
my $previous_count=0;
my $previous_line="";
my $line="";
my $A_variant=0;
my $T_variant=0;
my $G_variant=0;
my $C_variant=0;
my $total_tag_count=0;


open(IP, "$infile");
$line = <IP>;
while(substr($line,0,1) eq "#") {
	$line = <IP>;
}
chomp $line;
@line_split = split(/\t/,$line);
$previous_chrom=$line_split[0];
$previous_coordinate=$line_split[1];
$previous_count=$line_split[2];
$previous_ref=$line_split[3];
$previous_line=$line;
$total_tag_count= (split /,/, ((split /:/, $line_split[9])[2]))[0] + (split /,/, ((split /:/, $line_split[9])[2]))[1];
if($line_split[4] eq "A") {
	$A_variant=(split /,/, ((split /:/, $line_split[9])[2]))[1];
}
elsif($line_split[4] eq "T") {
	$T_variant=(split /,/, ((split /:/, $line_split[9])[2]))[1];
}
elsif($line_split[4] eq "G") {
	$G_variant=(split /,/, ((split /:/, $line_split[9])[2]))[1];
}
else {
	$C_variant=(split /,/, ((split /:/, $line_split[9])[2]))[1];
}
print "chr\tpos\tref\ttagcount\tA\tT\tG\tC\n";
while($line = <IP>){
	chomp $line;
	@line_split = split(/\t/,$line);
	if($line_split[1] eq $previous_coordinate){
		if($line_split[4] eq "A") {
			$A_variant=(split /,/, ((split /:/, $line_split[9])[2]))[1];
		}
		elsif($line_split[4] eq "T") {
			$T_variant=(split /,/, ((split /:/, $line_split[9])[2]))[1];
		}
		elsif($line_split[4] eq "G") {
			$G_variant=(split /,/, ((split /:/, $line_split[9])[2]))[1];
		}
		else {
			$C_variant=(split /,/, ((split /:/, $line_split[9])[2]))[1];
		}
		$previous_line=$line;
	}
	else {
		# print "$previous_line\n";
		print "$previous_chrom\t$previous_coordinate\t";
		print "$previous_ref\t$total_tag_count\t";
		print "$A_variant\t$T_variant\t$G_variant\t$C_variant\n";
		$A_variant=0;
		$T_variant=0;
		$G_variant=0;
		$C_variant=0;
		if($line_split[4] eq "A") {
			$A_variant=(split /,/, ((split /:/, $line_split[9])[2]))[1];
		}
		elsif($line_split[4] eq "T") {
			$T_variant=(split /,/, ((split /:/, $line_split[9])[2]))[1];
		}
		elsif($line_split[4] eq "G") {
			$G_variant=(split /,/, ((split /:/, $line_split[9])[2]))[1];
		}
		else {
			$C_variant=(split /,/, ((split /:/, $line_split[9])[2]))[1];
		}
		$previous_chrom=$line_split[0];
		$previous_coordinate=$line_split[1];
		$previous_count=$line_split[2];
		$previous_ref=$line_split[3];
		$previous_line=$line;
		$total_tag_count=(split /,/, ((split /:/, $line_split[9])[2]))[0] + (split /,/, ((split /:/, $line_split[9])[2]))[1];
	}
}

print "$previous_chrom\t$previous_coordinate\t";
print "$previous_ref\t$total_tag_count\t";
print "$A_variant\t$T_variant\t$G_variant\t$C_variant\n";

close IP;
exit;