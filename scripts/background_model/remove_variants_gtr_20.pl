#!/usr/bin/perl

## perl combine.pl <vcf>
## removes any variation >= 20% (real variation)

use warnings;
use strict;

my $infile =$ARGV[0];
my @line_split=();
my @percent_split=();
my $line="";

open(IP, "zcat $infile |") or die "could not zcat the infile";

while($line = <IP>){
	if (substr($line,0,1) eq "#") {
		print "$line";
		next;
	}
	chomp $line;
	@line_split = split(/\t/,$line);
	@percent_split = split(/:/,$line_split[9]);
	if ( $percent_split[4] <= "0.2" ) {
		print "$line\n";
	}
}

close IP;
exit;