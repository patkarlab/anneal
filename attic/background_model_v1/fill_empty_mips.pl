#!/usr/bin/perl

## perl fill_empty_mips.pl <mip list>  <multiple variant file>
## input is a non-overlapping list of all the mips and a multiple variant file
## the program will fill in all positions with no variants 

use warnings;
use strict;

my $mipfile = $ARGV[0];
my $variantfile = $ARGV[1]; 
my $outfile = $variantfile . "\.filled";
my @var_line_split=();
my $current_chrom="";
my $current_coordinate="";
my $var_line="";
my $end_file=0;

open(MIP, "$mipfile");
open(VAR, "$variantfile");
open(OUT, '>', "$outfile");
$var_line = <VAR>;
print OUT $var_line;
$var_line = <VAR>;
chomp $var_line;
@var_line_split = (split /\t/,$var_line);
while(my $mip_line = <MIP>) {
   my @mip_line_split = (); 
   chomp $mip_line;
   @mip_line_split = (split /\t/,$mip_line);
   for(my $i=$mip_line_split[1]; $i <= $mip_line_split[2]; $i++) {
      if($var_line_split[1] == $i) {
         print OUT "$var_line\t$mip_line_split[0]\n";
         $var_line = <VAR>;
         if(!(eof(VAR))) {
            chomp $var_line;
            @var_line_split = (split /\t/,$var_line);
         }
         elsif($end_file == 0) {
            chomp $var_line;
            @var_line_split = (split /\t/,$var_line);
            $end_file = 1;
         }
      }
      else {
         print OUT "$mip_line_split[0]\t$i\tX\t$var_line_split[3]\t0\t0\t0\t0\t$mip_line_split[0]\n";
      }
   } 
   
}

close MIP;
close VAR;
close OUT;

exit;

