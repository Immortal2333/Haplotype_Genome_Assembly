#!/usr/bin/perl
# -*- encoding: utf-8 -*-
'''
File    :   getFaLen.pl
Time    :   2022/09/01
Author  :   Qing Zhang, Xu Wang
Version :   1.0
Contact :   wangxu2018@cau.edu.cn
License :   (C)Copyright 2021-2023, CAAS ShenZhen
Desc    :   Get the length of chromosomes in the genome.fasta
'''

###This script is used to get length information for each sequence in a fasta file

use Getopt::Std;
getopts "i:o:";

$input = $opt_i;
$output = $opt_o;

if ((!defined $opt_i) || (!defined $opt_o) ) {
        die "************************************************************************
        Usage: perl getFalen.pl -i input.fasta -o output
          -h : help and usage.
          -i : input fasta
          -o : output
************************************************************************\n";
}
open(OUT, ">$output") or die"";
open(IN, $input) or die"";
$/='>';
<IN>;
while(<IN>){
	chomp;
	($gene,$seq) = split(/\n/,$_,2);
        $gene =~ s/\s+.*//g;
	$seq  =~ s/\s+//g;
	$len  = length $seq;
	print OUT "$gene	$len\n";
	}
close IN;
close OUT;

