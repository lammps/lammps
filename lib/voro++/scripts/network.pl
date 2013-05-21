#!/usr/bin/perl
# This script takes a .vol file that has neighbor information created with
# the custom output string "%i %q %v %n", and creates a map of the network, by
# drawing a line between all particles which are neighbors

open A,"@ARGV[0].vol" or die "Can't open input file";
open B,">@ARGV[0].net" or die "Can't open output file";

$maxn=0;

while(<A>) {
	@A=split;
	$n=@A[0];
	$x[$n]=@A[1];
	$y[$n]=@A[2];
	$z[$n]=@A[3];
	$c[$n]=$#A-5;
	$l[$n][$_-5]=@A[$_] foreach 5..$#A;
	$maxn=$n if $n>$maxn;
}

foreach $n (1..$maxn) {
	foreach (0..$c[$n]) {
		$j=$l[$n][$_];
		print B "$x[$j] $y[$j] $z[$j]\n$x[$n] $y[$n] $z[$n]\n\n\n" if $j>=$n;
	}
}
