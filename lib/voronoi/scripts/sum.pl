#!/usr/bin/perl
# This script will take a .vol file created by voro++, and it will sum up the
# volumes of all the Voronoi cells. In many cases, this will equal the volume
# of the container, and it can be a useful check.
#
# It assumes that the volume is stored in the final column, so it works with
# both regular and polydisperse calculations.

open A,"@ARGV[0]" or die "Can't open input file";
while(<A>) {
	@A=split;
	$c+=pop @A;
}
print "$c\n";
