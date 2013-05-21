#!/usr/bin/perl
# This script can be used to generate random lists of particles in the unit
# box 0<x<1, 0<y<1, 0<z<1, that can be used as test input to voro++. The -r
# switch adds extra radial information for use with the Voronoi radical
# tessellation.
use Getopt::Std;
getopts("hr:");

if ($opt_h) {
	print "Usage: random.pl [switches] <number of particles> <filename>\n";
	print "Switches:\n";
	print "-h                Print this information\n";
	print "-r <max radius>   Add extra radial column\n";
	exit 0;
};

@ARGV==2 or die "Exactly two command line arguments required\n";
open A,">@ARGV[1]" or die "Can't open output file";

foreach (1..@ARGV[0]) {
	printf A "$_ %f %f %f", rand(), rand(), rand();
	printf A " %f",$opt_r*rand() if $opt_r;
	printf A "\n";
}
