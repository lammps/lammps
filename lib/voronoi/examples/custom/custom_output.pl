#!/usr/bin/perl

# Open the custom output from the custom_output.cc program
open A,"packing.custom1" or die "Can't open file \"packing.custom1\"\n";

# Open the POV-Ray file
open B,">custom_output_p.pov" or die "Can't open output file\n";

# Loop over all lines in the packing.custom1 file
while(<A>) {

	# Use a regular expression to get the particle position and the number
	# of faces of the Voronoi cell. These will be stored in the variables
	# $1 and $2.
	m/pos=\((.*)\).*faces=(\d*)/;

	# Print a sphere to the POV-Ray file, giving it a different texture
	# depending on the number of faces of the Voronoi cell
	print B "sphere{<$1>,0.5 texture{t$2}}\n";
}

# Close the two files
close A;
close B;
