#!/usr/sbin/perl 

die "Usage: crd2pdb crdfile pdbfile\n" if $#ARGV!=1;

$occ=1.00;

open(CRD, $ARGV[0]) || die "couldn't open $ARGV[0]:$!\n";
open(PDB, ">$ARGV[1]") || die "couldn't open $ARGV[1]:$!\n";

while (<CRD>)
	{
	if (/^\*/)   # comment/header
		{
		s/\*\s*(.*)/REMARK $1/;
		chomp;
		print PDB "$_\n";
		}
	elsif (/^\s*(\d+)\s*$/) # number of atoms
		{
		$num_atoms=$1;
		}
	else    #coordinate line
		{
		($atom_num, $res_num, $res_name, $atom_name, $x,$y,$z, $segment,
		$res_in_segment, $bfactor)=split(' ');
#		$res_name =~ s/^(...).*/$1/;
		printf PDB
#			"ATOM  %5d %4s %4s  %4d   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n",
			"ATOM  %5d %4s %4s %4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n",
			$atom_num, $atom_name, $res_name, $res_in_segment, $x,$y,$z, 
			$occ, $bfactor, $segment;
		}
		
	}
print PDB "END\n";
