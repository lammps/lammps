This directory contains scripts used for converting AMBER parameter files
into moltemplate (.LT) format.  When a newer version of the AMBER parameters
is eventually published, you can use these scripts to convert the new files
again.  (Some tinkering may be necessary.)

The main bash script is a wrapper which simply splits up the parameter (".dat")
file into fragments which (it thinks) correspond to the mass, pair, bond,
angle, dihedral, and improper section of the original .dat file.
(However sometimes it gets this wrong and you have to split it up manually!)

Then this bash script invokes the relevant python script to convert 
each section into .LT format:
amberparm_to_mass.py
amberparm_to_pair.py
amberparm_to_bond.py
amberparm_to_angle.py
amberparm_to_dihedral.py
amberparm_to_improper.py
In case this goes wrong, you may have to run these scripts manaully.


Find out how to run this bash script by invoking it without any arguments:

./amberparm2lt.sh

------------ IMPORTANT ------------

BEFORE YOU RUN THIS SCRIPT, BE SURE TO CHANGE THE ORDER OF THE IMPROPER DIHEDRAL
PARAMETERS SO THAT THE "SPECIFIC" IMPROPER DIHEDRALS APPEAR LAST, AND THE
"GENERIC" IMPROPER DIHEDRALS APPEAR FIRST.

For example replace these two lines: 

X -o -c -o          1.1          180.          2.           JCC,7,(1986),230
X -X -c -o          10.5         180.          2.           JCC,7,(1986),230

with these two lines:

X -X -c -o          10.5         180.          2.           JCC,7,(1986),230
X -o -c -o          1.1          180.          2.           JCC,7,(1986),230

Why:
This is the order that moltemplate expects: generic first. specific last.
So far only the improper dihedral parameters in the gaff.dat file seem
to violate this order.  The bonds, angles and dihedrals seem to obey this,
but check to make sure.


There is a discussion of these parameters here:
http://structbio.vanderbilt.edu/archives/amber-archive/2005/3444.php

excerpt:

> > In the parm99 file (for example), sometimes the wild-card is used, as it
> > is done in the following example:
> >
> > X -X -C -O 10.5 180. 2. JCC,7,(1986),230
> >
> > The first example is the specific case while the second one is the generic
> > case. In page # 257 of the AMBER Manual, it is talking about Dihedral
> > Angle, and how these dihedral parameters are used to calculate the
> > energies. I am wondering what the difference between generic and specific
> > case is for improper torsions.
>
> "specific" torsions are search for first, and used if a match is found. If
> no match is found, then a search is made to see if a "generic" (aka wild-card)
> torsion with match. 
> ...good luck...dac 

Good luck

-Andrew
2014-4-19

