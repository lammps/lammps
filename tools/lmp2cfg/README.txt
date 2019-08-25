README.txt

lmp2cfg

Version 1.0
Ara Kooser
8/1/04

Contents
I.   Before you start
II.  Introduction
III. Compiling the program
IV.  Running the program


I. Before you start

1) Read the READMEFIRST file
2) You will need either f77 or g77 compiler
3) A text editor (VI or emacs)
4) AtomEye( http://164.107.79.177/Archive/Graphics/A/) good as of 8/1/04
5) A sense of adventure and humor


II. Introduction

The lmp2cfg program converts a LAMMPS atom dump file into a series of .cfg
files.(000001.cfg, 000002.cfg,...). These files can be read into AtomEye
and made into a movie or viewed indivdually.


III. Compiling the program

The program should comply with either command f77 lmp2cfg.f or
g77 lmp2cfg.f

IV. Running the program

You need to create a user input file. In the examples folder you will find
an example input script, LAMMPS dump file, and the correct .cfg output. 
The input script reads like this:

2898           #total number of atoms in system (may be more than in dump)
7              #number of atom types in your LAMMPS file       
'dump.atom'    #name of the LAMMPS dump file, you need the ' '
1              #first frame
10             #last frame
1              #first atom type
26.98154       #atomic weight
'ao'           #atom name
2              #second atom type
15.9994        #atomic weight
'oh'           #atom name
 

Make sure the input file and the atom dump file are in the same folder.
On the command line type

lmp2cfg < input_file.txt > out

You should get several .cfg files. For reading into AtomEye and making
movies see the AtomEye homepage.


If you get an error like:

open: 'new' file exists
apparent state: unit 21 named 00001.cfg
lately writing sequential formatted external IO
Abort

you need to first remove old copies of the cfg files.

If you get an error like:

open: No such file or directory
apparent state: unit 9 named dump.atom
last format: list io
lately reading direct formatted external IO
Abort

you need to check that the name of the dump file matches the name
in the input file, and that you enclosed it in single quotes.
