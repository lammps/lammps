README.txt

lmp2traj

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
4) LAMMPS
5) A graphing program that can do contours
6) A sense of adventure and humor

II. Introduction

This program will take a LAMMPS atom dump file and provide the following three
files.
1) data for making contour maps
2) density profile
3) dipole information

III. Compiling the program

To compile the program run either ./f77 lmp2traj.f or
./g77 traj.f

IV. Running the program

First you need to create a user input file. There is an example of an input file
in the examples folder.

The input file reads like this:

'atom'                        # dump file name, needs the ' '
1                             # first frame
60                            # last frame
38.26119                      # x dimension of the box
44.26119                      # y dimension of the box
48.33150                      # z dimension of the box
90.                           # angles of the box, always 90
90.                           # angles of the box, always 90                      
90.                           # angles of the box, always 90
82844.6                       # volumne of the box in cubic Angstroms
5                             # water oxygen atom type from LAMMPS (#)
6                             # water hydrogen atom type from LAMMPS(#)
0.                            # leave at 0 
5.                            # number of atom types
0.                            # z shift leave at 0
5.                            # number of density maps
'Surface (1) ho'              # Enter name/description of atom
2                             # atom type number from LAMMPS
'ho'                          # Column name for data
0                             # Defines inner sphere, in A
48.3                          # Defines out sphere, in A


Make sure you have the input file and the LAMMPS atom dump file in the same
directory.

To run the program type
   lmp2traj < inp_file.txt > out

This should give you three files like in the examples/output folder. 





