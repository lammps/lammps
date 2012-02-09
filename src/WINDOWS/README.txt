INSTRUCTIONS FOR COMPILING LAMMPS WITH VISUAL STUDIO 2005

Note: These instructions and this framework for building LAMMPS under
Windows via Visual Studio, were created by Ilya Valuev (JIHT), valuev
at physik.hu-berlin.de.  Please contact him for questions or updates.

There are 3 VS projects provided:
LAMMPS.vcproj            -- minimal package set: KSPACE, MANYBODY, MOLECULE
LAMMPS-std.vcproj        -- standard package set (except for REPLICA, REAX, GPU)
LAMMPS-std-user.vcproj   -- standard set with some user packages
 
Each of the projects has configurations to compile either with MPI
support or with MPI stubs.


To compile with MPI:

1.  Install MPICH for Windows, specify the corresponding include and
    lib directories in MSVS/Tools/Options/Projects and Solutions/VC++
    Directories

2.  Compile LAMMPS using Debug or Release configurations from the
provided projects

To compile with MPI STUBS
   
1.  Compile the STUBS.vcproj 

2.  Compile LAMMPS using Debug_STUBS or Release_STUBS configurations
from the provided projects

To run the code you will need the mpich and fftw213 dlls accessible by
the system search (they may be copied to Windows/system32 directory).
The fftw213 ddlls may be found in vs9/extra/fftw213 or downloaded from
the fftw site

To include additional packages into LAMMPS projects, you may follow
the pattern of LAMMPS-std-user.vcproj:

1. Add the appropriate *.cpp files to the project (for example to USER-* filter)

2. Add corresponding .h files from the project directory to the
   style_* headers listed in the Settings filter of the project.  For
   example, if there is a pair_*.h file in the project directory, it
   should be added to settings/style_pair.h aggregate header.
   

  
