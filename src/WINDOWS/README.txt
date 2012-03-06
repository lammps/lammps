INSTRUCTIONS FOR COMPILING LAMMPS WITH VISUAL STUDIO 2005

The provided project

LAMMPS.vcproj                  

includes the minimal package set: KSPACE, MANYBODY, MOLECULE.

The package set may be reconfiured with the help of the supplied VS
macro (see below).

The project has configurations to compile either with MPI support or
with MPI stubs.

To compile with MPI:

1.  Install MPICH for Windows, specify the corresponding include and
    lib directories in MSVS/Tools/Options/Projects and Solutions/VC++
    Directories

2.  Compile LAMMPS using Debug or Release configurations from the
    provided projects

To compile with MPI STUBS
   
1.  Compile STUBS.vcproj 

2.  Compile LAMMPS using Debug_STUBS or Release_STUBS configurations
from the provided project

To run the code you may need mpich and fftw213 dlls accessible by the
system search (they may be copied to Windows/system32 directory).  The
fftw213 dlls can be found in vs9/extra/fftw213 or downloaded from the
fftw site

To customise the packages via a Visual Basic macro:

1. Load LAMMPS solution in Visual Studio IDE
2. Select in the main menu "Tools/Macros/Load Macro Project..."
   and load the file src/WINDOWS/LAMMPS.vsmacros
3. In the "Macro Explorer" on the right panel open LAMMPS and LAMMPS_settings
4. Double click on "ManagePackages" to run the configuration
   macro. Please note that the window for running macro sometimes
   opens in the background, so use Alt-TAB to locate it.
5. Configure a custom set of packages and press Ok. Wait till the
   macro completes.
6. Rebuild the LAMMPS project

Before the first build or after an update from LAMMPS src repository
it is recommended to run "ManagePackages" macro an check "Refresh file
list in src filter" to get an up to date list of source files in the
"src" project filter. This may be needed as the file composition in
src may change between LAMMPS releases.

Some of the packages were not tested to be compatible with VS compiler
or require additional libraries. They are marked with asterics in the
package list displayed when the macro is running. If you wish to try
those packages, install them using the macro and then change the
project properties (libraries, includes, etc.) manually.

Please note that "ManagePackages" macro works on the project named
LAMMPS.  So if you rename the LAMMPS project and still want to use
automatic package configuration, please specify the new name in the
line "Dim LAMMPS_project_name As String =" at the beginning of the
macro code.

The default package options such as the path to include and library
files, description, etc can also be changed by editing the
"ManagePackages" macro code. To do this right click on
"ManagePackages" in the "Macro Explorer" and select Edit. Then go to
the section named 

"===================== Custom Package options ========================",

find the required package and change its properties by modyfing the
corresponding PKG_OPTS(...) entry.


