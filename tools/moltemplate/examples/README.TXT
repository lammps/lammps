These are examples for the "moltemplate" molecule builder for LAMMPS.
http://www.moltemplate.org

Each directory contains one or more examples.

Each example directory contains:

   images/              This folder has pictures of the molecules in the system
   moltemplate_files/   This folder contains LT files and other auxiliary files
   README_setup.sh      Instructions for how to use moltemplate (executable)
   README_visualize.txt Instructions for viewing in DATA/DUMP files in VMD

     ...and one or more LAMMPS input scripts with names like

   run.in.min
   run.in.npt
   run.in.nvt         
                      
You can run these scripts using
      lmp_linux -i run.in.npt
(The name of your lammps binary, "lmp_linux" in this example, may vary. 
 Sometimes, these scripts must be run in a certain order.  For example
 it may be necessary to run run.in.min to minimize the system before you can use run.in.npt, and later run.in.nvt.  The README_run.sh file in each subdirectory
 specifies indicates the order.  These files have not been optimized.)
