#!/usr/bin/env bash

  # cleanup_moltemplate.sh
  # This script attempts to remove irrelevant information from LAMMPS 
  # input scripts and data files (such as extra, unneeded atom types 
  # and force field parameters).
  #
  # Unfortunately, the system.data and system.in.settings file which are
  # created by moltemplate.sh often contain a lot of irrelevant information,
  # such as definition of parameters for atom types defined in some force field
  # file that the user is using, but not present in the system they are building
  #
  # In my experience, this extra information appears to be mostly harmless.
  # (Loading this information does not seem to slow down LAMMPS significantly.)
  # 
  # However it can make visualization difficult in VMD.  (Technically, this
  # extra information can take up megabytes of space in the system.data
  # and system.in.settings files.  Additionally, when you run LAMMPS, an O(n^2)
  # sized table is allocated to store the parameters for interactions between 
  # every possible pair of atom types (n atom types), and this occupies
  # significantly more memory if n is large.  For example, the "oplsaa.lt" file
  # and "oplsaa.prm" (TINKER-format) file both define almost 1000 atom types.)
  #
  # Usage: Invoke this script with no arguments, from a directory
  #        containing these files:
  #          system.data, system.in.init, system.in.settings, system.in.charges
  #        It will modify these files to remove unnecessary atoms and 
  #        parameters.  (If your files have other names, you must rename 
  #        them to match moltemplate file name conventions.)
  #
  # DO NOT USE THIS SCRIPT ON SIMULATIONS CONTAINING MANY-BODY PAIR STYLES,
  # DREIDING-STYLE HYDROGEN BONDS, OR SIMS NEEDING NON-STANDARD AUXILIARY FILES.
  # (This script relies on ltemplify.py and inherits its limitations.)

  PATH_TO_DATA_FILE="."

  pushd "$PATH_TO_DATA_FILE"

  mkdir new_lt_file_TMP
  cd new_lt_file_TMP

  # now run ltemplify.py

  ltemplify.py ../system.in.* ../system.data > system.lt

  # This creates a new .LT file named "system.lt" in the local directory.

  # The ltemplify.py script also does not copy the boundary dimensions.
  # We must do this manually.
  # Extract the header of the data file, reverse the order, and read lines
  # until you have 
  # If you did NOT throw away the "Data Boundary" file usually located in
  # "moltemplate_files/output_ttree/Data Boundary"
  # then you can copy that information from this file into system.lt


  # oops. looks like we don't need this after all
  #function _reverse_lines {
  #    # The following function reverses the order of lines in a file:
  #    # (Neither "tac", nor "tail -r" have cross-platform support.)
  #    awk '{a[i++]=$0} END {for (j=i-1; j>=0;) print a[j--] }'
  #}

  echo "" >> system.lt
  echo "write_once(\"Data Boundary\") {" >> system.lt
  # Extract the periodic boundary box dimensions from the 
  # end of the header section of the LAMMPS data file:
  extract_lammps_data.py Header < ../system.data | awk '{if (($3=="xlo") && ($4=="xhi")) {xl=$0} if (($3=="ylo") && ($4=="yhi")) {yl=$0} if (($3=="zlo") && ($4=="zhi")) {zl=$0} if (($4=="xy") && ($5=="xz") && ($6=="yz")) {xtr=$0}} END{print xl; print yl; print zl; if (xtr!="") {print xtr}}' >> system.lt
  echo "}" >> system.lt
  echo "" >> system.lt

  # Now, run moltemplate on this new .LT file.
  moltemplate.sh system.lt
  # This will create: "system.data" "system.in.init" "system.in.settings."

  # That's it.  The new "system.data" and system.in.settings files should
  # be ready to run in LAMMPS.

  # Special case: "set" commands
  # Typically "set type" or "set atom" commands are used to assign atom charge
  # If there is a "system.in.charges" file, then it contains these commands
  # however the atom type numbers will be wrong, so we must rewrite it.
  # Replace it with the corresponding commands from the system.in.settings
  # (whose atom type numbers are correct)
  if [ -f "../system.in.charges" ]; then
    awk '{ if ((NF >= 5) && ($1 == "set") && ($4 == "charge")) print $0}' \
        < system.in.settings > system.in.charges
    # There is no need to remove these lines from "system.in.settings",
    # (because there's no harm to invoke the "set" command twice)
    # ...but if you want to do that, try using a command similar to:
    #sed '/set type/,+1 d' < system.in.settings > system.in.settings.tmp
    #mv -f system.in.settings.tmp system.in.settings
  fi


  # Now move the system.data and system.in.* files to their original location:
  mv -f system.data system.in.* ../
  cd ../

  # Finally, delete all of the temporary files we generated
  rm -rf new_lt_file_TMP
  popd

