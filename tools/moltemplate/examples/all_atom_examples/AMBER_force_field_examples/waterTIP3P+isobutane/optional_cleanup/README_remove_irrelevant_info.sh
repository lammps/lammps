  # MOST USERS CAN IGNORE THIS FILE
  #
  # Unfortunately, as of 2014-4-19, the system.data and system.in.settings file
  # which are created by moltemplate.sh contain a lot of irrelevant information,
  # such as definition of parameters for atom types not present in the current
  # system.  This extra information takes up about 1 MB.
  #
  # This appears to be harmless.
  # (Loading this extra information does not seem to slow down LAMMPS.)
  #
  # --------- OPTIONAL STEPS FOR STRIPPING OUT JUNK ---------
  #
  # However if you want to eliminate this junk from these files
  # For now, we can strip this out using ltemplify.py to build a new .lt file.
  #
  # I suggest you do this in a temporary_directory

  mkdir new_lt_file
  cd new_lt_file/

  # now run ltemplify.py

  ltemplify.py ../system.in.init ../system.in.settings ../system.data > system.lt
  rm -rf ../system.data ../system.in*  # these old lammps files no longer needed

  # This creates a new .LT file named "system.lt" in the local directory.
  # Unfortunately, it may be missing some information because ltemplify.py
  # does not understand all the commands present in a LAMMPS input script.
  # If you define groups or use constraints, you must define them again. In this
  # case, we must add the SHAKE constraint for the "TIP3P_2004" water molecule.
  # So we have to remember the original name of the bond types and angle types.
  # (For this example, SHAKE is applied to the water molecule, which is defined
  #  in "tip3p_2004.lt" file in the "common/" directory.  Check this file.)

  ATOMTYPENUM_ow=`awk '{if ($1 == "@/atom:TIP3P_2004/ow") print $2}' < ../moltemplate_files/output_ttree/ttree_assignments.txt`
  ATOMTYPENUM_hw=`awk '{if ($1 == "@/atom:TIP3P_2004/hw") print $2}' < ../moltemplate_files/output_ttree/ttree_assignments.txt`
  BONDTYPENUM=`awk '{if ($1 == "@/bond:TIP3P_2004/OH") print $2}' < ../moltemplate_files/output_ttree/ttree_assignments.txt`
  ANGLETYPENUM=`awk '{if ($1 == "@/angle:TIP3P_2004/HOH") print $2}' < ../moltemplate_files/output_ttree/ttree_assignments.txt`
  echo "" >> system.lt
  echo "write_once(\"In Settings\") {" >> system.lt
  echo "  group tip3p type  @atom:type$ATOMTYPENUM_ow  @atom:type$ATOMTYPENUM_hw" >> system.lt
  echo "  fix fShakeTIP3P tip3p shake 0.0001 10 100 b @bond:type$BONDTYPENUM a @angle:type$ANGLETYPENUM" >> system.lt
  echo "}" >> system.lt
  echo "" >> system.lt

  # The ltemplify.py script also does not copy the boundary dimensions.
  # We must do this manually.
  # If you did NOT throw away the "Data Boundary" file usually located in
  # "moltemplate_files/output_ttree/Data Boundary"
  # then you can copy that information from this file into system.lt
  echo "write_once(\"Data Boundary\") {" >> system.lt
  cat "../moltemplate_files/output_ttree/Data Boundary" >> system.lt
  echo "}" >> system.lt
  echo "" >> system.lt

  # Now, run moltemplate on this new .LT file.
  moltemplate.sh system.lt
  # This will create: "system.data" "system.in.init" "system.in.settings."

  # That's it.  The new "system.data" and system.in.* files should
  # be ready to run in LAMMPS.

  # Now copy the system.data and system.in.* files to the place where
  # you plan to run moltemplate
  mv -f system.data system.in.* ../
  cd ../

  # Now delete all of the temporary files we generated
  rm -rf new_lt_file/
