# -------- REQUIREMENTS: ---------
#  You must define your MOLTEMPLATE_PATH environment variable
#  and set it to the "common" subdirectory of your moltemplate distribution.
#  (See the "Installation" section in the moltemplate manual.)


# Create LAMMPS input files this way:
cd moltemplate_files

  # run moltemplate

  moltemplate.sh system.lt

  # This will generate various files with names ending in *.in* and *.data. 

  # FIX THE PAIR STYLES
  #    (Sorry, this is messy)
  #
  # I was forced to change the default pair-style for AMBER-force-fields (GAFF)
  # from lj/charmm/coul/long to lj/charmm/coul/charmm.  (This is because
  # LAMMPS crashes when using lj/charmm/coul/long on a system without any
  # charged particles, and users were complaining it was moltemplate's fault.
  # I wish LAMMPS would not do this.)
  #
  # Unfortunately, this means that the "Isobutane" molecule (which uses
  # AMBER's GAFF), and the "TIP3P_2004" molecule now use different pair styles.
  #
  # The cleanest way to fix this is to force the two molecules to use
  # the same pair style.
  # (Using a hybrid pair_style is not practical because that disables mixing
  #  rules.  This would force us to add a huge list of pair_coeff commands to
  #  explain how TIP3P_2004 atoms interact with all of the various GAFF atoms.)

  # Add a line to systems.in.init to override the pair_style.
  # Change the pair_style to "lj/charmm/coul/long 10.0 10.5 10.5".

  echo "pair_style   hybrid lj/charmm/coul/long 10.0 10.5 10.5" >> system.in.init

  # Then use "sed" to replace "lj/charmm/coul/charmm" with "lj/charmm/coul/long"
  sed -i 's/lj\/charmm\/coul\/charmm/lj\/charmm\/coul\/long/g' system.in.settings

  # These files are the input files directly read by LAMMPS.  Move them to 
  # the parent directory (or wherever you plan to run the simulation).
  #cp -f system.data system.in* ../


  # --------- OPTIONAL STEPS FOR STRIPPING OUT JUNK ---------
  # --------- edit 2013-8-28 ---------
  echo "-----------------------------------------------------------------" >&2
  echo "OPTIONAL STEP: PRUNING THE RESULTING MOLTEMPLATE OUTPUT TO" >&2
  echo "               INCLUDE ONLY ATOMS AND TYPES WE ARE ACTUALLY USING." >&2
  # Unfortunately, as of 2013-8-28, these files contain a lot of irrelevant
  # information (for atom types not present in the current system).
  # For now, we can strip this out using ltemplify.py to build a new .lt file.
  # THIS IS AN UGLY WORKAROUND. HOPEFULLY IN THE FUTURE, WE CAN SKIP THESE STEPS

  # do this in a temporary_directory
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
  ATOMTYPENUM_ow=`awk '{if ($1 == "@/atom:TIP3P_2004/ow") print $2}' < ../output_ttree/ttree_assignments.txt`
  ATOMTYPENUM_hw=`awk '{if ($1 == "@/atom:TIP3P_2004/hw") print $2}' < ../output_ttree/ttree_assignments.txt`
  BONDTYPENUM=`awk '{if ($1 == "@/bond:TIP3P_2004/OH") print $2}' < ../output_ttree/ttree_assignments.txt`
  ANGLETYPENUM=`awk '{if ($1 == "@/angle:TIP3P_2004/HOH") print $2}' < ../output_ttree/ttree_assignments.txt`
  echo "" >> system.lt
  echo "write_once(\"In Settings\") {" >> system.lt
  echo "  group tip3p type  @atom:type$ATOMTYPENUM_ow  @atom:type$ATOMTYPENUM_hw" >> system.lt
  echo "  fix fShakeTIP3P tip3p shake 0.0001 10 100 b @bond:type$BONDTYPENUM a @angle:type$ANGLETYPENUM" >> system.lt
  echo "}" >> system.lt
  echo "" >> system.lt

  # The ltemplify.py script also does not copy the boundary dimensions.
  # We must do this manually as well.
  echo "write_once(\"Data Boundary\") {" >> system.lt
  cat "../output_ttree/Data Boundary" >> system.lt
  echo "}" >> system.lt
  echo "" >> system.lt
  # Now, run moltemplate on this new .LT file.
  moltemplate.sh system.lt
  # This will create: "system.data" "system.in.init" "system.in.settings."

  # move the final DATA and INput scripts to the desired location,
  mv -f system.data system.in* ../../

  # and clean up the mess
  rm -rf output_ttree/
  cd ..
  rm -rf new_lt_file/
  echo "---------------- DONE PRUNING MOLTEMPLATE OUTPUT ----------------" >&2
  echo "-----------------------------------------------------------------" >&2
  # --------- END OF OPTIONAL STEPS FOR STRIPPING OUT JUNK ---------




  # Optional:
  # The "./output_ttree/" directory is full of temporary files generated by 
  # moltemplate.  They can be useful for debugging, but are usually thrown away.
  rm -rf output_ttree/



cd ../
