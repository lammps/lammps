# Make.csh = update Makefile.lib or Makefile.list or style_*.h files
# use current list of *.cpp and *.h files in src dir

if ($1 == "Makefile.lib") then

  set list1 = `ls [abcde]*.cpp`
  set list2 = `ls [fghij]*.cpp`
  set list3 = `ls [klmno]*.cpp | sed s/^main\.cpp//`
  set list4 = `ls [pqrst]*.cpp`
  set list5 = `ls [uvwxyz]*.cpp`
  sed -i -e "s/SRC =	.*/SRC =	$list1/" Makefile.lib
  sed -i -e "s/SRC =	\(.*\)/SRC =	\1 $list2/" Makefile.lib
  sed -i -e "s/SRC =	\(.*\)/SRC =	\1 $list3/" Makefile.lib
  sed -i -e "s/SRC =	\(.*\)/SRC =	\1 $list4/" Makefile.lib
  sed -i -e "s/SRC =	\(.*\)/SRC =	\1 $list5/" Makefile.lib

  set list1 = `ls [abcde]*.h`
  set list2 = `ls [fghij]*.h`
  set list3 = `ls [klmno]*.h`
  set list4 = `ls [pqrst]*.h`
  set list5 = `ls [uvwxyz]*.h`
  sed -i -e "s/INC =	.*/INC =	$list1/" Makefile.lib
  sed -i -e "s/INC =	\(.*\)/INC =	\1 $list2/" Makefile.lib
  sed -i -e "s/INC =	\(.*\)/INC =	\1 $list3/" Makefile.lib
  sed -i -e "s/INC =	\(.*\)/INC =	\1 $list4/" Makefile.lib
  sed -i -e "s/INC =	\(.*\)/INC =	\1 $list5/" Makefile.lib

else if ($1 == "Makefile.list") then

  set list1 = `ls [abcde]*.cpp`
  set list2 = `ls [fghij]*.cpp`
  set list3 = `ls [klmno]*.cpp`
  set list4 = `ls [pqrst]*.cpp`
  set list5 = `ls [uvwxyz]*.cpp`
  sed -i -e "s/SRC =	.*/SRC =	$list1/" Makefile.list
  sed -i -e "s/SRC =	\(.*\)/SRC =	\1 $list2/" Makefile.list
  sed -i -e "s/SRC =	\(.*\)/SRC =	\1 $list3/" Makefile.list
  sed -i -e "s/SRC =	\(.*\)/SRC =	\1 $list4/" Makefile.list
  sed -i -e "s/SRC =	\(.*\)/SRC =	\1 $list5/" Makefile.list

  set list1 = `ls [abcde]*.h`
  set list2 = `ls [fghij]*.h`
  set list3 = `ls [klmno]*.h`
  set list4 = `ls [pqrst]*.h`
  set list5 = `ls [uvwxyz]*.h`
  sed -i -e "s/INC =	.*/INC =	$list1/" Makefile.list
  sed -i -e "s/INC =	\(.*\)/INC =	\1 $list2/" Makefile.list
  sed -i -e "s/INC =	\(.*\)/INC =	\1 $list3/" Makefile.list
  sed -i -e "s/INC =	\(.*\)/INC =	\1 $list4/" Makefile.list
  sed -i -e "s/INC =	\(.*\)/INC =	\1 $list5/" Makefile.list

else if ($1 == "style") then

  set list = `grep -l ANGLE_CLASS angle_*.h`
  if (-e style_angle.tmp) then
    rm style_angle.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_angle.tmp
  end
  if (! -e style_angle.tmp) then
     rm style_angle.h
     touch style_angle.h
  else if (! -e style_angle.h) then
     mv style_angle.tmp style_angle.h
     rm Obj_*/force.d
  else if (`diff style_angle.h style_angle.tmp` != "") then
     mv style_angle.tmp style_angle.h
     rm Obj_*/force.d
  else
     rm style_angle.tmp
  endif

  set list = `grep -l ATOM_CLASS atom_vec_*.h`
  if (-e style_atom.tmp) then
    rm style_atom.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_atom.tmp
  end
  if (! -e style_atom.tmp) then
     rm style_atom.h
     touch style_atom.h
  else if (! -e style_atom.h) then
     mv style_atom.tmp style_atom.h
     rm Obj_*/atom.d
  else if (`diff style_atom.h style_atom.tmp` != "") then
     mv style_atom.tmp style_atom.h
     rm Obj_*/atom.d
  else
     rm style_atom.tmp
  endif

  set list = `grep -l BOND_CLASS bond_*.h`
  if (-e style_bond.tmp) then
    rm style_bond.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_bond.tmp
  end
  if (! -e style_bond.tmp) then
     rm style_bond.h
     touch style_bond.h
  else if (! -e style_bond.h) then
     mv style_bond.tmp style_bond.h
     rm Obj_*/force.d
  else if (`diff style_bond.h style_bond.tmp` != "") then
     mv style_bond.tmp style_bond.h
     rm Obj_*/force.d
  else
     rm style_bond.tmp
  endif

  set list = `grep -l COMMAND_CLASS *.h`
  if (-e style_command.tmp) then
    rm style_command.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_command.tmp
  end
  if (! -e style_command.tmp) then
     rm style_command.h
     touch style_command.h
  else if (! -e style_command.h) then
     mv style_command.tmp style_command.h
     rm Obj_*/input.d
  else if (`diff style_command.h style_command.tmp` != "") then
     mv style_command.tmp style_command.h
     rm Obj_*/input.d
  else
     rm style_command.tmp
  endif

  set list = `grep -l COMPUTE_CLASS compute_*.h`
  if (-e style_compute.tmp) then
    rm style_compute.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_compute.tmp
  end
  if (! -e style_compute.tmp) then
     rm style_compute.h
     touch style_compute.h
  else if (! -e style_compute.h) then
     mv style_compute.tmp style_compute.h
     rm Obj_*/modify.d
  else if (`diff style_compute.h style_compute.tmp` != "") then
     mv style_compute.tmp style_compute.h
     rm Obj_*/modify.d
  else
     rm style_compute.tmp
  endif

  set list = `grep -l DIHEDRAL_CLASS dihedral_*.h`
  if (-e style_dihedral.tmp) then
    rm style_dihedral.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_dihedral.tmp
  end
  if (! -e style_dihedral.tmp) then
     rm style_dihedral.h
     touch style_dihedral.h
  else if (! -e style_dihedral.h) then
     mv style_dihedral.tmp style_dihedral.h
     rm Obj_*/force.d
  else if (`diff style_dihedral.h style_dihedral.tmp` != "") then
     mv style_dihedral.tmp style_dihedral.h
     rm Obj_*/force.d
  else
     rm style_dihedral.tmp
  endif

  set list = `grep -l DUMP_CLASS dump_*.h`
  if (-e style_dump.tmp) then
    rm style_dump.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_dump.tmp
  end
  if (! -e style_dump.tmp) then
     rm style_dump.h
     touch style_dump.h
  else if (! -e style_dump.h) then
     mv style_dump.tmp style_dump.h
     rm Obj_*/output.d
  else if (`diff style_dump.h style_dump.tmp` != "") then
     mv style_dump.tmp style_dump.h
     rm Obj_*/output.d
  else
     rm style_dump.tmp
  endif

  set list = `grep -l FIX_CLASS fix_*.h`
  if (-e style_fix.tmp) then
    rm style_fix.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_fix.tmp
  end
  if (! -e style_fix.tmp) then
     rm style_fix.h
     touch style_fix.h
  else if (! -e style_fix.h) then
     mv style_fix.tmp style_fix.h
     rm Obj_*/modify.d
  else if (`diff style_fix.h style_fix.tmp` != "") then
     mv style_fix.tmp style_fix.h
     rm Obj_*/modify.d
  else
     rm style_fix.tmp
  endif

  set list = `grep -l IMPROPER_CLASS improper_*.h`
  if (-e style_improper.tmp) then
    rm style_improper.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_improper.tmp
  end
  if (! -e style_improper.tmp) then
     rm style_improper.h
     touch style_improper.h
  else if (! -e style_improper.h) then
     mv style_improper.tmp style_improper.h
     rm Obj_*/force.d
  else if (`diff style_improper.h style_improper.tmp` != "") then
     mv style_improper.tmp style_improper.h
     rm Obj_*/force.d
  else
     rm style_improper.tmp
  endif

  set list = `grep -l INTEGRATE_CLASS *.h`
  if (-e style_integrate.tmp) then
    rm style_integrate.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_integrate.tmp
  end
  if (! -e style_integrate.tmp) then
     rm style_integrate.h
     touch style_integrate.h
  else if (! -e style_integrate.h) then
     mv style_integrate.tmp style_integrate.h
     rm Obj_*/update.d
  else if (`diff style_integrate.h style_integrate.tmp` != "") then
     mv style_integrate.tmp style_integrate.h
     rm Obj_*/update.d
  else
     rm style_integrate.tmp
  endif

  set list = `grep -l KSPACE_CLASS *.h`
  if (-e style_kspace.tmp) then
    rm style_kspace.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_kspace.tmp
  end
  if (! -e style_kspace.tmp) then
     rm style_kspace.h
     touch style_kspace.h
  else if (! -e style_kspace.h) then
     mv style_kspace.tmp style_kspace.h
     rm Obj_*/force.d
  else if (`diff style_kspace.h style_kspace.tmp` != "") then
     mv style_kspace.tmp style_kspace.h
     rm Obj_*/force.d
  else
     rm style_kspace.tmp
  endif

  set list = `grep -l MINIMIZE_CLASS min_*.h`
  if (-e style_minimize.tmp) then
    rm style_minimize.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_minimize.tmp
  end
  if (! -e style_minimize.tmp) then
     rm style_minimize.h
     touch style_minimize.h
  else if (! -e style_minimize.h) then
     mv style_minimize.tmp style_minimize.h
     rm Obj_*/update.d
  else if (`diff style_minimize.h style_minimize.tmp` != "") then
     mv style_minimize.tmp style_minimize.h
     rm Obj_*/update.d
  else
     rm style_minimize.tmp
  endif

  set list = `grep -l PAIR_CLASS pair_*.h`
  if (-e style_pair.tmp) then
    rm style_pair.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_pair.tmp
  end
  if (! -e style_pair.tmp) then
     rm style_pair.h
     touch style_pair.h
  else if (! -e style_pair.h) then
     mv style_pair.tmp style_pair.h
     rm Obj_*/force.d
  else if (`diff style_pair.h style_pair.tmp` != "") then
     mv style_pair.tmp style_pair.h
     rm Obj_*/force.d
  else
     rm style_pair.tmp
  endif

  set list = `grep -l REGION_CLASS region_*.h`
  if (-e style_region.tmp) then
    rm style_region.tmp
  endif
  foreach file ($list)
    set qfile = \"$file\"
    echo "#include $qfile" >>! style_region.tmp
  end
  if (! -e style_region.tmp) then
     rm style_region.h
     touch style_region.h
  else if (! -e style_region.h) then
     mv style_region.tmp style_region.h
     rm Obj_*/domain.d
  else if (`diff style_region.h style_region.tmp` != "") then
     mv style_region.tmp style_region.h
     rm Obj_*/domain.d
  else
     rm style_region.tmp
  endif

endif
