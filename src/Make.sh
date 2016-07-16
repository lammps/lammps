# Make.sh = update Makefile.lib, Makefile.shlib, Makefile.list 
#           or style_*.h files
# Syntax: sh Make.sh style
#         sh Make.sh Makefile.lib
#         sh Make.sh Makefile.shlib
#         sh Make.sh Makefile.list

# function to create one style_*.h file
# must whack *.d files that depend on style_*.h file,
# else Make will not recreate them

style () {
  list=`grep -sl $1 $2*.h`
  if (test -e style_$3.tmp) then
    rm -f style_$3.tmp
  fi
  for file in $list; do
    qfile="\"$file\""
    echo "#include $qfile" >> style_$3.tmp
  done
  if (test ! -e style_$3.tmp) then
    if (test ! -e style_$3.h) then
      touch style_$3.h
    elif (test "`cat style_$3.h`" != "") then
      rm -f style_$3.h
      touch style_$3.h
      rm -f Obj_*/$4.d
      if (test $5) then 
        rm -f Obj_*/$5.d
      fi
      rm -f Obj_*/lammps.d
    fi
  elif (test ! -e style_$3.h) then
    mv style_$3.tmp style_$3.h
    rm -f Obj_*/$4.d
    if (test $5) then 
      rm -f Obj_*/$5.d
    fi
    rm -f Obj_*/lammps.d
  elif (test "`diff --brief style_$3.h style_$3.tmp`" != "") then
    mv style_$3.tmp style_$3.h
    rm -f Obj_*/$4.d
    if (test $5) then 
      rm -f Obj_*/$5.d
    fi
    rm -f Obj_*/lammps.d
  else
    rm -f style_$3.tmp
  fi
}

# create individual style files
# called by "make machine"
# col 1 = string to search for
# col 2 = search in *.h files starting with this name
# col 3 = prefix of style file
# col 4 

if (test $1 = "style") then

  style ANGLE_CLASS     angle_      angle      force
  style ATOM_CLASS      atom_vec_   atom       atom      atom_vec_hybrid
  style BODY_CLASS      body_       body       atom_vec_body
  style BOND_CLASS      bond_       bond       force
  style COMMAND_CLASS   ""          command    input
  style COMPUTE_CLASS   compute_    compute    modify    modify_cuda
  style DIHEDRAL_CLASS  dihedral_   dihedral   force
  style DUMP_CLASS      dump_       dump       output    write_dump
  style FIX_CLASS       fix_        fix        modify
  style IMPROPER_CLASS  improper_   improper   force
  style INTEGRATE_CLASS ""          integrate  update
  style KSPACE_CLASS    ""          kspace     force
  style MINIMIZE_CLASS  min_        minimize   update
  style PAIR_CLASS      pair_       pair       force
  style READER_CLASS    reader_     reader     read_dump
  style REGION_CLASS    region_     region     domain

# edit Makefile.lib, for creating non-shared lib
# called by "make makelib"
# use current list of *.cpp and *.h files in src dir w/out main.cpp

elif (test $1 = "Makefile.lib") then

  list=`ls -1 *.cpp | sed s/^main\.cpp// | tr "[:cntrl:]" " "`
  sed -i -e "s/SRC =	.*/SRC =	$list/" Makefile.lib
  list=`ls -1 *.h | tr "[:cntrl:]" " "`
  sed -i -e "s/INC =	.*/INC =	$list/" Makefile.lib

# edit Makefile.shlib, for creating shared lib
# called by "make makeshlib"
# use current list of *.cpp and *.h files in src dir w/out main.cpp

elif (test $1 = "Makefile.shlib") then

  list=`ls -1 *.cpp | sed s/^main\.cpp// | tr "[:cntrl:]" " "`
  sed -i -e "s/SRC =	.*/SRC =	$list/" Makefile.shlib
  list=`ls -1 *.h | tr "[:cntrl:]" " "`
  sed -i -e "s/INC =	.*/INC =	$list/" Makefile.shlib

# edit Makefile.list
# called by "make makelist"
# use current list of *.cpp and *.h files in src dir

elif (test $1 = "Makefile.list") then

  list=`ls -1 *.cpp | tr "[:cntrl:]" " "`
  sed -i -e "s/SRC =	.*/SRC =	$list/" Makefile.list
  list=`ls -1 *.h | tr "[:cntrl:]" " "`
  sed -i -e "s/INC =	.*/INC =	$list/" Makefile.list

fi
