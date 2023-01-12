# Make.sh = update Makefile.lib, Makefile.shlib, Makefile.list
#           or style_*.h files
# Syntax: sh Make.sh style
#         sh Make.sh Makefile.lib
#         sh Make.sh Makefile.shlib
#         sh Make.sh Makefile.list

# turn off enforced customizations
GREP_OPTIONS=
# enforce using portable C locale
LC_ALL=C
export LC_ALL GREP_OPTIONS

# function to create one style_*.h file

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
    fi
  elif (test ! -e style_$3.h) then
    mv style_$3.tmp style_$3.h
  elif (test "`diff --brief style_$3.h style_$3.tmp`" != "") then
    mv style_$3.tmp style_$3.h
  else
    rm -f style_$3.tmp
  fi
}

packages () {
  list=`grep -sl $1 */$2*.h`
  if (test -e packages_$3.tmp) then
    rm -f packages_$3.tmp
  fi
  for file in $list; do
    dir="\"`dirname $file`\""
    echo "#undef PACKAGE" >> packages_$3.tmp
    echo "#define PACKAGE $dir" >> packages_$3.tmp
    qfile="\"$file\""
    echo "#include $qfile" >> packages_$3.tmp
  done
  if (test ! -e packages_$3.tmp) then
    if (test ! -e packages_$3.h) then
      touch packages_$3.h
    elif (test "`cat packages_$3.h`" != "") then
      rm -f packages_$3.h
      touch packages_$3.h
    fi
  elif (test ! -e packages_$3.h) then
    mv packages_$3.tmp packages_$3.h
  elif (test "`diff --brief packages_$3.h packages_$3.tmp`" != "") then
    mv packages_$3.tmp packages_$3.h
  else
    rm -f packages_$3.tmp
  fi
}

# create individual style files
# called by "make machine"
# col 1 = string to search for
# col 2 = search in *.h files starting with this name
# col 3 = name of style file
# col 4 = file that includes the style file
# col 5 = optional 2nd file that includes the style file

cmd=$1

if (test $cmd = "style") || (test $cmd = "packages") then

  $cmd ANGLE_CLASS     angle_      angle      force
  $cmd ATOM_CLASS      atom_vec_   atom       atom      atom_vec_hybrid
  $cmd BODY_CLASS      body_       body       atom_vec_body
  $cmd BOND_CLASS      bond_       bond       force
  $cmd COMMAND_CLASS   ""          command    input
  $cmd COMPUTE_CLASS   compute_    compute    modify
  $cmd DIHEDRAL_CLASS  dihedral_   dihedral   force
  $cmd DUMP_CLASS      dump_       dump       output    write_dump
  $cmd FIX_CLASS       fix_        fix        modify
  $cmd IMPROPER_CLASS  improper_   improper   force
  $cmd INTEGRATE_CLASS ""          integrate  update
  $cmd KSPACE_CLASS    ""          kspace     force
  $cmd MINIMIZE_CLASS  min_        minimize   update
  $cmd NBIN_CLASS      nbin_       nbin       neighbor
  $cmd NPAIR_CLASS     npair_      npair      neighbor
  $cmd NSTENCIL_CLASS  nstencil_   nstencil   neighbor
  $cmd NTOPO_CLASS     ntopo_      ntopo      neighbor
  $cmd PAIR_CLASS      pair_       pair       force
  $cmd READER_CLASS    reader_     reader     read_dump
  $cmd REGION_CLASS    region_     region     domain

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
