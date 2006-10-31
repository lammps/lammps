# Make.csh = update Makefile.lib and Makefile.list
# use current list of *.cpp and *.h files in src dir

if ($1 == "Makefile.lib") then

  set list = `ls *.cpp | sed s/^main\.cpp//`
  sed -i -e "s/SRC =\t.*/SRC =\t$list/" Makefile.lib
  set list = `ls *.h`
  sed -i -e "s/INC =\t.*/INC =\t$list/" Makefile.lib

else if ($1 == "Makefile.list") then

  set list = `ls *.cpp`
  sed -i -e "s/SRC =\t.*/SRC =\t$list/" Makefile.list
  set list = `ls *.h`
  sed -i -e "s/INC =\t.*/INC =\t$list/" Makefile.list

endif
