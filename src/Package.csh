# Package.csh = package management, called from Makefile
# Syntax: csh Package.csh DIR status/update/overwrite

# if last arg = "status":
#   print installation status of each package
#   if package not installed (0-length src/style file), do nothing besides
#     check that no package files are in src (except style file)
#   flag src files that do not exist or are 0-length
#   list package files that are different than src version

# if last arg = "update":
#   if 0-length src/style file doesn't exist, create src/style file
#     via touch command, since LAMMPS needs it to do a build,
#     but don't copy any more files since don't want to install package
#   if 0-length src/style file already exists, do nothing
#   if non-0-length src/style file exists,
#     copy any package file into src if file doesn't exist or is different

# if last arg = "overwrite":
#   if package not installed (0-length src/style file), do nothing besides
#     check that no package files are in src (except style file)
#   flag src files that do not exist or are 0-length
#   overwrite package files that are different than src version

# use diff to compare files
#   tried using cmp, but it doesn't satisfy if test if one file is
#   just longer than the other (has new stuff added)

set glob
set style = `echo $1 | sed 'y/-ABCDEFGHIJKLMNOPQRSTUVWXYZ/_abcdefghijklmnopqrstuvwxyz/'`
cd $1

if ($2 == "status") then

  if (-z ../style_$style.h) then
    echo "Installed  NO: package $1"
    foreach file (*.cpp *.h)
      if (-e ../$file && $file != "style_$style.h") then
        echo "  src/$file exists but should not"
      endif
    end
  else
    echo "Installed YES: package $1"
    foreach file (*.cpp *.h)
      if (! -e ../$file) then
        echo "  src/$file does not exist"
      else if (-z ../$file) then
        echo "  src/$file is empty file"
      else if (`diff --brief $file ../$file` != "") then
        echo "  src/$file and $1/$file are different"
      endif
    end
  endif

else if ($2 == "update") then

  echo "Updating src from $1 package"

  if (! -e ../style_$style.h) then
    echo "  $1 package is not installed, but creating dummy style file"
    touch ../style_$style.h
  else if (-z ../style_$style.h) then
    echo "  $1 package is not installed, no action"
  else
    foreach file (*.cpp *.h)
      if (! -e ../$file) then
        echo "  updating src/$file"
        cp $file ..
      else if (`diff --brief $file ../$file` != "") then
        echo "  updating src/$file"
        cp $file ..
      endif
    end
  endif

else if ($2 == "overwrite") then

  echo "Overwriting $1 package with src"

  if (-z ../style_$style.h) then
    echo "  $1 package is not installed, no action"
    foreach file (*.cpp *.h)
      if (-e ../$file && $file != "style_$style.h") then
        echo "  src/$file exists but should not"
      endif
    end
  else
    foreach file (*.cpp *.h)
      if (! -e ../$file) then
        echo "  src/$file does not exist"
      else if (-z ../$file) then
        echo "  src/$file is empty file"
      else if (`diff --brief $file ../$file` != "") then
        echo "  overwriting $1/$file"
        cp ../$file .
      endif
    end
  endif

endif
