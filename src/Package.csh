# Package.csh = package management, called from Makefile
# Syntax: csh Package.csh DIR status/update/overwrite

# package installed if any *.cpp or *.h file is in src

# if last arg = "status":
#   print installation status of each package
#   flag package files that do not exist in src
#   flag package files that are different than src version

# if last arg = "update":
#   if package installed, copy its files that are different from package to src

# if last arg = "overwrite":
#   if package installed, copy its files that are different from src to package

# use diff to compare files
#   tried using cmp, but it doesn't satisfy if test if one file is
#   just longer than the other (has new stuff added)

set glob
set style = `echo $1 | sed 'y/-ABCDEFGHIJKLMNOPQRSTUVWXYZ/_abcdefghijklmnopqrstuvwxyz/'`
cd $1

set installed = 0
foreach file (*.cpp *.h)
  if (-e ../$file) then
    set installed = 1
  endif
end

if ($2 == "status") then

  if ($installed) then
    echo "Installed YES: package $1"
    foreach file (*.cpp *.h)
      if (! -e ../$file) then
        echo "  src/$file does not exist"
      else if (`diff --brief $file ../$file` != "") then
        echo "  src/$file and $1/$file are different"
      endif
    end
  else
    echo "Installed  NO: package $1"
  endif

else if ($2 == "update") then

  echo "Updating src files from $1 package files"

  if ($installed) then
    foreach file (*.cpp *.h)
      if (! -e ../$file) then
        echo "  creating src/$file"
        cp $file ..
      else if (`diff --brief $file ../$file` != "") then
        echo "  updating src/$file"
        cp $file ..
      endif
    end
  else
    echo "  $1 package is not installed, no action"
  endif

else if ($2 == "overwrite") then

  echo "Overwriting $1 package files with src files"

  if ($installed) then
    foreach file (*.cpp *.h)
      if (! -e ../$file) then
        echo "  src/$file does not exist"
      else if (`diff --brief $file ../$file` != "") then
        echo "  overwriting $1/$file"
        cp ../$file .
      endif
    end
  else
    echo "  $1 package is not installed, no action"
  endif

endif
