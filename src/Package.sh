# Package.sh = package management, called from Makefile
# Syntax: sh Package.sh DIR status/update/overwrite/diff

# enforce using portable C locale
LC_ALL=C
export LC_ALL

# package is already installed if any package *.cpp or *.h file is in src
# else not installed

cd $1

installed=0
for file in *.cpp *.h; do
  if (test -e ../$file) then
    installed=1
  fi
done

# status, only if installed
# issue warning for any package file not in src or that is different

if (test $2 = "status") then
  if (test $installed = 1) then
    echo "Installed YES: package $1"
    for file in *.cpp *.h; do
      if (test ! -e ../$file) then
        echo "  src/$file does not exist"
      elif (! cmp -s $file ../$file) then
        echo "  src/$file and $1/$file are different"
      fi
    done
  else
    echo "Installed  NO: package $1"
  fi

# installed, list only if installed

elif (test $2 = "installed") then
  if (test $installed = 1) then
    echo "Installed YES: package $1"
  fi

# update, only if installed
# perform a re-install, but only if the package is already installed

elif (test $2 = "update") then
  echo "Updating src files from $1 package files"
  if (test $installed = 1) then
    echo "  updating package $1"
    if (test -e Install.sh) then
      /bin/sh Install.sh 2
    else
      /bin/sh ../Install.sh 2
    fi
    cd ..
    /bin/sh Depend.sh $1
  else
    echo "  $1 package is not installed"
  fi

# overwrite, only if installed
# overwrite package file with src file, if the two are different

elif (test $2 = "overwrite") then
  echo "Overwriting $1 package files with src files"
  if (test $installed = 1) then
    for file in *.cpp *.h; do
      if (test ! -e ../$file) then
        continue
      elif (! cmp -s $file ../$file) then
        echo "  overwriting $1/$file"
        cp ../$file .
      fi
    done
  else
    echo "  $1 package is not installed"
  fi

# diff
# if installed:
# show any differences between src files and package files

elif (test $2 = "diff") then
  if (test $installed = 1) then
    echo "Installed YES: package $1"
    for file in *.cpp *.h; do
      if (test ! -e ../$file) then
        echo "  src/$file does not exist"
      elif (! cmp -s $file ../$file) then
        echo "************************************************"
        echo "diff -u $1/$file src/$file "
        echo "************************************************"
	diff -u $file  ../$file 
      fi
    done
  fi
fi
