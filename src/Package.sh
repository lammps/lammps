# Package.sh = package management, called from Makefile
# Syntax: sh Package.sh DIR status/update/overwrite/diff

# style used to translate dir name to package name

style=`echo $1 | sed 'y/-ABCDEFGHIJKLMNOPQRSTUVWXYZ/_abcdefghijklmnopqrstuvwxyz/'`

# package is already installed if any package *.cpp or *.h file is in src
# else not installed

cd $1

installed=0
for file in *.cpp *.h; do
  if (test -e ../$file) then
    installed=1
  fi
done

# status
# if installed:
# issue warning if any package file is not in src or is different

if (test $2 = "status") then
  if (test $installed = 1) then
    echo "Installed YES: package $1"
    for file in *.cpp *.h; do
      if (test ! -e ../$file) then
        echo "  src/$file does not exist"
      elif (test "`diff --brief $file ../$file`" != "") then
        echo "  src/$file and $1/$file are different"
      fi
    done
  else
    echo "Installed  NO: package $1"
  fi

# update
# if installed:
# cp package file to src if doesn't exist or is different

elif (test $2 = "update") then
  echo "Updating src files from $1 package files"
  if (test $installed = 1) then
    if (test ! -e Package.sh) then
      for file in *.cpp *.h; do
        if (test ! -e ../$file) then
          echo "  creating src/$file"
          cp $file ..
        elif (test "`diff --brief $file ../$file`" != "") then
          echo "  updating src/$file"
          cp $file ..
        fi
      done
    else
      /bin/sh Package.sh
    fi
  else
    echo "  $1 package is not installed, no action"
  fi

# overwrite
# if installed:
# if package file not in src, issue warning
# if src file different than package file, overwrite package file

elif (test $2 = "overwrite") then
  echo "Overwriting $1 package files with src files"
  if (test $installed = 1) then
    for file in *.cpp *.h; do
      if (test ! -e ../$file) then
        echo "  src/$file does not exist"
      elif (test "`diff --brief $file ../$file`" != "") then
        echo "  overwriting $1/$file"
        cp ../$file .
      fi
    done
  else
    echo "  $1 package is not installed, no action"
  fi

# diff
# if installed:
# show any differences between src files and package files

elif (test $2 = "diff") then
  if (test $installed = 1) then
    echo "Installed YES: package $1"
    for file in *.cpp *.h; do
      if (test ! -e ../$file) then
        echo "************************************************************************"
        echo "  src/$file does not exist"
        echo "************************************************************************"
      elif (test "`diff --brief $file ../$file`" != "") then
        echo "************************************************************************"
        echo "diff $1/$file src/$file "
        echo "************************************************************************"
	diff $file  ../$file 
      fi
    done
  fi
fi



