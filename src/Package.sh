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

# status, only if installed
# issue warning for any package file not in src or that is different

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

# update, only if installed
# if package dir has its own Package.sh, use it
# cp package file to src if it exists and is different
# set installflag if any package file is not in src and do full install
# this is needed when a patch has added a new file to the package

elif (test $2 = "update") then
  echo "Updating src files from $1 package files"
  if (test $installed = 1) then
    if (test ! -e Package.sh) then
      installflag=0
      for file in *.cpp *.h; do
        if (test ! -e ../$file) then
          installflag=1
        elif (test "`diff --brief $file ../$file`" != "") then
          echo "  updating src/$file"
          cp $file ..
        fi
      done
      if (test $installflag = 1) then
         echo "  reinstalling package $1"
	 /bin/sh Install.sh 1
         /bin/sh ../Depend.sh $1 1
      fi
    else
      /bin/sh Package.sh
    fi
  else
    echo "  $1 package is not installed, no action"
  fi

# overwrite, only if installed
# overwrite package file with src file, if the two are different

elif (test $2 = "overwrite") then
  echo "Overwriting $1 package files with src files"
  if (test $installed = 1) then
    for file in *.cpp *.h; do
      if (test ! -e ../$file) then
        continue
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



