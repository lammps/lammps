# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

# this is default Install.sh for all packages
# if package has an auxiliary library or a file with a dependency,
# then package dir has its own customized Install.sh

mode=$1

# arg1 = file, arg2 = file it depends on

action () {
  currdir=${PWD##*/} 
  if (test $mode = 0) then #uninstall    
    if (test -e ../backup/$currdir/$1) then
		mv ../backup/$currdir/$1 ..
	elif (cmp -s $1 ../$1) then
		rm -f ../$1
	fi
  elif (! cmp -s $1 ../$1) then #install
	if (test -e ../$1) then
		if (test ! -e ../backup/$currdir/$1) then
			mv ../$1 ../backup/$currdir/
		fi
	fi
	cp $1 ..
	if (test $mode = 2) then
        echo "  updating src/$1"
    fi
  elif (test -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}


# all package files with no dependencies
currdir=${PWD##*/}       
#creating backup folders if they aren't already there
if (! test -d ../backup) then
   mkdir ../backup
fi 
if (! test -d ../backup/$currdir) then
   mkdir ../backup/$currdir
fi
for file in *.cpp *.h; do
  action $file
done
