# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# enforce using portable C locale
LC_ALL=C
export LC_ALL

cat <<EOF
WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING
WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING

 Support for building the ML-POD package with the legacy build system using GNU
 make will be removed in Summer 2025.  Please switch to using CMake to build
 LAMMPS as soon as possible and report any problems to developers@lammps.org
 or post a bug report issue at https://github.com/lammps/lammps/issues

WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING
WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING
EOF

# arg1 = file, arg2 = file it depends on

action () {
    if (test $mode = 0) then
        rm -f ../$1
    elif (! cmp -s $1 ../$1) then
        if (test -z "$2" || test -e ../$2) then
            cp $1 ..
            if (test $mode = 2) then
                echo "  updating src/$1"
            fi
        fi
    elif (test -n "$2") then
        if (test ! -e ../$2) then
            rm -f ../$1
        fi
    fi
}

# all package files with no dependencies

for file in *.cpp *.h; do
    test -f ${file} && action $file
done

if (test $1 = 1) then
    if (test -e ../Makefile.package) then
        sed -i -e 's/[^ \t]*mlpod[^ \t]* //g' ../Makefile.package
        sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/mlpod |' ../Makefile.package
        sed -i -e 's|^PKG_LIB =[ \t]*|&-lmlpod |' ../Makefile.package
        sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(mlpod_SYSPATH) |' ../Makefile.package
        sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(mlpod_SYSLIB) |' ../Makefile.package
    fi
    if (test -e ../Makefile.package.settings) then
        sed -i -e '/^[ \t]*include.*mlpod.*$/d' ../Makefile.package.settings
        # multiline form needed for BSD sed on Macs
        sed -i -e '4 i \
include ..\/..\/lib\/mlpod\/Makefile.lammps
' ../Makefile.package.settings
    fi

elif (test $1 = 0) then
    if (test -e ../Makefile.package) then
        sed -i -e 's/[^ \t]*mlpod[^ \t]* //g' ../Makefile.package
    fi
    if (test -e ../Makefile.package.settings) then
        sed -i -e '/^[ \t]*include.*mlpod.*$/d' ../Makefile.package.settings
    fi
fi

