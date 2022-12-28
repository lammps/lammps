#!/bin/bash

has_f2c=$(type f2c > /dev/null 2>&1 && echo 1 || echo 0)
if test ${has_f2c} -eq 0
then
    echo "Must have f2c installed to run this script"
    exit 1
fi

# cleanup
rm -f *.c *.cpp *.P *~ *.orig *.bak *.rej

# translate files directly, skip those for which we have replacements.
for f in fortran/*.f
do \
    b=$(basename $f .f)
    if test $b == dgetrf2 || test $b == disnan || test $b == dlaisnan ||    \
            test $b == dlamch || test $b == dlarft || test $b == dpotrf2 || \
            test $b == lsame || test $b == xerbla || test $b == zlarft
    then
        echo Skipping $b
    else
        f2c -C++ -a -f $f && mv $b.c $b.cpp || exit 2
        # silence c++ compiler warnings about string constants
        sed -i -e 's/\("[^"]\+"\)/(char *)\1/g' -e 's/^extern.*"C"/extern "C"/' \
            -e 's/^#include.*"f2c.h"/#include "lmp_f2c.h"/' $b.cpp
    fi
done

# translate modified versions
for f in static/*.f
do \
    b=$(basename $f .f)
    f2c -C++ -a -f $f && mv $b.c $b.cpp || exit 2
    # silence c++ compiler warnings about string constants
    sed -i -e 's/\("[^"]\+"\)/(char *)\1/g' -e 's/^extern.*"C"/extern "C"/' \
        -e 's/^#include.*"f2c.h"/#include "lmp_f2c.h"/' $b.cpp
done

# copy direct C++ alternatives
for c in static/*.cpp
do \
    cp -v $c .
done
