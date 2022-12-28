#!/bin/bash

has_f2c=$(type f2c > /dev/null 2>&1 && echo 1 || echo 0)
if test ${has_f2c} -eq 0
then
    echo "Must have f2c installed to run this script"
    exit 1
fi

# cleanup
rm -f *.c *.cpp *.P *~ *.orig *.bak *.rej

# translate original files directly
for f in fortran/*.f
do \
    b=$(basename $f .f)
    # skip files for which we have replacements
    if test $b == dgetrf2 || test $b == disnan || test $b == dlaisnan ||    \
            test $b == dlamch || test $b == dlarft || test $b == dpotrf2 || \
            test $b == lsame || test $b == xerbla || test $b == zlarft
    then
        echo Skipping $b
    else
        # convert to C++ with f2c, make local variables dynamic,
        # strip comments, and reindent with clang-format.
        f2c -C++ -a -f < $f \
            | g++ -fpreprocessed -dD -P -E - \
            | clang-format -style=file:static/.clang-format > $b.cpp || exit 2
        # silence c++ compiler warnings about string constants, use custom f2c header
        sed -i -e 's/\("[^"]\+"\)/(char *)\1/g' -e 's/^extern.*"C"/extern "C"/' \
            -e 's/^#include.*"f2c.h"/#include "lmp_f2c.h"/' $b.cpp
        # replace libf2c functions with local versions under different names
        sed -i -e 's/s_\(cat\|cmp\|copy\)(/s_lmp_\1(/g' \
             -e 's/d_\(sign\|cnjg\|imag\|lg10\)(/d_lmp_\1(/g' \
             -e 's/z_\(abs\|div\)(/z_lmp_\1(/g'         \
             -e 's/i_\(len\|nint\|dnnt\)(/i_lmp_\1(/g'  \
             -e 's/pow_\(dd\|di\|ii\)(/pow_lmp_\1(/g' $b.cpp
    fi
done

# translate modified versions
for f in static/*.f
do \
    b=$(basename $f .f)
    # convert to C++ with f2c, make local variables dynamic,
    # strip comments, and reindent with clang-format.
    f2c -C++ -a -f < $f \
        | g++ -fpreprocessed -dD -P -E - \
        | clang-format -style=file:static/.clang-format > $b.cpp || exit 2
    # silence c++ compiler warnings about string constants, use custom f2c header
    sed -i -e 's/\("[^"]\+"\)/(char *)\1/g' -e 's/^extern.*"C"/extern "C"/' \
        -e 's/^#include.*"f2c.h"/#include "lmp_f2c.h"/' $b.cpp
    # replace libf2c functions with local versions under different names
    sed -i -e 's/s_\(cat\|cmp\|copy\)(/s_lmp_\1(/g' \
        -e 's/d_\(sign\|cnjg\|imag\|lg10\)(/d_lmp_\1(/g'  \
        -e 's/z_\(abs\|div\)(/z_lmp_\1(/g'          \
        -e 's/i_\(len\|nint\|dnnt\)(/i_lmp_\1(/g'   \
        -e 's/pow_\(dd\|di\|ii\)(/pow_lmp_\1(/g' $b.cpp
done

# copy direct C++ alternatives
for c in static/*.cpp
do \
    cp -v $c .
done
