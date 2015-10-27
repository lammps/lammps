#!/bin/sh
# do some cleaning up to fix permissions and remove cruft

# delete backup and patch rejection files
find ./ -name \*.orig -print -or -name \*~ -print -or -name \*.rej -print \
  | xargs rm -vf

# remove executable permissions from sources
find ./ -name \*.cpp -print -or -name \*.c -print \
    -or -name \*.h -print -or -name Makefile\* -print \
  | xargs chmod -x $f

# remove training whitespace from C/C++ source files
find ./ -name \*.cpp -print -or -name \*.c -print  -or -name \*.h -print \
  | xargs sed -i -e 's,[	 ]\+$,,'

# change #include "..." to #include <...> for system headers
find ./ -name \*.cpp -print -or -name \*.c -print  -or -name \*.h -print \
  | xargs sed -i -e 's,^#include \+"\(Python\|assert\|ctype\|direct\|dirent\|errno\|float\|inttypes\|limits\|math\|mpi\|omp\|rpc/.\*\|stdint\|stdio\|stdlib\|string\|sys/\*\|time\|unistd\)\.h",#include <\1.h>,'
