#!/bin/sh

# convert ptx assembly output into 
# a c-style string constant written
# in portable posix shell script.
# requires: sed, rm, mv
#
# Author: Axel Kohlmeyer, Temple University
 
num_args=$#

# we write to a scratch file, since
# we know the real file name only at
# the very end.
output=geryon.tmp.$$
: > $output

# remove temporary file in case we're interrupted. 
cleanup () {
  rm -f geryon.tmp.$$
}
trap cleanup INT QUIT TERM

# loop over arguments and convert to 
# string constants. 
i=1
while [ $i -lt $num_args ]
do \
  src=$1
  krn=${src##*/}
  krn=${krn%.*}
  echo "Converting kernel $krn from $src to a c-style string"
  echo "const char * $krn = " >> $output
  sed -e 's/\\/\\\\/g'   \
      -e 's/"/\\"/g'     \
      -e 's/ *\/\/.*$//' \
      -e '/\.file/D'     \
      -e '/^[ 	]*$/D'   \
      -e 's/^\(.*\)$/"\1\\n"/' $src >> $output
  echo ';' >> $output
  shift
  i=`expr $i + 1`
done

# $1 holds now the real output file name
mv $output $1

