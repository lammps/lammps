#!/bin/sh

# convert ptx assembly output into
# a c-style string constant written
# in portable posix shell script.
# requires: sed, rm, mv
#
# Author: Axel Kohlmeyer, Temple University

num_args=$#

# Check command-line arguments
if [ $num_args -gt 9 ]; then
  echo "$0 can only take 9 arguments; not $num_args"
  exit 1
fi

if [ $num_args -lt 3 ]; then
  echo "Not enough arguments."
  echo "$0 name_for_string input_file1 input_file2 ... output"
  exit 1
fi

# Name is first arg, output file is last argument
string_name=$1
eval output=\${$num_args}
shift

# remove temporary file in case we're interrupted.
cleanup () {
  rm -f $output
}
trap cleanup INT QUIT TERM

# loop over arguments and convert to
# string constant.
i=2
echo "const char * $string_name = " > $output
while [ $i -lt $num_args ]
do \
  src=$1
  krn=${src##*/}
  krn=${krn%.*}
  echo "Converting $src to a c-style string"
  sed -e 's/\\/\\\\/g'   \
      -e 's/"/\\"/g'     \
      -e 's/ *\/\/.*$//' \
      -e '/\.file/D'     \
      -e '/^[ 	]*$/D'   \
      -e 's/^\(.*\)$/"\1\\n"/' $src >> $output
  shift
  i=`expr $i + 1`
done
echo ';' >> $output
