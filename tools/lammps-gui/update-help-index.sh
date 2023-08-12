#!/bin/sh
# this updates the help index table

grep '\.\. index::' ../../doc/src/*.rst | sort  | sed -e 's/^.*src\/\([^/]\+\)\.rst:/\1.html /' -e 's/\.\. \+index:: \+//' > help_index.table
