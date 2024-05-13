#!/bin/sh
# this updates the help index table

mv help_index.table help_index.oldtable
grep '\.\. index::' ../../doc/src/*.rst | sort  | sed -e 's/^.*src\/\([^/]\+\)\.rst:/\1.html /' -e 's/\.\. \+index:: \+//' > help_index.table
cmp help_index.table help_index.oldtable > /dev/null || touch lammpsgui.qrc
