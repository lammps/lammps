#!/bin/tcsh

set execs="asphere_vis"
set execdir="../bin"

if ( -e ../doc ) then
  echo Manpage directory exists...
else
  echo Creating directory 'manpages'
  mkdir ../doc
endif

cd ../doc

foreach exec ($execs)
  $execdir/$exec -h > $exec.manpage
  eqn $exec.manpage > $exec.1
  man -t -p eqn ./$exec.manpage > $exec.ps
  ps2pdf $exec.ps $exec.pdf
end

cd ../src


