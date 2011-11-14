#!/bin/csh
set total = `/bin/ls *.cpu.in | wc -l`
set command = "$argv[1]"

set runs = `echo *.cpu.in`
set counter = 1
foreach run ( $runs )
  echo "Running $run on CPU ($counter of $total)"
  $command -screen none -log $run:r.log -in $run
  @ counter = $counter + 1
end
