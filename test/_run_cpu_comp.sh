#!/bin/csh
set runs = `/bin/ls *.in | grep -v cpu`
foreach run ( $runs )
  lmp_diff numbers $run:r:r.cpu.log /tmp/$run:r.log 
end
