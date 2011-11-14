#!/bin/csh
set runs = `/bin/ls *.in | grep -v cpu`
foreach run ( $runs )
  lmp_diff numbers $run:r.log /tmp/$run:r.log
end
