#!/bin/csh
date > test_results.txt
set total = `/bin/ls *.in | grep -v cpu | wc -l`
set command = "$argv[1]"

set runs = `/bin/ls *.in | grep -v cpu`
set counter = 1
foreach run ( $runs )
  rm -f /tmp/$run:r.log
  echo "Running $run on CPU ($counter of $total)"
  $command -screen none -log /tmp/$run:r.log -in $run
  @ counter = $counter + 1
end

echo "\nGPU Comparison" >> test_results.txt
echo "--------------------------------------------------" >> test_results.txt
source _run_gpu_comp.sh | sort -g -k 5 -r >>& test_results.txt
echo "\n\n\n" >> test_results.txt

echo "CPU Comparison" >> test_results.txt
echo "--------------------------------------------------" >> test_results.txt
source _run_cpu_comp.sh | sort -g -k 5 -r >>& test_results.txt
echo "\n\n\n" >> test_results.txt

echo "Time Comparison" >> test_results.txt
echo "--------------------------------------------------" >> test_results.txt
source _run_time_comp.sh | sort -g -k 7 -r >>& test_results.txt
echo "\n\n\n" >> test_results.txt

date >> test_results.txt

set error_count = `grep ERROR test_results.txt | wc -l`
if ( $error_count == "0" ) then
  set mail_sub = "Passed: LAMMPS Regression"
else
  set mail_sub = "Failed: LAMMPS Regression"
endif

mail -s "$mail_sub" brownw@ornl.gov < test_results.txt
