now="$(date)" 
printf "Started at %s\n" "$now"

# Always generate a clean folder logs that helps us see just how far the tests ran
rm -r logs
mkdir logs

cd cgheptaneFileSet
# Test error conditions/messages
/mnt/c/openSource/eagunn/lammps/build/lmp -in cgheptane.unrecognizedCGBasisType.lmp > ../logs/cgheptane.unrecognizedCGBasisType.log
/mnt/c/openSource/eagunn/lammps/build/lmp -in cgheptane.badAnalyticCoefficients.lmp > ../logs/cgheptane.badAnalyticCoefficients.log
/mnt/c/openSource/eagunn/lammps/build/lmp -in cgheptane.badSplineDataFileName.lmp > ../logs/cgheptane.badSplineDataFileName.log
/mnt/c/openSource/eagunn/lammps/build/lmp -in cgheptane.badFvData.lmp > ../logs/cgheptane.badFvData.log
/mnt/c/openSource/eagunn/lammps/build/lmp -in cgheptane.unrecognizedKeyword.lmp > ../logs/cgheptane.unrecognizedKeyword.log

# Test valid input

# cubic splines, 200 steps
/mnt/c/openSource/eagunn/lammps/build/lmp -in cgheptane.cubic200.lmp > ../logs/cgheptane.cubic200.log
# cubic splines, 1000 steps -- original sample file received from Michael
/mnt/c/openSource/eagunn/lammps/build/lmp -in cgheptane.lmp > ../logs/cgheptane.log

# linear splines
/mnt/c/openSource/eagunn/lammps/build/lmp -in cgheptane.linear200.lmp > ../logs/cgheptane.linear200.log

now="$(date)" 
printf "cgheptane tests finished at %s\n" "$now"

# We have the example file set, uses analytic mode
cd ../exampleMethaneFileSet
/mnt/c/openSource/eagunn/lammps/build/lmp -in in.methanol > ../logs/methanol.log

#return home
cd ..

now="$(date)" 
printf "All tests finished at %s\n" "$now"
 