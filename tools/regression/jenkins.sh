# Sample Jenkins run script

### Build the voro++ library

cd lammps_testing/voro++
make
make install

### Build the KIM library

cd ../kim-api
find . -name Makefile.KIM_Config -type f -print0 | xargs -0 sed -i 's^$(HOME)^$(WORKSPACE)^g'
make add-ex_model_Ne_P_fastLJ || :
make all

cd ../../

### Build the LAMMPS libraries

cd lib

# Wipe out /lib directory
svn status | grep ^\? | cut -c9- | xargs -d \\n rm -r || :

cd awpmd
rm libawpmd.a || :
make -j8 -f Makefile.mpicc
cd ../colvars
rm libcolvars.a || :
make -j8 -f Makefile.g++
cd ../meam
make -f Makefile.gfortran clean
make -f Makefile.gfortran
cd ../poems
rm libpoems.a || :
make -j8 -f Makefile.g++
cd ../qmmm
rm libqmmm.a || :
make -j8 -f Makefile.gfortran
cd ../reax
rm libreax.a || :
make -j8 -f Makefile.gfortran

cd ../../

### Build LAMMPS

cd src

# Wipe out /src directory
svn revert * --recursive
svn status | grep ^\? | cut -c9- | xargs -d \\n rm -r || :

cp ../lammps_testing/MAKE/Makefile.voronoi ../lib/voronoi/Makefile.lammps
#cp ../lammps_testing/src/* . || :

# Set VORONOI library paths
find ../lib/voronoi/ -name Makefile.lammps -path '*/.svn' -prune -o -type f -print0 | xargs -0 sed -i 's^voronoi_SYSINC = ^voronoi_SYSINC = -I${WORKSPACE}/lammps_testing/voronoi/include/voro++^g'
find ../lib/voronoi/ -name Makefile.lammps -path '*/.svn' -prune -o -type f -print0 | xargs -0 sed -i 's^voronoi_SYSLIB = ^voronoi_SYSLIB = -lvoro++^g'
find ../lib/voronoi/ -name Makefile.lammps -path '*/.svn' -prune -o -type f -print0 | xargs -0 sed -i 's^voronoi_SYSPATH = ^voronoi_SYSPATH = -L${WORKSPACE}/lammps_testing/voronoi/lib^g'

# Set KIM library paths
find ../lib/kim/ -name Makefile.lammps -path '*/.svn' -prune -o -type f -print0 | xargs -0 sed -i 's^$(shell kim-api-build-config --includes)^-I${WORKSPACE}/lammps_testing/kim-api/src^g'
find ../lib/kim/ -name Makefile.lammps -path '*/.svn' -prune -o -type f -print0 | xargs -0 sed -i 's^$(shell kim-api-build-config --ldlibs)^-lkim-api-v1^g'
find ../lib/kim/ -name Makefile.lammps -path '*/.svn' -prune -o -type f -print0 | xargs -0 sed -i 's^$(shell kim-api-build-config --ldflags)^-L${WORKSPACE}/lammps_testing/kim-api/src^g'

# Install packages
make yes-asphere
make yes-body
make yes-class2
make yes-colloid
make yes-compress
make yes-coreshell
make yes-dipole
make yes-fld
make yes-granular
#make yes-kim
make yes-kspace
make yes-manybody
make yes-mc
make yes-meam
make yes-misc
make yes-molecule
make yes-mpiio
make yes-opt
make yes-peri
make yes-poems
make yes-python
make yes-qeq
make yes-reax
make yes-replica
make yes-rigid
make yes-shock
make yes-snap
make yes-srd
make yes-voronoi
make yes-xtc
make yes-user-awpmd
make yes-user-cg-cmm
make yes-user-colvars
make yes-user-diffraction
make yes-user-drude
make yes-user-eff
make yes-user-fep
make yes-user-lb
make yes-user-misc
make yes-user-molfile
make yes-user-phonon
make yes-user-qmmm
make yes-user-qtb
make yes-user-reaxc
make yes-user-sph
make yes-user-tally
make yes-user-dpd

make -j 8 mpi

cd ${WORKSPACE}

### Run the tests

cd lammps_testing

# Rebless tests as necessary
#rm ./examples/comb/log.*.linux.*comb.Cu2O.elastic || :
#rm ./examples/USER/eff/fixed-core/C2H6/log.*.linux.*bohr || :
#find ./examples/min -path '*/.svn' -prune -o -name "log.*.linux.*" -print0 | xargs -0 rm -rf
#rm ./examples/rigid/log.*.linux.*rigid.tnr || :
#rm ./examples/USER/tally/log.08Sep15.linux.8* || :
#rm ./examples/min/log.*.linux.*.min || :
#find ./examples/reax -path '*/.svn' -prune -o -name "*linux.8*" -print0 | xargs -0 rm -rf

# These tests are currently broken
rm -rf ./examples/ELASTIC/ || :

cd regress

# Rebless tests as necessary
#python regression.py 8 "mpiexec -np 8 ${WORKSPACE}/src/lmp_face -v CORES 8" ${WORKSPACE}/lammps_testing/examples -only USER/tally -auto-rebless True -min-same-rows 0 2>&1 |tee test.out
#python regression.py 8 "mpiexec -np 8 ${WORKSPACE}/src/lmp_face -v CORES 8" ${WORKSPACE}/lammps_testing/examples -only ASPHERE/box ASPHERE/star rigid deposit -auto-rebless True -min-same-rows 0 2>&1 |tee test.out
#python regression.py 8 "mpiexec -np 8 ${WORKSPACE}/src/lmp_face -v CORES 8" ${WORKSPACE}/lammps_testing/examples -only KAPPA peptide pour VISCOSITY USER bench -auto-rebless True -min-same-rows 0 2>&1 |tee test.out

# Run example tests

python regression.py 8 "mpiexec -np 8 ${WORKSPACE}/src/lmp_face -v CORES 8" ${WORKSPACE}/lammps_testing/examples -exclude kim 2>&1 |tee test.out

grep "*** no failures ***" test.out
error_code_examples=$?

exit ${error_code_examples}