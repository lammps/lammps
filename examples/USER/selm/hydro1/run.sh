#!/bin/sh 

OUT_DIR="$PWD/output/run_sh/test_001"
MODEL_DIR="$PWD/Model1"

#rm -ri $OUT_DIR
rm -rf $OUT_DIR

mkdir -p $OUT_DIR
mkdir -p $OUT_DIR/vtk

cp -r $MODEL_DIR/* $OUT_DIR

echo ""
echo "Setting up symbolic links"
echo $PWD

# May need to adjust these paths and links to point to 
# the compiled lammps binaries
PATH_LAMMPS_SRC=$PWD/../../../../src/

PATH_LAMMPS_BIN=$PATH_LAMMPS_SRC
PATH_LIBSELM=$PATH_LAMMPS_SRC/../lib/selm

# Give the symbolic link pointing to the compiled binary.
#ln -s (path-to-lammps/lmp) lmp_selm_lammps
ln -s $PATH_LAMMPS_BIN/lmp_atz_selm_serial $OUT_DIR/lmp_selm_lammps

echo ""
echo "Run the simulation:"

# Give the paths pointing to the library libselm codes
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PATH_LAMMPS_BIN
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PATH_LIBSELM

echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH 

echo ""
echo "Changing to output directory to run simulations."
cd $OUT_DIR

echo ""
echo "List files:"
ls -lah 

# run the simulation
./lmp_selm_lammps -in Model.LAMMPS_script


