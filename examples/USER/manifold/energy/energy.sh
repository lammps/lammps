#!/bin/bash
#
# Runs the various tests.
#

# Set this to your lammps executable
LMP="mpirun -np 2 lmp"

NO_GNUPLOT=$( command -v gnuplot >/dev/null 2>&1  )
NSTEPS=100000000
if [ $# -ge 1 ]; then
	NSTEPS=$1
fi
for MANIFOLD in torus sphere cylinder
do
	echo "+=== Running example on "$MANIFOLD" ===+"
	DATE=$(date | awk '{ print $3$2substr($6,3) }')
	LOG=$MANIFOLD".energy."$DATE".log"
	$LMP -in energy.in -var manifold $MANIFOLD -log $MANIFOLD".energy.log" -var STEPS $NSTEPS
	echo "+======================================+"
	echo ""
done
# 2D plane for verification:
$LMP -in energy.plane.in -log energy.plane.log -var STEPS $NSTEPS


if [ $NO_GNUPLOT ]
then
	echo "No Gnuplot found, not plotting."
	exit 0
fi

E0S=$( head -n 2 thermo.sphere.dat    | tail -n 1 | awk '{ print $4 }' )
E0C=$( head -n 2 thermo.cylinder.dat  | tail -n 1 | awk '{ print $4 }' )
E0T=$( head -n 2 thermo.torus.dat     | tail -n 1 | awk '{ print $4 }' )
E0P=$( head -n 2 thermo.plane.dat     | tail -n 1 | awk '{ print $4 }' )

echo "Plotting using Gnuplot"
gnuplot -e "E0S="$E0S -e "E0C="$E0C -e "E0T="$E0T -e "E0P="$E0P plot_energies.gpl
