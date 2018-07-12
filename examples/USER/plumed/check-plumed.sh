#  Needs make yes-molecule
#        make yes-rigid
#        make yes-kspace
#        make yes-user-plumed    ( obviously )


#!/bin/bash

LAMMPS=../../../src/lmp_mpi

# Run first LAMMPS calculation
$LAMMPS < in.peptide-plumed

# Check PLUMED positions
nlines=`paste plumed.xyz lammps.xyz | awk '{if( $2<$6-0.0001 || $2>$6+0.0001 || $3<$7-0.0001 || $3>$7+0.0001 || $4<$8-0.0002 || $4>$8+0.0002 ) if( $5!="Timestep:" && $1!=2004 ) print $0}' | wc -l` 
if [ "$nlines" -gt 0 ] ; then
   echo ERROR passing positions from LAMMPS to PLUMED
   paste plumed.xyz lammps.xyz | awk '{if( $2<$6-0.0001 || $2>$6+0.0001 || $3<$7-0.0001 || $3>$7+0.0001 || $4<$8-0.0002 || $4>$8+0.0002 ) if( $5!="Timestep:" && $1!=2004 ) print $0, $2-$6, $3-$7, $4-$8}'
fi

# CHECK PLUMED timestep
tstep=`grep timestep in.peptide-plumed | awk '{print $2}'`
tstep=`echo $tstep \* 0.001 \* 10 | bc -l`
nlines=`wc -l plmd_energy | awk '{print $1}'`

for ((i=3;i<$nlines;i++)); do
   told=`head -n $(($i-1)) plmd_energy | tail -n 1 | awk '{print $1}'`
   tnew=`head -n $i plmd_energy | tail -n 1 | awk '{print $1}'`
   tdiff=`echo \( $tnew - $told - $tstep \) \> 0 | bc -l`
   if [ $tdiff -gt 0 ] ; then
      echo ERROR passing timestep from LAMMPS to PLUMED
   fi
done

# Check PLUMED energy
tail -n +2 plmd_energy > plmd_energy2 
nlines=`paste lammps_energy plmd_energy2 | tail -n +2 | awk '{if( $2<$4-0.0001 || $2>$4+0.0001 ) print $0}' | wc -l` 
if [ "$nlines" -gt 0 ] ; then
   echo ERROR passing potential energy from LAMMPS to PLUMED
   paste lammps_energy plmd_energy2 | tail -n +2 | awk '{if( $2<$4-0.0001 || $2>$4+0.0001 ) print $0}'
fi
rm -f plmd_energy2

# Check PLMD mass and charge
nlines=`wc -l mq_lammps | awk '{print $1}'`
sline=`grep -n "mass q" mq_lammps | awk '{print $1}' | sed -e s/:ITEM://`
for ((i=$sline+1;i<$nlines;i++)); do
    # Mass and charge from LAMMPS
    index=`head -n $i mq_lammps | tail -n 1 | awk '{print $1}'`
    l_mass=`head -n $i mq_lammps | tail -n 1 | awk '{print $2}'`
    l_charge=`head -n $i mq_lammps | tail -n 1 | awk '{print $3}'`
    # Mass and charge from PLUMED
    p_mass=`head -n $(($index+1)) mq_plumed | tail -n 1 | awk '{print $2}'`
    p_charge=`head -n $(($index+1)) mq_plumed | tail -n 1 | awk '{print $3}'`
    # Check PLUMED mass is same as lammps mass
    mdiff=`echo \( $l_mass - $p_mass \) \> 0 | bc -l` 
    if [ "$mdiff" -gt 0 ] ; then
         echo ERROR passing masses from LAMMPS to PLUMED
    fi 
    # Check PLUMED charge is same as lammps charge
    qdiff=`echo \( $l_charge - $p_charge \) \> 0 | bc -l` 
    if [ "$qdiff" -gt 0 ] ; then
         echo ERROR passing charges from LAMMPS to PLUMED
    fi 
done

# Run calculations to test adding restraint on bond
$LAMMPS < in.peptide-plumed-plumed-restraint
$LAMMPS < in.peptide-plumed-lammps-restraint
# Now compare value of distance when lammps and plumed restraint the distance
nlines=`paste lammps_restraint plumed_restraint | tail -n +2 | awk '{if( $2<$4-0.0001 || $2>$4+0.0001 ) print $0}' | wc -l`
if [ "$nlines" -gt 0 ] ; then
      echo ERROR passing forces from PLUMED back to LAMMPS
fi

# Now run calculations to test virial
$LAMMPS < in.peptide-plumed-npt
$LAMMPS < in.peptide-plumed-npt2	

nlines=`paste plmd_volume_with_restraint plmd_volume_without_restraint | tail -n +2 | awk '{if( $2<$4-0.0001 || $2>$4+0.0001 ) print $0}' | wc -l`
if [ "$nlines" -gt 0 ] ; then
      echo ERROR passing virial from PLUMED back to LAMMPS
fi

# Nothing from here works

# Now run calculations to check forces on energy
$LAMMPS < in.peptide-plumed-engforce-ref
$LAMMPS < in.peptide-plumed-eng-force-plumed
