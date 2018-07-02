#!/bin/bash

# Data from first set of checks
rm bck.* plmd_energy lammps_energy mq_plumed mq_lammps lammps.xyz plumed.xyz p.log
# Data from checks on restraints
rm bck.* p.log lammps_restraint plumed_restraint
# Data from checks on virial
rm bck.* lammps_energy lammps.xyz log.lammps plmd_volume p.log plmd_volume_without_restraint plmd_volume_with_restraint
