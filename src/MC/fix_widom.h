/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(widom,FixWidom);
// clang-format on
#else

#ifndef LMP_FIX_WIDOM_H
#define LMP_FIX_WIDOM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWidom : public Fix {
 public:
  FixWidom(class LAMMPS *, int, char **);
  ~FixWidom() override;
  int setmask() override;
  void init() override;
  void pre_exchange() override;

  void attempt_atomic_insertion();
  void attempt_molecule_insertion();

  void attempt_atomic_insertion_full();
  void attempt_molecule_insertion_full();
  double energy(int, int, tagint, double *);
  double molecule_energy(tagint);
  double energy_full();
  double compute_vector(int) override;
  double memory_usage() override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  void grow_molecule_arrays(int);

 private:
  int molecule_group, molecule_group_bit;
  int molecule_group_inversebit;
  int exclusion_group, exclusion_group_bit;
  int nwidom_type, nevery, seed;
  int ninsertions;
  int exchmode;            // exchange ATOM or MOLECULE
  class Region *region;    // widom region
  char *idregion;          // widom region id
  bool charge_flag;        // true if user specified atomic charge
  bool full_flag;          // true if doing full system energy calculations

  int natoms_per_molecule;    // number of atoms in each inserted molecule
  int nmaxmolatoms;           // number of atoms allocated for molecule arrays

  double ave_widom_chemical_potential;

  int widom_nmax;
  int max_region_attempts;
  double gas_mass;
  double insertion_temperature;
  double beta, volume;
  double charge;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double region_xlo, region_xhi, region_ylo, region_yhi, region_zlo, region_zhi;
  double region_volume;
  double energy_stored;    // full energy of old/current configuration
  double *sublo, *subhi;
  double **cutsq;
  double **molcoords;
  double *molq;
  imageint *molimage;
  imageint imagezero;

  double energy_intra;

  class Pair *pair;

  class RanPark *random_equal;

  class Atom *model_atom;

  class Molecule *onemol;
  int triclinic;    // 0 = orthog box, 1 = triclinic

  class Compute *c_pe;

  void options(int, char **);
};

}    // namespace LAMMPS_NS

#endif
#endif
