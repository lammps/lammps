/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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
  ~FixWidom();
  int setmask();
  void init();
  void pre_exchange();

  void attempt_atomic_insertion();
  void attempt_molecule_insertion();

  void attempt_atomic_insertion_full();
  void attempt_molecule_insertion_full();
  double energy(int, int, tagint, double *);
  double molecule_energy(tagint);
  double energy_full();
  void update_gas_atoms_list();
  double compute_vector(int);
  double memory_usage();
  void write_restart(FILE *);
  void restart(char *);
  void grow_molecule_arrays(int);

 private:
  int molecule_group, molecule_group_bit;
  int molecule_group_inversebit;
  int exclusion_group, exclusion_group_bit;
  int nwidom_type, nevery, seed;
  int ninsertions;
  int ngas;            // # of gas atoms on all procs
  int ngas_local;      // # of gas atoms on this proc
  int exchmode;        // exchange ATOM or MOLECULE
  int movemode;        // move ATOM or MOLECULE
  int regionflag;      // 0 = anywhere in box, 1 = specific region
  int iregion;         // widom region
  char *idregion;      // widom region id
  bool charge_flag;    // true if user specified atomic charge
  bool full_flag;      // true if doing full system energy calculations

  int natoms_per_molecule;    // number of atoms in each inserted molecule
  int nmaxmolatoms;           // number of atoms allocated for molecule arrays

  double ave_widom_chemical_potential;

  int widom_nmax;
  int max_region_attempts;
  double gas_mass;
  double insertion_temperature;
  double beta, sigma, volume;
  double charge;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double region_xlo, region_xhi, region_ylo, region_yhi, region_zlo, region_zhi;
  double region_volume;
  double energy_stored;    // full energy of old/current configuration
  double *sublo, *subhi;
  int *local_gas_list;
  double **cutsq;
  double **molcoords;
  double *molq;
  imageint *molimage;
  imageint imagezero;

  double energy_intra;

  class Pair *pair;

  class RanPark *random_equal;

  class Atom *model_atom;

  class Molecule **onemols;
  int imol, nmol;
  int triclinic;    // 0 = orthog box, 1 = triclinic

  class Compute *c_pe;

  void options(int, char **);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix Widom does not (yet) work with atom_style template

Self-explanatory.

E: Fix Widom region does not support a bounding box

Not all regions represent bounded volumes.  You cannot use
such a region with the fix Widom command.

E: Fix Widom region cannot be dynamic

Only static regions can be used with fix Widom.

E: Fix Widom region extends outside simulation box

Self-explanatory.

E: Fix Widom molecule must have coordinates

The defined molecule does not specify coordinates.

E: Fix Widom molecule must have atom types

The defined molecule does not specify atom types.

E: Atom type must be zero in fix Widom mol command

Self-explanatory.

E: Fix Widom molecule has charges, but atom style does not

Self-explanatory.

E: Fix Widom molecule template ID must be same as atom_style template ID

When using atom_style template, you cannot insert molecules that are
not in that template.

E: Fix Widom atom has charge, but atom style does not

Self-explanatory.

UNDOCUMENTED

E: Molecule template ID for fix Widom does not exist

Self-explanatory.

W: Molecule template for fix Widom has multiple molecules

The fix Widom command will only create molecules of a single type,
i.e. the first molecule in the template.

E: Region ID for fix Widom does not exist

Self-explanatory.

W: Fix Widom using full_energy option

Fix Widom has automatically turned on the full_energy option since it
is required for systems like the one specified by the user. User input
included one or more of the following: kspace, a hybrid
pair style, an eam pair style, tail correction,
or no "single" function for the pair style.

E: Invalid atom type in fix Widom command

The atom type specified in the Widom command does not exist.

E: Fix Widom cannot exchange individual atoms belonging to a molecule

This is an error since you should not delete only one atom of a
molecule.  The user has specified atomic (non-molecular) gas
exchanges, but an atom belonging to a molecule could be deleted.

E: All mol IDs should be set for fix Widom group atoms

The molecule flag is on, yet not all molecule ids in the fix group
have been set to non-zero positive values by the user. This is an
error since all atoms in the fix Widom group are eligible for deletion,
rotation, and translation and therefore must have valid molecule ids.

E: Fix Widom molecule command requires that atoms have molecule attributes

Should not choose the Widom molecule feature if no molecules are being
simulated. The general molecule flag is off, but Widom's molecule flag
is on.

E: Cannot use fix Widom in a 2d simulation

Fix Widom is set up to run in 3d only. No 2d simulations with fix Widom
are allowed.

E: Could not find fix Widom exclusion group ID

Self-explanatory.

E: Illegal fix Widom gas mass <= 0

The computed mass of the designated gas molecule or atom type was less
than or equal to zero.

E: Cannot do Widom on atoms in atom_modify first group

This is a restriction due to the way atoms are organized in a list to
enable the atom_modify first command.

W: Fix Widom is being applied to the default group all

This is allowed, but it will result in Monte Carlo moves
being performed on all the atoms in the system, which is
often not what is intended.

E: Could not find specified fix Widom group ID

Self-explanatory.

E: fix Widom does currently not support full_energy option with molecules on more than 1 MPI process.

UNDOCUMENTED

E: Fix Widom put atom outside box

This should not normally happen.  Contact the developers.

E: Fix Widom ran out of available molecule IDs

See the setting for tagint in the src/lmptype.h file.

E: Fix Widom ran out of available atom IDs

See the setting for tagint in the src/lmptype.h file.

E: Too many total atoms

See the setting for bigint in the src/lmptype.h file.

*/
