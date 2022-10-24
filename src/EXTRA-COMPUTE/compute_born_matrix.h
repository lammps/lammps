/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */

/*------------------------------------------------------------------------
  Contributing Authors : Germain Clavier (TUe), Aidan Thompson (Sandia)
  --------------------------------------------------------------------------*/

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(born/matrix,ComputeBornMatrix);
// clang-format on
#else

#ifndef LMP_COMPUTE_BORN_MATRIX_H
#define LMP_COMPUTE_BORN_MATRIX_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeBornMatrix : public Compute {
 public:
  ComputeBornMatrix(class LAMMPS *, int, char **);
  ~ComputeBornMatrix() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_vector() override;
  double memory_usage() override;

 private:
  // Born matrix contributions

  void compute_pairs();                     // pair and manybody
  void compute_bonds();                     // bonds
  void compute_angles();                    // angles
  void compute_dihedrals();                 // dihedrals
  void compute_numdiff();                   // stress virial finite differences
  void displace_atoms(int, int, double);    // displace atoms
  void force_clear(int);                    // zero out force array
  void update_virial();                     // recalculate the virial
  void restore_atoms(int, int);             // restore atom positions
  void virial_addon();                      // restore atom positions
  void reallocate();                        // grow the atom arrays

  int nvalues;        // length of elastic tensor
  int numflag;        // 1 if using finite differences
  double numdelta;    // size of finite strain
  int maxatom;        // allocated size of atom arrays

  int pairflag, bondflag, angleflag;
  int dihedflag, impflag;

  double *values_local, *values_global;
  class NeighList *list;

  char *id_virial;                  // name of virial compute
  class Compute *compute_virial;    // pointer to virial compute

  static constexpr int NDIR_VIRIAL = 6;    // dimension of virial and strain vectors
  static constexpr int NXYZ_VIRIAL = 3;    // number of Cartesian coordinates
  int revalbe[NDIR_VIRIAL][NDIR_VIRIAL];
  int virialVtoV[NDIR_VIRIAL];
  double **temp_x;                   // original coords
  double **temp_f;                   // original forces
  double fixedpoint[NXYZ_VIRIAL];    // displacement field origin
  int dirlist[NDIR_VIRIAL][2];       // strain cartesian indices
};
}    // namespace LAMMPS_NS

#endif
#endif
