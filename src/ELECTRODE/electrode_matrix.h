/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meißner (TUHH)
------------------------------------------------------------------------- */

#include <map>

#include "electrode_kspace.h"
#include "pointers.h"

namespace LAMMPS_NS {

class ElectrodeMatrix : protected Pointers {
 public:
  ElectrodeMatrix(class LAMMPS *, int, double);
  ~ElectrodeMatrix() {}
  void setup(std::map<tagint, int>, class Pair *, class NeighList *);
  void compute_array(double **);
  int igroup;

 private:
  int groupbit;
  bigint ngroup;
  double **cutsq;
  double g_ewald, eta;
  std::map<tagint, int> tag_to_iele;
  bool assigned;
  std::vector<bigint> mpos;
  class Pair *pair;
  class NeighList *list;
  class ElectrodeKSpace *electrode_kspace;

  void update_mpos();
  void pair_contribution(double **);
  void self_contribution(double **);
  double calc_erfc(double);
};

}    // namespace LAMMPS_NS

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute group/group group ID does not exist

Self-explanatory.

E: Compute group/group molecule requires molecule IDs

UNDOCUMENTED

E: No pair style defined for compute group/group

Cannot calculate group interactions without a pair style defined.

E: Pair style does not support compute group/group

The pair_style does not have a single() function, so it cannot be
invoked by the compute group/group command.

E: No KSpace style defined for compute group/group

Self-explanatory.

E: KSpace style does not support compute group/group

Self-explanatory.

W: Both groups in compute group/group have a net charge; the KSpace boundary
correction to energy will be non-zero

Self-explanatory.

*/
