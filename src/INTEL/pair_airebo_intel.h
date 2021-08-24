// clang-format off
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

/* ----------------------------------------------------------------------
   Contributing author: Markus Hohnerbach (RWTH)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(airebo/intel,PairAIREBOIntel);
// clang-format on
#else

#ifndef LMP_PAIR_AIREBO_INTEL_H
#define LMP_PAIR_AIREBO_INTEL_H

#include "fix_intel.h"
#include "pair.h"
#include "pair_airebo.h"

namespace LAMMPS_NS {

template <class flt_t, class acc_t> struct PairAIREBOIntelParam;

class PairAIREBOIntel : public PairAIREBO {
 public:
  PairAIREBOIntel(class LAMMPS *);
  virtual ~PairAIREBOIntel();
  virtual void compute(int, int);
  virtual void init_style();

 protected:
  template <class flt_t, class acc_t>
  void compute(int eflag, int vflag, IntelBuffers<flt_t, acc_t> *buffers);

  template <int EVFLAG, int EFLAG, class flt_t, class acc_t>
  void eval(const int offload, const int vflag, IntelBuffers<flt_t, acc_t> *buffers,
            const int astart, const int aend);

  template <class flt_t, class acc_t> void pack_force_const(IntelBuffers<flt_t, acc_t> *buffers);

  template <class flt_t, class acc_t> PairAIREBOIntelParam<flt_t, acc_t> get_param();

  FixIntel *fix;
  int _cop;

  int *REBO_cnumneigh;
  int *REBO_num_skin;
  int *REBO_list_data;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style AIREBO requires atom IDs

This is a requirement to use the AIREBO potential.

E: Pair style AIREBO requires newton pair on

See the newton command.  This is a restriction to use the AIREBO
potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Neighbor list overflow, boost neigh_modify one

There are too many neighbors of a single atom.  Use the neigh_modify
command to increase the max number of neighbors allowed for one atom.
You may also want to boost the page size.

E: Cannot open AIREBO potential file %s

The specified AIREBO potential file cannot be opened.  Check that the
path and name are correct.

E: Cannot yet use airebo/intel with hybrid.

Pair style airebo/intel cannot currently be used as part of a hybrid
pair style (with the exception of hybrid/overlay).


*/
