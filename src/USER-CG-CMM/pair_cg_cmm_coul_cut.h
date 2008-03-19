/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef PAIR_CG_CMM_COUL_CUT_H
#define PAIR_CG_CMM_COUL_CUT_H

#include "pair_cmm_common.h"

namespace LAMMPS_NS {

    class PairCGCMMCoulCut : public PairCMMCommon {
    public:
    PairCGCMMCoulCut(class LAMMPS *);
    ~PairCGCMMCoulCut();

    void compute(int, int);
    void compute_inner();
    void compute_middle();
    void compute_outer(int, int);

    void init_style();
    double init_one(int, int);

    void write_restart(FILE *);
    void read_restart(FILE *);

    double memory_usage();

    double single(int, int, int, int, double, double, double, double &);

    protected:
    void allocate();
  };
}

#endif
