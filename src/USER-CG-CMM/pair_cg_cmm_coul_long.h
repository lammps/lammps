/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef PAIR_CG_CMM_COUL_LONG_H
#define PAIR_CG_CMM_COUL_LONG_H

#include "pair_cmm_common.h"

namespace LAMMPS_NS {

  class PairCGCMMCoulLong : public PairCMMCommon {
    public:
    PairCGCMMCoulLong(class LAMMPS *);
    ~PairCGCMMCoulLong();

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
    void *extract(char *str);

    protected:
    void allocate();
    void init_tables();
    void free_tables();
  };
}

#endif
