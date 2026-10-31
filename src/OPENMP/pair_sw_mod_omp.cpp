// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
                    and Wengen Ouyang (Wuhan U)
------------------------------------------------------------------------- */

#include "pair_sw_mod_omp.h"

#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

PairSWMODOMP::PairSWMODOMP(LAMMPS *lmp) : PairSWOMP(lmp)
{
  sw_mod_flag = 1;
  delta1 = 0.25;
  delta2 = 0.35;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSWMODOMP::settings(int narg, char **arg)
{
  // process optional keywords

  int iarg = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"maxdelcs") == 0) {
      if (iarg+3 > narg) utils::missing_cmd_args(FLERR,"pair_style sw/mod/omp", error);
      delta1 = utils::numeric(FLERR,arg[iarg+1],false,Pointers::lmp);
      delta2 = utils::numeric(FLERR,arg[iarg+2],false,Pointers::lmp);
      iarg += 3;
      if ((delta1 < 0.0) || (delta1 > 1.0) || (delta2 < 0.0) || (delta2 > 1.0) || (delta1 > delta2))
        error->all(FLERR, Error::NOPOINTER, "Illegal values for maxdelcs keyword");
    } else error->all(FLERR, Error::NOPOINTER,
                      "Unkown pair_style sw/mod/omp keyword {}", arg[iarg]);
  }
  PairSWOMP::settings(narg-iarg,arg+iarg);
}
