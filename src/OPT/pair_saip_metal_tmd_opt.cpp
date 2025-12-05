/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   This is an optimized version of ilp/tmd based on the contribution of:
     author: Wengen Ouyang (Wuhan University)
     e-mail: w.g.ouyang at gmail dot com

   Optimizations are done by:
     author1: Xiaohui Duan (National Supercomputing Center in Wuxi, China)
     e-mail: sunrise_duan at 126 dot com

     author2: Ping Gao (National Supercomputing Center in Wuxi, China)
     e-mail: qdgaoping at gmail dot com

   Optimizations are described in:
     Gao, Ping and Duan, Xiaohui, et al.:
       LMFF: Efficient and Scalable Layered Materials Force Field on Heterogeneous Many-Core Processors
     DOI: 10.1145/3458817.3476137

   Potential is described by:
     [Ouyang et al, J. Chem. Theory Comput. 17, 7237 (2021).]
*/
#include "pair_saip_metal_tmd_opt.h"
#include "atom.h"
#include "interlayer_taper.h"
#include "memory.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace InterLayer;

PairSAIPMETALTMDOpt::PairSAIPMETALTMDOpt(LAMMPS *lmp) :
    PairILPGrapheneHBN(lmp), PairILPTMD(lmp), PairSAIPMETALTMD(lmp), PairILPGrapheneHBNOpt(lmp),
    PairILPTMDOpt(lmp)
{
}

void PairSAIPMETALTMDOpt::coeff(int narg, char **args)
{
  PairILPTMDOpt::coeff(narg, args);
  for (int i = 1; i <= atom->ntypes; i++) {
    int itype = map[i];
    if (strcmp(elements[itype], "Mo") == 0 || strcmp(elements[itype], "W") == 0 ||
        strcmp(elements[itype], "S") == 0 || strcmp(elements[itype], "Se") == 0 ||
        strcmp(elements[itype], "Te") == 0) {
      special_type[i] = TMD_METAL;
    }
  }
  for (int i = 1; i <= atom->ntypes; i++) {
    int itype = map[i];
    if (strcmp(elements[itype], "Au") == 0 || strcmp(elements[itype], "Cu") == 0 ||
        strcmp(elements[itype], "Ag") == 0 || strcmp(elements[itype], "Ru") == 0 ||
        strcmp(elements[itype], "Pt") == 0 || strcmp(elements[itype], "Ni") == 0) {
      special_type[i] = SAIP_METAL;
    }
  }
}
