/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   This is an optimized version of saip/metal based on the contribution of:
     author: Wengen Ouyang (Wuhan University)
     e-mail: w.g.ouyang at gmail dot com

   Optimizations are done by:
     author1: Xiaohui Duan (National Supercomputing Center in Wuxi, China)
     e-mail: sunrise_duan at 126 dot com

     author2: Ping Gao (National Supercomputing Center in Wuxi, China)
     e-mail: qdgaoping at gmail dot com

   Optimizations are described in:
     Gao, Ping and Duan, Xiaohui, et al:
       LMFF: Efficient and Scalable Layered Materials Force Field on Heterogeneous Many-Core Processors
     DOI: 10.1145/3458817.3476137

   Potential is described by:
     [Ouyang et al, J. Chem. Theory Comput. 17, 7215-7223 (2021)]
*/

#include "pair_saip_metal_opt.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "interlayer_taper.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace InterLayer;

PairSAIPMETALOpt::PairSAIPMETALOpt(LAMMPS *lmp) :
    PairILPGrapheneHBN(lmp), PairSAIPMETAL(lmp), PairILPGrapheneHBNOpt(lmp)
{
}

void PairSAIPMETALOpt::coeff(int narg, char **args)
{
  PairSAIPMETAL::coeff(narg, args);
  memory->create(special_type, atom->ntypes + 1, "PairSAIPMETALOpt:special_type");
  for (int i = 1; i <= atom->ntypes; i++) {
    int itype = map[i];
    if (strcmp(elements[itype], "C") != 0 && strcmp(elements[itype], "H") != 0 &&
        strcmp(elements[itype], "B") != 0 && strcmp(elements[itype], "N") != 0) {
      special_type[i] = true;
    } else {
      special_type[i] = false;
    }
  }
}
