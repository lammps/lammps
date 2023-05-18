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
   Contributing author: Wengen Ouyang (Wuhan University)
   e-mail: w.g.ouyang at gmail dot com

   This is a full version of the potential described in
   [Feng and Ouyang et al, J. Phys. Chem. C 127, 8704-8713 (2023).]
------------------------------------------------------------------------- */

#include "pair_ilp_water_2dm.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "interlayer_taper.h"
#include "memory.h"
#include "my_page.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace InterLayer;

#define MAXLINE 1024
#define DELTA 4
#define PGDELTA 1

static const char cite_ilp_water[] = "ilp/water/2dm potential doi/10.1021/acs.jpcc.2c08464\n"
                                   "@Article{Feng2023\n"
                                   " author = {Z. Feng, Y. Yao, J. Liu, B. Wu, Z. Liu, and W. Ouyang},\n"
                                   " title = {Registry-Dependent Potential for Interfaces of Water with Graphene},\n"
                                   " journal = {J. Phys. Chem. C},\n"
                                   " volume =  127,\n"
                                   " pages =   {8704-8713}\n"
                                   " year =    2023,\n"
                                   "}\n\n";

/* ---------------------------------------------------------------------- */

PairILPWATER2DM::PairILPWATER2DM(LAMMPS *lmp) : PairILPGrapheneHBN(lmp), PairILPTMD(lmp)
{
  variant = ILP_WATER_2DM;
  single_enable = 0;

  // for TMD, each atom have six neighbors
  Nnei = 6;

  if (lmp->citeme) lmp->citeme->add(cite_ilp_water);
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairILPWATER2DM::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2) error->all(FLERR, "Illegal pair_style command");
  if (!utils::strmatch(force->pair_style, "^hybrid/overlay"))
    error->all(FLERR, "Pair style ilp/water/2dm must be used as sub-style with hybrid/overlay");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);
  if (narg == 2) tap_flag = utils::numeric(FLERR, arg[1], false, lmp);
}
