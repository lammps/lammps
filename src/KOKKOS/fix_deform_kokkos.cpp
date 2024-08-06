// clang-format off
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

/* ----------------------------------------------------------------------
   Contributing author:   Pieter in 't Veld (SNL)
   Refactoring (2024/08): Mitch Murphy (alphataubio@gmail.com)
------------------------------------------------------------------------- */

#include "fix_deform_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain_kokkos.h"
#include "force.h"
#include "input.h"
#include "irregular.h"
#include "kspace.h"
#include "math_const.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{NONE=0,FINAL,DELTA,SCALE,VEL,ERATE,TRATE,VOLUME,WIGGLE,VARIABLE};
enum{ONE_FROM_ONE,ONE_FROM_TWO,TWO_FROM_ONE};

/* ---------------------------------------------------------------------- */

FixDeformKokkos::FixDeformKokkos(LAMMPS *lmp, int narg, char **arg) : FixDeform(lmp, narg, arg)
{
  kokkosable = 1;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ----------------------------------------------------------------------
  box flipped on previous step
  reset box tilts for flipped config and create new box in domain
  image_flip() adjusts image flags due to box shape change induced by flip
  remap() puts atoms outside the new box back into the new box
  perform irregular on atoms in lamda coords to migrate atoms to new procs
  important that image_flip comes before remap, since remap may change
    image flags to new values, making eqs in doc of Domain:image_flip incorrect
------------------------------------------------------------------------- */

void FixDeformKokkos::pre_exchange()
{
  atomKK->sync(Host,ALL_MASK);
  FixDeform::pre_exchange();
  atomKK->modified(Host,ALL_MASK);
}

/* ---------------------------------------------------------------------- */

void FixDeformKokkos::end_of_step()
{
  atomKK->sync(Host,ALL_MASK);
  FixDeform::end_of_step();
  atomKK->modified(Host,ALL_MASK);
}
