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

#include "compute_extra_property_atom.h"
#include "compute.h"
#include "atom.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "pair_mliap.h"
#include "pair_hybrid.h"
#include <cstring>
using namespace LAMMPS_NS;


ComputeExtraPropertyAtom::ComputeExtraPropertyAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), data(nullptr)
{
  peratom_flag = 1;
  hybridIndex = -1;
  //Check for args here
  int curArg = 3;
  int nameFlag = 0, dimFlag = 0;
  while (curArg < narg - 1) //All args come in pairs of two argName argValue
  //If this changes you need to do bounds checking in each argName
  {
    if(strcmp(arg[curArg], "name") == 0) { //Parse name
      extra_property_name = arg[curArg+1];
      nameFlag = 1; 
    }
    else if(strcmp(arg[curArg], "dim") == 0) { //Parse dim
      size_peratom_cols = atoi(arg[curArg+1]);
      dimFlag = 1;
      if (size_peratom_cols <= 0) {
        error->all(FLERR, "Compute extraProperty/atom requires dim > 0.");
      }
    }
    else if(strcmp(arg[curArg], "hybridIndex") == 0) {
      hybridIndex = atoi(arg[curArg+1]);
      if (hybridIndex <= 0) {
        error->all(FLERR, "Compute extraProperty/atom requires hybridIndex > 0.");
      }
    }
    curArg += 2; //Must put inside ifs if variable length args
  }

  //If not all essential properties were set. Error
  if (nameFlag == 0 || dimFlag == 0)
  {
    error->all(FLERR, "Compute extraProperty/atom missing args (Expecting name and dim).");
  }
}

/* ---------------------------------------------------------------------- */

ComputeExtraPropertyAtom::~ComputeExtraPropertyAtom()
{
  //data is managed by PairMLIAP so we should only set it to nullptr
  data = nullptr;
}

/* ---------------------------------------------------------------------- */

void ComputeExtraPropertyAtom::init()
{
  // Check if a pair style has been defined
  if (force->pair == nullptr)
    error->all(FLERR,"Compute extraProperty/atom requires a pair style be defined");

  // Check if it is safe downcast to PairMLIAP pair style
  
  PairMLIAP* castedPair;
  if (hybridIndex == -1) {
    castedPair = dynamic_cast<PairMLIAP *>(force->pair);
    if (castedPair == nullptr)
      error->all(FLERR, "Compute extraProperty/atom requires pair mliap to be active");
  } else {
    PairHybrid* castedPairHybrid = dynamic_cast<PairHybrid *>(force->pair);
    if (castedPairHybrid == nullptr)
      error->all(FLERR, "Compute extraPorperty/atom with hybridIndex requires pair hybrid to be active");
    if (hybridIndex > castedPairHybrid->nstyles)
      error->all(FLERR, "Compute extraProperty/atom requires hybridIndex <= number of styles in pair hybird");
    castedPair = dynamic_cast<PairMLIAP *>(castedPairHybrid->styles[hybridIndex-1]);
    if (castedPair == nullptr)
      error->all(FLERR, "Compute extraProperty/atom requires pair style at hybridIndex to be a pair mliap");
  }

  //Get the pointers
  data = castedPair->data;

  //Register the extra_property compute with the data class
  data->register_extra_property(extra_property_name, size_peratom_cols);
}

/* ---------------------------------------------------------------------- */

void ComputeExtraPropertyAtom::compute_peratom()
{
  array_atom = data->extra_properties.get_pointer(extra_property_name);
}
