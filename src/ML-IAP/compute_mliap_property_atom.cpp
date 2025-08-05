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

#include "compute_mliap_property_atom.h"
#include "compute.h"
#include "atom.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "pair_mliap.h"
#include "pair_hybrid.h"
#include <cstring>
using namespace LAMMPS_NS;


ComputeMLIAPPropertyAtom::ComputeMLIAPPropertyAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), data(nullptr)
{
  peratom_flag = 1;
  hybrid_index = -1;
  //Check for args here
  int curArg = 3;
  int nameFlag = 0, dimFlag = 0;
  while (curArg < narg - 1) //All args come in pairs of two argName argValue
  //If this changes you need to do bounds checking in each argName
  {
    if(strcmp(arg[curArg], "name") == 0) { //Parse name
      property_name = arg[curArg+1];
      nameFlag = 1; 
    }
    else if(strcmp(arg[curArg], "dim") == 0) { //Parse dim
      size_peratom_cols = atoi(arg[curArg+1]);
      dimFlag = 1;
      if (size_peratom_cols <= 0) {
        error->all(FLERR, "Compute mliap/property/atom requires dim > 0.");
      }
    }
    else if(strcmp(arg[curArg], "hybrid_index") == 0) {
      hybrid_index = atoi(arg[curArg+1]);
      if (hybrid_index <= 0) {
        error->all(FLERR, "Compute mliap/property/atom requires hybrid_index > 0.");
      }
    }
    curArg += 2; //Must put inside ifs if variable length args
  }

  //If not all essential properties were set. Error
  if (nameFlag == 0 || dimFlag == 0)
  {
    error->all(FLERR, "Compute mliap/property/atom missing args (Expecting name and dim).");
  }
}

/* ---------------------------------------------------------------------- */

ComputeMLIAPPropertyAtom::~ComputeMLIAPPropertyAtom()
{
  //data is managed by PairMLIAP so we should only set it to nullptr
  data = nullptr;
}

/* ---------------------------------------------------------------------- */

void ComputeMLIAPPropertyAtom::init()
{
  // Check if a pair style has been defined
  if (force->pair == nullptr)
    error->all(FLERR,"Compute mliap/property/atom requires a pair style be defined");

  // Check if it is safe downcast to PairMLIAP pair style
  
  PairMLIAP* castedPair;
  if (hybrid_index == -1) {
    castedPair = dynamic_cast<PairMLIAP *>(force->pair);
    if (castedPair == nullptr)
      error->all(FLERR, "Compute mliap/property/atom requires pair mliap to be active");
  } else {
    castedPair = dynamic_cast<PairMLIAP *>(force->pair_match("mliap", 0, hybrid_index));
    if (castedPair == nullptr)
      error->all(FLERR, "Compute mliap/property/atom requires a mliap pair style at hybrid_index to be a pair mliap");
  }

  //Get the pointers
  data = castedPair->data;

  //Register the extra_property compute with the data class
  data->register_extra_property(property_name, size_peratom_cols);
}

/* ---------------------------------------------------------------------- */

void ComputeMLIAPPropertyAtom::compute_peratom()
{
  array_atom = data->extra_properties.get_pointer(property_name);
}
