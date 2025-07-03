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
#include <cstring>
using namespace LAMMPS_NS;


ComputeExtraPropertyAtom::ComputeExtraPropertyAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), extra_property_data(nullptr), data(nullptr), descriptor(nullptr)
{
  peratom_flag = 1;
  nmax = 0;
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
    if(strcmp(arg[curArg], "dim") == 0) { //Parse dim
      size_peratom_cols = atoi(arg[curArg+1]);
      dimFlag = 1;
      if (size_peratom_cols <= 0) {
        error->all(FLERR, "Compute extra_property/atom requires dim > 0.");
      }
    }
    curArg += 2; //Must put inside ifs if variable length args
  }

  //If not all essential properties were set. Error
  if (nameFlag == 0 || dimFlag == 0)
  {
    error->all(FLERR, "Compute extra_property/atom missing args (Expecting name and dim).");
  }
}

/* ---------------------------------------------------------------------- */

ComputeExtraPropertyAtom::~ComputeExtraPropertyAtom()
{
  memory->destroy(extra_property_data);

  //data is managed by PairMLIAP so we should only set it to nullptr
  data = nullptr;
  descriptor = nullptr;
}

/* ---------------------------------------------------------------------- */

void ComputeExtraPropertyAtom::init()
{
  // Check if a pair style has been defined
  if (force->pair == nullptr)
    error->all(FLERR,"Compute extra_property/atom requires a pair style be defined");

  // Check if it is safe downcast to PairMLIAP pair style
  PairMLIAP* castedPair = dynamic_cast<PairMLIAP *>(force->pair);
  if (castedPair == nullptr)
    error->all(FLERR, "Compute extra_property/atom requires pair mliap to be active");  

  //Get the pointers
  data = castedPair->data;
  descriptor = castedPair->descriptor;

  //Register the extra_property compute with the data class
  extra_property_index = data->register_extra_property(extra_property_name, size_peratom_cols);
}

/* ---------------------------------------------------------------------- */

void ComputeExtraPropertyAtom::compute_peratom()
{

  //Resize or create eatom_stdev array if needed
  if (extra_property_data == nullptr || atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->create(extra_property_data, nmax, size_peratom_cols, "compute_extra_property_atom:extra_property_data");
    array_atom = extra_property_data;
  }

  //Copy the values in data to eatom_stdev
  for (int i = 0; i < atom->nlocal; i++) {
    for (int j = 0; j < size_peratom_cols; j++) {
      extra_property_data[i][j] = data->extra_properties[extra_property_index][i][j];
    }
  }
}

/* ---------------------------------------------------------------------- */

double ComputeExtraPropertyAtom::memory_usage()
{
  double bytes = nmax * size_peratom_cols * sizeof(double); //extra_property_data

  return bytes;
}
