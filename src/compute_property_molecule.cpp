/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "compute_property_molecule.h"
#include "atom.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePropertyMolecule::
ComputePropertyMolecule(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute property/molecule command");

  if (atom->molecular == 0)
    error->all(FLERR,"Compute property/molecule requires molecular atom style");

  nvalues = narg - 3;

  pack_choice = new FnPtrPack[nvalues];

  int i;
  for (int iarg = 3; iarg < narg; iarg++) {
    i = iarg-3;

    if (strcmp(arg[iarg],"mol") == 0)
      pack_choice[i] = &ComputePropertyMolecule::pack_mol;
    else if (strcmp(arg[iarg],"count") == 0)
      pack_choice[i] = &ComputePropertyMolecule::pack_count;
    else error->all(FLERR,
                    "Invalid keyword in compute property/molecule command");
  }

  // setup molecule-based data

  nmolecules = molecules_in_group(idlo,idhi);

  vector = NULL;
  array = NULL;

  if (nvalues == 1) {
    memory->create(vector,nmolecules,"property/molecule:vector");
    vector_flag = 1;
    size_vector = nmolecules;
    extvector = 0;
  } else {
    memory->create(array,nmolecules,nvalues,"property/molecule:array");
    array_flag = 1;
    size_array_rows = nmolecules;
    size_array_cols = nvalues;
    extarray = 0;
  }

  // fill vector or array with molecule values

  if (nvalues == 1) {
    buf = vector;
    (this->*pack_choice[0])(0);
  } else {
    if (array) buf = &array[0][0];
    for (int n = 0; n < nvalues; n++)
      (this->*pack_choice[n])(n);
  }
}

/* ---------------------------------------------------------------------- */

ComputePropertyMolecule::~ComputePropertyMolecule()
{
  delete [] pack_choice;
  memory->destroy(vector);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void ComputePropertyMolecule::init()
{
  int ntmp = molecules_in_group(idlo,idhi);
  if (ntmp != nmolecules)
    error->all(FLERR,"Molecule count changed in compute property/molecule");
}

/* ---------------------------------------------------------------------- */

void ComputePropertyMolecule::compute_vector()
{
  invoked_vector = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void ComputePropertyMolecule::compute_array()
{
  invoked_array = update->ntimestep;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputePropertyMolecule::memory_usage()
{
  double bytes = (bigint) nmolecules * nvalues * sizeof(double);
  if (molmap) bytes += (idhi-idlo+1) * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   one method for every keyword compute property/molecule can output
   the atom property is packed into buf starting at n with stride nvalues
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

void ComputePropertyMolecule::pack_mol(int n)
{
  for (tagint m = idlo; m <= idhi; m++)
    if (molmap == NULL || molmap[m-idlo] >= 0) {
      buf[n] = m;
      n += nvalues;
    }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyMolecule::pack_count(int n)
{
  int i,m;
  tagint imol;

  int *count_one = new int[nmolecules];
  for (m = 0; m < nmolecules; m++) count_one[m] = 0;

  tagint *molecule = atom->molecule;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      count_one[imol]++;
    }

  int *count_all = new int[nmolecules];
  MPI_Allreduce(count_one,count_all,nmolecules,MPI_INT,MPI_SUM,world);

  for (m = 0; m < nmolecules; m++)
    if (molmap == NULL || molmap[m] >= 0) {
      buf[n] = count_all[m];
      n += nvalues;
    }

  delete [] count_one;
  delete [] count_all;
}
