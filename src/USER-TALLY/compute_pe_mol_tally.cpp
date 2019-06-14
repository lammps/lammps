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

#include <cstring>
#include "compute_pe_mol_tally.h"
#include "atom.h"
#include "group.h"
#include "pair.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "comm.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePEMolTally::ComputePEMolTally(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute pe/mol/tally command");

  igroup2 = group->find(arg[3]);
  if (igroup2 == -1)
    error->all(FLERR,"Could not find compute pe/mol/tally second group ID");
  groupbit2 = group->bitmask[igroup2];

  vector_flag = 1;
  size_vector = 4;
  timeflag = 1;

  extvector = 1;
  peflag = 1;                   // we need Pair::ev_tally() to be run

  did_setup = invoked_vector = -1;
  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputePEMolTally::~ComputePEMolTally()
{
  if (force && force->pair) force->pair->del_tally_callback(this);
  delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ComputePEMolTally::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Trying to use compute pe/mol/tally without pair style");
  else
    force->pair->add_tally_callback(this);

  if (atom->molecule_flag == 0)
    error->all(FLERR,"Compute pe/mol/tally requires molecule IDs");

  if (comm->me == 0) {
    if (force->pair->single_enable == 0 || force->pair->manybody_flag)
      error->warning(FLERR,"Compute pe/mol/tally used with incompatible pair style");

    if (force->bond || force->angle || force->dihedral
                    || force->improper || force->kspace)
      error->warning(FLERR,"Compute pe/mol/tally only called from pair style");
  }
  did_setup = -1;
}

/* ---------------------------------------------------------------------- */

void ComputePEMolTally::pair_setup_callback(int, int)
{
  // run setup only once per time step.
  // we may be called from multiple pair styles

  if (did_setup == update->ntimestep) return;

  etotal[0] = etotal[1] = etotal[2] = etotal[3] = 0.0;
  did_setup = update->ntimestep;
}

/* ---------------------------------------------------------------------- */
void ComputePEMolTally::pair_tally_callback(int i, int j, int nlocal, int newton,
                                         double evdwl, double ecoul, double,
                                         double, double, double)
{
  const int * const mask = atom->mask;
  const tagint * const molid = atom->molecule;

  if ( ((mask[i] & groupbit) && (mask[j] & groupbit2))
     || ((mask[i] & groupbit2) && (mask[j] & groupbit)) ){

    evdwl *= 0.5; ecoul *= 0.5;
    if (newton || i < nlocal) {
      if (molid[i] == molid[j]) {
        etotal[0] += evdwl; etotal[1] += ecoul;
      } else {
        etotal[2] += evdwl; etotal[3] += ecoul;
      }
    }
    if (newton || j < nlocal) {
      if (molid[i] == molid[j]) {
        etotal[0] += evdwl; etotal[1] += ecoul;
      } else {
        etotal[2] += evdwl; etotal[3] += ecoul;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePEMolTally::compute_vector()
{
  invoked_vector = update->ntimestep;
  if ((did_setup != invoked_vector) || (update->eflag_global != invoked_vector))
    error->all(FLERR,"Energy was not tallied on needed timestep");

  // sum accumulated energies across procs

  MPI_Allreduce(etotal,vector,size_vector,MPI_DOUBLE,MPI_SUM,world);
}

