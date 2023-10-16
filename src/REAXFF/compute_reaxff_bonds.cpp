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
   Contributing author: Richard Berger (LANL)
------------------------------------------------------------------------- */

#include "compute_reaxff_bonds.h"
#include "atom.h"
#include "molecule.h"
#include "update.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "neigh_list.h"

#include "pair_reaxff.h"
#include "reaxff_api.h"

using namespace LAMMPS_NS;
using namespace ReaxFF;

/* ---------------------------------------------------------------------- */

ComputeReaxFFBonds::ComputeReaxFFBonds(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  abo(nullptr), neighid(nullptr), numneigh(nullptr), reaxff(nullptr)
{
  if (atom->tag_consecutive() == 0)
    error->all(FLERR,"Atom IDs must be consecutive for compute reaxff/bonds");

  local_flag = 1;
  peratom_flag = 1;

  // initialize output

  nlocal = -1;
  prev_nbonds = -1;

  size_peratom_cols = 7;

  size_local_rows = 0;
  size_local_cols = 3;

  invoked_bonds = -1;
}


/* ---------------------------------------------------------------------- */

ComputeReaxFFBonds::~ComputeReaxFFBonds()
{
  memory->destroy(array_local);
  memory->destroy(array_atom);
  memory->destroy(abo);
  memory->destroy(neighid);
  memory->destroy(numneigh);
}

/* ---------------------------------------------------------------------- */

void ComputeReaxFFBonds::init()
{
  reaxff = dynamic_cast<PairReaxFF *>(force->pair_match("^reax..",0));
  if (reaxff == nullptr) error->all(FLERR,"Cannot use compute reaxff/bonds without "
                                          "pair_style reaxff, reaxff/kk, or reaxff/omp");
}

/* ---------------------------------------------------------------------- */

int ComputeReaxFFBonds::FindBond()
{
  int *ilist, i, ii, inum;
  int j, pj, nj;
  tagint jtag;
  double bo_tmp,bo_cut;

  inum = reaxff->list->inum;
  ilist = reaxff->list->ilist;
  bond_data *bo_ij;
  bo_cut = reaxff->api->control->bg_cut;

  tagint *tag = atom->tag;
  int numbonds = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    nj = 0;

    for (pj = Start_Index(i, reaxff->api->lists); pj < End_Index(i, reaxff->api->lists); ++pj) {
      bo_ij = &(reaxff->api->lists->select.bond_list[pj]);
      j = bo_ij->nbr;
      jtag = tag[j];
      bo_tmp = bo_ij->bo_data.BO;

      if (bo_tmp > bo_cut) {
        neighid[i][nj] = jtag;
        abo[i][nj] = bo_tmp;
        nj++;
      }
    }
    numneigh[i] = nj;
    numbonds += nj;
  }
  return numbonds;
}


/* ---------------------------------------------------------------------- */

void ComputeReaxFFBonds::compute_bonds()
{
  invoked_bonds = update->ntimestep;

  if (atom->nlocal > nlocal) {
    memory->destroy(abo);
    memory->destroy(neighid);
    memory->destroy(numneigh);
    memory->destroy(array_atom);
    nlocal = atom->nlocal;
    memory->create(abo, nlocal, MAXREAXBOND, "reaxff/bonds:abo");
    memory->create(neighid, nlocal, MAXREAXBOND, "reaxff/bonds:neighid");
    memory->create(numneigh, nlocal, "reaxff/bonds:numneigh");
    memory->create(array_atom, nlocal, 7, "reaxff/bonds:array_atom");
  }

  for (int i = 0; i < nlocal; i++) {
    numneigh[i] = 0;
    for (int j = 0; j < MAXREAXBOND; j++) {
      neighid[i][j] = 0;
      abo[i][j] = 0.0;
    }
  }

  nbonds = FindBond();
}

/* ---------------------------------------------------------------------- */

void ComputeReaxFFBonds::compute_local()
{
  invoked_local = update->ntimestep;

  if(invoked_bonds < update->ntimestep) {
    compute_bonds();
  }

  if(nbonds > prev_nbonds) {
    // grow array_local
    memory->destroy(array_local);
    memory->create(array_local, nbonds, 3, "reaxff/bonds:array_local");
    prev_nbonds = nbonds;
  }

  size_local_rows = nbonds;

  int b = 0;

  for (int i = 0; i < nlocal; ++i) {
    const int numbonds = numneigh[i];

    for (int k = 0; k < numbonds; k++) {
      auto bond = array_local[b++];
      bond[0] = i;
      bond[1] = neighid[i][k];
      bond[2] = abo[i][k];
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeReaxFFBonds::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  if(invoked_bonds < update->ntimestep) {
    compute_bonds();
  }

  for (int i = 0; i < nlocal; ++i) {
    auto ptr = array_atom[i];
    ptr[0] = atom->tag[i];
    ptr[1] = atom->type[i];
    ptr[2] = numneigh[i];
    ptr[3] = (atom->molecule == nullptr)  ? 0.0 : atom->molecule[i];
    ptr[4] = reaxff->api->workspace->total_bond_order[i];
    ptr[5] = reaxff->api->workspace->nlp[i];
    ptr[6] = atom->q[i];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeReaxFFBonds::memory_usage()
{
  double bytes = (double)(nbonds*3) * sizeof(double);
  bytes += (double)(nlocal*7) * sizeof(double);
  bytes += (double)(2*nlocal*MAXREAXBOND) * sizeof(double);
  bytes += (double)(nlocal) * sizeof(int);
  return bytes;
}
