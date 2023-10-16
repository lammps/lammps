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

#include "compute_reaxff_bonds_local.h"
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

ComputeReaxFFBondsLocal::ComputeReaxFFBondsLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  alocal(nullptr), abo(nullptr), neighid(nullptr), numneigh(nullptr), reaxff(nullptr)
{
  if (atom->tag_consecutive() == 0)
    error->all(FLERR,"Atom IDs must be consecutive for compute reaxff/bonds/local");

  local_flag = 1;

  nvalues = 7 + 2*MAXREAXBOND;
  prev_nvalues = 0;

  // initialize output

  nlocal = atom->nlocal;

  size_local_rows = atom->nlocal;
  size_local_cols = 7 + 2*MAXREAXBOND;

  allocate(nlocal);
}


/* ---------------------------------------------------------------------- */

ComputeReaxFFBondsLocal::~ComputeReaxFFBondsLocal()
{
  memory->destroy(alocal);
  destroy();
}

void ComputeReaxFFBondsLocal::destroy()
{
  memory->destroy(abo);
  memory->destroy(neighid);
  memory->destroy(numneigh);
}

/* ---------------------------------------------------------------------- */

void ComputeReaxFFBondsLocal::allocate(int n)
{
  memory->create(abo,n,MAXREAXBOND,"reaxff/bonds/local:abo");
  memory->create(neighid,n,MAXREAXBOND,"reaxff/bonds/local:neighid");
  memory->create(numneigh,n,"reaxff/bonds/local:numneigh");
}

/* ---------------------------------------------------------------------- */

void ComputeReaxFFBondsLocal::init()
{
  reaxff = dynamic_cast<PairReaxFF *>(force->pair_match("^reax..",0));
  if (reaxff == nullptr) error->all(FLERR,"Cannot use compute reaxff/bonds/local without "
                                "pair_style reaxff, reaxff/kk, or reaxff/omp");
}

/* ---------------------------------------------------------------------- */

int ComputeReaxFFBondsLocal::FindBond()
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
        nj ++;
      }
    }
    numneigh[i] = nj;
    if (nj > numbonds) numbonds = nj;
  }
  return numbonds;
}

/* ---------------------------------------------------------------------- */

void ComputeReaxFFBondsLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and compute bond info
  if (atom->nlocal > nlocal) {
    destroy();
    allocate(atom->nlocal);
  }

  {
    const int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      numneigh[i] = 0;
      for (int j = 0; j < MAXREAXBOND; j++) {
        neighid[i][j] = 0;
        abo[i][j] = 0.0;
      }
    }
  }

  int maxnumbonds = FindBond();
  nvalues = 7+2*maxnumbonds;

  if(atom->nlocal > nlocal || nvalues > prev_nvalues) {
    reallocate();
  }

  size_local_rows = nlocal;
  size_local_cols = nvalues;

  for (int i = 0; i < nlocal; ++i) {
    auto ptr = alocal[i];
    int numbonds = numneigh[i];
    ptr[0] = atom->tag[i];
    ptr[1] = atom->type[i];
    ptr[2] = numbonds;

    int j = 3;

    for (int k = 0; k < numbonds; k++) {
      ptr[j++] = neighid[i][k];
    }

    ptr[j++] = (atom->molecule == nullptr)  ? 0.0 : atom->molecule[i];

    for (int k = 0; k < numbonds; k++) {
      ptr[j++] = abo[i][k];
    }

    ptr[j++] = reaxff->api->workspace->total_bond_order[i];
    ptr[j++] = reaxff->api->workspace->nlp[i];
    ptr[j++] = atom->q[i];

    // clear any remaining
    for(; j < nvalues; ++j) {
      ptr[j] = 0.0;
    }
  }
}

void ComputeReaxFFBondsLocal::reallocate()
{
  nlocal = atom->nlocal;

  // grow array_local
  memory->destroy(alocal);
  memory->create(alocal,nlocal,nvalues,"reaxff/bonds/local:array_local");
  array_local = alocal;

  prev_nvalues = nvalues;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeReaxFFBondsLocal::memory_usage()
{
  double bytes = (double)nlocal*nvalues * sizeof(double);
  bytes += (double)(2*nlocal*MAXREAXBOND) * sizeof(double);
  bytes += (double)(nlocal) * sizeof(int);
  return bytes;
}
