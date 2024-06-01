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

#include "compute_reaxff_atom.h"

#include "atom.h"
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

ComputeReaxFFAtom::ComputeReaxFFAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), neighid(nullptr), abo(nullptr), bondcount(nullptr), reaxff(nullptr)
{
  if (atom->tag_consecutive() == 0)
    error->all(FLERR, "Atom IDs must be consecutive for compute reaxff/atom");

  peratom_flag = 1;

  // initialize output

  nmax = -1;
  nbonds = 0;
  prev_nbonds = -1;

  size_peratom_cols = 3;

  size_local_rows = 0;
  size_local_cols = 3;

  invoked_bonds = -1;

  store_bonds = false;
  nsub = 0;

  int iarg = 3;
  while (iarg<narg) {
    if (strcmp(arg[iarg], "pair") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "compute reaxff/atom pair", error);
      ++iarg;

      if (isdigit(arg[iarg][0])) {
        nsub = utils::inumeric(FLERR, arg[iarg], false, lmp);
        ++iarg;
        if (nsub > 0) continue;
      }
      error->all(FLERR, "Illegal compute reaxff/atom command");
    } else if (strcmp(arg[iarg], "bonds") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "compute reaxff/atom bonds", error);
      store_bonds = utils::logical(FLERR, arg[iarg+1], false, lmp);
      iarg += 2;
    } else error->all(FLERR,"Illegal compute reaxff/atom command");
  }

  local_flag = store_bonds;
}

/* ---------------------------------------------------------------------- */

ComputeReaxFFAtom::~ComputeReaxFFAtom()
{
  memory->destroy(array_local);
  memory->destroy(array_atom);
  memory->destroy(abo);
  memory->destroy(neighid);
  memory->destroy(bondcount);
}

/* ---------------------------------------------------------------------- */

void ComputeReaxFFAtom::init()
{
  if (lmp->suffix_enable) {
    if (lmp->suffix)
      reaxff = dynamic_cast<PairReaxFF *>(force->pair_match(fmt::format("^reax../{}", lmp->suffix), 0, nsub));
    if (!reaxff && lmp->suffix2)
      reaxff = dynamic_cast<PairReaxFF *>(force->pair_match(fmt::format("^reax../{}", lmp->suffix2), 0, nsub));
  }

  if (!reaxff) reaxff = dynamic_cast<PairReaxFF *>(force->pair_match("^reax..", 0, nsub));

  if (!reaxff) error->all(FLERR,"Cannot use compute reaxff/atom without "
                                "pair_style reaxff or reaxff/omp");

  if (reaxff->kokkosable && !kokkosable)
    error->all(FLERR,"Cannot use compute reaxff/atom with pair_style reaxff/kk. Use reaxff/atom/kk.");
}

/* ---------------------------------------------------------------------- */

int ComputeReaxFFAtom::FindBond()
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
  int * mask = atom->mask;
  int numbonds = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      nj = 0;

      for (pj = Start_Index(i, reaxff->api->lists); pj < End_Index(i, reaxff->api->lists); ++pj) {
        bo_ij = &(reaxff->api->lists->select.bond_list[pj]);
        j = bo_ij->nbr;
        if (mask[j] & groupbit) {
          jtag = tag[j];
          bo_tmp = bo_ij->bo_data.BO;

          if (bo_tmp > bo_cut) {
            if (store_bonds) {
              neighid[i][nj] = jtag;
              abo[i][nj] = bo_tmp;
            }
            nj++;
          }
        }
      }
      bondcount[i] = nj;
      numbonds += nj;
    }
  }
  return numbonds;
}

/* ---------------------------------------------------------------------- */

void ComputeReaxFFAtom::compute_bonds()
{
  invoked_bonds = update->ntimestep;

  if (atom->nmax > nmax) {
    memory->destroy(abo);
    memory->destroy(neighid);
    memory->destroy(bondcount);
    memory->destroy(array_atom);
    nmax = atom->nmax;
    if (store_bonds) {
      memory->create(abo, nmax, MAXREAXBOND, "reaxff/atom:abo");
      memory->create(neighid, nmax, MAXREAXBOND, "reaxff/atom:neighid");
    }
    memory->create(bondcount, nmax, "reaxff/atom:bondcount");
    memory->create(array_atom, nmax, 3, "reaxff/atom:array_atom");
  }

  const int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    bondcount[i] = 0;
    for (int j = 0; store_bonds && j < MAXREAXBOND; j++) {
      neighid[i][j] = 0;
      abo[i][j] = 0.0;
    }
  }

  nbonds = FindBond();
}

/* ---------------------------------------------------------------------- */

void ComputeReaxFFAtom::compute_local()
{
  invoked_local = update->ntimestep;

  if (invoked_bonds < update->ntimestep)
    compute_bonds();

  if (nbonds > prev_nbonds) {
    // grow array_local
    memory->destroy(array_local);
    memory->create(array_local, nbonds, 3, "reaxff/atom:array_local");
    prev_nbonds = nbonds;
  }

  size_local_rows = nbonds;
  auto tag = atom->tag;

  int b = 0;

  const int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; ++i) {
    const int numbonds = bondcount[i];

    for (int k = 0; k < numbonds; k++) {
      auto bond = array_local[b++];
      bond[0] = tag[i];
      bond[1] = neighid[i][k];
      bond[2] = abo[i][k];
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeReaxFFAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  if (invoked_bonds < update->ntimestep) {
    compute_bonds();
  }

  const int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; ++i) {
    auto ptr = array_atom[i];
    ptr[0] = reaxff->api->workspace->total_bond_order[i];
    ptr[1] = reaxff->api->workspace->nlp[i];
    ptr[2] = bondcount[i];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeReaxFFAtom::memory_usage()
{
  double bytes = (double)(nmax*3) * sizeof(double);
  bytes += (double)(nmax) * sizeof(int);
  if (store_bonds) {
    bytes += (double)(2*nmax*MAXREAXBOND) * sizeof(double);
    bytes += (double)(nbonds*3) * sizeof(double);
  }
  return bytes;
}
