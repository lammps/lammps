// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_bond_history.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "string.h"

#include <cstring>
#include <mpi.h>

using namespace LAMMPS_NS;
using namespace FixConst;

#define LB_FACTOR 1.5
#define DELTA 8192

/* ---------------------------------------------------------------------- */

FixBondHistory::FixBondHistory(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal fix bond/history command");
  update_flag = utils::inumeric(FLERR,arg[3],false,lmp);
  ndata = utils::inumeric(FLERR,arg[4],false,lmp);
  nbond = atom->bond_per_atom;

  if (nbond == 0)
    error->all(FLERR, "Cannot store bond variables without any bonds");

  stored_flag = false;
  restart_global = 1;
  create_attribute = 1;

  bondstore = nullptr;
  maxbond = 0;
  allocate();

  new_fix_id = nullptr;
  array_id = nullptr;
}

/* ---------------------------------------------------------------------- */

FixBondHistory::~FixBondHistory()
{
  if (new_fix_id && modify->nfix) modify->delete_fix(new_fix_id);
  delete [] new_fix_id;
  delete [] array_id;

  memory->destroy(bondstore);
}

/* ---------------------------------------------------------------------- */

int FixBondHistory::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= POST_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBondHistory::post_constructor()
{
  // Store saved bond quantities for each atom using fix property atom

  char *id_fix = utils::strdup(id + std::string("_FIX_PROP_ATOM"));
  char *id_array = utils::strdup(std::string("d2_") + id);
  modify->add_fix(fmt::format("{} {} property/atom {} {}",
                                id_fix, group->names[igroup], id_array, nbond*ndata));
  int tmp1, tmp2;
  index = atom->find_custom(&id_array[3],tmp1,tmp2);
}

/* ---------------------------------------------------------------------- */

void FixBondHistory::update_atom_value(int i, int m, int idata, double value)
{
  if (idata >= ndata || m > nbond) error->all(FLERR, "Index exceeded in fix bond history");
  atom->darray[index][i][m*ndata+idata] = value;
}

/* ---------------------------------------------------------------------- */

double FixBondHistory::get_atom_value(int i, int m, int idata)
{
   if (idata >= ndata || m > nbond) error->all(FLERR, "Index exceeded in fix bond history");
   return atom->darray[index][i][m*ndata+idata];
}

/* ----------------------------------------------------------------------
   If stored values are updated, need to copy to atom arrays before exchanging
   If bondstore array has been allocated, call before nlist rebuild
   Also call prior to irregular communication in other fixes (e.g. deform)
------------------------------------------------------------------------- */

void FixBondHistory::pre_exchange()
{
  if (!update_flag) return;
  if (!stored_flag) return;
  if (!bondstore) return;

  int i1, i2, n, m, idata;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  double **stored = atom->darray[index];

  int nlocal = atom->nlocal;
  tagint **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;
  tagint *tag = atom->tag;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];

    // skip bond if already broken
    if (bondlist[n][2] <= 0) {
      continue;
    }

    if (i1 < nlocal) {
      for (m = 0; m < num_bond[i1]; m++) {
        if (bond_atom[i1][m] == tag[i2]) {
          for (idata = 0; idata < ndata; idata++) {
            stored[i1][m*ndata+idata] = bondstore[n][idata];
          }
        }
      }
    }

    if (i2 < nlocal) {
      for (m = 0; m < num_bond[i2]; m++) {
        if (bond_atom[i2][m] == tag[i1]) {
          for (idata = 0; idata < ndata; idata++) {
            stored[i2][m*ndata+idata] = bondstore[n][idata];
          }
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBondHistory::allocate()
{
  //Ideally would just ask ntopo for maxbond, but protected
  if (comm->nprocs == 1) maxbond = atom->nbonds;
  else maxbond = static_cast<int> (LB_FACTOR * atom->nbonds / comm->nprocs);
  memory->create(bondstore,maxbond,ndata,"fix_bond_store:bondstore");
}

/* ---------------------------------------------------------------------- */

void FixBondHistory::setup_post_neighbor()
{
  //Grow array if number of bonds has increased
  while (neighbor->nbondlist >= maxbond) {
    maxbond += DELTA;
    memory->grow(bondstore,maxbond,ndata,"fix_bond_store:bondstore");
  }
  
  pre_exchange();
  post_neighbor();
}

/* ----------------------------------------------------------------------
   called after neighbor list is build
   build array of stored bond quantities from fix property atom
------------------------------------------------------------------------- */

void FixBondHistory::post_neighbor()
{
  //Grow array if number of bonds has increased
  while (neighbor->nbondlist >= maxbond) {
    maxbond += DELTA;
    memory->grow(bondstore,maxbond,ndata,"fix_bond_store:bondstore");
  }

  int i1, i2, n, m, idata;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  double **stored = atom->darray[index];

  int nlocal = atom->nlocal;
  tagint **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;
  tagint *tag = atom->tag;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];

    // skip bond if already broken
    if (bondlist[n][2] <= 0) {
      continue;
    }

    if (i1 < nlocal) {
      for (m = 0; m < num_bond[i1]; m++) {
        if (bond_atom[i1][m] == tag[i2]) {
          for (idata = 0; idata < ndata; idata++) {
            bondstore[n][idata] = stored[i1][m*ndata+idata];
          }
        }
      }
    }

    if (i2 < nlocal){
      for (m = 0; m < num_bond[i2]; m++) {
        if (bond_atom[i2][m] == tag[i1]) {
          for (idata = 0; idata < ndata; idata++) {
            bondstore[n][idata] = stored[i2][m*ndata+idata];
          }
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixBondHistory::memory_usage()
{
  return maxbond * ndata * sizeof(double);
}

/* ---------------------------------------------------------------------- */

void FixBondHistory::write_restart(FILE *fp)
{
  int n = 0;
  double list[1];
  list[n++] = stored_flag;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ---------------------------------------------------------------------- */

void FixBondHistory::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  stored_flag = static_cast<int> (list[n++]);
}

/* ----------------------------------------------------------------------
   initialize bond values to zero, called when atom is created
------------------------------------------------------------------------- */

void FixBondHistory::set_arrays(int i)
{
  double **stored = atom->darray[index];
  for (int m = 0; m < nbond; m++)
    for (int idata = 0; idata < ndata; idata++)
      stored[i][m*ndata+idata] = 0.0;
}

/* ----------------------------------------------------------------------
   Remove all data for row by compressing data - moving last element
------------------------------------------------------------------------- */

void FixBondHistory::delete_bond(int i, int m)
{
  double **stored = atom->darray[index];
  int n = atom->num_bond[i];
  if (m != n-1) shift_bond(i, m, n-1);

  for (int idata = 0; idata < ndata; idata ++)
    stored[i][(n-1)*ndata+idata] = 0.0;
}

/* ----------------------------------------------------------------------
   Shift bond data to a new location
   Used by some fixes that delete bonds which don't move last element
------------------------------------------------------------------------- */

void FixBondHistory::shift_bond(int i, int m, int k)
{
  double **stored = atom->darray[index];
  int n = atom->num_bond[i];
  for (int idata = 0; idata < ndata; idata ++)
    stored[i][m*ndata+idata] = stored[i][k*ndata+idata];
}