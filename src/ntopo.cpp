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

#include <mpi.h>
#include "ntopo.h"
#include "atom.h"
#include "neighbor.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define LB_FACTOR 1.5

/* ---------------------------------------------------------------------- */

NTopo::NTopo(LAMMPS *lmp) : Pointers(lmp)
{
  me = comm->me;
  nprocs = comm->nprocs;

  nbondlist = nanglelist = ndihedrallist = nimproperlist = 0;
  maxbond = maxangle = maxdihedral = maximproper = 0;
  bondlist = anglelist = dihedrallist = improperlist = NULL;

  cluster_check = neighbor->cluster_check;
}

/* ---------------------------------------------------------------------- */

NTopo::~NTopo()
{
  memory->destroy(bondlist);
  memory->destroy(anglelist);
  memory->destroy(dihedrallist);
  memory->destroy(improperlist);
}

/* ---------------------------------------------------------------------- */

void NTopo::allocate_bond()
{
  if (nprocs == 1) maxbond = atom->nbonds;
  else maxbond = static_cast<int> (LB_FACTOR * atom->nbonds / nprocs);
  memory->create(bondlist,maxbond,3,"neigh_topo:bondlist");
}

/* ---------------------------------------------------------------------- */

void NTopo::allocate_angle()
{
  if (nprocs == 1) maxangle = atom->nangles;
  else maxangle = static_cast<int> (LB_FACTOR * atom->nangles / nprocs);
  memory->create(anglelist,maxangle,4,"neigh_topo:anglelist");
}

/* ---------------------------------------------------------------------- */

void NTopo::allocate_dihedral()
{
  if (nprocs == 1) maxdihedral = atom->ndihedrals;
  else maxdihedral = static_cast<int> (LB_FACTOR * atom->ndihedrals / nprocs);
  memory->create(dihedrallist,maxdihedral,5,"neigh_topo:dihedrallist");
}

/* ---------------------------------------------------------------------- */

void NTopo::allocate_improper()
{
  if (nprocs == 1) maximproper = atom->nimpropers;
  else maximproper = static_cast<int> (LB_FACTOR * atom->nimpropers / nprocs);
  memory->create(improperlist,maximproper,5,"neigh_topo:improperlist");
}

/* ---------------------------------------------------------------------- */

void NTopo::bond_check()
{
  int i,j;
  double dx,dy,dz,dxstart,dystart,dzstart;

  double **x = atom->x;
  int flag = 0;

  for (int m = 0; m < nbondlist; m++) {
    i = bondlist[m][0];
    j = bondlist[m][1];
    dxstart = dx = x[i][0] - x[j][0];
    dystart = dy = x[i][1] - x[j][1];
    dzstart = dz = x[i][2] - x[j][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
  }

  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all(FLERR,"Bond extent > half of periodic box length");
}

/* ---------------------------------------------------------------------- */

void NTopo::angle_check()
{
  int i,j,k;
  double dx,dy,dz,dxstart,dystart,dzstart;

  double **x = atom->x;
  int flag = 0;

  // check all 3 distances
  // in case angle potential computes any of them

  for (int m = 0; m < nanglelist; m++) {
    i = anglelist[m][0];
    j = anglelist[m][1];
    k = anglelist[m][2];
    dxstart = dx = x[i][0] - x[j][0];
    dystart = dy = x[i][1] - x[j][1];
    dzstart = dz = x[i][2] - x[j][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[i][0] - x[k][0];
    dystart = dy = x[i][1] - x[k][1];
    dzstart = dz = x[i][2] - x[k][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[j][0] - x[k][0];
    dystart = dy = x[j][1] - x[k][1];
    dzstart = dz = x[j][2] - x[k][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
  }

  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all(FLERR,"Angle extent > half of periodic box length");
}

/* ---------------------------------------------------------------------- */

void NTopo::dihedral_check(int nlist, int **list)
{
  int i,j,k,l;
  double dx,dy,dz,dxstart,dystart,dzstart;

  double **x = atom->x;
  int flag = 0;

  // check all 6 distances
  // in case dihedral/improper potential computes any of them

  for (int m = 0; m < nlist; m++) {
    i = list[m][0];
    j = list[m][1];
    k = list[m][2];
    l = list[m][3];
    dxstart = dx = x[i][0] - x[j][0];
    dystart = dy = x[i][1] - x[j][1];
    dzstart = dz = x[i][2] - x[j][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[i][0] - x[k][0];
    dystart = dy = x[i][1] - x[k][1];
    dzstart = dz = x[i][2] - x[k][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[i][0] - x[l][0];
    dystart = dy = x[i][1] - x[l][1];
    dzstart = dz = x[i][2] - x[l][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[j][0] - x[k][0];
    dystart = dy = x[j][1] - x[k][1];
    dzstart = dz = x[j][2] - x[k][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[j][0] - x[l][0];
    dystart = dy = x[j][1] - x[l][1];
    dzstart = dz = x[j][2] - x[l][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[k][0] - x[l][0];
    dystart = dy = x[k][1] - x[l][1];
    dzstart = dz = x[k][2] - x[l][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
  }

  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all)
    error->all(FLERR,"Dihedral/improper extent > half of periodic box length");
}

/* ---------------------------------------------------------------------- */

bigint NTopo::memory_usage()
{
  bigint bytes = 0;
  bytes += 3*maxbond * sizeof(int);
  bytes += 4*maxangle * sizeof(int);
  bytes += 5*maxdihedral * sizeof(int);
  bytes += 5*maximproper * sizeof(int);
  return bytes;
}

