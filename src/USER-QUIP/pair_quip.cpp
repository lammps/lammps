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

/* ----------------------------------------------------------------------
   Contributing authors: Albert Bartok (Cambridge University)
                         Aidan Thompson (Sandia, athomps@sandia.gov)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_quip.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairQUIP::PairQUIP(LAMMPS *lmp) : Pair(lmp)
{
   single_enable = 0;
   restartinfo = 0;
   one_coeff = 1;
   no_virial_fdotr_compute = 1;
   manybody_flag = 1;

   map = NULL;
   quip_potential = NULL;
   quip_file = NULL;
   quip_string = NULL;
}

PairQUIP::~PairQUIP()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
  }
  delete [] quip_potential;
  delete [] quip_file;
  delete [] quip_string;
}

void PairQUIP::compute(int eflag, int vflag)
{
  int inum, jnum, sum_num_neigh, ii, jj, i, iquip;
  int *ilist;
  int *jlist;
  int *numneigh, **firstneigh;
  int *quip_num_neigh, *quip_neigh, *atomic_numbers;

  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int ntotal = nlocal + nghost;
  int *type = atom->type;
  tagint *tag = atom->tag;

  double **x = atom->x;
  double **f = atom->f;

  double *quip_local_e, *quip_force, *quip_local_virial, *quip_virial, quip_energy, *lattice;

  ev_init(eflag,vflag);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  sum_num_neigh = 0;
  quip_num_neigh = new int [inum];

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    quip_num_neigh[ii] = numneigh[i];
    sum_num_neigh += numneigh[i];
  }

  quip_neigh = new int [sum_num_neigh];
  iquip = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
       quip_neigh[iquip] = (jlist[jj] & NEIGHMASK) + 1;
       iquip++;
    }
  }

  atomic_numbers = new int[ntotal];
  for (ii = 0; ii < ntotal; ii++)
     atomic_numbers[ii] = map[type[ii]-1];

  quip_local_e = new double [ntotal];
  quip_force = new double [ntotal*3];
  quip_local_virial = new double [ntotal*9];
  quip_virial = new double [9];

  lattice = new double [9];
  lattice[0] = domain->xprd;
  lattice[1] = 0.0;
  lattice[2] = 0.0;
  lattice[3] = domain->xy;
  lattice[4] = domain->yprd;
  lattice[5] = 0.0;
  lattice[6] = domain->xz;
  lattice[7] = domain->yz;
  lattice[8] = domain->zprd;

#if defined(LAMMPS_BIGBIG)
  int *tmptag = new int[ntotal];
  int tmplarge = 0, toolarge = 0;
  for (ii = 0; ii < ntotal; ++ii) {
    tmptag[ii] = tag[ii];
    if (tag[ii] > MAXSMALLINT) tmplarge=1;
  }
  MPI_Allreduce(&tmplarge,&toolarge,1,MPI_INT,MPI_MAX,world);
  if (toolarge > 0)
    error->all(FLERR,"Pair style quip does not support 64-bit atom IDs");

  quip_lammps_wrapper(&nlocal,&nghost,atomic_numbers,tmptag,
                      &inum,&sum_num_neigh,ilist,
                      quip_num_neigh,quip_neigh,lattice,
                      quip_potential,&n_quip_potential,&x[0][0],
                      &quip_energy,quip_local_e,quip_virial,
                      quip_local_virial,quip_force);

  delete[] tmptag;
#else
  quip_lammps_wrapper(&nlocal,&nghost,atomic_numbers,tag,
                      &inum,&sum_num_neigh,ilist,
                      quip_num_neigh,quip_neigh,lattice,
                      quip_potential,&n_quip_potential,&x[0][0],
                      &quip_energy,quip_local_e,quip_virial,
                      quip_local_virial,quip_force);
#endif

  iquip = 0;
  for (ii = 0; ii < ntotal; ii++) {
     for( jj = 0; jj < 3; jj++ ) {
        f[ii][jj] += quip_force[iquip];
        iquip++;
     }
  }

  if(eflag_global) {
    eng_vdwl = quip_energy;
  }

  if(eflag_atom) {
    for (ii = 0; ii < ntotal; ii++) {
      eatom[ii] = quip_local_e[ii];
    }
  }

  if (vflag_global) {
      virial[0] = quip_virial[0];
      virial[1] = quip_virial[4];
      virial[2] = quip_virial[8];
      virial[3] = (quip_virial[3] + quip_virial[1])*0.5;
      virial[4] = (quip_virial[2] + quip_virial[6])*0.5;
      virial[5] = (quip_virial[5] + quip_virial[7])*0.5;
  }

  if(vflag_atom) {
    int iatom = 0;
     for(ii = 0; ii < ntotal; ii++) {
       vatom[ii][0] += quip_local_virial[iatom+0];
       vatom[ii][1] += quip_local_virial[iatom+4];
       vatom[ii][2] += quip_local_virial[iatom+8];
       vatom[ii][3] += (quip_local_virial[iatom+3] +
                        quip_local_virial[iatom+1])*0.5;
       vatom[ii][4] += (quip_local_virial[iatom+2] +
                        quip_local_virial[iatom+6])*0.5;
       vatom[ii][5] += (quip_local_virial[iatom+5] +
                        quip_local_virial[iatom+7])*0.5;
       iatom += 9;
     }
  }

  delete [] atomic_numbers;
  delete [] quip_num_neigh;
  delete [] quip_neigh;
  delete [] quip_local_e;
  delete [] quip_force;
  delete [] quip_virial;
  delete [] quip_local_virial;
  delete [] lattice;

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairQUIP::settings(int narg, char ** /* arg */)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

  // check if linked to the correct QUIP library API version
  // as of 2017-07-19 this is API_VERSION 1
  if (quip_lammps_api_version() != 1)
    error->all(FLERR,"QUIP LAMMPS wrapper API version is not compatible "
        "with this version of LAMMPS");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */
void PairQUIP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create(setflag,n+1,n+1,"pair:setflag");
  cutsq = memory->create(cutsq,n+1,n+1,"pair:cutsq");
  map = new int[n];
}

void PairQUIP::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  int n = atom->ntypes;
  if (narg != (4+n)) {
    char str[1024];
    sprintf(str,"Number of arguments %d is not correct, it should be %d", narg, 4+n);
    error->all(FLERR,str);
  }

  // ensure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  n_quip_file = strlen(arg[2]);
  quip_file = new char[n_quip_file+1];
  strcpy(quip_file,arg[2]);

  n_quip_string = strlen(arg[3]);
  quip_string = new char[n_quip_string+1];
  strcpy(quip_string,arg[3]);

  for (int i = 4; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0)
      map[i-3] = -1;
    else
      map[i-3] = force->inumeric(FLERR,arg[i]);
  }

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

  // Initialise potential
  // First call initialises potential via the fortran code in memory, and returns the necessary size
  // of quip_potential. This behaviour is invoked by setting n_potential_quip to 0.
  n_quip_potential = 0;
  quip_potential = new int[0];
  quip_lammps_potential_initialise(quip_potential,&n_quip_potential,&cutoff,quip_file,&n_quip_file,quip_string,&n_quip_string);
  delete [] quip_potential;

  // Allocate quip_potential integer array. This initialise call will transfer the location of the
  // previously initialised potential to the quip_potential variable, and we will use it as a handle
  // when calling the actual calculation routine. We return the cutoff as well.
  quip_potential = new int[n_quip_potential];
  quip_lammps_potential_initialise(quip_potential,&n_quip_potential,&cutoff,quip_file,&n_quip_file,quip_string,&n_quip_string);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairQUIP::init_style()
{
  // Require newton pair on

  if (force->newton_pair != 1)
    error->all(FLERR,"Pair style quip requires newton pair on");

  // Initialise neighbor list
  int irequest_full = neighbor->request(this);
  neighbor->requests[irequest_full]->id = 1;
  neighbor->requests[irequest_full]->half = 0;
  neighbor->requests[irequest_full]->full = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairQUIP::init_one(int /*i*/, int /*j*/)
{
  return cutoff;
}
