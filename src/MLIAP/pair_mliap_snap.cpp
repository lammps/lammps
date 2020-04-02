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
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "mliap_model_snap.h"
#include "mliap_descriptor_snap.h"
#include "pair_mliap_snap.h"
#include "atom.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairMLIAPSNAP::PairMLIAPSNAP(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  beta_max = 0;
  beta = NULL;
  bispectrum = NULL;
  atomenergy = NULL;

  model = NULL;
  descriptor = NULL;
}

/* ---------------------------------------------------------------------- */

PairMLIAPSNAP::~PairMLIAPSNAP()
{
  if (copymode) return;

  memory->destroy(beta);
  memory->destroy(bispectrum);
  memory->destroy(atomenergy);

  delete model;
  delete descriptor;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }

}

/* ----------------------------------------------------------------------
   This version is a straightforward implementation
   ---------------------------------------------------------------------- */

void PairMLIAPSNAP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  // resize lists

  if (beta_max < list->inum) {
    memory->grow(beta,list->inum,ncoeff,"PairMLIAPSNAP:beta");
    memory->grow(bispectrum,list->inum,ncoeff,"PairMLIAPSNAP:bispectrum");
    memory->grow(atomenergy,list->inum,"PairMLIAPSNAP:atomenergy");
    beta_max = list->inum;
  }

  // compute descriptors

  if (model->quadraticflag || eflag)
    descriptor->forward(list, bispectrum);

  // compute E_i and beta_i = dE_i/dB_i for all i in list

  model->gradient(list, bispectrum, atomenergy, beta, eflag);

  // calculate force contributions beta_i*dB_i/dR_j
 
  descriptor->backward(list, beta);

  // tally energy contributions

  if (eflag) {
    for (int ii = 0; ii < list->inum; ii++) {
      int i = list->ilist[ii];
      double evdwl = atomenergy[ii];
      ev_tally_full(i,2.0*evdwl,0.0,0.0,0.0,0.0,0.0);
    }

  }

  // calculate stress 

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMLIAPSNAP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMLIAPSNAP::settings(int narg, char ** /* arg */)
{
  if (narg > 0)
    error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMLIAPSNAP::coeff(int narg, char **arg)
{
  if (narg < 5) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  model = new MLIAPModelSNAP(lmp);
  descriptor = new MLIAPDescriptorSNAP(lmp);

  model->init(narg, arg);
  descriptor->init(narg, arg);

  ncoeff = model->ncoeff;
  if (ncoeff != descriptor->ncoeff)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // clear setflag since coeff() called once with I,J = * *

  int n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (descriptor->map[i] >= 0 && descriptor->map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairMLIAPSNAP::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style SNAP requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMLIAPSNAP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  return sqrt(descriptor->cutsq[i][j]);
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double PairMLIAPSNAP::memory_usage()
{
  double bytes = Pair::memory_usage();

  int n = atom->ntypes+1;
  bytes += n*n*sizeof(int);      // setflag
  bytes += n*n*sizeof(double);   // cutsq
  bytes += beta_max*ncoeff*sizeof(double); // bispectrum
  bytes += beta_max*ncoeff*sizeof(double); // beta
  bytes += beta_max*ncoeff*sizeof(double); // energyatom

  bytes += descriptor->memory_usage(); // Descriptor object
  bytes += model->memory_usage();      // Model object

  return bytes;
}

