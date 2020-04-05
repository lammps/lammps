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
#include "mliap_model_linear.h"
#include "mliap_model_quadratic.h"
#include "mliap_descriptor_snap.h"
#include "pair_mliap.h"
#include "atom.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairMLIAP::PairMLIAP(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  beta_max = 0;
  beta = NULL;
  descriptors = NULL;

  model = NULL;
  descriptor = NULL;
}

/* ---------------------------------------------------------------------- */

PairMLIAP::~PairMLIAP()
{
  if (copymode) return;

  memory->destroy(beta);
  memory->destroy(descriptors);

  delete model;
  delete descriptor;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(map);
  }

}

/* ----------------------------------------------------------------------
   This version is a straightforward implementation
   ---------------------------------------------------------------------- */

void PairMLIAP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  // resize lists

  if (beta_max < list->inum) {
    memory->grow(beta,list->inum,ndescriptors,"PairMLIAP:beta");
    memory->grow(descriptors,list->inum,ndescriptors,"PairMLIAP:descriptors");
    beta_max = list->inum;
  }

  // compute descriptors, if needed

  if (model->nonlinearflag || eflag)
    descriptor->forward(list, descriptors);

  // compute E_i and beta_i = dE_i/dB_i for all i in list

  model->gradient(list, descriptors, beta, eflag);

  // calculate force contributions beta_i*dB_i/dR_j
 
  descriptor->backward(list, beta, vflag);

  // calculate stress 

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMLIAP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(map,n+1,"pair:map");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMLIAP::settings(int narg, char ** arg)
{
  if (narg < 4)
    error->all(FLERR,"Illegal pair_style command");

  // process keywords

  int iarg = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"model") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style mliap command");
      if (strcmp(arg[iarg+1],"linear") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal pair_style mliap command");
        model = new MLIAPModelLinear(lmp,arg[iarg+2],this);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"quadratic") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal pair_style mliap command");
        model = new MLIAPModelQuadratic(lmp,arg[iarg+2],this);
        iarg += 3;
      } else error->all(FLERR,"Illegal pair_style mliap command");
    } else if (strcmp(arg[iarg],"descriptor") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style mliap command");
      if (strcmp(arg[iarg+1],"sna") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal pair_style mliap command");
        descriptor = new MLIAPDescriptorSNAP(lmp,arg[iarg+2],this);
        iarg += 3;
      } else error->all(FLERR,"Illegal pair_style mliap command");
    }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMLIAP::coeff(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  char* type1 = arg[0];
  char* type2 = arg[1];
  char** elemtypes = &arg[2];

  // insure I,J args are * *

  if (strcmp(type1,"*") != 0 || strcmp(type2,"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom types to elements
  // map[i] = which element the Ith atom type is, -1 if not mapped
  // map[0] is not used

  for (int i = 1; i <= atom->ntypes; i++) {
    char* elemname = elemtypes[i-1];
    int jelem;
    for (jelem = 0; jelem < descriptor->nelements; jelem++)
      if (strcmp(elemname,descriptor->elements[jelem]) == 0)
        break;

    if (jelem < descriptor->nelements)
      map[i] = jelem;
    else if (strcmp(elemname,"NULL") == 0) map[i] = -1;
    else error->all(FLERR,"Incorrect args for pair coefficients");
  }

  // clear setflag since coeff() called once with I,J = * *

  int n = atom->ntypes;
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

  model->init();
  descriptor->init();

  // consistency checks

  ndescriptors = descriptor->ndescriptors;
  if (ndescriptors != model->ndescriptors)
    error->all(FLERR,"Incompatible model and descriptor definitions");
  if (descriptor->nelements != model->nelements)
    error->all(FLERR,"Incompatible model and descriptor definitions");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairMLIAP::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style MLIAP requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMLIAP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  return descriptor->get_cutoff(map[i],map[j]);
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double PairMLIAP::memory_usage()
{
  double bytes = Pair::memory_usage();

  int n = atom->ntypes+1;
  bytes += n*n*sizeof(int);      // setflag
  bytes += n*n*sizeof(double);   // cutsq
  bytes += beta_max*ndescriptors*sizeof(double); // descriptors
  bytes += beta_max*ndescriptors*sizeof(double); // beta

  bytes += descriptor->memory_usage(); // Descriptor object
  bytes += model->memory_usage();      // Model object

  return bytes;
}

