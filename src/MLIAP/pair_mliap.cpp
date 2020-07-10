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

  beta = NULL;
  descriptors = NULL;

  model = NULL;
  descriptor = NULL;
  map = NULL;

  natomdesc_max = 0;
  natomneigh_max = 0;
  nneigh_max = 0;
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
   MLIAP force calculation
   ---------------------------------------------------------------------- */

void PairMLIAP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);
    
  // grow atom arrays if necessary

  const int natomdesc = list->inum;
  if (natomdesc_max < natomdesc) {
    memory->grow(beta,natomdesc,ndescriptors,"PairMLIAP:beta");
    memory->grow(descriptors,natomdesc,ndescriptors,"PairMLIAP:descriptors");
    natomdesc_max = natomdesc;
  }

  generate_neigharrays();
    
  // compute descriptors, if needed

  if (model->nonlinearflag || eflag)
    descriptor->compute_descriptors(natomdesc, iatommliap, ielemmliap, numneighmliap, 
                                    jatommliap, jelemmliap, descriptors);

  // compute E_i and beta_i = dE_i/dB_i for all i in list

  model->gradient(natomdesc, iatommliap, ielemmliap, descriptors, beta, this, eflag);

  // calculate force contributions beta_i*dB_i/dR_j

  descriptor->compute_forces(natomdesc, iatommliap, ielemmliap, numneighmliap, 
                             jatommliap, jelemmliap, beta, this, vflag);

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

  // set flags for required keywords

  int modelflag = 0;
  int descriptorflag = 0;

  // process keywords

  int iarg = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"model") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style mliap command");
      if (strcmp(arg[iarg+1],"linear") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal pair_style mliap command");
        model = new MLIAPModelLinear(lmp,arg[iarg+2]);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"quadratic") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal pair_style mliap command");
        model = new MLIAPModelQuadratic(lmp,arg[iarg+2]);
        iarg += 3;
      } else error->all(FLERR,"Illegal pair_style mliap command");
      modelflag = 1;
    } else if (strcmp(arg[iarg],"descriptor") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style mliap command");
      if (strcmp(arg[iarg+1],"sna") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal pair_style mliap command");
        descriptor = new MLIAPDescriptorSNAP(lmp,arg[iarg+2]);
        iarg += 3;
      } else error->all(FLERR,"Illegal pair_style mliap command");
      descriptorflag = 1;
    } else
      error->all(FLERR,"Illegal pair_style mliap command");
  }

  if (modelflag == 0 || descriptorflag == 0)
    error->all(FLERR,"Illegal pair_style command");

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
   generate neighbor arrays
------------------------------------------------------------------------- */

void PairMLIAP::generate_neigharrays()
{
  double **x = atom->x;
  int *type = atom->type;
  
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  
  grow_neigharrays();
  
  int ij = 0;
  for (int ii = 0; ii < list->inum; ii++) {
    const int i = list->ilist[ii];
    
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int itype = type[i];
    const int ielem = map[itype];
    
    int *jlist = firstneigh[i];
    const int jnum = numneigh[i];
    
    int ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      const double delx = x[j][0] - xtmp;
      const double dely = x[j][1] - ytmp;
      const double delz = x[j][2] - ztmp;
      const double rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type[j];
      const int jelem = map[jtype];
      
      if (rsq < descriptor->cutsq[ielem][jelem]) {
        jatommliap[ij] = j;
        jelemmliap[ij] = jelem;
        ij++;
        ninside++;
      }
    }
    iatommliap[ii] = i;
    ielemmliap[ii] = ielem;
    numneighmliap[ii] = ninside;
  }
}

/* ----------------------------------------------------------------------
   grow neighbor arrays to handle all neighbors
------------------------------------------------------------------------- */

void PairMLIAP::grow_neigharrays()
{

  // grow neighbor atom arrays if necessary
    
  const int natomneigh = list->inum;
  if (natomneigh_max < natomneigh) {
    memory->grow(iatommliap,natomneigh,"ComputeMLIAP:iatommliap");
    memory->grow(ielemmliap,natomneigh,"ComputeMLIAP:ielemmliap");
    memory->grow(numneighmliap,natomneigh,"ComputeMLIAP:numneighmliap");
    natomneigh_max = natomneigh;
  }

  // grow neighbor arrays if necessary

  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  
  int iilast = list->inum-1;
  int ilast = list->ilist[iilast];
  int upperbound = firstneigh[ilast] - firstneigh[0] + numneigh[ilast];
  if (nneigh_max >= upperbound) return;

  double **x = atom->x;
  int *type = atom->type;
  
  int nneigh = 0;
  for (int ii = 0; ii < list->inum; ii++) {
    const int i = list->ilist[ii];
    
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int itype = type[i];
    const int ielem = map[itype];
    
    int *jlist = firstneigh[i];
    const int jnum = numneigh[i];
    
    int ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      const double delx = x[j][0] - xtmp;
      const double dely = x[j][1] - ytmp;
      const double delz = x[j][2] - ztmp;
      const double rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type[j];
      const int jelem = map[jtype];
      if (rsq < descriptor->cutsq[ielem][jelem]) ninside++;
    }
    nneigh += ninside;
  }
  
  if (nneigh_max < nneigh) {
    memory->grow(jatommliap,nneigh,"ComputeMLIAP:jatommliap");
    memory->grow(jelemmliap,nneigh,"ComputeMLIAP:jelemmliap");
    memory->grow(graddesc,nneigh,ndescriptors,3,"ComputeMLIAP:graddesc");
    nneigh_max = nneigh;
  }
}

/* ----------------------------------------------------------------------
   add energy of atom i to global and per-atom energy
------------------------------------------------------------------------- */

void PairMLIAP::e_tally(int i, double ei)
{
  if (eflag_global) eng_vdwl += ei;
  if (eflag_atom) eatom[i] += ei;
}

/* ----------------------------------------------------------------------
   add virial contribution into global and per-atom accumulators
------------------------------------------------------------------------- */

void PairMLIAP::v_tally(int i, int j,
                        double *fij, double *rij)
{
  double v[6];

  if (vflag_either) {
    v[0] = -rij[0]*fij[0];
    v[1] = -rij[1]*fij[1];
    v[2] = -rij[2]*fij[2];
    v[3] = -rij[0]*fij[1];
    v[4] = -rij[0]*fij[2];
    v[5] = -rij[1]*fij[2];

    if (vflag_global) {
      virial[0] += v[0];
      virial[1] += v[1];
      virial[2] += v[2];
      virial[3] += v[3];
      virial[4] += v[4];
      virial[5] += v[5];
    }

    if (vflag_atom) {
      vatom[i][0] += 0.5*v[0];
      vatom[i][1] += 0.5*v[1];
      vatom[i][2] += 0.5*v[2];
      vatom[i][3] += 0.5*v[3];
      vatom[i][4] += 0.5*v[4];
      vatom[i][5] += 0.5*v[5];

      vatom[j][0] += 0.5*v[0];
      vatom[j][1] += 0.5*v[1];
      vatom[j][2] += 0.5*v[2];
      vatom[j][3] += 0.5*v[3];
      vatom[j][4] += 0.5*v[4];
      vatom[j][5] += 0.5*v[5];
    }
  }
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
  return sqrt(descriptor->cutsq[map[i]][map[j]]);
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double PairMLIAP::memory_usage()
{
  double bytes = Pair::memory_usage();

  int n = atom->ntypes+1;
  bytes += n*n*sizeof(int);      // setflag
  bytes += natomdesc_max*ndescriptors*sizeof(double); // descriptors
  bytes += natomdesc_max*ndescriptors*sizeof(double); // beta

  bytes += descriptor->memory_usage(); // Descriptor object
  bytes += model->memory_usage();      // Model object

  return bytes;
}

