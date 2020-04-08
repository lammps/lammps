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
#include <cstdlib>
#include "mliap_model_linear.h"
#include "mliap_model_quadratic.h"
#include "mliap_descriptor_snap.h"
#include "compute_mliap.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{SCALAR,VECTOR,ARRAY};

ComputeMLIAP::ComputeMLIAP(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), cutsq(NULL), list(NULL), mliap(NULL),
  mliap_peratom(NULL), mliapall(NULL)
{

  array_flag = 1;
  extarray = 0;

  int ntypes = atom->ntypes;

  if (narg < 4)
    error->all(FLERR,"Illegal compute mliap command");

  // process keywords

  int iarg = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"model") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute mliap command");
      if (strcmp(arg[iarg+1],"linear") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal compute mliap command");
        model = new MLIAPModelLinear(lmp,arg[iarg+2]);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"quadratic") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal compute mliap command");
        model = new MLIAPModelQuadratic(lmp,arg[iarg+2]);
        iarg += 3;
      } else error->all(FLERR,"Illegal compute mliap command");
    } else if (strcmp(arg[iarg],"descriptor") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute mliap command");
      if (strcmp(arg[iarg+1],"sna") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal compute mliap command");
        descriptor = new MLIAPDescriptorSNAP(lmp,arg[iarg+2]);
        iarg += 3;
      } else error->all(FLERR,"Illegal compute mliap command");
    }
  }

  nparams = model->nparams;
  nperdim = nparams;
  ndims_force = 3;
  ndims_virial = 6;
  yoffset = nperdim;
  zoffset = 2*nperdim;
  natoms = atom->natoms;
  size_array_rows = 1+ndims_force*natoms+ndims_virial;
  size_array_cols = nperdim*atom->ntypes+1;
  lastcol = size_array_cols-1;

  ndims_peratom = ndims_force;
  size_peratom = ndims_peratom*nperdim*atom->ntypes;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeMLIAP::~ComputeMLIAP()
{
  memory->destroy(mliap);
  memory->destroy(mliapall);
  memory->destroy(mliap_peratom);
  memory->destroy(cutsq);

  memory->destroy(map);
}

/* ---------------------------------------------------------------------- */

void ComputeMLIAP::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute mliap requires a pair style be defined");

  if (descriptor->get_cutmax() > force->pair->cutforce)
    error->all(FLERR,"Compute mliap cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"mliap") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute mliap");

  // allocate memory for global array

  memory->create(mliap,size_array_rows,size_array_cols,
                 "mliap:mliap");
  memory->create(mliapall,size_array_rows,size_array_cols,
                 "mliap:mliapall");
  array = mliapall;

  // find compute for reference energy

  char *id_pe = (char *) "thermo_pe";
  int ipe = modify->find_compute(id_pe);
  if (ipe == -1)
    error->all(FLERR,"compute thermo_pe does not exist.");
  c_pe = modify->compute[ipe];

  // add compute for reference virial tensor

  char *id_virial = (char *) "mliap_press";
  char **newarg = new char*[5];
  newarg[0] = id_virial;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pressure";
  newarg[3] = (char *) "NULL";
  newarg[4] = (char *) "virial";
  modify->add_compute(5,newarg);
  delete [] newarg;

  int ivirial = modify->find_compute(id_virial);
  if (ivirial == -1)
    error->all(FLERR,"compute mliap_press does not exist.");
  c_virial = modify->compute[ivirial];

}


/* ---------------------------------------------------------------------- */

void ComputeMLIAP::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeMLIAP::compute_array()
{
  int ntotal = atom->nlocal + atom->nghost;

  invoked_array = update->ntimestep;

  // grow mliap_peratom array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(mliap_peratom);
    nmax = atom->nmax;
    memory->create(mliap_peratom,nmax,size_peratom,
                   "mliap:mliap_peratom");
  }

  if (gamma_max < list->inum) {
    memory->grow(descriptors,list->inum,ndescriptors,"PairMLIAP:descriptors");
    memory->grow(gamma,nparams,list->inum,ndescriptors,"PairMLIAP:gamma");
    gamma_max = list->inum;
  }

  // clear global array

  for (int irow = 0; irow < size_array_rows; irow++)
    for (int icoeff = 0; icoeff < size_array_cols; icoeff++)
      mliap[irow][icoeff] = 0.0;

  // clear local peratom array

  for (int i = 0; i < ntotal; i++)
    for (int icoeff = 0; icoeff < size_peratom; icoeff++) {
      mliap_peratom[i][icoeff] = 0.0;
    }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  // compute descriptors, if needed

  if (model->nonlinearflag)
    descriptor->forward(map, list, descriptors);

  // calculate descriptor contributions to parameter gradients
  // and gamma = double gradient w.r.t. parameters and descriptors

  // i.e. gamma = d2E/d\sigma.dB_i
  // sigma is a parameter and B_i is a descriptor of atom i
  // for SNAP, this is a sparse nparams*natoms*ndescriptors matrix,
  // but in general it could be fully dense. 
 
  // *******Not implemented yet*****************
  // This should populate the energy row and gamma
  // For the linear model energy row will look just like the Bi accumulation
  // in ComputeSNAP i.e. accumulating the intput descriptors vector, 
  // while gamma will be just 1's and 0's
  // For the quadratic model, the energy row will be similar,
  // while gamma will be 1's, 0's and Bi's  

  //   model->param_gradient(list, descriptors, mliap[0], gamma);

  // calculate descriptor gradient contributions to parameter gradients

  // *******Not implemented yet*****************
  // This will just take gamma and multiply it with
  // descriptor gradient contributions i.e. dblist
  // this will resemble snadi accumualation in ComputeSNAP
  // descriptor->param_backward(list, gamma, snadi);

  // accumulate descriptor gradient contributions to global array

  for (int itype = 0; itype < atom->ntypes; itype++) {
    const int typeoffset_local = ndims_peratom*nperdim*itype;
    const int typeoffset_global = nperdim*itype;
    for (int icoeff = 0; icoeff < nperdim; icoeff++) {
      int irow = 1;
      for (int i = 0; i < ntotal; i++) {
        double *snadi = mliap_peratom[i]+typeoffset_local;
        int iglobal = atom->tag[i];
        int irow = 3*(iglobal-1)+1;
        mliap[irow][icoeff+typeoffset_global] += snadi[icoeff];
        mliap[irow+1][icoeff+typeoffset_global] += snadi[icoeff+yoffset];
        mliap[irow+2][icoeff+typeoffset_global] += snadi[icoeff+zoffset];
      }
    }
  }

 // accumulate forces to global array

  for (int i = 0; i < atom->nlocal; i++) {
    int iglobal = atom->tag[i];
    int irow = 3*(iglobal-1)+1;
    mliap[irow][lastcol] = atom->f[i][0];
    mliap[irow+1][lastcol] = atom->f[i][1];
    mliap[irow+2][lastcol] = atom->f[i][2];
  }

  // accumulate bispectrum virial contributions to global array

  dbdotr_compute();

  // sum up over all processes

  MPI_Allreduce(&mliap[0][0],&mliapall[0][0],size_array_rows*size_array_cols,MPI_DOUBLE,MPI_SUM,world);

  // assign energy to last column

  int irow = 0;
  double reference_energy = c_pe->compute_scalar();
  mliapall[irow++][lastcol] = reference_energy;

  // assign virial stress to last column
  // switch to Voigt notation

  c_virial->compute_vector();
  irow += 3*natoms;
  mliapall[irow++][lastcol] = c_virial->vector[0];
  mliapall[irow++][lastcol] = c_virial->vector[1];
  mliapall[irow++][lastcol] = c_virial->vector[2];
  mliapall[irow++][lastcol] = c_virial->vector[5];
  mliapall[irow++][lastcol] = c_virial->vector[4];
  mliapall[irow++][lastcol] = c_virial->vector[3];

}

/* ----------------------------------------------------------------------
   compute global virial contributions via summing r_i.dB^j/dr_i over
   own & ghost atoms
------------------------------------------------------------------------- */

void ComputeMLIAP::dbdotr_compute()
{
  double **x = atom->x;
  int irow0 = 1+ndims_force*natoms;

  // sum over bispectrum contributions to forces
  // on all particles including ghosts

  int nall = atom->nlocal + atom->nghost;
  for (int i = 0; i < nall; i++)
    for (int itype = 0; itype < atom->ntypes; itype++) {
      const int typeoffset_local = ndims_peratom*nperdim*itype;
      const int typeoffset_global = nperdim*itype;
      double *snadi = mliap_peratom[i]+typeoffset_local;
      for (int icoeff = 0; icoeff < nperdim; icoeff++) {
        double dbdx = snadi[icoeff];
        double dbdy = snadi[icoeff+yoffset];
        double dbdz = snadi[icoeff+zoffset];
        int irow = irow0;
        mliap[irow++][icoeff+typeoffset_global] += dbdx*x[i][0];
        mliap[irow++][icoeff+typeoffset_global] += dbdy*x[i][1];
        mliap[irow++][icoeff+typeoffset_global] += dbdz*x[i][2];
        mliap[irow++][icoeff+typeoffset_global] += dbdz*x[i][1];
        mliap[irow++][icoeff+typeoffset_global] += dbdz*x[i][0];
        mliap[irow++][icoeff+typeoffset_global] += dbdy*x[i][0];
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ComputeMLIAP::memory_usage()
{

  double bytes = size_array_rows*size_array_cols *
    sizeof(double);                                     // mliap
  bytes += size_array_rows*size_array_cols *
    sizeof(double);                                     // mliapall
  bytes += nmax*size_peratom * sizeof(double);          // mliap_peratom
  int n = atom->ntypes+1;
  bytes += n*sizeof(int);        // map

  return bytes;
}
