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
#include "mliap.h"
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
  Compute(lmp, narg, arg), mliaparray(NULL), 
  mliaparrayall(NULL), map(NULL)
{
  array_flag = 1;
  extarray = 0;

  if (narg < 4)
    error->all(FLERR,"Illegal compute mliap command");

  // default values

  gradgradflag = 1;

  // set flags for required keywords

  int modelflag = 0;
  int descriptorflag = 0;

  // process keywords

  int iarg = 3;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"model") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute mliap command");
      if (strcmp(arg[iarg+1],"linear") == 0) {
	if (iarg+4 > narg) error->all(FLERR,"Illegal compute mliap command");
	int ntmp1 = atoi(arg[iarg+2]);
	int ntmp2 = atoi(arg[iarg+3]);
        model = new MLIAPModelLinear(lmp,ntmp1,ntmp2);
        iarg += 4;
      } else if (strcmp(arg[iarg+1],"quadratic") == 0) {
	if (iarg+4 > narg) error->all(FLERR,"Illegal compute mliap command");
	int ntmp1 = atoi(arg[iarg+2]);
	int ntmp2 = atoi(arg[iarg+3]);
        model = new MLIAPModelQuadratic(lmp,ntmp1,ntmp2);
        iarg += 4;
      } else error->all(FLERR,"Illegal compute mliap command");
      modelflag = 1;
    } else if (strcmp(arg[iarg],"descriptor") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute mliap command");
      if (strcmp(arg[iarg+1],"sna") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal compute mliap command");
        descriptor = new MLIAPDescriptorSNAP(lmp,arg[iarg+2]);
        iarg += 3;
      } else error->all(FLERR,"Illegal compute mliap command");
      descriptorflag = 1;
    } else if (strcmp(arg[iarg],"gradgradflag") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal compute mliap command");
      gradgradflag = atoi(arg[iarg+1]);
      if (gradgradflag != 0 && gradgradflag != 1) 
        error->all(FLERR,"Illegal compute mliap command");
      iarg += 2;
    } else
      error->all(FLERR,"Illegal compute mliap command");
  }

  if (modelflag == 0 || descriptorflag == 0)
    error->all(FLERR,"Illegal compute_style command");

  ndescriptors = descriptor->ndescriptors;
  nparams = model->nparams;
  nelements = model->nelements;

  // create a minimal map, placeholder for more general map

  memory->create(map,atom->ntypes+1,"compute_mliap:map");

  for (int i = 1; i <= atom->ntypes; i++)
    map[i] = i-1;

  mliap = new MLIAP(lmp, ndescriptors, nparams, nelements, gradgradflag, map, model, descriptor);

  size_array_rows = mliap->size_array_rows;
  size_array_cols = mliap->size_array_cols;
  lastcol = size_array_cols-1;
  nmax = 0;
  natomgamma_max = 0;
  printf("Made it to here\n");
}

/* ---------------------------------------------------------------------- */

ComputeMLIAP::~ComputeMLIAP()
{

  modify->delete_compute(id_virial);

  memory->destroy(mliaparray);
  memory->destroy(mliaparrayall);
  memory->destroy(mliap->gradforce);
  memory->destroy(map);

  delete mliap;
  delete model;
  delete descriptor;
}

/* ---------------------------------------------------------------------- */

void ComputeMLIAP::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute mliap requires a pair style be defined");

  if (descriptor->cutmax > force->pair->cutforce)
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

  // initialize model and descriptor

  model->init();
  descriptor->init();
  mliap->init();

  // consistency checks

  if (descriptor->ndescriptors != model->ndescriptors)
    error->all(FLERR,"Incompatible model and descriptor definitions");
  if (descriptor->nelements != model->nelements)
    error->all(FLERR,"Incompatible model and descriptor definitions");
  if (nelements != atom->ntypes)
    error->all(FLERR,"nelements must equal ntypes");

  // allocate memory for global array

  memory->create(mliaparray,size_array_rows,size_array_cols,
                 "mliap:mliaparray");
  memory->create(mliaparrayall,size_array_rows,size_array_cols,
                 "mliap:mliaparrayall");
  array = mliaparrayall;

  // find compute for reference energy

  std::string id_pe = std::string("thermo_pe");
  int ipe = modify->find_compute(id_pe);
  if (ipe == -1)
    error->all(FLERR,"compute thermo_pe does not exist.");
  c_pe = modify->compute[ipe];

  // add compute for reference virial tensor

  id_virial = id + std::string("_press");
  std::string pcmd = id_virial + " all pressure NULL virial";
  modify->add_compute(pcmd);

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

  // clear global array

  for (int irow = 0; irow < size_array_rows; irow++)
    for (int jcol = 0; jcol < size_array_cols; jcol++)
      mliaparray[irow][jcol] = 0.0;

  // grow nmax gradforce array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(mliap->gradforce);
    nmax = atom->nmax;
    memory->create(mliap->gradforce,nmax,mliap->size_gradforce,
                   "mliap:gradforce");
  }

  // clear gradforce array
  
  for (int i = 0; i < ntotal; i++)
    for (int j = 0; j < mliap->size_gradforce; j++) {
      mliap->gradforce[i][j] = 0.0;
    }
  
  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  mliap->generate_neigharrays(list);
    
  // compute descriptors
  
  descriptor->compute_descriptors(mliap->natomdesc, mliap->iatommliap, mliap->ielemmliap, mliap->numneighmliap, 
                                  mliap->jatommliap, mliap->jelemmliap, mliap->descriptors);
    
  if (gradgradflag) {

  // grow gamma arrays if necessary

    const int natomgamma = list->inum;
    if (natomgamma_max < natomgamma) {
      memory->grow(mliap->gamma_row_index,natomgamma,mliap->gamma_nnz,"ComputeMLIAP:gamma_row_index");
      memory->grow(mliap->gamma_col_index,natomgamma,mliap->gamma_nnz,"ComputeMLIAP:gamma_col_index");
      memory->grow(mliap->gamma,natomgamma,mliap->gamma_nnz,"ComputeMLIAP:gamma");
      natomgamma_max = natomgamma;
    }

    // calculate double gradient w.r.t. parameters and descriptors
    
    model->param_gradient(mliap->natomdesc, mliap->iatommliap, mliap->ielemmliap, mliap->descriptors, mliap->gamma_row_index, 
                          mliap->gamma_col_index, mliap->gamma, mliap->egradient);
    
    // calculate descriptor gradient contributions to parameter gradients
    
    descriptor->compute_gradients(mliap->natomdesc, mliap->iatommliap, mliap->ielemmliap, mliap->numneighmliap, 
                                  mliap->jatommliap, mliap->jelemmliap, 
                                  mliap->gamma_nnz, mliap->gamma_row_index, 
                                  mliap->gamma_col_index, mliap->gamma, mliap->gradforce,
                                  mliap->yoffset, mliap->zoffset);
    
  } else {

    // calculate descriptor gradients
    
    descriptor->compute_descriptor_gradients(mliap->natomdesc, mliap->iatommliap, mliap->ielemmliap, mliap->numneighmliap, 
                                             mliap->jatommliap, mliap->jelemmliap, mliap->graddesc);

    // calculate force gradients w.r.t. parameters
    
    model->compute_force_gradients(mliap->descriptors, mliap->natomdesc, mliap->iatommliap, mliap->ielemmliap, 
                                   mliap->numneighmliap, mliap->jatommliap, mliap->jelemmliap, mliap->graddesc, 
                                   mliap->yoffset, mliap->zoffset, mliap->gradforce, mliap->egradient);
    
  }

  // accumulate descriptor gradient contributions to global array

  for (int ielem = 0; ielem < nelements; ielem++) {
    const int elemoffset = nparams*ielem;
    for (int jparam = 0; jparam < nparams; jparam++) {
      int irow = 1;
      for (int i = 0; i < ntotal; i++) {
        double *gradforcei = mliap->gradforce[i]+elemoffset;
        int iglobal = atom->tag[i];
        int irow = 3*(iglobal-1)+1;
        mliaparray[irow][jparam+elemoffset] += gradforcei[jparam];
        mliaparray[irow+1][jparam+elemoffset] += gradforcei[jparam+mliap->yoffset];
        mliaparray[irow+2][jparam+elemoffset] += gradforcei[jparam+mliap->zoffset];
      }
    }
  }

 // copy forces to global array

  for (int i = 0; i < atom->nlocal; i++) {
    int iglobal = atom->tag[i];
    int irow = 3*(iglobal-1)+1;
    mliaparray[irow][mliap->lastcol] = atom->f[i][0];
    mliaparray[irow+1][mliap->lastcol] = atom->f[i][1];
    mliaparray[irow+2][mliap->lastcol] = atom->f[i][2];
  }

  // accumulate bispectrum virial contributions to global array

  dbdotr_compute();

  // copy energy gradient contributions to global array

  for (int ielem = 0; ielem < nelements; ielem++) {
    const int elemoffset = nparams*ielem;
    for (int jparam = 0; jparam < nparams; jparam++)
      mliaparray[0][jparam+elemoffset] = mliap->egradient[jparam+elemoffset];
  }

  // sum up over all processes

  MPI_Allreduce(&mliaparray[0][0],&mliaparrayall[0][0],size_array_rows*size_array_cols,MPI_DOUBLE,MPI_SUM,world);

  // copy energy to last column

  int irow = 0;
  double reference_energy = c_pe->compute_scalar();
  mliaparrayall[irow++][mliap->lastcol] = reference_energy;

  // copy virial stress to last column
  // switch to Voigt notation

  c_virial->compute_vector();
  irow += 3*mliap->natoms;
  mliaparrayall[irow++][mliap->lastcol] = c_virial->vector[0];
  mliaparrayall[irow++][mliap->lastcol] = c_virial->vector[1];
  mliaparrayall[irow++][mliap->lastcol] = c_virial->vector[2];
  mliaparrayall[irow++][mliap->lastcol] = c_virial->vector[5];
  mliaparrayall[irow++][mliap->lastcol] = c_virial->vector[4];
  mliaparrayall[irow++][mliap->lastcol] = c_virial->vector[3];

}

/* ----------------------------------------------------------------------
   compute global virial contributions via summing r_i.dB^j/dr_i over
   own & ghost atoms
------------------------------------------------------------------------- */

  void ComputeMLIAP::dbdotr_compute()
{
  double **x = atom->x;
  int irow0 = 1+mliap->ndims_force*mliap->natoms;

  // sum over bispectrum contributions to forces
  // on all particles including ghosts

  int nall = atom->nlocal + atom->nghost;
  for (int i = 0; i < nall; i++)
    for (int ielem = 0; ielem < nelements; ielem++) {
      const int elemoffset = nparams*ielem;
      double *gradforcei = mliap->gradforce[i]+elemoffset;
      for (int jparam = 0; jparam < nparams; jparam++) {
        double dbdx = gradforcei[jparam];
        double dbdy = gradforcei[jparam+mliap->yoffset];
        double dbdz = gradforcei[jparam+mliap->zoffset];
        int irow = irow0;
        mliaparray[irow++][jparam+elemoffset] += dbdx*x[i][0];
        mliaparray[irow++][jparam+elemoffset] += dbdy*x[i][1];
        mliaparray[irow++][jparam+elemoffset] += dbdz*x[i][2];
        mliaparray[irow++][jparam+elemoffset] += dbdz*x[i][1];
        mliaparray[irow++][jparam+elemoffset] += dbdz*x[i][0];
        mliaparray[irow++][jparam+elemoffset] += dbdy*x[i][0];
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ComputeMLIAP::memory_usage()
{

  double bytes = size_array_rows*size_array_cols *
    sizeof(double);                                     // mliaparray
  bytes += size_array_rows*size_array_cols *
    sizeof(double);                                     // mliaparrayall
  bytes += nmax*mliap->size_gradforce * sizeof(double);          // gradforce
  int n = atom->ntypes+1;
  bytes += n*sizeof(int);        // map

  return bytes;
}
