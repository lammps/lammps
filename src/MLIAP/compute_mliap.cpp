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
  Compute(lmp, narg, arg), list(NULL), mliap(NULL),
  gradforce(NULL), mliapall(NULL), map(NULL), 
  descriptors(NULL), gamma_row_index(NULL), gamma_col_index(NULL),
  gamma(NULL), egradient(NULL), model(NULL), descriptor(NULL),
  iatommliap(NULL), ielemmliap(NULL), numneighmliap(NULL),  
  jatommliap(NULL), jelemmliap(NULL), graddesc(NULL)
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
  gamma_nnz = model->get_gamma_nnz();
  ndims_force = 3;
  ndims_virial = 6;
  yoffset = nparams*nelements;
  zoffset = 2*yoffset;
  natoms = atom->natoms;
  size_array_rows = 1+ndims_force*natoms+ndims_virial;
  size_array_cols = nparams*nelements+1;
  lastcol = size_array_cols-1;

  size_gradforce = ndims_force*nparams*nelements;

  nmax = 0;
  natomdesc_max = 0;
  natomgamma_max = 0;
  natomneigh_max = 0;
  nneigh_max = 0;

  // create a minimal map, placeholder for more general map

  memory->create(map,atom->ntypes+1,"compute_mliap:map");

  for (int i = 1; i <= atom->ntypes; i++)
    map[i] = i-1;

}

/* ---------------------------------------------------------------------- */

ComputeMLIAP::~ComputeMLIAP()
{

  modify->delete_compute(id_virial);

  memory->destroy(mliap);
  memory->destroy(mliapall);
  memory->destroy(gradforce);

  memory->destroy(map);

  memory->destroy(descriptors);
  memory->destroy(gamma_row_index);
  memory->destroy(gamma_col_index);
  memory->destroy(gamma);
  memory->destroy(egradient);
  memory->destroy(graddesc);

  memory->destroy(iatommliap);
  memory->destroy(ielemmliap);
  memory->destroy(numneighmliap);
  memory->destroy(jatommliap);
  memory->destroy(jelemmliap);

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

  // consistency checks

  if (descriptor->ndescriptors != model->ndescriptors)
    error->all(FLERR,"Incompatible model and descriptor definitions");
  if (descriptor->nelements != model->nelements)
    error->all(FLERR,"Incompatible model and descriptor definitions");
  if (nelements != atom->ntypes)
    error->all(FLERR,"nelements must equal ntypes");

  // allocate memory for global array

  memory->create(mliap,size_array_rows,size_array_cols,
                 "mliap:mliap");
  memory->create(mliapall,size_array_rows,size_array_cols,
                 "mliap:mliapall");
  array = mliapall;

  memory->create(egradient,nelements*nparams,"ComputeMLIAP:egradient");

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

  // grow nmax gradforce array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(gradforce);
    nmax = atom->nmax;
    memory->create(gradforce,nmax,size_gradforce,
                   "mliap:gradforce");
  }

  // clear gradforce array

  for (int i = 0; i < ntotal; i++)
    for (int j = 0; j < size_gradforce; j++) {
      gradforce[i][j] = 0.0;
    }

  // clear global array

  for (int irow = 0; irow < size_array_rows; irow++)
    for (int jcol = 0; jcol < size_array_cols; jcol++)
      mliap[irow][jcol] = 0.0;

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  // grow arrays if necessary

  const int natomdesc = list->inum;
  if (natomdesc_max < natomdesc) {
    memory->grow(descriptors,natomdesc,ndescriptors,"ComputeMLIAP:descriptors");
    natomdesc_max = natomdesc;
  }

  generate_neigharrays();
    
  // compute descriptors
  
  descriptor->compute_descriptors(natomdesc, iatommliap, ielemmliap, numneighmliap, 
                                  jatommliap, jelemmliap, descriptors);
    
  if (gradgradflag) {

  // grow gamma arrays if necessary

    const int natomgamma = list->inum;
    if (natomgamma_max < natomgamma) {
      memory->grow(gamma_row_index,natomgamma,gamma_nnz,"ComputeMLIAP:gamma_row_index");
      memory->grow(gamma_col_index,natomgamma,gamma_nnz,"ComputeMLIAP:gamma_col_index");
      memory->grow(gamma,natomgamma,gamma_nnz,"ComputeMLIAP:gamma");
      natomgamma_max = natomgamma;
    }

    // calculate double gradient w.r.t. parameters and descriptors
    
    model->param_gradient(natomdesc, iatommliap, ielemmliap, descriptors, gamma_row_index, 
                          gamma_col_index, gamma, egradient);
    
    // calculate descriptor gradient contributions to parameter gradients
    
    descriptor->compute_gradients(natomdesc, iatommliap, ielemmliap, numneighmliap, 
                                  jatommliap, jelemmliap, 
                                  gamma_nnz, gamma_row_index, 
                                  gamma_col_index, gamma, gradforce,
                                  yoffset, zoffset);
    
  } else {

    // calculate descriptor gradients
    
    descriptor->compute_descriptor_gradients(natomdesc, iatommliap, ielemmliap, numneighmliap, 
                                             jatommliap, jelemmliap, graddesc);

    // calculate force gradients w.r.t. parameters
    
    model->compute_force_gradients(descriptors, natomdesc, iatommliap, ielemmliap, 
                                   numneighmliap, jatommliap, jelemmliap, graddesc, 
                                   yoffset, zoffset, gradforce, egradient);
    
  }

  // accumulate descriptor gradient contributions to global array

  for (int ielem = 0; ielem < nelements; ielem++) {
    const int elemoffset = nparams*ielem;
    for (int jparam = 0; jparam < nparams; jparam++) {
      int irow = 1;
      for (int i = 0; i < ntotal; i++) {
        double *gradforcei = gradforce[i]+elemoffset;
        int iglobal = atom->tag[i];
        int irow = 3*(iglobal-1)+1;
        mliap[irow][jparam+elemoffset] += gradforcei[jparam];
        mliap[irow+1][jparam+elemoffset] += gradforcei[jparam+yoffset];
        mliap[irow+2][jparam+elemoffset] += gradforcei[jparam+zoffset];
      }
    }
  }

 // copy forces to global array

  for (int i = 0; i < atom->nlocal; i++) {
    int iglobal = atom->tag[i];
    int irow = 3*(iglobal-1)+1;
    mliap[irow][lastcol] = atom->f[i][0];
    mliap[irow+1][lastcol] = atom->f[i][1];
    mliap[irow+2][lastcol] = atom->f[i][2];
  }

  // accumulate bispectrum virial contributions to global array

  dbdotr_compute();

  // copy energy gradient contributions to global array

  for (int ielem = 0; ielem < nelements; ielem++) {
    const int elemoffset = nparams*ielem;
    for (int jparam = 0; jparam < nparams; jparam++)
      mliap[0][jparam+elemoffset] = egradient[jparam+elemoffset];
  }

  // sum up over all processes

  MPI_Allreduce(&mliap[0][0],&mliapall[0][0],size_array_rows*size_array_cols,MPI_DOUBLE,MPI_SUM,world);

  // copy energy to last column

  int irow = 0;
  double reference_energy = c_pe->compute_scalar();
  mliapall[irow++][lastcol] = reference_energy;

  // copy virial stress to last column
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
   generate neighbor arrays
------------------------------------------------------------------------- */

void ComputeMLIAP::generate_neigharrays()
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

void ComputeMLIAP::grow_neigharrays()
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
    for (int ielem = 0; ielem < nelements; ielem++) {
      const int elemoffset = nparams*ielem;
      double *gradforcei = gradforce[i]+elemoffset;
      for (int jparam = 0; jparam < nparams; jparam++) {
        double dbdx = gradforcei[jparam];
        double dbdy = gradforcei[jparam+yoffset];
        double dbdz = gradforcei[jparam+zoffset];
        int irow = irow0;
        mliap[irow++][jparam+elemoffset] += dbdx*x[i][0];
        mliap[irow++][jparam+elemoffset] += dbdy*x[i][1];
        mliap[irow++][jparam+elemoffset] += dbdz*x[i][2];
        mliap[irow++][jparam+elemoffset] += dbdz*x[i][1];
        mliap[irow++][jparam+elemoffset] += dbdz*x[i][0];
        mliap[irow++][jparam+elemoffset] += dbdy*x[i][0];
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
  bytes += nmax*size_gradforce * sizeof(double);          // gradforce
  int n = atom->ntypes+1;
  bytes += n*sizeof(int);        // map

  return bytes;
}
