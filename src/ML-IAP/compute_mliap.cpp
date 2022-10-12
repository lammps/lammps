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

/* ----------------------------------------------------------------------
   Contributing author: Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "compute_mliap.h"

#include "mliap_data.h"
#include "mliap_model_linear.h"
#include "mliap_model_quadratic.h"
#include "mliap_descriptor_snap.h"
#include "mliap_descriptor_so3.h"
#ifdef MLIAP_PYTHON
#include "mliap_model_python.h"
#endif

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

enum{SCALAR,VECTOR,ARRAY};

ComputeMLIAP::ComputeMLIAP(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), mliaparray(nullptr),
  mliaparrayall(nullptr), map(nullptr)
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
        model = new MLIAPModelLinear(lmp);
        iarg += 2;
      } else if (strcmp(arg[iarg+1],"quadratic") == 0) {
        model = new MLIAPModelQuadratic(lmp);
        iarg += 2;
      }
#ifdef MLIAP_PYTHON
      else if (strcmp(arg[iarg+1],"mliappy") == 0) {
      model = new MLIAPModelPython(lmp);
      iarg += 2;
      }
#endif
      else error->all(FLERR,"Illegal compute mliap command");
      modelflag = 1;
    } else if (strcmp(arg[iarg],"descriptor") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute mliap command");
      if (strcmp(arg[iarg+1],"sna") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal compute mliap command");
        descriptor = new MLIAPDescriptorSNAP(lmp,arg[iarg+2]);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"so3") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal pair_style mliap command");
        descriptor = new MLIAPDescriptorSO3(lmp,arg[iarg+2]);
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

  // need to tell model how many descriptors
  // so it can figure out how many parameters

  model->set_ndescriptors(descriptor->ndescriptors);

  // create a minimal map, placeholder for more general map

  memory->create(map,atom->ntypes+1,"compute_mliap:map");

  for (int i = 1; i <= atom->ntypes; i++)
    map[i] = i-1;

  data = new MLIAPData(lmp, gradgradflag, map, model, descriptor);

  size_array_rows = data->size_array_rows;
  size_array_cols = data->size_array_cols;
  lastcol = size_array_cols-1;
}

/* ---------------------------------------------------------------------- */

ComputeMLIAP::~ComputeMLIAP()
{

  modify->delete_compute(id_virial);

  memory->destroy(mliaparray);
  memory->destroy(mliaparrayall);
  memory->destroy(map);

  delete data;
  delete model;
  delete descriptor;
}

/* ---------------------------------------------------------------------- */

void ComputeMLIAP::init()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Compute mliap requires a pair style be defined");

  if (descriptor->cutmax > force->pair->cutforce)
    error->all(FLERR,"Compute mliap cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

  if (modify->get_compute_by_style("mliap").size() > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute mliap");

  // initialize model and descriptor

  model->init();
  descriptor->init();
  data->init();

  // consistency checks

  if (data->nelements != atom->ntypes)
    error->all(FLERR,"nelements must equal ntypes");

  // allocate memory for global array

  memory->create(mliaparray,size_array_rows,size_array_cols,
                 "compute_mliap:mliaparray");
  memory->create(mliaparrayall,size_array_rows,size_array_cols,
                 "compute_mliap:mliaparrayall");
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
  int nall = atom->nlocal + atom->nghost;
  invoked_array = update->ntimestep;

  // clear global array

  for (int irow = 0; irow < size_array_rows; irow++)
    for (int jcol = 0; jcol < size_array_cols; jcol++)
      mliaparray[irow][jcol] = 0.0;

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  data->generate_neighdata(list);

  // compute descriptors

  descriptor->compute_descriptors(data);

  if (gradgradflag == 1) {

    // calculate double gradient w.r.t. parameters and descriptors

    model->compute_gradgrads(data);

    // calculate gradients of forces w.r.t. parameters

    descriptor->compute_force_gradients(data);

  } else if (gradgradflag == 0) {

    // calculate descriptor gradients

    descriptor->compute_descriptor_gradients(data);

    // calculate gradients of forces w.r.t. parameters

    model->compute_force_gradients(data);

  } else error->all(FLERR,"Invalid value for gradgradflag");

  // accumulate descriptor gradient contributions to global array

  for (int ielem = 0; ielem < data->nelements; ielem++) {
    const int elemoffset = data->nparams*ielem;
    for (int jparam = 0; jparam < data->nparams; jparam++) {
      for (int i = 0; i < nall; i++) {
        double *gradforcei = data->gradforce[i]+elemoffset;
        tagint irow = 3*(atom->tag[i]-1)+1;
        mliaparray[irow][jparam+elemoffset] += gradforcei[jparam];
        mliaparray[irow+1][jparam+elemoffset] += gradforcei[jparam+data->yoffset];
        mliaparray[irow+2][jparam+elemoffset] += gradforcei[jparam+data->zoffset];
      }
    }
  }

 // copy forces to global array

  for (int i = 0; i < atom->nlocal; i++) {
    int iglobal = atom->tag[i];
    int irow = 3*(iglobal-1)+1;
    mliaparray[irow][lastcol] = atom->f[i][0];
    mliaparray[irow+1][lastcol] = atom->f[i][1];
    mliaparray[irow+2][lastcol] = atom->f[i][2];
  }

  // accumulate bispectrum virial contributions to global array

  dbdotr_compute();

  // copy energy gradient contributions to global array

  for (int ielem = 0; ielem < data->nelements; ielem++) {
    const int elemoffset = data->nparams*ielem;
    for (int jparam = 0; jparam < data->nparams; jparam++)
      mliaparray[0][jparam+elemoffset] = data->egradient[jparam+elemoffset];
  }

  // sum up over all processes

  MPI_Allreduce(&mliaparray[0][0],&mliaparrayall[0][0],size_array_rows*size_array_cols,MPI_DOUBLE,MPI_SUM,world);

  // copy energy to last column

  int irow = 0;
  double reference_energy = c_pe->compute_scalar();
  mliaparrayall[irow++][lastcol] = reference_energy;

  // copy virial stress to last column
  // switch to Voigt notation

  c_virial->compute_vector();
  irow += 3*data->natoms;
  mliaparrayall[irow++][lastcol] = c_virial->vector[0];
  mliaparrayall[irow++][lastcol] = c_virial->vector[1];
  mliaparrayall[irow++][lastcol] = c_virial->vector[2];
  mliaparrayall[irow++][lastcol] = c_virial->vector[5];
  mliaparrayall[irow++][lastcol] = c_virial->vector[4];
  mliaparrayall[irow++][lastcol] = c_virial->vector[3];

}

/* ----------------------------------------------------------------------
   compute global virial contributions via summing r_i.dB^j/dr_i over
   own & ghost atoms
------------------------------------------------------------------------- */

  void ComputeMLIAP::dbdotr_compute()
{
  double **x = atom->x;
  int irow0 = 1+data->ndims_force*data->natoms;

  // sum over bispectrum contributions to forces
  // on all particles including ghosts

  int nall = atom->nlocal + atom->nghost;
  for (int i = 0; i < nall; i++)
    for (int ielem = 0; ielem < data->nelements; ielem++) {
      const int elemoffset = data->nparams*ielem;
      double *gradforcei = data->gradforce[i]+elemoffset;
      for (int jparam = 0; jparam < data->nparams; jparam++) {
        double dbdx = gradforcei[jparam];
        double dbdy = gradforcei[jparam+data->yoffset];
        double dbdz = gradforcei[jparam+data->zoffset];
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

  double bytes = (double)size_array_rows*size_array_cols *
    sizeof(double);                                   // mliaparray
  bytes += (double)size_array_rows*size_array_cols *
    sizeof(double);                                   // mliaparrayall
  int n = atom->ntypes+1;
  bytes += (double)n*sizeof(int);                             // map

  bytes += descriptor->memory_usage(); // Descriptor object
  bytes += model->memory_usage();      // Model object
  bytes += data->memory_usage();       // Data object

  return bytes;
}
