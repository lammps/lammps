// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "pair_mliap.h"

#include "mliap_data.h"
#include "mliap_descriptor_snap.h"
#include "mliap_descriptor_so3.h"
#include "mliap_model_linear.h"
#include "mliap_model_nn.h"
#include "mliap_model_quadratic.h"
#ifdef MLIAP_PYTHON
#include "mliap_unified.h"
#include "mliap_model_python.h"
#endif

#include "atom.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairMLIAP::PairMLIAP(LAMMPS *lmp) :
    Pair(lmp), map(nullptr), model(nullptr), descriptor(nullptr), data(nullptr)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  is_child = false;
  centroidstressflag = CENTROID_NOTAVAIL;
  model=nullptr;
  descriptor=nullptr;
}

/* ---------------------------------------------------------------------- */

PairMLIAP::~PairMLIAP()
{
  if (copymode) return;

  delete model;
  delete descriptor;
  delete data;
  model=nullptr;
  descriptor=nullptr;
  data=nullptr;
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cutghost);
    memory->destroy(map);
  }
}

/* ----------------------------------------------------------------------
   MLIAP force calculation
   ---------------------------------------------------------------------- */

void PairMLIAP::compute(int eflag, int vflag)
{

  // consistency checks
  if (data->ndescriptors != model->ndescriptors)
    error->all(FLERR, "Inconsistent model and descriptor descriptor count: {} vs {}",
               model->ndescriptors, data->ndescriptors);

  if (data->nelements != model->nelements)
    error->all(FLERR, "Inconsistent model and descriptor element count: {} vs {}",
               model->nelements, data->nelements);

  ev_init(eflag, vflag);
  data->generate_neighdata(list, eflag, vflag);

  // compute descriptors, if needed

  if (model->nonlinearflag || eflag) descriptor->compute_descriptors(data);

  // compute E_i and beta_i = dE_i/dB_i for all i in list

  model->compute_gradients(data);

  // calculate force contributions beta_i*dB_i/dR_j

  descriptor->compute_forces(data);
  e_tally(data);

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
  memory->create(cutghost,n+1,n+1,"pair:cutghost");
  memory->create(map,n+1,"pair:map");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMLIAP::settings(int narg, char ** arg)
{

  // This is needed because the unit test calls settings twice
  if (!is_child) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "pair_style mliap", error);
    delete model;
    model = nullptr;
    delete descriptor;
    descriptor = nullptr;
  }

  // process keywords
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"model") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "pair_style mliap model", error);
      if (model != nullptr) error->all(FLERR,"Illegal multiple pair_style mliap model definition");
      if (strcmp(arg[iarg+1],"linear") == 0) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "pair_style mliap model linear", error);
        model = new MLIAPModelLinear(lmp,arg[iarg+2]);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"quadratic") == 0) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "pair_style mliap model quadratic", error);
        model = new MLIAPModelQuadratic(lmp,arg[iarg+2]);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"nn") == 0) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "pair_style mliap model nn", error);
        model = new MLIAPModelNN(lmp,arg[iarg+2]);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"mliappy") == 0) {
#ifdef MLIAP_PYTHON
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "pair_style mliap mliappy", error);
        model = new MLIAPModelPython(lmp,arg[iarg+2]);
        iarg += 3;
#else
        error->all(FLERR,"Using pair_style mliap model mliappy requires ML-IAP with python support");
#endif
      } else error->all(FLERR,"Unknown pair_style mliap model keyword: {}", arg[iarg]);
    } else if (strcmp(arg[iarg],"descriptor") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "pair_style mliap descriptor", error);
      if (descriptor != nullptr) error->all(FLERR,"Illegal multiple pair_style mliap descriptor definition");
      if (strcmp(arg[iarg+1],"sna") == 0) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "pair_style mliap descriptor sna", error);
        descriptor = new MLIAPDescriptorSNAP(lmp,arg[iarg+2]);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"so3") == 0) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "pair_style mliap descriptor so3", error);
        descriptor = new MLIAPDescriptorSO3(lmp,arg[iarg+2]);
        iarg += 3;

      } else error->all(FLERR,"Illegal pair_style mliap command");
    } else if (strcmp(arg[iarg], "unified") == 0) {
#ifdef MLIAP_PYTHON
      if (model != nullptr) error->all(FLERR,"Illegal multiple pair_style mliap model definitions");
      if (descriptor != nullptr) error->all(FLERR,"Illegal multiple pair_style mliap descriptor definitions");
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "pair_style mliap unified", error);
      MLIAPBuildUnified_t build = build_unified(arg[iarg+1], data, lmp);
      if (iarg+3 > narg) {
        ghostneigh = 0;
      } else {
        ghostneigh = utils::logical(FLERR, arg[iarg+2], false, lmp);
      }

      iarg += 3;
      model = build.model;
      descriptor = build.descriptor;
#else
      error->all(FLERR,"Using pair_style mliap unified requires ML-IAP with python support");
#endif
    } else
      error->all(FLERR,"Unknown pair_style mliap keyword: {}", arg[iarg]);
  }

  if (model == nullptr || descriptor == nullptr)
    error->all(FLERR,"Incomplete pair_style mliap setup: need model and descriptor, or unified");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMLIAP::coeff(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  char** elemtypes = &arg[2];

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

  // set up model, descriptor, and mliap data structures

  model->init();
  descriptor->init();
  constexpr int gradgradflag = -1;
  delete data;
  data = new MLIAPData(lmp, gradgradflag, map, model, descriptor, this);
  data->init();
}

/* ----------------------------------------------------------------------
   add energies to eng_vdwl and per-atom energy
------------------------------------------------------------------------- */

void PairMLIAP::e_tally(MLIAPData* data)
{
  if (eflag_global) eng_vdwl += data->energy;
  if (eflag_atom)
    for (int ii = 0; ii < data->nlistatoms; ii++)
      eatom[data->iatoms[ii]] += data->eatoms[ii];
}

/* ----------------------------------------------------------------------
   add virial contribution into global and per-atom accumulators
------------------------------------------------------------------------- */

void PairMLIAP::v_tally(int i, int j, double *fij, double *rij)
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

  if (ghostneigh == 1) {
    neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_GHOST);
  } else {
    neighbor->add_request(this, NeighConst::REQ_FULL);
  }
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMLIAP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  double cutmax = sqrt(descriptor->cutsq[map[i]][map[j]]);
  cutghost[i][j] = cutghost[j][i] = 2.0 * cutmax + neighbor->skin;
  return cutmax;
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double PairMLIAP::memory_usage()
{
  double bytes = Pair::memory_usage();

  int n = atom->ntypes+1;
  bytes += (double)n*n*sizeof(int);            // setflag
  bytes += (double)n*n*sizeof(int);            // cutsq
  bytes += (double)n*n*sizeof(int);            // cutghost
  bytes += (double)n*sizeof(int);              // map
  bytes += descriptor->memory_usage(); // Descriptor object
  bytes += model->memory_usage();      // Model object
  bytes += data->memory_usage();       // Data object

  return bytes;
}

