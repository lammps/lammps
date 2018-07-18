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
   Contributing author: Rene Halver (JSC)
------------------------------------------------------------------------- */

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "scafacos.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "memory.h"
#include "error.h"

// ScaFaCoS library

#include <string>
#include <sstream>
#include "fcs.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Scafacos::Scafacos(LAMMPS *lmp, int narg, char **arg) : KSpace(lmp, narg, arg)
{
  if (narg < 2) error->all(FLERR,"Illegal scafacos command");

  int n = strlen(arg[0]) + 1;
  method = new char[n];
  strcpy(method,arg[0]);

  // additional args

  int tflag = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"tolerance") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal scafacos command");
      tflag = 1;
      if (strcmp(arg[iarg+1],"energy") == 0)
        tolerance_type = FCS_TOLERANCE_TYPE_ENERGY;     
      else if (strcmp(arg[iarg+1],"energy_rel") == 0)
        tolerance_type = FCS_TOLERANCE_TYPE_ENERGY_REL;     
      else if (strcmp(arg[iarg+1],"field") == 0)
        tolerance_type = FCS_TOLERANCE_TYPE_FIELD;     
      else if (strcmp(arg[iarg+1],"field_rel") == 0)
        tolerance_type = FCS_TOLERANCE_TYPE_FIELD_REL;     
      else if (strcmp(arg[iarg+1],"potential") == 0)
        tolerance_type = FCS_TOLERANCE_TYPE_POTENTIAL;     
      else if (strcmp(arg[iarg+1],"potential_rel") == 0)
        tolerance_type = FCS_TOLERANCE_TYPE_POTENTIAL_REL;     
      else error->all(FLERR,"Illegal scafacos command");
      tolerance = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else error->all(FLERR,"Illegal scafacos command");
  }

  if (!tflag) error->all(FLERR,"Must set tolerance in scafacos command");

  // initializations

  me = comm->me;
  initialized = 0;

  maxatom = 0;
  epot = NULL;
  efield = NULL;
}

/* ---------------------------------------------------------------------- */

Scafacos::~Scafacos()
{
  delete [] method;

  memory->destroy(epot);
  memory->destroy(efield);

  // NOTE: any clean-up call to ScaFaCoS needed?
}

/* ---------------------------------------------------------------------- */

void Scafacos::init()
{
  // error checks
  // NOTE: allow triclinic at some point

  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use ScaFaCoS with 2d simulation");

  if (domain->triclinic)
    error->all(FLERR,"Cannot use ScaFaCoS with triclinic domain yet");

  // one-time initialization of ScaFaCoS

  scale = 1.0;
  qqrd2e = force->qqrd2e;
  //qsum_qsq();

  if (!initialized) {
    result = fcs_init(&fcs,method,world);
    check_result(result);

    setup_handle();

    result = fcs_set_tolerance(fcs,tolerance_type,tolerance);
    check_result(result);

    if (me == 0) fcs_print_parameters(fcs);

    double **x = atom->x;
    double *q = atom->q; 
    int nlocal = atom->nlocal;

    result = fcs_tune(fcs,nlocal,&x[0][0],q);
  }

  initialized = 1;
}

/* ---------------------------------------------------------------------- */

void Scafacos::setup() {}

/* ---------------------------------------------------------------------- */

void Scafacos::compute(int eflag, int vflag)
{
  double **x = atom->x;
  double *q = atom->q;
  int nlocal = atom->nlocal;

  const double qscale = qqrd2e * scale;

  // if simluation box has changed, call fcs_tune()

  if (eflag || vflag) ev_setup(eflag,vflag);
  else
    eflag_atom = 0;

  if (box_has_changed()) {
    setup_handle();
    result = fcs_tune(fcs,nlocal,&x[0][0],q);
    check_result(result);
  }
  
  // grow epot & efield if necessary

  if (nlocal > maxatom) {
    memory->destroy(epot);
    memory->destroy(efield);
    maxatom = atom->nmax;
    memory->create(epot,maxatom,"scafacos:epot");
    memory->create(efield,maxatom,3,"scafacos:efield");
     
  }

  // initialize epot & efield
  // NOTE: is this necessary?

  for (int i = 0; i < nlocal; i++) {
    epot[i] = 0.0;
    efield[i][0] = efield[i][1] = efield[i][2] = 0.0;
  }
 
  // ScaFaCoS calculation of full Coulombics

  result = fcs_run(fcs,nlocal,&x[0][0],q,&efield[0][0],epot);
  check_result(result);

  // apply Efield to each particle
  // accumulate total energy

  double **f = atom->f;

  double qone;
  double myeng = 0.0;

  for (int i = 0; i < nlocal; i++) {
    qone = q[i] * qscale;
    f[i][0] += qone * efield[i][0];
    f[i][1] += qone * efield[i][1];
    f[i][2] += qone * efield[i][2];
    myeng += 0.5 * qone * epot[i];
  } 

  if (eflag_atom)
  {
    for (int i = 0; i < nlocal; i++)
    {
      eatom[i] = qscale * epot[i];
    }
  }

  MPI_Allreduce(&myeng,&energy,1,MPI_DOUBLE,MPI_SUM,world);
}

/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

double Scafacos::memory_usage()
{
  double bytes = 0.0;
  bytes += maxatom * sizeof(double);
  bytes += 3*maxatom * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
    setup of ScaFaCoS handle with common parameters 
------------------------------------------------------------------------- */

void Scafacos::setup_handle()
{
  // store simulation box params
  // NOTE: this assumes orthogonal box
  // NOTE: total particles may not be a 4-byte int
  // NOTE: what does SCFCS mean by offset?
  //       it's an integer flag in LAMMPS

  old_periodicity[0] = domain->xperiodic;
  old_periodicity[1] = domain->yperiodic;
  old_periodicity[2] = domain->zperiodic;

  old_offset[0] = domain->boundary[0][0];
  old_offset[1] = domain->boundary[1][0];
  old_offset[2] = domain->boundary[2][0];

  old_box_x[0] = domain->prd[0];
  old_box_x[1] = old_box_x[2] = 0.0;
  old_box_y[1] = domain->prd[1];
  old_box_y[0] = old_box_y[2] = 0.0;
  old_box_z[2] = domain->prd[2];
  old_box_z[1] = old_box_z[0] = 0.0;

  old_natoms = atom->natoms;

  // set all required ScaFaCoS params

  result = fcs_set_box_a(fcs,old_box_x);
  check_result(result);

  result = fcs_set_box_b(fcs,old_box_y);
  check_result(result);

  result = fcs_set_box_c(fcs,old_box_z);
  check_result(result);

  result = fcs_set_box_origin(fcs,old_offset);
  check_result(result);

  result = fcs_set_periodicity(fcs,old_periodicity);
  check_result(result);

  result = fcs_set_total_particles(fcs,old_natoms);
  check_result(result);

  // NOTE: for now disable short-range calculations within LAMMPS

  int near_field_flag = 0;
  result = fcs_set_near_field_flag(fcs,near_field_flag);
  check_result(result);
}

/* ----------------------------------------------------------------------
    check if box parameters changed, requiring a new call to fcs_tune
------------------------------------------------------------------------- */

bool Scafacos::box_has_changed()
{
  int *periodicity = domain->periodicity;
  double *prd = domain->prd;

  bool changed = 
    (periodicity[0] != old_periodicity[0]) ||
    (periodicity[1] != old_periodicity[1]) ||
    (periodicity[2] != old_periodicity[2]) ||
    (domain->boundary[0][0] != old_offset[0]) ||
    (domain->boundary[1][0] != old_offset[1]) ||
    (domain->boundary[2][0] != old_offset[2]) ||
    (prd[0] != old_box_x[0]) ||
    (prd[1] != old_box_y[1]) ||
    (prd[2] != old_box_z[2]) ||
    (atom->natoms != old_natoms);

  return changed;
}

/* ----------------------------------------------------------------------
   check ScaFaCoS result for error condition
   NOTE: will all procs have same error?
------------------------------------------------------------------------- */

void Scafacos::check_result(FCSResult result) 
{
  if (!result) return;

  std::stringstream ss;
  ss << "ScaFaCoS: " << fcs_result_get_function(result) << "\n" 
     << fcs_result_get_message(result) << "\n";
  fcs_result_destroy(result);
  std::string err_msg = ss.str();
  const char *str = err_msg.c_str();

  error->one(FLERR,str);
}
 
