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
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "compute_temp_asphere.h"

#include <cstring>
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{ROTATE,ALL};

static constexpr double INERTIA = 0.2;          // moment of inertia prefactor for ellipsoid

/* ---------------------------------------------------------------------- */

ComputeTempAsphere::ComputeTempAsphere(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  id_bias(nullptr), tbias(nullptr), avec(nullptr)
{
  if (narg < 3) error->all(FLERR,"Illegal compute temp/asphere command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;

  tempbias = 0;
  id_bias = nullptr;
  mode = ALL;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"bias") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute temp/asphere command");
      tempbias = 1;
      id_bias = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dof") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute temp/asphere command");
      if (strcmp(arg[iarg+1],"rotate") == 0) mode = ROTATE;
      else if (strcmp(arg[iarg+1],"all") == 0) mode = ALL;
      else error->all(FLERR,"Illegal compute temp/asphere command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute temp/asphere command");
  }

  // when computing only the rotational temperature,
  // do not remove DOFs for translation as set by default

  if (mode == ROTATE) extra_dof = 0;

  vector = new double[size_vector];

}

/* ---------------------------------------------------------------------- */

ComputeTempAsphere::~ComputeTempAsphere()
{
  delete [] id_bias;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempAsphere::init()
{
  // error check

  avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  if (!avec)
    error->all(FLERR,"Compute temp/asphere requires atom style ellipsoid");

  // check that all particles are finite-size, no point particles allowed

  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (ellipsoid[i] < 0)
        error->one(FLERR,"Compute temp/asphere requires extended particles");

  if (tempbias) {
    int i = modify->find_compute(id_bias);
    if (i < 0)
      error->all(FLERR,"Could not find compute ID for temperature bias");
    tbias = modify->compute[i];
    if (tbias->tempflag == 0)
      error->all(FLERR,"Bias compute does not calculate temperature");
    if (tbias->tempbias == 0)
      error->all(FLERR,"Bias compute does not calculate a velocity bias");
    if (tbias->igroup != igroup)
      error->all(FLERR,"Bias compute group does not match compute group");
    if (strcmp(tbias->style,"temp/region") == 0) tempbias = 2;
    else tempbias = 1;

    // init and setup bias compute because
    // this compute's setup()->dof_compute() may be called first

    tbias->init();
    tbias->setup();
  }
}

/* ---------------------------------------------------------------------- */

void ComputeTempAsphere::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTempAsphere::dof_compute()
{
  adjust_dof_fix();

  // 6 dof for 3d, 3 dof for 2d
  // which dof are included also depends on mode
  // assume full rotation of extended particles
  // user should correct this via compute_modify if needed

  natoms_temp = group->count(igroup);
  int nper;
  if (domain->dimension == 3) {
    if (mode == ALL) nper = 6;
    else nper = 3;
  } else {
    if (mode == ALL) nper = 3;
    else nper = 1;
  }
  dof = nper*natoms_temp;

  // additional adjustments to dof

  if (tempbias == 1) {
    if (mode == ALL) dof -= tbias->dof_remove(-1) * natoms_temp;

  } else if (tempbias == 2) {
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    tbias->dof_remove_pre();

    int count = 0;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        if (tbias->dof_remove(i)) count++;
    int count_all;
    MPI_Allreduce(&count,&count_all,1,MPI_INT,MPI_SUM,world);
    dof -= nper*count_all;
  }

  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempAsphere::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  if (tempbias) {
    if (tbias->invoked_scalar != update->ntimestep) tbias->compute_scalar();
    tbias->remove_bias_all();
  }

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  double **v = atom->v;
  double **angmom = atom->angmom;
  double *rmass = atom->rmass;
  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *shape,*quat;
  double wbody[3],inertia[3];
  double rot[3][3];

  // sum translational and rotational energy for each particle
  // no point particles since divide by inertia

  double t = 0.0;

  if (mode == ALL) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * rmass[i];

        // principal moments of inertia

        shape = bonus[ellipsoid[i]].shape;
        quat = bonus[ellipsoid[i]].quat;

        inertia[0] = INERTIA*rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]);
        inertia[1] = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]);
        inertia[2] = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]);

        // wbody = angular velocity in body frame

        MathExtra::quat_to_mat(quat,rot);
        MathExtra::transpose_matvec(rot,angmom[i],wbody);
        wbody[0] /= inertia[0];
        wbody[1] /= inertia[1];
        wbody[2] /= inertia[2];

        t += inertia[0]*wbody[0]*wbody[0] +
          inertia[1]*wbody[1]*wbody[1] + inertia[2]*wbody[2]*wbody[2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {

        // principal moments of inertia

        shape = bonus[ellipsoid[i]].shape;
        quat = bonus[ellipsoid[i]].quat;

        inertia[0] = INERTIA*rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]);
        inertia[1] = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]);
        inertia[2] = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]);

        // wbody = angular velocity in body frame

        MathExtra::quat_to_mat(quat,rot);
        MathExtra::transpose_matvec(rot,angmom[i],wbody);
        wbody[0] /= inertia[0];
        wbody[1] /= inertia[1];
        wbody[2] /= inertia[2];

        t += inertia[0]*wbody[0]*wbody[0] +
          inertia[1]*wbody[1]*wbody[1] + inertia[2]*wbody[2]*wbody[2];
      }
  }

  if (tempbias) tbias->restore_bias_all();

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic || tempbias == 2) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempAsphere::compute_vector()
{
  int i;

  invoked_vector = update->ntimestep;

  if (tempbias) {
    if (tbias->invoked_vector != update->ntimestep) tbias->compute_vector();
    tbias->remove_bias_all();
  }

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  double **v = atom->v;
  double **angmom = atom->angmom;
  double *rmass = atom->rmass;
  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *shape,*quat;
  double wbody[3],inertia[3],t[6];
  double rot[3][3];
  double massone;

  // sum translational and rotational energy for each particle
  // no point particles since divide by inertia

  for (i = 0; i < 6; i++) t[i] = 0.0;

  if (mode == ALL) {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = rmass[i];
        t[0] += massone * v[i][0]*v[i][0];
        t[1] += massone * v[i][1]*v[i][1];
        t[2] += massone * v[i][2]*v[i][2];
        t[3] += massone * v[i][0]*v[i][1];
        t[4] += massone * v[i][0]*v[i][2];
        t[5] += massone * v[i][1]*v[i][2];

        // principal moments of inertia

        shape = bonus[ellipsoid[i]].shape;
        quat = bonus[ellipsoid[i]].quat;

        inertia[0] = INERTIA*massone * (shape[1]*shape[1]+shape[2]*shape[2]);
        inertia[1] = INERTIA*massone * (shape[0]*shape[0]+shape[2]*shape[2]);
        inertia[2] = INERTIA*massone * (shape[0]*shape[0]+shape[1]*shape[1]);

        // wbody = angular velocity in body frame

        MathExtra::quat_to_mat(quat,rot);
        MathExtra::transpose_matvec(rot,angmom[i],wbody);
        wbody[0] /= inertia[0];
        wbody[1] /= inertia[1];
        wbody[2] /= inertia[2];

        // rotational kinetic energy

        t[0] += inertia[0]*wbody[0]*wbody[0];
        t[1] += inertia[1]*wbody[1]*wbody[1];
        t[2] += inertia[2]*wbody[2]*wbody[2];
        t[3] += inertia[0]*wbody[0]*wbody[1];
        t[4] += inertia[1]*wbody[0]*wbody[2];
        t[5] += inertia[2]*wbody[1]*wbody[2];
      }

  } else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {

        // principal moments of inertia

        shape = bonus[ellipsoid[i]].shape;
        quat = bonus[ellipsoid[i]].quat;
        massone = rmass[i];

        inertia[0] = INERTIA*massone * (shape[1]*shape[1]+shape[2]*shape[2]);
        inertia[1] = INERTIA*massone * (shape[0]*shape[0]+shape[2]*shape[2]);
        inertia[2] = INERTIA*massone * (shape[0]*shape[0]+shape[1]*shape[1]);

        // wbody = angular velocity in body frame

        MathExtra::quat_to_mat(quat,rot);
        MathExtra::transpose_matvec(rot,angmom[i],wbody);
        wbody[0] /= inertia[0];
        wbody[1] /= inertia[1];
        wbody[2] /= inertia[2];

        // rotational kinetic energy

        t[0] += inertia[0]*wbody[0]*wbody[0];
        t[1] += inertia[1]*wbody[1]*wbody[1];
        t[2] += inertia[2]*wbody[2]*wbody[2];
        t[3] += inertia[0]*wbody[0]*wbody[1];
        t[4] += inertia[1]*wbody[0]*wbody[2];
        t[5] += inertia[2]*wbody[1]*wbody[2];
      }
  }

  if (tempbias) tbias->restore_bias_all();

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempAsphere::remove_bias(int i, double *v)
{
  if (tbias) tbias->remove_bias(i,v);
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempAsphere::remove_bias_thr(int i, double *v, double *b)
{
  if (tbias) tbias->remove_bias_thr(i,v,b);
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempAsphere::restore_bias(int i, double *v)
{
  if (tbias) tbias->restore_bias(i,v);
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias_thr()
   assume remove_bias_thr() was previously called with the same buffer b
------------------------------------------------------------------------- */

void ComputeTempAsphere::restore_bias_thr(int i, double *v, double *b)
{
  if (tbias) tbias->restore_bias_thr(i,v,b);
}
