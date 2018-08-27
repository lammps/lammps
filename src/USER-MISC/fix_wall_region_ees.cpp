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
   Contributing author:  Abdoreza Ershadinia, a.ershadinia at gmail.com
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "fix_wall_region_ees.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_ellipsoid.h"
#include "domain.h"
#include "region.h"
#include "force.h"
#include "lattice.h"
#include "update.h"
#include "output.h"
#include "respa.h"
#include "error.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixWallRegionEES::FixWallRegionEES(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 7) error->all(FLERR,"Illegal fix wall/region/ees command");
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  // parse args

  iregion = domain->find_region(arg[3]);
  if (iregion == -1)
    error->all(FLERR,"Region ID for fix wall/region/ees does not exist");
  int n = strlen(arg[3]) + 1;
  idregion = new char[n];
  strcpy(idregion,arg[3]);

  epsilon = force->numeric(FLERR,arg[4]);
  sigma = force->numeric(FLERR,arg[5]);
  cutoff = force->numeric(FLERR,arg[6]);

  if (cutoff <= 0.0) error->all(FLERR,"Fix wall/region/ees cutoff <= 0.0");

  eflag = 0;
  ewall[0] = ewall[1] = ewall[2] = ewall[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixWallRegionEES::~FixWallRegionEES()
{
  delete [] idregion;
}

/* ---------------------------------------------------------------------- */

int FixWallRegionEES::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallRegionEES::init()
{
  // set index and check validity of region

  iregion = domain->find_region(idregion);
  if (iregion == -1)
    error->all(FLERR,"Region ID for fix wall/region/ees does not exist");

  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec)
    error->all(FLERR,"Fix wall/region/ees requires atom style ellipsoid");

  // check that all particles are finite-size ellipsoids
  // no point particles allowed, spherical is OK

  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (ellipsoid[i] < 0)
        error->one(FLERR,"Fix wall/region/ees requires extended particles");

  // setup coefficients

  coeff1 = ( 2. / 4725. ) * epsilon * pow(sigma,12.0);
  coeff2 = ( 1. / 24. ) * epsilon * pow(sigma,6.0);
  coeff3 = ( 2. / 315. ) * epsilon * pow(sigma,12.0);
  coeff4 = ( 1. / 3. ) * epsilon * pow(sigma,6.0);
  coeff5 = ( 4. / 315. ) * epsilon * pow(sigma,12.0);
  coeff6 = ( 1. / 12. )  * epsilon * pow(sigma,6.0);
  offset = 0;


  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixWallRegionEES::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallRegionEES::min_setup(int vflag)
{
    post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallRegionEES::post_force(int /*vflag*/)
{
  //sth is needed here, but I dont know what
  //that is calculation of sn

  int i,m,n;
  double rinv,fx,fy,fz,sn,tooclose[3];

  eflag = 0;
  ewall[0] = ewall[1] = ewall[2] = ewall[3] = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double **tor = atom->torque;

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  Region *region = domain->regions[iregion];
  region->prematch();

  int onflag = 0;

  // region->match() insures particle is in region or on surface, else error
  // if returned contact dist r = 0, is on surface, also an error
  // in COLLOID case, r <= radius is an error

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (!region->match(x[i][0],x[i][1],x[i][2])) {
        onflag = 1;
        continue;
      }

      double A[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
      double tempvec[3]= {0,0,0};
      double sn2 = 0.0;
      double nhat[3] = {0,0,0};
      double* shape = bonus[ellipsoid[i]].shape;;
      MathExtra::quat_to_mat(bonus[ellipsoid[i]].quat,A);

      for(int which = 0 ; which < 3; which ++){//me
        nhat[which]=1;
        nhat[(which+1)%3] = 0 ;
        nhat[(which+2)%3] = 0 ;
        sn2 = 0 ;
        MathExtra::transpose_matvec(A,nhat,tempvec);
        for(int k = 0; k<3; k++) tempvec[k] *= shape[k];
        for(int k = 0; k<3 ; k++) sn2 += tempvec[k]*tempvec[k];
        sn = sqrt(sn2);
        tooclose[which] = sn;
      }

      n = region->surface(x[i][0],x[i][1],x[i][2],cutoff);

      for (m = 0; m < n; m++) {

        if (region->contact[m].delx != 0 && region->contact[m].r <= tooclose[0]){
          onflag = 1;
          continue;
        } else if (region->contact[m].dely != 0 && region->contact[m].r <= tooclose[1]){
          onflag = 1;
          continue;
        } else if (region->contact[m].delz !=0 && region->contact[m].r <= tooclose[2]){
          onflag = 1;
          continue;
        } else rinv = 1.0/region->contact[m].r;

        ees(m,i);

        ewall[0] += eng;
        fx = fwall * region->contact[m].delx * rinv;
        fy = fwall * region->contact[m].dely * rinv;
        fz = fwall * region->contact[m].delz * rinv;
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        ewall[1] -= fx;
        ewall[2] -= fy;
        ewall[3] -= fz;

        tor[i][0] += torque[0];
        tor[i][1] += torque[1];
        tor[i][2] += torque[2];
      }
    }

  if (onflag) error->one(FLERR,"Particle on or inside surface of region "
                         "used in fix wall/region/ees");
}

/* ---------------------------------------------------------------------- */

void FixWallRegionEES::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallRegionEES::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of wall interaction
------------------------------------------------------------------------- */

double FixWallRegionEES::compute_scalar()
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewall,ewall_all,4,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return ewall_all[0];
}

/* ----------------------------------------------------------------------
   components of force on wall
------------------------------------------------------------------------- */

double FixWallRegionEES::compute_vector(int n)
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewall,ewall_all,4,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return ewall_all[n+1];
}

/* ----------------------------------------------------------------------
   EES interaction for ellipsoid particle with wall
   compute eng and fwall and twall = magnitude of wall force and torque
------------------------------------------------------------------------- */

void FixWallRegionEES::ees(int m, int i)
{
  Region *region = domain->regions[iregion];
  region->prematch();

  double delta, delta2, delta3, delta4, delta5, delta6;
  double sigman, sigman2 , sigman3, sigman4, sigman5, sigman6;
  double hhss, hhss2, hhss4, hhss7, hhss8; //h^2 - s_n^2
  double hps;                              //h+s_n
  double hms;                              //h-s_n
  double twall;

  double A[3][3], nhat[3], SAn[3], that[3];

  double tempvec[3]= {0,0,0};
  double tempvec2[3]= {0,0,0};

  double Lx[3][3] = {{0,0,0},{0,0,-1},{0,1,0}};
  double Ly[3][3] = {{0,0,1},{0,0,0},{-1,0,0}};
  double Lz[3][3] = {{0,-1,0},{1,0,0},{0,0,0}};

  nhat[0] = region->contact[m].delx / region->contact[m].r;
  nhat[1] = region->contact[m].dely / region->contact[m].r;
  nhat[2] = region->contact[m].delz / region->contact[m].r;

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;

  double* shape = bonus[ellipsoid[i]].shape;;
  MathExtra::quat_to_mat(bonus[ellipsoid[i]].quat,A);

  sigman2 = 0.0;
  MathExtra::transpose_matvec(A,nhat,tempvec);
  for(int k = 0; k<3; k++) tempvec[k] *= shape[k];
  for(int k = 0; k<3; k++) sigman2 += tempvec[k]*tempvec[k];
  for(int k = 0; k<3; k++) SAn[k] = tempvec[k];

  sigman = sqrt(sigman2);
  delta = fabs(region->contact[m].r);

  sigman3 = sigman2 * sigman;
  sigman4 = sigman2 * sigman2;
  sigman5 = sigman4 * sigman;
  sigman6 = sigman3 * sigman3;

  delta2 = delta  * delta;
  delta3 = delta2 * delta;
  delta4 = delta2 * delta2;
  delta5 = delta3 * delta2;
  delta6 = delta3 * delta3;

  hhss = delta2 - sigman2;
  hhss2 = hhss  * hhss;
  hhss4 = hhss2 * hhss2;
  hhss8 = hhss4 * hhss4;
  hhss7 = hhss4 * hhss2 * hhss;

  hps = delta + sigman;
  hms = delta - sigman;

  fwall =  -1*coeff4/hhss2 + coeff3
    * (21*delta6 + 63*delta4*sigman2 + 27*delta2*sigman4 + sigman6) / hhss8;

  eng = -1*coeff2 * (4*delta/sigman2/hhss + 2*log(hms/hps)/sigman3) +
    coeff1 * (35*delta5 + 70*delta3*sigman2 + 15*delta*sigman4) / hhss7;

  twall = coeff6 * (6*delta3/sigman4/hhss2 - 10*delta/sigman2/hhss2
                    + 3*log(hms/hps)/sigman5)
    + coeff5 * (21.*delta5 + 30.*delta3*sigman2 + 5.*delta*sigman4) / hhss8;

  MathExtra::matvec(Lx,nhat,tempvec);
  MathExtra::transpose_matvec(A,tempvec,tempvec2);
  for(int k = 0; k<3; k++) tempvec2[k] *= shape[k];
  that[0] = MathExtra::dot3(SAn,tempvec2);

  MathExtra::matvec(Ly,nhat,tempvec);
  MathExtra::transpose_matvec(A,tempvec,tempvec2);
  for(int k = 0; k<3; k++) tempvec2[k] *= shape[k];
  that[1] = MathExtra::dot3(SAn,tempvec2);

  MathExtra::matvec(Lz,nhat,tempvec);
  MathExtra::transpose_matvec(A,tempvec,tempvec2);
  for(int k = 0; k < 3; k++) tempvec2[k] *= shape[k];
  that[2] = MathExtra::dot3(SAn,tempvec2);

  for(int j = 0; j<3 ; j++)
    torque[j] = twall * that[j];
}
