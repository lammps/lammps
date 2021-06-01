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
   Contributing author:  Abdoreza Ershadinia, a.ershadinia at gmail.com
------------------------------------------------------------------------- */

#include "fix_wall_ees.h"

#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "error.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixWallEES::FixWallEES(LAMMPS *lmp, int narg, char **arg) :
  FixWall(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

void FixWallEES::precompute(int m)
{
  coeff1[m] = ( 2. / 4725. ) * epsilon[m] * pow(sigma[m],12.0);
  coeff2[m] = ( 1. / 24. ) * epsilon[m] * pow(sigma[m],6.0);

  coeff3[m] = ( 2. / 315. ) * epsilon[m] * pow(sigma[m],12.0);
  coeff4[m] = ( 1. / 3. ) * epsilon[m] * pow(sigma[m],6.0);

  coeff5[m] = ( 4. / 315. ) * epsilon[m] * pow(sigma[m],12.0);
  coeff6[m] = ( 1. / 12. )  * epsilon[m] * pow(sigma[m],6.0);
}

/* ---------------------------------------------------------------------- */
void FixWallEES::init()
{
  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec)
    error->all(FLERR,"Fix wall/ees requires atom style ellipsoid");

  // check that all particles are finite-size ellipsoids
  // no point particles allowed, spherical is OK

  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (ellipsoid[i] < 0)
        error->one(FLERR,"Fix wall/ees requires extended particles");

  FixWall::init();
}


/* ----------------------------------------------------------------------
   interaction of all particles in group with a wall
   m = index of wall coeffs
   which = xlo,xhi,ylo,yhi,zlo,zhi
   error if any particle is on or behind wall
------------------------------------------------------------------------- */

void FixWallEES::wall_particle(int m, int which, double coord)
{
  double delta;

  double **x = atom->x;
  double **f = atom->f;
  double **tor = atom->torque;

  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int dim = which / 2;
  int side = which % 2;
  if (side == 0) side = -1;

  int onflag = 0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      if (side < 0)
        delta = x[i][dim] - coord;
      else
        delta = coord - x[i][dim];

      if (delta >= cutoff[m])
        continue;

      double A[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
      double tempvec[3]= {0,0,0};
      double sigman = 0.0, sigman2 = 0.0;
      double nhat[3] = {0,0,0};

      nhat[dim]=-1*side;
      nhat[(dim+1)%3] = 0 ;
      nhat[(dim+2)%3] = 0 ;


      double* shape = bonus[ellipsoid[i]].shape;;
      MathExtra::quat_to_mat(bonus[ellipsoid[i]].quat,A);
      MathExtra::transpose_matvec(A,nhat,tempvec);
      for (int k = 0; k<3; k++) tempvec[k] *= shape[k];
      for (int k = 0; k<3 ; k++) sigman2 += tempvec[k]*tempvec[k];
      sigman = sqrt(sigman2);

      if (delta <= sigman) {
        onflag = 1;
        continue;
      }


      double fwall = 0.0, twall = 0.0;
      double delta2 = 0.0, delta3 = 0.0, delta4 = 0.0, delta5 = 0.0, delta6 = 0.0;
      double sigman3 = 0.0, sigman4 = 0.0, sigman5 = 0.0, sigman6 = 0.0;
      double hhss = 0.0, hhss2 = 0.0, hhss4 = 0.0, hhss7 = 0.0, hhss8 = 0.0;
      double hps = 0.0;
      double hms = 0.0;

      double tempvec2[3]= {0,0,0};

      double SAn[3] = {0,0,0};
      double that[3] = {0,0,0};

      double Lx[3][3] =  {{0,0,0},{0,0,-1},{0,1,0}};
      double Ly[3][3] =  {{0,0,1},{0,0,0},{-1,0,0}};
      double Lz[3][3] =  {{0,-1,0},{1,0,0},{0,0,0}};


      for (int k = 0; k<3; k++) SAn[k] = tempvec[k];

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

      fwall = side*(
        -1*coeff4[m]/hhss2 +
        coeff3[m] * (  21*delta6 + 63*delta4*sigman2 + 27*delta2*sigman4 + sigman6  ) / hhss8
        );
      f[i][dim] -= fwall;

      ewall[0] += -1*coeff2[m] * (  4*delta/sigman2/hhss + 2*log(hms/hps)/sigman3  ) +
        coeff1[m] * (  35*delta5 + 70*delta3*sigman2 + 15*delta*sigman4  ) / hhss7;

      ewall[m+1] += fwall;

      twall = coeff6[m] * (  6.*delta3/sigman4/hhss2  - 10.*delta/sigman2/hhss2 + 3.*log(hms/hps)/sigman5  )  +
        coeff5[m] * (  21.*delta5 + 30.*delta3*sigman2 + 5.*delta*sigman4  ) / hhss8    ;


      MathExtra::matvec(Lx,nhat,tempvec);
      MathExtra::transpose_matvec(A,tempvec,tempvec2);
      for (int k = 0; k<3; k++) tempvec2[k] *= shape[k];
      that[0] = MathExtra::dot3(SAn,tempvec2);

      MathExtra::matvec(Ly,nhat,tempvec);
      MathExtra::transpose_matvec(A,tempvec,tempvec2);
      for (int k = 0; k<3; k++) tempvec2[k] *= shape[k];
      that[1] = MathExtra::dot3(SAn,tempvec2);

      MathExtra::matvec(Lz,nhat,tempvec);
      MathExtra::transpose_matvec(A,tempvec,tempvec2);
      for (int k = 0; k < 3; k++) tempvec2[k] *= shape[k];
      that[2] = MathExtra::dot3(SAn,tempvec2);


      for (int j = 0; j<3 ; j++)
        tor[i][j] += twall * that[j];

    }

  if (onflag) error->one(FLERR,"Particle on or inside fix wall surface");
}
