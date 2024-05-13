// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
                    and Wengen Ouyang (Wuhan U)
------------------------------------------------------------------------- */

#include "pair_sw_mod_omp.h"

#include "error.h"
#include "math_const.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

PairSWMODOMP::PairSWMODOMP(LAMMPS *lmp) : PairSWOMP(lmp)
{
  delta1 = 0.25;
  delta2 = 0.35;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSWMODOMP::settings(int narg, char **arg)
{
  // process optional keywords

  int iarg = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"maxdelcs") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal pair_style command");
      delta1 = utils::numeric(FLERR,arg[iarg+1],false,Pointers::lmp);
      delta2 = utils::numeric(FLERR,arg[iarg+2],false,Pointers::lmp);
      iarg += 3;
      if ((delta1 < 0.0) || (delta1 > 1.0) || (delta2 < 0.0) || (delta2 > 1.0) || (delta1 > delta2))
        error->all(FLERR,"Illegal values for maxdelcs keyword");
    } else error->all(FLERR,"Illegal pair_style command");
  }
  PairSWOMP::settings(narg-iarg,arg+iarg);
}

/* ---------------------------------------------------------------------- */

void PairSWMODOMP::threebody(Param *paramij, Param *paramik, Param *paramijk,
                            double rsq1, double rsq2,
                            double *delr1, double *delr2,
                            double *fj, double *fk, int eflag, double &eng)
{
  double r1,rinvsq1,rainv1,gsrainv1,gsrainvsq1,expgsrainv1;
  double r2,rinvsq2,rainv2,gsrainv2,gsrainvsq2,expgsrainv2;
  double rinv12,cs,delcs,delcssq,facexp,facrad,frad1,frad2;
  double facang,facang12,csfacang,csfac1,csfac2,factor;

  r1 = sqrt(rsq1);
  rinvsq1 = 1.0/rsq1;
  rainv1 = 1.0/(r1 - paramij->cut);
  gsrainv1 = paramij->sigma_gamma * rainv1;
  gsrainvsq1 = gsrainv1*rainv1/r1;
  expgsrainv1 = exp(gsrainv1);

  r2 = sqrt(rsq2);
  rinvsq2 = 1.0/rsq2;
  rainv2 = 1.0/(r2 - paramik->cut);
  gsrainv2 = paramik->sigma_gamma * rainv2;
  gsrainvsq2 = gsrainv2*rainv2/r2;
  expgsrainv2 = exp(gsrainv2);

  rinv12 = 1.0/(r1*r2);
  cs = (delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2]) * rinv12;
  delcs = cs - paramijk->costheta;

  // Modification to delcs
  if(fabs(delcs) >= delta2) delcs = 0.0;
  else if(fabs(delcs) < delta2 && fabs(delcs) > delta1) {
    factor = 0.5 + 0.5*cos(MY_PI*(fabs(delcs) - delta1)/(delta2 - delta1));
    delcs *= factor;
  }
  delcssq = delcs*delcs;

  facexp = expgsrainv1*expgsrainv2;

  // facrad = sqrt(paramij->lambda_epsilon*paramik->lambda_epsilon) *
  //          facexp*delcssq;

  facrad = paramijk->lambda_epsilon * facexp*delcssq;
  frad1 = facrad*gsrainvsq1;
  frad2 = facrad*gsrainvsq2;
  facang = paramijk->lambda_epsilon2 * facexp*delcs;
  facang12 = rinv12*facang;
  csfacang = cs*facang;
  csfac1 = rinvsq1*csfacang;

  fj[0] = delr1[0]*(frad1+csfac1)-delr2[0]*facang12;
  fj[1] = delr1[1]*(frad1+csfac1)-delr2[1]*facang12;
  fj[2] = delr1[2]*(frad1+csfac1)-delr2[2]*facang12;

  csfac2 = rinvsq2*csfacang;

  fk[0] = delr2[0]*(frad2+csfac2)-delr1[0]*facang12;
  fk[1] = delr2[1]*(frad2+csfac2)-delr1[1]*facang12;
  fk[2] = delr2[2]*(frad2+csfac2)-delr1[2]*facang12;

  if (eflag) eng = facrad;
}
