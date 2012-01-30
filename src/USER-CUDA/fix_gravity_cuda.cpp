/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator 

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov 

   See the README file in the top-level LAMMPS directory. 

   ----------------------------------------------------------------------- 

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/ 

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany 

   See the README file in the USER-CUDA directory. 

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "fix_gravity_cuda.h"
#include "fix_gravity_cuda_cu.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "cuda.h"
#include "cuda_modify_flags.h"
#include "math_const.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace FixConstCuda;
using namespace MathConst;

enum{CHUTE,SPHERICAL,GRADIENT,VECTOR};

/* ---------------------------------------------------------------------- */

FixGravityCuda::FixGravityCuda(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

  if (narg < 5) error->all(FLERR,"Illegal fix gravity command");

  time_depend = 1;

  magnitude = atof(arg[3]);

  if (strcmp(arg[4],"chute") == 0) {
    if (narg != 6) error->all(FLERR,"Illegal fix gravity command");
    style = CHUTE;
    phi = 0.0;
    theta = 180.0 - atof(arg[5]);
  } else if (strcmp(arg[4],"spherical") == 0) {
    if (narg != 7) error->all(FLERR,"Illegal fix gravity command");
    style = SPHERICAL;
    phi = atof(arg[5]);
    theta = atof(arg[6]);
  } else if (strcmp(arg[4],"gradient") == 0) {
    if (narg != 9) error->all(FLERR,"Illegal fix gravity command");
    style = GRADIENT;
    phi = atof(arg[5]);
    theta = atof(arg[6]);
    phigrad = atof(arg[7]);
    thetagrad = atof(arg[8]);
  } else if (strcmp(arg[4],"vector") == 0) {
    if (narg != 8) error->all(FLERR,"Illegal fix gravity command");
    style = VECTOR;
    xdir = atof(arg[5]);
    ydir = atof(arg[6]);
    zdir = atof(arg[7]);
  } else error->all(FLERR,"Illegal fix gravity command");

  degree2rad = MY_PI/180.0;

  if (style == CHUTE || style == SPHERICAL || style == GRADIENT) {
    if (domain->dimension == 3) {
      xgrav = sin(degree2rad * theta) * cos(degree2rad * phi);
      ygrav = sin(degree2rad * theta) * sin(degree2rad * phi);
      zgrav = cos(degree2rad * theta);
    } else {
      xgrav = sin(degree2rad * theta);
      ygrav = cos(degree2rad * theta);
      zgrav = 0.0;
    }
  } else if (style == VECTOR) {
    if (domain->dimension == 3) {
      double length = sqrt(xdir*xdir + ydir*ydir + zdir*zdir);
      xgrav = xdir/length;
      ygrav = ydir/length;
      zgrav = zdir/length;
    } else {
      double length = sqrt(xdir*xdir + ydir*ydir);
      xgrav = xdir/length;
      ygrav = ydir/length;
      zgrav = 0.0;
    }
  }

  time_origin = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

int FixGravityCuda::setmask()
{
  int mask = 0;
  mask |= POST_FORCE_CUDA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGravityCuda::init()
{
  dt = update->dt;

  xacc = magnitude*xgrav;
  yacc = magnitude*ygrav;
  zacc = magnitude*zgrav;
}

/* ---------------------------------------------------------------------- */

void FixGravityCuda::setup(int vflag)
{
  MYDBG( printf("# CUDA: FixGravityCuda::setup\n"); )
	
  if (strstr(update->integrate_style,"verlet"))
  {
    Cuda_FixGravityCuda_Init(&cuda->shared_data);
    cuda->cu_f->upload();
    post_force(vflag);
    cuda->cu_f->download();
    
  }
  else {
  }
  MYDBG( printf("# CUDA: FixGravityCuda::setup done\n"); )
}

/* ---------------------------------------------------------------------- */

void FixGravityCuda::post_force(int vflag)
{
  // update direction of gravity vector if gradient style

  if (style == GRADIENT) {
    if (domain->dimension == 3) {
      double phi_current = degree2rad * 
	(phi + (update->ntimestep - time_origin)*dt*phigrad*360.0);
      double theta_current = degree2rad * 
	(theta + (update->ntimestep - time_origin)*dt*thetagrad*360.0);
      xgrav = sin(theta_current) * cos(phi_current);
      ygrav = sin(theta_current) * sin(phi_current);
      zgrav = cos(theta_current);
    } else {
      double theta_current = degree2rad * 
	(theta + (update->ntimestep - time_origin)*dt*thetagrad*360.0);
      xgrav = sin(theta_current);
      ygrav = cos(theta_current);
    }
    xacc = magnitude*xgrav;
    yacc = magnitude*ygrav;
    zacc = magnitude*zgrav;
  }

  MYDBG( printf("# CUDA: FixGravityCuda::postforce start\n"); )
  Cuda_FixGravityCuda_PostForce(&cuda->shared_data, groupbit, xacc,yacc,zacc);
}


