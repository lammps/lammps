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
   Contributing author: Carsten Svaneborg
   (Max Planck Institute for Complex Systems, Dresden, Germany)
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "fix_uniaxial.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "comm.h"
#include "kspace.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixUniaxial::FixUniaxial(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all("Illegal fix uniaxial command");
  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all("Illegal fix uniaxial command");

  if (strcmp(arg[4],"x") == 0) dir = 0;
  else if (strcmp(arg[4],"y") == 0) dir = 1;
  else if (strcmp(arg[4],"z") == 0) dir = 2;
  else error->all("Illegal fix uniaxial command");
  lambda_final = atof(arg[5]);

  if (lambda_final <= 0) error->all("Illegal fix uniaxial command");
  if (domain->nonperiodic) 
    error->all("Cannot fix uniaxial on non-periodic system");

  nrigid = 0;
  rfix = NULL;
}

/* ---------------------------------------------------------------------- */

FixUniaxial::~FixUniaxial()
{
  delete [] rfix;
}

/* ---------------------------------------------------------------------- */

int FixUniaxial::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixUniaxial::init()
{
  // store pointers to domain variable so can loop over dimensions

  domainlo[0] = &domain->boxxlo;
  domainlo[1] = &domain->boxylo;
  domainlo[2] = &domain->boxzlo;
  
  domainhi[0] = &domain->boxxhi;
  domainhi[1] = &domain->boxyhi;
  domainhi[2] = &domain->boxzhi;

  domainprd[0] = &domain->xprd;
  domainprd[1] = &domain->yprd;
  domainprd[2] = &domain->zprd;

  double L = pow((domain->boxxhi-domain->boxxlo)*
                 (domain->boxyhi-domain->boxylo)*
                 (domain->boxzhi-domain->boxzlo) ,1.0/3.0);
   
  // save box sizes for coordinate rescaling
  // calculate strains and asymmetry parameter
  // alpha=lampdai[first]/lampbdai[second] for the two perp directions

  alpha0 = 1.0;
  for (int m = 0; m < 3; m++) {
    domainloi[m] = *domainlo[m];
    domainhii[m] = *domainhi[m];
    lambdai[m] = (*domainhi[m] - *domainlo[m])/L;    
    lambdaf[m] = ( m==dir ? lambda_final : 1.0/sqrt(lambda_final) ) ;    
    if (m != dir) {
      alpha0*= lambdai[m];
      alpha0=1.0/alpha0;
    }
  } 
  
  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Initial strain = %g %g %g\n",
	      lambdai[0],lambdai[1],lambdai[2]);
      fprintf(screen,"Target strain = %g %g %g\n",
	      lambdaf[0],lambdaf[1],lambdaf[2]);
    }
    if (logfile) {
      fprintf(logfile,"Initial strain = %g %g %g\n",
	      lambdai[0],lambdai[1],lambdai[2]);
      fprintf(logfile,"Target strain = %g %g %g\n",
	      lambdaf[0],lambdaf[1],lambdaf[2]);
    }
  }

  if (force->kspace) kspace_flag = 1;
  else kspace_flag = 0;

  // detect if any fix rigid exist so rigid bodies can be re-scaled
  // rfix[] = indices to each fix rigid

  delete [] rfix;
  nrigid = 0;
  rfix = NULL;

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"rigid") == 0 ||
	strcmp(modify->fix[i]->style,"poems") == 0) nrigid++;
  if (nrigid) {
    rfix = new int[nrigid];
    nrigid = 0;
    for (int i = 0; i < modify->nfix; i++)
      if (strcmp(modify->fix[i]->style,"rigid") == 0 ||
	  strcmp(modify->fix[i]->style,"poems") == 0) rfix[nrigid++] = i;
  }
}

/* ---------------------------------------------------------------------- */

void FixUniaxial::end_of_step()
{
  int i,m;
  double oldlo,oldhi,newlo,newhi,ratio;

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
 
  double lvalue[3];

  // linear interpolation of strain in specified direction
  
  lvalue[dir] = lambdai[dir]*(1.0-delta) + lambdaf[dir]*delta;

  // linear interpolation of asymmetry parameter in the perp direction

  double alpha = alpha0*(1-delta) + delta;

  // calculate strains perpendicular to dir

  for (m = 0; m < 3; m++)
    if (m != dir) {
      lvalue[m] = sqrt(alpha/lvalue[dir]);
      alpha=1.0/alpha;
    }

  // apply appropriate rescaling in each dimension

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (m = 0; m < 3; m++) {
    oldlo = *domainlo[m];
    oldhi = *domainhi[m];
    
    newlo = domainloi[m] * lvalue[m]/lambdai[m];
    newhi = domainhii[m] * lvalue[m]/lambdai[m];
    ratio = (newhi - newlo) / *domainprd[m];

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	x[i][m] = newlo + (x[i][m] - oldlo) * ratio;
    
    *domainlo[m] = newlo;
    *domainhi[m] = newhi;
    *domainprd[m] = newhi - newlo;
  
    if (nrigid)
      for (i = 0; i < nrigid; i++)
	modify->fix[rfix[i]]->dilate(m,oldlo,oldhi,newlo,newhi);
  }

  // redo KSpace coeffs since volume has changed

  if (kspace_flag) force->kspace->setup();
}

