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

#include "fix_wall_reflect_stochastic.h"
#include <cstring>
#include <cstdlib>
#include "atom.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "force.h"
#include "lattice.h"
#include "input.h"
#include "variable.h"
#include "random_mars.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE=0,DIFFUSIVE=1,MAXWELL=2,CERCIGNANILAMPIS=3};

/* ---------------------------------------------------------------------- */

FixWallReflectStochastic::
FixWallReflectStochastic(LAMMPS *lmp, int narg, char **arg) :
  FixWallReflect(lmp, narg, arg), random(NULL)
{
  if (narg < 8) error->all(FLERR,"Illegal fix wall/reflect/stochastic command");
  
  if (domain->triclinic != 0) 
    error->all(FLERR, "Fix wall/reflect/stochastic cannot be used with "
               "triclinic simulation box");

  dynamic_group_allow = 1;

  // parse args

  int arginc;

  nwall = 0;
  int scaleflag = 1;
  reflectionstyle = NONE;
  
  if (strcmp(arg[3],"diffusive") == 0) {
    reflectionstyle = DIFFUSIVE;
    arginc = 6;
  } else if (strcmp(arg[3],"maxwell") == 0) {
    reflectionstyle = MAXWELL;
    arginc = 7;
  } else if (strcmp(arg[3],"cercignanilampis") == 0) {
    reflectionstyle = CERCIGNANILAMPIS;
    arginc = 9;
  } else error->all(FLERR,"Illegal fix wall/reflect/stochastic command");


  seedfix = force->inumeric(FLERR,arg[4]);
  if (seedfix <= 0) error->all(FLERR,"Random seed must be a postive number");

  int iarg = 5;
  while (iarg < narg) {
    if ((strcmp(arg[iarg],"xlo") == 0) || (strcmp(arg[iarg],"xhi") == 0) ||
        (strcmp(arg[iarg],"ylo") == 0) || (strcmp(arg[iarg],"yhi") == 0) ||
        (strcmp(arg[iarg],"zlo") == 0) || (strcmp(arg[iarg],"zhi") == 0)) {
      if (iarg+2 > narg) 
        error->all(FLERR,"Illegal fix wall/reflect/stochastic command");

      int newwall;
      if (strcmp(arg[iarg],"xlo") == 0) newwall = XLO;
      else if (strcmp(arg[iarg],"xhi") == 0) newwall = XHI;
      else if (strcmp(arg[iarg],"ylo") == 0) newwall = YLO;
      else if (strcmp(arg[iarg],"yhi") == 0) newwall = YHI;
      else if (strcmp(arg[iarg],"zlo") == 0) newwall = ZLO;
      else if (strcmp(arg[iarg],"zhi") == 0) newwall = ZHI;
      
      for (int m = 0; (m < nwall) && (m < 6); m++)
        if (newwall == wallwhich[m])
          error->all(FLERR,
                     "Fix wall/reflect/stochastic command wall defined twice");

      wallwhich[nwall] = newwall;
      if (strcmp(arg[iarg+1],"EDGE") == 0) {
        wallstyle[nwall] = EDGE;
        int dim = wallwhich[nwall] / 2;
        int side = wallwhich[nwall] % 2;
        if (side == 0) coord0[nwall] = domain->boxlo[dim];
        else coord0[nwall] = domain->boxhi[dim];
      } else {
        wallstyle[nwall] = CONSTANT;
        coord0[nwall] = force->numeric(FLERR,arg[iarg+1]);
      }

      walltemp[nwall]= force->numeric(FLERR,arg[iarg+2]);

      for (int dir = 0; dir < 3; dir++) {
        wallvel[nwall][dir]= force->numeric(FLERR,arg[iarg+dir+3]);
        int dim = wallwhich[nwall] / 2;
        if ((wallvel[nwall][dir] !=0) & (dir == dim)) 
          error->all(FLERR,"The wall velocity must be tangential");
        
        if (reflectionstyle == CERCIGNANILAMPIS) {
          wallaccom[nwall][dir]= force->numeric(FLERR,arg[iarg+dir+6]);
        } else if (reflectionstyle == MAXWELL) {
          // one accommodation coefficient for all directions
          wallaccom[nwall][dir]= force->numeric(FLERR,arg[iarg+6]);
        }
      }
      
      nwall++;
      iarg += arginc;
      
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) 
        error->all(FLERR,"Illegal wall/reflect/stochastic command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix wall/reflect/stochastic command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix wall/reflect/stochastic command");
  }

  // error check

  if (nwall == 0) error->all(FLERR,"Illegal fix wall command");

  for (int m = 0; m < nwall; m++) {
    if ((wallwhich[m] == XLO || wallwhich[m] == XHI) && domain->xperiodic)
      error->all(FLERR,"Cannot use fix wall/reflect/stochastic "
                 "in periodic dimension");
    if ((wallwhich[m] == YLO || wallwhich[m] == YHI) && domain->yperiodic)
      error->all(FLERR,"Cannot use fix wall/reflect/stochastic "
                 "in periodic dimension");
    if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->zperiodic)
      error->all(FLERR,"Cannot use fix wall/reflect/stochastic "
                 "in periodic dimension");
  }
  
  for (int m = 0; m < nwall; m++)
    if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->dimension == 2)
      error->all(FLERR,
                 "Cannot use fix wall/reflect/stochastic zlo/zhi "
                 "for a 2d simulation");
  
  // scale factors for CONSTANT walls, VARIABLE wall is not allowed
  
  int flag = 0;
  for (int m = 0; m < nwall; m++)
    if (wallstyle[m] != EDGE) flag = 1;
  
  if (flag) {
    if (scaleflag) {
      xscale = domain->lattice->xlattice;
      yscale = domain->lattice->ylattice;
      zscale = domain->lattice->zlattice;
    }
    else xscale = yscale = zscale = 1.0;
    
    for (int m = 0; m < nwall; m++) {
      if (wallstyle[m] != CONSTANT) continue;
      if (wallwhich[m] < YLO) coord0[m] *= xscale;
      else if (wallwhich[m] < ZLO) coord0[m] *= yscale;
      else coord0[m] *= zscale;
    }
  }

  // random number generator

  random = new RanMars(lmp,seedfix + comm->me);
}

/* ---------------------------------------------------------------------- */

FixWallReflectStochastic::~FixWallReflectStochastic()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

void FixWallReflectStochastic::wall_particle(int m, int which, double coord)
{
  int i, dir, dim, side, sign;
  double vsave,factor,timecol,difftest,theta;

  double *rmass;
  double *mass = atom->mass;

  // coord = current position of wall

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  dim = which / 2;
  side = which % 2;

  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
      
    sign = 0;
    if ((side == 0) & (x[i][dim] < coord)) sign = 1;
    else if ((side == 1) & (x[i][dim] > coord)) sign = -1;
    if (sign == 0) continue;

    // theta = kT/m

    if (rmass) theta = force->boltz*walltemp[m]/(rmass[i]*force->mvv2e);
    else theta = force->boltz*walltemp[m]/(mass[type[i]]*force->mvv2e);
    factor = sqrt(theta);
      
    // time travelling after collision

    timecol = (x[i][dim]-coord) / v[i][dim]; 

    // only needed for Maxwell model
    
    if (reflectionstyle == MAXWELL) difftest = random->uniform();
    
    for (dir = 0; dir < 3; dir++) {

      // pushing back atoms to wall position, assign new reflected velocity
      
      x[i][dir] -= v[i][dir]*timecol; 
      
      // diffusive reflection
      
      if (reflectionstyle  == DIFFUSIVE) {
        if (dir != dim) 
          v[i][dir] = wallvel[m][dir] + random->gaussian(0,factor);
        else v[i][dir] =  sign*random->rayleigh(factor);
        
      // Maxwell reflection
          
      } else if (reflectionstyle  == MAXWELL) {
        if (difftest < wallaccom[m][dir]) {
          if (dir != dim) 
            v[i][dir] = wallvel[m][dir] + random->gaussian(0,factor);
          else v[i][dir] =  sign*random->rayleigh(factor);
        } else {
          if (dir == dim) v[i][dir] = -v[i][dir];
        }

      // Cercignani Lampis reflection

      } else if (reflectionstyle  == CERCIGNANILAMPIS) {
        if (dir != dim)
          v[i][dir] = wallvel[m][dir] + 
            random->gaussian((1-wallaccom[m][dir]) * 
                             (v[i][dir] - wallvel[m][dir]),
                             sqrt((2-wallaccom[m][dir]) * 
                                  wallaccom[m][dir]*theta));
        else v[i][dir] = random->besselexp(theta,wallaccom[m][dir],v[i][dir]);
      }

      // update new position due to the new velocity
      
      if (dir != dim) x[i][dir] += v[i][dir]*timecol;
      else x[i][dir] = coord + v[i][dir]*timecol;
    }
  }
}
