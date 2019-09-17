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

#include "fix_wall_cl.h"
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "lattice.h"
#include "input.h"
#include "variable.h"
#include "error.h"
#include "force.h"
#include "random_mars.h"
#include <stdlib.h> 
/*#include <iostream> */

using namespace LAMMPS_NS;
using namespace FixConst;
/* using namespace std; */

enum{XLO=0,XHI=1,YLO=2,YHI=3,ZLO=4,ZHI=5};
enum{NONE=0,EDGE,CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */

FixWallCL::FixWallCL(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), random(NULL),
  nwall(0)
{
  int dir;
  if (narg < 4) error->all(FLERR,"Illegal fix wall/cl command");

  dynamic_group_allow = 1;

  // parse args

  nwall = 0;
  int scaleflag = 1;
   
  do seedfix = rand(); while (seedfix > 900000000);
  random = new RanMars(lmp,seedfix);

  int iarg = 3;
  while (iarg < narg) {
    if ((strcmp(arg[iarg],"xlo") == 0) || (strcmp(arg[iarg],"xhi") == 0) ||
        (strcmp(arg[iarg],"ylo") == 0) || (strcmp(arg[iarg],"yhi") == 0) ||
        (strcmp(arg[iarg],"zlo") == 0) || (strcmp(arg[iarg],"zhi") == 0)) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall/cl command");

      int newwall;
      if (strcmp(arg[iarg],"xlo") == 0) newwall = XLO;
      else if (strcmp(arg[iarg],"xhi") == 0) newwall = XHI;
      else if (strcmp(arg[iarg],"ylo") == 0) newwall = YLO;
      else if (strcmp(arg[iarg],"yhi") == 0) newwall = YHI;
      else if (strcmp(arg[iarg],"zlo") == 0) newwall = ZLO;
      else if (strcmp(arg[iarg],"zhi") == 0) newwall = ZHI;

      for (int m = 0; (m < nwall) && (m < 6); m++)
        if (newwall == wallwhich[m])
          error->all(FLERR,"Wall defined twice in fix wall/cl command");

      wallwhich[nwall] = newwall;
      if (strcmp(arg[iarg+1],"EDGE") == 0) {
        wallstyle[nwall] = EDGE;
        int dim = wallwhich[nwall] / 2;
        int side = wallwhich[nwall] % 2;
        if (side == 0) coord0[nwall] = domain->boxlo[dim];
        else coord0[nwall] = domain->boxhi[dim];
      } else if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        wallstyle[nwall] = VARIABLE;
        int n = strlen(&arg[iarg+1][2]) + 1;
        varstr[nwall] = new char[n];
        strcpy(varstr[nwall],&arg[iarg+1][2]);
      } else {
        wallstyle[nwall] = CONSTANT;
        coord0[nwall] = force->numeric(FLERR,arg[iarg+1]);
      }
      walltemp[nwall]= force->numeric(FLERR,arg[iarg+2]);
      for (dir = 0; dir < 3; dir++) {
      wallvel[nwall][dir]= force->numeric(FLERR,arg[iarg+dir+3]);
      wallaccom[nwall][dir]= force->numeric(FLERR,arg[iarg+dir+6]);

      }  
      nwall++;
      iarg += 9;

    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal wall/cl command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix wall/cl command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix wall/cl command");
  }

  // error check

  if (nwall == 0) error->all(FLERR,"Illegal fix wall command");

  for (int m = 0; m < nwall; m++) {
    if ((wallwhich[m] == XLO || wallwhich[m] == XHI) && domain->xperiodic)
      error->all(FLERR,"Cannot use fix wall/cl in periodic dimension");
    if ((wallwhich[m] == YLO || wallwhich[m] == YHI) && domain->yperiodic)
      error->all(FLERR,"Cannot use fix wall/cl in periodic dimension");
    if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->zperiodic)
      error->all(FLERR,"Cannot use fix wall/cl in periodic dimension");
  }

  for (int m = 0; m < nwall; m++)
    if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->dimension == 2)
      error->all(FLERR,
                 "Cannot use fix wall/cl zlo/zhi for a 2d simulation");

  // scale factors for CONSTANT and VARIABLE walls

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

  // set varflag if any wall positions are variable

  varflag = 0;
  for (int m = 0; m < nwall; m++)
    if (wallstyle[m] == VARIABLE) varflag = 1;
}

/* ---------------------------------------------------------------------- */

FixWallCL::~FixWallCL()
{
  if (copymode) return;

  for (int m = 0; m < nwall; m++)
    if (wallstyle[m] == VARIABLE) delete [] varstr[m];
  delete random;
}

/* ---------------------------------------------------------------------- */

int FixWallCL::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallCL::init()
{
  for (int m = 0; m < nwall; m++) {
    if (wallstyle[m] != VARIABLE) continue;
    varindex[m] = input->variable->find(varstr[m]);
    if (varindex[m] < 0)
      error->all(FLERR,"Variable name for fix wall/cl does not exist");
    if (!input->variable->equalstyle(varindex[m]))
      error->all(FLERR,"Variable for fix wall/cl is invalid style");
  }

  int nrigid = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) nrigid++;

  if (nrigid && comm->me == 0)
    error->warning(FLERR,"Should not allow rigid bodies to bounce off "
                   "relecting walls");
}


/* ---------------------------------------------------------------------- */

void FixWallCL::post_integrate()
{
  int i,dim,side,dir;
  double coord,vsave,factor,theta,time1,time2;

  double *rmass;
  double *mass = atom->mass;


  // coord = current position of wall
  // evaluate variable if necessary, wrap with clear/add


  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  

  if (varflag) modify->clearstep_compute();

  for (int m = 0; m < nwall; m++) {
    if (wallstyle[m] == VARIABLE) {
      coord = input->variable->compute_equal(varindex[m]);
      if (wallwhich[m] < YLO) coord *= xscale;
      else if (wallwhich[m] < ZLO) coord *= yscale;
      else coord *= zscale;
    } else coord = coord0[m];

    dim = wallwhich[m] / 2;
    side = wallwhich[m] % 2;
  
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (side == 0) {
          if (x[i][dim] < coord) {
           if (rmass) theta = force->boltz*walltemp[m]/(rmass[i]*force->mvv2e);// theta = kT/m
           else theta = force->boltz*walltemp[m]/(mass[type[i]]*force->mvv2e);

             time2= (x[i][dim] - coord)/(v[i][dim]-wallvel[m][dim]);// after collision


            for (dir = 0; dir < 3; dir++) 
            if (dir != dim) {
            x[i][dir] -=v[i][dir]*time2;
            v[i][dir] = wallvel[m][dir] + random->gaussian((1-wallaccom[m][dir])*(v[i][dir]-wallvel[m][dir]),sqrt((2-wallaccom[m][dir])*wallaccom[m][dir]*theta)) ;
            x[i][dir] +=v[i][dir]*time2;
            } else {
            v[i][dir] = wallvel[m][dir] + random->besselexp(theta, wallaccom[m][dir], v[i][dir]-wallvel[m][dir]) ; 
            x[i][dir] = coord + v[i][dir]*time2;
            }
            

          }
        } else {
          if (x[i][dim] > coord) {
           if (rmass) theta = force->boltz*walltemp[m]/(rmass[i]*force->mvv2e);
           else theta = force->boltz*walltemp[m]/(mass[type[i]]*force->mvv2e);


             time2= (x[i][dim] - coord)/(v[i][dim]-wallvel[m][dim]);// after collision


            for (dir = 0; dir < 3; dir++) 
            if (dir != dim) {
            x[i][dir] -=v[i][dir]*time2;
            v[i][dir] = wallvel[m][dir] + random->gaussian((1-wallaccom[m][dir])*(v[i][dir]-wallvel[m][dir]),sqrt((2-wallaccom[m][dir])*wallaccom[m][dir]*theta)) ;
            x[i][dir] +=v[i][dir]*time2;
            } else {
            v[i][dir] = wallvel[m][dir] + random->besselexp(theta, wallaccom[m][dir], v[i][dir]-wallvel[m][dir]) ;  
            x[i][dir] = coord + v[i][dir]*time2;
            }



          }
        }
      }
  }

  if (varflag) modify->addstep_compute(update->ntimestep + 1);
}
