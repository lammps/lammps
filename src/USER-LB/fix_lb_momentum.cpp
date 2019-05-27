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
   Contributing authors: Frances Mackay, Santtu Ollila, Colin Denniston (UWO)

   Based on fix_momentum,
   Contributing author: Naveen Michaud-Agrawal (Johns Hopkins U)
------------------------------------------------------------------------- */

#include <cstdlib>
#include <cstring>
#include "fix_lb_momentum.h"
#include "atom.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "fix_lb_fluid.h"
#include "modify.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;


/* ---------------------------------------------------------------------- */

FixLbMomentum::FixLbMomentum(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix lb/momentum command");
  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix lb/momentum command");

  linear = 1;
  xflag = 1;
  yflag = 1;
  zflag = 1;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"linear") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix lb/momentum command");
      linear = 1;
      xflag = atoi(arg[iarg+1]);
      yflag = atoi(arg[iarg+2]);
      zflag = atoi(arg[iarg+3]);
      iarg += 4;
    } else error->all(FLERR,"Illegal fix lb/momentum command");
  }

  if (linear == 0)
    error->all(FLERR,"Illegal fix lb/momentum command");

  if (linear)
    if (xflag < 0 || xflag > 1 || yflag < 0 || yflag > 1 ||
        zflag < 0 || zflag > 1) error->all(FLERR,"Illegal fix lb/momentum command");

  // cannot have 0 atoms in group

  if (group->count(igroup) == 0.0)
    error->all(FLERR,"Fix lb/momentum group has no atoms");

  for(int ifix=0; ifix<modify->nfix; ifix++)
    if(strcmp(modify->fix[ifix]->style,"lb/fluid")==0)
      fix_lb_fluid = (FixLbFluid *)modify->fix[ifix];



}

/* ---------------------------------------------------------------------- */

int FixLbMomentum::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLbMomentum::init()
{
  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

void FixLbMomentum::end_of_step()
{
  double masslb,masslbloc;
  double momentumlbloc[3],momentumlb[3];
  double vcmtotal[3];
  const int numvel = fix_lb_fluid->numvel;
  double etacov[19]; // = double etacov[numvel]; i.e. 15 or 19
  double rho;

  if (linear) {
    double vcm[3];
    group->vcm(igroup,masstotal,vcm);

    // adjust velocities by vcm to zero linear momentum
    // only adjust a component if flag is set

    double **v = atom->v;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    int typeLB = fix_lb_fluid->typeLB;
    int subNbx = fix_lb_fluid->subNbx;
    int subNby = fix_lb_fluid->subNby;
    int subNbz = fix_lb_fluid->subNbz;
    double ***density_lb = fix_lb_fluid->density_lb;
    double ****f_lb = fix_lb_fluid->f_lb;
    double ****u_lb = fix_lb_fluid->u_lb;
    double dm_lb = fix_lb_fluid->dm_lb;
    double dx_lb = fix_lb_fluid->dx_lb;
    double dt_lb = fix_lb_fluid->dt_lb;
    double *Ng_lb = fix_lb_fluid->Ng_lb;
    double *w_lb = fix_lb_fluid->w_lb;
    double **mg_lb = fix_lb_fluid->mg_lb;
    double ucmx,ucmy,ucmz;

    //Calculate the total fluid mass and momentum.
    masslbloc = 0.0;
    momentumlbloc[0] = momentumlbloc[1] = momentumlbloc[2] = 0.0;

     for(int i = 1; i<subNbx-1; i++)
      for(int j = 1; j<subNby-1; j++)
        for(int k = 1; k<subNbz-1; k++){
          masslbloc += density_lb[i][j][k];
          momentumlbloc[0] += density_lb[i][j][k]*u_lb[i][j][k][0];
          momentumlbloc[1] += density_lb[i][j][k]*u_lb[i][j][k][1];
          momentumlbloc[2] += density_lb[i][j][k]*u_lb[i][j][k][2];
        }
    MPI_Allreduce(&masslbloc,&masslb,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&momentumlbloc[0],&momentumlb[0],3,MPI_DOUBLE,MPI_SUM,world);

    momentumlb[0] *= dm_lb*dx_lb/dt_lb;
    momentumlb[1] *= dm_lb*dx_lb/dt_lb;
    momentumlb[2] *= dm_lb*dx_lb/dt_lb;
    masslb *= dm_lb;

    //Calculate the center of mass velocity of the entire system.
    vcmtotal[0] = (masstotal*vcm[0] + momentumlb[0])/(masslb + masstotal);
    vcmtotal[1] = (masstotal*vcm[1] + momentumlb[1])/(masslb + masstotal);
    vcmtotal[2] = (masstotal*vcm[2] + momentumlb[2])/(masslb + masstotal);

    //Subtract vcm from the particles.
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (xflag) v[i][0] -= vcmtotal[0];
        if (yflag) v[i][1] -= vcmtotal[1];
        if (zflag) v[i][2] -= vcmtotal[2];
      }

    vcmtotal[0] *= dt_lb/dx_lb;
    vcmtotal[1] *= dt_lb/dx_lb;
    vcmtotal[2] *= dt_lb/dx_lb;

    ucmx = ucmy = ucmz = 0.0;
    double density_old;
    double u_old[3];
    int **e = fix_lb_fluid->e;

    //Subtract vcm from the fluid.
    for(int i=0; i<subNbx; i++)
      for(int j=0; j<subNby; j++)
        for(int k=0; k<subNbz; k++){
          rho = density_lb[i][j][k];
          if(xflag) ucmx = vcmtotal[0];
          if(yflag) ucmy = vcmtotal[1];
          if(zflag) ucmz = vcmtotal[2];
          if(numvel==15){
            etacov[0]=0.0;
            etacov[1]=rho*ucmx;
            etacov[2]=rho*ucmy;
            etacov[3]=rho*ucmz;
            etacov[4]=rho*(2.*u_lb[i][j][k][0]*ucmx-ucmx*ucmx);
            etacov[5]=rho*(2.*u_lb[i][j][k][1]*ucmy-ucmy*ucmy);
            etacov[6]=rho*(2.*u_lb[i][j][k][2]*ucmz-ucmz*ucmz);
            etacov[7]=rho*(u_lb[i][j][k][0]*ucmy+u_lb[i][j][k][1]*ucmx-ucmx*ucmy);
            etacov[8]=rho*(u_lb[i][j][k][1]*ucmz+u_lb[i][j][k][2]*ucmy-ucmy*ucmz);
            etacov[9]=rho*(u_lb[i][j][k][0]*ucmz+u_lb[i][j][k][2]*ucmx-ucmx*ucmz);
            etacov[10]=0.0;
            etacov[11]=0.0;
            etacov[12]=0.0;
            etacov[13]=rho*(u_lb[i][j][k][0]*u_lb[i][j][k][1]*ucmz+u_lb[i][j][k][0]*ucmy*u_lb[i][j][k][2]-
                            u_lb[i][j][k][0]*ucmy*ucmz+ucmx*u_lb[i][j][k][1]*u_lb[i][j][k][2]-
                            ucmx*u_lb[i][j][k][1]*ucmz-ucmx*ucmy*u_lb[i][j][k][2]+
                            ucmx*ucmy*ucmz);
            etacov[14]=0.0;
          } else {
            etacov[0] = 0.0;
            etacov[1] = rho*ucmx;
            etacov[2] = rho*ucmy;
            etacov[3] = rho*ucmz;
            etacov[4]=rho*(2.*u_lb[i][j][k][0]*ucmx-ucmx*ucmx);
            etacov[5]=rho*(2.*u_lb[i][j][k][1]*ucmy-ucmy*ucmy);
            etacov[6]=rho*(2.*u_lb[i][j][k][2]*ucmz-ucmz*ucmz);
            etacov[7]=rho*(u_lb[i][j][k][0]*ucmy+u_lb[i][j][k][1]*ucmx-ucmx*ucmy);
            etacov[8]=rho*(u_lb[i][j][k][0]*ucmz+u_lb[i][j][k][2]*ucmx-ucmx*ucmz);
            etacov[9]=rho*(u_lb[i][j][k][1]*ucmz+u_lb[i][j][k][2]*ucmy-ucmy*ucmz);
            etacov[10] = 0.0;
            etacov[11] = 0.0;
            etacov[12] = 0.0;
            etacov[13] = 0.0;
            etacov[14] = 0.0;
            etacov[15] = 0.0;
            etacov[16] = 0.0;
            etacov[17] = 0.0;
            etacov[18] = 0.0;
          }

          for(int l=0; l<numvel; l++)
            for(int ii=0; ii<numvel; ii++){
              f_lb[i][j][k][l] -= w_lb[l]*mg_lb[ii][l]*etacov[ii]*Ng_lb[ii];
            }


          if(typeLB == 2){
            double ****feqold = fix_lb_fluid->feqold;
            double ****feqoldn = fix_lb_fluid->feqoldn;
            density_old = 0.0;
            u_old[0] = u_old[1] = u_old[2] = 0.0;
            for(int l=0; l<numvel; l++){
              density_old += feqold[i][j][k][l];
              u_old[0] += feqold[i][j][k][l]*e[l][0];
              u_old[1] += feqold[i][j][k][l]*e[l][1];
              u_old[2] += feqold[i][j][k][l]*e[l][2];
            }
            u_old[0] = u_old[0]/density_old;
            u_old[1] = u_old[1]/density_old;
            u_old[2] = u_old[2]/density_old;

            if(numvel==15){
              etacov[0]=0.0;
              etacov[1]=density_old*ucmx;
              etacov[2]=density_old*ucmy;
              etacov[3]=density_old*ucmz;
              etacov[4]=density_old*(2.*u_old[0]*ucmx-ucmx*ucmx);
              etacov[5]=density_old*(2.*u_old[1]*ucmy-ucmy*ucmy);
              etacov[6]=density_old*(2.*u_old[2]*ucmz-ucmz*ucmz);
              etacov[7]=density_old*(u_old[0]*ucmy+u_old[1]*ucmx-ucmx*ucmy);
              etacov[8]=density_old*(u_old[1]*ucmz+u_old[2]*ucmy-ucmy*ucmz);
              etacov[9]=density_old*(u_old[0]*ucmz+u_old[2]*ucmx-ucmx*ucmz);
              etacov[10]=0.0;
              etacov[11]=0.0;
              etacov[12]=0.0;
              etacov[13]=density_old*(u_old[0]*u_old[1]*ucmz+u_old[0]*ucmy*u_old[2]-
                                      u_old[0]*ucmy*ucmz+ucmx*u_old[1]*u_old[2]-
                                      ucmx*u_old[1]*ucmz-ucmx*ucmy*u_old[2]+
                                      ucmx*ucmy*ucmz);
              etacov[14]=0.0;
            } else {
              etacov[0] = 0.0;
              etacov[1] = density_old*ucmx;
              etacov[2] = density_old*ucmy;
              etacov[3] = density_old*ucmz;
              etacov[4] = density_old*(2.*u_lb[i][j][k][0]*ucmx-ucmx*ucmx);
              etacov[5] = density_old*(2.*u_lb[i][j][k][1]*ucmy-ucmy*ucmy);
              etacov[6] = density_old*(2.*u_lb[i][j][k][2]*ucmz-ucmz*ucmz);
              etacov[7] = density_old*(u_lb[i][j][k][0]*ucmy+u_lb[i][j][k][1]*ucmx-ucmx*ucmy);
              etacov[8] = density_old*(u_lb[i][j][k][0]*ucmz+u_lb[i][j][k][2]*ucmx-ucmx*ucmz);
              etacov[9] = density_old*(u_lb[i][j][k][1]*ucmz+u_lb[i][j][k][2]*ucmy-ucmy*ucmz);
              etacov[10] = 0.0;
              etacov[11] = 0.0;
              etacov[12] = 0.0;
              etacov[13] = 0.0;
              etacov[14] = 0.0;
              etacov[15] = 0.0;
              etacov[16] = 0.0;
              etacov[17] = 0.0;
              etacov[18] = 0.0;
            }

            for(int l=0; l<numvel; l++)
              for(int ii=0; ii<numvel; ii++){
                feqold[i][j][k][l] -= w_lb[l]*mg_lb[ii][l]*etacov[ii]*Ng_lb[ii];
                feqoldn[i][j][k][l] -= w_lb[l]*mg_lb[ii][l]*etacov[ii]*Ng_lb[ii];
              }
          }
        }
  }

  fix_lb_fluid->parametercalc_full();

}
