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

#include <mpi.h>
#include <cmath>
#include "min_spinmin.h"
#include "universe.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "output.h"
#include "timer.h"
#include "error.h"

#include <cstdlib>
#include <cstring>
#include "modify.h"
#include "math_special.h"
#include "math_const.h"
#include "fix_neb_spin.h"

using namespace LAMMPS_NS;
using namespace MathConst;

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8

#define DELAYSTEP 5

/* ---------------------------------------------------------------------- */

MinSpinMin::MinSpinMin(LAMMPS *lmp) : Min(lmp) {}

/* ---------------------------------------------------------------------- */

void MinSpinMin::init()
{
  Min::init();

  dt = update->dt;
  last_negative = update->ntimestep;

  // test dts
  dts = dt;

}

/* ---------------------------------------------------------------------- */

void MinSpinMin::setup_style()
{
  double **v = atom->v;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    v[i][0] = v[i][1] = v[i][2] = 0.0;
}

/* ----------------------------------------------------------------------
   set current vector lengths and pointers
   called after atoms have migrated
------------------------------------------------------------------------- */

void MinSpinMin::reset_vectors()
{
  // atomic dof

  // not really good size => sp is 4N vector
  nvec = 4 * atom->nlocal;
  if (nvec) spvec = atom->sp[0];
  
  nvec = 3 * atom->nlocal;
  if (nvec) fmvec = atom->fm[0];
  
  if (nvec) xvec = atom->x[0];
  if (nvec) fvec = atom->f[0];
}

/* ----------------------------------------------------------------------
   minimization via QuickMin damped dynamics
------------------------------------------------------------------------- */

int MinSpinMin::iterate(int maxiter)
{
  bigint ntimestep;
  //double vmax,vdotf,vdotfall,fdotf,fdotfall,scale;
  //double dtvone,dtv,dtf,dtfm;
  int flag,flagall;

  //alpha_final = 0.0;
  
  // search for and allocate neb_spin fix
  
  //int ineb;
  //for (ineb = 0; ineb < modify->nfix; ineb++)
  //  if (strcmp(modify->fix[ineb]->style,"neb/spin") == 0) break;
  //if (ineb == modify->nfix) error->all(FLERR,"spinmin requires use of fix neb/spin");
  //fneb = (FixNEB_spin *) modify->fix[ineb];

  for (int iter = 0; iter < maxiter; iter++) {

    if (timer->check_timeout(niter))
      return TIMEOUT;

    ntimestep = ++update->ntimestep;
    niter++;

    // optimize timestep accross processes / replicas
    
    //dts = fneb->evaluate_dt();
    dts = evaluate_dt();
   
    // apply damped precessional dynamics to the spins
      
    //fneb->advance_spins(dts);
    advance_spins(dts);


    //// zero velocity if anti-parallel to force
    //// else project velocity in direction of force

    //double **v = atom->v;
    //double **f = atom->f;
    //int nlocal = atom->nlocal;

    //vdotf = 0.0;
    //for (int i = 0; i < nlocal; i++)
    //  vdotf += v[i][0]*f[i][0] + v[i][1]*f[i][1] + v[i][2]*f[i][2];
    //MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,world);

    // sum vdotf over replicas, if necessary
    // this communicator would be invalid for multiprocess replicas

    //if (update->multireplica == 1) {
    //  vdotf = vdotfall;
    //  MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
    //}

    //if (vdotfall < 0.0) {
    //  last_negative = ntimestep;
    //  for (int i = 0; i < nlocal; i++)
    //    v[i][0] = v[i][1] = v[i][2] = 0.0;

    //} else {
    //  fdotf = 0.0;
    //  for (int i = 0; i < nlocal; i++)
    //    fdotf += f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2];
    //  MPI_Allreduce(&fdotf,&fdotfall,1,MPI_DOUBLE,MPI_SUM,world);

    //  // sum fdotf over replicas, if necessary
    //  // this communicator would be invalid for multiprocess replicas

    //  if (update->multireplica == 1) {
    //    fdotf = fdotfall;
    //    MPI_Allreduce(&fdotf,&fdotfall,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
    //  }

    //  if (fdotfall == 0.0) scale = 0.0;
    //  else scale = vdotfall/fdotfall;
    //  for (int i = 0; i < nlocal; i++) {
    //    v[i][0] = scale * f[i][0];
    //    v[i][1] = scale * f[i][1];
    //    v[i][2] = scale * f[i][2];
    //  }
    //}

    //// limit timestep so no particle moves further than dmax

    //double *rmass = atom->rmass;
    //double *mass = atom->mass;
    //int *type = atom->type;

    //dtvone = dt;

    //for (int i = 0; i < nlocal; i++) {
    //  vmax = MAX(fabs(v[i][0]),fabs(v[i][1]));
    //  vmax = MAX(vmax,fabs(v[i][2]));
    //  if (dtvone*vmax > dmax) dtvone = dmax/vmax;
    //}
    //MPI_Allreduce(&dtvone,&dtv,1,MPI_DOUBLE,MPI_MIN,world);

    //// min dtv over replicas, if necessary
    //// this communicator would be invalid for multiprocess replicas

    //if (update->multireplica == 1) {
    //  dtvone = dtv;
    //  MPI_Allreduce(&dtvone,&dtv,1,MPI_DOUBLE,MPI_MIN,universe->uworld);
    //}

    //dtf = dtv * force->ftm2v;

    //// Euler integration step

    //double **x = atom->x;

    //if (rmass) {
    //  for (int i = 0; i < nlocal; i++) {
    //    dtfm = dtf / rmass[i];
    //    x[i][0] += dtv * v[i][0];
    //    x[i][1] += dtv * v[i][1];
    //    x[i][2] += dtv * v[i][2];
    //    v[i][0] += dtfm * f[i][0];
    //    v[i][1] += dtfm * f[i][1];
    //    v[i][2] += dtfm * f[i][2];
    //  }
    //} else {
    //  for (int i = 0; i < nlocal; i++) {
    //    dtfm = dtf / mass[type[i]];
    //    x[i][0] += dtv * v[i][0];
    //    x[i][1] += dtv * v[i][1];
    //    x[i][2] += dtv * v[i][2];
    //    v[i][0] += dtfm * f[i][0];
    //    v[i][1] += dtfm * f[i][1];
    //    v[i][2] += dtfm * f[i][2];
    //  }
    //}

    eprevious = ecurrent;
    ecurrent = energy_force(0);
    neval++;

    //// energy tolerance criterion
    //// only check after DELAYSTEP elapsed since velocties reset to 0
    //// sync across replicas if running multi-replica minimization

    if (update->etol > 0.0 && ntimestep-last_negative > DELAYSTEP) {
      if (update->multireplica == 0) {
        if (fabs(ecurrent-eprevious) <
            update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
          return ETOL;
      } else {
        if (fabs(ecurrent-eprevious) <
            update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
          flag = 0;
        else flag = 1;
        MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,universe->uworld);
        if (flagall == 0) return ETOL;
      }
    }

    //// magnetic force tolerance criterion
    //// sync across replicas if running multi-replica minimization

    //if (update->fmtol > 0.0) {
    //  fmdotfm = fmnorm_sqr();
    //  if (update->multireplica == 0) {
    //    if (fmdotfm < update->fmtol*update->fmtol) return FTOL;
    //  } else {
    //    if (fmdotfm < update->fmtol*update->fmtol) flag = 0;
    //    else flag = 1;
    //    MPI_Allreduce(&fmlag,&fmlagall,1,MPI_INT,MPI_SUM,universe->uworld);
    //    if (fmlagall == 0) return FTOL;
    //  }
    //}
    
    //// force tolerance criterion
    //// sync across replicas if running multi-replica minimization

    //if (update->ftol > 0.0) {
    //  fdotf = fnorm_sqr();
    //  if (update->multireplica == 0) {
    //    if (fdotf < update->ftol*update->ftol) return FTOL;
    //  } else {
    //    if (fdotf < update->ftol*update->ftol) flag = 0;
    //    else flag = 1;
    //    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,universe->uworld);
    //    if (flagall == 0) return FTOL;
    //  }
    //}

    //// output for thermo, dump, restart files

    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }

  return MAXITER;
}

/* ----------------------------------------------------------------------
   evaluate max timestep
---------------------------------------------------------------------- */

double MinSpinMin::evaluate_dt()
{
  double dtmax;
  double fmsq;
  double fmaxsqone,fmaxsqloc,fmaxsqall;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **fm = atom->fm;

  // finding max fm on this proc. 
  
  fmsq = fmaxsqone = fmaxsqloc = fmaxsqall = 0.0;
  for (int i = 0; i < nlocal; i++) {
    fmsq = fm[i][0]*fm[i][0]+fm[i][1]*fm[i][1]+fm[i][2]*fm[i][2];
    fmaxsqone = MAX(fmaxsqone,fmsq);
  }
  //printf("test inside evaluate dt, fmaxsqone = %g \n",fmaxsqone);
  //for (int i = 0; i < nlocal; i++)
  //  if (mask[i] & groupbit) {
  //    fmsq = fm[i][0]*fm[i][0]+fm[i][1]*fm[i][1]+fm[i][2]*fm[i][2];
  //    fmaxsqone = MAX(fmaxsqone,fmsq);
  //  }

  // finding max fm on this replica
 
  fmaxsqloc = fmaxsqone; 
  MPI_Allreduce(&fmaxsqone,&fmaxsqloc,1,MPI_DOUBLE,MPI_MAX,world); 
  
  // finding max fm over all replicas, if necessary
  // this communicator would be invalid for multiprocess replicas

  fmaxsqall = fmaxsqloc;
  if (update->multireplica == 1) {
    fmaxsqall = fmaxsqloc;
    MPI_Allreduce(&fmaxsqloc,&fmaxsqall,1,MPI_DOUBLE,MPI_MAX,universe->uworld);
  }

  //if (fmaxsqall < fmaxsqloc)
  //  error->all(FLERR,"Incorrect fmaxall calc.");

  // define max timestep
  // dividing by 10 the inverse of max frequency
  //printf("test inside evaluate dt, fmaxsqall = %g \n",fmaxsqall);

  dtmax = MY_2PI/(10.0*sqrt(fmaxsqall));
  //printf("test inside evaluate dt, dtmax = %g \n",dtmax);

  return dtmax;
}

/* ----------------------------------------------------------------------
   geometric damped advance of spins
---------------------------------------------------------------------- */

void MinSpinMin::advance_spins(double dts)
{
  //int j=0;
  //int *sametag = atom->sametag;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **sp = atom->sp;
  double **fm = atom->fm;
  double tdampx,tdampy,tdampz;
  double msq,scale,fm2,energy,dts2;
  double alpha;
  //double spi[3],fmi[3];
  double cp[3],g[3];

  //cp[0] = cp[1] = cp[2] = 0.0;
  //g[0] = g[1] = g[2] = 0.0;
  dts2 = dts*dts;		

  // fictitious Gilbert damping of 1
  alpha = 1.0;

  //printf("test inside spinmin, dts %g \n",dts);
  //printf("test inside spinmin, fmi i=%d, %g %g %g \n",1,fm[1][0],fm[1][1],fm[1][2]);
  
  // loop on all spins on proc.

  //if (ireplica != nreplica-1 && ireplica != 0)
  //  for (int i = 0; i < nlocal; i++)
  //    if (mask[i] & groupbit) {
  for (int i = 0; i < nlocal; i++) {

    //spi[0] = sp[i][0];
    //spi[1] = sp[i][1];
    //spi[2] = sp[i][2];
    //
    //fmi[0] = fm[i][0];
    //fmi[1] = fm[i][1];
    //fmi[2] = fm[i][2];
    
    // calc. damping torque
    
    tdampx = -alpha*(fm[i][1]*sp[i][2] - fm[i][2]*sp[i][1]);
    tdampy = -alpha*(fm[i][2]*sp[i][0] - fm[i][0]*sp[i][2]);
    tdampz = -alpha*(fm[i][0]*sp[i][1] - fm[i][1]*sp[i][0]);
    
    //printf("for %d, test tdamp: %g %g %g \n",i,tdampx,tdampy,tdampz);
    
    // apply advance algorithm (geometric, norm preserving)
    
    fm2 = (tdampx*tdampx+tdampy*tdampy+tdampz*tdampz);
    energy = (sp[i][0]*tdampx)+(sp[i][1]*tdampy)+(sp[i][2]*tdampz);
    
    cp[0] = tdampy*sp[i][2]-tdampz*sp[i][1];
    cp[1] = tdampz*sp[i][0]-tdampx*sp[i][2];
    cp[2] = tdampx*sp[i][1]-tdampy*sp[i][0];
    
    //printf("for %d, test cp: %g %g %g \n",i,cp[0],cp[1],cp[2]);
    
    g[0] = sp[i][0]+cp[0]*dts;
    g[1] = sp[i][1]+cp[1]*dts;
    g[2] = sp[i][2]+cp[2]*dts;
    		
    //g[0] += (fm[i][0]*energy-0.5*sp[i][0]*fm2)*0.5*dts2;
    //g[1] += (fm[i][1]*energy-0.5*sp[i][1]*fm2)*0.5*dts2;
    //g[2] += (fm[i][2]*energy-0.5*sp[i][2]*fm2)*0.5*dts2;
    g[0] += (tdampx*energy-0.5*sp[i][0]*fm2)*0.5*dts2;
    g[1] += (tdampy*energy-0.5*sp[i][1]*fm2)*0.5*dts2;
    g[2] += (tdampz*energy-0.5*sp[i][2]*fm2)*0.5*dts2;
    		
    g[0] /= (1+0.25*fm2*dts2);
    g[1] /= (1+0.25*fm2*dts2);
    g[2] /= (1+0.25*fm2*dts2);
   
    //printf("test inside spinmin, spi i=%d, %g %g %g \n",i,sp[i][0],sp[i][1],sp[i][2]);
    //printf("test inside spinmin, fmi i=%d, %g %g %g \n",i,fm[i][0],fm[i][1],fm[i][2]);
    //printf("for %d, test g: %g %g %g \n",i,g[0],g[1],g[2]);

    sp[i][0] = g[0];
    sp[i][1] = g[1];
    sp[i][2] = g[2];			
    
    // renormalization (check if necessary)
    
    msq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
    scale = 1.0/sqrt(msq);
    sp[i][0] *= scale;
    sp[i][1] *= scale;
    sp[i][2] *= scale;
    
    // comm. sp[i] to atoms with same tag (for serial algo)
    
    // no need for simplecticity
    //if (sector_flag == 0) {
    //  if (sametag[i] >= 0) {
    //    j = sametag[i];
    //    while (j >= 0) {
    //      sp[j][0] = sp[i][0];
    //      sp[j][1] = sp[i][1];
    //      sp[j][2] = sp[i][2];
    //      j = sametag[j];
    //    }
    //  }
    //}
    //
  }

  //printf("test inside spinmin, dts = %g \n",dts);
  //printf("test inside spinmin, fmi i=%d, %g %g %g \n",1,fm[1][0],fm[1][1],fm[1][2]);
  //printf("test inside spinmin, spi i=%d, %g %g %g \n",1,sp[1][0],sp[1][1],sp[1][2]);
}

