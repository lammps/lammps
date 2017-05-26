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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "fix_nve_spin.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "error.h"
#include "math_vector.h"
#include "math_extra.h"
#include "math_const.h"
#include "modify.h" 

#include "pair.h"
#include "timer.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace MathExtra;

enum{NONE,SPIN};

/* ---------------------------------------------------------------------- */

FixNVESpin::FixNVESpin(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
	
  if (narg < 3) error->all(FLERR,"Illegal fix nve/spin command");	

  time_integrate = 1;
  
  extra = NONE;

  int iarg = 2;
  if (strcmp(arg[iarg],"nve/spin") == 0) {
	  if (iarg+1 > narg) error->all(FLERR,"Illegal fix nve/spin command");
	  extra = SPIN;
  }
  
  // error checks
  if (extra == SPIN && !atom->mumag_flag)
    error->all(FLERR,"Fix nve/spin requires spin attribute mumag");
   
}

/* ---------------------------------------------------------------------- */

void FixNVESpin::init()
{
  FixNVE::init();       
  
  dts = update->dt;

  /*int idamp;
  for (idamp = 0; idamp < modify->nfix; idamp++)
    if (strstr(modify->fix[idamp]->style,"damping/spin")) break;
  if (idamp == modify->nfix)
    error->all(FLERR,"Integration of spin systems requires use of fix damping (set damping to 0.0 for NVE)");
  
  lockspindamping = (FixSpinDamping *) modify->fix[idamp]; 
  alpha_t = lockspindamping->get_damping(0); 
  */
}

/* ---------------------------------------------------------------------- */

void FixNVESpin::initial_integrate(int vflag)
{
  double dtfm,msq,scale,fm2,fmsq,sp2,spsq,energy;
  double cp[3],g[3]; 	
	
  double **x = atom->x;	
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass; 
  double *mass = atom->mass;  
  double **sp = atom->sp;
  double **fm = atom->fm;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;  
  int *type = atom->type;
  int *mask = atom->mask;  
  
  // update half v all particles
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass) dtfm = dtf / rmass[i];
      else dtfm = dtf / mass[type[i]]; 
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      }
  }

  // update half x for all particles
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i][0] += 0.5 * dtv * v[i][0];
      x[i][1] += 0.5 * dtv * v[i][1];
      x[i][2] += 0.5 * dtv * v[i][2];
      }
  }  


  size_t nbytes;
  nbytes = sizeof(double) * nlocal;
  int eflag = 3;
  // update sp for all particles
  if (extra == SPIN) {
      // Advance spins
      //See Omelyan et al., PRL 86, 2001 and P.W. Ma et al, PRB 83, 2011

#define CONC
#if defined CONC
      for (int i = 0; i < nlocal; i++){
         AdvanceSingleSpin(i,dts,sp,fm);
      }
#endif

//#define SEQ
#if defined SEQ
      //advance N-1 spins to half
      for (int i = 0; i < nlocal-1; i++){
         //Recomp field
         atom->avec->force_clear(0,nbytes);
         timer->stamp();
         modify->pre_force(vflag);
         timer->stamp(Timer::PAIR);
         force->pair->compute(eflag,vflag);
         timer->stamp(Timer::PAIR);
         modify->pre_reverse(eflag,vflag);
         timer->stamp(Timer::MODIFY);    
         comm->reverse_comm();
         timer->stamp(Timer::COMM); 
         modify->post_force(vflag);

         AdvanceSingleSpin(i,0.5*dts,sp,fm);
      }
      
      //advance N spin
      //Recomp field
      atom->avec->force_clear(0,nbytes);
      timer->stamp();
      modify->pre_force(vflag);
      timer->stamp(Timer::PAIR);
      force->pair->compute(eflag,vflag);
      timer->stamp(Timer::PAIR);
      modify->pre_reverse(eflag,vflag);
      timer->stamp(Timer::MODIFY);    
      comm->reverse_comm();
      timer->stamp(Timer::COMM); 
      modify->post_force(vflag);

      AdvanceSingleSpin(nlocal-1,dts,sp,fm);


      //advance N-1 spins to half
      for (int i = nlocal-2; i >= 0; i--){
         //Recomp field
         atom->avec->force_clear(0,nbytes);
         timer->stamp();
         modify->pre_force(vflag);
         timer->stamp(Timer::PAIR);
         force->pair->compute(eflag,vflag);
         timer->stamp(Timer::PAIR);
         modify->pre_reverse(eflag,vflag);
         timer->stamp(Timer::MODIFY);    
         comm->reverse_comm();
         timer->stamp(Timer::COMM); 
         modify->post_force(vflag);

         AdvanceSingleSpin(i,0.5*dts,sp,fm); 
      }

#endif

  }

#define FORCE_PRINT
#if defined FORCE_PRINT
     FILE* file_force=NULL;
     file_force=fopen("spin_force_Lammps.dat","a");
     fprintf(file_force,"---------------------------------- \n");
     for(int i=0;i<nlocal;i++){
        fprintf(file_force,"%d %lf %lf %lf \n",i,fm[i][0],fm[i][1],fm[i][2]);
     }
     if (file_force!=NULL) fclose(file_force);
#endif

}



/* ---------------------------------------------------------------------- */

void FixNVESpin::AdvanceSingleSpin(int i, double dts, double **sp, double **fm)
{
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;  
  int *type = atom->type;
  int *mask = atom->mask; 
  
  double dtfm,msq,scale,fm2,fmsq,sp2,spsq,energy;
  double cp[3],g[3]; 	

  cp[0] = cp[1] = cp[2] = 0.0;
  g[0] = g[1] = g[2] = 0.0;
  fm2 = (fm[i][0]*fm[i][0])+(fm[i][1]*fm[i][1])+(fm[i][2]*fm[i][2]);
  fmsq = sqrt(fm2);
  energy = (sp[i][0]*fm[i][0])+(sp[i][1]*fm[i][1])+(sp[i][2]*fm[i][2]);
			    
  cp[0] = fm[i][1]*sp[i][2]-fm[i][2]*sp[i][1];
  cp[1] = fm[i][2]*sp[i][0]-fm[i][0]*sp[i][2];
  cp[2] = fm[i][0]*sp[i][1]-fm[i][1]*sp[i][0];

//#define ALG_MA //Ma algo
#if defined ALG_MA
  double A[3];
  A[0]= A[1] = A[2] = 0.0;
  double xi=fmsq*dts;
  double zeta=0.0; // this is because omega contains already the transverse damping
  double chi=energy/fmsq;
  double expo=exp(2.0*zeta);
  double K1=1.0+expo+chi*(1.0-expo);
  double K=2.0*exp(zeta)/K1;
  double Ktrigo=1.0+0.25*xi*xi;
  double cosinus=(1.0-0.25*xi*xi)/Ktrigo; //cos(xi)
  double sinus=xi/Ktrigo; //sin(xi)
  A[0]=K*cosinus;
  A[1]=K*sinus;
  A[2]=(1.0-expo+chi*(1.0+expo-2.0*exp(zeta)*cosinus))/K1;

  g[0] = A[0]*sp[i][0]+A[1]*cp[0]/fmsq+A[2]*fm[i][0]/fmsq;
  g[1] = A[0]*sp[i][1]+A[1]*cp[1]/fmsq+A[2]*fm[i][0]/fmsq;
  g[2] = A[0]*sp[i][2]+A[1]*cp[2]/fmsq+A[2]*fm[i][0]/fmsq;
#endif

#define ALG_OM  //Omelyan algo
#if defined ALG_OM		  
  g[0] = sp[i][0]+cp[0]*dts;
  g[1] = sp[i][1]+cp[1]*dts;
  g[2] = sp[i][2]+cp[2]*dts;
			  
  g[0] += (fm[i][0]*energy-0.5*sp[i][0]*fm2)*0.5*dts*dts;
  g[1] += (fm[i][1]*energy-0.5*sp[i][1]*fm2)*0.5*dts*dts;
  g[2] += (fm[i][2]*energy-0.5*sp[i][2]*fm2)*0.5*dts*dts;
			  
  g[0] /= (1+0.25*fm2*dts*dts);
  g[1] /= (1+0.25*fm2*dts*dts);
  g[2] /= (1+0.25*fm2*dts*dts);
#endif
			  
  sp[i][0] = g[0];
  sp[i][1] = g[1];
  sp[i][2] = g[2];			  
			  
  //Renormalization (may not be necessary)
  msq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
  scale = 1.0/sqrt(msq);
  sp[i][0] *= scale;
  sp[i][1] *= scale;
  sp[i][2] *= scale;
 
}

/* ---------------------------------------------------------------------- */

void FixNVESpin::final_integrate()
{	
  double dtfm,msq,scale,fm2,fmsq,energy;
  double cp[3],g[3]; 	
	
  double **x = atom->x;	
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;  
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;  
  int *type = atom->type;
  int *mask = atom->mask; 

  // update half x for all particles
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i][0] += 0.5 * dtv * v[i][0];
      x[i][1] += 0.5 * dtv * v[i][1];
      x[i][2] += 0.5 * dtv * v[i][2];
      }
  }   	
  
  // update half v for all particles
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass) dtfm = dtf / rmass[i];
      else dtfm = dtf / mass[type[i]]; 
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];  
      }
  }
}
