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

  // update sp for all particles
  if (extra == SPIN) {
      // Advance spins
      //See Omelyan et al., PRL 86, 2001 and P.W. Ma et al, PRB 83, 2011
      for (int i = 0; i < nlocal; i++)
	      if (mask[i] & groupbit) {
                          cp[0] = cp[1] = cp[2] = 0.0;
                          g[0] = g[1] = g[2] = 0.0;
			  fm2 = (fm[i][0]*fm[i][0])+(fm[i][1]*fm[i][1])+(fm[i][2]*fm[i][2]);
			  fmsq = sqrt(fm2);
			  energy = (sp[i][0]*fm[i][0])+(sp[i][1]*fm[i][1])+(sp[i][2]*fm[i][2]);
			    
			  cp[0] = fm[i][1]*sp[i][2]-fm[i][2]*sp[i][1];
			  cp[1] = fm[i][2]*sp[i][0]-fm[i][0]*sp[i][2];
			  cp[2] = fm[i][0]*sp[i][1]-fm[i][1]*sp[i][0];
			  
			  g[0] = sp[i][0]+cp[0]*dts;
			  g[1] = sp[i][1]+cp[1]*dts;
			  g[2] = sp[i][2]+cp[2]*dts;
			  
			  g[0] += (fm[i][0]*energy-0.5*sp[i][0]*fm2)*0.5*dts*dts;
			  g[1] += (fm[i][1]*energy-0.5*sp[i][1]*fm2)*0.5*dts*dts;
			  g[2] += (fm[i][2]*energy-0.5*sp[i][2]*fm2)*0.5*dts*dts;
			  
			  g[0] /= (1+0.25*fm2*dts*dts);
			  g[1] /= (1+0.25*fm2*dts*dts);
			  g[2] /= (1+0.25*fm2*dts*dts);
			  
			  sp[i][0] = g[0];
			  sp[i][1] = g[1];
			  sp[i][2] = g[2];			  
			  
			  //Renormalization (may not be necessary)
                          msq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
                          scale = 1.0/sqrt(msq);
                          sp[i][0] *= scale;
                          sp[i][1] *= scale;
                          sp[i][2] *= scale;
                          
                          //printf("test fix integ. 1;i=%d, fx=%g, fy=%g, fz=%g \n",i,fm[i][0],fm[i][1],fm[i][2]);	  

                          //printf("test fix integ.; i=%d, sx=%g, sy=%g, sz=%g, norm=%g \n",i,sp[i][0],sp[i][1],sp[i][2],scale);	  
		  }
	  }

                          //printf("test fix integ. 1;i=0, fx=%g, fy=%g, fz=%g \n",fm[0][0],fm[0][1],fm[0][2]);	  
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
