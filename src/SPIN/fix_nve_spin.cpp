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

//Add headers (see delete later)
#include "neighbor.h"
#include "neigh_list.h"
#include "pair.h"
#include "pair_spin.h"
#include "memory.h"
#include "fix_force_spin.h"
#include "fix_langevin_spin.h"

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

  exch_flag = dmi_flag = me_flag = 0;
  zeeman_flag = aniso_flag = 0;
  tdamp_flag = temp_flag = 0;

  lockpairspin = NULL;
  lockforcespin = NULL;
  locklangevinspin = NULL;
}

/* ---------------------------------------------------------------------- */
FixNVESpin::~FixNVESpin(){
  //delete lockpairspin;
  //delete lockforcespin;
  memory->destroy(xi);
#if defined SECTORING  
  memory->destroy(sec);
  memory->destroy(rsec);
  memory->destroy(seci);
#endif
#if defined SECTOR_PRINT
  fclose(file_sect);
#endif
  memory->destroy(spi);
  memory->destroy(spj);
  memory->destroy(fmi);
  memory->destroy(fmj);
}

/* ---------------------------------------------------------------------- */

void FixNVESpin::init()
{
  FixNVE::init();       
  
  dts = update->dt;
  memory->create(xi,3,"nves:xi");
#if defined SECTORING
  memory->create(sec,3,"nves:sec");
  memory->create(rsec,3,"nves:rsec");
  memory->create(seci,3,"nves:seci");
#endif
  memory->create(spi,3,"nves:spi");
  memory->create(spj,3,"nves:spj");
  memory->create(fmi,3,"nves:fmi");
  memory->create(fmj,3,"nves:fmj");

  lockpairspin = (PairSpin *) force->pair;

  int iforce;
  for (iforce = 0; iforce < modify->nfix; iforce++)
    if (strstr(modify->fix[iforce]->style,"force/spin")) break;
  lockforcespin = (FixForceSpin *) modify->fix[iforce]; 

  for (iforce = 0; iforce < modify->nfix; iforce++)
    if (strstr(modify->fix[iforce]->style,"langevin/spin")) break;
  locklangevinspin = (FixLangevinSpin *) modify->fix[iforce]; 

  exch_flag = lockpairspin->exch_flag;
  dmi_flag = lockpairspin->dmi_flag;
  me_flag = lockpairspin->me_flag; 

  zeeman_flag = lockforcespin->zeeman_flag;
  aniso_flag = lockforcespin->aniso_flag;

  tdamp_flag = locklangevinspin->tdamp_flag;
  temp_flag = locklangevinspin->temp_flag;


#if defined SECTORING 
  sectoring();
#endif

#if defined SECTOR_PRINT
  file_sect=fopen("sectoring.lammpstrj", "w");
  fprintf(file_sect,"ITEM: TIMESTEP\n");
  fprintf(file_sect,"%g\n",0.0);
  fprintf(file_sect,"ITEM: NUMBER OF ATOMS\n");
  //int natoms = atom->natoms;
  int natoms = atom->nlocal;
  fprintf(file_sect,"%d\n",natoms);
  fprintf(file_sect,"ITEM: BOX BOUNDS\n");
  for(int d=0; d<3; d++) fprintf(file_sect,"%lf %lf\n",domain->boxlo[d],domain->boxhi[d]);
  fprintf(file_sect,"ITEM: ATOMS type x y z vx vy vz\n"); 
#endif

}

/* ---------------------------------------------------------------------- */

void FixNVESpin::initial_integrate(int vflag)
{
  double dtfm,msq,scale,fm2,fmsq,sp2,spsq,energy;
  double cp[3],g[3]; 	
	
  double **x = atom->x;	
  double **v = atom->v;
  double **f = atom->f;
  double **sp = atom->sp;
  double **fm = atom->fm;
  double *rmass = atom->rmass; 
  double *mass = atom->mass;  
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

#if defined SECTORING
  int nseci;
  // Seq. update spins for all particles
  if (extra == SPIN) {
    for (int j = 0; j < nsectors; j++) { 
      comm->forward_comm();
      for (int i = 0; i < nlocal; i++) {
        xi[0] = x[i][0];
        xi[1] = x[i][1];
        xi[2] = x[i][2];
        nseci = coords2sector(xi);
	if (j != nseci) continue;
	ComputeSpinInteractionsNeigh(i);
        AdvanceSingleSpin(i,0.5*dts,sp,fm);
      }    
    }
    for (int j = nsectors-1; j >= 0; j--) { 
      comm->forward_comm();
      for (int i = nlocal-1; i >= 0; i--) {
        xi[0] = x[i][0];
        xi[1] = x[i][1];
        xi[2] = x[i][2];
        nseci = coords2sector(xi);
        if (j != nseci) continue;
        ComputeSpinInteractionsNeigh(i);
        AdvanceSingleSpin(i,0.5*dts,sp,fm);
      }    
    }
  }

#else 
  // Seq. update spins for all particles
  if (extra == SPIN) {
    for (int i = 0; i < nlocal; i++){
      ComputeSpinInteractionsNeigh(i);
      AdvanceSingleSpin(i,0.5*dts,sp,fm);
    }

    for (int i = nlocal-1; i >= 0; i--){
      ComputeSpinInteractionsNeigh(i);
      AdvanceSingleSpin(i,0.5*dts,sp,fm);
    }
  }
#endif
  
  // update x for all particles
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i][0] += 0.5 * dtv * v[i][0];
      x[i][1] += 0.5 * dtv * v[i][1];
      x[i][2] += 0.5 * dtv * v[i][2];
      }
  }


#if defined SECTOR_PRINT
  int my_rank;
  MPI_Comm_rank(world, &my_rank);
  if (my_rank == 0) { 
    for (int j = 0; j < nsectors; j++) { 
      for (int i = 0; i < nlocal; i++) {
        xi[0] = x[i][0];
        xi[1] = x[i][1];
        xi[2] = x[i][2];
        nseci = coords2sector(xi);
        if (j != nseci) continue;
        fprintf(file_sect,"%d %lf %lf %lf %lf %lf %lf\n",j,xi[0],xi[1],xi[2],0.0,0.0,1.0);
      }    
    }
  }
#endif

}

/* ---------------------------------------------------------------------- */
void FixNVESpin::ComputeSpinInteractionsNeigh(int ii)
{
  const int nlocal = atom->nlocal;

  //Force compute quantities
  int i,j,jj,inum,jnum,itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double **x = atom->x;
  double **sp = atom->sp;
  double **fm = atom->fm;
  int *type = atom->type;
  const int newton_pair = force->newton_pair;

  inum = lockpairspin->list->inum;
  ilist = lockpairspin->list->ilist;
  numneigh = lockpairspin->list->numneigh;
  firstneigh = lockpairspin->list->firstneigh;
 
  double xtmp,ytmp,ztmp;
  double rsq,rd,delx,dely,delz;
  double cut_ex_2, cut_dmi_2, cut_me_2;
  cut_ex_2 = cut_dmi_2 = cut_me_2 = 0.0;  

  int eflag = 1;
  int vflag = 0;
  int pair_compute_flag = 1;

#if !defined(SECTORING)
  comm->forward_comm();
#endif

  ///////Force computation for spin i/////////////
  i = ilist[ii];
  
  //Clear atom i 
  fm[i][0] = fm[i][1] = fm[i][2] = 0.0;
  
  spi[0] = sp[i][0];
  spi[1] = sp[i][1];
  spi[2] = sp[i][2];
 
  xi[0] = x[i][0];
  xi[1] = x[i][1];
  xi[2] = x[i][2];
  fmi[0] = fmi[1] = fmi[2] = 0.0;
  fmj[0] = fmj[1] = fmj[2] = 0.0;
  jlist = firstneigh[i];
  jnum = numneigh[i];

  //Pair interaction
  for (int jj = 0; jj < jnum; jj++) {
    j = jlist[jj];
    j &= NEIGHMASK;
    spj[0] = sp[j][0];
    spj[1] = sp[j][1];
    spj[2] = sp[j][2];

    delx = xi[0] - x[j][0];
    dely = xi[1] - x[j][1];
    delz = xi[2] - x[j][2];
    rsq = delx*delx + dely*dely + delz*delz;
    itype = type[ii];
    jtype = type[j];

    if (exch_flag) {
      cut_ex_2 = (lockpairspin->cut_spin_exchange[itype][jtype])*(lockpairspin->cut_spin_exchange[itype][jtype]);
      if (rsq <= cut_ex_2) {
        lockpairspin->compute_exchange(i,j,rsq,fmi,fmj,spi,spj);
      }  
    }

    if (dmi_flag) {
      cut_dmi_2 = (lockpairspin->cut_spin_dmi[itype][jtype])*(lockpairspin->cut_spin_dmi[itype][jtype]);
      if (rsq <= cut_dmi_2) {
        lockpairspin->compute_dmi(i,j,fmi,fmj,spi,spj);
      }  
    }

    if (me_flag) {
      cut_me_2 = (lockpairspin->cut_spin_me[itype][jtype])*(lockpairspin->cut_spin_me[itype][jtype]);
      if (rsq <= cut_me_2) {
        lockpairspin->compute_me(i,j,fmi,fmj,spi,spj);
      }  
    }
  }

  //post force
  if (zeeman_flag) {
    lockforcespin->compute_zeeman(i,fmi);
  }

  if (aniso_flag) {
    spi[0] = sp[i][0];
    spi[1] = sp[i][1];                                       
    spi[2] = sp[i][2]; 
    lockforcespin->compute_anisotropy(i,spi,fmi);
  }

  if (tdamp_flag) {
    locklangevinspin->add_tdamping(spi,fmi);   
  }

  if (temp_flag) {
    locklangevinspin->add_temperature(fmi);
  } 
 
  //Replace the force by its new value
  fm[i][0] = fmi[0];
  fm[i][1] = fmi[1];
  fm[i][2] = fmi[2];

}

#if defined SECTORING
/* ---------------------------------------------------------------------- */
void FixNVESpin::sectoring()
{
  double sublo[3],subhi[3];
  double* sublotmp = domain->sublo;
  double* subhitmp = domain->subhi;
  for (int dim = 0 ; dim<3 ; dim++) {
    sublo[dim]=sublotmp[dim];
    subhi[dim]=subhitmp[dim];
  }

  const double rsx = subhi[0] - sublo[0];  
  const double rsy = subhi[1] - sublo[1];  
  const double rsz = subhi[2] - sublo[2];  

  const double rv = lockpairspin->cut_spin_pair_global;

  double rax = rsx/rv;  
  double ray = rsy/rv;  
  double raz = rsz/rv;  
 
  sec[0] = 1;
  sec[1] = 1;
  sec[2] = 1;
  if (rax >= 2.0) sec[0] = 2;
  if (ray >= 2.0) sec[1] = 2;
  if (raz >= 2.0) sec[2] = 2;

  nsectors = sec[0]*sec[1]*sec[2];

  rsec[0] = rsx;
  rsec[1] = rsy;
  rsec[2] = rsz;
  if (sec[0] == 2) rsec[0] = rsx/2.0;
  if (sec[1] == 2) rsec[1] = rsy/2.0;
  if (sec[2] == 2) rsec[2] = rsz/2.0;

  if (2.0 * rv >= rsx && sec[0] >= 2)
    error->all(FLERR,"Illegal number of sectors"); 

  if (2.0 * rv >= rsy && sec[1] >= 2)
    error->all(FLERR,"Illegal number of sectors"); 

  if (2.0 * rv >= rsz && sec[2] >= 2)
    error->all(FLERR,"Illegal number of sectors"); 

}

/* ---------------------------------------------------------------------- */
int FixNVESpin::coords2sector(double *xi)
{
  int nseci;
  double sublo[3];
  double* sublotmp = domain->sublo;
  for (int dim = 0 ; dim<3 ; dim++) {
    sublo[dim]=sublotmp[dim];
  }

  double rix = (xi[0] - sublo[0])/rsec[0];
  double riy = (xi[1] - sublo[1])/rsec[1];
  double riz = (xi[2] - sublo[2])/rsec[2];

  seci[0] = (int)rix;
  seci[1] = (int)riy;
  seci[2] = (int)riz;

  if (nsectors == 1) {
    nseci = 0;
  } else if (nsectors == 2) {
    nseci = seci[0] + seci[1] + seci[2];
  } else if (nsectors == 4) {
    if (sec[1]*sec[2] == 4) { //plane normal to x
      nseci = (seci[1] + 2*seci[2]);
    } else if (sec[0]*sec[2] == 4) { //plane normal to y
      nseci = (seci[0] + 2*seci[2]);
    } else if (sec[0]*sec[1] == 4) { //plane normal to z
      nseci = (seci[0] + 2*seci[1]);
    }
  } else if (nsectors == 8) {
    nseci = (seci[0] + 2*seci[1] + 4*seci[2]); 
  }

  return nseci;
}

#endif

/* ---------------------------------------------------------------------- */

void FixNVESpin::AdvanceSingleSpin(int i, double dts, double **sp, double **fm)
{
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
