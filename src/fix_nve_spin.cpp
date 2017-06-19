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
#include "pair.h"
#include "timer.h"
#include "integrate.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "pair_spin.h"
#include "memory.h"
#include "fix_force_spin.h"

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

#if defined SEQNEI
  lockpairspin = NULL;
  lockforcespin = NULL;
  exch_flag = dmi_flag = me_flag = 0;
  zeeman_flag = aniso_flag = 0;
#endif 
}

/* ---------------------------------------------------------------------- */
FixNVESpin::~FixNVESpin(){
#if defined SEQNEI
  //delete lockpairspin;
  //delete lockforcespin;
  memory->destroy(spi);
  memory->destroy(fmi);
  memory->destroy(fmj);
#endif
}

/* ---------------------------------------------------------------------- */

void FixNVESpin::init()
{
  FixNVE::init();       
  
  dts = update->dt;

  #if defined SEQNEI
  lockpairspin = (PairSpin *) force->pair;

  memory->create(spi,3,"nves:spi");
  memory->create(fmi,3,"nves:fmi");
  memory->create(fmj,3,"nves:fmj");

  int iforce;
  for (iforce = 0; iforce < modify->nfix; iforce++)
    if (strstr(modify->fix[iforce]->style,"force/spin")) break;
  lockforcespin = (FixForceSpin *) modify->fix[iforce]; 

  exch_flag = lockpairspin->exch_flag;
  dmi_flag = lockpairspin->dmi_flag;
  me_flag = lockpairspin->me_flag; 

  zeeman_flag = lockforcespin->zeeman_flag;
  aniso_flag = lockforcespin->aniso_flag;
  #endif
  
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
  double **sp = atom->sp;
  double **fm = atom->fm;
  double *rmass = atom->rmass; 
  double *mass = atom->mass;  
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;  
  int *type = atom->type;
  int *mask = atom->mask;  

  // Advance half spins all particles
  //See Omelyan et al., PRL 86, 2001 and P.W. Ma et al, PRB 83, 2011
  if (extra == SPIN) {
#if defined SEQNEI
    for (int i = 0; i < nlocal; i++){
      ComputeSpinInteractionsNei(i);
      AdvanceSingleSpin(i,0.5*dts,sp,fm);
    }
#endif
#if defined SEQ
    for (int i = 0; i < nlocal; i++){
      AdvanceSingleSpin(i,0.5*dts,sp,fm);
      ComputeSpinInteractions();
    }
#endif 
  }
  
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

  // update x for all particles
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i][0] += 0.5 * dtv * v[i][0];
      x[i][1] += 0.5 * dtv * v[i][1];
      x[i][2] += 0.5 * dtv * v[i][2];
      }
  }  
}

#if defined SEQNEI
/* ---------------------------------------------------------------------- */
void FixNVESpin::ComputeSpinInteractionsNei(int ii)
{
  int nflag,sortflag;
  int nlocal = atom->nlocal;

  int n_post_integrate = modify->n_post_integrate;
  int n_pre_exchange = modify->n_pre_exchange;
  int n_pre_neighbor = modify->n_pre_neighbor;
  int n_pre_force = modify->n_pre_force;
  int n_pre_reverse = modify->n_pre_reverse;
  int n_post_force = modify->n_post_force;
  int n_end_of_step = modify->n_end_of_step;

  bigint ntimestep;
  ntimestep = update->ntimestep;

  //Force compute quantities
  int i,j,jj,inum,jnum,itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double **x = atom->x;
  double **sp = atom->sp;
  double **fm = atom->fm;
  int *type = atom->type;
  int newton_pair = force->newton_pair;

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

  if (atom->sortfreq > 0) sortflag = 1;
  else sortflag = 0;
  if (n_post_integrate) modify->post_integrate();
  timer->stamp(Timer::MODIFY);

  // regular communication vs neighbor list rebuild
  nflag = neighbor->decide();
  if (nflag == 0) {
  timer->stamp();
  comm->forward_comm();
  timer->stamp(Timer::COMM);
  } else {
    if (n_pre_exchange) {
      timer->stamp();
      modify->pre_exchange();
      timer->stamp(Timer::MODIFY);
    }
    //if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    if (domain->box_change) {
      domain->reset_box();
      comm->setup();
      if (neighbor->style) neighbor->setup_bins();
    }
    timer->stamp();
    comm->exchange();
    if (sortflag && ntimestep >= atom->nextsort) atom->sort();
    comm->borders();
    //if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    timer->stamp(Timer::COMM);
    if (n_pre_neighbor) {
      modify->pre_neighbor();
      timer->stamp(Timer::MODIFY);
    }
    neighbor->build();
    timer->stamp(Timer::NEIGH);
  }

  ///////Force computation for spin i/////////////
  i = ilist[ii];
  //Clear atom i 
  fm[i][0] = fm[i][1] = fm[i][2] = 0.0;

  timer->stamp();

  xtmp = x[i][0];
  ytmp = x[i][1];
  ztmp = x[i][2];
  fmi[0] = fmi[1] = fmi[2] = 0.0;
  fmj[0] = fmj[1] = fmj[2] = 0.0;
  jlist = firstneigh[i];
  jnum = numneigh[i];

//  printf("Test inum: %g \n",inum);
/*
  //Pair interaction
  for (int jj = 0; jj < inum; jj++) {
    j = jlist[jj];
    j &= NEIGHMASK;

    delx = xtmp - x[j][0];
    dely = ytmp - x[j][1];
    delz = ztmp - x[j][2];
    rsq = delx*delx + dely*dely + delz*delz;
    itype = type[ii];
    jtype = type[j];

    if (exch_flag) {
      cut_ex_2 = (lockpairspin->cut_spin_exchange[itype][jtype])*(lockpairspin->cut_spin_exchange[itype][jtype]);
      if (rsq <= cut_ex_2) {
        lockpairspin->compute_exchange(i,j,rsq,fmi,fmj);
      }  
    }
    if (dmi_flag) {
      cut_dmi_2 = (lockpairspin->cut_spin_dmi[itype][jtype])*(lockpairspin->cut_spin_dmi[itype][jtype]);
      if (rsq <= cut_dmi_2) {
        lockpairspin->compute_dmi(i,j,fmi,fmj);
      }  
    }
    if (me_flag) {
      cut_me_2 = (lockpairspin->cut_spin_me[itype][jtype])*(lockpairspin->cut_spin_me[itype][jtype]);
      if (rsq <= cut_me_2) {
        lockpairspin->compute_me(i,j,fmi,fmj);
      }  
    }
  }
*/
  
  //Pair interaction
  int natom = nlocal + atom->nghost;
  for (int k = 0; k < natom; k++) {
    delx = xtmp - x[k][0];
    dely = ytmp - x[k][1];
    delz = ztmp - x[k][2];
    rsq = delx*delx + dely*dely + delz*delz;
    itype = type[ii];
    jtype = type[k];

    if (exch_flag) {
      cut_ex_2 = (lockpairspin->cut_spin_exchange[itype][jtype])*(lockpairspin->cut_spin_exchange[itype][jtype]);
      if (rsq <= cut_ex_2) {
        lockpairspin->compute_exchange(i,k,rsq,fmi,fmj);
      }  
    }
    if (dmi_flag) {
      cut_dmi_2 = (lockpairspin->cut_spin_dmi[itype][jtype])*(lockpairspin->cut_spin_dmi[itype][jtype]);
      if (rsq <= cut_dmi_2) {
        lockpairspin->compute_dmi(i,k,fmi,fmj);
      }  
    }
    if (me_flag) {
      cut_me_2 = (lockpairspin->cut_spin_me[itype][jtype])*(lockpairspin->cut_spin_me[itype][jtype]);
      if (rsq <= cut_me_2) {
        lockpairspin->compute_me(i,k,fmi,fmj);
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
    
  //Replace the force by its new value
  fm[i][0] = fmi[0];
  fm[i][1] = fmi[1];
  fm[i][2] = fmi[2];

}
#endif

/* ---------------------------------------------------------------------- */
void FixNVESpin::ComputeSpinInteractions()
{
  int nflag,sortflag;
  int nlocal = atom->nlocal;

  int n_post_integrate = modify->n_post_integrate;
  int n_pre_exchange = modify->n_pre_exchange;
  int n_pre_neighbor = modify->n_pre_neighbor;
  int n_pre_force = modify->n_pre_force;
  int n_pre_reverse = modify->n_pre_reverse;
  int n_post_force = modify->n_post_force;
  int n_end_of_step = modify->n_end_of_step;

  bigint ntimestep;
  ntimestep = update->ntimestep; 

  //int eflag = update->integrate->eflag;
  //int vflag = update->integrate->vflag;
  int eflag = 1;
  int vflag = 0;
  int pair_compute_flag = 1;

  if (atom->sortfreq > 0) sortflag = 1;
  else sortflag = 0;
  if (n_post_integrate) modify->post_integrate();
  timer->stamp(Timer::MODIFY);

  // regular communication vs neighbor list rebuild
  nflag = neighbor->decide();
  if (nflag == 0) {
  timer->stamp();
  comm->forward_comm();
  timer->stamp(Timer::COMM);
  } else {
    if (n_pre_exchange) {
      timer->stamp();
      modify->pre_exchange();
      timer->stamp(Timer::MODIFY);
    }
    //if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    if (domain->box_change) {
      domain->reset_box();
      comm->setup();
      if (neighbor->style) neighbor->setup_bins();
    }
    timer->stamp();
    comm->exchange();
    if (sortflag && ntimestep >= atom->nextsort) atom->sort();
    comm->borders();
    //if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    timer->stamp(Timer::COMM);
    if (n_pre_neighbor) {
      modify->pre_neighbor();
      timer->stamp(Timer::MODIFY);
    }
    neighbor->build();
    timer->stamp(Timer::NEIGH);
  }

  // force computations
  // important for pair to come before bonded contributions
  // since some bonded potentials tally pairwise energy/virial
  // and Pair:ev_tally() needs to be called before any tallying

  size_t nbytes;
  nbytes = sizeof(double) * nlocal;
  if (force->newton) nbytes += sizeof(double) * atom->nghost;

  atom->avec->force_clear(0,nbytes); 

  timer->stamp();

  if (n_pre_force) {
    modify->pre_force(vflag);
    timer->stamp(Timer::MODIFY);
  }

  if (pair_compute_flag) {
    force->pair->compute(eflag,vflag);
    timer->stamp(Timer::PAIR);
 }
 
  /*if (kspace_compute_flag) {
    force->kspace->compute(eflag,vflag);
    timer->stamp(Timer::KSPACE);
  }*/
 
  if (n_pre_reverse) {
    modify->pre_reverse(eflag,vflag);
    timer->stamp(Timer::MODIFY);
  }
  
  // reverse communication of forces
 
  if (force->newton) {
    comm->reverse_comm();
    timer->stamp(Timer::COMM);
  }
 
  // force modifications
 
  if (n_post_force) modify->post_force(vflag);
  timer->stamp(Timer::MODIFY);

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
  double **sp = atom->sp;
  double **fm = atom->fm;
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

  // Advance half spins all particles
  //See Omelyan et al., PRL 86, 2001 and P.W. Ma et al, PRB 83, 2011
  if (extra == SPIN) {
#if defined SEQNEI
    for (int i = nlocal-1; i >= 0; i--){
      ComputeSpinInteractionsNei(i);
      AdvanceSingleSpin(i,0.5*dts,sp,fm);
    }
#endif
#if defined SEQ
    for (int i = nlocal-1; i >= 0; i--){
      AdvanceSingleSpin(i,0.5*dts,sp,fm);
      ComputeSpinInteractions();
    }
#endif 
  }
}
