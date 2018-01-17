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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "fix_force_spin.h"
#include "fix_integration_spin.h"
#include "fix_langevin_spin.h"
#include "force.h"
#include "math_vector.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h" 
#include "neighbor.h"
#include "neigh_list.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "pair_spin_exchange.h"
#include "pair_spin_me.h"
#include "pair_spin_soc_dmi.h"
#include "pair_spin_soc_neel.h"
#include "pair_spin_soc_landau.h"
#include "respa.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace MathExtra;

enum{NONE,SPIN};

/* ---------------------------------------------------------------------- */

FixIntegrationSpin::FixIntegrationSpin(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), 
  rsec(NULL), stack_head(NULL), stack_foot(NULL), 
  backward_stacks(NULL), forward_stacks(NULL),
  lockpairspinexchange(NULL), lockpairspinsocneel(NULL), lockforcespin(NULL),
  locklangevinspin(NULL)
{
  
  if (narg < 4) error->all(FLERR,"Illegal fix/integration/spin command");  
  
  time_integrate = 1;
  
  extra = NONE;
  mpi_flag = NONE;
  mech_flag = 1;
 
  if (strcmp(arg[2],"integration/spin") == 0) {
    extra = SPIN;
  } else error->all(FLERR,"Illegal fix integration/spin command");
  
  int iarg = 3;

  while (iarg < narg) { 
    if (strcmp(arg[iarg],"serial") == 0){
      mpi_flag = 0;
      iarg += 1;
    } else if (strcmp(arg[iarg],"mpi") == 0) {
      mpi_flag = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"lattice") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix/integration/spin command");
      if (strcmp(arg[iarg+1],"no") == 0) mech_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) mech_flag = 1;
      else error->all(FLERR,"Illegal fix/integration/spin command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix langevin command");
  }

  if (extra == SPIN && !atom->mumag_flag)
    error->all(FLERR,"Fix integration/spin requires spin attribute mumag");

  if (mpi_flag == NONE)
    error->all(FLERR,"Illegal fix/integration/spin command");

  magpair_flag = 0;
  exch_flag = 0;
  soc_neel_flag = soc_dmi_flag = 0;
  me_flag = 0;
  magforce_flag = 0;
  zeeman_flag = aniso_flag = 0;
  maglangevin_flag = 0;
  tdamp_flag = temp_flag = 0;

}

/* ---------------------------------------------------------------------- */

FixIntegrationSpin::~FixIntegrationSpin()
{
  memory->destroy(rsec);
  memory->destroy(stack_head);
  memory->destroy(stack_foot);
  memory->destroy(forward_stacks);
  memory->destroy(backward_stacks);
}

/* ---------------------------------------------------------------------- */

int FixIntegrationSpin::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_NEIGHBOR;
  mask |= FINAL_INTEGRATE;
  return mask;  
}

/* ---------------------------------------------------------------------- */

void FixIntegrationSpin::init()
{
  // set timesteps

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dts = 0.25 * update->dt;

  // set necessary pointers to interaction classes

  if (strstr(force->pair_style,"pair/spin/exchange")) {
    exch_flag = 1;
    lockpairspinexchange = (PairSpinExchange *) force->pair;
  } else if (strstr(force->pair_style,"pair/spin/soc/neel")) {
    soc_neel_flag = 1;
    lockpairspinsocneel = (PairSpinSocNeel *) force->pair;
  } else if (strstr(force->pair_style,"pair/spin/soc/dmi")) {
    soc_dmi_flag = 1;
    lockpairspinsocdmi = (PairSpinSocDmi *) force->pair;
  } else if (strstr(force->pair_style,"pair/spin/me")) {
    me_flag = 1;
    lockpairspinme = (PairSpinMe *) force->pair;
  } else if (strstr(force->pair_style,"hybrid/overlay")) {
    PairHybrid *lockhybrid = (PairHybrid *) force->pair;
    int nhybrid_styles = lockhybrid->nstyles;
    int ipair;
    for (ipair = 0; ipair < nhybrid_styles; ipair++) {
      if (strstr(lockhybrid->keywords[ipair],"pair/spin/exchange")) {
        exch_flag = 1;
	lockpairspinexchange = (PairSpinExchange *) lockhybrid->styles[ipair];
      } else if (strstr(lockhybrid->keywords[ipair],"pair/spin/soc/neel")) {
	soc_neel_flag = 1;
	lockpairspinsocneel = (PairSpinSocNeel *) lockhybrid->styles[ipair];
      } else if (strstr(lockhybrid->keywords[ipair],"pair/spin/soc/dmi")) {
	soc_dmi_flag = 1;
	lockpairspinsocdmi = (PairSpinSocDmi *) lockhybrid->styles[ipair];
      } else if (strstr(lockhybrid->keywords[ipair],"pair/spin/me")) {
	me_flag = 1;
	lockpairspinme = (PairSpinMe *) lockhybrid->styles[ipair];
      }
    }
  } else error->all(FLERR,"Illegal fix integration/spin command");

  // check errors, and handle simple hybrid (not overlay), and no pair/spin interaction

  int iforce;
  for (iforce = 0; iforce < modify->nfix; iforce++) {
    if (strstr(modify->fix[iforce]->style,"force/spin")) {
      magforce_flag = 1;
      lockforcespin = (FixForceSpin *) modify->fix[iforce]; 
    }
  }

  for (iforce = 0; iforce < modify->nfix; iforce++) {
    if (strstr(modify->fix[iforce]->style,"langevin/spin")) {
      maglangevin_flag = 1;
      locklangevinspin = (FixLangevinSpin *) modify->fix[iforce]; 
    } 
  }

  if (magforce_flag) { 
    if (lockforcespin->zeeman_flag == 1) zeeman_flag = 1;
    if (lockforcespin->aniso_flag == 1) aniso_flag = 1;
  }

  if (maglangevin_flag) {
   if (locklangevinspin->tdamp_flag == 1) tdamp_flag = 1;
   if (locklangevinspin->temp_flag == 1) temp_flag = 1;
  }
  

  memory->create(rsec,3,"integration/spin:rsec");
  
  // perform the sectoring if mpi integration

  if (mpi_flag) {
    sectoring();
    
    // grow tables of stacking variables
    
    stack_head = memory->grow(stack_head,nsectors,"integration/spin:stack_head");
    stack_foot = memory->grow(stack_foot,nsectors,"integration/spin:stack_foot");
    forward_stacks = memory->grow(forward_stacks,atom->nmax,"integration/spin:forward_stacks");
    backward_stacks = memory->grow(backward_stacks,atom->nmax,"integration/spin:backward_stacks");
  }

  // grow tables of stacking variables
 /* 
  stack_head = memory->grow(stack_head,nsectors,"integration/spin:stack_head");
  stack_foot = memory->grow(stack_foot,nsectors,"integration/spin:stack_foot");
  forward_stacks = memory->grow(forward_stacks,atom->nmax,"integration/spin:forward_stacks");
  backward_stacks = memory->grow(backward_stacks,atom->nmax,"integration/spin:backward_stacks");
*/
}

/* ---------------------------------------------------------------------- */

void FixIntegrationSpin::initial_integrate(int vflag)
{
  double dtfm,msq,scale,fm2,fmsq,sp2,spsq,energy;
  double spi[3], fmi[3];
	
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

  // update half v for all atoms
  
  if (mech_flag) {
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

  // update half s for all atoms

  if (extra == SPIN) {
    if (mpi_flag) {				// mpi seq. update
      for (int j = 0; j < nsectors; j++) {	// advance quarter s for nlocal
        comm->forward_comm();
	int i = stack_foot[j];
        while (i >= 0) {
	  ComputeInteractionsSpin(i);
    	  AdvanceSingleSpin(i,dts);
	  i = forward_stacks[i];
	}
      }
      for (int j = nsectors-1; j >= 0; j--) {	// advance quarter s for nlocal 
        comm->forward_comm();
	int i = stack_head[j];
	while (i >= 0) {
	  ComputeInteractionsSpin(i);
    	  AdvanceSingleSpin(i,dts);
	  i = backward_stacks[i];
	}
      }
    } else if (mpi_flag == 0) {			// serial seq. update
      comm->forward_comm();			// comm. positions of ghost atoms
      for (int i = 0; i < nlocal-1; i++){	// advance quarter s for nlocal-1
        ComputeInteractionsSpin(i);
    	AdvanceSingleSpin(i,dts);
      }
      ComputeInteractionsSpin(nlocal-1);
      AdvanceSingleSpin(nlocal-1,2.0*dts);	// advance half s for 1
      for (int i = nlocal-2; i >= 0; i--){	// advance quarter s for nlocal-1
        ComputeInteractionsSpin(i);
        AdvanceSingleSpin(i,dts);
      }
    } else error->all(FLERR,"Illegal fix integration/spin command");
  }


  // update x for all particles

  if (mech_flag) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
    }
  }

  // update half s for all particles

  if (extra == SPIN) {
    if (mpi_flag) {				// mpi seq. update
      for (int j = 0; j < nsectors; j++) {	// advance quarter s for nlocal
        comm->forward_comm();
	int i = stack_foot[j];
        while (i >= 0) {
	  ComputeInteractionsSpin(i);
    	  AdvanceSingleSpin(i,dts);
	  i = forward_stacks[i];
	}
      }
      for (int j = nsectors-1; j >= 0; j--) {	// advance quarter s for nlocal 
        comm->forward_comm();
	int i = stack_head[j];
	while (i >= 0) {
	  ComputeInteractionsSpin(i);
    	  AdvanceSingleSpin(i,dts);
	  i = backward_stacks[i];
	}
      }
    } else if (mpi_flag == 0) {			// serial seq. update
      comm->forward_comm();			// comm. positions of ghost atoms
      for (int i = 0; i < nlocal-1; i++){	// advance quarter s for nlocal-1
        ComputeInteractionsSpin(i);
        AdvanceSingleSpin(i,dts);
      }
      ComputeInteractionsSpin(nlocal-1);
      AdvanceSingleSpin(nlocal-1,2.0*dts);	// advance half s for 1
      for (int i = nlocal-2; i >= 0; i--){	// advance quarter s for nlocal-1
        ComputeInteractionsSpin(i);
        AdvanceSingleSpin(i,dts);
      }
    } else error->all(FLERR,"Illegal fix integration/spin command");
  }

}

/* ---------------------------------------------------------------------- 
   setup pre_neighbor()
---------------------------------------------------------------------- */

void FixIntegrationSpin::setup_pre_neighbor()
{
  pre_neighbor();
}

/* ---------------------------------------------------------------------- 
   store in two linked lists the advance order of the spins
---------------------------------------------------------------------- */

void FixIntegrationSpin::pre_neighbor()
{
  double **x = atom->x;	
  int nlocal = atom->nlocal;

  for (int j = 0; j < nsectors; j++) {
    stack_head[j] = -1;
    stack_foot[j] = -1;
  }

  int nseci;
  for (int j = 0; j < nsectors; j++) {		// stacking backward order
    for (int i = 0; i < nlocal; i++) {
      nseci = coords2sector(x[i]);
      if (j != nseci) continue;
      backward_stacks[i] = stack_head[j];
      stack_head[j] = i;
    }
  }
  for (int j = nsectors-1; j >= 0; j--) {	// stacking forward order
    for (int i = nlocal-1; i >= 0; i--) {
      nseci = coords2sector(x[i]);
      if (j != nseci) continue;
      forward_stacks[i] = stack_foot[j];
      stack_foot[j] = i;
    }
  }

}

/* ---------------------------------------------------------------------- 
   compute the magnetic force for the spin ii
---------------------------------------------------------------------- */

void FixIntegrationSpin::ComputeInteractionsSpin(int ii)
{
  const int nlocal = atom->nlocal;
  double xi[3], rij[3], eij[3];
  double spi[3], spj[3];
  double fmi[3];

  int i,j,jj,inum,jnum,itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double **x = atom->x;
  double **sp = atom->sp;
  double **fm = atom->fm;
  int *type = atom->type;
  const int newton_pair = force->newton_pair;

  // add test here (only exchange quantities here)
  if (exch_flag) { 
    inum = lockpairspinexchange->list->inum;
    ilist = lockpairspinexchange->list->ilist;
    numneigh = lockpairspinexchange->list->numneigh;
    firstneigh = lockpairspinexchange->list->firstneigh;
  }
 
  double rsq, rd;
  double temp_cut, cut_2, inorm;
  temp_cut = cut_2 = inorm = 0.0;  

  int eflag = 1;
  int vflag = 0;
  int pair_compute_flag = 1;

  // force computation for spin i

  i = ilist[ii];
  
  spi[0] = sp[i][0];
  spi[1] = sp[i][1];
  spi[2] = sp[i][2];
 
  xi[0] = x[i][0];
  xi[1] = x[i][1];
  xi[2] = x[i][2];
  fmi[0] = fmi[1] = fmi[2] = 0.0;
  jlist = firstneigh[i];
  jnum = numneigh[i];

  // loop on neighbors for pair interactions

  for (int jj = 0; jj < jnum; jj++) {

    j = jlist[jj];
    j &= NEIGHMASK;
    itype = type[ii];
    jtype = type[j];

    spj[0] = sp[j][0];
    spj[1] = sp[j][1];
    spj[2] = sp[j][2];

    rij[0] = x[j][0] - xi[0];
    rij[1] = x[j][1] - xi[1];
    rij[2] = x[j][2] - xi[2];
    rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
    inorm = 1.0/sqrt(rsq);

    temp_cut = 0.0;

    if (exch_flag) {		// exchange
      temp_cut = lockpairspinexchange->cut_spin_exchange[itype][jtype];
      cut_2 = temp_cut*temp_cut;
      if (rsq <= cut_2) {
	lockpairspinexchange->compute_exchange(i,j,rsq,fmi,spi,spj);
      }  
    }

    if (soc_neel_flag) {		// soc_neel
      temp_cut = lockpairspinsocneel->cut_soc_neel[itype][jtype];
      cut_2 = temp_cut*temp_cut;
      if (rsq <= cut_2) {
	eij[0] = rij[0]*inorm;
	eij[1] = rij[1]*inorm;
	eij[2] = rij[2]*inorm;
	lockpairspinsocneel->compute_soc_neel(i,j,rsq,eij,fmi,spi,spj);
      }
    }
    
    if (soc_dmi_flag) {		// soc_dmi
      temp_cut = lockpairspinsocdmi->cut_soc_dmi[itype][jtype];
      cut_2 = temp_cut*temp_cut;
      if (rsq <= cut_2) {
	eij[0] = rij[0]*inorm;
	eij[1] = rij[1]*inorm;
	eij[2] = rij[2]*inorm;
	lockpairspinsocdmi->compute_soc_dmi(i,j,fmi,spi,spj);
      }
    }

    if (me_flag) {		// me
      temp_cut = lockpairspinme->cut_spin_me[itype][jtype];
      cut_2 = temp_cut*temp_cut;
      if (rsq <= cut_2) {
	eij[0] = rij[0]*inorm;
	eij[1] = rij[1]*inorm;
	eij[2] = rij[2]*inorm;
	lockpairspinme->compute_me(i,j,eij,fmi,spi,spj);
      }
    }

  }  
  
  if (magforce_flag) {		// mag. forces
    if (zeeman_flag) {		// zeeman
      lockforcespin->compute_zeeman(i,fmi);
    }
    if (aniso_flag) {		// aniso
      lockforcespin->compute_anisotropy(i,spi,fmi);
    }
  }

  if (maglangevin_flag) {	// mag. langevin
    if (tdamp_flag) {		// transverse damping
      locklangevinspin->add_tdamping(spi,fmi);   
    }
    if (temp_flag) { 		// spin temperature
      locklangevinspin->add_temperature(fmi);
    } 
  }

  // replace the magnetic force fm[i] by its new value

  fm[i][0] = fmi[0];
  fm[i][1] = fmi[1];
  fm[i][2] = fmi[2];

}

/* ---------------------------------------------------------------------- 
   divide each domain into sectors
---------------------------------------------------------------------- */

void FixIntegrationSpin::sectoring()
{
  int sec[3];
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

  const double rv = lockpairspinexchange->cut_spin_exchange_global;

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

  if (mpi_flag == 1 && nsectors != 8) 
    error->all(FLERR,"Illegal sectoring operation");

  rsec[0] = rsx;
  rsec[1] = rsy;
  rsec[2] = rsz;
  if (sec[0] == 2) rsec[0] = rsx/2.0;
  if (sec[1] == 2) rsec[1] = rsy/2.0;
  if (sec[2] == 2) rsec[2] = rsz/2.0;

}

/* ---------------------------------------------------------------------- 
   define sector for an atom at a position x[i]
---------------------------------------------------------------------- */

int FixIntegrationSpin::coords2sector(double *x)
{
  int nseci;
  int seci[3];
  double sublo[3];
  double* sublotmp = domain->sublo;
  for (int dim = 0 ; dim<3 ; dim++) {
    sublo[dim]=sublotmp[dim];
  }

//#define M1
#if defined M1
  double rix = (x[0] - sublo[0])/rsec[0];
  double riy = (x[1] - sublo[1])/rsec[1];
  double riz = (x[2] - sublo[2])/rsec[2];

  seci[0] = static_cast<int>(rix);
  seci[1] = static_cast<int>(riy);
  seci[2] = static_cast<int>(riz);
#endif

#define M2
#if defined M2
  seci[0] = x[0] > (sublo[0] + rsec[0]);
  seci[1] = x[1] > (sublo[1] + rsec[1]);
  seci[2] = x[2] > (sublo[2] + rsec[2]);
#endif

  nseci = (seci[0] + 2*seci[1] + 4*seci[2]); 

  return nseci;
}

/* ---------------------------------------------------------------------- 
   advance the spin i of a timestep dtl
---------------------------------------------------------------------- */

void FixIntegrationSpin::AdvanceSingleSpin(int i, double dtl)
{
  int j=0;
  int *sametag = atom->sametag;
  double **sp = atom->sp;
  double **fm = atom->fm;
  double dtfm,msq,scale,fm2,fmsq,sp2,spsq,energy,dts2;
  double cp[3],g[3]; 	

  cp[0] = cp[1] = cp[2] = 0.0;
  g[0] = g[1] = g[2] = 0.0;
  fm2 = (fm[i][0]*fm[i][0])+(fm[i][1]*fm[i][1])+(fm[i][2]*fm[i][2]);
  fmsq = sqrt(fm2);
  energy = (sp[i][0]*fm[i][0])+(sp[i][1]*fm[i][1])+(sp[i][2]*fm[i][2]);
  dts2 = dtl*dtl;		

  cp[0] = fm[i][1]*sp[i][2]-fm[i][2]*sp[i][1];
  cp[1] = fm[i][2]*sp[i][0]-fm[i][0]*sp[i][2];
  cp[2] = fm[i][0]*sp[i][1]-fm[i][1]*sp[i][0];

  g[0] = sp[i][0]+cp[0]*dtl;
  g[1] = sp[i][1]+cp[1]*dtl;
  g[2] = sp[i][2]+cp[2]*dtl;
			  
  g[0] += (fm[i][0]*energy-0.5*sp[i][0]*fm2)*0.5*dts2;
  g[1] += (fm[i][1]*energy-0.5*sp[i][1]*fm2)*0.5*dts2;
  g[2] += (fm[i][2]*energy-0.5*sp[i][2]*fm2)*0.5*dts2;
			  
  g[0] /= (1+0.25*fm2*dts2);
  g[1] /= (1+0.25*fm2*dts2);
  g[2] /= (1+0.25*fm2*dts2);
			  
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

  if (mpi_flag == 0) {
    if (sametag[i] >= 0) {
      j = sametag[i];
      while (j >= 0) {
        sp[j][0] = sp[i][0];
        sp[j][1] = sp[i][1];
        sp[j][2] = sp[i][2];
        j = sametag[j];
      }
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixIntegrationSpin::final_integrate()
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

  if (mech_flag) {
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

}
