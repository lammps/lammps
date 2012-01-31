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

#include <cstdio>
#include <cstring>
#include "modify_cuda.h"
#include "style_compute.h"
#include "style_fix.h"
#include "atom.h"
#include "comm.h"
#include "fix.h"
#include "compute.h"
#include "group.h"
#include "update.h"
#include "domain.h"
#include "cuda.h"
#include "memory.h"
#include "error.h"

#include "cuda_modify_flags.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace FixConstCuda;

#define DELTA 4

#define BIG 1.0e20


/* ---------------------------------------------------------------------- */

ModifyCuda::ModifyCuda(LAMMPS *lmp) : Modify(lmp)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

  n_initial_integrate_cuda = 0;
  n_post_integrate_cuda = 0;
  n_pre_exchange = 0;
  n_pre_neighbor_cuda = 0;
  n_pre_force_cuda = 0;
  n_post_force_cuda = 0;
  n_final_integrate_cuda = 0;
  n_end_of_step_cuda = 0;
  n_thermo_energy_cuda = 0;

  n_initial_integrate_host = 0;
  n_post_integrate_host = 0;
  n_pre_exchange = 0;
  n_pre_neighbor_host = 0;
  n_pre_force_host = 0;
  n_post_force_host = 0;
  n_final_integrate_host = 0;
  n_end_of_step_host = 0;
  n_thermo_energy_host = 0;

  list_initial_integrate_cuda = NULL;
  list_post_integrate_cuda = NULL;
  list_pre_exchange_cuda = NULL;
  list_pre_neighbor_cuda = NULL;
  list_pre_force_cuda = NULL;
  list_post_force_cuda = NULL;
  list_final_integrate_cuda = NULL;
  list_end_of_step_cuda = NULL;
  list_thermo_energy_cuda = NULL;
  end_of_step_every_cuda = NULL;
}

/* ---------------------------------------------------------------------- */

ModifyCuda::~ModifyCuda()
{
  delete [] list_initial_integrate_cuda;
  delete [] list_post_integrate_cuda;
  delete [] list_pre_exchange_cuda;
  delete [] list_pre_neighbor_cuda;
  delete [] list_pre_force_cuda;
  delete [] list_post_force_cuda;
  delete [] list_final_integrate_cuda;
  delete [] list_end_of_step_cuda;
  delete [] list_thermo_energy_cuda;
  delete [] end_of_step_every_cuda;
}

/* ----------------------------------------------------------------------
   initialize all fixes and computes
------------------------------------------------------------------------- */

void ModifyCuda::init()
{
  int i,j;

  // delete storage of restart info since it is not valid after 1st run

  restart_deallocate();

  // create lists of fixes to call at each stage of run

  list_init(INITIAL_INTEGRATE,n_initial_integrate,list_initial_integrate);
  list_init(POST_INTEGRATE,n_post_integrate,list_post_integrate);
  list_init(PRE_EXCHANGE,n_pre_exchange,list_pre_exchange);
  list_init(PRE_NEIGHBOR,n_pre_neighbor,list_pre_neighbor);
  list_init(PRE_FORCE,n_pre_force,list_pre_force);
  list_init(POST_FORCE,n_post_force,list_post_force);
  list_init(FINAL_INTEGRATE,n_final_integrate,list_final_integrate);
  list_init_end_of_step(END_OF_STEP,n_end_of_step,list_end_of_step);
  list_init_thermo_energy(THERMO_ENERGY,n_thermo_energy,list_thermo_energy);

  list_init(INITIAL_INTEGRATE_CUDA, n_initial_integrate_cuda, list_initial_integrate_cuda);
  list_init(POST_INTEGRATE_CUDA, n_post_integrate_cuda, list_post_integrate_cuda);
  list_init(PRE_EXCHANGE_CUDA, n_pre_exchange_cuda, list_pre_exchange_cuda);
  list_init(PRE_NEIGHBOR_CUDA, n_pre_neighbor_cuda, list_pre_neighbor_cuda);
  list_init(PRE_FORCE_CUDA, n_pre_force_cuda, list_pre_force_cuda);
  list_init(POST_FORCE_CUDA, n_post_force_cuda, list_post_force_cuda);
  list_init(FINAL_INTEGRATE_CUDA, n_final_integrate_cuda, list_final_integrate_cuda);
  list_init_end_of_step_cuda(END_OF_STEP_CUDA, n_end_of_step_cuda, list_end_of_step_cuda);
  list_init_thermo_energy(THERMO_ENERGY_CUDA, n_thermo_energy_cuda, list_thermo_energy_cuda);

  n_initial_integrate_host = n_initial_integrate;
  n_post_integrate_host = n_post_integrate;
  n_pre_exchange_host = n_pre_exchange;
  n_pre_neighbor_host = n_pre_neighbor;
  n_pre_force_host = n_pre_force;
  n_post_force_host = n_post_force;
  n_final_integrate_host = n_final_integrate;
  n_end_of_step_host = n_end_of_step;
  n_thermo_energy_host = n_thermo_energy;
  
  n_initial_integrate = n_initial_integrate_cuda+n_initial_integrate_host;
  n_post_integrate = n_post_integrate_cuda+n_post_integrate_host;
  n_pre_exchange = n_pre_exchange_cuda+n_pre_exchange_host;
  n_pre_neighbor = n_pre_neighbor_cuda+n_pre_neighbor_host;
  n_pre_force = n_pre_force_cuda+n_pre_force_host;
  n_post_force = n_post_force_cuda+n_post_force_host;
  n_final_integrate = n_final_integrate_cuda+n_final_integrate_host;
  n_end_of_step = n_end_of_step_cuda+n_end_of_step_host;
  n_thermo_energy = n_thermo_energy_cuda+n_thermo_energy_host;
  
  list_init(INITIAL_INTEGRATE_RESPA,
	    n_initial_integrate_respa,list_initial_integrate_respa);
  list_init(POST_INTEGRATE_RESPA,
	    n_post_integrate_respa,list_post_integrate_respa);
  list_init(POST_FORCE_RESPA,
	    n_post_force_respa,list_post_force_respa);
  list_init(PRE_FORCE_RESPA,
	    n_pre_force_respa,list_pre_force_respa);
  list_init(FINAL_INTEGRATE_RESPA,
	    n_final_integrate_respa,list_final_integrate_respa);

  list_init(MIN_PRE_EXCHANGE,n_min_pre_exchange,list_min_pre_exchange);
  list_init(MIN_POST_FORCE,n_min_post_force,list_min_post_force);
  list_init(MIN_ENERGY,n_min_energy,list_min_energy);

  // init each fix
  // needs to come before compute init
  // this is b/c some computes call fix->dof()
  // FixRigid::dof() depends on its own init having been called

  for (i = 0; i < nfix; i++) fix[i]->init();

  // set global flag if any fix has its restart_pbc flag set

  restart_pbc_any = 0;
  for (i = 0; i < nfix; i++)
    if (fix[i]->restart_pbc) restart_pbc_any = 1;

  // create list of computes that store invocation times

  list_init_compute();

  // init each compute
  // set invoked_scalar,vector,etc to -1 to force new run to re-compute them
  // add initial timestep to all computes that store invocation times
  //   since any of them may be invoked by initial thermo
  // do not clear out invocation times stored within a compute,
  //   b/c some may be holdovers from previous run, like for ave fixes

  for (i = 0; i < ncompute; i++) {
    compute[i]->init();
    compute[i]->invoked_scalar = -1;
    compute[i]->invoked_vector = -1;
    compute[i]->invoked_array = -1;
    compute[i]->invoked_peratom = -1;
    compute[i]->invoked_local = -1;
  }
  addstep_compute_all(update->ntimestep);

  // warn if any particle is time integrated more than once

  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  int *flag = new int[nlocal];
  for (i = 0; i < nlocal; i++) flag[i] = 0;

  int groupbit;
  for (i = 0; i < nfix; i++) {
    if (fix[i]->time_integrate == 0) continue;
    groupbit = fix[i]->groupbit;
    for (j = 0; j < nlocal; j++)
      if (mask[j] & groupbit) flag[j]++;
  }

  int check = 0;
  for (i = 0; i < nlocal; i++)
    if (flag[i] > 1) check = 1;

  delete [] flag;

  int checkall;
  MPI_Allreduce(&check,&checkall,1,MPI_INT,MPI_SUM,world);
  if (comm->me == 0 && checkall)
    error->warning(FLERR,"One or more atoms are time integrated more than once");
}

/* ----------------------------------------------------------------------
   1st half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyCuda::initial_integrate(int vflag)
{
	for(int i = 0; i < n_initial_integrate_cuda; i++)
		fix[list_initial_integrate_cuda[i]]->initial_integrate(vflag);

	if(n_initial_integrate_host != 0)
	{
		cuda->downloadAll(); cuda->oncpu = true;
		for (int i = 0; i < n_initial_integrate_host; i++)
			fix[list_initial_integrate[i]]->initial_integrate(vflag);
		cuda->uploadAll(); cuda->oncpu = false;
	}
}

/* ----------------------------------------------------------------------
   post_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyCuda::post_integrate()
{
	for(int i = 0; i < n_post_integrate_cuda; i++)
		fix[list_post_integrate_cuda[i]]->post_integrate();
	
	if(n_post_integrate_host != 0)
	{
		cuda->downloadAll(); cuda->oncpu = true;
		for (int i = 0; i < n_post_integrate_host; i++)
			fix[list_post_integrate[i]]->post_integrate();
		cuda->uploadAll(); cuda->oncpu = false;
	}
}

/* ----------------------------------------------------------------------
   pre_exchange call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyCuda::pre_exchange()
{
	for(int i = 0; i < n_pre_exchange_cuda; i++)
		fix[list_pre_exchange_cuda[i]]->pre_exchange();

	if(n_pre_exchange_host != 0)
	{
		cuda->downloadAll(); cuda->oncpu = true;
		for (int i = 0; i < n_pre_exchange_host; i++)
			fix[list_pre_exchange[i]]->pre_exchange();
		cuda->uploadAll(); cuda->oncpu = false;
	}
}

/* ----------------------------------------------------------------------
   pre_neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyCuda::pre_neighbor()
{
	for(int i = 0; i < n_pre_neighbor_cuda; i++)
		fix[list_pre_neighbor_cuda[i]]->pre_neighbor();
	
	if(n_pre_neighbor_host != 0)
	{
		cuda->downloadAll(); cuda->oncpu = true;
		for (int i = 0; i < n_pre_neighbor_host; i++)
			fix[list_pre_neighbor[i]]->pre_neighbor();
		cuda->uploadAll(); cuda->oncpu = false;
	}
}

/* ----------------------------------------------------------------------
   pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyCuda::setup_pre_force(int vflag)
{
	for(int i = 0; i < n_pre_force_cuda; i++)
		fix[list_pre_force_cuda[i]]->pre_force(vflag);

	if(n_pre_force_host != 0)
	{
		cuda->downloadAll(); cuda->oncpu = true;
		for (int i = 0; i < n_pre_force_host; i++)
			fix[list_pre_force[i]]->pre_force(vflag);
		cuda->uploadAll(); cuda->oncpu = false;
	}
}

void ModifyCuda::pre_force(int vflag)
{
	for(int i = 0; i < n_pre_force_cuda; i++)
		fix[list_pre_force_cuda[i]]->pre_force(vflag);

	if(n_pre_force_host != 0)
	{
		cuda->downloadAll(); cuda->oncpu = true;
		for (int i = 0; i < n_pre_force_host; i++)
			fix[list_pre_force[i]]->pre_force(vflag);
		cuda->uploadAll(); cuda->oncpu = false;
	}
}

/* ----------------------------------------------------------------------
   post_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyCuda::post_force(int vflag)
{
	for(int i = 0; i < n_post_force_cuda; i++)
			fix[list_post_force_cuda[i]]->post_force(vflag);

	if(n_post_force_host != 0)
	{
		cuda->downloadAll(); cuda->oncpu = true;
		for (int i = 0; i < n_post_force_host; i++)
			fix[list_post_force[i]]->post_force(vflag);
		cuda->uploadAll(); cuda->oncpu = false;
	}
}

/* ----------------------------------------------------------------------
   2nd half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyCuda::final_integrate()
{
	for (int i = 0; i < n_final_integrate_cuda; i++)
		fix[list_final_integrate_cuda[i]]->final_integrate();
	
	if(n_final_integrate_host != 0)
	{
		cuda->downloadAll(); cuda->oncpu = true;
		for (int i = 0; i < n_final_integrate_host; i++)
			fix[list_final_integrate[i]]->final_integrate();
		cuda->uploadAll(); cuda->oncpu = false;
	}
}

/* ----------------------------------------------------------------------
   end-of-timestep call, only for relevant fixes
   only call fix->end_of_step() on timesteps that are multiples of nevery
------------------------------------------------------------------------- */

void ModifyCuda::end_of_step()
{
	for (int i = 0; i < n_end_of_step_cuda; i++)
		if (update->ntimestep % end_of_step_every_cuda[i] == 0)
			fix[list_end_of_step_cuda[i]]->end_of_step();
	
	if(n_end_of_step_host != 0)
	{
		int do_thisstep=0;
		for (int i = 0; i < n_end_of_step_host; i++)
			if (update->ntimestep % end_of_step_every[i] == 0) do_thisstep=1;
		if(do_thisstep)
		{
		  cuda->downloadAll(); cuda->oncpu = true;
		  for (int i = 0; i < n_end_of_step_host; i++)
			 if (update->ntimestep % end_of_step_every[i] == 0)
				fix[list_end_of_step[i]]->end_of_step();
		  cuda->uploadAll(); cuda->oncpu = false;
		}
	}
}

/* ----------------------------------------------------------------------
   thermo energy call, only for relevant fixes
   called by Thermo class
   compute_scalar() is fix call to return energy
------------------------------------------------------------------------- */

double ModifyCuda::thermo_energy()
{
	double energy = 0.0;
	
	for (int i = 0; i < n_thermo_energy_cuda; i++)
		energy += fix[list_thermo_energy_cuda[i]]->compute_scalar();
	
	if(n_thermo_energy_host != 0)
	{
		cuda->downloadAll(); cuda->oncpu = true;
		for (int i = 0; i < n_thermo_energy_host; i++)
			energy += fix[list_thermo_energy[i]]->compute_scalar();
		cuda->uploadAll(); cuda->oncpu = false;
	}
	
	return energy;
}



void ModifyCuda::list_init_end_of_step_cuda(int mask, int &n, int *&list)
{
  delete [] list;
  delete [] end_of_step_every_cuda;

  n = 0;
  for (int i = 0; i < nfix; i++) if (fmask[i] & mask) n++;
  list = new int[n];
  end_of_step_every_cuda = new int[n];

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask) {
      list[n] = i;
      end_of_step_every_cuda[n++] = fix[i]->nevery;
    }
}
