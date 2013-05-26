/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_gran_hooke_cuda.h"
#include "pair_gran_hooke_cuda_cu.h"
#include "cuda_data.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "modify.h"
#include "fix_pour.h"
#include "cuda_neigh_list.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "cuda.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGranHookeCuda::PairGranHookeCuda(LAMMPS *lmp) : PairGranHooke(lmp)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

        allocated2 = false;
        cuda->shared_data.pair.cudable_force = 1;
        cuda->setSystemParams();
}

/* ----------------------------------------------------------------------
   remember pointer to arrays in cuda shared data
------------------------------------------------------------------------- */

void PairGranHookeCuda::allocate()
{
        if(! allocated) PairGranHooke::allocate();
        if(! allocated2)
        {
                allocated2 = true;
                 int n = atom->ntypes;
                cuda->shared_data.pair.cutsq     = cutsq;
                memory->create(cuda->shared_data.pair.coeff1,n+1,n+1,
                               "pair:cuda_coeff1");
                memory->create(cuda->shared_data.pair.coeff2,
                               n+1,n+1,"pair:cuda_coeff2");
                cuda->shared_data.pair.coeff1[0][0]=kn;
                cuda->shared_data.pair.coeff1[0][1]=kt;
                cuda->shared_data.pair.coeff1[1][0]=gamman;
                cuda->shared_data.pair.coeff1[1][1]=gammat;
                cuda->shared_data.pair.coeff2[0][0]=xmu;
                cuda->shared_data.pair.coeff2[0][1]=dampflag;
        }
}

/* ---------------------------------------------------------------------- */

void PairGranHookeCuda::compute(int eflag, int vflag)
{
             cuda->shared_data.pair.use_block_per_atom = 0;
        //cuda->cu_debugdata->memset_device(0);
        if (eflag || vflag) ev_setup(eflag,vflag);
        if(eflag) cuda->cu_eng_vdwl->upload();
        if(vflag) cuda->cu_virial->upload();

        Cuda_PairGranHookeCuda(& cuda->shared_data, & cuda_neigh_list->sneighlist, eflag, vflag, eflag_atom, vflag_atom);

    if(not cuda->shared_data.pair.collect_forces_later)
    {
          if(eflag) cuda->cu_eng_vdwl->download();
          if(vflag) cuda->cu_virial->download();
    }
        //cuda->cu_debugdata->download();
        //printf("%lf %lf %lf %lf %lf %lf\n",1.0e-6*cuda->debugdata[0],1.0e-6*cuda->debugdata[1],1.0e-6*cuda->debugdata[2],1.0e-6*cuda->debugdata[3],1.0e-6*cuda->debugdata[4],1.0e-6*cuda->debugdata[5]);

}

/* ---------------------------------------------------------------------- */

void PairGranHookeCuda::settings(int narg, char **arg)
{
        PairGranHooke::settings(narg, arg);
 }

/* ---------------------------------------------------------------------- */

void PairGranHookeCuda::coeff(int narg, char **arg)
{
        PairGranHooke::coeff(narg, arg);
        allocate();
}

void PairGranHookeCuda::init_style()
{
        int i;
        MYDBG(printf("# CUDA PairGranHookeCuda::init_style start\n"); )
  // request regular or rRESPA neighbor lists

  int irequest;

  if (update->whichflag == 0 && strstr(update->integrate_style,"respa")) {

  }
  else
  {
          irequest = neighbor->request(this);
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->gran = 1;
    neighbor->requests[irequest]->cudable = 1;
    //neighbor->style=0; //0=NSQ neighboring
  }

  if (!atom->radius_flag || !atom->omega_flag || !atom->torque_flag)
    error->all(FLERR,"Pair granular requires atom attributes radius, omega, torque");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair granular requires ghost atoms store velocity");

  // need a half neigh list and optionally a granular history neigh list

  dt = update->dt;



  // check for Fix freeze and set freeze_group_bit

  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"freeze") == 0) break;
  if (i < modify->nfix) freeze_group_bit = modify->fix[i]->groupbit;
  else freeze_group_bit = 0;

  cuda->shared_data.pair.freeze_group_bit=freeze_group_bit;
  // check for Fix pour and set pour_type and pour_maxdiam

  int pour_type = 0;
  double pour_maxrad = 0.0;
  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"pour") == 0) break;
  if (i < modify->nfix) {
    pour_type = ((FixPour *) modify->fix[i])->ntype;
    pour_maxrad = ((FixPour *) modify->fix[i])->radius_max;
  }

  // set maxrad_dynamic and maxrad_frozen for each type
  // include future Fix pour particles as dynamic

  for (i = 1; i <= atom->ntypes; i++)
    onerad_dynamic[i] = onerad_frozen[i] = 0.0;
  if (pour_type) onerad_dynamic[pour_type] = pour_maxrad;

  double *radius = atom->radius;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++){
    if (mask[i] & freeze_group_bit)
      onerad_frozen[type[i]] = MAX(onerad_frozen[type[i]],radius[i]);
    else
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]],radius[i]);
  }

  MPI_Allreduce(&onerad_dynamic[1],&maxrad_dynamic[1],atom->ntypes,
                MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&onerad_frozen[1],&maxrad_frozen[1],atom->ntypes,
                MPI_DOUBLE,MPI_MAX,world);

  MYDBG(printf("# CUDA PairGranHookeCuda::init_style end\n"); )
}

void PairGranHookeCuda::init_list(int id, NeighList *ptr)
{
        MYDBG(printf("# CUDA PairGranHookeCuda::init_list\n");)
        PairGranHooke::init_list(id, ptr);
        #ifndef CUDA_USE_BINNING
        // right now we can only handle verlet (id 0), not respa
        if(id == 0) cuda_neigh_list = cuda->registerNeighborList(ptr);
        // see Neighbor::init() for details on lammps lists' logic
        #endif
        MYDBG(printf("# CUDA PairGranHookeCuda::init_list end\n");)
}

void PairGranHookeCuda::ev_setup(int eflag, int vflag)
{
        int maxeatomold=maxeatom;
        PairGranHooke::ev_setup(eflag,vflag);

  if (eflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_eatom; cuda->cu_eatom = new cCudaData<double, ENERGY_FLOAT, x > ((double*)eatom, & cuda->shared_data.atom.eatom , atom->nmax  );}

  if (eflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_vatom; cuda->cu_vatom = new cCudaData<double, ENERGY_FLOAT, yx > ((double*)vatom, & cuda->shared_data.atom.eatom , atom->nmax, 6  );}

}
