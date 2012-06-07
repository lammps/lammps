/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   Contributing author: Paul Crozier (SNL)
   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_lj_charmm_coul_charmm_cuda.h"
#include "pair_lj_charmm_coul_charmm_cuda_cu.h"
#include "cuda_data.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "cuda_neigh_list.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "cuda.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCharmmCoulCharmmCuda::PairLJCharmmCoulCharmmCuda(LAMMPS *lmp) : PairLJCharmmCoulCharmm(lmp)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

        allocated2 = false;
        cuda->shared_data.pair.cudable_force = 1;
        cuda->shared_data.pair.use_block_per_atom = 0;
        cuda->setSystemParams();
}

/* ----------------------------------------------------------------------
   remember pointer to arrays in cuda shared data
------------------------------------------------------------------------- */

void PairLJCharmmCoulCharmmCuda::allocate()
{
        if(! allocated) PairLJCharmmCoulCharmm::allocate();
        if(! allocated2)
        {
                cuda->accelerator(0,NULL);
                allocated2 = true;
                cuda->shared_data.pair.coeff1  = lj1;
                cuda->shared_data.pair.coeff2  = lj2;
                cuda->shared_data.pair.coeff3  = lj3;
                cuda->shared_data.pair.coeff4  = lj4;
                cuda->shared_data.pair.special_lj  = force->special_lj;
                cuda->shared_data.pair.special_coul  = force->special_coul;
            cu_lj1_gm = new cCudaData<double, F_FLOAT, x> ((double*)lj1, &cuda->shared_data.pair.coeff1_gm, (atom->ntypes+1)*(atom->ntypes+1));
            cu_lj2_gm = new cCudaData<double, F_FLOAT, x> ((double*)lj2, &cuda->shared_data.pair.coeff2_gm, (atom->ntypes+1)*(atom->ntypes+1));
            cu_lj3_gm = new cCudaData<double, F_FLOAT, x> ((double*)lj3, &cuda->shared_data.pair.coeff3_gm, (atom->ntypes+1)*(atom->ntypes+1));
            cu_lj4_gm = new cCudaData<double, F_FLOAT, x> ((double*)lj4, &cuda->shared_data.pair.coeff4_gm, (atom->ntypes+1)*(atom->ntypes+1));
        }
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulCharmmCuda::compute(int eflag, int vflag)
{
          if (eflag || vflag) ev_setup(eflag,vflag);
        if(not cuda->shared_data.pair.collect_forces_later)
        {
          if(eflag) cuda->cu_eng_vdwl->upload();
          if(eflag) cuda->cu_eng_coul->upload();
          if(vflag) cuda->cu_virial->upload();
        }

        Cuda_PairLJCharmmCoulCharmmCuda(& cuda->shared_data, & cuda_neigh_list->sneighlist, eflag, vflag, eflag_atom, vflag_atom,denom_lj,cut_coul_innersq,denom_coul);

        if(not cuda->shared_data.pair.collect_forces_later)
        {
          if(eflag) cuda->cu_eng_vdwl->download();
          if(eflag) cuda->cu_eng_coul->download();
          if(vflag) cuda->cu_virial->download();
        }
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulCharmmCuda::settings(int narg, char **arg)
{
        PairLJCharmmCoulCharmm::settings(narg, arg);
        cuda->shared_data.pair.cut_global = (X_FLOAT) cut_lj;
        cuda->shared_data.pair.cut_coulsq_global = (X_FLOAT) cut_coulsq;
        cuda->shared_data.pair.cut_inner_global = (F_FLOAT) cut_lj_inner;
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulCharmmCuda::coeff(int narg, char **arg)
{
        PairLJCharmmCoulCharmm::coeff(narg, arg);
        allocate();
}

void PairLJCharmmCoulCharmmCuda::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style lj/charmm/coul/long requires atom attribute q");
  // request regular or rRESPA neighbor lists

        if(atom->molecular)
        {
          cuda->shared_data.pair.collect_forces_later = 1;
        }

  int irequest;

           irequest = neighbor->request(this);
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->cudable = 1;

   if (cut_lj_inner >= cut_lj || cut_coul_inner >= cut_coul)
    error->all(FLERR,"Pair inner cutoff >= Pair outer cutoff");

  cut_lj_innersq = cut_lj_inner * cut_lj_inner;
  cut_ljsq = cut_lj * cut_lj;
  cut_coul_innersq = cut_coul_inner * cut_coul_inner;
  cut_coulsq = cut_coul * cut_coul;
  cut_bothsq = MAX(cut_ljsq,cut_coulsq);

  denom_lj = (cut_ljsq-cut_lj_innersq) * (cut_ljsq-cut_lj_innersq) *
    (cut_ljsq-cut_lj_innersq);
  denom_coul = (cut_coulsq-cut_coul_innersq) * (cut_coulsq-cut_coul_innersq) *
    (cut_coulsq-cut_coul_innersq);

  cut_coulsq = cut_coul * cut_coul;

  cuda->shared_data.pair.cut_coulsq_global=cut_coulsq;

  cuda->shared_data.pppm.qqrd2e=force->qqrd2e;
}

void PairLJCharmmCoulCharmmCuda::init_list(int id, NeighList *ptr)
{
        MYDBG(printf("# CUDA PairLJCharmmCoulCharmmCuda::init_list\n");)
        PairLJCharmmCoulCharmm::init_list(id, ptr);
        #ifndef CUDA_USE_BINNING
        // right now we can only handle verlet (id 0), not respa
        if(id == 0) cuda_neigh_list = cuda->registerNeighborList(ptr);
        // see Neighbor::init() for details on lammps lists' logic
        #endif
        MYDBG(printf("# CUDA PairLJCharmmCoulCharmmCuda::init_list end\n");)
}

void PairLJCharmmCoulCharmmCuda::ev_setup(int eflag, int vflag)
{
        int maxeatomold=maxeatom;
        PairLJCharmmCoulCharmm::ev_setup(eflag,vflag);

  if (eflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_eatom; cuda->cu_eatom = new cCudaData<double, ENERGY_FLOAT, x > ((double*)eatom, & cuda->shared_data.atom.eatom , atom->nmax  );}

  if (vflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_vatom; cuda->cu_vatom = new cCudaData<double, ENERGY_FLOAT, yx > ((double*)vatom, & cuda->shared_data.atom.vatom , atom->nmax, 6  );}

}
