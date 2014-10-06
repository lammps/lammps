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
#include "pair_lj_class2_coul_long_cuda.h"
#include "pair_lj_class2_coul_long_cuda_cu.h"
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
#include "user_cuda.h"

using namespace LAMMPS_NS;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429
/* ---------------------------------------------------------------------- */

PairLJClass2CoulLongCuda::PairLJClass2CoulLongCuda(LAMMPS *lmp) : PairLJClass2CoulLong(lmp)
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

void PairLJClass2CoulLongCuda::allocate()
{
        if(! allocated) PairLJClass2CoulLong::allocate();
        if(! allocated2)
        {
                allocated2 = true;
                cuda->shared_data.pair.cut     = cut_lj;
                cuda->shared_data.pair.coeff1  = lj1;
                cuda->shared_data.pair.coeff2  = lj2;
                cuda->shared_data.pair.coeff3  = lj3;
                cuda->shared_data.pair.coeff4  = lj4;
                cuda->shared_data.pair.offset  = offset;
                cuda->shared_data.pair.special_lj  = force->special_lj;
                cuda->shared_data.pair.special_coul  = force->special_coul;
        }
}

/* ---------------------------------------------------------------------- */

void PairLJClass2CoulLongCuda::compute(int eflag, int vflag)
{
        if (eflag || vflag) ev_setup(eflag,vflag);
        if(eflag) cuda->cu_eng_vdwl->upload();
        if(eflag) cuda->cu_eng_coul->upload();
        if(vflag) cuda->cu_virial->upload();

        Cuda_PairLJClass2CoulLongCuda(& cuda->shared_data, & cuda_neigh_list->sneighlist, eflag, vflag, eflag_atom, vflag_atom);

    if(not cuda->shared_data.pair.collect_forces_later)
    {
          if(eflag) cuda->cu_eng_vdwl->download();
          if(eflag) cuda->cu_eng_coul->download();
          if(vflag) cuda->cu_virial->download();
    }
}

/* ---------------------------------------------------------------------- */

void PairLJClass2CoulLongCuda::settings(int narg, char **arg)
{
        PairLJClass2CoulLong::settings(narg, arg);
        cuda->shared_data.pair.cut_global = (F_CFLOAT) cut_lj_global;
}

/* ---------------------------------------------------------------------- */

void PairLJClass2CoulLongCuda::coeff(int narg, char **arg)
{
        PairLJClass2CoulLong::coeff(narg, arg);
        allocate();
}

void PairLJClass2CoulLongCuda::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style lj/cut/coul/long requires atom attribute q");
  // request regular or rRESPA neighbor lists

  int irequest;


          irequest = neighbor->request(this);
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->cudable = 1;

  cut_coulsq = cut_coul * cut_coul;
  cuda->shared_data.pair.cut_coul_global=cut_coul;
  cuda->shared_data.pair.cut_coulsq_global=cut_coulsq;
  // set rRESPA cutoffs

  if (force->newton) error->warning(FLERR,"Pair style uses does not use \"newton\" setting. You might test if \"newton off\" makes the simulation run faster.");
  if (force->kspace == NULL)
    error->all(FLERR,"Pair style is incompatible with KSpace style");
  g_ewald = force->kspace->g_ewald;
  cuda->shared_data.pair.g_ewald=g_ewald;
  cuda->shared_data.pppm.qqrd2e=force->qqrd2e;


  if(ncoultablebits) error->warning(FLERR,"# CUDA: You asked for the usage of Coulomb Tables. This is not supported in CUDA Pair forces. Setting is ignored.\n");
}

void PairLJClass2CoulLongCuda::init_list(int id, NeighList *ptr)
{
        MYDBG(printf("# CUDA PairLJClass2CoulLongCuda::init_list\n");)
        PairLJClass2CoulLong::init_list(id, ptr);
        if(id == 0) cuda_neigh_list = cuda->registerNeighborList(ptr);
        // see Neighbor::init() for details on lammps lists' logic
        MYDBG(printf("# CUDA PairLJClass2CoulLongCuda::init_list end\n");)
}

void PairLJClass2CoulLongCuda::ev_setup(int eflag, int vflag)
{
        int maxeatomold=maxeatom;
        PairLJClass2CoulLong::ev_setup(eflag,vflag);

  if (eflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_eatom; cuda->cu_eatom = new cCudaData<double, ENERGY_CFLOAT, x > ((double*)eatom, & cuda->shared_data.atom.eatom , atom->nmax  );}

  if (vflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_vatom; cuda->cu_vatom = new cCudaData<double, ENERGY_CFLOAT, yx > ((double*)vatom, & cuda->shared_data.atom.vatom , atom->nmax, 6  );}

}
