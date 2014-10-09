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
#include "pair_buck_coul_cut_cuda.h"
#include "pair_buck_coul_cut_cuda_cu.h"
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

/* ---------------------------------------------------------------------- */

PairBuckCoulCutCuda::PairBuckCoulCutCuda(LAMMPS *lmp) : PairBuckCoulCut(lmp)
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

void PairBuckCoulCutCuda::allocate()
{
        if(! allocated) PairBuckCoulCut::allocate();
        if(! allocated2)
        {
                allocated2 = true;
                cuda->shared_data.pair.cut_coul     = cut_coul;
                cuda->shared_data.pair.cut     = cut_lj;
                cuda->shared_data.pair.coeff1  = rhoinv;
                cuda->shared_data.pair.coeff2  = buck1;
                cuda->shared_data.pair.coeff3  = buck2;
                cuda->shared_data.pair.coeff4  = a;
                cuda->shared_data.pair.coeff5  = c;
                cuda->shared_data.pair.offset  = offset;
                cuda->shared_data.pair.special_lj  = force->special_lj;
                cuda->shared_data.pair.special_coul  = force->special_coul;
        }
}

/* ---------------------------------------------------------------------- */

void PairBuckCoulCutCuda::compute(int eflag, int vflag)
{
        MYDBG( printf("PairBuckCoulCutCuda compute start\n"); fflush(stdout);)
        if (eflag || vflag) ev_setup(eflag,vflag);
        if(eflag) cuda->cu_eng_vdwl->upload();
        if(eflag) cuda->cu_eng_coul->upload();
        if(vflag) cuda->cu_virial->upload();

        Cuda_PairBuckCoulCutCuda(& cuda->shared_data, & cuda_neigh_list->sneighlist, eflag, vflag, eflag_atom, vflag_atom);

    if(not cuda->shared_data.pair.collect_forces_later)
    {
          if(eflag) cuda->cu_eng_vdwl->download();
          if(eflag) cuda->cu_eng_coul->download();
          if(vflag) cuda->cu_virial->download();
    }
        MYDBG( printf("PairBuckCoulCutCuda compute end\n"); fflush(stdout);)
}

/* ---------------------------------------------------------------------- */

void PairBuckCoulCutCuda::settings(int narg, char **arg)
{
        PairBuckCoulCut::settings(narg, arg);
        cuda->shared_data.pair.cut_coul_global = (F_CFLOAT) cut_coul_global;
        cuda->shared_data.pair.cut_global = (F_CFLOAT) cut_lj_global;
}

/* ---------------------------------------------------------------------- */

void PairBuckCoulCutCuda::coeff(int narg, char **arg)
{
        PairBuckCoulCut::coeff(narg, arg);
        allocate();
}

void PairBuckCoulCutCuda::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style buck/coul/long requires atom attribute q");
  // request regular or rRESPA neighbor lists

  int irequest;

  if (strstr(update->integrate_style,"respa")) error->all(FLERR,"Integrate Style Respa is not supported by pair style buck/coul/long/cuda");

          irequest = neighbor->request(this);
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->cudable = 1;


  cuda->shared_data.pppm.qqrd2e=force->qqrd2e;

   cuda->shared_data.pair.cut_coulsq_global=cut_coul_global * cut_coul_global;

  if(ncoultablebits) error->warning(FLERR,"# CUDA: You asked for the usage of Coulomb Tables. This is not supported in CUDA Pair forces. Setting is ignored.\n");
}

void PairBuckCoulCutCuda::init_list(int id, NeighList *ptr)
{
        MYDBG(printf("# CUDA PairBuckCoulCutCuda::init_list\n");)
        PairBuckCoulCut::init_list(id, ptr);
        #ifndef CUDA_USE_BINNING
        // right now we can only handle verlet (id 0), not respa
        if(id == 0) cuda_neigh_list = cuda->registerNeighborList(ptr);
        // see Neighbor::init() for details on lammps lists' logic
        #endif
        MYDBG(printf("# CUDA PairBuckCoulCutCuda::init_list end\n");)
}

void PairBuckCoulCutCuda::ev_setup(int eflag, int vflag)
{
        int maxeatomold=maxeatom;
        PairBuckCoulCut::ev_setup(eflag,vflag);

  if (eflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_eatom; cuda->cu_eatom = new cCudaData<double, ENERGY_CFLOAT, x > ((double*)eatom, & cuda->shared_data.atom.eatom , atom->nmax  );}

  if (vflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_vatom; cuda->cu_vatom = new cCudaData<double, ENERGY_CFLOAT, yx > ((double*)vatom, & cuda->shared_data.atom.vatom , atom->nmax, 6  );}

}
