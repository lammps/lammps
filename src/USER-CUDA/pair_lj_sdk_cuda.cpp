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
#include "pair_lj_sdk_cuda.h"
#include "pair_lj_sdk_cuda_cu.h"
#include "cuda_data.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
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

PairLJSDKCuda::PairLJSDKCuda(LAMMPS *lmp) : PairLJSDK(lmp)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

        allocated2 = false;
        lj_type_double = NULL;
        cuda->shared_data.pair.cudable_force = 1;
        cuda->setSystemParams();
}

/* ----------------------------------------------------------------------
   remember pointer to arrays in cuda shared data
------------------------------------------------------------------------- */

void PairLJSDKCuda::allocate()
{
        if(! allocated) PairLJSDK::allocate();
        int n = atom->ntypes;
        if(! allocated2)
        {
                allocated2 = true;


                  memory->create(lj_type_double,n+1,n+1,"pairlj:ljtypedouble");

                cuda->shared_data.pair.cut     = cut;
                cuda->shared_data.pair.coeff1  = lj1;
                cuda->shared_data.pair.coeff2  = lj2;
                cuda->shared_data.pair.coeff3  = lj3;
                cuda->shared_data.pair.coeff4  = lj4;
                cuda->shared_data.pair.coeff5  = lj_type_double;
            /*cu_lj1_gm = new cCudaData<double, F_FLOAT, x> ((double*)lj1, &cuda->shared_data.pair.coeff1_gm, (atom->ntypes+1)*(atom->ntypes+1));
            cu_lj2_gm = new cCudaData<double, F_FLOAT, x> ((double*)lj2, &cuda->shared_data.pair.coeff2_gm, (atom->ntypes+1)*(atom->ntypes+1));
            cu_lj3_gm = new cCudaData<double, F_FLOAT, x> ((double*)lj3, &cuda->shared_data.pair.coeff3_gm, (atom->ntypes+1)*(atom->ntypes+1));
            cu_lj4_gm = new cCudaData<double, F_FLOAT, x> ((double*)lj4, &cuda->shared_data.pair.coeff4_gm, (atom->ntypes+1)*(atom->ntypes+1));
            cu_lj_type_double_gm = new cCudaData<double, F_FLOAT, x> ((double*)lj_type_double, &cuda->shared_data.pair.coeff5_gm, (atom->ntypes+1)*(atom->ntypes+1));*/
                cuda->shared_data.pair.offset  = offset;
                cuda->shared_data.pair.special_lj  = force->special_lj;
        }
          for (int i = 1; i <= n; i++) {
      for (int j = i; j <= n; j++) {
        lj_type_double[i][j] = lj_type[i][j];
        lj_type_double[j][i] = lj_type[i][j];
      }
    }
}

/* ---------------------------------------------------------------------- */

void PairLJSDKCuda::compute(int eflag, int vflag)
{
        if (eflag || vflag) ev_setup(eflag,vflag);
        if(eflag) cuda->cu_eng_vdwl->upload();
        if(vflag) cuda->cu_virial->upload();

        Cuda_PairLJSDKCuda(& cuda->shared_data, & cuda_neigh_list->sneighlist, eflag, vflag, eflag_atom, vflag_atom);

    if(not cuda->shared_data.pair.collect_forces_later)
    {
          if(eflag) cuda->cu_eng_vdwl->download();
          if(vflag) cuda->cu_virial->download();
    }

}

/* ---------------------------------------------------------------------- */

void PairLJSDKCuda::settings(int narg, char **arg)
{
        PairLJSDK::settings(narg, arg);
        cuda->shared_data.pair.cut_global = (F_FLOAT) cut_global;
}

/* ---------------------------------------------------------------------- */

void PairLJSDKCuda::coeff(int narg, char **arg)
{
        PairLJSDK::coeff(narg, arg);
        allocate();
}

void PairLJSDKCuda::init_style()
{
  MYDBG(printf("# CUDA PairLJSDKCuda::init_style start\n"); )

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->cudable = 1;

  MYDBG(printf("# CUDA PairLJSDKCuda::init_style end\n"); )
}

void PairLJSDKCuda::init_list(int id, NeighList *ptr)
{
        MYDBG(printf("# CUDA PairLJSDKCuda::init_list\n");)
        PairLJSDK::init_list(id, ptr);
        #ifndef CUDA_USE_BINNING
        // right now we can only handle verlet (id 0), not respa
        if(id == 0) cuda_neigh_list = cuda->registerNeighborList(ptr);
        // see Neighbor::init() for details on lammps lists' logic
        #endif
        MYDBG(printf("# CUDA PairLJSDKCuda::init_list end\n");)
}

void PairLJSDKCuda::ev_setup(int eflag, int vflag)
{
        int maxeatomold=maxeatom;
        PairLJSDK::ev_setup(eflag,vflag);

  if (eflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_eatom; cuda->cu_eatom = new cCudaData<double, ENERGY_FLOAT, x > ((double*)eatom, & cuda->shared_data.atom.eatom , atom->nmax  );}

  if (vflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_vatom; cuda->cu_vatom = new cCudaData<double, ENERGY_FLOAT, yx > ((double*)vatom, & cuda->shared_data.atom.vatom , atom->nmax, 6  );}

}
