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
#include "pair_lj_class2_cuda.h"
#include "pair_lj_class2_cuda_cu.h"
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
#include "cuda.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJClass2Cuda::PairLJClass2Cuda(LAMMPS *lmp) : PairLJClass2(lmp)
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

void PairLJClass2Cuda::allocate()
{
        if(! allocated) PairLJClass2::allocate();
        if(! allocated2)
        {
                allocated2 = true;
                cuda->shared_data.pair.cut     = cut;
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

void PairLJClass2Cuda::compute(int eflag, int vflag)
{
        if (eflag || vflag) ev_setup(eflag,vflag);
        if(eflag) cuda->cu_eng_vdwl->upload();
        if(vflag) cuda->cu_virial->upload();

        Cuda_PairLJClass2Cuda(& cuda->shared_data, & cuda_neigh_list->sneighlist, eflag, vflag, eflag_atom, vflag_atom);

    if(not cuda->shared_data.pair.collect_forces_later)
    {
          if(eflag) cuda->cu_eng_vdwl->download();
          if(vflag) cuda->cu_virial->download();
    }

}

/* ---------------------------------------------------------------------- */

void PairLJClass2Cuda::settings(int narg, char **arg)
{
        PairLJClass2::settings(narg, arg);
        cuda->shared_data.pair.cut_global = (F_FLOAT) cut_global;
}

/* ---------------------------------------------------------------------- */

void PairLJClass2Cuda::coeff(int narg, char **arg)
{
        PairLJClass2::coeff(narg, arg);
        allocate();
}

void PairLJClass2Cuda::init_style()
{
        MYDBG(printf("# CUDA PairLJClass2Cuda::init_style start\n"); )
  // request regular or rRESPA neighbor lists

  int irequest;

          irequest = neighbor->request(this);
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->cudable = 1;
    //neighbor->style=0; //0=NSQ neighboring
  MYDBG(printf("# CUDA PairLJClass2Cuda::init_style end\n"); )
}

void PairLJClass2Cuda::init_list(int id, NeighList *ptr)
{
        MYDBG(printf("# CUDA PairLJClass2Cuda::init_list\n");)
        PairLJClass2::init_list(id, ptr);
        if(id == 0) cuda_neigh_list = cuda->registerNeighborList(ptr);
        // see Neighbor::init() for details on lammps lists' logic
        MYDBG(printf("# CUDA PairLJClass2Cuda::init_list end\n");)
}

void PairLJClass2Cuda::ev_setup(int eflag, int vflag)
{
        int maxeatomold=maxeatom;
        PairLJClass2::ev_setup(eflag,vflag);

  if (eflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_eatom; cuda->cu_eatom = new cCudaData<double, ENERGY_FLOAT, x > ((double*)eatom, & cuda->shared_data.atom.eatom , atom->nmax  );}

  if (vflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_vatom; cuda->cu_vatom = new cCudaData<double, ENERGY_FLOAT, yx > ((double*)vatom, & cuda->shared_data.atom.vatom , atom->nmax, 6  );}

}
