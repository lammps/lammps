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
#include "pair_morse_cuda.h"
#include "pair_morse_cuda_cu.h"
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

PairMorseCuda::PairMorseCuda(LAMMPS *lmp) : PairMorse(lmp)
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

void PairMorseCuda::allocate()
{
        if(! allocated) PairMorse::allocate();
        if(! allocated2)
        {
                allocated2 = true;
                cuda->shared_data.pair.cut     = cut;
                cuda->shared_data.pair.coeff1  = r0;
                cuda->shared_data.pair.coeff2  = alpha;
                cuda->shared_data.pair.coeff3  = morse1;
                cuda->shared_data.pair.coeff4  = d0;
                cuda->shared_data.pair.offset  = offset;
                cuda->shared_data.pair.special_lj  = force->special_lj;
        }
}

/* ---------------------------------------------------------------------- */

void PairMorseCuda::compute(int eflag, int vflag)
{
        if (eflag || vflag) ev_setup(eflag,vflag);
        if(eflag) cuda->cu_eng_vdwl->upload();
        if(vflag) cuda->cu_virial->upload();

        Cuda_PairMorseCuda(& cuda->shared_data, & cuda_neigh_list->sneighlist, eflag, vflag, eflag_atom, vflag_atom);

    if(not cuda->shared_data.pair.collect_forces_later)
    {
          if(eflag) cuda->cu_eng_vdwl->download();
          if(vflag) cuda->cu_virial->download();
    }

}

/* ---------------------------------------------------------------------- */

void PairMorseCuda::settings(int narg, char **arg)
{
        PairMorse::settings(narg, arg);
        cuda->shared_data.pair.cut_global = (F_FLOAT) cut_global;
}

/* ---------------------------------------------------------------------- */

void PairMorseCuda::coeff(int narg, char **arg)
{
        PairMorse::coeff(narg, arg);
        allocate();
}

void PairMorseCuda::init_style()
{
        MYDBG(printf("# CUDA PairMorseCuda::init_style start\n"); )
  // request regular or rRESPA neighbor lists

  int irequest;

  if (update->whichflag == 0 && strstr(update->integrate_style,"respa")) {

  }
  else
  {
          irequest = neighbor->request(this);
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->cudable = 1;
    //neighbor->style=0; //0=NSQ neighboring
  }


  MYDBG(printf("# CUDA PairMorseCuda::init_style end\n"); )
}

void PairMorseCuda::init_list(int id, NeighList *ptr)
{
        MYDBG(printf("# CUDA PairMorseCuda::init_list\n");)
        PairMorse::init_list(id, ptr);
        #ifndef CUDA_USE_BINNING
        // right now we can only handle verlet (id 0), not respa
        if(id == 0) cuda_neigh_list = cuda->registerNeighborList(ptr);
        // see Neighbor::init() for details on lammps lists' logic
        #endif
        MYDBG(printf("# CUDA PairMorseCuda::init_list end\n");)
}

void PairMorseCuda::ev_setup(int eflag, int vflag)
{
        int maxeatomold=maxeatom;
        PairMorse::ev_setup(eflag,vflag);

  if (eflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_eatom; cuda->cu_eatom = new cCudaData<double, ENERGY_FLOAT, x > ((double*)eatom, & cuda->shared_data.atom.eatom , atom->nmax  );}

  if (vflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_vatom; cuda->cu_vatom = new cCudaData<double, ENERGY_FLOAT, yx > ((double*)vatom, & cuda->shared_data.atom.vatom , atom->nmax, 6  );}

}
