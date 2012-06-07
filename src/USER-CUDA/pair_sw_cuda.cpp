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
#include "pair_sw_cuda.h"
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

PairSWCuda::PairSWCuda(LAMMPS *lmp) : PairSW(lmp)
{
  cuda = lmp->cuda;
  if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

        allocated2 = false;
        params_f = NULL;
        cuda->setSystemParams();
  cuda->shared_data.pair.cudable_force = 1;
  cuda->shared_data.pair.override_block_per_atom = 0;
  cuda->shared_data.pair.neighall = true;
  init = false;
}

/* ----------------------------------------------------------------------
   remember pointer to arrays in cuda shared data
------------------------------------------------------------------------- */

void PairSWCuda::allocate()
{
        if(! allocated) PairSW::allocate();
        if(! allocated2)
        {
                allocated2 = true;
                cuda->shared_data.pair.cutsq   = cutsq;
                cuda->shared_data.pair.special_lj  = force->special_lj;
                cuda->shared_data.pair.special_coul  = force->special_coul;
        }
}

/* ---------------------------------------------------------------------- */

void PairSWCuda::compute(int eflag, int vflag)
{
  if(!init) {Cuda_PairSWCuda_Init(&cuda->shared_data,params_f,map, &elem2param[0][0][0],nelements); init=true;}
        if (eflag || vflag) ev_setup(eflag,vflag);
        if(eflag) cuda->cu_eng_vdwl->upload();
        if(vflag) cuda->cu_virial->upload();

        Cuda_PairSWCuda(& cuda->shared_data, & cuda_neigh_list->sneighlist, eflag, vflag, eflag_atom, vflag_atom);//,&elem2param[0][0][0],map
  if(not cuda->shared_data.pair.collect_forces_later)
  {
          if(eflag) cuda->cu_eng_vdwl->download();
          if(vflag) cuda->cu_virial->download();
  }
}

/* ---------------------------------------------------------------------- */

void PairSWCuda::settings(int narg, char **arg)
{
        PairSW::settings(narg, arg);
}

/* ---------------------------------------------------------------------- */

void PairSWCuda::coeff(int narg, char **arg)
{
        PairSW::coeff(narg, arg);
        allocate();
  params_f = (ParamSW_Float *) memory->srealloc(params_f,maxparam*sizeof(ParamSW_Float),
        "pair:params_f");
  for(int i=0;i<maxparam;i++)
  {
    printf("%e %e\n",params[i].cut,params[i].cutsq);
    params_f[i].cut = params[i].cut;
    params_f[i].cutsq = params[i].cutsq;
    params_f[i].c1 = params[i].c1;
    params_f[i].c2 = params[i].c2;
    params_f[i].c3 = params[i].c3;
    params_f[i].c4 = params[i].c4;
    params_f[i].c5 = params[i].c5;
    params_f[i].c6 = params[i].c6;
    params_f[i].ielement = params[i].ielement;
    params_f[i].jelement = params[i].jelement;
    params_f[i].kelement = params[i].kelement;
    params_f[i].epsilon = params[i].epsilon;
    params_f[i].sigma = params[i].sigma;
    params_f[i].littlea = params[i].littlea;
    params_f[i].lambda = params[i].lambda;
    params_f[i].costheta = params[i].costheta;
    params_f[i].tol = params[i].tol;
    params_f[i].sigma_gamma = params[i].sigma_gamma;
    params_f[i].lambda_epsilon = params[i].lambda_epsilon;
    params_f[i].lambda_epsilon2 = params[i].lambda_epsilon2;
    params_f[i].gamma = params[i].gamma;

    params_f[i].biga = params[i].biga;
    params_f[i].bigb = params[i].bigb;
    params_f[i].gamma = params[i].gamma;
    params_f[i].powerp = params[i].powerp;
    params_f[i].powerq = params[i].powerq;
  }
  cuda->shared_data.pair.cut_global = cutmax;
}

void PairSWCuda::init_style()
{
        MYDBG(printf("# CUDA PairSWCuda::init_style start\n"); )

  int irequest;

        irequest = neighbor->request(this);
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->cudable = 1;
  neighbor->requests[irequest]->ghost = 1;


  MYDBG(printf("# CUDA PairSWCuda::init_style end\n"); )
}

void PairSWCuda::init_list(int id, NeighList *ptr)
{
        MYDBG(printf("# CUDA PairSWCuda::init_list\n");)
        PairSW::init_list(id, ptr);
        // right now we can only handle verlet (id 0), not respa
        if(id == 0) cuda_neigh_list = cuda->registerNeighborList(ptr);
        // see Neighbor::init() for details on lammps lists' logic
        MYDBG(printf("# CUDA PairSWCuda::init_list end\n");)
  cu_params_f = (ParamSW_Float*) CudaWrapper_AllocCudaData(sizeof(ParamSW_Float)*maxparam);
  CudaWrapper_UploadCudaData((void*) params_f,(void*) cu_params_f,sizeof(ParamSW_Float)*maxparam);
  cu_elem2param = new cCudaData<int, int, xyz > ((int*) elem2param, nelements,nelements,nelements);
  cu_elem2param->upload();
  cu_map = new cCudaData<int, int, x > ( map,atom->ntypes+1 );
  cu_map->upload();
}

void PairSWCuda::ev_setup(int eflag, int vflag)
{
        int maxeatomold=maxeatom;
        PairSW::ev_setup(eflag,vflag);

  if (eflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_eatom; cuda->cu_eatom = new cCudaData<double, ENERGY_FLOAT, x > ((double*)eatom, & cuda->shared_data.atom.eatom , atom->nmax  );}

  if (vflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_vatom; cuda->cu_vatom = new cCudaData<double, ENERGY_FLOAT, yx > ((double*)vatom, & cuda->shared_data.atom.vatom , atom->nmax, 6  );}
}
