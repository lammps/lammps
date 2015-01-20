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
#include "pair_eam_cuda.h"
#include "pair_eam_cuda_cu.h"
#include "pair_virial_compute_cu.h"
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

PairEAMCuda::PairEAMCuda(LAMMPS* lmp) : PairEAM(lmp)
{
  cuda = lmp->cuda;

  if(cuda == NULL)
    error->all(FLERR, "You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

  allocated2 = false;
  cuda->shared_data.pair.cudable_force = 1;
  cuda->shared_data.pair.override_block_per_atom = 0;

  cuda->setSystemParams();
  cu_rho = NULL;
  cu_fp = NULL;
  cu_frho_spline = NULL;
  cu_z2r_spline = NULL;
  cu_rhor_spline = NULL;
}

/* ----------------------------------------------------------------------
   remember pointer to arrays in cuda shared data
------------------------------------------------------------------------- */

void PairEAMCuda::allocate()
{
  if(! allocated) PairEAM::allocate();

  cuda->shared_data.pair.cutsq     = cutsq;
  cuda->shared_data.pair.cut_global = (F_CFLOAT) cutforcesq;
}

/* ---------------------------------------------------------------------- */

void PairEAMCuda::compute(int eflag, int vflag)
{
  cuda->shared_data.pair.cut_global = (F_CFLOAT) cutforcesq;
  cuda->shared_data.pair.use_block_per_atom = 0;
  cuda->shared_data.pair.collect_forces_later = 0;

  if(atom->nmax > nmax || cuda->finished_setup == false) {
    memory->destroy(rho);
    memory->destroy(fp);
    nmax = atom->nmax;
    memory->create(rho, nmax, "pair:rho");
    memory->create(fp, nmax, "pair:fp");
    delete cu_rho;
    delete cu_fp;
    cu_rho = new cCudaData<double, F_CFLOAT, x> (rho, atom->nmax);
    cu_fp  = new cCudaData<double, F_CFLOAT, x> (fp, atom->nmax);
    Cuda_PairEAMCuda_Init(&cuda->shared_data, rdr, rdrho, nfrho, nrhor, nr, nrho, nz2r,
                          cu_frho_spline->dev_data(), cu_rhor_spline->dev_data(), cu_z2r_spline->dev_data(),
                          cu_rho->dev_data(), cu_fp->dev_data(), type2frho, type2z2r, type2rhor);
  }



  if(eflag || vflag) ev_setup(eflag, vflag);

  if(eflag) cuda->cu_eng_vdwl->upload();

  if(vflag) cuda->cu_virial->upload();

  Cuda_PairEAM1Cuda(& cuda->shared_data, & cuda_neigh_list->sneighlist, eflag, vflag, eflag_atom, vflag_atom);

  comm->forward_comm_pair(this);

  Cuda_PairEAM2Cuda(& cuda->shared_data, & cuda_neigh_list->sneighlist, eflag, vflag, eflag_atom, vflag_atom);

  if(eflag) cuda->cu_eng_vdwl->download();

  if(vflag) cuda->cu_virial->download();
}

/* ---------------------------------------------------------------------- */

void PairEAMCuda::settings(int narg, char** arg)
{
  PairEAM::settings(narg, arg);
  cuda->shared_data.pair.cut_global = (F_CFLOAT) cutforcesq;
}

/* ---------------------------------------------------------------------- */

void PairEAMCuda::coeff(int narg, char** arg)
{
  PairEAM::coeff(narg, arg);
  allocate();
}

void PairEAMCuda::init_style()
{
  MYDBG(printf("# CUDA PairEAMCuda::init_style start\n");)
  // request regular or rRESPA neighbor lists
  file2array();
  array2spline();
  int irequest;


  irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->cudable = 1;

  delete cu_rhor_spline;
  delete cu_z2r_spline;
  delete cu_frho_spline;

  cu_rhor_spline = new cCudaData<double, F_CFLOAT, xyz>((double*)rhor_spline, nrhor, nr + 1, EAM_COEFF_LENGTH);
  cu_z2r_spline = new cCudaData<double, F_CFLOAT, xyz>((double*)z2r_spline, nz2r, nr + 1, EAM_COEFF_LENGTH);
  cu_frho_spline = new cCudaData<double, F_CFLOAT, xyz>((double*)frho_spline, nfrho, nrho + 1, EAM_COEFF_LENGTH);

  cu_rhor_spline->upload();
  cu_z2r_spline->upload();
  cu_frho_spline->upload();

  MYDBG(printf("# CUDA PairEAMCuda::init_style end\n");)
}

void PairEAMCuda::init_list(int id, NeighList* ptr)
{
  MYDBG(printf("# CUDA PairEAMCuda::init_list\n");)
  PairEAM::init_list(id, ptr);

  // right now we can only handle verlet (id 0), not respa
  if(id == 0) cuda_neigh_list = cuda->registerNeighborList(ptr);

  // see Neighbor::init() for details on lammps lists' logic
  MYDBG(printf("# CUDA PairEAMCuda::init_list end\n");)
}

void PairEAMCuda::array2spline()
{
  rdr = 1.0 / dr;
  rdrho = 1.0 / drho;

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);

  memory->create(frho_spline, nfrho, nrho + 1, 8, "pair:frho");
  memory->create(rhor_spline, nrhor, nr + 1, 8, "pair:rhor");
  memory->create(z2r_spline, nz2r, nr + 1, 8, "pair:z2r");

  for(int i = 0; i < nfrho; i++) {
    interpolate(nrho, drho, frho[i], frho_spline[i]);

    for(int j = 0; j < nrho + 1; j++)
      frho_spline[i][j][7] = frho_spline[i][j][3];
  }

  for(int i = 0; i < nrhor; i++) {
    interpolate(nr, dr, rhor[i], rhor_spline[i]);

    for(int j = 0; j < nr + 1; j++)
      rhor_spline[i][j][7] = rhor_spline[i][j][3];
  }

  for(int i = 0; i < nz2r; i++) {
    interpolate(nr, dr, z2r[i], z2r_spline[i]);

    for(int j = 0; j < nr + 1; j++)
      z2r_spline[i][j][7] = z2r_spline[i][j][3];
  }
}

/* ---------------------------------------------------------------------- */

int PairEAMCuda::pack_forward_comm(int n, int* iswap, double* buf, 
                                   int pbc_flag, int* pbc)
{
  Cuda_PairEAMCuda_PackComm(&cuda->shared_data, n, *iswap, buf);

  if(sizeof(F_CFLOAT) < sizeof(double)) return n;
  else return n;
}

/* ---------------------------------------------------------------------- */

void PairEAMCuda::unpack_forward_comm(int n, int first, double* buf)
{
  Cuda_PairEAMCuda_UnpackComm(&cuda->shared_data, n, first, buf, cu_fp->dev_data());
}

void PairEAMCuda::ev_setup(int eflag, int vflag)
{
  int maxeatomold = maxeatom;
  PairEAM::ev_setup(eflag, vflag);

  if(eflag_atom && atom->nmax > maxeatomold) {
    delete cuda->cu_eatom;
    cuda->cu_eatom = new cCudaData<double, ENERGY_CFLOAT, x > ((double*)eatom, & cuda->shared_data.atom.eatom , atom->nmax);
  }

  if(vflag_atom && atom->nmax > maxeatomold) {
    delete cuda->cu_vatom;
    cuda->cu_vatom = new cCudaData<double, ENERGY_CFLOAT, yx > ((double*)vatom, & cuda->shared_data.atom.vatom , atom->nmax, 6);
  }

}
