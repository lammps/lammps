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
   Contributing authors: Mike Brown (ORNL), Axel Kohlmeyer (Temple)
------------------------------------------------------------------------- */

#include "pppm_gpu.h"
#include <mpi.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "atom.h"
#include "comm.h"
#include "gridcomm.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "domain.h"
#include "fft3d_wrap.h"
#include "remap_wrap.h"
#include "gpu_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "universe.h"
#include "fix.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXORDER 7
#define OFFSET 16384
#define SMALL 0.00001
#define LARGE 10000.0
#define EPS_HOC 1.0e-7

enum{REVERSE_RHO_GPU,REVERSE_RHO};
enum{FORWARD_IK,FORWARD_AD,FORWARD_IK_PERATOM,FORWARD_AD_PERATOM};

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

// external functions from cuda library for atom decomposition

#ifdef FFT_SINGLE
#define PPPM_GPU_API(api)  pppm_gpu_ ## api ## _f
#else
#define PPPM_GPU_API(api)  pppm_gpu_ ## api ## _d
#endif

FFT_SCALAR* PPPM_GPU_API(init)(const int nlocal, const int nall, FILE *screen,
                               const int order, const int nxlo_out,
                               const int nylo_out, const int nzlo_out,
                               const int nxhi_out, const int nyhi_out,
                               const int nzhi_out, FFT_SCALAR **rho_coeff,
                               FFT_SCALAR **_vd_brick,
                               const double slab_volfactor,
                               const int nx_pppm, const int ny_pppm,
                               const int nz_pppm, const bool split,
                               const bool respa, int &success);
void PPPM_GPU_API(clear)(const double poisson_time);
int PPPM_GPU_API(spread)(const int ago, const int nlocal, const int nall,
                      double **host_x, int *host_type, bool &success,
                      double *host_q, double *boxlo, const double delxinv,
                      const double delyinv, const double delzinv);
void PPPM_GPU_API(interp)(const FFT_SCALAR qqrd2e_scale);
double PPPM_GPU_API(bytes)();
void PPPM_GPU_API(forces)(double **f);

/* ---------------------------------------------------------------------- */

PPPMGPU::PPPMGPU(LAMMPS *lmp) : PPPM(lmp)
{
  triclinic_support = 0;
  density_brick_gpu = vd_brick = NULL;
  kspace_split = false;
  im_real_space = false;

  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

PPPMGPU::~PPPMGPU()
{
  PPPM_GPU_API(clear)(poisson_time);
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

void PPPMGPU::init()
{
  // PPPM init manages all arrays except density_brick_gpu and vd_brick
  //      thru its deallocate(), allocate()
  // NOTE: could free density_brick and vdxyz_brick after PPPM allocates them,
  //       before allocating db_gpu and vd_brick down below, if don't need,
  //       if do this, make sure to set them to NULL

  destroy_3d_offset(density_brick_gpu,nzlo_out,nylo_out);
  destroy_3d_offset(vd_brick,nzlo_out,nylo_out);
  density_brick_gpu = vd_brick = NULL;

  PPPM::init();

  // insure no conflict with fix balance

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"balance") == 0)
      error->all(FLERR,"Cannot currently use pppm/gpu with fix balance.");

  // unsupported option

  if (differentiation_flag == 1)
    error->all(FLERR,"Cannot (yet) do analytic differentiation with pppm/gpu");

  if (strcmp(update->integrate_style,"verlet/split") == 0) {
    kspace_split=true;
    old_nlocal = 0;
  }

  if (kspace_split && universe->iworld == 0) {
    im_real_space = true;
    return;
  }

  // GPU precision specific init

  bool respa_value=false;
  if (strstr(update->integrate_style,"respa"))
    respa_value=true;

  if (order>8)
    error->all(FLERR,"Cannot use order greater than 8 with pppm/gpu.");
  PPPM_GPU_API(clear)(poisson_time);

  int success;
  FFT_SCALAR *data, *h_brick;
  h_brick = PPPM_GPU_API(init)(atom->nlocal, atom->nlocal+atom->nghost, screen,
                               order, nxlo_out, nylo_out, nzlo_out, nxhi_out,
                               nyhi_out, nzhi_out, rho_coeff, &data,
                               slab_volfactor,nx_pppm,ny_pppm,nz_pppm,
                               kspace_split,respa_value,success);

  GPU_EXTRA::check_flag(success,error,world);

  // allocate density_brick_gpu and vd_brick

  density_brick_gpu =
    create_3d_offset(nzlo_out,nzhi_out,nylo_out,nyhi_out,
                     nxlo_out,nxhi_out,"pppm:density_brick_gpu",h_brick,1);
  vd_brick =
    create_3d_offset(nzlo_out,nzhi_out,nylo_out,nyhi_out,
                     nxlo_out,nxhi_out,"pppm:vd_brick",data,4);

  poisson_time = 0.0;
}

/* ----------------------------------------------------------------------
   compute the PPPMGPU long-range force, energy, virial
------------------------------------------------------------------------- */

void PPPMGPU::compute(int eflag, int vflag)
{
  int i,j;

  int nago;
  if (kspace_split) {
    if (im_real_space) return;
    if (atom->nlocal > old_nlocal) {
      nago=0;
      old_nlocal = atom->nlocal;
    } else nago = 1;
  } else nago = neighbor->ago;

  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  ev_init(eflag,vflag);

  // If need per-atom energies/virials, allocate per-atom arrays here
  // so that particle map on host can be done concurrently with GPU calculations

  if (evflag_atom && !peratom_allocate_flag) {
    allocate_peratom();
    cg_peratom->ghost_notify();
    cg_peratom->setup();
  }

  bool success = true;
  int flag=PPPM_GPU_API(spread)(nago, atom->nlocal, atom->nlocal +
                             atom->nghost, atom->x, atom->type, success,
                             atom->q, domain->boxlo, delxinv, delyinv,
                             delzinv);
  if (!success)
    error->one(FLERR,"Insufficient memory on accelerator");
  if (flag != 0)
    error->one(FLERR,"Out of range atoms - cannot compute PPPM");

  // convert atoms from box to lamda coords

  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  // If need per-atom energies/virials, also do particle map on host
  // concurrently with GPU calculations

  if (evflag_atom) {

    // extend size of per-atom arrays if necessary

    if (atom->nmax > nmax) {
      memory->destroy(part2grid);
      nmax = atom->nmax;
      memory->create(part2grid,nmax,3,"pppm:part2grid");
    }

    particle_map();
  }

  double t3 = MPI_Wtime();

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  cg->reverse_comm(this,REVERSE_RHO_GPU);
  brick2fft_gpu();

  // compute potential gradient on my FFT grid and
  //   portion of e_long on this proc's FFT grid
  // return gradients (electric fields) in 3d brick decomposition

  poisson();

  // all procs communicate E-field values
  // to fill ghost cells surrounding their 3d bricks

  if (differentiation_flag == 1) cg->forward_comm(this,FORWARD_AD);
  else cg->forward_comm(this,FORWARD_IK);

  // extra per-atom energy/virial communication

  if (evflag_atom) {
    if (differentiation_flag == 1 && vflag_atom)
      cg_peratom->forward_comm(this,FORWARD_AD_PERATOM);
    else if (differentiation_flag == 0)
      cg_peratom->forward_comm(this,FORWARD_IK_PERATOM);
  }

  poisson_time += MPI_Wtime()-t3;

  // calculate the force on my particles

  FFT_SCALAR qscale = force->qqrd2e * scale;
  PPPM_GPU_API(interp)(qscale);

  // per-atom energy/virial
  // energy includes self-energy correction

  if (evflag_atom) fieldforce_peratom();

  // update qsum and qsqsum, if atom count has changed and energy needed

  if ((eflag_global || eflag_atom) && atom->natoms != natoms_original) {
    qsum_qsq();
    natoms_original = atom->natoms;
  }

  // sum energy across procs and add in volume-dependent term

  if (eflag_global) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;

    energy *= 0.5*volume;
    energy -= g_ewald*qsqsum/MY_PIS +
      MY_PI2*qsum*qsum / (g_ewald*g_ewald*volume);
    energy *= qscale;
  }

  // sum virial across procs

  if (vflag_global) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*qscale*volume*virial_all[i];
  }

  // per-atom energy/virial
  // energy includes self-energy correction

  if (evflag_atom) {
    double *q = atom->q;
    int nlocal = atom->nlocal;

    if (eflag_atom) {
      for (i = 0; i < nlocal; i++) {
        eatom[i] *= 0.5;
        eatom[i] -= g_ewald*q[i]*q[i]/MY_PIS + MY_PI2*q[i]*qsum /
          (g_ewald*g_ewald*volume);
        eatom[i] *= qscale;
      }
    }

    if (vflag_atom) {
      for (i = 0; i < nlocal; i++)
        for (j = 0; j < 6; j++) vatom[i][j] *= 0.5*qscale;
    }
  }

  // 2d slab correction

  if (slabflag) slabcorr();

  // convert atoms back from lamda to box coords

  if (triclinic) domain->lamda2x(atom->nlocal);

  if (kspace_split) PPPM_GPU_API(forces)(atom->f);
}

/* ----------------------------------------------------------------------
   remap density from 3d brick decomposition to FFT decomposition
------------------------------------------------------------------------- */

void PPPMGPU::brick2fft_gpu()
{
  int n,ix,iy,iz;

  // copy grabs inner portion of density from 3d brick
  // remap could be done as pre-stage of FFT,
  //   but this works optimally on only double values, not complex values

  n = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        density_fft[n++] = density_brick_gpu[iz][iy][ix];

  remap->perform(density_fft,density_fft,work1);
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver
------------------------------------------------------------------------- */

void PPPMGPU::poisson_ik()
{
  int i,j,k,n;
  double eng;

  // transform charge density (r -> k)

  n = 0;
  for (i = 0; i < nfft; i++) {
    work1[n++] = density_fft[i];
    work1[n++] = ZEROF;
  }

  fft1->compute(work1,work1,1);

  // if requested, compute energy and virial contribution

  double scaleinv = 1.0/(nx_pppm*ny_pppm*nz_pppm);
  double s2 = scaleinv*scaleinv;

  if (eflag_global || vflag_global) {
    if (vflag_global) {
      n = 0;
      for (i = 0; i < nfft; i++) {
        eng = s2 * greensfn[i] * (work1[n]*work1[n] + work1[n+1]*work1[n+1]);
        for (j = 0; j < 6; j++) virial[j] += eng*vg[i][j];
        if (eflag_global) energy += eng;
        n += 2;
      }
    } else {
      n = 0;
      for (i = 0; i < nfft; i++) {
        energy +=
          s2 * greensfn[i] * (work1[n]*work1[n] + work1[n+1]*work1[n+1]);
        n += 2;
      }
    }
  }

  // scale by 1/total-grid-pts to get rho(k)
  // multiply by Green's function to get V(k)

  n = 0;
  for (i = 0; i < nfft; i++) {
    work1[n++] *= scaleinv * greensfn[i];
    work1[n++] *= scaleinv * greensfn[i];
  }

  // extra FFTs for per-atom energy/virial

  if (evflag_atom) poisson_peratom();

  // compute gradients of V(r) in each of 3 dims by transformimg -ik*V(k)
  // FFT leaves data in 3d brick decomposition
  // copy it into inner portion of vdx,vdy,vdz arrays

  // x direction gradient

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work2[n] = fkx[i]*work1[n+1];
        work2[n+1] = -fkx[i]*work1[n];
        n += 2;
      }

  fft2->compute(work2,work2,-1);

  n = 0;
  int x_hi = nxhi_in * 4 + 3;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in * 4; i < x_hi; i+=4) {
        vd_brick[k][j][i] = work2[n];
        n += 2;
      }

  // y direction gradient

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work2[n] = fky[j]*work1[n+1];
        work2[n+1] = -fky[j]*work1[n];
        n += 2;
      }

  fft2->compute(work2,work2,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in * 4 + 1; i < x_hi; i+=4) {
        vd_brick[k][j][i] = work2[n];
        n += 2;
      }

  // z direction gradient

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work2[n] = fkz[k]*work1[n+1];
        work2[n+1] = -fkz[k]*work1[n];
        n += 2;
      }

  fft2->compute(work2,work2,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in * 4 + 2; i < x_hi; i+=4) {
        vd_brick[k][j][i] = work2[n];
        n += 2;
      }
}

/* ----------------------------------------------------------------------
   pack own values to buf to send to another proc
------------------------------------------------------------------------- */

void PPPMGPU::pack_forward(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;

  if (flag == FORWARD_IK) {
    int offset;
    FFT_SCALAR *src = &vd_brick[nzlo_out][nylo_out][4*nxlo_out];
    for (int i = 0; i < nlist; i++) {
      offset = 4*list[i];
      buf[n++] = src[offset++];
      buf[n++] = src[offset++];
      buf[n++] = src[offset];
    }
  } else if (flag == FORWARD_AD) {
    FFT_SCALAR *src = &u_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (flag == FORWARD_IK_PERATOM) {
    FFT_SCALAR *esrc = &u_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v0src = &v0_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1src = &v1_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2src = &v2_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3src = &v3_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4src = &v4_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5src = &v5_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      if (eflag_atom) buf[n++] = esrc[list[i]];
      if (vflag_atom) {
        buf[n++] = v0src[list[i]];
        buf[n++] = v1src[list[i]];
        buf[n++] = v2src[list[i]];
        buf[n++] = v3src[list[i]];
        buf[n++] = v4src[list[i]];
        buf[n++] = v5src[list[i]];
      }
    }
  } else if (flag == FORWARD_AD_PERATOM) {
    FFT_SCALAR *v0src = &v0_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1src = &v1_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2src = &v2_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3src = &v3_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4src = &v4_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5src = &v5_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      buf[n++] = v0src[list[i]];
      buf[n++] = v1src[list[i]];
      buf[n++] = v2src[list[i]];
      buf[n++] = v3src[list[i]];
      buf[n++] = v4src[list[i]];
      buf[n++] = v5src[list[i]];
    }
  }
}

/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void PPPMGPU::unpack_forward(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;

  if (flag == FORWARD_IK) {
    int offset;
    FFT_SCALAR *dest = &vd_brick[nzlo_out][nylo_out][4*nxlo_out];
    for (int i = 0; i < nlist; i++) {
      offset = 4*list[i];
      dest[offset++] = buf[n++];
      dest[offset++] = buf[n++];
      dest[offset] = buf[n++];
    }
  } else if (flag == FORWARD_AD) {
    FFT_SCALAR *dest = &u_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      dest[list[i]] = buf[i];
  } else if (flag == FORWARD_IK_PERATOM) {
    FFT_SCALAR *esrc = &u_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v0src = &v0_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1src = &v1_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2src = &v2_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3src = &v3_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4src = &v4_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5src = &v5_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      if (eflag_atom) esrc[list[i]] = buf[n++];
      if (vflag_atom) {
        v0src[list[i]] = buf[n++];
        v1src[list[i]] = buf[n++];
        v2src[list[i]] = buf[n++];
        v3src[list[i]] = buf[n++];
        v4src[list[i]] = buf[n++];
        v5src[list[i]] = buf[n++];
      }
    }
  } else if (flag == FORWARD_AD_PERATOM) {
    FFT_SCALAR *v0src = &v0_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1src = &v1_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2src = &v2_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3src = &v3_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4src = &v4_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5src = &v5_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      v0src[list[i]] = buf[n++];
      v1src[list[i]] = buf[n++];
      v2src[list[i]] = buf[n++];
      v3src[list[i]] = buf[n++];
      v4src[list[i]] = buf[n++];
      v5src[list[i]] = buf[n++];
    }
  }
}

/* ----------------------------------------------------------------------
   pack ghost values into buf to send to another proc
------------------------------------------------------------------------- */

void PPPMGPU::pack_reverse(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  if (flag == REVERSE_RHO_GPU) {
    FFT_SCALAR *src = &density_brick_gpu[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (flag == REVERSE_RHO) {
    FFT_SCALAR *src = &density_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  }
}

/* ----------------------------------------------------------------------
   unpack another proc's ghost values from buf and add to own values
------------------------------------------------------------------------- */

void PPPMGPU::unpack_reverse(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  if (flag == REVERSE_RHO_GPU) {
    FFT_SCALAR *dest = &density_brick_gpu[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      dest[list[i]] += buf[i];
  } else if (flag == REVERSE_RHO) {
    FFT_SCALAR *dest = &density_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      dest[list[i]] += buf[i];
  }
}

/* ----------------------------------------------------------------------
   create array using offsets from pinned memory allocation
------------------------------------------------------------------------- */

FFT_SCALAR ***PPPMGPU::create_3d_offset(int n1lo, int n1hi, int n2lo, int n2hi,
                                        int n3lo, int n3hi, const char *name,
                                        FFT_SCALAR *data, int vec_length)
{
  int i,j;
  int n1 = n1hi - n1lo + 1;
  int n2 = n2hi - n2lo + 1;
  int n3 = n3hi - n3lo + 1;

  FFT_SCALAR **plane = (FFT_SCALAR **)
    memory->smalloc(n1*n2*sizeof(FFT_SCALAR *),name);
  FFT_SCALAR ***array = (FFT_SCALAR ***)
    memory->smalloc(n1*sizeof(FFT_SCALAR **),name);

  int n = 0;
  for (i = 0; i < n1; i++) {
    array[i] = &plane[i*n2];
    for (j = 0; j < n2; j++) {
      plane[i*n2+j] = &data[n];
      n += n3*vec_length;
    }
  }

  for (i = 0; i < n1*n2; i++) array[0][i] -= n3lo*vec_length;
  for (i = 0; i < n1; i++) array[i] -= n2lo;
  return array-n1lo;
}

/* ----------------------------------------------------------------------
   3d memory offsets
------------------------------------------------------------------------- */

void PPPMGPU::destroy_3d_offset(FFT_SCALAR ***array, int n1_offset,
                                 int n2_offset)
{
  if (array == NULL) return;
  memory->sfree(&array[n1_offset][n2_offset]);
  memory->sfree(array + n1_offset);
}


/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

double PPPMGPU::memory_usage()
{
  double bytes = PPPM::memory_usage();

  // NOTE: add tallying here for density_brick_gpu and vd_brick
  //       could subtract out density_brick and vdxyz_brick if freed them above
  //       it the net efffect is zero, do nothing

  return bytes + PPPM_GPU_API(bytes)();
}

/* ----------------------------------------------------------------------
   perform and time the 1d FFTs required for N timesteps
------------------------------------------------------------------------- */

int PPPMGPU::timing_1d(int n, double &time1d)
{
  if (im_real_space) {
    time1d = 1.0;
    return 4;
  }
  PPPM::timing_1d(n,time1d);
  return 4;
}

/* ----------------------------------------------------------------------
   perform and time the 3d FFTs required for N timesteps
------------------------------------------------------------------------- */

int PPPMGPU::timing_3d(int n, double &time3d)
{
  if (im_real_space) {
    time3d = 1.0;
    return 4;
  }
  PPPM::timing_3d(n,time3d);
  return 4;
}

/* ----------------------------------------------------------------------
   adjust PPPM coeffs, called initially and whenever volume has changed
------------------------------------------------------------------------- */

void PPPMGPU::setup()
{
  if (im_real_space) return;
  PPPM::setup();
}

/* ----------------------------------------------------------------------
   group-group interactions
 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   compute the PPPM total long-range force and energy for groups A and B
 ------------------------------------------------------------------------- */

void PPPMGPU::compute_group_group(int groupbit_A, int groupbit_B, int AA_flag)
{
  if (slabflag && triclinic)
    error->all(FLERR,"Cannot (yet) use K-space slab "
               "correction with compute group/group for triclinic systems");

  if (differentiation_flag)
    error->all(FLERR,"Cannot (yet) use kspace_modify "
               "diff ad with compute group/group");

  if (!group_allocate_flag) allocate_groups();

  // convert atoms from box to lamda coords

  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  // extend size of per-atom arrays if necessary
  // part2grid needs to be allocated

  if (atom->nmax > nmax || part2grid == NULL) {
    memory->destroy(part2grid);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"pppm:part2grid");
  }

  particle_map();

  e2group = 0.0; //energy
  f2group[0] = 0.0; //force in x-direction
  f2group[1] = 0.0; //force in y-direction
  f2group[2] = 0.0; //force in z-direction

  // map my particle charge onto my local 3d density grid

  make_rho_groups(groupbit_A,groupbit_B,AA_flag);

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  // temporarily store and switch pointers so we can
  //  use brick2fft() for groups A and B (without
  //  writing an additional function)

  FFT_SCALAR ***density_brick_real = density_brick;
  FFT_SCALAR *density_fft_real = density_fft;

  // group A

  density_brick = density_A_brick;
  density_fft = density_A_fft;

  cg->reverse_comm(this,REVERSE_RHO);
  brick2fft();

  // group B

  density_brick = density_B_brick;
  density_fft = density_B_fft;

  cg->reverse_comm(this,REVERSE_RHO);
  brick2fft();

  // switch back pointers

  density_brick = density_brick_real;
  density_fft = density_fft_real;

  // compute potential gradient on my FFT grid and
  //   portion of group-group energy/force on this proc's FFT grid

  poisson_groups(AA_flag);

  const double qscale = qqrd2e * scale;

  // total group A <--> group B energy
  // self and boundary correction terms are in compute_group_group.cpp

  double e2group_all;
  MPI_Allreduce(&e2group,&e2group_all,1,MPI_DOUBLE,MPI_SUM,world);
  e2group = e2group_all;

  e2group *= qscale*0.5*volume;

  // total group A <--> group B force

  double f2group_all[3];
  MPI_Allreduce(f2group,f2group_all,3,MPI_DOUBLE,MPI_SUM,world);

  f2group[0] = qscale*volume*f2group_all[0];
  f2group[1] = qscale*volume*f2group_all[1];
  if (slabflag != 2) f2group[2] = qscale*volume*f2group_all[2];

  // convert atoms back from lamda to box coords

  if (triclinic) domain->lamda2x(atom->nlocal);

  if (slabflag == 1)
    slabcorr_groups(groupbit_A, groupbit_B, AA_flag);
}

