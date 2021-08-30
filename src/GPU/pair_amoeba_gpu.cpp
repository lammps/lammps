// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#include "pair_amoeba_gpu.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int amoeba_gpu_init(const int ntypes, const int max_amtype,
                    const double *host_pdamp, const double *host_thole,
                    const double *host_special_polar_wscale,
                    const double *host_special_polar_piscale,
                    const double *host_special_polar_pscale,
                    const int nlocal, const int nall, const int max_nbors,
                    const int maxspecial, const int maxspecial15,
                    const double cell_size, int &gpu_mode, FILE *screen,
                    const double aewald, const double felec,
                    const double off2, const double polar_dscale,
                    const double polar_uscale, int& tep_size);
void amoeba_gpu_clear();

int ** amoeba_gpu_compute_n(const int ago, const int inum, const int nall,
                            double **host_x, int *host_type, int *host_amtype, int *host_amgroup,
                            double **host_rpole, double **host_uind, double **host_uinp,
                            double *sublo, double *subhi, tagint *tag, int **nspecial,
                            tagint **special, int* nspecial15, tagint** special15,
                            const bool eflag, const bool vflag,
                            const bool eatom, const bool vatom, int &host_start,
                            int **ilist, int **jnum, const double cpu_time,
                            bool &success, double *host_q, double *boxlo, double *prd,
                            void **tep_ptr);
void amoeba_gpu_compute(const int ago, const int inum,
                        const int nall, double **host_x, int *host_type,
                        int *host_amtype, int *host_amgroup,
                        double **host_rpole, double **host_uind, double **host_uinp,
                        int *ilist, int *numj, int **firstneigh,
                        const bool eflag, const bool vflag, const bool eatom,
                        const bool vatom, int &host_start, const double cpu_time,
                        bool &success, double *host_q, const int nlocal,
                        double *boxlo, double *prd, void **tep_ptr);

double amoeba_gpu_bytes();

enum{VDWL,REPULSE,QFER,DISP,MPOLE,POLAR,USOLV,DISP_LONG,MPOLE_LONG,POLAR_LONG};

/* ---------------------------------------------------------------------- */

PairAmoebaGPU::PairAmoebaGPU(LAMMPS *lmp) : PairAmoeba(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  reinitflag = 0;
  cpu_time = 0.0;
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairAmoebaGPU::~PairAmoebaGPU()
{
  amoeba_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairAmoebaGPU::polar_real()
{
  int eflag=1, vflag=1;
  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    double sublo[3],subhi[3];
    if (domain->triclinic == 0) {
      sublo[0] = domain->sublo[0];
      sublo[1] = domain->sublo[1];
      sublo[2] = domain->sublo[2];
      subhi[0] = domain->subhi[0];
      subhi[1] = domain->subhi[1];
      subhi[2] = domain->subhi[2];
    } else {
      domain->bbox(domain->sublo_lamda,domain->subhi_lamda,sublo,subhi);
    }
    inum = atom->nlocal;

    firstneigh = amoeba_gpu_compute_n(neighbor->ago, inum, nall, atom->x,
                                      atom->type, amtype, amgroup,
                                      rpole, uind, uinp, sublo, subhi,
                                      atom->tag, atom->nspecial, atom->special,
                                      atom->nspecial15, atom->special15,
                                      eflag, vflag, eflag_atom, vflag_atom,
                                      host_start, &ilist, &numneigh, cpu_time,
                                      success, atom->q, domain->boxlo,
                                      domain->prd, &tep_pinned);

  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    
    amoeba_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type,
                       amtype, amgroup, rpole, uind, uinp,
                       ilist, numneigh, firstneigh, eflag, vflag, eflag_atom,
                       vflag_atom, host_start, cpu_time, success, atom->q,
                       atom->nlocal, domain->boxlo, domain->prd, &tep_pinned);
  }
  if (!success)
    error->one(FLERR,"Insufficient memory on accelerator");

  // reference to the tep array from GPU lib

  if (tep_single) {
    float *tep_ptr = (float *)tep_pinned;
    compute_force_from_tep<float>(tep_ptr);
  } else {
    double *tep_ptr = (double *)tep_pinned;
    compute_force_from_tep<double>(tep_ptr);
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template <class numtyp>
void PairAmoebaGPU::compute_force_from_tep(const numtyp* tep_ptr)
{
  int i,ix,iy,iz;
  double ci,dix,diy,diz;
  double qixx,qixy,qixz;
  double qiyy,qiyz,qizz;
  double xix,yix,zix;
  double xiy,yiy,ziy;
  double xiz,yiz,ziz;
  double vxx,vyy,vzz;
  double vxy,vxz,vyz;
  double fix[3],fiy[3],fiz[3],tep[4];

  double** x = atom->x;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    dix = rpole[i][1];
    diy = rpole[i][2];
    diz = rpole[i][3];
    qixx = rpole[i][4];
    qixy = rpole[i][5];
    qixz = rpole[i][6];
    qiyy = rpole[i][8];
    qiyz = rpole[i][9];
    qizz = rpole[i][12];
    
    tep[0] = tep_ptr[4*i];
    tep[1] = tep_ptr[4*i+1];
    tep[2] = tep_ptr[4*i+2];
    torque2force(i,tep,fix,fiy,fiz,fpolar);

    iz = zaxis2local[i];
    ix = xaxis2local[i];
    iy = yaxis2local[i];

    xiz = x[iz][0] - x[i][0];
    yiz = x[iz][1] - x[i][1];
    ziz = x[iz][2] - x[i][2];
    xix = x[ix][0] - x[i][0];
    yix = x[ix][1] - x[i][1];
    zix = x[ix][2] - x[i][2];
    xiy = x[iy][0] - x[i][0];
    yiy = x[iy][1] - x[i][1];
    ziy = x[iy][2] - x[i][2];

    vxx = xix*fix[0] + xiy*fiy[0] + xiz*fiz[0];
    vyy = yix*fix[1] + yiy*fiy[1] + yiz*fiz[1];
    vzz = zix*fix[2] + ziy*fiy[2] + ziz*fiz[2];
    vxy = 0.5 * (yix*fix[0] + yiy*fiy[0] + yiz*fiz[0] + 
                xix*fix[1] + xiy*fiy[1] + xiz*fiz[1]);
    vxz = 0.5 * (zix*fix[0] + ziy*fiy[0] + ziz*fiz[0] + 
                xix*fix[2] + xiy*fiy[2] + xiz*fiz[2]);
    vyz = 0.5 * (zix*fix[1] + ziy*fiy[1] + ziz*fiz[1] + 
                yix*fix[2] + yiy*fiy[2] + yiz*fiz[2]);

    virpolar[0] += vxx;
    virpolar[1] += vyy;
    virpolar[2] += vzz;
    virpolar[3] += vxy;
    virpolar[4] += vxz;
    virpolar[5] += vyz;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairAmoebaGPU::init_style()
{
  PairAmoeba::init_style();

  // Repeat cutsq calculation because done after call to init_style

  double maxcut = -1.0;
  double cut;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        cut = init_one(i,j);
        cut *= cut;
        if (cut > maxcut)
          maxcut = cut;
        cutsq[i][j] = cutsq[j][i] = cut;
      } else
        cutsq[i][j] = cutsq[j][i] = 0.0;
    }
  }

  // select the cutoff (off2) for neighbor list builds (the polar term for now)

  if (use_ewald) choose(POLAR_LONG);
  else choose(POLAR);

  double cell_size = sqrt(off2) + neighbor->skin;

  int maxspecial=0;
  int maxspecial15=0;
  if (atom->molecular != Atom::ATOMIC) {
    maxspecial=atom->maxspecial;
    maxspecial15=atom->maxspecial15;
  }
    
  int tep_size;
  int mnf = 5e-2 * neighbor->oneatom;

  // set the energy unit conversion factor for polar real-space calculation

  double felec = 0.5 * electric / am_dielectric;
  
  int success = amoeba_gpu_init(atom->ntypes+1, max_amtype, pdamp, thole,
                                special_polar_wscale, special_polar_piscale,
                                special_polar_pscale, atom->nlocal,
                                atom->nlocal+atom->nghost, mnf, maxspecial,
                                maxspecial15, cell_size, gpu_mode, screen,
                                aewald, felec, off2, polar_dscale, polar_uscale,
                                tep_size);
  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE)
    error->all(FLERR,"Pair style amoeba/gpu does not support neigh no for now");

  if (tep_size == sizeof(double))
    tep_single = false;
  else
    tep_single = true;
}

/* ---------------------------------------------------------------------- */

double PairAmoebaGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + amoeba_gpu_bytes();
}
