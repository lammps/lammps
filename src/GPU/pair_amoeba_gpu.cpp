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
#include "math_const.h"
#include "my_page.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "suffix.h"
#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

enum{MUTUAL,OPT,TCG,DIRECT};

// External functions from cuda library for atom decomposition

int amoeba_gpu_init(const int ntypes, const int max_amtype,
                    const double *host_pdamp, const double *host_thole,
                    const double *host_dirdamp,
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

int ** amoeba_gpu_compute_udirect2b(const int ago, const int inum, const int nall,
                            double **host_x, int *host_type, int *host_amtype, int *host_amgroup,
                            double **host_rpole, double *sublo, double *subhi, tagint *tag, int **nspecial,
                            tagint **special, int* nspecial15, tagint** special15,
                            const bool eflag, const bool vflag,
                            const bool eatom, const bool vatom, int &host_start,
                            int **ilist, int **jnum, const double cpu_time,
                            bool &success, double *host_q, double *boxlo, double *prd,
                            void **fieldp_ptr);

int ** amoeba_gpu_compute_polar_real(const int ago, const int inum, const int nall,
                            double **host_x, int *host_type, int *host_amtype, int *host_amgroup,
                            double **host_rpole, double **host_uind, double **host_uinp,
                            double *sublo, double *subhi, tagint *tag, int **nspecial,
                            tagint **special, int* nspecial15, tagint** special15,
                            const bool eflag, const bool vflag,
                            const bool eatom, const bool vatom, int &host_start,
                            int **ilist, int **jnum, const double cpu_time,
                            bool &success, double *host_q, double *boxlo, double *prd,
                            void **tep_ptr);

double amoeba_gpu_bytes();

enum{VDWL,REPULSE,QFER,DISP,MPOLE,POLAR,USOLV,DISP_LONG,MPOLE_LONG,POLAR_LONG};

/* ---------------------------------------------------------------------- */

PairAmoebaGPU::PairAmoebaGPU(LAMMPS *lmp) : PairAmoeba(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  reinitflag = 0;
  cpu_time = 0.0;
  suffix_flag |= Suffix::GPU;
  fieldp_pinned = nullptr;
  tep_pinned = nullptr;
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
  bool gpu_polar_real_ready = true;
  if (!gpu_polar_real_ready) {
    PairAmoeba::polar_real();
    return;
  }

  int eflag=1, vflag=1;
  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  
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

  firstneigh = amoeba_gpu_compute_polar_real(neighbor->ago, inum, nall, atom->x,
                                        atom->type, amtype, amgroup,
                                        rpole, uind, uinp, sublo, subhi,
                                        atom->tag, atom->nspecial, atom->special,
                                        atom->nspecial15, atom->special15,
                                        eflag, vflag, eflag_atom, vflag_atom,
                                        host_start, &ilist, &numneigh, cpu_time,
                                        success, atom->q, domain->boxlo,
                                        domain->prd, &tep_pinned);

  
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
  // NOTE: induce and polar terms are using the same flags here

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
  
  int success = amoeba_gpu_init(atom->ntypes+1, max_amtype, pdamp, thole, dirdamp,
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

/* ----------------------------------------------------------------------
   udirect2b = Ewald real direct field via list
   udirect2b computes the real space contribution of the permanent
   atomic multipole moments to the field via a neighbor list
------------------------------------------------------------------------- */

void PairAmoebaGPU::udirect2b(double **field, double **fieldp)
{
  bool gpu_udirect2b_ready = true;
  if (!gpu_udirect2b_ready) {
    PairAmoeba::udirect2b(field, fieldp);
    return;
  }
   
  int eflag=1, vflag=1;
  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  
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

  firstneigh = amoeba_gpu_compute_udirect2b(neighbor->ago, inum, nall, atom->x,
                                        atom->type, amtype, amgroup, rpole, sublo,
                                        subhi, atom->tag, atom->nspecial, atom->special,
                                        atom->nspecial15, atom->special15,
                                        eflag, vflag, eflag_atom, vflag_atom,
                                        host_start, &ilist, &numneigh, cpu_time,
                                        success, atom->q, domain->boxlo,
                                        domain->prd, &fieldp_pinned);
  if (!success)
    error->one(FLERR,"Insufficient memory on accelerator");

  // get field and fieldp values from the GPU lib

  int nlocal = atom->nlocal;
  double *field_ptr = (double *)fieldp_pinned;

  for (int i = 0; i < nlocal; i++) {
    int idx = 4*i;
    field[i][0] = field_ptr[idx];
    field[i][1] = field_ptr[idx+1];
    field[i][2] = field_ptr[idx+2]; 
  }

  double* fieldp_ptr = (double *)fieldp_pinned;
  fieldp_ptr += 4*inum;
  for (int i = 0; i < nlocal; i++) {
    int idx = 4*i;
    fieldp[i][0] = fieldp_ptr[idx];
    fieldp[i][1] = fieldp_ptr[idx+1];
    fieldp[i][2] = fieldp_ptr[idx+2];
  }

  // rebuild dipole-dipole pair list and store pairwise dipole matrices
  // done one atom at a time in real-space double loop over atoms & neighs

  udirect2b_cpu();
}

/* ----------------------------------------------------------------------
   udirect2b = Ewald real direct field via list
   udirect2b computes the real space contribution of the permanent
     atomic multipole moments to the field via a neighbor list
------------------------------------------------------------------------- */

void PairAmoebaGPU::udirect2b_cpu()
{
  int i,j,k,m,n,ii,jj,kk,kkk,jextra,ndip,itype,jtype,igroup,jgroup;
  double xr,yr,zr,r,r2;
  double rr1,rr2,rr3,rr5;
  double bfac,exp2a;
  double ralpha,aefac;
  double aesq2,aesq2n;
  double pdi,pti,ddi;
  double pgamma;
  double damp,expdamp;
  double scale3,scale5;
  double scale7,scalek;
  double bn[4],bcn[3];
  double factor_dscale,factor_pscale,factor_uscale,factor_wscale;

  int inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;

  // neigh list

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // NOTE: doesn't this have a problem if aewald is tiny ??
  
  aesq2 = 2.0 * aewald * aewald;
  aesq2n = 0.0;
  if (aewald > 0.0) aesq2n = 1.0 / (MY_PIS*aewald);

  // rebuild dipole-dipole pair list and store pairwise dipole matrices
  // done one atom at a time in real-space double loop over atoms & neighs

  int *neighptr;
  double *tdipdip;

  // compute the real space portion of the Ewald summation

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = amtype[i];
    igroup = amgroup[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    n = ndip = 0;
    neighptr = ipage_dipole->vget();
    tdipdip = dpage_dipdip->vget();

    pdi = pdamp[itype];
    pti = thole[itype];
    ddi = dirdamp[itype];
    
    // evaluate all sites within the cutoff distance

    for (jj = 0; jj < jnum; jj++) {
      jextra = jlist[jj];
      j = jextra & NEIGHMASK15;
      
      xr = x[j][0] - x[i][0];
      yr = x[j][1] - x[i][1];
      zr = x[j][2] - x[i][2];
      r2 = xr*xr + yr* yr + zr*zr;
      if (r2 > off2) continue;

      jtype = amtype[j];
      jgroup = amgroup[j];
      
      factor_wscale = special_polar_wscale[sbmask15(jextra)];
      if (igroup == jgroup) {
        factor_pscale = special_polar_piscale[sbmask15(jextra)];
        factor_dscale = polar_dscale;
        factor_uscale = polar_uscale;
      } else {
        factor_pscale = special_polar_pscale[sbmask15(jextra)];
        factor_dscale = factor_uscale = 1.0;
      }

      r = sqrt(r2);
      rr1 = 1.0 / r;
      rr2 = rr1 * rr1;
      rr3 = rr2 * rr1;
      rr5 = 3.0 * rr2 * rr3;

      // calculate the real space Ewald error function terms

      ralpha = aewald * r;
      bn[0] = erfc(ralpha) * rr1;
      exp2a = exp(-ralpha*ralpha);
      aefac = aesq2n;
      for (m = 1; m <= 3; m++) {
        bfac = m+m-1;
        aefac = aesq2 * aefac;
        bn[m] = (bfac*bn[m-1]+aefac*exp2a) * rr2;
      }
      
      // find terms needed later to compute mutual polarization

      if (poltyp != DIRECT) {
        scale3 = 1.0;
        scale5 = 1.0;
        damp = pdi * pdamp[jtype];
        if (damp != 0.0) {
          pgamma = MIN(pti,thole[jtype]);
          damp = pgamma * pow(r/damp,3.0);
          if (damp < 50.0) {
            expdamp = exp(-damp);
            scale3 = 1.0 - expdamp;
            scale5 = 1.0 - expdamp*(1.0+damp);
          }
        }
        scalek = factor_uscale;
        bcn[0] = bn[1] - (1.0-scalek*scale3)*rr3;
        bcn[1] = bn[2] - (1.0-scalek*scale5)*rr5;
        
        neighptr[n++] = j;
        tdipdip[ndip++] = -bcn[0] + bcn[1]*xr*xr;
        tdipdip[ndip++] = bcn[1]*xr*yr;
        tdipdip[ndip++] = bcn[1]*xr*zr;
        tdipdip[ndip++] = -bcn[0] + bcn[1]*yr*yr;
        tdipdip[ndip++] = bcn[1]*yr*zr;
        tdipdip[ndip++] = -bcn[0] + bcn[1]*zr*zr;
      }
      
    } // jj

    firstneigh_dipole[i] = neighptr;
    firstneigh_dipdip[i] = tdipdip;
    numneigh_dipole[i] = n;
    ipage_dipole->vgot(n);
    dpage_dipdip->vgot(ndip);
  }
}

/* ---------------------------------------------------------------------- */

double PairAmoebaGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + amoeba_gpu_bytes();
}
