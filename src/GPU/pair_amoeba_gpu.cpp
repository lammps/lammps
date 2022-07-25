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

#include "amoeba_convolution.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_store.h"
#include "force.h"
#include "gpu_extra.h"
#include "math_const.h"
#include "memory.h"
#include "my_page.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "suffix.h"
#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

enum{INDUCE,RSD,SETUP_AMOEBA,SETUP_HIPPO,KMPOLE,AMGROUP};   // forward comm
enum{FIELD,ZRSD,TORQUE,UFLD};                               // reverse comm
enum{VDWL,REPULSE,QFER,DISP,MPOLE,POLAR,USOLV,DISP_LONG,MPOLE_LONG,POLAR_LONG};
enum{MUTUAL,OPT,TCG,DIRECT};
enum{GEAR,ASPC,LSQR};
enum{BUILD,APPLY};
enum{GORDON1,GORDON2};

#define DEBYE 4.80321    // conversion factor from q-Angs (real units) to Debye

// External functions from cuda library for atom decomposition

int amoeba_gpu_init(const int ntypes, const int max_amtype, const int max_amclass,
                    const double *host_pdamp, const double *host_thole,
                    const double *host_dirdamp, const int* host_amtype2class,
                    const double *host_special_hal, const double *host_special_repel,
                    const double *host_special_disp, const double *host_special_mpole,
                    const double *host_special_polar_wscale,
                    const double *host_special_polar_piscale,
                    const double *host_special_polar_pscale,
                    const double *host_csix, const double *host_adisp,
                    const int nlocal, const int nall, const int max_nbors,
                    const int maxspecial, const int maxspecial15,
                    const double cell_size, int &gpu_mode, FILE *screen,
                    const double polar_dscale, const double polar_uscale, int& tq_size);
void amoeba_gpu_clear();

int ** amoeba_gpu_compute_multipole_real(const int ago, const int inum, const int nall,
              double **host_x, int *host_type, int *host_amtype, int *host_amgroup,
              double **host_rpole, double *sublo, double *subhi, tagint *tag,
              int **nspecial, tagint **special, int* nspecial15, tagint** special15,
              const bool eflag, const bool vflag, const bool eatom, const bool vatom,
              int &host_start, int **ilist, int **jnum, const double cpu_time,
              bool &success, const double aewald, const double felec, const double off2,
              double *host_q, double *boxlo, double *prd, void **tq_ptr);

void amoeba_gpu_compute_udirect2b(int *host_amtype, int *host_amgroup,
              double **host_rpole, double **host_uind, double **host_uinp,
              const double aewald, const double off2, void **fieldp_ptr);

void amoeba_gpu_compute_umutual2b(int *host_amtype, int *host_amgroup,
              double **host_rpole, double **host_uind, double **host_uinp,
              const double aewald, const double off2, void **fieldp_ptr);

void amoeba_gpu_update_fieldp(void **fieldp_ptr);

void amoeba_gpu_compute_polar_real(int *host_amtype, int *host_amgroup,
              double **host_rpole, double **host_uind, double **host_uinp,
              const bool eflag, const bool vflag, const bool eatom, const bool vatom,
              const double aewald, const double felec, const double off2,
              void **tq_ptr);

double amoeba_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairAmoebaGPU::PairAmoebaGPU(LAMMPS *lmp) : PairAmoeba(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  reinitflag = 0;
  cpu_time = 0.0;
  suffix_flag |= Suffix::GPU;
  fieldp_pinned = nullptr;
  tq_pinned = nullptr;

  gpu_hal_ready = false;               // true for AMOEBA when ready
  gpu_repulsion_ready = false;         // always false for AMOEBA
  gpu_dispersion_real_ready = false;   // always false for AMOEBA
  gpu_multipole_real_ready = true;     // need to be true for precompute()
  gpu_udirect2b_ready = true;
  gpu_umutual2b_ready = true;
  gpu_polar_real_ready = true;         // need to be true for copying data from device back to host

  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairAmoebaGPU::~PairAmoebaGPU()
{
  amoeba_gpu_clear();
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

  double cell_size = sqrt(maxcut) + neighbor->skin;

  int maxspecial=0;
  int maxspecial15=0;
  if (atom->molecular != Atom::ATOMIC) {
    maxspecial=atom->maxspecial;
    maxspecial15=atom->maxspecial15;
  }

  int tq_size;
  int mnf = 5e-2 * neighbor->oneatom;
  int success = amoeba_gpu_init(atom->ntypes+1, max_amtype, max_amclass,
                                pdamp, thole, dirdamp, amtype2class, special_hal,
                                special_repel, special_disp, special_mpole,
                                special_polar_wscale, special_polar_piscale,
                                special_polar_pscale, csix, adisp, atom->nlocal,
                                atom->nlocal+atom->nghost, mnf, maxspecial,
                                maxspecial15, cell_size, gpu_mode, screen,
                                polar_dscale, polar_uscale, tq_size);
  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE)
    error->all(FLERR,"Pair style amoeba/gpu does not support neigh no for now");

  if (tq_size == sizeof(double))
    tq_single = false;
  else
    tq_single = true;
}

/* ----------------------------------------------------------------------
   multipole_real = real-space portion of mulipole interactions
   adapted from Tinker emreal1d() routine
------------------------------------------------------------------------- */

void PairAmoebaGPU::multipole_real()
{
  if (!gpu_multipole_real_ready) {
    PairAmoeba::multipole_real();
    return;
  }

  int eflag=1, vflag=1;
  double **f = atom->f;
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

  // select the correct cutoff for the term

  if (use_ewald) choose(MPOLE_LONG);
  else choose(MPOLE);

  // set the energy unit conversion factor for multipolar real-space calculation

  double felec = electric / am_dielectric;

  firstneigh = amoeba_gpu_compute_multipole_real(neighbor->ago, inum, nall, atom->x,
                                                 atom->type, amtype, amgroup, rpole,
                                                 sublo, subhi, atom->tag,
                                                 atom->nspecial, atom->special,
                                                 atom->nspecial15, atom->special15,
                                                 eflag, vflag, eflag_atom, vflag_atom,
                                                 host_start, &ilist, &numneigh, cpu_time,
                                                 success, aewald, felec, off2, atom->q,
                                                 domain->boxlo, domain->prd, &tq_pinned);

  if (!success)
    error->one(FLERR,"Insufficient memory on accelerator");

  // reference to the tep array from GPU lib

  if (tq_single) {
    float *tq_ptr = (float *)tq_pinned;
    compute_force_from_torque<float>(tq_ptr, f, virmpole); // fmpole
  } else {
    double *tq_ptr = (double *)tq_pinned;
    compute_force_from_torque<double>(tq_ptr, f, virmpole); // fmpole
  }
}

/* ----------------------------------------------------------------------
   induce = induced dipole moments via pre-conditioned CG solver
   adapted from Tinker induce0a() routine
   NOTE: Almost the same in the CPU version, except that there is no need
      to call reverse_comm() for crstyle = FIELD;
------------------------------------------------------------------------- */

void PairAmoebaGPU::induce()
{
  bool done;
  int i,j,m,ii,itype;
  int iter,maxiter;
  double polmin;
  double eps,epsold;
  double epsd,epsp;
  double udsum,upsum;
  double a,ap,b,bp;
  double sum,sump,term;
  double reduce[4],allreduce[4];

  int debug = 1;

  // set cutoffs, taper coeffs, and PME params
  // create qfac here, free at end of polar()

  if (use_ewald) choose(POLAR_LONG);
  else choose(POLAR);

  // owned atoms

  int nlocal = atom->nlocal;

  // zero out the induced dipoles at each site

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < 3; j++) {
      uind[i][j] = 0.0;
      uinp[i][j] = 0.0;
    }
  }

  // get the electrostatic field due to permanent multipoles

  dfield0c(field,fieldp);

  // need reverse_comm if dfield0c (i.e. udirect2b) is CPU-only

  if (!gpu_udirect2b_ready) {
    crstyle = FIELD;
    comm->reverse_comm(this);
  }

  // set induced dipoles to polarizability times direct field

  for (i = 0; i < nlocal; i++) {
    itype = amtype[i];
    for (j = 0; j < 3; j++) {
      udir[i][j] = polarity[itype] * field[i][j];
      udirp[i][j] = polarity[itype] * fieldp[i][j];
      if (pcgguess) {
        uind[i][j] = udir[i][j];
        uinp[i][j] = udirp[i][j];
      }
    }
  }

  // get induced dipoles via the OPT extrapolation method
  // NOTE: any way to rewrite these loops to avoid allocating
  //       uopt,uoptp with a optorder+1 dimension, just optorder ??
  //       since no need to store optorder+1 values after these loops

  if (poltyp == OPT) {
    for (i = 0; i < nlocal; i++) {
      for (j = 0; j < 3; j++) {
        uopt[i][0][j] = udir[i][j];
        uoptp[i][0][j] = udirp[i][j];
      }
    }

    for (m = 1; m <= optorder; m++) {
      optlevel = m - 1;     // used in umutual1() for fopt,foptp

      cfstyle = INDUCE;
      comm->forward_comm(this);

      ufield0c(field,fieldp);

      if (!gpu_umutual2b_ready) {
        crstyle = FIELD;
        comm->reverse_comm(this);
      }

      for (i = 0; i < nlocal; i++) {
	      itype = amtype[i];
        for (j = 0; j < 3; j++) {
          uopt[i][m][j] = polarity[itype] * field[i][j];
          uoptp[i][m][j] = polarity[itype] * fieldp[i][j];
          uind[i][j] = uopt[i][m][j];
          uinp[i][j] = uoptp[i][m][j];
        }
      }
    }

    for (i = 0; i < nlocal; i++) {
      for (j = 0; j < 3; j++) {
        uind[i][j] = 0.0;
        uinp[i][j] = 0.0;
        usum[i][j] = 0.0;
        usump[i][j] = 0.0;
        for (m = 0; m <= optorder; m++) {
          usum[i][j] += uopt[i][m][j];
          usump[i][j] += uoptp[i][m][j];
          uind[i][j] += copt[m]*usum[i][j];
          uinp[i][j] += copt[m]*usump[i][j];
        }
      }
    }
  }

  // set tolerances for computation of mutual induced dipoles

  if (poltyp == MUTUAL) {
    done = false;
    maxiter = 100;
    iter = 0;
    polmin = 0.00000001;
    eps = 100.0;

    // estimate induced dipoles using a polynomial predictor

    if (use_pred && nualt == maxualt) {
      ulspred();

      double ***udalt = fixudalt->tstore;
      double ***upalt = fixupalt->tstore;

      for (i = 0; i < nlocal; i++) {
        for (j = 0; j < 3; j++) {
          udsum = 0.0;
          upsum = 0.0;
          for (m = 0; m < nualt; m++) {
            udsum += bpred[m]*udalt[i][m][j];
            upsum += bpredp[m]*upalt[i][m][j];
          }
          uind[i][j] = udsum;
          uinp[i][j] = upsum;
        }
      }
    }

    // estimate induced dipoles via inertial extended Lagrangian
    // not supported for now
    // requires uaux,upaux to persist with each atom
    // also requires a velocity vector(s) to persist
    // also requires updating uaux,upaux in the Verlet integration

    /*
    if (use_ielscf) {
      for (i = 0; i < nlocal; i++) {
        for (j = 0; j < 3; j++) {
          uind[i][j] = uaux[i][j];
          uinp[i][j] = upaux[i][j];
        }
      }
    }
    */

    // get the electrostatic field due to induced dipoles

    cfstyle = INDUCE;
    comm->forward_comm(this);

    ufield0c(field,fieldp);

    if (!gpu_umutual2b_ready) {
      crstyle = FIELD;
      comm->reverse_comm(this);
    }

    //error->all(FLERR,"STOP GPU");

    // set initial conjugate gradient residual and conjugate vector

    for (i = 0; i < nlocal; i++) {
      itype = amtype[i];

      poli[i] = MAX(polmin,polarity[itype]);
      for (j = 0; j < 3; j++) {
        if (pcgguess) {
          rsd[i][j] = (udir[i][j]-uind[i][j])/poli[i] + field[i][j];
          rsdp[i][j] = (udirp[i][j]-uinp[i][j])/poli[i] + fieldp[i][j];
        } else {
          rsd[i][j] = udir[i][j] / poli[i];
          rsdp[i][j] = udirp[i][j] / poli[i];
        }
        zrsd[i][j] = rsd[i][j];
        zrsdp[i][j] = rsdp[i][j];
      }
    }

    if (pcgprec) {
      cfstyle = RSD;
      comm->forward_comm(this);
      uscale0b(BUILD,rsd,rsdp,zrsd,zrsdp);
      uscale0b(APPLY,rsd,rsdp,zrsd,zrsdp);
      crstyle = ZRSD;
      comm->reverse_comm(this);
   }

    for (i = 0; i < nlocal; i++) {
      for (j = 0; j < 3; j++) {
        conj[i][j] = zrsd[i][j];
        conjp[i][j] = zrsdp[i][j];
      }
    }

    // conjugate gradient iteration of the mutual induced dipoles

    while (!done) {
      iter++;

      for (i = 0; i < nlocal; i++) {
        for (j = 0; j < 3; j++) {
          vec[i][j] = uind[i][j];
          vecp[i][j] = uinp[i][j];
          uind[i][j] = conj[i][j];
          uinp[i][j] = conjp[i][j];
        }
      }

      cfstyle = INDUCE;
      comm->forward_comm(this);

      ufield0c(field,fieldp);

      if (!gpu_umutual2b_ready) {
        crstyle = FIELD;
        comm->reverse_comm(this);
      }

      //error->all(FLERR,"STOP");

      for (i = 0; i < nlocal; i++) {
        for (j = 0; j < 3; j++) {
          uind[i][j] = vec[i][j];
          uinp[i][j] = vecp[i][j];
          vec[i][j] = conj[i][j]/poli[i] - field[i][j];
          vecp[i][j] = conjp[i][j]/poli[i] - fieldp[i][j];
        }
      }

      a = 0.0;
      ap = 0.0;
      sum = 0.0;
      sump = 0.0;

      for (i = 0; i < nlocal; i++) {
        for (j = 0; j < 3; j++) {
          a += conj[i][j]*vec[i][j];
          ap += conjp[i][j]*vecp[i][j];
          sum += rsd[i][j]*zrsd[i][j];
          sump += rsdp[i][j]*zrsdp[i][j];
        }
      }

      reduce[0] = a;
      reduce[1] = ap;
      reduce[2] = sum;
      reduce[3] = sump;
      MPI_Allreduce(reduce,allreduce,4,MPI_DOUBLE,MPI_SUM,world);
      a = allreduce[0];
      ap = allreduce[1];
      sum = allreduce[2];
      sump = allreduce[3];

      if (a != 0.0) a = sum / a;
      if (ap != 0.0) ap = sump / ap;

      for (i = 0; i < nlocal; i++) {
        for (j = 0; j < 3; j++) {
          uind[i][j] = uind[i][j] + a*conj[i][j];
          uinp[i][j] = uinp[i][j] + ap*conjp[i][j];
          rsd[i][j] = rsd[i][j] - a*vec[i][j];
          rsdp[i][j] = rsdp[i][j] - ap*vecp[i][j];
          zrsd[i][j] = rsd[i][j];
          zrsdp[i][j] = rsdp[i][j];
        }
      }

      if (pcgprec) {
        cfstyle = RSD;
        comm->forward_comm(this);
        uscale0b(APPLY,rsd,rsdp,zrsd,zrsdp);
        crstyle = ZRSD;
        comm->reverse_comm(this);
      }

      b = 0.0;
      bp = 0.0;

      for (i = 0; i < nlocal; i++) {
        for (j = 0; j < 3; j++) {
          b += rsd[i][j]*zrsd[i][j];
          bp += rsdp[i][j]*zrsdp[i][j];
        }
      }

      reduce[0] = b;
      reduce[1] = bp;
      MPI_Allreduce(reduce,allreduce,4,MPI_DOUBLE,MPI_SUM,world);
      b = allreduce[0];
      bp = allreduce[1];

      if (sum != 0.0) b /= sum;
      if (sump != 0.0) bp /= sump;

      epsd = 0.0;
      epsp = 0.0;

      for (i = 0; i < nlocal; i++) {
        for (j = 0; j < 3; j++) {
          conj[i][j] = zrsd[i][j] + b*conj[i][j];
          conjp[i][j] = zrsdp[i][j] + bp*conjp[i][j];
          epsd += rsd[i][j]*rsd[i][j];
          epsp += rsdp[i][j]*rsdp[i][j];
        }
      }

      reduce[0] = epsd;
      reduce[1] = epsp;
      MPI_Allreduce(reduce,allreduce,4,MPI_DOUBLE,MPI_SUM,world);
      epsd = allreduce[0];
      epsp = allreduce[1];

      // check the convergence of the mutual induced dipoles

      epsold = eps;
      eps = MAX(epsd,epsp);
      eps = DEBYE * sqrt(eps/atom->natoms);

      if (eps < poleps) done = true;
      if (eps > epsold) done = true;
      if (iter >= politer) done = true;

      //  apply a "peek" iteration to the mutual induced dipoles

      if (done) {
        for (i = 0; i < nlocal; i++) {
          term = pcgpeek * poli[i];
          for (j = 0; j < 3; j++) {
            uind[i][j] += term*rsd[i][j];
            uinp[i][j] += term*rsdp[i][j];
          }
        }
      }

    }

    // terminate the calculation if dipoles failed to converge
    // NOTE: could make this an error

    if (iter >= maxiter || eps > epsold)
      if (comm->me == 0)
	      error->warning(FLERR,"AMOEBA induced dipoles did not converge");
  }

  // update the lists of previous induced dipole values
  // shift previous m values up to m+1, add new values at m = 0
  // only when preconditioner is used

  if (use_pred) {
    double ***udalt = fixudalt->tstore;
    double ***upalt = fixupalt->tstore;

    nualt = MIN(nualt+1,maxualt);
    for (i = 0; i < nlocal; i++) {
      for (j = 0; j < 3; j++) {
        for (m = nualt-1; m > 0; m--) {
          udalt[i][m][j] = udalt[i][m-1][j];
          upalt[i][m][j] = upalt[i][m-1][j];
        }
        udalt[i][0][j] = uind[i][j];
        upalt[i][0][j] = uinp[i][j];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   udirect2b = Ewald real direct field via list
   udirect2b computes the real space contribution of the permanent
   atomic multipole moments to the field via a neighbor list
------------------------------------------------------------------------- */

void PairAmoebaGPU::udirect2b(double **field, double **fieldp)
{
  if (!gpu_udirect2b_ready) {
    PairAmoeba::udirect2b(field, fieldp);
    return;
  }

  int eflag=1, vflag=1;
  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

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

  // select the correct cutoff (off2) for the term

  if (use_ewald) choose(POLAR_LONG);
  else choose(POLAR);

  amoeba_gpu_compute_udirect2b(amtype, amgroup, rpole, uind, uinp,
                               aewald, off2, &fieldp_pinned);

  // rebuild dipole-dipole pair list and store pairwise dipole matrices
  // done one atom at a time in real-space double loop over atoms & neighs
  // NOTE: for the moment the tdipdip values are computed just in time in umutual2b()
  //   so no need to call ubdirect2b_cpu().
  // udirect2b_cpu();

  // accumulate the field and fieldp values from the GPU lib
  //   field and fieldp may already have some nonzero values from kspace (udirect1)

  int nlocal = atom->nlocal;
  double *field_ptr = (double *)fieldp_pinned;

  for (int i = 0; i < nlocal; i++) {
    int idx = 4*i;
    field[i][0] += field_ptr[idx];
    field[i][1] += field_ptr[idx+1];
    field[i][2] += field_ptr[idx+2];
  }

  double* fieldp_ptr = (double *)fieldp_pinned;
  fieldp_ptr += 4*inum;
  for (int i = 0; i < nlocal; i++) {
    int idx = 4*i;
    fieldp[i][0] += fieldp_ptr[idx];
    fieldp[i][1] += fieldp_ptr[idx+1];
    fieldp[i][2] += fieldp_ptr[idx+2];
  }

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
      } else {
        if (comm->me == 0) printf("i = %d: j = %d: poltyp == DIRECT\n", i, j);
      }

    } // jj

    firstneigh_dipole[i] = neighptr;
    firstneigh_dipdip[i] = tdipdip;
    numneigh_dipole[i] = n;
    ipage_dipole->vgot(n);
    dpage_dipdip->vgot(ndip);
  }
}

/* ----------------------------------------------------------------------
   ufield0c = mutual induction via Ewald sum
   ufield0c computes the mutual electrostatic field due to
   induced dipole moments via Ewald summation
------------------------------------------------------------------------- */

void PairAmoebaGPU::ufield0c(double **field, double **fieldp)
{
  int i,j;
  double term;

  // zero field,fieldp for owned and ghost atoms

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  for (i = 0; i < nall; i++) {
    for (j = 0; j < 3; j++) {
      field[i][j] = 0.0;
      fieldp[i][j] = 0.0;
    }
  }

  // get the real space portion of the mutual field first

  if (polar_rspace_flag) umutual2b(field,fieldp);

  // get the reciprocal space part of the mutual field

  if (polar_kspace_flag) umutual1(field,fieldp);

  // add the self-energy portion of the mutual field

  term = (4.0/3.0) * aewald*aewald*aewald / MY_PIS;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < 3; j++) {
      field[i][j] += term*uind[i][j];
      fieldp[i][j] += term*uinp[i][j];
    }
  }

  // accumulate the field and fieldp values from the real space portion from umutual2b() on the GPU
  //   field and fieldp may already have some nonzero values from kspace (umutual1 and self)

  amoeba_gpu_update_fieldp(&fieldp_pinned);
  
  int inum = atom->nlocal;
  double *field_ptr = (double *)fieldp_pinned;

  for (int i = 0; i < nlocal; i++) {
    int idx = 4*i;
    field[i][0] += field_ptr[idx];
    field[i][1] += field_ptr[idx+1];
    field[i][2] += field_ptr[idx+2];
  }

  double* fieldp_ptr = (double *)fieldp_pinned;
  fieldp_ptr += 4*inum;
  for (int i = 0; i < nlocal; i++) {
    int idx = 4*i;
    fieldp[i][0] += fieldp_ptr[idx];
    fieldp[i][1] += fieldp_ptr[idx+1];
    fieldp[i][2] += fieldp_ptr[idx+2];
  }
}

/* ----------------------------------------------------------------------
   umutual2b = Ewald real mutual field via list
   umutual2b computes the real space contribution of the induced
   atomic dipole moments to the field via a neighbor list
------------------------------------------------------------------------- */

void PairAmoebaGPU::umutual2b(double **field, double **fieldp)
{
  if (!gpu_umutual2b_ready) {
    PairAmoeba::umutual2b(field, fieldp);
    return;
  }

  int eflag=1, vflag=1;
  int nall = atom->nlocal + atom->nghost;
  int inum;

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

  // select the correct cutoff (off2) for the term

  if (use_ewald) choose(POLAR_LONG);
  else choose(POLAR);

  amoeba_gpu_compute_umutual2b(amtype, amgroup, rpole, uind, uinp,
                               aewald, off2, &fieldp_pinned);
}

/* ----------------------------------------------------------------------
   polar_real = real-space portion of induced dipole polarization
   adapted from Tinker epreal1d() routine
------------------------------------------------------------------------- */

void PairAmoebaGPU::polar_real()
{
  if (!gpu_polar_real_ready) {
    PairAmoeba::polar_real();
    return;
  }

  int eflag=1, vflag=1;
  double **f = atom->f;
  int nall = atom->nlocal + atom->nghost;
  int inum;

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

  // select the correct cutoff and aewald for the term

  if (use_ewald) choose(POLAR_LONG);
  else choose(POLAR);

  // set the energy unit conversion factor for polar real-space calculation

  double felec = 0.5 * electric / am_dielectric;

  amoeba_gpu_compute_polar_real(amtype, amgroup, rpole, uind, uinp,
                                eflag, vflag, eflag_atom, vflag_atom,
                                aewald, felec, off2, &tq_pinned);

  // reference to the tep array from GPU lib

  if (tq_single) {
    float *tep_ptr = (float *)tq_pinned;
    compute_force_from_torque<float>(tep_ptr, f, virpolar); // fpolar
  } else {
    double *tep_ptr = (double *)tq_pinned;
    compute_force_from_torque<double>(tep_ptr, f, virpolar); // fpolar
  }
}

/* ----------------------------------------------------------------------
   compute atom forces from torques
------------------------------------------------------------------------- */

template <class numtyp>
void PairAmoebaGPU::compute_force_from_torque(const numtyp* tq_ptr,
                                              double** force_comp,
                                              double* virial_comp)
{
  int i,ix,iy,iz;
  double xix,yix,zix;
  double xiy,yiy,ziy;
  double xiz,yiz,ziz;
  double vxx,vyy,vzz;
  double vxy,vxz,vyz;
  double fix[3],fiy[3],fiz[3],_tq[4];

  double** x = atom->x;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    _tq[0] = tq_ptr[4*i];
    _tq[1] = tq_ptr[4*i+1];
    _tq[2] = tq_ptr[4*i+2];
    torque2force(i,_tq,fix,fiy,fiz,force_comp);

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

    virial_comp[0] -= vxx;
    virial_comp[1] -= vyy;
    virial_comp[2] -= vzz;
    virial_comp[3] -= vxy;
    virial_comp[4] -= vxz;
    virial_comp[5] -= vyz;
  }
}

/* ---------------------------------------------------------------------- */

double PairAmoebaGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + amoeba_gpu_bytes();
}
