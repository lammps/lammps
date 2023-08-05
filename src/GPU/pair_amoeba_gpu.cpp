// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Nguyen (Northwestern/UChicago)
------------------------------------------------------------------------- */

#include "pair_amoeba_gpu.h"

#include "amoeba_convolution_gpu.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_store_atom.h"
#include "force.h"
#include "gpu_extra.h"
#include "info.h"
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

// same as in amoeba_induce.cpp
enum{INDUCE,RSD,SETUP_AMOEBA,SETUP_HIPPO,KMPOLE,AMGROUP};   // forward comm
enum{FIELD,ZRSD,TORQUE,UFLD};                               // reverse comm
enum{VDWL,REPULSE,QFER,DISP,MPOLE,POLAR,USOLV,DISP_LONG,MPOLE_LONG,POLAR_LONG};
enum{MUTUAL,OPT,TCG,DIRECT};
enum{GEAR,ASPC,LSQR};
enum{BUILD,APPLY};
enum{GORDON1,GORDON2};

// same as in pair_amoeba.cpp
enum{MPOLE_GRID,POLAR_GRID,POLAR_GRIDC,DISP_GRID,INDUCE_GRID,INDUCE_GRIDC};

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
                    const double polar_dscale, const double polar_uscale);
void amoeba_gpu_clear();

int** amoeba_gpu_precompute(const int ago, const int inum_full, const int nall,
                            double **host_x, int *host_type, int *host_amtype,
                            int *host_amgroup, double **host_rpole,
                            double **host_uind, double **host_uinp, double *host_pval,
                            double *sublo, double *subhi, tagint *tag,
                            int **nspecial, tagint **special,
                            int *nspecial15, tagint **special15,
                            const bool eflag_in, const bool vflag_in,
                            const bool eatom, const bool vatom, int &host_start,
                            int **ilist, int **jnum, const double cpu_time,
                            bool &success, double *host_q, double *boxlo, double *prd);

void amoeba_gpu_compute_multipole_real(const int ago, const int inum, const int nall,
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

void amoeba_gpu_precompute_kspace(const int inum_full, const int bsorder,
              double ***host_thetai1, double ***host_thetai2,
              double ***host_thetai3, int** igrid,
              const int nzlo_out, const int nzhi_out,
              const int nylo_out, const int nyhi_out,
              const int nxlo_out, const int nxhi_out);

void amoeba_gpu_fphi_uind(double ****host_grid_brick, void **host_fdip_phi1,
                          void **host_fdip_phi2, void **host_fdip_sum_phi);

void amoeba_gpu_fphi_mpole(double ***host_grid_brick, void **host_fdip_sum_phi,
                           const double felec);

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
  gpu_umutual1_ready = true;
  gpu_fphi_uind_ready = true;
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

/* ---------------------------------------------------------------------- */

void PairAmoebaGPU::compute(int eflag, int vflag)
{
  if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0)
    neighbor->build_topology();
  PairAmoeba::compute(eflag, vflag);
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

  int mnf = 5e-2 * neighbor->oneatom;
  int success = amoeba_gpu_init(atom->ntypes+1, max_amtype, max_amclass,
                                pdamp, thole, dirdamp, amtype2class, special_hal,
                                special_repel, special_disp, special_mpole,
                                special_polar_wscale, special_polar_piscale,
                                special_polar_pscale, csix, adisp, atom->nlocal,
                                atom->nlocal+atom->nghost, mnf, maxspecial,
                                maxspecial15, cell_size, gpu_mode, screen,
                                polar_dscale, polar_uscale);
  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE)
    error->all(FLERR,"Pair style amoeba/gpu does not support neigh no for now");

  acc_float = Info::has_accelerator_feature("GPU", "precision", "single");

  // replace with the gpu counterpart

  if (gpu_umutual1_ready) {
    if (use_ewald && ic_kspace) {
      delete ic_kspace;
      ic_kspace =
        new AmoebaConvolutionGPU(lmp,this,nefft1,nefft2,nefft3,bsporder,INDUCE_GRIDC);
    }
  }
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
  int *ilist, *numneigh;

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

  amoeba_gpu_precompute(neighbor->ago, inum, nall, atom->x,
                        atom->type, amtype, amgroup, rpole,
                        nullptr, nullptr, nullptr,
                        sublo, subhi, atom->tag,
                        atom->nspecial, atom->special,
                        atom->nspecial15, atom->special15,
                        eflag, vflag, eflag_atom, vflag_atom,
                        host_start, &ilist, &numneigh, cpu_time,
                        success, atom->q, domain->boxlo, domain->prd);
  if (!success)
    error->one(FLERR,"Insufficient memory on accelerator");

  // select the correct cutoff for the term

  if (use_ewald) choose(MPOLE_LONG);
  else choose(MPOLE);

  // set the energy unit conversion factor for multipolar real-space calculation

  double felec = electric / am_dielectric;

  amoeba_gpu_compute_multipole_real(neighbor->ago, inum, nall, atom->x,
                                    atom->type, amtype, amgroup, rpole,
                                    sublo, subhi, atom->tag,
                                    atom->nspecial, atom->special,
                                    atom->nspecial15, atom->special15,
                                    eflag, vflag, eflag_atom, vflag_atom,
                                    host_start, &ilist, &numneigh, cpu_time,
                                    success, aewald, felec, off2, atom->q,
                                    domain->boxlo, domain->prd, &tq_pinned);



  // reference to the tep array from GPU lib

  if (acc_float) {
    auto *tq_ptr = (float *)tq_pinned;
    compute_force_from_torque<float>(tq_ptr, f, virmpole); // fmpole
  } else {
    auto *tq_ptr = (double *)tq_pinned;
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
  int i,j,m,itype;
  int iter,maxiter;
  double polmin;
  double eps,epsold;
  double epsd,epsp;
  double udsum,upsum;
  double a,ap,b,bp;
  double sum,sump,term;
  double reduce[4],allreduce[4];

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

  // allocate memory and make early host-device transfers
  // must be done before the first ufield0c
  // NOTE: this is for ic_kspace, and thetai[1-3]

  if (ic_kspace)
    amoeba_gpu_precompute_kspace(atom->nlocal, bsorder, thetai1, thetai2,
                                 thetai3, igrid,
                                 ic_kspace->nzlo_out, ic_kspace->nzhi_out,
                                 ic_kspace->nylo_out, ic_kspace->nyhi_out,
                                 ic_kspace->nxlo_out, ic_kspace->nxhi_out);

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
  if (acc_float) {
    auto field_ptr = (float *)fieldp_pinned;

    for (int i = 0; i < nlocal; i++) {
      int idx = 3*i;
      field[i][0] += field_ptr[idx];
      field[i][1] += field_ptr[idx+1];
      field[i][2] += field_ptr[idx+2];
    }

    field_ptr += 3*inum;
    for (int i = 0; i < nlocal; i++) {
      int idx = 3*i;
      fieldp[i][0] += field_ptr[idx];
      fieldp[i][1] += field_ptr[idx+1];
      fieldp[i][2] += field_ptr[idx+2];
    }
  } else {
    auto field_ptr = (double *)fieldp_pinned;

    for (int i = 0; i < nlocal; i++) {
      int idx = 3*i;
      field[i][0] += field_ptr[idx];
      field[i][1] += field_ptr[idx+1];
      field[i][2] += field_ptr[idx+2];
    }

    field_ptr += 3*inum;
    for (int i = 0; i < nlocal; i++) {
      int idx = 3*i;
      fieldp[i][0] += field_ptr[idx];
      fieldp[i][1] += field_ptr[idx+1];
      fieldp[i][2] += field_ptr[idx+2];
    }
  }


}

/* ----------------------------------------------------------------------
   udirect2b = Ewald real direct field via list
   udirect2b computes the real space contribution of the permanent
     atomic multipole moments to the field via a neighbor list
------------------------------------------------------------------------- */

void PairAmoebaGPU::udirect2b_cpu()
{
  int i,j,m,n,ii,jj,jextra,ndip,itype,jtype,igroup,jgroup;
  double xr,yr,zr,r,r2;
  double rr1,rr2,rr3,rr5;
  double bfac,exp2a;
  double ralpha,aefac;
  double aesq2,aesq2n;
  double pdi,pti;
  double pgamma;
  double damp,expdamp;
  double scale3,scale5;
  double scalek;
  double bn[4],bcn[3];
  double factor_uscale;

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

      if (igroup == jgroup) factor_uscale = polar_uscale;
      else factor_uscale = 1.0;

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
  double term;

  // zero field,fieldp for owned and ghost atoms

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  memset(&field[0][0], 0, 3*nall *sizeof(double));
  memset(&fieldp[0][0], 0, 3*nall *sizeof(double));

  // get the real space portion of the mutual field first

  double time0, time1, time2;

  MPI_Barrier(world);
  time0 = platform::walltime();

  if (polar_rspace_flag) umutual2b(field,fieldp);
  time1 = platform::walltime();

  // get the reciprocal space part of the mutual field

  if (polar_kspace_flag) umutual1(field,fieldp);
  time2 = platform::walltime();

  // add the self-energy portion of the mutual field

  term = (4.0/3.0) * aewald*aewald*aewald / MY_PIS;
  for (int i = 0; i < nlocal; i++) {
    field[i][0] += term*uind[i][0];
    field[i][1] += term*uind[i][1];
    field[i][2] += term*uind[i][2];
  }

  for (int i = 0; i < nlocal; i++) {
    fieldp[i][0] += term*uinp[i][0];
    fieldp[i][1] += term*uinp[i][1];
    fieldp[i][2] += term*uinp[i][2];
  }

  // accumulate the field and fieldp values from the real-space portion from umutual2b() on the GPU
  //   field and fieldp may already have some nonzero values from kspace (umutual1 and self)

  amoeba_gpu_update_fieldp(&fieldp_pinned);

  int inum = atom->nlocal;
  if (acc_float) {
    auto field_ptr = (float *)fieldp_pinned;

    for (int i = 0; i < nlocal; i++) {
      int idx = 3*i;
      field[i][0] += field_ptr[idx];
      field[i][1] += field_ptr[idx+1];
      field[i][2] += field_ptr[idx+2];
    }

    field_ptr += 3*inum;
    for (int i = 0; i < nlocal; i++) {
      int idx = 3*i;
      fieldp[i][0] += field_ptr[idx];
      fieldp[i][1] += field_ptr[idx+1];
      fieldp[i][2] += field_ptr[idx+2];
    }
  } else {
    auto field_ptr = (double *)fieldp_pinned;

    for (int i = 0; i < nlocal; i++) {
      int idx = 3*i;
      field[i][0] += field_ptr[idx];
      field[i][1] += field_ptr[idx+1];
      field[i][2] += field_ptr[idx+2];
    }

    field_ptr += 3*inum;
    for (int i = 0; i < nlocal; i++) {
      int idx = 3*i;
      fieldp[i][0] += field_ptr[idx];
      fieldp[i][1] += field_ptr[idx+1];
      fieldp[i][2] += field_ptr[idx+2];
    }
  }


  // accumulate timing information

  time_mutual_rspace += time1 - time0;
  time_mutual_kspace += time2 - time1;
}

/* ----------------------------------------------------------------------
   umutual1 = Ewald recip mutual induced field
   umutual1 computes the reciprocal space contribution of the
   induced atomic dipole moments to the field
------------------------------------------------------------------------- */

void PairAmoebaGPU::umutual1(double **field, double **fieldp)
{
  int m,n;
  int nxlo,nxhi,nylo,nyhi,nzlo,nzhi;
  double term;
  double a[3][3];  // indices not flipped vs Fortran

  // return if the Ewald coefficient is zero

  if (aewald < 1.0e-6) return;

  // convert Cartesian dipoles to fractional coordinates

  for (int j = 0; j < 3; j++) {
    a[0][j] = nfft1 * recip[0][j];
    a[1][j] = nfft2 * recip[1][j];
    a[2][j] = nfft3 * recip[2][j];
  }

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    fuind[i][0] = a[0][0]*uind[i][0] + a[0][1]*uind[i][1] + a[0][2]*uind[i][2];
    fuind[i][1] = a[1][0]*uind[i][0] + a[1][1]*uind[i][1] + a[1][2]*uind[i][2];
    fuind[i][2] = a[2][0]*uind[i][0] + a[2][1]*uind[i][1] + a[2][2]*uind[i][2];
  }

  for (int i = 0; i < nlocal; i++) {
    fuinp[i][0] = a[0][0]*uinp[i][0] + a[0][1]*uinp[i][1] + a[0][2]*uinp[i][2];
    fuinp[i][1] = a[1][0]*uinp[i][0] + a[1][1]*uinp[i][1] + a[1][2]*uinp[i][2];
    fuinp[i][2] = a[2][0]*uinp[i][0] + a[2][1]*uinp[i][1] + a[2][2]*uinp[i][2];
  }

  // gridpre = my portion of 4d grid in brick decomp w/ ghost values

  FFT_SCALAR ****gridpre = (FFT_SCALAR ****) ic_kspace->zero();

  // map 2 values to grid

  double time0, time1;
  MPI_Barrier(world);
  time0 = platform::walltime();

  grid_uind(fuind,fuinp,gridpre);

  time1 = platform::walltime();
  time_grid_uind += (time1 - time0);

  // pre-convolution operations including forward FFT
  // gridfft = my portion of complex 3d grid in FFT decomposition

  FFT_SCALAR *gridfft = ic_kspace->pre_convolution();

  // ---------------------
  // convolution operation
  // ---------------------

  nxlo = ic_kspace->nxlo_fft;
  nxhi = ic_kspace->nxhi_fft;
  nylo = ic_kspace->nylo_fft;
  nyhi = ic_kspace->nyhi_fft;
  nzlo = ic_kspace->nzlo_fft;
  nzhi = ic_kspace->nzhi_fft;

  // use qfac values stored in udirect1()

  m = n = 0;
  for (int k = nzlo; k <= nzhi; k++) {
    for (int j = nylo; j <= nyhi; j++) {
      for (int i = nxlo; i <= nxhi; i++) {
        term = qfac[m++];
        gridfft[n] *= term;
        gridfft[n+1] *= term;
        n += 2;
      }
    }
  }

  // post-convolution operations including backward FFT
  // gridppost = my portion of 4d grid in brick decomp w/ ghost values

  FFT_SCALAR ****gridpost = (FFT_SCALAR ****) ic_kspace->post_convolution();

  // get potential

  MPI_Barrier(world);
  time0 = platform::walltime();

  fphi_uind(gridpost,fdip_phi1,fdip_phi2,fdip_sum_phi);

  time1 = platform::walltime();
  time_fphi_uind += (time1 - time0);

  // store fractional reciprocal potentials for OPT method

  if (poltyp == OPT) {
    for (int i = 0; i < nlocal; i++) {
      for (int j = 0; j < 10; j++) {
        fopt[i][optlevel][j] = fdip_phi1[i][j];
        foptp[i][optlevel][j] = fdip_phi2[i][j];
      }
    }
  }

  for (int i = 0; i < nlocal; i++) {
    double dfx = a[0][0]*fdip_phi1[i][1] +
      a[0][1]*fdip_phi1[i][2] + a[0][2]*fdip_phi1[i][3];
    double dfy = a[1][0]*fdip_phi1[i][1] +
      a[1][1]*fdip_phi1[i][2] + a[1][2]*fdip_phi1[i][3];
    double dfz = a[2][0]*fdip_phi1[i][1] +
      a[2][1]*fdip_phi1[i][2] + a[2][2]*fdip_phi1[i][3];
    field[i][0] -= dfx;
    field[i][1] -= dfy;
    field[i][2] -= dfz;
  }

  for (int i = 0; i < nlocal; i++) {
    double dfx = a[0][0]*fdip_phi2[i][1] +
      a[0][1]*fdip_phi2[i][2] + a[0][2]*fdip_phi2[i][3];
    double dfy = a[1][0]*fdip_phi2[i][1] +
      a[1][1]*fdip_phi2[i][2] + a[1][2]*fdip_phi2[i][3];
    double dfz = a[2][0]*fdip_phi2[i][1] +
      a[2][1]*fdip_phi2[i][2] + a[2][2]*fdip_phi2[i][3];
    fieldp[i][0] -= dfx;
    fieldp[i][1] -= dfy;
    fieldp[i][2] -= dfz;
  }

}

/* ----------------------------------------------------------------------
   fphi_uind = induced potential from grid
   fphi_uind extracts the induced dipole potential from the particle mesh Ewald grid
------------------------------------------------------------------------- */

void PairAmoebaGPU::fphi_uind(FFT_SCALAR ****grid, double **fdip_phi1,
                              double **fdip_phi2, double **fdip_sum_phi)
{
  if (!gpu_fphi_uind_ready) {
    PairAmoeba::fphi_uind(grid, fdip_phi1, fdip_phi2, fdip_sum_phi);
    return;
  }

  void* fdip_phi1_pinned = nullptr;
  void* fdip_phi2_pinned = nullptr;
  void* fdip_sum_phi_pinned = nullptr;
  amoeba_gpu_fphi_uind(grid, &fdip_phi1_pinned, &fdip_phi2_pinned,
                       &fdip_sum_phi_pinned);

  int nlocal = atom->nlocal;
  if (acc_float) {
    auto _fdip_phi1_ptr = (float *)fdip_phi1_pinned;
    for (int i = 0; i < nlocal; i++) {
      int n = i;
      for (int m = 0; m < 10; m++) {
        fdip_phi1[i][m] = _fdip_phi1_ptr[n];
        n += nlocal;
      }
    }

    auto _fdip_phi2_ptr = (float *)fdip_phi2_pinned;
    for (int i = 0; i < nlocal; i++) {
      int n = i;
      for (int m = 0; m < 10; m++) {
        fdip_phi2[i][m] = _fdip_phi2_ptr[n];
        n += nlocal;
      }
    }

    auto _fdip_sum_phi_ptr = (float *)fdip_sum_phi_pinned;
    for (int i = 0; i < nlocal; i++) {
      int n = i;
      for (int m = 0; m < 20; m++) {
        fdip_sum_phi[i][m] = _fdip_sum_phi_ptr[n];
        n += nlocal;
      }
    }

  } else {
    auto _fdip_phi1_ptr = (double *)fdip_phi1_pinned;
    for (int i = 0; i < nlocal; i++) {
      int n = i;
      for (int m = 0; m < 10; m++) {
        fdip_phi1[i][m] = _fdip_phi1_ptr[n];
        n += nlocal;
      }
    }

    auto _fdip_phi2_ptr = (double *)fdip_phi2_pinned;
    for (int i = 0; i < nlocal; i++) {
      int n = i;
      for (int m = 0; m < 10; m++) {
        fdip_phi2[i][m] = _fdip_phi2_ptr[n];
        n += nlocal;
      }
    }

    auto _fdip_sum_phi_ptr = (double *)fdip_sum_phi_pinned;
    for (int i = 0; i < nlocal; i++) {
      int n = i;
      for (int m = 0; m < 20; m++) {
        fdip_sum_phi[i][m] = _fdip_sum_phi_ptr[n];
        n += nlocal;
      }
    }
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

  // select the correct cutoff and aewald for the term

  if (use_ewald) choose(POLAR_LONG);
  else choose(POLAR);

  // set the energy unit conversion factor for polar real-space calculation

  double felec = 0.5 * electric / am_dielectric;

  amoeba_gpu_compute_polar_real(amtype, amgroup, rpole, uind, uinp,
                                eflag, vflag, eflag_atom, vflag_atom,
                                aewald, felec, off2, &tq_pinned);

  // reference to the tep array from GPU lib

  if (acc_float) {
    auto *tep_ptr = (float *)tq_pinned;
    compute_force_from_torque<float>(tep_ptr, f, virpolar); // fpolar
  } else {
    auto *tep_ptr = (double *)tq_pinned;
    compute_force_from_torque<double>(tep_ptr, f, virpolar); // fpolar
  }
}

/* ----------------------------------------------------------------------
   polar_kspace = KSpace portion of induced dipole polarization
   adapted from Tinker eprecip1() routine
   same as PairAmoeba, except that fphi_uind() is reimplemented here
 ------------------------------------------------------------------------- */

void PairAmoebaGPU::polar_kspace()
{
  int i,j,k,m,n;
  int nhalf1,nhalf2,nhalf3;
  int nxlo,nxhi,nylo,nyhi,nzlo,nzhi;
  int j1,j2,j3;
  int ix,iy,iz;
  double eterm,felec;
  double r1,r2,r3;
  double h1,h2,h3;
  double f1,f2,f3;
  double xix,yix,zix;
  double xiy,yiy,ziy;
  double xiz,yiz,ziz;
  double vxx,vyy,vzz;
  double vxy,vxz,vyz;
  double volterm,denom;
  double hsq,expterm;
  double term,pterm;
  double vterm,struc2;
  double tep[3];
  double fix[3],fiy[3],fiz[3];
  double cphid[4],cphip[4];
  double a[3][3];    // indices not flipped vs Fortran

  bool gpu_fphi_mpole_ready = true;

  // indices into the electrostatic field array
  // decremented by 1 versus Fortran

  int deriv1[10] = {1, 4, 7, 8, 10, 15, 17, 13, 14, 19};
  int deriv2[10] = {2, 7, 5, 9, 13, 11, 18, 15, 19, 16};
  int deriv3[10] = {3, 8, 9, 6, 14, 16, 12, 19, 17, 18};

  // return if the Ewald coefficient is zero

  if (aewald < 1.0e-6) return;

  // owned atoms

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;

  double volbox = domain->prd[0] * domain->prd[1] * domain->prd[2];
  pterm = pow((MY_PI/aewald),2.0);
  volterm = MY_PI * volbox;

  // initialize variables required for the scalar summation

  felec = electric / am_dielectric;

  // remove scalar sum virial from prior multipole FFT
  // can only do this if multipoles were computed with same aeewald = apewald
  // else need to re-compute it via new long-range solve

  nfft1 = p_kspace->nx;
  nfft2 = p_kspace->ny;
  nfft3 = p_kspace->nz;
  bsorder = p_kspace->order;

  nhalf1 = (nfft1+1) / 2;
  nhalf2 = (nfft2+1) / 2;
  nhalf3 = (nfft3+1) / 2;

  nxlo = p_kspace->nxlo_fft;
  nxhi = p_kspace->nxhi_fft;
  nylo = p_kspace->nylo_fft;
  nyhi = p_kspace->nyhi_fft;
  nzlo = p_kspace->nzlo_fft;
  nzhi = p_kspace->nzhi_fft;

  // use previous results or compute new qfac and convolution

  if (aewald == aeewald) {
    vxx = -vmsave[0];
    vyy = -vmsave[1];
    vzz = -vmsave[2];
    vxy = -vmsave[3];
    vxz = -vmsave[4];
    vyz = -vmsave[5];

  } else {

    // setup stencil size and B-spline coefficients

    moduli();
    bspline_fill();

    // allocate memory and make early host-device transfers

    // NOTE: this is for p_kspace, and igrid and thetai[1-3] are filled by bpsline_fill
    if (gpu_fphi_mpole_ready) {
       amoeba_gpu_precompute_kspace(atom->nlocal, bsorder,
                                    thetai1, thetai2, thetai3, igrid,
                                    p_kspace->nzlo_out, p_kspace->nzhi_out,
                                    p_kspace->nylo_out, p_kspace->nyhi_out,
                                    p_kspace->nxlo_out, p_kspace->nxhi_out);
    }


    // convert Cartesian multipoles to fractional coordinates

    cmp_to_fmp(cmp,fmp);

    // gridpre = my portion of 3d grid in brick decomp w/ ghost values

    FFT_SCALAR ***gridpre = (FFT_SCALAR ***) p_kspace->zero();

    // map atoms to grid

    grid_mpole(fmp,gridpre);

    // pre-convolution operations including forward FFT
    // gridfft = my portion of complex 3d grid in FFT decomp as 1d vector

    FFT_SCALAR *gridfft = p_kspace->pre_convolution();

    // ---------------------
    // convolution operation
    // ---------------------

    // zero virial accumulation variables

    vxx = vyy = vzz = vxy = vxz = vyz = 0.0;

    // perform convolution on K-space points I own

    m = n = 0;
    for (k = nzlo; k <= nzhi; k++) {
      for (j = nylo; j <= nyhi; j++) {
        for (i = nxlo; i <= nxhi; i++) {
          r1 = (i >= nhalf1) ? i-nfft1 : i;
          r2 = (j >= nhalf2) ? j-nfft2 : j;
          r3 = (k >= nhalf3) ? k-nfft3 : k;
          h1 = recip[0][0]*r1 + recip[0][1]*r2 + recip[0][2]*r3;  // matvec
          h2 = recip[1][0]*r1 + recip[1][1]*r2 + recip[1][2]*r3;
          h3 = recip[2][0]*r1 + recip[2][1]*r2 + recip[2][2]*r3;
          hsq = h1*h1 + h2*h2 + h3*h3;
          term = -pterm * hsq;
          expterm = 0.0;
          if (term > -50.0 && hsq != 0.0) {
            denom = volterm*hsq*bsmod1[i]*bsmod2[j]*bsmod3[k];
            if (hsq) expterm = exp(term) / denom;
            struc2 = gridfft[n]*gridfft[n] + gridfft[n+1]*gridfft[n+1];
            eterm = 0.5 * felec * expterm * struc2;
            vterm = (2.0/hsq) * (1.0-term) * eterm;
            vxx -= h1*h1*vterm - eterm;
            vyy -= h2*h2*vterm - eterm;
            vzz -= h3*h3*vterm - eterm;
            vxy -= h1*h2*vterm;
            vxz -= h1*h3*vterm;
            vyz -= h2*h3*vterm;
          }

          expterm = qfac[m++];
          gridfft[n] *= expterm;
          gridfft[n+1] *= expterm;
          n += 2;
        }
      }
    }

    // post-convolution operations including backward FFT
    // gridppost = my portion of 3d grid in brick decomp w/ ghost values

    FFT_SCALAR ***gridpost = (FFT_SCALAR ***) p_kspace->post_convolution();

    // get potential

    if (!gpu_fphi_mpole_ready) {
      fphi_mpole(gridpost,fphi);

      for (i = 0; i < nlocal; i++) {
        for (k = 0; k < 20; k++)
          fphi[i][k] *= felec;
      }

    } else {
      void* fphi_pinned = nullptr;
      amoeba_gpu_fphi_mpole(gridpost, &fphi_pinned, felec);
      if (acc_float) {
        auto _fphi_ptr = (float *)fphi_pinned;
        for (int i = 0; i < nlocal; i++) {
          int idx = i;
          for (int m = 0; m < 20; m++) {
            fphi[i][m] = _fphi_ptr[idx];
            idx += nlocal;
          }
        }
      } else {
        auto _fphi_ptr = (double *)fphi_pinned;
        for (int i = 0; i < nlocal; i++) {
          int idx = i;
          for (int m = 0; m < 20; m++) {
            fphi[i][m] = _fphi_ptr[idx];
            idx += nlocal;
          }
        }
      }
    }

    // convert field from fractional to Cartesian

    fphi_to_cphi(fphi,cphi);
  }

  // convert Cartesian induced dipoles to fractional coordinates

  for (i = 0; i < 3; i++) {
    a[0][i] = nfft1 * recip[0][i];
    a[1][i] = nfft2 * recip[1][i];
    a[2][i] = nfft3 * recip[2][i];
  }

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < 3; j++) {
      fuind[i][j] = a[j][0]*uind[i][0] + a[j][1]*uind[i][1] + a[j][2]*uind[i][2];
      fuinp[i][j] = a[j][0]*uinp[i][0] + a[j][1]*uinp[i][1] + a[j][2]*uinp[i][2];
    }
  }

  // gridpre2 = my portion of 4d grid in brick decomp w/ ghost values

  FFT_SCALAR ****gridpre2 = (FFT_SCALAR ****) pc_kspace->zero();

  // map 2 values to grid

  grid_uind(fuind,fuinp,gridpre2);

  // pre-convolution operations including forward FFT
  // gridfft = my portion of complex 3d grid in FFT decomposition

  FFT_SCALAR *gridfft = pc_kspace->pre_convolution();

  // ---------------------
  // convolution operation
  // ---------------------

  // use qfac values from above or from induce()

  m = n = 0;
  for (k = nzlo; k <= nzhi; k++) {
    for (j = nylo; j <= nyhi; j++) {
      for (i = nxlo; i <= nxhi; i++) {
        term = qfac[m++];
        gridfft[n] *= term;
        gridfft[n+1] *= term;
        n += 2;
      }
    }
  }

  // post-convolution operations including backward FFT
  // gridppost = my portion of 4d grid in brick decomp w/ ghost values

  FFT_SCALAR ****gridpost = (FFT_SCALAR ****) pc_kspace->post_convolution();

  // get potential

  fphi_uind(gridpost,fphid,fphip,fphidp);

  // TODO: port the remaining loops to the GPU

  for (i = 0; i < nlocal; i++) {
    for (j = 1; j < 10; j++) {
      fphid[i][j] = felec * fphid[i][j];
      fphip[i][j] = felec * fphip[i][j];
    }
    for (j = 0; j < 20; j++)
      fphidp[i][j] = felec * fphidp[i][j];
  }

  // increment the dipole polarization gradient contributions

  for (i = 0; i < nlocal; i++) {
    f1 = 0.0;
    f2 = 0.0;
    f3 = 0.0;
    for (k = 0; k < 3; k++) {
      j1 = deriv1[k+1];
      j2 = deriv2[k+1];
      j3 = deriv3[k+1];
      f1 += (fuind[i][k]+fuinp[i][k])*fphi[i][j1];
      f2 += (fuind[i][k]+fuinp[i][k])*fphi[i][j2];
      f3 += (fuind[i][k]+fuinp[i][k])*fphi[i][j3];
      if (poltyp == MUTUAL) {
        f1 += fuind[i][k]*fphip[i][j1] + fuinp[i][k]*fphid[i][j1];
        f2 += fuind[i][k]*fphip[i][j2] + fuinp[i][k]*fphid[i][j2];
        f3 += fuind[i][k]*fphip[i][j3] + fuinp[i][k]*fphid[i][j3];
      }
    }
    for (k = 0; k < 10; k++) {
      f1 += fmp[i][k]*fphidp[i][deriv1[k]];
      f2 += fmp[i][k]*fphidp[i][deriv2[k]];
      f3 += fmp[i][k]*fphidp[i][deriv3[k]];
    }
    f1 *= 0.5 * nfft1;
    f2 *= 0.5 * nfft2;
    f3 *= 0.5 * nfft3;
    h1 = recip[0][0]*f1 + recip[0][1]*f2 + recip[0][2]*f3;
    h2 = recip[1][0]*f1 + recip[1][1]*f2 + recip[1][2]*f3;
    h3 = recip[2][0]*f1 + recip[2][1]*f2 + recip[2][2]*f3;
    f[i][0] -= h1;
    f[i][1] -= h2;
    f[i][2] -= h3;
  }

  // set the potential to be the induced dipole average

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < 10; j++)
      fphidp[i][j] *= 0.5;
  }

  fphi_to_cphi(fphidp,cphidp);

  // get the fractional to Cartesian transformation matrix

  //frac_to_cart();

  // increment the dipole polarization virial contributions

  for (i = 0; i < nlocal; i++) {
    for (j = 1; j < 4; j++) {
      cphid[j] = 0.0;
      cphip[j] = 0.0;
      for (k = 1; k < 4; k++) {
        cphid[j] += ftc[j][k]*fphid[i][k];
        cphip[j] += ftc[j][k]*fphip[i][k];
      }
    }

    vxx -= cmp[i][1]*cphidp[i][1] +
      0.5*((uind[i][0]+uinp[i][0])*cphi[i][1]);
    vyy -= cmp[i][2]*cphidp[i][2] +
      0.5*((uind[i][1]+uinp[i][1])*cphi[i][2]);
    vzz -= cmp[i][3]*cphidp[i][3] +
      0.5*((uind[i][2]+uinp[i][2])*cphi[i][3]);
    vxy -= 0.5*(cphidp[i][1]*cmp[i][2]+cphidp[i][2]*cmp[i][1]) +
      0.25*((uind[i][1]+uinp[i][1])*cphi[i][1] +
            (uind[i][0]+uinp[i][0])*cphi[i][2]);
    vyz -= 0.5*(cphidp[i][2]*cmp[i][3]+cphidp[i][3]*cmp[i][2]) +
      0.25*((uind[i][2]+uinp[i][2])*cphi[i][2] +
            (uind[i][1]+uinp[i][1])*cphi[i][3]);
    vxz -= 0.5*(cphidp[i][1]*cmp[i][3]+cphidp[i][3]*cmp[i][1]) +
      0.25*((uind[i][2]+uinp[i][2])*cphi[i][1] +
            (uind[i][0]+uinp[i][0])*cphi[i][3]);

    vxx -= 2.0*cmp[i][4]*cphidp[i][4] + cmp[i][7]*cphidp[i][7] +
      cmp[i][8]*cphidp[i][8];
    vyy -= 2.0*cmp[i][5]*cphidp[i][5] + cmp[i][7]*cphidp[i][7] +
      cmp[i][9]*cphidp[i][9];
    vzz -= 2.0*cmp[i][6]*cphidp[i][6] + cmp[i][8]*cphidp[i][8] +
      cmp[i][9]*cphidp[i][9];
    vxy -= (cmp[i][4]+cmp[i][5])*cphidp[i][7] +
      0.5*(cmp[i][7]*(cphidp[i][5]+cphidp[i][4]) +
           cmp[i][8]*cphidp[i][9]+cmp[i][9]*cphidp[i][8]);
    vyz -= (cmp[i][5]+cmp[i][6])*cphidp[i][9] +
      0.5*(cmp[i][9]*(cphidp[i][5]+cphidp[i][6]) +
           cmp[i][7]*cphidp[i][8]+cmp[i][8]*cphidp[i][7]);
    vxz -= (cmp[i][4]+cmp[i][6])*cphidp[i][8] +
      0.5*(cmp[i][8]*(cphidp[i][4]+cphidp[i][6]) +
           cmp[i][7]*cphidp[i][9]+cmp[i][9]*cphidp[i][7]);

    if (poltyp == MUTUAL) {
      vxx -= 0.5 * (cphid[1]*uinp[i][0]+cphip[1]*uind[i][0]);
      vyy -= 0.5 * (cphid[2]*uinp[i][1]+cphip[2]*uind[i][1]);
      vzz -= 0.5 * (cphid[3]*uinp[i][2]+cphip[3]*uind[i][2]);
      vxy -= 0.25 * (cphid[1]*uinp[i][1]+cphip[1]*uind[i][1] +
                     cphid[2]*uinp[i][0]+cphip[2]*uind[i][0]);
      vyz -= 0.25 * (cphid[2]*uinp[i][2]+cphip[2]*uind[i][2] +
                     cphid[3]*uinp[i][1]+cphip[3]*uind[i][1]);
      vxz -= 0.25 * (cphid[1]*uinp[i][2]+cphip[1]*uind[i][2] +
                     cphid[3]*uinp[i][0]+cphip[3]*uind[i][0]);
    }
  }


  // resolve site torques then increment forces and virial

  for (i = 0; i < nlocal; i++) {
    tep[0] = cmp[i][3]*cphidp[i][2] - cmp[i][2]*cphidp[i][3] +
      2.0*(cmp[i][6]-cmp[i][5])*cphidp[i][9] + cmp[i][8]*cphidp[i][7] +
      cmp[i][9]*cphidp[i][5]- cmp[i][7]*cphidp[i][8] - cmp[i][9]*cphidp[i][6];
    tep[1] = cmp[i][1]*cphidp[i][3] - cmp[i][3]*cphidp[i][1] +
      2.0*(cmp[i][4]-cmp[i][6])*cphidp[i][8] + cmp[i][7]*cphidp[i][9] +
      cmp[i][8]*cphidp[i][6] - cmp[i][8]*cphidp[i][4] - cmp[i][9]*cphidp[i][7];
    tep[2] = cmp[i][2]*cphidp[i][1] - cmp[i][1]*cphidp[i][2] +
      2.0*(cmp[i][5]-cmp[i][4])*cphidp[i][7] + cmp[i][7]*cphidp[i][4] +
      cmp[i][9]*cphidp[i][8] - cmp[i][7]*cphidp[i][5] - cmp[i][8]*cphidp[i][9];

    torque2force(i,tep,fix,fiy,fiz,f);

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

    vxx += xix*fix[0] + xiy*fiy[0] + xiz*fiz[0];
    vyy += yix*fix[1] + yiy*fiy[1] + yiz*fiz[1];
    vzz += zix*fix[2] + ziy*fiy[2] + ziz*fiz[2];
    vxy += 0.5*(yix*fix[0] + yiy*fiy[0] + yiz*fiz[0] +
                xix*fix[1] + xiy*fiy[1] + xiz*fiz[1]);
    vyz += 0.5*(zix*fix[1] + ziy*fiy[1] + ziz*fiz[1] +
                yix*fix[2] + yiy*fiy[2] + yiz*fiz[2]);
    vxz += 0.5*(zix*fix[0] + ziy*fiy[0] + ziz*fiz[0] +
                xix*fix[2] + xiy*fiy[2] + xiz*fiz[2]);
  }

  // account for dipole response terms in the OPT method

  if (poltyp == OPT) {
    for (i = 0; i < nlocal; i++) {
      for (k = 0; k < optorder; k++) {
        for (j = 1; j < 10; j++) {
          fphid[i][j] = felec * fopt[i][k][j];
          fphip[i][j] = felec * foptp[i][k][j];
        }

        for (m = 0; m < optorder-k; m++) {
          for (j = 0; j < 3; j++) {
            fuind[i][j] = a[0][j]*uopt[i][m][0] + a[1][j]*uopt[i][m][1] +
              a[2][j]*uopt[i][m][2];
            fuinp[i][j] = a[0][j]*uoptp[i][m][0] + a[1][j]*uoptp[i][m][1] +
              a[2][j]*uoptp[i][m][2];
          }

          f1 = 0.0;
          f2 = 0.0;
          f3 = 0.0;

          for (j = 0; j < 3; j++) {
            j1 = deriv1[j+1];
            j2 = deriv2[j+1];
            j3 = deriv3[j+1];
            f1 += fuind[i][j]*fphip[i][j1] + fuinp[i][j]*fphid[i][j1];
            f2 += fuind[i][j]*fphip[i][j2] + fuinp[i][j]*fphid[i][j2];
            f3 += fuind[i][j]*fphip[i][j3] + fuinp[i][j]*fphid[i][j3];
          }

          f1 *= 0.5 * nfft1;
          f2 *= 0.5 * nfft2;
          f3 *= 0.5 * nfft3;
          h1 = recip[0][0]*f1 + recip[0][1]*f2 + recip[0][2]*f3;
          h2 = recip[1][0]*f1 + recip[1][1]*f2 + recip[1][2]*f3;
          h3 = recip[2][0]*f1 + recip[2][1]*f2 + recip[2][2]*f3;

          f[i][0] -= copm[k+m+1]*h1;
          f[i][1] -= copm[k+m+1]*h2;
          f[i][2] -= copm[k+m+1]*h3;

          for (j = 1; j < 4; j++) {
            cphid[j] = 0.0;
            cphip[j] = 0.0;
            for (j1 = 1; j1 < 4; j1++) {
              cphid[j] += ftc[j][j1]*fphid[i][j1];
              cphip[j] += ftc[j][j1]*fphip[i][j1];
            }
          }

          vxx -= 0.5*copm[k+m+1] *
            (cphid[1]*uoptp[i][m][0] + cphip[1]*uopt[i][m][0]);
          vyy -= 0.5*copm[k+m+1] *
            (cphid[2]*uoptp[i][m][1]+ cphip[2]*uopt[i][m][1]);
          vzz -= 0.5*copm[k+m+1] *
            (cphid[3]*uoptp[i][m][2]+ cphip[3]*uopt[i][m][2]);
          vxy -= 0.25*copm[k+m+1] *
            (cphid[1]*uoptp[i][m][1]+ cphip[1]*uopt[i][m][1]+
             cphid[2]*uoptp[i][m][0]+ cphip[2]*uopt[i][m][0]);
          vyz -= 0.25*copm[k+m+1] *
            (cphid[1]*uoptp[i][m][2]+ cphip[1]*uopt[i][m][2]+
             cphid[3]*uoptp[i][m][0]+ cphip[3]*uopt[i][m][0]);
          vxz -= 0.25*copm[k+m+1] *
            (cphid[2]*uoptp[i][m][2]+ cphip[2]*uopt[i][m][2]+
             cphid[3]*uoptp[i][m][1]+ cphip[3]*uopt[i][m][1]);
        }
      }
    }
  }

  // assign permanent and induced multipoles to the PME grid

  for (i = 0; i < nlocal; i++) {
    for (j = 1; j < 4; j++)
      cmp[i][j] += uinp[i][j-1];
  }

  // convert Cartesian multipoles to fractional multipoles

  cmp_to_fmp(cmp,fmp);

  // gridpre = my portion of 3d grid in brick decomp w/ ghost values
  // zeroed by zero()

  FFT_SCALAR ***gridpre = (FFT_SCALAR ***) p_kspace->zero();

  // map atoms to grid

  grid_mpole(fmp,gridpre);

  // pre-convolution operations including forward FFT
  // gridfft = my portion of complex 3d grid in FFT decomp as 1d vector

  gridfft = p_kspace->pre_convolution();

  // gridfft1 = copy of first FFT

  int nfft_owned = p_kspace->nfft_owned;
  memcpy(gridfft1,gridfft,2*nfft_owned*sizeof(FFT_SCALAR));

  // assign induced dipoles to the PME grid

  for (i = 0; i < nlocal; i++) {
    for (j = 1; j < 4; j++)
      cmp[i][j] += uind[i][j-1] - uinp[i][j-1];
  }

  // convert Cartesian multipoles to fractional multipoles

  cmp_to_fmp(cmp,fmp);

  // gridpre = my portion of 3d grid in brick decomp w/ ghost values
  // zeroed by zero()

  gridpre = (FFT_SCALAR ***) p_kspace->zero();

  // map atoms to grid

  grid_mpole(fmp,gridpre);

  // pre-convolution operations including forward FFT
  // gridfft1/2 = my portions of complex 3d grid in FFT decomp as 1d vectors

  FFT_SCALAR *gridfft2 = p_kspace->pre_convolution();

  // ---------------------
  // convolution operation
  // ---------------------

  m = n = 0;
  for (k = nzlo; k <= nzhi; k++) {
    for (j = nylo; j <= nyhi; j++) {
      for (i = nxlo; i <= nxhi; i++) {
        r1 = (i >= nhalf1) ? i-nfft1 : i;
        r2 = (j >= nhalf2) ? j-nfft2 : j;
        r3 = (k >= nhalf3) ? k-nfft3 : k;
        h1 = recip[0][0]*r1 + recip[0][1]*r2 + recip[0][2]*r3;  // matvec
        h2 = recip[1][0]*r1 + recip[1][1]*r2 + recip[1][2]*r3;
        h3 = recip[2][0]*r1 + recip[2][1]*r2 + recip[2][2]*r3;
        hsq = h1*h1 + h2*h2 + h3*h3;
        term = -pterm * hsq;
        expterm = 0.0;
        if (term > -50.0 && hsq != 0.0) {
          denom = volterm*hsq*bsmod1[i]*bsmod2[j]*bsmod3[k];
          expterm = exp(term) / denom;
          struc2 = gridfft1[n]*gridfft2[n] + gridfft1[n+1]*gridfft2[n+1];
          eterm = 0.5 * felec * expterm * struc2;
          vterm = (2.0/hsq) * (1.0-term) * eterm;
          vxx += h1*h1*vterm - eterm;
          vyy += h2*h2*vterm - eterm;
          vzz += h3*h3*vterm - eterm;
          vxy += h1*h2*vterm;
          vyz += h2*h3*vterm;
          vxz += h1*h3*vterm;
        }
        n += 2;
      }
    }
  }

  // assign only the induced dipoles to the PME grid
  // and perform the 3-D FFT forward transformation
  // NOTE: why is there no inverse FFT in this section?

  if (poltyp == DIRECT || poltyp == TCG) {

    for (i = 0; i < nlocal; i++) {
      for (j = 0; j < 10; j++)
        cmp[i][j] = 0.0;
      for (j = 1; j < 4; j++)
        cmp[i][j] = uinp[i][j-1];
    }

    // convert Cartesian multipoles to fractional multipoles

    cmp_to_fmp(cmp,fmp);

    // gridpre = my portion of 3d grid in brick decomp w/ ghost values
    // zeroed by zero()

    FFT_SCALAR ***gridpre = (FFT_SCALAR ***) p_kspace->zero();

    // map atoms to grid

    grid_mpole(fmp,gridpre);

    // pre-convolution operations including forward FFT
    // gridfft = my portion of complex 3d grid in FFT decomp as 1d vector

    FFT_SCALAR *gridfft = p_kspace->pre_convolution();

    // gridfft1 = copy of first FFT

    int nfft_owned = p_kspace->nfft_owned;
    memcpy(gridfft1,gridfft,2*nfft_owned*sizeof(FFT_SCALAR));

    // assign ??? to the PME grid

    for (i = 0; i < nlocal; i++) {
      for (j = 1; j < 4; j++)
        cmp[i][j] = uind[i][j-1];
    }

    // convert Cartesian multipoles to fractional multipoles

    cmp_to_fmp(cmp,fmp);

    // gridpre = my portion of 3d grid in brick decomp w/ ghost values

    gridpre = (FFT_SCALAR ***) p_kspace->zero();

    // map atoms to grid

    grid_mpole(fmp,gridpre);

    // pre-convolution operations including forward FFT
    // gridfft = my portion of complex 3d grid in FFT decomp as 1d vector

    FFT_SCALAR *gridfft2 = p_kspace->pre_convolution();

    // ---------------------
    // convolution operation
    // ---------------------

    m = n = 0;
    for (k = nzlo; k <= nzhi; k++) {
      for (j = nylo; j <= nyhi; j++) {
        for (i = nxlo; i <= nxhi; i++) {
          r1 = (i >= nhalf1) ? i-nfft1 : i;
          r2 = (j >= nhalf2) ? j-nfft2 : j;
          r3 = (k >= nhalf3) ? k-nfft3 : k;
          h1 = recip[0][0]*r1 + recip[0][1]*r2 + recip[0][2]*r3;  // matvec
          h2 = recip[1][0]*r1 + recip[1][1]*r2 + recip[1][2]*r3;
          h3 = recip[2][0]*r1 + recip[2][1]*r2 + recip[2][2]*r3;
          hsq = h1*h1 + h2*h2 + h3*h3;
          term = -pterm * hsq;
          expterm = 0.0;
          if (term > -50.0 && hsq != 0.0) {
            denom = volterm*hsq*bsmod1[i]*bsmod2[j]*bsmod3[k];
            expterm = exp(term) / denom;
            struc2 = gridfft1[n]*gridfft2[n] + gridfft1[n+1]*gridfft2[n+1];
            eterm = 0.5 * felec * expterm * struc2;
            vterm = (2.0/hsq) * (1.0-term) * eterm;
            vxx += h1*h1*vterm - eterm;
            vyy += h2*h2*vterm - eterm;
            vzz += h3*h3*vterm - eterm;
            vxy += h1*h2*vterm;
            vyz += h2*h3*vterm;
            vxz += h1*h3*vterm;
          }
          n += 2;
        }
      }
    }
  }

  // increment the total internal virial tensor components

  if (vflag_global) {
    virpolar[0] -= vxx;
    virpolar[1] -= vyy;
    virpolar[2] -= vzz;
    virpolar[3] -= vxy;
    virpolar[4] -= vxz;
    virpolar[5] -= vyz;
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
    _tq[0] = tq_ptr[3*i];
    _tq[1] = tq_ptr[3*i+1];
    _tq[2] = tq_ptr[3*i+2];
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
