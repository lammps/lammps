// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_amoeba.h"

#include "amoeba_convolution.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_store_atom.h"
#include "math_const.h"
#include "math_special.h"
#include "my_page.h"
#include "neigh_list.h"
#include "timer.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

using MathSpecial::cube;

enum{INDUCE,RSD,SETUP_AMOEBA,SETUP_HIPPO,KMPOLE,AMGROUP};   // forward comm
enum{FIELD,ZRSD,TORQUE,UFLD};                               // reverse comm
enum{VDWL,REPULSE,QFER,DISP,MPOLE,POLAR,USOLV,DISP_LONG,MPOLE_LONG,POLAR_LONG};
enum{MUTUAL,OPT,TCG,DIRECT};
enum{GEAR,ASPC,LSQR};
enum{BUILD,APPLY};
enum{GORDON1,GORDON2};

#define DEBYE 4.80321    // conversion factor from q-Angs (real units) to Debye

/* ----------------------------------------------------------------------
   induce = induced dipole moments via pre-conditioned CG solver
   adapted from Tinker induce0a() routine
------------------------------------------------------------------------- */

void PairAmoeba::induce()
{
  bool done;
  int i,j,m,itype,iter;
  double polmin;
  double eps,epsold;
  double epsd,epsp;
  double udsum,upsum;
  double a,ap,b,bp;
  double sum,sump,term;
  double reduce[4],allreduce[4];

  // set cutoffs, taper coeffs, and PME params

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

  // reverse comm to sum field,fieldp from ghost atoms to owned atoms

  crstyle = FIELD;
  comm->reverse_comm(this);

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
  // NOTE: could rewrite these loops to avoid allocating
  //       uopt,uoptp with a optorder+1 dimension, just optorder
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

      crstyle = FIELD;
      comm->reverse_comm(this);

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

    crstyle = FIELD;
    comm->reverse_comm(this);

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

      crstyle = FIELD;
      comm->reverse_comm(this);

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
      // also commented out in induce.f of Tinker
      // if (eps > epsold) done = true;
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

    if (iter >= politer || eps > epsold)
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
   ulspred = induced dipole prediction coeffs

   ulspred uses standard extrapolation or a least squares fit
   to set coefficients of an induced dipole predictor polynomial
   literature references:

   J. Kolafa, "Time-Reversible Always Stable Predictor-Corrector
   Method for Molecular Dynamics of Polarizable Molecules", Journal
   of Computational Chemistry, 25, 335-342 (2004)

   W. Wang and R. D. Skeel, "Fast Evaluation of Polarizable Forces",
   Journal of Chemical Physics, 123, 164107 (2005)
------------------------------------------------------------------------- */

void PairAmoeba::ulspred()
{
  int i,j,k,m;
  double coeff,udk,upk;
  double amax,apmax;

  // set the Gear predictor binomial coefficients

  if (polpred == GEAR) {
    for (i = 0; i < nualt; i++) {
      coeff = gear[i];
      bpred[i] = coeff;
      bpredp[i] = coeff;
      bpreds[i] = coeff;
      bpredps[i] = coeff;
    }

  // set always stable predictor-corrector (ASPC) coefficients

  } else if (polpred == ASPC) {
    for (i = 0; i < nualt; i++) {
      coeff = aspc[i];
      bpred[i] = coeff;
      bpredp[i] = coeff;
      bpreds[i] = coeff;
      bpredps[i] = coeff;
    }

  // derive normal equations corresponding to least squares fit

  } else if (polpred == LSQR) {
    double ***udalt = fixudalt->tstore;
    double ***upalt = fixupalt->tstore;

    for (k = 0; k < nualt; k++) {
      b_ualt[k] = 0.0;
      bp_ualt[k] = 0.0;
      for (m = k; m < nualt; m++) {
        c_ualt[k][m] = 0.0;
        cp_ualt[k][m] = 0.0;
      }
    }

    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++) {
      for (j = 0; j < 3; j++) {
        for (k = 0; k < nualt; k++) {
          udk = udalt[i][k][j];
          upk = upalt[i][k][j];
          for (m = k; m < nualt; m++) {
            c_ualt[k][m] += udk*udalt[i][m][j];
            cp_ualt[k][m] += upk*upalt[i][m][j];
          }
        }
      }
    }

    i = 0;
    for (k = 1; k < nualt; k++) {
      b_ualt[k-1] = c_ualt[0][k];
      bp_ualt[k-1] = cp_ualt[0][k];
      for (m = k; m < nualt; m++) {
        a_ualt[i] = c_ualt[k][m];
        ap_ualt[i] = cp_ualt[k][m];
        i++;
      }
    }

    // check for nonzero coefficients and solve normal equations

    k = nualt - 1;
    amax = 0.0;
    apmax = 0.0;
    for (i = 0; i < k*(k+1)/2; i++) {
      amax = MAX(amax,a_ualt[i]);
      apmax = MAX(apmax,ap_ualt[i]);
    }
    if (amax != 0.0) cholesky(nualt-1,a_ualt,b_ualt);
    if (apmax != 0.0) cholesky(nualt-1,ap_ualt,bp_ualt);

    // transfer the final solution to the coefficient vector

    for (k = 0; k < nualt-1; k++) {
      bpred[k] = b_ualt[k];
      bpredp[k] = bp_ualt[k];
      bpreds[k] = b_ualt[k];
      bpredps[k] = bp_ualt[k];
    }
    bpred[nualt-1] = 0.0;
    bpredp[nualt-1] = 0.0;
    bpreds[nualt-1] = 0.0;
    bpredps[nualt-1] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   ufield0c = mutual induction via Ewald sum
   ufield0c computes the mutual electrostatic field due to
   induced dipole moments via Ewald summation
------------------------------------------------------------------------- */

void PairAmoeba::ufield0c(double **field, double **fieldp)
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

  double time0, time1, time2;
  if (timer->has_sync()) MPI_Barrier(world);
  time0 = platform::walltime();

  // get the real space portion of the mutual field

  if (polar_rspace_flag) umutual2b(field,fieldp);
  time1 = platform::walltime();

  // get the reciprocal space part of the mutual field

  if (polar_kspace_flag) umutual1(field,fieldp);
  time2 = platform::walltime();

  // add the self-energy portion of the mutual field

  term = (4.0/3.0) * aewald*aewald*aewald / MY_PIS;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < 3; j++) {
      field[i][j] += term*uind[i][j];
      fieldp[i][j] += term*uinp[i][j];
    }
  }

  // accumulate timing information

  time_mutual_rspace += time1 - time0;
  time_mutual_kspace += time2 - time1;
}

/* ----------------------------------------------------------------------
   uscale0b = dipole preconditioner via neigh list
   uscale0b builds and applies a preconditioner for the conjugate
   gradient induced dipole solver using a neighbor pair list
------------------------------------------------------------------------- */

void PairAmoeba::uscale0b(int mode, double **rsd, double **rsdp,
                          double **zrsd, double **zrsdp)
{
  int i,j,itype,jtype,iclass,jclass,igroup,jgroup;
  int ii,jj;
  double xi,yi,zi;
  double xr,yr,zr;
  double r,r2,rr3,rr5;
  double pdi,pti;
  double polmin;
  double poli,polik;
  double alphai,alphak;
  double damp,expdamp;
  double pgamma;
  double scale3,scale5;
  double m1,m2,m3;
  double m4,m5,m6;
  double factor_uscale,factor_wscale;
  double dmpik[5];

  // owned atoms

  double **x = atom->x;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  // neighbor list info

  int inum,jnum;
  int *ilist,*jlist;
  double *pclist;

  inum = list->inum;
  ilist = list->ilist;

  // ------------------------------------------------
  // apply the preconditioning matrix to the current residual
  // ------------------------------------------------

  if (mode == APPLY) {

    // use diagonal preconditioner elements as first approximation

    polmin = 0.00000001;
    for (i = 0; i < nlocal; i++) {
      itype = amtype[i];
      poli = udiag * MAX(polmin,polarity[itype]);
      for (j = 0; j < 3; j++) {
        zrsd[i][j] = poli * rsd[i][j];
        zrsdp[i][j] = poli * rsdp[i][j];
      }
    }

    // zero zrsd,zrsdp for ghost atoms only

    for (i = nlocal; i < nall; i++) {
      for (j = 0; j < 3; j++) {
        zrsd[i][j] = 0.0;
        zrsdp[i][j] = 0.0;
      }
    }

    // use the off-diagonal preconditioner elements in second phase

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      jlist = firstneigh_precond[i];
      jnum = numneigh_precond[i];
      pclist = firstneigh_pcpc[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK15;

        m1 = pclist[0];
        m2 = pclist[1];
        m3 = pclist[2];
        m4 = pclist[3];
        m5 = pclist[4];
        m6 = pclist[5];

        zrsd[i][0] += m1*rsd[j][0] + m2*rsd[j][1] + m3*rsd[j][2];
        zrsd[i][1] += m2*rsd[j][0] + m4*rsd[j][1] + m5*rsd[j][2];
        zrsd[i][2] += m3*rsd[j][0] + m5*rsd[j][1] + m6*rsd[j][2];
        zrsd[j][0] += m1*rsd[i][0] + m2*rsd[i][1] + m3*rsd[i][2];
        zrsd[j][1] += m2*rsd[i][0] + m4*rsd[i][1] + m5*rsd[i][2];
        zrsd[j][2] += m3*rsd[i][0] + m5*rsd[i][1] + m6*rsd[i][2];
        zrsdp[i][0] += m1*rsdp[j][0] + m2*rsdp[j][1] + m3*rsdp[j][2];
        zrsdp[i][1] += m2*rsdp[j][0] + m4*rsdp[j][1] + m5*rsdp[j][2];
        zrsdp[i][2] += m3*rsdp[j][0] + m5*rsdp[j][1] + m6*rsdp[j][2];
        zrsdp[j][0] += m1*rsdp[i][0] + m2*rsdp[i][1] + m3*rsdp[i][2];
        zrsdp[j][1] += m2*rsdp[i][0] + m4*rsdp[i][1] + m5*rsdp[i][2];
        zrsdp[j][2] += m3*rsdp[i][0] + m5*rsdp[i][1] + m6*rsdp[i][2];

        pclist += 6;
      }
    }

    return;
  }

  // ------------------------------------------------
  // build the off-diagonal elements of preconditioning matrix
  // ------------------------------------------------

  dpage_pcpc->reset();

  // determine the off-diagonal elements of the preconditioner

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = amtype[i];
    iclass = amtype2class[itype];
    igroup = amgroup[i];

    jlist = firstneigh_precond[i];
    jnum = numneigh_precond[i];

    xi = x[i][0];
    yi = x[i][1];
    zi = x[i][2];

    poli = polarity[itype];
    if (amoeba) {
      pdi = pdamp[itype];
      pti = thole[itype];
    } else {
      alphai = palpha[iclass];
    }

    // evaluate all sites in induce neigh list, no cutoff
    // store results in plist for re-use in APPLY

    pclist = dpage_pcpc->get(6*jnum);
    firstneigh_pcpc[i] = pclist;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_wscale = special_polar_wscale[sbmask15(j)];
      j &= NEIGHMASK15;

      xr = x[j][0] - xi;
      yr = x[j][1] - yi;
      zr = x[j][2] - zi;
      r2 = xr*xr + yr* yr + zr*zr;
      r = sqrt(r2);

      jtype = amtype[j];
      jclass = amtype2class[jtype];
      jgroup = amgroup[j];

      if (igroup == jgroup) factor_uscale = polar_uscale;
      else factor_uscale = 1.0;

      if (amoeba) {
        scale3 = factor_uscale;
        scale5 = factor_uscale;
        damp = pdi * pdamp[jtype];
        if (damp != 0.0) {
          pgamma = MIN(pti,thole[jtype]);
          damp = -pgamma * cube(r/damp);
          if (damp > -50.0) {
            expdamp = exp(damp);
            scale3 *= 1.0 - expdamp;
            scale5 *= 1.0 - expdamp*(1.0-damp);
          }
        }
      } else {
        alphak = palpha[jclass];
        dampmut(r,alphai,alphak,dmpik);
        scale3 = factor_wscale * dmpik[2];
        scale5 = factor_wscale * dmpik[4];
      }

      polik = poli * polarity[jtype];
      rr3 = scale3 * polik / (r*r2);
      rr5 = 3.0 * scale5 * polik / (r*r2*r2);

      pclist[0] = rr5*xr*xr - rr3;
      pclist[1] = rr5*xr*yr;
      pclist[2] = rr5*xr*zr;
      pclist[3] = rr5*yr*yr - rr3;
      pclist[4] = rr5*yr*zr;
      pclist[5] = rr5*zr*zr - rr3;
      pclist += 6;
    }
  }
}

/* ----------------------------------------------------------------------
   dfield0c =  direct induction via Ewald sum
   dfield0c computes the mutual electrostatic field due to
   permanent multipole moments via Ewald summation
------------------------------------------------------------------------- */

void PairAmoeba::dfield0c(double **field, double **fieldp)
{
  int i,j;
  double term;

  // zero out field,fieldp for owned and ghost atoms

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  for (i = 0; i < nall; i++) {
    for (j = 0; j < 3; j++) {
      field[i][j] = 0.0;
      fieldp[i][j] = 0.0;
    }
  }

  // get the reciprocal space part of the permanent field

  double time0, time1, time2;
  if (timer->has_sync()) MPI_Barrier(world);
  time0 = platform::walltime();

  if (polar_kspace_flag) udirect1(field);
  time1 = platform::walltime();

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < 3; j++) {
      fieldp[i][j] = field[i][j];
    }
  }

  // get the real space portion of the permanent field

  if (polar_rspace_flag) udirect2b(field,fieldp);
  time2 = platform::walltime();

  // get the self-energy portion of the permanent field

  term = (4.0/3.0) * aewald*aewald*aewald / MY_PIS;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < 3; j++) {
      field[i][j] += term*rpole[i][j+1];
      fieldp[i][j] += term*rpole[i][j+1];
    }
  }

  // accumulate timing information

  time_direct_kspace += time1 - time0;
  time_direct_rspace += time2 - time1;
}

/* ----------------------------------------------------------------------
   umutual1 = Ewald recip mutual induced field
   umutual1 computes the reciprocal space contribution of the
   induced atomic dipole moments to the field
------------------------------------------------------------------------- */

void PairAmoeba::umutual1(double **field, double **fieldp)
{
  int i,j,k,m,n;
  int nxlo,nxhi,nylo,nyhi,nzlo,nzhi;
  double term;
  double a[3][3];  // indices not flipped vs Fortran

  // return if the Ewald coefficient is zero

  if (aewald < 1.0e-6) return;

  // convert Cartesian dipoles to fractional coordinates

  for (j = 0; j < 3; j++) {
    a[0][j] = nfft1 * recip[0][j];
    a[1][j] = nfft2 * recip[1][j];
    a[2][j] = nfft3 * recip[2][j];
  }

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < 3; j++) {
      fuind[i][j] = a[j][0]*uind[i][0] + a[j][1]*uind[i][1] + a[j][2]*uind[i][2];
      fuinp[i][j] = a[j][0]*uinp[i][0] + a[j][1]*uinp[i][1] + a[j][2]*uinp[i][2];
    }
  }

  double time0, time1;

  // gridpre = my portion of 4d grid in brick decomp w/ ghost values

  FFT_SCALAR ****gridpre = (FFT_SCALAR ****) ic_kspace->zero();

  // map 2 values to grid

  if (timer->has_sync()) MPI_Barrier(world);
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

  FFT_SCALAR ****gridpost = (FFT_SCALAR ****) ic_kspace->post_convolution();

  // get potential

  if (timer->has_sync()) MPI_Barrier(world);
  time0 = platform::walltime();

  fphi_uind(gridpost,fdip_phi1,fdip_phi2,fdip_sum_phi);

  time1 = platform::walltime();
  time_fphi_uind += (time1 - time0);

  // store fractional reciprocal potentials for OPT method

  if (poltyp == OPT) {
    for (i = 0; i < nlocal; i++) {
      for (j = 0; j < 10; j++) {
        fopt[i][optlevel][j] = fdip_phi1[i][j];
        foptp[i][optlevel][j] = fdip_phi2[i][j];
      }
    }
  }

  // convert the dipole fields from fractional to Cartesian

  for (i = 0; i < 3; i++) {
    a[0][i] = nfft1 * recip[0][i];
    a[1][i] = nfft2 * recip[1][i];
    a[2][i] = nfft3 * recip[2][i];
  }

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < 3; j++) {
      dipfield1[i][j] = a[j][0]*fdip_phi1[i][1] +
        a[j][1]*fdip_phi1[i][2] + a[j][2]*fdip_phi1[i][3];
      dipfield2[i][j] = a[j][0]*fdip_phi2[i][1] +
        a[j][1]*fdip_phi2[i][2] + a[j][2]*fdip_phi2[i][3];
    }
  }

  // increment the field at each multipole site

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < 3; j++) {
      field[i][j] -= dipfield1[i][j];
      fieldp[i][j] -= dipfield2[i][j];
    }
  }
}

/* ----------------------------------------------------------------------
   umutual2b = Ewald real mutual field via list
   umutual2b computes the real space contribution of the induced
     atomic dipole moments to the field via a neighbor list
------------------------------------------------------------------------- */

void PairAmoeba::umutual2b(double **field, double **fieldp)
{
  int i,j,m,ii,jj,jnum;
  double fid[3],fkd[3];
  double fip[3],fkp[3];
  double *uindi,*uindj,*uinpi,*uinpj;

  // neigh list

  int inum = list->inum;
  int *ilist = list->ilist;
  int *jlist;
  double *tdipdip;

  // loop over owned atoms and neighs
  // compute field terms for each pairwise interaction
  // using tdipdip values stored by udirect2b()

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    uindi = uind[i];
    uinpi = uinp[i];
    jlist = firstneigh_dipole[i];
    tdipdip = firstneigh_dipdip[i];
    jnum = numneigh_dipole[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      uindj = uind[j];
      uinpj = uinp[j];

      fid[0] = tdipdip[0]*uindj[0] + tdipdip[1]*uindj[1] + tdipdip[2]*uindj[2];
      fid[1] = tdipdip[1]*uindj[0] + tdipdip[3]*uindj[1] + tdipdip[4]*uindj[2];
      fid[2] = tdipdip[2]*uindj[0] + tdipdip[4]*uindj[1] + tdipdip[5]*uindj[2];

      fkd[0] = tdipdip[0]*uindi[0] + tdipdip[1]*uindi[1] + tdipdip[2]*uindi[2];
      fkd[1] = tdipdip[1]*uindi[0] + tdipdip[3]*uindi[1] + tdipdip[4]*uindi[2];
      fkd[2] = tdipdip[2]*uindi[0] + tdipdip[4]*uindi[1] + tdipdip[5]*uindi[2];

      fip[0] = tdipdip[0]*uinpj[0] + tdipdip[1]*uinpj[1] + tdipdip[2]*uinpj[2];
      fip[1] = tdipdip[1]*uinpj[0] + tdipdip[3]*uinpj[1] + tdipdip[4]*uinpj[2];
      fip[2] = tdipdip[2]*uinpj[0] + tdipdip[4]*uinpj[1] + tdipdip[5]*uinpj[2];

      fkp[0] = tdipdip[0]*uinpi[0] + tdipdip[1]*uinpi[1] + tdipdip[2]*uinpi[2];
      fkp[1] = tdipdip[1]*uinpi[0] + tdipdip[3]*uinpi[1] + tdipdip[4]*uinpi[2];
      fkp[2] = tdipdip[2]*uinpi[0] + tdipdip[4]*uinpi[1] + tdipdip[5]*uinpi[2];

      tdipdip += 6;

      // increment the field at each site due to this interaction

      for (m = 0; m < 3; m++) {
        field[i][m] += fid[m];
        field[j][m] += fkd[m];
        fieldp[i][m] += fip[m];
        fieldp[j][m] += fkp[m];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   udirect1 = Ewald recip direct induced field
   udirect1 computes the reciprocal space contribution of the
   permanent atomic multipole moments to the field
   since corresponding values in empole and epolar are different
------------------------------------------------------------------------- */

void PairAmoeba::udirect1(double **field)
{
  int i,j,k,m,n;
  int nhalf1,nhalf2,nhalf3;
  int nxlo,nxhi,nylo,nyhi,nzlo,nzhi;
  double r1,r2,r3;
  double h1,h2,h3;
  double volterm,denom;
  double hsq,expterm;
  double term,pterm;

  // return if the Ewald coefficient is zero

  if (aewald < 1.0e-6) return;

  pterm = (MY_PI/aewald) * (MY_PI/aewald);
  double volbox = domain->prd[0] * domain->prd[1] * domain->prd[2];
  volterm = MY_PI * volbox;

  // FFT moduli pre-computations
  // set igrid for each atom and its B-spline coeffs

  nfft1 = i_kspace->nx;
  nfft2 = i_kspace->ny;
  nfft3 = i_kspace->nz;
  bsorder = i_kspace->order;

  moduli();
  bspline_fill();

  // copy the multipole moments into local storage areas

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    cmp[i][0] = rpole[i][0];
    cmp[i][1] = rpole[i][1];
    cmp[i][2] = rpole[i][2];
    cmp[i][3] = rpole[i][3];
    cmp[i][4] = rpole[i][4];
    cmp[i][5] = rpole[i][8];
    cmp[i][6] = rpole[i][12];
    cmp[i][7] = 2.0 * rpole[i][5];
    cmp[i][8] = 2.0 * rpole[i][6];
    cmp[i][9] = 2.0 * rpole[i][9];
  }

  // convert Cartesian multipoles to fractional coordinates

  cmp_to_fmp(cmp,fmp);

  // gridpre = my portion of 3d grid in brick decomp w/ ghost values
  // zeroed by setup()

  FFT_SCALAR ***gridpre = (FFT_SCALAR ***) i_kspace->zero();

  // map multipole moments to grid

  grid_mpole(fmp,gridpre);

  // pre-convolution operations including forward FFT
  // gridfft = my 1d portion of complex 3d grid in FFT decomp

  FFT_SCALAR *gridfft = i_kspace->pre_convolution();

  // ---------------------
  // convolution operation
  // ---------------------

  nhalf1 = (nfft1+1) / 2;
  nhalf2 = (nfft2+1) / 2;
  nhalf3 = (nfft3+1) / 2;

  nxlo = i_kspace->nxlo_fft;
  nxhi = i_kspace->nxhi_fft;
  nylo = i_kspace->nylo_fft;
  nyhi = i_kspace->nyhi_fft;
  nzlo = i_kspace->nzlo_fft;
  nzhi = i_kspace->nzhi_fft;

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
        }
        qfac[m++] = expterm;
        gridfft[n] *= expterm;
        gridfft[n+1] *= expterm;
        n += 2;
      }
    }
  }

  // post-convolution operations including backward FFT
  // gridppost = my portion of 3d grid in brick decomp w/ ghost values

  FFT_SCALAR ***gridpost = (FFT_SCALAR ***) i_kspace->post_convolution();

  // get potential

  fphi_mpole(gridpost,fphi);

  // convert the field from fractional to Cartesian

  fphi_to_cphi(fphi,cphi);

  // increment the field at each multipole site

  for (i = 0; i < nlocal; i++) {
    field[i][0] -= cphi[i][1];
    field[i][1] -= cphi[i][2];
    field[i][2] -= cphi[i][3];
  }
}

/* ----------------------------------------------------------------------
   udirect2b = Ewald real direct field via list
   udirect2b computes the real space contribution of the permanent
   atomic multipole moments to the field via a neighbor list
------------------------------------------------------------------------- */

void PairAmoeba::udirect2b(double **field, double **fieldp)
{
  int i,j,m,n,ii,jj,jextra,ndip,itype,jtype,iclass,jclass,igroup,jgroup;
  double xr,yr,zr,r,r2;
  double rr1,rr2,rr3;
  double rr5,rr7;
  double rr3i,rr5i,rr7i;
  double rr3k,rr5k,rr7k;
  double rr3ik,rr5ik;
  double bfac,exp2a;
  double ci,dix,diy,diz;
  double qixx,qiyy,qizz;
  double qixy,qixz,qiyz;
  double ck,dkx,dky,dkz;
  double qkxx,qkyy,qkzz;
  double qkxy,qkxz,qkyz;
  double dir,dkr;
  double qix,qiy,qiz,qir;
  double qkx,qky,qkz,qkr;
  double corei,corek;
  double vali,valk;
  double alphai,alphak;
  double ralpha,aefac;
  double aesq2,aesq2n;
  double pdi,pti,ddi;
  double pgamma;
  double damp,expdamp;
  double scale3,scale5;
  double scale7,scalek;
  double bn[4],bcn[3];
  double fid[3],fkd[3];
  double fip[3],fkp[3];
  double dmpi[7],dmpk[7];
  double dmpik[5];
  double factor_dscale,factor_pscale,factor_uscale,factor_wscale;

  int inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;

  // owned atoms

  double **x = atom->x;
  double *pval = atom->dvector[index_pval];

  // neigh list

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  aesq2 = 2.0 * aewald * aewald;
  aesq2n = 0.0;
  if (aewald > 0.0) aesq2n = 1.0 / (MY_PIS*aewald);

  // rebuild dipole-dipole pair list and store pairwise dipole matrices
  // done one atom at a time in real-space double loop over atoms & neighs

  int *neighptr;
  double *tdipdip;

  // compute the real space portion of the Ewald summation

  ipage_dipole->reset();
  dpage_dipdip->reset();

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = amtype[i];
    iclass = amtype2class[itype];
    igroup = amgroup[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    n = ndip = 0;
    neighptr = ipage_dipole->vget();
    tdipdip = dpage_dipdip->vget();

    ci = rpole[i][0];
    dix = rpole[i][1];
    diy = rpole[i][2];
    diz = rpole[i][3];
    qixx = rpole[i][4];
    qixy = rpole[i][5];
    qixz = rpole[i][6];
    qiyy = rpole[i][8];
    qiyz = rpole[i][9];
    qizz = rpole[i][12];
    if (amoeba) {
      pdi = pdamp[itype];
      pti = thole[itype];
      ddi = dirdamp[itype];
    } else {
      corei = pcore[iclass];
      alphai = palpha[iclass];
      vali = pval[i];
    }

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
      jclass = amtype2class[jtype];
      jgroup = amgroup[j];

      if (amoeba) {
        factor_wscale = special_polar_wscale[sbmask15(jextra)];
        if (igroup == jgroup) {
          factor_pscale = special_polar_piscale[sbmask15(jextra)];
          factor_dscale = polar_dscale;
          factor_uscale = polar_uscale;
        } else {
          factor_pscale = special_polar_pscale[sbmask15(jextra)];
          factor_dscale = factor_uscale = 1.0;
        }

      } else {
        factor_wscale = special_polar_wscale[sbmask15(jextra)];
        if (igroup == jgroup) {
          factor_dscale = factor_pscale = special_polar_piscale[sbmask15(jextra)];
          factor_uscale = polar_uscale;
        } else {
          factor_dscale = factor_pscale = special_polar_pscale[sbmask15(jextra)];
          factor_uscale = 1.0;
        }
      }

      r = sqrt(r2);
      rr1 = 1.0 / r;
      rr2 = rr1 * rr1;
      rr3 = rr2 * rr1;
      rr5 = 3.0 * rr2 * rr3;
      rr7 = 5.0 * rr2 * rr5;
      ck = rpole[j][0];
      dkx = rpole[j][1];
      dky = rpole[j][2];
      dkz = rpole[j][3];
      qkxx = rpole[j][4];
      qkxy = rpole[j][5];
      qkxz = rpole[j][6];
      qkyy = rpole[j][8];
      qkyz = rpole[j][9];
      qkzz = rpole[j][12];

      // intermediates involving moments and separation distance

      dir = dix*xr + diy*yr + diz*zr;
      qix = qixx*xr + qixy*yr + qixz*zr;
      qiy = qixy*xr + qiyy*yr + qiyz*zr;
      qiz = qixz*xr + qiyz*yr + qizz*zr;
      qir = qix*xr + qiy*yr + qiz*zr;
      dkr = dkx*xr + dky*yr + dkz*zr;
      qkx = qkxx*xr + qkxy*yr + qkxz*zr;
      qky = qkxy*xr + qkyy*yr + qkyz*zr;
      qkz = qkxz*xr + qkyz*yr + qkzz*zr;
      qkr = qkx*xr + qky*yr + qkz*zr;

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

      // find the field components for Thole polarization damping

      if (amoeba) {
        scale3 = 1.0;
        scale5 = 1.0;
        scale7 = 1.0;
        damp = pdi * pdamp[jtype];
        if (damp != 0.0) {
          pgamma = MIN(ddi,dirdamp[jtype]);
          if (pgamma != 0.0) {
            damp = pgamma * pow(r/damp,1.5);
            if (damp < 50.0) {
              expdamp = exp(-damp) ;
              scale3 = 1.0 - expdamp ;
              scale5 = 1.0 - expdamp*(1.0+0.5*damp);
              scale7 = 1.0 - expdamp*(1.0+0.65*damp + 0.15*damp*damp);
            }
          } else {
            pgamma = MIN(pti,thole[jtype]);
            damp = pgamma * cube(r/damp);
            if (damp < 50.0) {
              expdamp = exp(-damp);
              scale3 = 1.0 - expdamp;
              scale5 = 1.0 - expdamp*(1.0+damp);
              scale7 = 1.0 - expdamp*(1.0+damp + 0.6*damp*damp);
            }
          }
        }

        scalek = factor_dscale;
        bcn[0] = bn[1] - (1.0-scalek*scale3)*rr3;
        bcn[1] = bn[2] - (1.0-scalek*scale5)*rr5;
        bcn[2] = bn[3] - (1.0-scalek*scale7)*rr7;
        fid[0] = -xr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) -
          bcn[0]*dkx + 2.0*bcn[1]*qkx;
        fid[1] = -yr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) -
          bcn[0]*dky + 2.0*bcn[1]*qky;
        fid[2] = -zr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) -
          bcn[0]*dkz + 2.0*bcn[1]*qkz;
        fkd[0] = xr*(bcn[0]*ci+bcn[1]*dir+bcn[2]*qir) -
          bcn[0]*dix - 2.0*bcn[1]*qix;
        fkd[1] = yr*(bcn[0]*ci+bcn[1]*dir+bcn[2]*qir) -
          bcn[0]*diy - 2.0*bcn[1]*qiy;
        fkd[2] = zr*(bcn[0]*ci+bcn[1]*dir+bcn[2]*qir) -
          bcn[0]*diz - 2.0*bcn[1]*qiz;

        scalek = factor_pscale;
        bcn[0] = bn[1] - (1.0-scalek*scale3)*rr3;
        bcn[1] = bn[2] - (1.0-scalek*scale5)*rr5;
        bcn[2] = bn[3] - (1.0-scalek*scale7)*rr7;
        fip[0] = -xr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) -
          bcn[0]*dkx + 2.0*bcn[1]*qkx;
        fip[1] = -yr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) -
          bcn[0]*dky + 2.0*bcn[1]*qky;
        fip[2] = -zr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) -
          bcn[0]*dkz + 2.0*bcn[1]*qkz;
        fkp[0] = xr*(bcn[0]*ci+bcn[1]*dir+bcn[2]*qir) -
          bcn[0]*dix - 2.0*bcn[1]*qix;
        fkp[1] = yr*(bcn[0]*ci+bcn[1]*dir+bcn[2]*qir) -
          bcn[0]*diy - 2.0*bcn[1]*qiy;
        fkp[2] = zr*(bcn[0]*ci+bcn[1]*dir+bcn[2]*qir) -
          bcn[0]*diz - 2.0*bcn[1]*qiz;

        // find terms needed later to compute mutual polarization

        if (poltyp != DIRECT) {
          scale3 = 1.0;
          scale5 = 1.0;
          damp = pdi * pdamp[jtype];
          if (damp != 0.0) {
            pgamma = MIN(pti,thole[jtype]);
            damp = pgamma * cube(r/damp);
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

      // find the field components for charge penetration damping

      } else {
        corek = pcore[jclass];
        alphak = palpha[jclass];
        valk = pval[j];
        dampdir(r,alphai,alphak,dmpi,dmpk);

        scalek = factor_dscale;
        rr3i = bn[1] - (1.0-scalek*dmpi[2])*rr3;
        rr5i = bn[2] - (1.0-scalek*dmpi[4])*rr5;
        rr7i = bn[3] - (1.0-scalek*dmpi[6])*rr7;
        rr3k = bn[1] - (1.0-scalek*dmpk[2])*rr3;
        rr5k = bn[2] - (1.0-scalek*dmpk[4])*rr5;
        rr7k = bn[3] - (1.0-scalek*dmpk[6])*rr7;
        rr3 = bn[1] - (1.0-scalek)*rr3;
        fid[0] = -xr*(rr3*corek + rr3k*valk - rr5k*dkr + rr7k*qkr) -
          rr3k*dkx + 2.0*rr5k*qkx;
        fid[1] = -yr*(rr3*corek + rr3k*valk - rr5k*dkr + rr7k*qkr) -
          rr3k*dky + 2.0*rr5k*qky;
        fid[2] = -zr*(rr3*corek + rr3k*valk - rr5k*dkr + rr7k*qkr) -
          rr3k*dkz + 2.0*rr5k*qkz;
        fkd[0] = xr*(rr3*corei + rr3i*vali + rr5i*dir + rr7i*qir) -
          rr3i*dix - 2.0*rr5i*qix;
        fkd[1] = yr*(rr3*corei + rr3i*vali + rr5i*dir + rr7i*qir) -
          rr3i*diy - 2.0*rr5i*qiy;
        fkd[2] = zr*(rr3*corei + rr3i*vali + rr5i*dir + rr7i*qir) -
          rr3i*diz - 2.0*rr5i*qiz;

        scalek = factor_pscale;
        rr3 = rr2 * rr1;
        rr3i = bn[1] - (1.0-scalek*dmpi[2])*rr3;
        rr5i = bn[2] - (1.0-scalek*dmpi[4])*rr5;
        rr7i = bn[3] - (1.0-scalek*dmpi[6])*rr7;
        rr3k = bn[1] - (1.0-scalek*dmpk[2])*rr3;
        rr5k = bn[2] - (1.0-scalek*dmpk[4])*rr5;
        rr7k = bn[3] - (1.0-scalek*dmpk[6])*rr7;
        rr3 = bn[1] - (1.0-scalek)*rr3;
        fip[0] = -xr*(rr3*corek + rr3k*valk - rr5k*dkr + rr7k*qkr) -
          rr3k*dkx + 2.0*rr5k*qkx;
        fip[1] = -yr*(rr3*corek + rr3k*valk - rr5k*dkr + rr7k*qkr) -
          rr3k*dky + 2.0*rr5k*qky;
        fip[2] = -zr*(rr3*corek + rr3k*valk - rr5k*dkr + rr7k*qkr) -
          rr3k*dkz + 2.0*rr5k*qkz;
        fkp[0] = xr*(rr3*corei + rr3i*vali + rr5i*dir + rr7i*qir) -
          rr3i*dix - 2.0*rr5i*qix;
        fkp[1] = yr*(rr3*corei + rr3i*vali + rr5i*dir + rr7i*qir) -
          rr3i*diy - 2.0*rr5i*qiy;
        fkp[2] = zr*(rr3*corei + rr3i*vali + rr5i*dir + rr7i*qir) -
          rr3i*diz - 2.0*rr5i*qiz;

        // find terms needed later to compute mutual polarization

        if (poltyp != DIRECT) {
          dampmut(r,alphai,alphak,dmpik);
          scalek = factor_wscale;
          rr3 = rr2 * rr1;
          rr3ik = bn[1] - (1.0-scalek*dmpik[2])*rr3;
          rr5ik = bn[2] - (1.0-scalek*dmpik[4])*rr5;

          neighptr[n++] = j;
          tdipdip[ndip++] = -rr3ik + rr5ik*xr*xr;
          tdipdip[ndip++] = rr5ik*xr*yr;
          tdipdip[ndip++] = rr5ik*xr*zr;
          tdipdip[ndip++] = -rr3ik + rr5ik*yr*yr;
          tdipdip[ndip++] = rr5ik*yr*zr;
          tdipdip[ndip++] = -rr3ik + rr5ik*zr*zr;
        }
      }

      // increment the field at each site due to this interaction

      for (m = 0; m < 3; m++) {
        field[i][m] += fid[m];
        field[j][m] += fkd[m];
        fieldp[i][m] += fip[m];
        fieldp[j][m] += fkp[m];
      }
    }

    firstneigh_dipole[i] = neighptr;
    firstneigh_dipdip[i] = tdipdip;
    numneigh_dipole[i] = n;
    ipage_dipole->vgot(n);
    dpage_dipdip->vgot(ndip);
  }
}

/* ----------------------------------------------------------------------
   dampmut = mutual field damping coefficents
   dampmut generates coefficients for the mutual field damping
   function for powers of the interatomic distance
------------------------------------------------------------------------- */

void PairAmoeba::dampmut(double r, double alphai, double alphak, double *dmpik)
{
  double termi,termk;
  double termi2,termk2;
  double alphai2,alphak2;
  double eps,diff;
  double expi,expk;
  double dampi,dampk;
  double dampi2,dampi3;
  double dampi4,dampi5;
  double dampk2,dampk3;

  // compute tolerance and exponential damping factors

  eps = 0.001;
  diff = fabs(alphai-alphak);
  dampi = alphai * r;
  dampk = alphak * r;
  expi = exp(-dampi);
  expk = exp(-dampk);

  // valence-valence charge penetration damping for Gordon f1 (HIPPO)

  dampi2 = dampi * dampi;
  dampi3 = dampi * dampi2;
  if (diff < eps) {
    dampi4 = dampi2 * dampi2;
    dampi5 = dampi2 * dampi3;
    dmpik[2] = 1.0 - (1.0 + dampi + 0.5*dampi2 +
                      7.0*dampi3/48.0 + dampi4/48.0)*expi;
    dmpik[4] = 1.0 - (1.0 + dampi + 0.5*dampi2 + dampi3/6.0 +
                      dampi4/24.0 + dampi5/144.0)*expi;
  } else {
    dampk2 = dampk * dampk;
    dampk3 = dampk * dampk2;
    alphai2 = alphai * alphai;
    alphak2 = alphak * alphak;
    termi = alphak2 / (alphak2-alphai2);
    termk = alphai2 / (alphai2-alphak2);
    termi2 = termi * termi;
    termk2 = termk * termk;
    dmpik[2] = 1.0 - termi2*(1.0+dampi+0.5*dampi2)*expi -
      termk2*(1.0+dampk+0.5*dampk2)*expk -
      2.0*termi2*termk*(1.0+dampi)*expi - 2.0*termk2*termi*(1.0+dampk)*expk;
    dmpik[4] = 1.0 - termi2*(1.0+dampi+0.5*dampi2 + dampi3/6.0)*expi -
      termk2*(1.0+dampk+0.5*dampk2 + dampk3/6.00)*expk -
      2.0*termi2*termk *(1.0+dampi+dampi2/3.0)*expi -
      2.0*termk2*termi *(1.0+dampk+dampk2/3.0)*expk;
  }
}

/* ----------------------------------------------------------------------
   dampdir = direct field damping coefficents
   dampdir generates coefficients for the direct field damping
   function for powers of the interatomic distance
------------------------------------------------------------------------- */

void PairAmoeba::dampdir(double r, double alphai, double alphak,
                         double *dmpi, double *dmpk)
{
  double eps,diff;
  double expi,expk;
  double dampi,dampk;
  double dampi2,dampk2;
  double dampi3,dampk3;
  double dampi4,dampk4;

  // compute tolerance and exponential damping factors

  eps = 0.001;
  diff = fabs(alphai-alphak);
  dampi = alphai * r;
  dampk = alphak * r;
  expi = exp(-dampi);
  expk = exp(-dampk);

  // core-valence charge penetration damping for Gordon f1 (HIPPO)

  dampi2 = dampi * dampi;
  dampi3 = dampi * dampi2;
  dampi4 = dampi2 * dampi2;
  dmpi[2] = 1.0 - (1.0 + dampi + 0.5*dampi2)*expi;
  dmpi[4] = 1.0 - (1.0 + dampi + 0.5*dampi2 + dampi3/6.0)*expi;
  dmpi[6] = 1.0 - (1.0 + dampi + 0.5*dampi2 + dampi3/6.0 + dampi4/30.0)*expi;
  if (diff < eps) {
    dmpk[2] = dmpi[2];
    dmpk[4] = dmpi[4];
    dmpk[6] = dmpi[6];
  } else {
    dampk2 = dampk * dampk;
    dampk3 = dampk * dampk2;
    dampk4 = dampk2 * dampk2;
    dmpk[2] = 1.0 - (1.0 + dampk + 0.5*dampk2)*expk;
    dmpk[4] = 1.0 - (1.0 + dampk + 0.5*dampk2 + dampk3/6.0)*expk;
    dmpk[6] = 1.0 - (1.0 + dampk + 0.5*dampk2 + dampk3/6.0 + dampk4/30.0)*expk;
  }
}

/* ----------------------------------------------------------------------
   cholesky = modified Cholesky linear solver
   cholesky uses a modified Cholesky method to solve the linear
     system Ax = b, returning "x" in "b"; "A" is a real symmetric
     positive definite matrix with its upper triangle (including the
     diagonal) stored by rows
   literature reference:
     R. S. Martin, G. Peters and J. H. Wilkinson, "Symmetric
     Decomposition of a Positive Definite Matrix", Numerische
     Mathematik, 7, 362-383 (1965)
------------------------------------------------------------------------- */

void PairAmoeba::cholesky(int nvar, double *a, double *b)
{
  int i,j,k;
  int ii,ij,ik,ki,kk;
  int im,jk,jm;
  double r,s,t;

  // all code in this method is exact translation from Fortran version
  // decrement pointers so can access vectors like Fortran does from 1:N

  a--;
  b--;

  // Cholesky factorization to reduce A to (L)(D)(L transpose)
  // L has a unit diagonal; store 1.0/D on the diagonal of A

  ii = 1;
  for (i = 1; i <= nvar; i++) {
    im = i - 1;
    if (i != 1) {
      ij = i;
      for (j = 1; j <= im; j++) {
        r = a[ij];
        if (j != 1) {
          ik = i;
          jk = j;
          jm = j - 1;
          for (k = 1; k <= jm; k++) {
            r = r - a[ik]*a[jk];
            ik = nvar - k + ik;
            jk = nvar - k + jk;
          }
        }
        a[ij] = r;
        ij = nvar - j + ij;
      }
    }

    r = a[ii];
    if (i != 1) {
      kk = 1;
      ik = i;
      for (k = 1; k <= im; k++) {
        s = a[ik];
        t = s * a[kk];
        a[ik] = t;
        r = r - s*t;
        ik = nvar - k + ik;
        kk = nvar - k + 1 + kk;
      }
    }

    a[ii] = 1.0 / r;
    ii = nvar - i + 1 + ii;
  }

  // solve linear equations; first solve Ly = b for y

  for (i = 1; i <= nvar; i++) {
    if (i != 1) {
      ik = i;
      im = i - 1;
      r = b[i];
      for (k = 1; k <= im; k++) {
        r = r - b[k]*a[ik];
        ik = nvar - k + ik;
      }
      b[i] = r;
    }
  }

  // finally, solve (D)(L transpose)(x) = y for x

  ii = nvar*(nvar+1)/2;
  for (j = 1; j <= nvar; j++) {
    i = nvar + 1 - j;
    r = b[i] * a[ii];
    if (j != 1) {
      im = i + 1;
      ki = ii + 1;
      for (k = im; k <= nvar; k++) {
        r = r - a[ki]*b[k];
        ki = ki + 1;
      }
    }
    b[i] = r;
    ii = ii - j - 1;
  }
}
