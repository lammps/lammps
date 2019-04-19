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
   Contributing author: Andres Jaramillo-Botero
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_eff_cut.h"
#include "pair_eff_inline.h"
#include "atom.h"
#include "update.h"
#include "min.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "atom_vec_electron.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairEffCut::PairEffCut(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;

  nmax = 0;
  min_eradius = NULL;
  min_erforce = NULL;
  nextra = 4;
  pvector = new double[nextra];
}

/* ---------------------------------------------------------------------- */

PairEffCut::~PairEffCut()
{
  delete [] pvector;
  memory->destroy(min_eradius);
  memory->destroy(min_erforce);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
}

/* ---------------------------------------------------------------------- */

void PairEffCut::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,energy;
  double eke,ecoul,epauli,errestrain,halfcoul,halfpauli;
  double fpair,fx,fy,fz;
  double e1rforce,e2rforce,e1rvirial,e2rvirial;
  double s_fpair, s_e1rforce, s_e2rforce;
  double ecp_epauli, ecp_fpair, ecp_e1rforce, ecp_e2rforce;
  double rsq,rc;
  int *ilist,*jlist,*numneigh,**firstneigh;

  energy = eke = epauli = ecp_epauli = ecoul = errestrain = 0.0;
  // pvector = [KE, Pauli, ecoul, radial_restraint]
  for (i=0; i<4; i++) pvector[i] = 0.0;

  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  double *erforce = atom->erforce;
  double *eradius = atom->eradius;
  int *spin = atom->spin;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // add electron wavefuntion kinetic energy (not pairwise)

    if (abs(spin[i])==1 || spin[i]==2) {
      // reset energy and force temp variables
      eke = epauli = ecoul = 0.0;
      fpair = e1rforce = e2rforce = 0.0;
      s_fpair = 0.0;

      KinElec(eradius[i],&eke,&e1rforce);

      // Fixed-core
      if (spin[i] == 2) {
        // KE(2s)+Coul(1s-1s)+Coul(2s-nuclei)+Pauli(2s)
        eke *= 2;
        ElecNucElec(q[i],0.0,eradius[i],&ecoul,&fpair,&e1rforce);
        ElecNucElec(q[i],0.0,eradius[i],&ecoul,&fpair,&e1rforce);
        ElecElecElec(0.0,eradius[i],eradius[i],&ecoul,&fpair,&e1rforce,&e2rforce);

        // opposite spin electron interactions
        PauliElecElec(0,0.0,eradius[i],eradius[i],
            &epauli,&s_fpair,&e1rforce,&e2rforce);

        // fix core electron size, i.e. don't contribute to ervirial
        e2rforce = e1rforce = 0.0;
      }

      // apply unit conversion factors
      eke *= hhmss2e;
      ecoul *= qqrd2e;
      fpair *= qqrd2e;
      epauli *= hhmss2e;
      s_fpair *= hhmss2e;
      e1rforce *= hhmss2e;

      // Sum up contributions
      energy = eke + epauli + ecoul;
      fpair = fpair + s_fpair;

      erforce[i] += e1rforce;

      // Tally energy and compute radial atomic virial contribution
      if (evflag) {
        ev_tally_eff(i,i,nlocal,newton_pair,energy,0.0);
        if (pressure_with_evirials_flag) // iff flexible pressure flag on
          ev_tally_eff(i,i,nlocal,newton_pair,0.0,e1rforce*eradius[i]);
      }
      if (eflag_global) {
        pvector[0] += eke;
        pvector[1] += epauli;
        pvector[2] += ecoul;
      }
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      rc = sqrt(rsq);

      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {

        energy = ecoul = epauli = ecp_epauli = 0.0;
        fx = fy = fz = fpair = s_fpair = ecp_fpair = 0.0;

        double taper = sqrt(cutsq[itype][jtype]);
        double dist = rc / taper;
        double spline = cutoff(dist);
        double dspline = dcutoff(dist) / taper;

        // nucleus (i) - nucleus (j) Coul interaction

        if (spin[i] == 0 && spin[j] == 0) {
          double qxq = q[i]*q[j];

          ElecNucNuc(qxq, rc, &ecoul, &fpair);
        }

        // fixed-core (i) - nucleus (j) nuclear Coul interaction
        else if (spin[i] == 2 && spin[j] == 0) {
          double qxq = q[i]*q[j];
          e1rforce = 0.0;

          ElecNucNuc(qxq, rc, &ecoul, &fpair);
          ElecNucElec(q[j],rc,eradius[i],&ecoul,&fpair,&e1rforce);
          ElecNucElec(q[j],rc,eradius[i],&ecoul,&fpair,&e1rforce);
        }

        // nucleus (i) - fixed-core (j) nuclear Coul interaction
        else if (spin[i] == 0 && spin[j] == 2) {
          double qxq = q[i]*q[j];
          e1rforce = 0.0;

          ElecNucNuc(qxq, rc, &ecoul, &fpair);
          ElecNucElec(q[i],rc,eradius[j],&ecoul,&fpair,&e1rforce);
          ElecNucElec(q[i],rc,eradius[j],&ecoul,&fpair,&e1rforce);
        }

        // pseudo-core nucleus (i) - nucleus (j) interaction
        else if (spin[i] == 3 && spin[j] == 0) {
          double qxq = q[i]*q[j];

          ElecCoreNuc(qxq, rc, eradius[i], &ecoul, &fpair);
        }

        else if (spin[i] == 4 && spin[j] == 0) {
          double qxq = q[i]*q[j];

          ElecCoreNuc(qxq, rc, eradius[i], &ecoul, &fpair);
        }

        // nucleus (i) - pseudo-core nucleus (j) interaction
        else if (spin[i] == 0 && spin[j] == 3) {
          double qxq = q[i]*q[j];

          ElecCoreNuc(qxq, rc, eradius[j], &ecoul, &fpair);
        }

        else if (spin[i] == 0 && spin[j] == 4) {
          double qxq = q[i]*q[j];

          ElecCoreNuc(qxq, rc, eradius[j], &ecoul, &fpair);
        }

        // nucleus (i) - electron (j) Coul interaction

        else if  (spin[i] == 0 && abs(spin[j]) == 1) {
          e1rforce = 0.0;

          ElecNucElec(q[i],rc,eradius[j],&ecoul,&fpair,&e1rforce);

          e1rforce = spline * qqrd2e * e1rforce;
          erforce[j] += e1rforce;

          // Radial electron virial, iff flexible pressure flag set
          if (evflag && pressure_with_evirials_flag) {
            e1rvirial = eradius[j] * e1rforce;
            ev_tally_eff(j,j,nlocal,newton_pair,0.0,e1rvirial);
          }
        }

        // electron (i) - nucleus (j) Coul interaction

        else if (abs(spin[i]) == 1 && spin[j] == 0) {
          e1rforce = 0.0;

          ElecNucElec(q[j],rc,eradius[i],&ecoul,&fpair,&e1rforce);

          e1rforce = spline * qqrd2e * e1rforce;
          erforce[i] += e1rforce;

          // Radial electron virial, iff flexible pressure flag set
          if (evflag && pressure_with_evirials_flag) {
            e1rvirial = eradius[i] * e1rforce;
            ev_tally_eff(i,i,nlocal,newton_pair,0.0,e1rvirial);
          }
        }

        // electron (i) - electron (j) interactions

        else if (abs(spin[i]) == 1 && abs(spin[j]) == 1) {
          e1rforce = e2rforce = 0.0;
          s_e1rforce = s_e2rforce = 0.0;

          ElecElecElec(rc,eradius[i],eradius[j],&ecoul,&fpair,
                       &e1rforce,&e2rforce);
          PauliElecElec(spin[i] == spin[j],rc,eradius[i],eradius[j],
                       &epauli,&s_fpair,&s_e1rforce,&s_e2rforce);

          // Apply conversion factor
          epauli *= hhmss2e;
          s_fpair *= hhmss2e;

          e1rforce = spline * (qqrd2e * e1rforce + hhmss2e * s_e1rforce);
          erforce[i] += e1rforce;
          e2rforce = spline * (qqrd2e * e2rforce + hhmss2e * s_e2rforce);
          erforce[j] += e2rforce;

          // Radial electron virial, iff flexible pressure flag set
          if (evflag && pressure_with_evirials_flag) {
            e1rvirial = eradius[i] * e1rforce;
            e2rvirial = eradius[j] * e2rforce;
            ev_tally_eff(i,j,nlocal,newton_pair,0.0,e1rvirial+e2rvirial);
          }
        }

        // fixed-core (i) - electron (j) interactions

        else if (spin[i] == 2 && abs(spin[j]) == 1) {
          e1rforce = e2rforce = 0.0;
          s_e1rforce = s_e2rforce = 0.0;

          ElecNucElec(q[i],rc,eradius[j],&ecoul,&fpair,&e2rforce);
          ElecElecElec(rc,eradius[i],eradius[j],&ecoul,&fpair,
                         &e1rforce,&e2rforce);
          ElecElecElec(rc,eradius[i],eradius[j],&ecoul,&fpair,
                         &e1rforce,&e2rforce);
          PauliElecElec(0,rc,eradius[i],eradius[j],&epauli,
                       &s_fpair,&s_e1rforce,&s_e2rforce);
          PauliElecElec(1,rc,eradius[i],eradius[j],&epauli,
                       &s_fpair,&s_e1rforce,&s_e2rforce);

          // Apply conversion factor
          epauli *= hhmss2e;
          s_fpair *= hhmss2e;

          // only update virial for j electron
          e2rforce = spline * (qqrd2e * e2rforce + hhmss2e * s_e2rforce);
          erforce[j] += e2rforce;

          // Radial electron virial, iff flexible pressure flag set
          if (evflag && pressure_with_evirials_flag) {
            e2rvirial = eradius[j] * e2rforce;
            ev_tally_eff(j,j,nlocal,newton_pair,0.0,e2rvirial);
          }
        }

        // electron (i) - fixed-core (j) interactions

        else if (abs(spin[i]) == 1 && spin[j] == 2) {
          e1rforce = e2rforce = 0.0;
          s_e1rforce = s_e2rforce = 0.0;

          ElecNucElec(q[j],rc,eradius[i],&ecoul,&fpair,&e2rforce);
          ElecElecElec(rc,eradius[j],eradius[i],&ecoul,&fpair,
                         &e1rforce,&e2rforce);
          ElecElecElec(rc,eradius[j],eradius[i],&ecoul,&fpair,
                         &e1rforce,&e2rforce);

          PauliElecElec(0,rc,eradius[j],eradius[i],&epauli,
                       &s_fpair,&s_e1rforce,&s_e2rforce);
          PauliElecElec(1,rc,eradius[j],eradius[i],&epauli,
                       &s_fpair,&s_e1rforce,&s_e2rforce);

          // Apply conversion factor
          epauli *= hhmss2e;
          s_fpair *= hhmss2e;

          // only update virial for i electron
          e2rforce = spline * (qqrd2e * e2rforce + hhmss2e * s_e2rforce);
          erforce[i] += e2rforce;

          // add radial atomic virial, iff flexible pressure flag set
          if (evflag && pressure_with_evirials_flag) {
            e2rvirial = eradius[i] * e2rforce;
            ev_tally_eff(i,i,nlocal,newton_pair,0.0,e2rvirial);
          }
        }

        // fixed-core (i) - fixed-core (j) interactions

        else if (spin[i] == 2 && spin[j] == 2) {
          e1rforce = e2rforce = 0.0;
          s_e1rforce = s_e2rforce = 0.0;
          double qxq = q[i]*q[j];

          ElecNucNuc(qxq, rc, &ecoul, &fpair);
          ElecNucElec(q[i],rc,eradius[j],&ecoul,&fpair,&e1rforce);
          ElecNucElec(q[i],rc,eradius[j],&ecoul,&fpair,&e1rforce);
          ElecNucElec(q[j],rc,eradius[i],&ecoul,&fpair,&e1rforce);
          ElecNucElec(q[j],rc,eradius[i],&ecoul,&fpair,&e1rforce);
          ElecElecElec(rc,eradius[i],eradius[j],&ecoul,&fpair,
                         &e1rforce,&e2rforce);
          ElecElecElec(rc,eradius[i],eradius[j],&ecoul,&fpair,
                         &e1rforce,&e2rforce);
          ElecElecElec(rc,eradius[i],eradius[j],&ecoul,&fpair,
                         &e1rforce,&e2rforce);
          ElecElecElec(rc,eradius[i],eradius[j],&ecoul,&fpair,
                         &e1rforce,&e2rforce);

          PauliElecElec(0,rc,eradius[i],eradius[j],&epauli,
                       &s_fpair,&s_e1rforce,&s_e2rforce);
          PauliElecElec(1,rc,eradius[i],eradius[j],&epauli,
                       &s_fpair,&s_e1rforce,&s_e2rforce);
          epauli *= 2;
          s_fpair *= 2;

          // Apply conversion factor
          epauli *= hhmss2e;
          s_fpair *= hhmss2e;
        }

        // pseudo-core (i) - electron/fixed-core electrons (j) interactions

        else if (spin[i] == 3 && (abs(spin[j]) == 1 || spin[j] == 2)) {
          e2rforce = ecp_e2rforce = 0.0;

          if (((PAULI_CORE_D[ecp_type[itype]]) == 0.0) && ((PAULI_CORE_E[ecp_type[itype]]) == 0.0)) {
            if (abs(spin[j]) == 1) {
              ElecCoreElec(q[i],rc,eradius[i],eradius[j],&ecoul,
                          &fpair,&e2rforce);
              PauliCoreElec(rc,eradius[j],&ecp_epauli,&ecp_fpair,
                          &ecp_e2rforce,PAULI_CORE_A[ecp_type[itype]], PAULI_CORE_B[ecp_type[itype]],
                          PAULI_CORE_C[ecp_type[itype]]);
            } else { // add second s electron contribution from fixed-core
              double qxq = q[i]*q[j];
              ElecCoreNuc(qxq, rc, eradius[j], &ecoul, &fpair);
              ElecCoreElec(q[i],rc,eradius[i],eradius[j],&ecoul,
                          &fpair,&e2rforce);
              ElecCoreElec(q[i],rc,eradius[i],eradius[j],&ecoul,
                          &fpair,&e2rforce);
              PauliCoreElec(rc,eradius[j],&ecp_epauli,&ecp_fpair,
                          &ecp_e2rforce,PAULI_CORE_A[ecp_type[itype]], PAULI_CORE_B[ecp_type[itype]],
                          PAULI_CORE_C[ecp_type[itype]]);
              PauliCoreElec(rc,eradius[j],&ecp_epauli,&ecp_fpair,
                          &ecp_e2rforce,PAULI_CORE_A[ecp_type[itype]], PAULI_CORE_B[ecp_type[itype]],
                          PAULI_CORE_C[ecp_type[itype]]);
            }
          } else {
            if (abs(spin[j]) == 1) {
              ElecCoreElec(q[i],rc,eradius[i],eradius[j],&ecoul,
                          &fpair,&e2rforce);
              PauliCorePElec(rc,eradius[j],&ecp_epauli,&ecp_fpair,
                          &ecp_e2rforce,PAULI_CORE_A[ecp_type[itype]],PAULI_CORE_B[ecp_type[itype]],
                          PAULI_CORE_C[ecp_type[itype]],PAULI_CORE_D[ecp_type[itype]],PAULI_CORE_E[ecp_type[itype]]);
            } else { // add second s electron contribution from fixed-core
              double qxq = q[i]*q[j];
              ElecCoreNuc(qxq, rc, eradius[j], &ecoul, &fpair);
              ElecCoreElec(q[i],rc,eradius[i],eradius[j],&ecoul,
                          &fpair,&e2rforce);
              ElecCoreElec(q[i],rc,eradius[i],eradius[j],&ecoul,
                          &fpair,&e2rforce);
              PauliCorePElec(rc,eradius[j],&ecp_epauli,&ecp_fpair,
                          &ecp_e2rforce,PAULI_CORE_A[ecp_type[itype]], PAULI_CORE_B[ecp_type[itype]],
                          PAULI_CORE_C[ecp_type[itype]],PAULI_CORE_D[ecp_type[itype]],PAULI_CORE_E[ecp_type[itype]]);
              PauliCorePElec(rc,eradius[j],&ecp_epauli,&ecp_fpair,
                          &ecp_e2rforce,PAULI_CORE_A[ecp_type[itype]], PAULI_CORE_B[ecp_type[itype]],
                          PAULI_CORE_C[ecp_type[itype]],PAULI_CORE_D[ecp_type[itype]],PAULI_CORE_E[ecp_type[itype]]);
            }
          }

          // Apply conversion factor from Hartree to kcal/mol
          ecp_epauli *= h2e;
          ecp_fpair *= h2e;

          // only update virial for j electron
          e2rforce = spline * (qqrd2e * e2rforce + h2e * ecp_e2rforce);
          erforce[j] += e2rforce;

          // add radial atomic virial, iff flexible pressure flag set
          if (evflag && pressure_with_evirials_flag) {
            e2rvirial = eradius[j] * e2rforce;
            ev_tally_eff(j,j,nlocal,newton_pair,0.0,e2rvirial);
          }
        }

        // electron/fixed-core electrons (i) - pseudo-core (j) interactions

        else if ((abs(spin[i]) == 1 || spin[i] == 2) && spin[j] == 3) {
          e1rforce = ecp_e1rforce = 0.0;

          if (((PAULI_CORE_D[ecp_type[jtype]]) == 0.0) && ((PAULI_CORE_E[ecp_type[jtype]]) == 0.0)) {
            if (abs(spin[i]) == 1) {
              ElecCoreElec(q[j],rc,eradius[j],eradius[i],&ecoul,
                          &fpair,&e1rforce);
              PauliCoreElec(rc,eradius[i],&ecp_epauli,&ecp_fpair,
                          &ecp_e1rforce,PAULI_CORE_A[ecp_type[jtype]],PAULI_CORE_B[ecp_type[jtype]],
                          PAULI_CORE_C[ecp_type[jtype]]);
            } else {
              double qxq = q[i]*q[j];
              ElecCoreNuc(qxq,rc,eradius[i],&ecoul,&fpair);
              ElecCoreElec(q[j],rc,eradius[j],eradius[i],&ecoul,
                          &fpair,&e1rforce);
              ElecCoreElec(q[j],rc,eradius[j],eradius[i],&ecoul,
                          &fpair,&e1rforce);
              PauliCoreElec(rc,eradius[i],&ecp_epauli,&ecp_fpair,
                          &ecp_e1rforce,PAULI_CORE_A[ecp_type[jtype]], PAULI_CORE_B[ecp_type[jtype]],
                          PAULI_CORE_C[ecp_type[jtype]]);
              PauliCoreElec(rc,eradius[i],&ecp_epauli,&ecp_fpair,
                          &ecp_e1rforce,PAULI_CORE_A[ecp_type[jtype]], PAULI_CORE_B[ecp_type[jtype]],
                          PAULI_CORE_C[ecp_type[jtype]]);
            }
          } else {
            if (abs(spin[i]) == 1) {
              ElecCoreElec(q[j],rc,eradius[j],eradius[i],&ecoul,
                          &fpair,&e1rforce);
              PauliCorePElec(rc,eradius[i],&ecp_epauli,&ecp_fpair,
                          &ecp_e1rforce,PAULI_CORE_A[ecp_type[jtype]],PAULI_CORE_B[ecp_type[jtype]],
                          PAULI_CORE_C[ecp_type[jtype]],PAULI_CORE_D[ecp_type[jtype]],PAULI_CORE_E[ecp_type[jtype]]);
            } else {
              double qxq = q[i]*q[j];
              ElecCoreNuc(qxq,rc,eradius[i],&ecoul,&fpair);
              ElecCoreElec(q[j],rc,eradius[j],eradius[i],&ecoul,
                          &fpair,&e1rforce);
              ElecCoreElec(q[j],rc,eradius[j],eradius[i],&ecoul,
                          &fpair,&e1rforce);
              PauliCorePElec(rc,eradius[i],&ecp_epauli,&ecp_fpair,
                          &ecp_e1rforce,PAULI_CORE_A[ecp_type[jtype]], PAULI_CORE_B[ecp_type[jtype]],
                          PAULI_CORE_C[ecp_type[jtype]],PAULI_CORE_D[ecp_type[jtype]],PAULI_CORE_E[ecp_type[jtype]]);
              PauliCorePElec(rc,eradius[i],&ecp_epauli,&ecp_fpair,
                          &ecp_e1rforce,PAULI_CORE_A[ecp_type[jtype]], PAULI_CORE_B[ecp_type[jtype]],
                          PAULI_CORE_C[ecp_type[jtype]],PAULI_CORE_D[ecp_type[jtype]],PAULI_CORE_E[ecp_type[jtype]]);
            }
          }

          // Apply conversion factor from Hartree to kcal/mol
          ecp_epauli *= h2e;
          ecp_fpair *= h2e;

          // only update virial for j electron
          e1rforce = spline * (qqrd2e * e1rforce + h2e * ecp_e1rforce);
          erforce[i] += e1rforce;

          // add radial atomic virial, iff flexible pressure flag set
          if (evflag && pressure_with_evirials_flag) {
            e1rvirial = eradius[i] * e1rforce;
            ev_tally_eff(i,i,nlocal,newton_pair,0.0,e1rvirial);
          }
        }

        // pseudo-core (i) - pseudo-core (j) interactions

        else if (spin[i] == 3 && spin[j] == 3) {
          double qxq = q[i]*q[j];

          ElecCoreCore(qxq,rc,eradius[i],eradius[j],&ecoul,&fpair);
        }

        // Apply Coulomb conversion factor for all cases
        ecoul *= qqrd2e;
        fpair *= qqrd2e;

        // Sum up energy and force contributions
        epauli += ecp_epauli;
        energy = ecoul + epauli;
        fpair = fpair + s_fpair + ecp_fpair;

        // Apply cutoff spline
        fpair = fpair * spline - energy * dspline;
        energy = spline * energy;

        // Tally cartesian forces
        SmallRForce(delx,dely,delz,rc,fpair,&fx,&fy,&fz);
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        if (newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
        }

        // Tally energy (in ecoul) and compute normal pressure virials
        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,0.0,
                             energy,fx,fy,fz,delx,dely,delz);
        if (eflag_global) {
          if (newton_pair) {
            pvector[1] += spline * epauli;
            pvector[2] += spline * ecoul;
          }
          else {
            halfpauli = 0.5 * spline * epauli;
            halfcoul = 0.5 * spline * ecoul;
            if (i < nlocal) {
              pvector[1] += halfpauli;
              pvector[2] += halfcoul;
            }
            if (j < nlocal) {
              pvector[1] += halfpauli;
              pvector[2] += halfcoul;
            }
          }
        }

      }
    }

    // limit electron stifness (size) for periodic systems, to max=half-box-size

    if (abs(spin[i]) == 1 && limit_eradius_flag) {
      double half_box_length=0, dr, kfactor=hhmss2e*1.0;
      e1rforce = errestrain = 0.0;

      if (domain->xperiodic == 1 || domain->yperiodic == 1 ||
          domain->zperiodic == 1) {
        delx = domain->boxhi[0]-domain->boxlo[0];
        dely = domain->boxhi[1]-domain->boxlo[1];
        delz = domain->boxhi[2]-domain->boxlo[2];
        half_box_length = 0.5 * MIN(delx, MIN(dely, delz));
        if (eradius[i] > half_box_length) {
          dr = eradius[i]-half_box_length;
          errestrain=0.5*kfactor*dr*dr;
          e1rforce=-kfactor*dr;
          if (eflag_global) pvector[3] += errestrain;

          erforce[i] += e1rforce;

          // Tally radial restrain energy and add radial restrain virial
          if (evflag) {
            ev_tally_eff(i,i,nlocal,newton_pair,errestrain,0.0);
            if (pressure_with_evirials_flag)  // flexible electron pressure
              ev_tally_eff(i,i,nlocal,newton_pair,0.0,eradius[i]*e1rforce);
          }
        }
      }
    }

  }
  if (vflag_fdotr) {
    virial_fdotr_compute();
    if (pressure_with_evirials_flag) virial_eff_compute();
  }
}

/* ----------------------------------------------------------------------
   eff-specific contribution to global virial
------------------------------------------------------------------------- */

void PairEffCut::virial_eff_compute()
{
  double *eradius = atom->eradius;
  double *erforce = atom->erforce;
  double e_virial;
  int *spin = atom->spin;

  // sum over force on all particles including ghosts

  if (neighbor->includegroup == 0) {
    int nall = atom->nlocal + atom->nghost;
    for (int i = 0; i < nall; i++) {
      if (spin[i]) {
        e_virial = erforce[i]*eradius[i]/3;
        virial[0] += e_virial;
        virial[1] += e_virial;
        virial[2] += e_virial;
      }
    }

  // neighbor includegroup flag is set
  // sum over force on initial nfirst particles and ghosts

  } else {
    int nall = atom->nfirst;
    for (int i = 0; i < nall; i++) {
      if (spin[i]) {
        e_virial = erforce[i]*eradius[i]/3;
        virial[0] += e_virial;
        virial[1] += e_virial;
        virial[2] += e_virial;
      }
    }

    nall = atom->nlocal + atom->nghost;
    for (int i = atom->nlocal; i < nall; i++) {
      if (spin[i]) {
        e_virial = erforce[i]*eradius[i]/3;
        virial[0] += e_virial;
        virial[1] += e_virial;
        virial[2] += e_virial;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into per-atom accumulators
   for virial radial electronic contributions
------------------------------------------------------------------------- */

void PairEffCut::ev_tally_eff(int i, int j, int nlocal, int newton_pair,
                              double energy, double e_virial)
{
  double energyhalf;
  double partial_evirial = e_virial/3.0;
  double half_partial_evirial = partial_evirial/2;

  int *spin = atom->spin;

  if (eflag_either) {
    if (eflag_global) {
      if (newton_pair)
        eng_coul += energy;
      else {
        energyhalf = 0.5*energy;
        if (i < nlocal)
          eng_coul += energyhalf;
        if (j < nlocal)
          eng_coul += energyhalf;
      }
    }
    if (eflag_atom) {
      if (newton_pair || i < nlocal) eatom[i] += 0.5 * energy;
      if (newton_pair || j < nlocal) eatom[j] += 0.5 * energy;
    }
  }

  if (vflag_either) {
    if (vflag_global) {
      if (spin[i] && i < nlocal) {
        virial[0] += half_partial_evirial;
        virial[1] += half_partial_evirial;
        virial[2] += half_partial_evirial;
      }
      if (spin[j] && j < nlocal) {
        virial[0] += half_partial_evirial;
        virial[1] += half_partial_evirial;
        virial[2] += half_partial_evirial;
      }
    }
    if (vflag_atom) {
      if (spin[i]) {
        if (newton_pair || i < nlocal) {
          vatom[i][0] += half_partial_evirial;
          vatom[i][1] += half_partial_evirial;
          vatom[i][2] += half_partial_evirial;
        }
      }
      if (spin[j]) {
        if (newton_pair || j < nlocal) {
          vatom[j][0] += half_partial_evirial;
          vatom[j][1] += half_partial_evirial;
          vatom[j][2] += half_partial_evirial;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairEffCut::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
}

/* ---------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairEffCut::settings(int narg, char **arg)
{
  if (narg < 1)
    error->all(FLERR,"Illegal pair_style command");

  // Defaults ECP parameters for C (radius=0.154)
  PAULI_CORE_A[6] = 22.721015;
  PAULI_CORE_B[6] = 0.728733;
  PAULI_CORE_C[6] = 1.103199;
  PAULI_CORE_D[6] = 17.695345;
  PAULI_CORE_E[6] = 6.693621;

  // Defaults ECP parameters for N (radius=0.394732)
  PAULI_CORE_A[7] = 16.242367;
  PAULI_CORE_B[7] = 0.602818;
  PAULI_CORE_C[7] = 1.081856;
  PAULI_CORE_D[7] = 7.150803;
  PAULI_CORE_E[7] = 5.351936;

  // Defaults p-element ECP parameters for Oxygen (radius=0.15)
  PAULI_CORE_A[8] = 29.5185;
  PAULI_CORE_B[8] = 0.32995;
  PAULI_CORE_C[8] = 1.21676;
  PAULI_CORE_D[8] = 11.98757;
  PAULI_CORE_E[8] = 3.073417;

  // Defaults ECP parameters for Al (radius=1.660)
  PAULI_CORE_A[13] = 0.486;
  PAULI_CORE_B[13] = 1.049;
  PAULI_CORE_C[13] = 0.207;
  PAULI_CORE_D[13] = 0.0;
  PAULI_CORE_E[13] = 0.0;

  // Defaults ECP parameters for Si (radius=1.691)
  PAULI_CORE_A[14] = 0.320852;
  PAULI_CORE_B[14] = 2.283269;
  PAULI_CORE_C[14] = 0.814857;
  PAULI_CORE_D[14] = 0.0;
  PAULI_CORE_E[14] = 0.0;

  cut_global = force->numeric(FLERR,arg[0]);
  limit_eradius_flag = 0;
  pressure_with_evirials_flag = 0;

  int atype;
  int iarg = 1;
  int ecp_found = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"limit/eradius") == 0) {
      limit_eradius_flag = 1;
      iarg += 1;
    }
    else if (strcmp(arg[iarg],"pressure/evirials") == 0) {
      pressure_with_evirials_flag = 1;
      iarg += 1;
    }
    else if (strcmp(arg[iarg],"ecp") == 0) {
      iarg += 1;
      while (iarg < narg) {
        atype = force->inumeric(FLERR,arg[iarg]);
        if (strcmp(arg[iarg+1],"C") == 0) ecp_type[atype] = 6;
        else if (strcmp(arg[iarg+1],"N") == 0) ecp_type[atype] = 7;
        else if (strcmp(arg[iarg+1],"O") == 0) ecp_type[atype] = 8;
        else if (strcmp(arg[iarg+1],"Al") == 0) ecp_type[atype] = 13;
        else if (strcmp(arg[iarg+1],"Si") == 0) ecp_type[atype] = 14;
        else error->all(FLERR, "Note: there are no default parameters for this atom ECP\n");
        iarg += 2;
        ecp_found = 1;
      }
    }
  }

  if (!ecp_found && atom->ecp_flag)
    error->all(FLERR,"Need to specify ECP type on pair_style command");

  // Need to introduce 2 new constants w/out changing update.cpp
  if (force->qqr2e==332.06371) {        // i.e. Real units chosen
    h2e = 627.509;                      // hartree->kcal/mol
    hhmss2e = 175.72044219620075;       // hartree->kcal/mol * (Bohr->Angstrom)^2
  } else if (force->qqr2e==1.0) {        // electron units
    h2e = 1.0;
    hhmss2e = 1.0;
  } else error->all(FLERR,"Check your units");

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEffCut::init_style()
{
  // error and warning checks

  if (!atom->q_flag || !atom->spin_flag ||
      !atom->eradius_flag || !atom->erforce_flag)
    error->all(FLERR,"Pair eff/cut requires atom attributes "
               "q, spin, eradius, erforce");

  // add hook to minimizer for eradius and erforce

  if (update->whichflag == 2)
    update->minimize->request(this,1,0.01);

  // make sure to use the appropriate timestep when using real units

  if (update->whichflag == 1) {
    if (force->qqr2e == 332.06371 && update->dt == 1.0)
      error->all(FLERR,"You must lower the default real units timestep for pEFF ");
  }

  // need a half neigh list and optionally a granular history neigh list

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type electron pairs (ECP-only)
------------------------------------------------------------------------- */

void PairEffCut::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if ((strcmp(arg[0],"*") == 0) || (strcmp(arg[1],"*") == 0)) {
    int ilo,ihi,jlo,jhi;
    force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
    force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

    double cut_one = cut_global;
    if (narg == 3) cut_one = force->numeric(FLERR,arg[2]);

    int count = 0;
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo,i); j <= jhi; j++) {
        cut[i][j] = cut_one;
        setflag[i][j] = 1;
        count++;
      }
    }
    if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
  } else {
    int ecp;
    ecp = force->inumeric(FLERR,arg[0]);
    if (strcmp(arg[1],"s") ==0) {
      PAULI_CORE_A[ecp_type[ecp]] = force->numeric(FLERR,arg[2]);
      PAULI_CORE_B[ecp_type[ecp]] = force->numeric(FLERR,arg[3]);
      PAULI_CORE_C[ecp_type[ecp]] = force->numeric(FLERR,arg[4]);
      PAULI_CORE_D[ecp_type[ecp]] = 0.0;
      PAULI_CORE_E[ecp_type[ecp]] = 0.0;
    } else if (strcmp(arg[1],"p") ==0) {
      PAULI_CORE_A[ecp_type[ecp]] = force->numeric(FLERR,arg[2]);
      PAULI_CORE_B[ecp_type[ecp]] = force->numeric(FLERR,arg[3]);
      PAULI_CORE_C[ecp_type[ecp]] = force->numeric(FLERR,arg[4]);
      PAULI_CORE_D[ecp_type[ecp]] = force->numeric(FLERR,arg[5]);
      PAULI_CORE_E[ecp_type[ecp]] = force->numeric(FLERR,arg[6]);
    } else error->all(FLERR,"Illegal pair_coeff command");
  }
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairEffCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairEffCut::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) fwrite(&cut[i][j],sizeof(double),1,fp);
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairEffCut::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) fread(&cut[i][j],sizeof(double),1,fp);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairEffCut::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairEffCut::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   returns pointers to the log() of electron radius and corresponding force
   minimizer operates on log(radius) so radius never goes negative
   these arrays are stored locally by pair style
------------------------------------------------------------------------- */

void PairEffCut::min_xf_pointers(int /*ignore*/, double **xextra, double **fextra)
{
  // grow arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(min_eradius);
    memory->destroy(min_erforce);
    nmax = atom->nmax;
    memory->create(min_eradius,nmax,"pair:min_eradius");
    memory->create(min_erforce,nmax,"pair:min_erforce");
  }

  *xextra = min_eradius;
  *fextra = min_erforce;
}

/* ----------------------------------------------------------------------
   minimizer requests the log() of electron radius and corresponding force
   calculate and store in min_eradius and min_erforce
------------------------------------------------------------------------- */

void PairEffCut::min_xf_get(int /*ignore*/)
{
  double *eradius = atom->eradius;
  double *erforce = atom->erforce;
  int *spin = atom->spin;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (spin[i]) {
      min_eradius[i] = log(eradius[i]);
      min_erforce[i] = eradius[i]*erforce[i];
    } else min_eradius[i] = min_erforce[i] = 0.0;
}

/* ----------------------------------------------------------------------
   minimizer has changed the log() of electron radius
   propagate the change back to eradius
------------------------------------------------------------------------- */

void PairEffCut::min_x_set(int /*ignore*/)
{
  double *eradius = atom->eradius;
  int *spin = atom->spin;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (spin[i]) eradius[i] = exp(min_eradius[i]);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairEffCut::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  bytes += 2 * nmax * sizeof(double);
  return bytes;
}
