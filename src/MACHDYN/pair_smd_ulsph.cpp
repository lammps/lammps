// clang-format off
/* ----------------------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the MACHDYN package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */

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

#include "pair_smd_ulsph.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "smd_kernels.h"
#include "smd_material_models.h"
#include "smd_math.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace SMD_Kernels;
using namespace std;
using namespace LAMMPS_NS;
using namespace SMD_Math;

#include <Eigen/Eigen>
using namespace Eigen;

#define FORMAT1 "%60s : %g\n"
#define FORMAT2 "\n.............................. %s \n"

PairULSPH::PairULSPH(LAMMPS *lmp) :
                Pair(lmp) {

        // per-type arrays
        Q1 = nullptr;
        eos = viscosity = strength = nullptr;
        c0_type = nullptr;
        c0 = nullptr;
        Lookup = nullptr;
        artificial_stress = nullptr;
        artificial_pressure = nullptr;

        nmax = 0; // make sure no atom on this proc such that initial memory allocation is correct
        stressTensor = L = K = nullptr;
        shepardWeight = nullptr;
        smoothVel = nullptr;
        numNeighs = nullptr;
        F = nullptr;
        rho = nullptr;
        effm = nullptr;

        velocity_gradient_required = false; // turn off computation of velocity gradient by default
        density_summation = velocity_gradient = false;

        comm_forward = 18; // this pair style communicates 18 doubles to ghost atoms
        updateFlag = 0;
}

/* ---------------------------------------------------------------------- */

PairULSPH::~PairULSPH() {
        if (allocated) {
                memory->destroy(setflag);
                memory->destroy(cutsq);
                memory->destroy(Q1);
                memory->destroy(rho0);
                memory->destroy(eos);
                memory->destroy(viscosity);
                memory->destroy(strength);
                memory->destroy(c0_type);
                memory->destroy(Lookup);
                memory->destroy(artificial_pressure);
                memory->destroy(artificial_stress);

                delete[] onerad_dynamic;
                delete[] onerad_frozen;
                delete[] maxrad_dynamic;
                delete[] maxrad_frozen;

                delete[] K;
                delete[] shepardWeight;
                delete[] c0;
                delete[] smoothVel;
                delete[] stressTensor;
                delete[] L;
                delete[] numNeighs;
                delete[] F;
                delete[] rho;
                delete[] effm;

        }
}

/* ----------------------------------------------------------------------
 *
 * Re-compute mass density from scratch.
 * Only used for plain fluid SPH with no physical viscosity models.
 *
 ---------------------------------------------------------------------- */

void PairULSPH::PreCompute_DensitySummation() {
        double *radius = atom->radius;
        double **x = atom->x;
        double *rmass = atom->rmass;
        int *type = atom->type;
        int *ilist, *jlist, *numneigh;
        int **firstneigh;
        int nlocal = atom->nlocal;
        int inum, jnum, ii, jj, i, itype, jtype, j;
        double h, irad, hsq, rSq, wf;
        Vector3d dx, xi, xj;

        // set up neighbor list variables
        inum = list->inum;
        ilist = list->ilist;
        numneigh = list->numneigh;
        firstneigh = list->firstneigh;

        // zero accumulators
        for (i = 0; i < nlocal; i++) {
                rho[i] = 0.0;
                //shepardWeight[i] = 0.0;
        }

        /*
         * only recompute mass density if density summation is used.
         * otherwise, change in mass density is time-integrated
         */
        for (i = 0; i < nlocal; i++) {
                itype = type[i];
                if (setflag[itype][itype] == 1) {
                        // initialize particle density with self-contribution.
                        h = 2.0 * radius[i];
                        hsq = h * h;
                        Poly6Kernel(hsq, h, 0.0, domain->dimension, wf);
                        rho[i] = wf * rmass[i]; // / shepardWeight[i];
                        //printf("SIC to rho is %f\n", rho[i]);
                }
        }

        for (ii = 0; ii < inum; ii++) {
                i = ilist[ii];
                itype = type[i];
                jlist = firstneigh[i];
                jnum = numneigh[i];
                irad = radius[i];

                xi << x[i][0], x[i][1], x[i][2];

                for (jj = 0; jj < jnum; jj++) {
                        j = jlist[jj];
                        j &= NEIGHMASK;

                        xj << x[j][0], x[j][1], x[j][2];
                        dx = xj - xi;
                        rSq = dx.squaredNorm();
                        h = irad + radius[j];
                        hsq = h * h;
                        if (rSq < hsq) {

                                jtype = type[j];
                                Poly6Kernel(hsq, h, rSq, domain->dimension, wf);

                                if (setflag[itype][itype] == 1) {
                                        rho[i] += wf * rmass[j]; // / shepardWeight[i];
                                }

                                if (j < nlocal) {
                                        if (setflag[jtype][jtype] == 1) {
                                                rho[j] += wf * rmass[i]; // / shepardWeight[j];
                                        }
                                }
                        } // end if check distance
                } // end loop over j
        } // end loop over i
}

/* ----------------------------------------------------------------------
 *
 * Compute shape matrix for kernel gradient correction and velocity gradient.
 * This is used if material strength or viscosity models are employed.
 *
 ---------------------------------------------------------------------- */

void PairULSPH::PreCompute() {
        double **atom_data9 = atom->smd_data_9;
        double *radius = atom->radius;
        double **x = atom->x;
        double **x0 = atom->x0;
        double **v = atom->vest;
        double *vfrac = atom->vfrac;
        int *type = atom->type;
        int *ilist, *jlist, *numneigh;
        int **firstneigh;
        int nlocal = atom->nlocal;
        int inum, jnum, ii, jj, i, itype, j, idim;
        double wfd, h, irad, r, rSq, wf, ivol, jvol;
        Vector3d dx, dv, g, du;
        Matrix3d Ktmp, Ltmp, Ftmp, K3di, D;
        Vector3d xi, xj, vi, vj, x0i, x0j, dx0;
        Matrix2d K2di, K2d;

        // zero accumulators
        for (i = 0; i < nlocal; i++) {
                itype = type[i];
                if (setflag[itype][itype]) {
                        if (gradient_correction_flag) {
                                K[i].setZero();
                        } else {
                                K[i].setIdentity();
                        }
                        L[i].setZero();
                        F[i].setZero();
                }
        }

        // set up neighbor list variables
        inum = list->inum;
        ilist = list->ilist;
        numneigh = list->numneigh;
        firstneigh = list->firstneigh;

        for (ii = 0; ii < inum; ii++) {
                i = ilist[ii];
                itype = type[i];
                jlist = firstneigh[i];
                jnum = numneigh[i];
                irad = radius[i];
                ivol = vfrac[i];

                // initialize Eigen data structures from LAMMPS data structures
                for (idim = 0; idim < 3; idim++) {
                        x0i(idim) = x0[i][idim];
                        xi(idim) = x[i][idim];
                        vi(idim) = v[i][idim];
                }

                for (jj = 0; jj < jnum; jj++) {
                        j = jlist[jj];
                        j &= NEIGHMASK;

                        for (idim = 0; idim < 3; idim++) {
                                x0j(idim) = x0[j][idim];
                                xj(idim) = x[j][idim];
                                vj(idim) = v[j][idim];
                        }

                        dx = xj - xi;

                        rSq = dx.squaredNorm();
                        h = irad + radius[j];
                        if (rSq < h * h) {

                                r = sqrt(rSq);
                                jvol = vfrac[j];

                                // distance vectors in current and reference configuration, velocity difference
                                dv = vj - vi;
                                dx0 = x0j - x0i;

                                // kernel and derivative
                                spiky_kernel_and_derivative(h, r, domain->dimension, wf, wfd);
                                //barbara_kernel_and_derivative(h, r, domain->dimension, wf, wfd);

                                // uncorrected kernel gradient
                                g = (wfd / r) * dx;

                                /* build correction matrix for kernel derivatives */
                                if (gradient_correction_flag) {
                                        Ktmp = -g * dx.transpose();
                                        K[i] += jvol * Ktmp;
                                }

                                // velocity gradient L
                                Ltmp = -dv * g.transpose();
                                L[i] += jvol * Ltmp;

                                // deformation gradient F in Eulerian frame
                                du = dx - dx0;
                                Ftmp = dv * g.transpose();
                                F[i] += jvol * Ftmp;

                                if (j < nlocal) {

                                        if (gradient_correction_flag) {
                                                K[j] += ivol * Ktmp;
                                        }

                                        L[j] += ivol * Ltmp;
                                        F[j] += ivol * Ftmp;
                                }
                        } // end if check distance
                } // end loop over j
        } // end loop over i

        /*
         * invert shape matrix and compute corrected quantities
         */

        for (i = 0; i < nlocal; i++) {
                itype = type[i];
                if (setflag[itype][itype]) {
                        if (gradient_correction_flag) {
                                pseudo_inverse_SVD(K[i]);
                                K[i] = LimitEigenvalues(K[i], 2.0);
                                L[i] *= K[i];
                                F[i] *= K[i];
                        } // end if (gradient_correction[itype]) {

                        /*
                         * accumulate strain increments
                         * we abuse the atom array "atom_data_9" for this purpose, which was originally designed to hold the deformation gradient.
                         */
                        D = update->dt * 0.5 * (L[i] + L[i].transpose());
                        atom_data9[i][0] += D(0, 0); // xx
                        atom_data9[i][1] += D(1, 1); // yy
                        atom_data9[i][2] += D(2, 2); // zz
                        atom_data9[i][3] += D(0, 1); // xy
                        atom_data9[i][4] += D(0, 2); // xz
                        atom_data9[i][5] += D(1, 2); // yz

                } // end if (setflag[itype][itype])
        } // end loop over i = 0 to nlocal

}

/* ---------------------------------------------------------------------- */

void PairULSPH::compute(int eflag, int vflag) {
        double **x = atom->x;
        double **v = atom->vest;
        double **vint = atom->v; // Velocity-Verlet algorithm velocities
        double **f = atom->f;
        double *vfrac = atom->vfrac;
        double *desph = atom->desph;
        double *rmass = atom->rmass;
        double *radius = atom->radius;
        double *contact_radius = atom->contact_radius;
        double **atom_data9 = atom->smd_data_9;

        int *type = atom->type;
        int nlocal = atom->nlocal;
        int i, j, ii, jj, jnum, itype, jtype, iDim, inum;
        double r, wf, wfd, h, rSq, ivol, jvol;
        double mu_ij, c_ij, rho_ij;
        double delVdotDelR, visc_magnitude, deltaE;
        int *ilist, *jlist, *numneigh;
        int **firstneigh;
        Vector3d fi, fj, dx, dv, f_stress, g, vinti, vintj, dvint;
        Vector3d xi, xj, vi, vj, f_visc, sumForces, f_stress_new;
        Vector3d gamma, f_hg, dx0, du_est, du;
        double r_ref, weight, p;
        //int periodic = (domain->xperiodic || domain->yperiodic || domain->zperiodic);

        double ini_dist;
        Matrix3d S, D, V, eye;
        eye.setIdentity();
        int k;
        SelfAdjointEigenSolver < Matrix3d > es;

        ev_init(eflag, vflag);

        if (atom->nmax > nmax) {
//printf("... allocating in compute with nmax = %d\n", atom->nmax);
                nmax = atom->nmax;
                delete[] K;
                K = new Matrix3d[nmax];
                delete[] shepardWeight;
                shepardWeight = new double[nmax];
                delete[] c0;
                c0 = new double[nmax];
                delete[] smoothVel;
                smoothVel = new Vector3d[nmax];
                delete[] stressTensor;
                stressTensor = new Matrix3d[nmax];
                delete[] L;
                L = new Matrix3d[nmax];
                delete[] numNeighs;
                numNeighs = new int[nmax];
                delete[] F;
                F = new Matrix3d[nmax];
                delete[] rho;
                rho = new double[nmax];
                delete[] effm;
                effm = new double[nmax];
        }

// zero accumulators
        for (i = 0; i < nlocal; i++) {
                shepardWeight[i] = 0.0;
                smoothVel[i].setZero();
                numNeighs[i] = 0;

                h = 2.0 * radius[i];
                r = 0.0;
                spiky_kernel_and_derivative(h, r, domain->dimension, wf, wfd);
        }

        /*
         * if this is the very first step, zero the array which holds the accumulated strain
         */
        if (update->ntimestep == 0) {
                for (i = 0; i < nlocal; i++) {
                        itype = type[i];
                        if (setflag[itype][itype]) {
                                for (j = 0; j < 9; j++) {
                                        atom_data9[i][j] = 0.0;
                                }
                        }
                }
        }

        if (density_summation) {
                //printf("dens summ\n");
                PreCompute_DensitySummation();

                for (i = 0; i < nlocal; i++) { //compute volumes from rho
                        itype = type[i];
                        if (setflag[itype][itype]) {
                                vfrac[i] = rmass[i] / rho[i];
                        }
                }

        }

        if (velocity_gradient) {
                PairULSPH::PreCompute(); // get velocity gradient and kernel gradient correction
        }

        PairULSPH::AssembleStressTensor();

        /*
         * QUANTITIES ABOVE HAVE ONLY BEEN CALCULATED FOR NLOCAL PARTICLES.
         * NEED TO DO A FORWARD COMMUNICATION TO GHOST ATOMS NOW
         */
        comm->forward_comm(this);

        updateFlag = 0;

        /*
         * iterate over pairs of particles i, j and assign forces using pre-computed pressure
         */

// set up neighbor list variables
        inum = list->inum;
        ilist = list->ilist;
        numneigh = list->numneigh;
        firstneigh = list->firstneigh;

        for (ii = 0; ii < inum; ii++) {
                i = ilist[ii];
                itype = type[i];
                jlist = firstneigh[i];
                jnum = numneigh[i];
                ivol = vfrac[i];

                // initialize Eigen data structures from LAMMPS data structures
                for (iDim = 0; iDim < 3; iDim++) {
                        xi(iDim) = x[i][iDim];
                        vi(iDim) = v[i][iDim];
                        vinti(iDim) = vint[i][iDim];
                }

                for (jj = 0; jj < jnum; jj++) {
                        j = jlist[jj];
                        j &= NEIGHMASK;

                        xj(0) = x[j][0];
                        xj(1) = x[j][1];
                        xj(2) = x[j][2];

                        dx = xj - xi;
                        rSq = dx.squaredNorm();
                        h = radius[i] + radius[j];
                        if (rSq < h * h) {

                                // initialize Eigen data structures from LAMMPS data structures
                                for (iDim = 0; iDim < 3; iDim++) {
                                        vj(iDim) = v[j][iDim];
                                        vintj(iDim) = vint[j][iDim];
                                }

                                r = sqrt(rSq);
                                jtype = type[j];
                                jvol = vfrac[j];

                                // distance vectors in current and reference configuration, velocity difference
                                dv = vj - vi;
                                dvint = vintj - vinti;

                                // kernel and derivative
                                spiky_kernel_and_derivative(h, r, domain->dimension, wf, wfd);
                                //barbara_kernel_and_derivative(h, r, domain->dimension, wf, wfd);

                                // uncorrected kernel gradient
                                g = (wfd / r) * dx;

                                delVdotDelR = dx.dot(dv) / (r + 0.1 * h); // project relative velocity onto unit particle distance vector [m/s]

                                S = stressTensor[i] + stressTensor[j];

                                if (artificial_pressure[itype][jtype] > 0.0) {
                                        p = S.trace();
                                        if (p > 0.0) { // we are in tension
                                                r_ref = contact_radius[i] + contact_radius[j];
                                                weight = Kernel_Cubic_Spline(r, h) / Kernel_Cubic_Spline(r_ref, h);
                                                weight = pow(weight, 4.0);
                                                S -= artificial_pressure[itype][jtype] * weight * p * eye;
                                        }
                                }

                                /*
                                 * artificial stress to control tensile instability
                                 * Only works if particles are uniformly spaced initially.
                                 */
                                if (artificial_stress[itype][jtype] > 0.0) {
                                        ini_dist = contact_radius[i] + contact_radius[j];
                                        weight = Kernel_Cubic_Spline(r, h) / Kernel_Cubic_Spline(ini_dist, h);
                                        weight = pow(weight, 4.0);

                                        es.compute(S);
                                        D = es.eigenvalues().asDiagonal();
                                        for (k = 0; k < 3; k++) {
                                                if (D(k, k) > 0.0) {
                                                        D(k, k) -= weight * artificial_stress[itype][jtype] * D(k, k);
                                                }
                                        }
                                        V = es.eigenvectors();
                                        S = V * D * V.inverse();
                                }

                                // compute forces
                                f_stress = -ivol * jvol * S * g; // DO NOT TOUCH SIGN

                                /*
                                 * artificial viscosity -- alpha is dimensionless
                                 * MonaghanBalsara form of the artificial viscosity
                                 */

                                c_ij = 0.5 * (c0[i] + c0[j]);
                                LimitDoubleMagnitude(delVdotDelR, 1.1 * c_ij);

                                mu_ij = h * delVdotDelR / (r + 0.1 * h); // units: [m * m/s / m = m/s]
                                rho_ij = 0.5 * (rmass[i] / ivol + rmass[j] / jvol);
                                visc_magnitude = 0.5 * (Q1[itype] + Q1[jtype]) * c_ij * mu_ij / rho_ij;
                                f_visc = -rmass[i] * rmass[j] * visc_magnitude * g;

                                if ((Lookup[HOURGLASS_CONTROL_AMPLITUDE][itype] > 0.0) && (Lookup[HOURGLASS_CONTROL_AMPLITUDE][jtype] > 0.0)) {
                                        f_hg = ComputeHourglassForce(i, itype, j, jtype, dv, dx, g, c_ij, mu_ij, rho_ij);

                                } else {
                                        f_hg.setZero();
                                }

                                sumForces = f_stress + f_visc + f_hg;

                                // energy rate -- project velocity onto force vector
                                deltaE = sumForces.dot(dv);

                                // apply forces to pair of particles
                                f[i][0] += sumForces(0);
                                f[i][1] += sumForces(1);
                                f[i][2] += sumForces(2);
                                desph[i] += deltaE;

                                // accumulate smooth velocities
                                shepardWeight[i] += jvol * wf;
                                smoothVel[i] += jvol * wf * dvint;
                                numNeighs[i] += 1;

                                if (j < nlocal) {
                                        f[j][0] -= sumForces(0);
                                        f[j][1] -= sumForces(1);
                                        f[j][2] -= sumForces(2);
                                        desph[j] += deltaE;

                                        shepardWeight[j] += ivol * wf;
                                        smoothVel[j] -= ivol * wf * dvint;
                                        numNeighs[j] += 1;
                                }

                                // tally atomistic stress tensor
                                if (evflag) {
                                        ev_tally_xyz(i, j, nlocal, 0, 0.0, 0.0, sumForces(0), sumForces(1), sumForces(2), dx(0), dx(1), dx(2));
                                }
                        }

                }
        }

        for (i = 0; i < nlocal; i++) {
                itype = type[i];
                if (setflag[itype][itype] == 1) {
                        if (shepardWeight[i] != 0.0) {
                                smoothVel[i] /= shepardWeight[i];
                        } else {
                                smoothVel[i].setZero();
                        }
                } // end check if particle is SPH-type
        } // end loop over i = 0 to nlocal

        if (vflag_fdotr)
                virial_fdotr_compute();

}

/* ----------------------------------------------------------------------
 Assemble total stress tensor with pressure, material sterength, and
 viscosity contributions.
 ------------------------------------------------------------------------- */
void PairULSPH::AssembleStressTensor() {
        double *radius = atom->radius;
        double *vfrac = atom->vfrac;
        double *rmass = atom->rmass;
        double *eff_plastic_strain = atom->eff_plastic_strain;
        double **tlsph_stress = atom->smd_stress;
        double *esph = atom->esph;
        int *type = atom->type;
        int i, itype;
        int nlocal = atom->nlocal;
        Matrix3d D, Ddev, W, V, sigma_diag;
        Matrix3d eye, stressRate, StressRateDevJaumann;
        Matrix3d sigmaInitial_dev, d_dev, sigmaFinal_dev, stressRateDev, oldStressDeviator, newStressDeviator;
        double plastic_strain_increment, yieldStress;
        double dt = update->dt;
        double vol, newPressure;
        double G_eff = 0.0; // effective shear modulus
        double K_eff; // effective bulk modulus
        double M, p_wave_speed;
        double rho, effectiveViscosity;
        Matrix3d deltaStressDev;

        dtCFL = 1.0e22;
        eye.setIdentity();

        for (i = 0; i < nlocal; i++) {
                itype = type[i];
                if (setflag[itype][itype] == 1) {
                        newStressDeviator.setZero();
                        newPressure = 0.0;
                        stressTensor[i].setZero();
                        vol = vfrac[i];
                        rho = rmass[i] / vfrac[i];
                        effectiveViscosity = 0.0;
                        K_eff = 0.0;
                        G_eff = 0.0;

                        //printf("rho = %f\n", rho);

                        switch (eos[itype]) {
                        default:
                                error->one(FLERR, "unknown EOS.");
                                break;
                        case NONE:
                                c0[i] = 1.0;
                                break;
                        case EOS_TAIT:
                                TaitEOS_density(Lookup[EOS_TAIT_EXPONENT][itype], Lookup[REFERENCE_SOUNDSPEED][itype],
                                                Lookup[REFERENCE_DENSITY][itype], rho, newPressure, c0[i]);
                                //printf("new pressure =%f\n", newPressure);

                                break;
                        case EOS_PERFECT_GAS:
                                PerfectGasEOS(Lookup[EOS_PERFECT_GAS_GAMMA][itype], vol, rmass[i], esph[i], newPressure, c0[i]);
                                break;
                        case EOS_LINEAR:
                                newPressure = Lookup[BULK_MODULUS][itype] * (rho / Lookup[REFERENCE_DENSITY][itype] - 1.0);
                                //printf("p=%f, rho0=%f, rho=%f\n", newPressure, Lookup[REFERENCE_DENSITY][itype], rho);
                                c0[i] = Lookup[REFERENCE_SOUNDSPEED][itype];
                                break;
                        }

                        K_eff = c0[i] * c0[i] * rho; // effective bulk modulus

                        /*
                         * ******************************* STRENGTH MODELS ************************************************
                         */

                        if (strength[itype] != NONE) {
                                /*
                                 * initial stress state: given by the unrotateted Cauchy stress.
                                 * Assemble Eigen 3d matrix from stored stress state
                                 */
                                oldStressDeviator(0, 0) = tlsph_stress[i][0];
                                oldStressDeviator(0, 1) = tlsph_stress[i][1];
                                oldStressDeviator(0, 2) = tlsph_stress[i][2];
                                oldStressDeviator(1, 1) = tlsph_stress[i][3];
                                oldStressDeviator(1, 2) = tlsph_stress[i][4];
                                oldStressDeviator(2, 2) = tlsph_stress[i][5];
                                oldStressDeviator(1, 0) = oldStressDeviator(0, 1);
                                oldStressDeviator(2, 0) = oldStressDeviator(0, 2);
                                oldStressDeviator(2, 1) = oldStressDeviator(1, 2);

                                D = 0.5 * (L[i] + L[i].transpose());
                                W = 0.5 * (L[i] - L[i].transpose()); // spin tensor:: need this for Jaumann rate
                                d_dev = Deviator(D);

                                switch (strength[itype]) {
                                default:
                                        error->one(FLERR, "unknown strength model.");
                                        break;
                                case STRENGTH_LINEAR:

                                        // here in a version with pressure part
//                                      stressRateDev = Lookup[BULK_MODULUS][itype] * d_iso * eye + 2.0 * Lookup[SHEAR_MODULUS][itype] * d_dev;
//                                      c0[i] = Lookup[REFERENCE_SOUNDSPEED][itype];
//                                      newPressure = 0.0;

                                        // here only stress deviator
                                        stressRateDev = 2.0 * Lookup[SHEAR_MODULUS][itype] * d_dev;
                                        //cout << "stress rate deviator is " << endl << stressRateDev << endl;
                                        break;

                                case STRENGTH_LINEAR_PLASTIC:
                                        yieldStress = Lookup[YIELD_STRENGTH][itype] + Lookup[HARDENING_PARAMETER][itype] * eff_plastic_strain[i];
                                        LinearPlasticStrength(Lookup[SHEAR_MODULUS][itype], yieldStress, oldStressDeviator, d_dev, dt,
                                                        newStressDeviator, stressRateDev, plastic_strain_increment);
                                        eff_plastic_strain[i] += plastic_strain_increment;

                                        break;
                                }

                                //double m = effective_longitudinal_modulus(itype, dt, d_iso, p_rate, d_dev, stressRate_dev, damage);

                                StressRateDevJaumann = stressRateDev - W * oldStressDeviator + oldStressDeviator * W;
                                newStressDeviator = oldStressDeviator + dt * StressRateDevJaumann;

                                tlsph_stress[i][0] = newStressDeviator(0, 0);
                                tlsph_stress[i][1] = newStressDeviator(0, 1);
                                tlsph_stress[i][2] = newStressDeviator(0, 2);
                                tlsph_stress[i][3] = newStressDeviator(1, 1);
                                tlsph_stress[i][4] = newStressDeviator(1, 2);
                                tlsph_stress[i][5] = newStressDeviator(2, 2);

                                // estimate effective shear modulus for time step stability
                                deltaStressDev = oldStressDeviator - newStressDeviator;
                                G_eff = effective_shear_modulus(d_dev, deltaStressDev, dt, itype);

                        } // end if (strength[itype] != NONE)

                        if (viscosity[itype] != NONE) {
                                D = 0.5 * (L[i] + L[i].transpose());
                                d_dev = Deviator(D);

                                switch (viscosity[itype]) {
                                default:
                                        error->one(FLERR, "unknown viscosity model.");
                                        break;
                                case VISCOSITY_NEWTON:
                                        effectiveViscosity = Lookup[VISCOSITY_MU][itype];
//                                      double shear_rate = 2.0
//                                                      * sqrt(d_dev(0, 1) * d_dev(0, 1) + d_dev(0, 2) * d_dev(0, 2) + d_dev(1, 2) * d_dev(1, 2)); // 3d
                                        //cout << "shear rate: " << shear_rate << endl;
                                        //effectiveViscosity = PA6_270C(shear_rate);
                                        //if (effectiveViscosity > 178.062e-6) {
                                        //      printf("effective visc is %f\n", effectiveViscosity);
                                        //}
                                        newStressDeviator = 2.0 * effectiveViscosity * d_dev; // newton original
                                        //cout << "this is Ddev " << endl << d_dev << endl << endl;
                                        break;
                                }
                        } // end if (viscosity[itype] != NONE)

                        /*
                         * assemble stress Tensor from pressure and deviatoric parts
                         */

                        stressTensor[i] = -newPressure * eye + newStressDeviator;

                        /*
                         * stable timestep based on speed-of-sound
                         */

                        M = K_eff + 4.0 * G_eff / 3.0;
                        p_wave_speed = sqrt(M / rho);
                        effm[i] = G_eff;
                        dtCFL = MIN(2 * radius[i] / p_wave_speed, dtCFL);

                        /*
                         * stable timestep based on viscosity
                         */
                        if (viscosity[itype] != NONE) {
                                dtCFL = MIN(4 * radius[i] * radius[i] * rho / effectiveViscosity, dtCFL);
                        }

                        /*
                         * kernel gradient correction
                         */
                        if (gradient_correction_flag) {
                                stressTensor[i] *= K[i];
                        }
                }
                // end if (setflag[itype][itype] == 1)
        } // end loop over nlocal

//printf("stable timestep = %g\n", 0.1 * hMin * MaxBulkVelocity);
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairULSPH::allocate() {

        allocated = 1;
        int n = atom->ntypes;

        memory->create(setflag, n + 1, n + 1, "pair:setflag");

        memory->create(Q1, n + 1, "pair:Q1");
        memory->create(rho0, n + 1, "pair:Q2");
        memory->create(c0_type, n + 1, "pair:c0_type");
        memory->create(eos, n + 1, "pair:eosmodel");
        memory->create(viscosity, n + 1, "pair:viscositymodel");
        memory->create(strength, n + 1, "pair:strengthmodel");

        memory->create(Lookup, MAX_KEY_VALUE, n + 1, "pair:LookupTable");

        memory->create(artificial_pressure, n + 1, n + 1, "pair:artificial_pressure");
        memory->create(artificial_stress, n + 1, n + 1, "pair:artificial_stress");
        memory->create(cutsq, n + 1, n + 1, "pair:cutsq");              // always needs to be allocated, even with granular neighborlist

        /*
         * initialize arrays to default values
         */

        for (int i = 1; i <= n; i++) {
                for (int j = i; j <= n; j++) {
                        artificial_pressure[i][j] = 0.0;
                        artificial_stress[i][j] = 0.0;
                        setflag[i][j] = 0;
                }
        }

        onerad_dynamic = new double[n + 1];
        onerad_frozen = new double[n + 1];
        maxrad_dynamic = new double[n + 1];
        maxrad_frozen = new double[n + 1];

}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairULSPH::settings(int narg, char **arg) {
        if (narg != 3) {
                printf("narg = %d\n", narg);
                error->all(FLERR, "Illegal number of arguments for pair_style ulsph");
        }

        if (comm->me == 0) {
                printf("\n>>========>>========>>========>>========>>========>>========>>========>>========\n");
                printf("... SMD / ULSPH PROPERTIES\n\n");
        }

        if (strcmp(arg[0], "*DENSITY_SUMMATION") == 0) {
                density_summation = true;
                density_continuity = false;
                if (comm->me == 0)
                        printf("... density summation active\n");
        } else if (strcmp(arg[0], "*DENSITY_CONTINUITY") == 0) {
                density_continuity = true;
                density_summation = false;
                if (comm->me == 0)
                        printf("... density continuity active\n");
        } else {
                error->all(FLERR,
                                "Illegal settings keyword for first keyword of pair style ulsph. Must be either *DENSITY_SUMMATION or *DENSITY_CONTINUITY");
        }

        if (strcmp(arg[1], "*VELOCITY_GRADIENT") == 0) {
                velocity_gradient = true;
                if (comm->me == 0)
                        printf("... computation of velocity gradients active\n");
        } else if (strcmp(arg[1], "*NO_VELOCITY_GRADIENT") == 0) {
                velocity_gradient = false;
                if (comm->me == 0)
                        printf("... computation of velocity gradients NOT active\n");
        } else {
                error->all(FLERR,
                                "Illegal settings keyword for first keyword of pair style ulsph. Must be either *VELOCITY_GRADIENT or *NO_VELOCITY_GRADIENT");
        }

        if (strcmp(arg[2], "*GRADIENT_CORRECTION") == 0) {
                gradient_correction_flag = true;
                if (comm->me == 0)
                        printf("... first order correction of kernel gradients is active\n");
        } else if (strcmp(arg[2], "*NO_GRADIENT_CORRECTION") == 0) {
                gradient_correction_flag = false;
                if (comm->me == 0)
                        printf("... first order correction of kernel gradients is NOT active\n");
        } else {
                error->all(FLERR, "Illegal settings keyword for pair style ulsph");
        }

// error check
        //if ((gradient_correction_flag == true) && (density_summation)) {
        //      error->all(FLERR, "Cannot use *DENSITY_SUMMATION in combination with *YES_GRADIENT_CORRECTION");
        //}

        if (comm->me == 0)
                printf(">>========>>========>>========>>========>>========>>========>>========>>========\n");

}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairULSPH::coeff(int narg, char **arg) {
        int ioffset, iarg, iNextKwd, itype, jtype;
        std::string s, t;

        if (narg < 3) utils::missing_cmd_args(FLERR, "pair ulsph", error);

        if (!allocated) allocate();

        /*
         * if parameters are give in i,i form, i.e., no a cross interaction, set material parameters
         */

        if (utils::inumeric(FLERR, arg[0], false, lmp) == utils::inumeric(FLERR, arg[1], false, lmp)) {

                itype = utils::inumeric(FLERR, arg[0],false,lmp);
                eos[itype] = viscosity[itype] = strength[itype] = NONE;

                if (comm->me == 0) {
                        printf("\n>>========>>========>>========>>========>>========>>========>>========>>========\n");
                        printf("...SMD / ULSPH PROPERTIES OF PARTICLE TYPE %d\n\n", itype);
                }

                /*
                 * read parameters which are common -- regardless of material / eos model
                 */

                ioffset = 2;
                if (strcmp(arg[ioffset], "*COMMON") != 0) error->all(FLERR, "common keyword missing!");

                t = string("*");
                iNextKwd = -1;
                for (iarg = ioffset + 1; iarg < narg; iarg++) {
                        s = string(arg[iarg]);
                        if (s.compare(0, t.length(), t) == 0) {
                                iNextKwd = iarg;
                                break;
                        }
                }

                if (iNextKwd < 0) error->all(FLERR, "no *KEYWORD terminates *COMMON");

                if (iNextKwd - ioffset != 5 + 1)
                  error->all(FLERR, "expected 5 arguments following *COMMON but got {}\n", iNextKwd - ioffset - 1);

                Lookup[REFERENCE_DENSITY][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);
                Lookup[REFERENCE_SOUNDSPEED][itype] = utils::numeric(FLERR, arg[ioffset + 2],false,lmp);
                Q1[itype] = utils::numeric(FLERR, arg[ioffset + 3],false,lmp);
                Lookup[HEAT_CAPACITY][itype] = utils::numeric(FLERR, arg[ioffset + 4],false,lmp);
                Lookup[HOURGLASS_CONTROL_AMPLITUDE][itype] = utils::numeric(FLERR, arg[ioffset + 5],false,lmp);

                Lookup[BULK_MODULUS][itype] = Lookup[REFERENCE_SOUNDSPEED][itype] * Lookup[REFERENCE_SOUNDSPEED][itype]
                                * Lookup[REFERENCE_DENSITY][itype];

                if (comm->me == 0) {
                        printf("material unspecific properties for SMD/ULSPH definition of particle type %d:\n", itype);
                        printf(FORMAT1, "reference density", Lookup[REFERENCE_DENSITY][itype]);
                        printf(FORMAT1, "reference speed of sound", Lookup[REFERENCE_SOUNDSPEED][itype]);
                        printf(FORMAT1, "linear viscosity coefficient", Q1[itype]);
                        printf(FORMAT1, "heat capacity [energy / (mass * temperature)]", Lookup[HEAT_CAPACITY][itype]);
                        printf(FORMAT1, "bulk modulus", Lookup[BULK_MODULUS][itype]);
                        printf(FORMAT1, "hourglass control amplitude", Lookup[HOURGLASS_CONTROL_AMPLITUDE][itype]);
                }

                /*
                 * read following material cards
                 */

                while (true) {
                        if (strcmp(arg[iNextKwd], "*END") == 0) {
                                break;
                        }

                        ioffset = iNextKwd;
                        if (strcmp(arg[ioffset], "*EOS_TAIT") == 0) {

                                /*
                                 * Tait EOS
                                 */

                                eos[itype] = EOS_TAIT;

                                t = string("*");
                                iNextKwd = -1;
                                for (iarg = ioffset + 1; iarg < narg; iarg++) {
                                        s = string(arg[iarg]);
                                        if (s.compare(0, t.length(), t) == 0) {
                                                iNextKwd = iarg;
                                                break;
                                        }
                                }

                                if (iNextKwd < 0) error->all(FLERR, "no *KEYWORD terminates *EOS_TAIT");

                                if (iNextKwd - ioffset != 1 + 1)
                                        error->all(FLERR, "expected 1 arguments following *EOS_TAIT but got {}\n", iNextKwd - ioffset - 1);

                                Lookup[EOS_TAIT_EXPONENT][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);

                                if (comm->me == 0) {
                                        printf(FORMAT2, "Tait EOS");
                                        printf(FORMAT1, "Exponent", Lookup[EOS_TAIT_EXPONENT][itype]);
                                }
                        } // end Tait EOS

                        else if (strcmp(arg[ioffset], "*EOS_PERFECT_GAS") == 0) {

                                /*
                                 * Perfect Gas EOS
                                 */

                                eos[itype] = EOS_PERFECT_GAS;

                                t = string("*");
                                iNextKwd = -1;
                                for (iarg = ioffset + 1; iarg < narg; iarg++) {
                                        s = string(arg[iarg]);
                                        if (s.compare(0, t.length(), t) == 0) {
                                                iNextKwd = iarg;
                                                break;
                                        }
                                }

                                if (iNextKwd < 0) error->all(FLERR, "no *KEYWORD terminates *EOS_PERFECT_GAS");
                                if (iNextKwd - ioffset != 1 + 1)
                                  error->all(FLERR, "expected 1 arguments following *EOS_PERFECT_GAS but got {}\n", iNextKwd - ioffset - 1);

                                Lookup[EOS_PERFECT_GAS_GAMMA][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);

                                if (comm->me == 0) {
                                        printf(FORMAT2, "Perfect Gas EOS");
                                        printf(FORMAT1, "Heat Capacity Ratio Gamma", Lookup[EOS_PERFECT_GAS_GAMMA][itype]);
                                }
                        } // end Perfect Gas EOS
                        else if (strcmp(arg[ioffset], "*EOS_LINEAR") == 0) {

                                /*
                                 * Linear EOS
                                 */

                                eos[itype] = EOS_LINEAR;

                                t = string("*");
                                iNextKwd = -1;
                                for (iarg = ioffset + 1; iarg < narg; iarg++) {
                                        s = string(arg[iarg]);
                                        if (s.compare(0, t.length(), t) == 0) {
                                                iNextKwd = iarg;
                                                break;
                                        }
                                }

                                if (iNextKwd < 0) error->all(FLERR, "no *KEYWORD terminates *EOS_LINEAR");
                                if (iNextKwd - ioffset != 0 + 1)
                                        error->all(FLERR, "expected 0 arguments following *EOS_LINEAR but got {}\n", iNextKwd - ioffset - 1);

                                if (comm->me == 0) {
                                        printf(FORMAT2, "Linear EOS");
                                        printf(FORMAT1, "Bulk modulus", Lookup[BULK_MODULUS][itype]);
                                }
                        } // end Linear EOS
                        else if (strcmp(arg[ioffset], "*STRENGTH_LINEAR_PLASTIC") == 0) {

                                if (!velocity_gradient) {
                                        error->all(FLERR, "A strength model was requested but *VELOCITY_GRADIENT is not set");
                                }

                                /*
                                 * linear elastic / ideal plastic material model with strength
                                 */

                                strength[itype] = STRENGTH_LINEAR_PLASTIC;
                                velocity_gradient_required = true;

                                t = string("*");
                                iNextKwd = -1;
                                for (iarg = ioffset + 1; iarg < narg; iarg++) {
                                        s = string(arg[iarg]);
                                        if (s.compare(0, t.length(), t) == 0) {
                                                iNextKwd = iarg;
                                                break;
                                        }
                                }

                                if (iNextKwd < 0) error->all(FLERR, "no *KEYWORD terminates *STRENGTH_LINEAR_PLASTIC");
                                if (iNextKwd - ioffset != 3 + 1)
                                        error->all(FLERR, "expected 3 arguments following *STRENGTH_LINEAR_PLASTIC but got {}\n", iNextKwd - ioffset - 1);

                                Lookup[SHEAR_MODULUS][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);
                                Lookup[YIELD_STRENGTH][itype] = utils::numeric(FLERR, arg[ioffset + 2],false,lmp);
                                Lookup[HARDENING_PARAMETER][itype] = utils::numeric(FLERR, arg[ioffset + 3],false,lmp);

                                if (comm->me == 0) {
                                        printf(FORMAT2, "linear elastic / ideal plastic material mode");
                                        printf(FORMAT1, "yield_strength", Lookup[YIELD_STRENGTH][itype]);
                                        printf(FORMAT1, "constant hardening parameter", Lookup[HARDENING_PARAMETER][itype]);
                                        printf(FORMAT1, "shear modulus", Lookup[SHEAR_MODULUS][itype]);
                                }
                        } // end *STRENGTH_LINEAR_PLASTIC
                        else if (strcmp(arg[ioffset], "*STRENGTH_LINEAR") == 0) {

                                if (!velocity_gradient) {
                                        error->all(FLERR, "A strength model was requested but *VELOCITY_GRADIENT is not set");
                                }

                                /*
                                 * linear elastic / ideal plastic material model with strength
                                 */

                                strength[itype] = STRENGTH_LINEAR;
                                t = string("*");
                                iNextKwd = -1;
                                for (iarg = ioffset + 1; iarg < narg; iarg++) {
                                        s = string(arg[iarg]);
                                        if (s.compare(0, t.length(), t) == 0) {
                                                iNextKwd = iarg;
                                                break;
                                        }
                                }

                                if (iNextKwd < 0) error->all(FLERR, "no *KEYWORD terminates *STRENGTH_LINEAR");
                                if (iNextKwd - ioffset != 1 + 1)
                                        error->all(FLERR, "expected 1 arguments following *STRENGTH_LINEAR but got {}\n", iNextKwd - ioffset - 1);

                                Lookup[SHEAR_MODULUS][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);

                                if (comm->me == 0) {
                                        printf(FORMAT2, "linear elastic strength model");
                                        printf(FORMAT1, "shear modulus", Lookup[SHEAR_MODULUS][itype]);
                                }
                        } // end *STRENGTH_LINEAR
                        else if (strcmp(arg[ioffset], "*VISCOSITY_NEWTON") == 0) {

                                if (!velocity_gradient) {
                                        error->all(FLERR, "A viscosity model was requested but *VELOCITY_GRADIENT is not set");
                                }

                                /*
                                 * linear elastic / ideal plastic material model with strength
                                 */

                                viscosity[itype] = VISCOSITY_NEWTON;
                                t = string("*");
                                iNextKwd = -1;
                                for (iarg = ioffset + 1; iarg < narg; iarg++) {
                                        s = string(arg[iarg]);
                                        if (s.compare(0, t.length(), t) == 0) {
                                                iNextKwd = iarg;
                                                break;
                                        }
                                }

                                if (iNextKwd < 0) error->all(FLERR, "no *KEYWORD terminates *VISCOSITY_NEWTON");
                                if (iNextKwd - ioffset != 1 + 1)
                                        error->all(FLERR, "expected 1 arguments following *VISCOSITY_NEWTON but got {}\n", iNextKwd - ioffset - 1);

                                Lookup[VISCOSITY_MU][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);

                                if (comm->me == 0) {
                                        printf(FORMAT2, "Newton viscosity model");
                                        printf(FORMAT1, "viscosity mu", Lookup[VISCOSITY_MU][itype]);
                                }
                        } // end *STRENGTH_VISCOSITY_NEWTON

                        else if (strcmp(arg[ioffset], "*ARTIFICIAL_PRESSURE") == 0) {

                                /*
                                 * use Monaghan's artificial pressure to prevent particle clumping
                                 */

                                t = string("*");
                                iNextKwd = -1;
                                for (iarg = ioffset + 1; iarg < narg; iarg++) {
                                        s = string(arg[iarg]);
                                        if (s.compare(0, t.length(), t) == 0) {
                                                iNextKwd = iarg;
                                                break;
                                        }
                                }

                                if (iNextKwd < 0) error->all(FLERR, "no *KEYWORD terminates *ARTIFICIAL_PRESSURE");
                                if (iNextKwd - ioffset != 1 + 1)
                                        error->all(FLERR, "expected 1 arguments following *ARTIFICIAL_PRESSURE but got {}\n", iNextKwd - ioffset - 1);

                                artificial_pressure[itype][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);

                                if (comm->me == 0) {
                                        printf(FORMAT2, "Artificial Pressure is enabled.");
                                        printf(FORMAT1, "Artificial Pressure amplitude", artificial_pressure[itype][itype]);
                                }
                        } // end *ARTIFICIAL_PRESSURE

                        else if (strcmp(arg[ioffset], "*ARTIFICIAL_STRESS") == 0) {

                                /*
                                 * use Monaghan's artificial stress to prevent particle clumping
                                 */

                                t = string("*");
                                iNextKwd = -1;
                                for (iarg = ioffset + 1; iarg < narg; iarg++) {
                                        s = string(arg[iarg]);
                                        if (s.compare(0, t.length(), t) == 0) {
                                                iNextKwd = iarg;
                                                break;
                                        }
                                }

                                if (iNextKwd < 0) error->all(FLERR, "no *KEYWORD terminates *ARTIFICIAL_STRESS");
                                if (iNextKwd - ioffset != 1 + 1)
                                        error->all(FLERR, "expected 1 arguments following *ARTIFICIAL_STRESS but got {}\n", iNextKwd - ioffset - 1);

                                artificial_stress[itype][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);

                                if (comm->me == 0) {
                                        printf(FORMAT2, "Artificial Stress is enabled.");
                                        printf(FORMAT1, "Artificial Stress amplitude", artificial_stress[itype][itype]);
                                }
                        } // end *ARTIFICIAL_STRESS

                        else error->all(FLERR, "unknown *KEYWORD: {}", arg[ioffset]);
                }

                /*
                 * copy data which is looked up in inner pairwise loops from slow maps to fast arrays
                 */

                rho0[itype] = Lookup[REFERENCE_DENSITY][itype];
                c0_type[itype] = Lookup[REFERENCE_SOUNDSPEED][itype];
                setflag[itype][itype] = 1;

                /*
                 * error checks
                 */

                if ((viscosity[itype] != NONE) && (strength[itype] != NONE))
                  error->all(FLERR, "cannot have both a strength and viscosity model for particle type {}", itype);

                if (eos[itype] == NONE)
                  error->all(FLERR, "must specify an EOS for particle type {}", itype);

        } else {
                /*
                 * we are reading a cross-interaction line for particle types i, j
                 */

                itype = utils::inumeric(FLERR, arg[0],false,lmp);
                jtype = utils::inumeric(FLERR, arg[1],false,lmp);

                if (strcmp(arg[2], "*CROSS") != 0)
                        error->all(FLERR, "ulsph cross interaction between particle type {} and {} requested, however, *CROSS keyword is missing",
                                        itype, jtype);

                if (setflag[itype][itype] != 1)
                        error->all(FLERR, "ulsph cross interaction between particle type {} and {} requested, however, properties of type {}  have not yet been specified", itype, jtype, itype);

                if (setflag[jtype][jtype] != 1)
                        error->all(FLERR, "ulsph cross interaction between particle type {} and {} requested, however, properties of type {}  have not yet been specified", itype, jtype, jtype);

                setflag[itype][jtype] = 1;
                setflag[jtype][itype] = 1;

                if ((artificial_pressure[itype][itype] > 0.0) && (artificial_pressure[jtype][jtype] > 0.0)) {
                        artificial_pressure[itype][jtype] = 0.5 * (artificial_pressure[itype][itype] + artificial_pressure[jtype][jtype]);
                        artificial_pressure[jtype][itype] = artificial_pressure[itype][jtype];
                } else {
                        artificial_pressure[itype][jtype] = artificial_pressure[jtype][itype] = 0.0;
                }

                if ((artificial_stress[itype][itype] > 0.0) && (artificial_stress[jtype][jtype] > 0.0)) {
                        artificial_stress[itype][jtype] = 0.5 * (artificial_stress[itype][itype] + artificial_stress[jtype][jtype]);
                        artificial_stress[jtype][itype] = artificial_stress[itype][jtype];
                } else {
                        artificial_stress[itype][jtype] = artificial_stress[jtype][itype] = 0.0;
                }

                if (comm->me == 0) {
                        printf(">>========>>========>>========>>========>>========>>========>>========>>========\n");
                }

        }
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairULSPH::init_one(int i, int j) {

        if (!allocated) allocate();

        if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

// cutoff = sum of max I,J radii for
// dynamic/dynamic & dynamic/frozen interactions, but not frozen/frozen

        double cutoff = maxrad_dynamic[i] + maxrad_dynamic[j];
        cutoff = MAX(cutoff, maxrad_frozen[i] + maxrad_dynamic[j]);
        cutoff = MAX(cutoff, maxrad_dynamic[i] + maxrad_frozen[j]);
        return cutoff;
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairULSPH::init_style() {
        int i;

 // request a granular neighbor list
        neighbor->add_request(this, NeighConst::REQ_SIZE);

// set maxrad_dynamic and maxrad_frozen for each type
// include future Fix pour particles as dynamic

        for (i = 1; i <= atom->ntypes; i++)
                onerad_dynamic[i] = onerad_frozen[i] = 0.0;

        double *radius = atom->radius;
        int *type = atom->type;
        int nlocal = atom->nlocal;

        for (i = 0; i < nlocal; i++)
                onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]], radius[i]);

        MPI_Allreduce(&onerad_dynamic[1], &maxrad_dynamic[1], atom->ntypes, MPI_DOUBLE, MPI_MAX, world);
        MPI_Allreduce(&onerad_frozen[1], &maxrad_frozen[1], atom->ntypes, MPI_DOUBLE, MPI_MAX, world);

}

/* ----------------------------------------------------------------------
 neighbor callback to inform pair style of neighbor list to use
 optional granular history list
 ------------------------------------------------------------------------- */

void PairULSPH::init_list(int id, NeighList *ptr) {
        if (id == 0)
                list = ptr;
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based arrays
 ------------------------------------------------------------------------- */

double PairULSPH::memory_usage() {
        return 11 * nmax * sizeof(double);
}

/* ---------------------------------------------------------------------- */

int PairULSPH::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/) {
        double *vfrac = atom->vfrac;
        double *eff_plastic_strain = atom->eff_plastic_strain;
        int i, j, m;

//printf("packing comm\n");
        m = 0;
        for (i = 0; i < n; i++) {
                j = list[i];
                buf[m++] = vfrac[j];
                buf[m++] = c0[j]; //2

                buf[m++] = stressTensor[j](0, 0); // pack symmetric stress tensor
                buf[m++] = stressTensor[j](1, 1);
                buf[m++] = stressTensor[j](2, 2);
                buf[m++] = stressTensor[j](0, 1);
                buf[m++] = stressTensor[j](0, 2);
                buf[m++] = stressTensor[j](1, 2); // 2 + 6 = 8

                buf[m++] = F[j](0, 0); // F is not symmetric
                buf[m++] = F[j](0, 1);
                buf[m++] = F[j](0, 2);
                buf[m++] = F[j](1, 0);
                buf[m++] = F[j](1, 1);
                buf[m++] = F[j](1, 2);
                buf[m++] = F[j](2, 0);
                buf[m++] = F[j](2, 1);
                buf[m++] = F[j](2, 2); // 8 + 9 = 17

                buf[m++] = eff_plastic_strain[j]; // 18
        }
        return m;
}

/* ---------------------------------------------------------------------- */

void PairULSPH::unpack_forward_comm(int n, int first, double *buf) {
        double *vfrac = atom->vfrac;
        double *eff_plastic_strain = atom->eff_plastic_strain;
        int i, m, last;

        m = 0;
        last = first + n;
        for (i = first; i < last; i++) {
                vfrac[i] = buf[m++];
                c0[i] = buf[m++]; // 2

                stressTensor[i](0, 0) = buf[m++];
                stressTensor[i](1, 1) = buf[m++];
                stressTensor[i](2, 2) = buf[m++];
                stressTensor[i](0, 1) = buf[m++];
                stressTensor[i](0, 2) = buf[m++];
                stressTensor[i](1, 2) = buf[m++]; // 2 + 6 = 8
                stressTensor[i](1, 0) = stressTensor[i](0, 1);
                stressTensor[i](2, 0) = stressTensor[i](0, 2);
                stressTensor[i](2, 1) = stressTensor[i](1, 2);

                F[i](0, 0) = buf[m++];
                F[i](0, 1) = buf[m++];
                F[i](0, 2) = buf[m++];
                F[i](1, 0) = buf[m++];
                F[i](1, 1) = buf[m++];
                F[i](1, 2) = buf[m++];
                F[i](2, 0) = buf[m++];
                F[i](2, 1) = buf[m++];
                F[i](2, 2) = buf[m++]; // 8 + 9 = 17

                eff_plastic_strain[i] = buf[m++]; // 18
        }
}

/*
 * EXTRACT
 */

void *PairULSPH::extract(const char *str, int &/*i*/) {
//printf("in extract\n");
        if (strcmp(str, "smd/ulsph/smoothVel_ptr") == 0) {
                return (void *) smoothVel;
        } else if (strcmp(str, "smd/ulsph/stressTensor_ptr") == 0) {
                return (void *) stressTensor;
        } else if (strcmp(str, "smd/ulsph/velocityGradient_ptr") == 0) {
                return (void *) L;
        } else if (strcmp(str, "smd/ulsph/numNeighs_ptr") == 0) {
                return (void *) numNeighs;
        } else if (strcmp(str, "smd/ulsph/dtCFL_ptr") == 0) {
//printf("dtcfl = %f\n", dtCFL);
                return (void *) &dtCFL;
        } else if (strcmp(str, "smd/ulsph/updateFlag_ptr") == 0) {
                return (void *) &updateFlag;
        } else if (strcmp(str, "smd/ulsph/effective_modulus_ptr") == 0) {
                return (void *) effm;
        } else if (strcmp(str, "smd/ulsph/shape_matrix_ptr") == 0) {
                return (void *) K;
        }

        return nullptr;
}

/* ----------------------------------------------------------------------
 compute effective shear modulus by dividing rate of deviatoric stress with rate of shear deformation
 ------------------------------------------------------------------------- */

double PairULSPH::effective_shear_modulus(const Matrix3d& d_dev, const Matrix3d& deltaStressDev, const double dt, const int itype) {
        double G_eff; // effective shear modulus, see Pronto 2d eq. 3.4.7
        double deltaStressDevSum, shearRateSq, strain_increment;

        if (domain->dimension == 3) {
                deltaStressDevSum = deltaStressDev(0, 1) * deltaStressDev(0, 1) + deltaStressDev(0, 2) * deltaStressDev(0, 2)
                                + deltaStressDev(1, 2) * deltaStressDev(1, 2);
                shearRateSq = d_dev(0, 1) * d_dev(0, 1) + d_dev(0, 2) * d_dev(0, 2) + d_dev(1, 2) * d_dev(1, 2);
        } else {
                deltaStressDevSum = deltaStressDev(0, 1) * deltaStressDev(0, 1);
                shearRateSq = d_dev(0, 1) * d_dev(0, 1);
        }

        strain_increment = dt * dt * shearRateSq;

        if (strain_increment > 1.0e-12) {
                G_eff = 0.5 * sqrt(deltaStressDevSum / strain_increment);
        } else {
                if (strength[itype] != NONE) {
                        G_eff = Lookup[SHEAR_MODULUS][itype];
                } else {
                        G_eff = 0.0;
                }
        }

        return G_eff;

}

/* ----------------------------------------------------------------------
 hourglass force for updated Lagrangian SPH
 input: particles indices i, j, particle types ityep, jtype
 ------------------------------------------------------------------------- */

Vector3d PairULSPH::ComputeHourglassForce(const int i, const int itype, const int j, const int jtype, const Vector3d& dv,
                const Vector3d& xij, const Vector3d& g, const double c_ij, const double mu_ij, const double rho_ij) {

        double *rmass = atom->rmass;
        Vector3d dv_est, f_hg;
        double visc_magnitude;

        dv_est = -0.5 * (F[i] + F[j]) * xij;
        double hurz = dv_est.dot(dv) / (dv_est.norm() * dv.norm() + 1.0e-16);
        if (hurz < 0.0) {

                visc_magnitude = 0.5 * (Q1[itype] + Q1[jtype]) * c_ij * mu_ij / rho_ij;
                f_hg = -rmass[i] * rmass[j] * visc_magnitude * g;
//              printf(" f_hg   = %f %f %f\n", f_hg(0), f_hg(1), f_hg(2));
//              printf("\nnegative\n");
//              printf(" dv_est = %f %f %f\n", dv_est(0), dv_est(1), dv_est(2));
//              printf(" dv     = %f %f %f\n", dv(0), dv(1), dv(2));
        } else {
                f_hg.setZero();
        }

        return f_hg;

}
