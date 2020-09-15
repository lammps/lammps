/* ----------------------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the USER-SMD package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */

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
#include "smd_material_models.h"
#include <cmath>
#include <cstdlib>
#include <utility>
#include <iostream>
#include <cstdio>
#include "math_special.h"

#include <Eigen/Eigen>

using namespace LAMMPS_NS::MathSpecial;
using namespace std;
using namespace Eigen;

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

/* ----------------------------------------------------------------------
 linear EOS for use with linear elasticity
 input: initial pressure pInitial, isotropic part of the strain rate d, time-step dt
 output: final pressure pFinal, pressure rate p_rate
 ------------------------------------------------------------------------- */
void LinearEOS(double lambda, double pInitial, double d, double dt, double &pFinal, double &p_rate) {

        /*
         * pressure rate
         */
        p_rate = lambda * d;

        pFinal = pInitial + dt * p_rate; // increment pressure using pressure rate
        //cout << "hurz" << endl;

}

/* ----------------------------------------------------------------------
 shock EOS
 input:
 current density rho
 reference density rho0
 current energy density e
 reference energy density e0
 reference speed of sound c0
 shock Hugoniot parameter S
 Grueneisen parameter Gamma
 initial pressure pInitial
 time step dt

 output:
 pressure rate p_rate
 final pressure pFinal

 ------------------------------------------------------------------------- */
void ShockEOS(double rho, double rho0, double e, double e0, double c0, double S, double Gamma, double pInitial, double dt,
                double &pFinal, double &p_rate) {

        double mu = rho / rho0 - 1.0;
        double pH = rho0 * square(c0) * mu * (1.0 + mu) / square(1.0 - (S - 1.0) * mu);

        pFinal = (-pH + rho * Gamma * (e - e0));

        //printf("shock EOS: rho = %g, rho0 = %g, Gamma=%f, c0=%f, S=%f, e=%f, e0=%f\n", rho, rho0, Gamma, c0, S, e, e0);
        //printf("pFinal = %f\n", pFinal);
        p_rate = (pFinal - pInitial) / dt;

}

/* ----------------------------------------------------------------------
 polynomial EOS
 input:
 current density rho
 reference density rho0
 coefficients 0 .. 6
 initial pressure pInitial
 time step dt

 output:
 pressure rate p_rate
 final pressure pFinal

 ------------------------------------------------------------------------- */
void polynomialEOS(double rho, double rho0, double /*e*/, double C0, double C1, double C2, double C3, double /*C4*/, double /*C5*/, double /*C6*/,
                double pInitial, double dt, double &pFinal, double &p_rate) {

        double mu = rho / rho0 - 1.0;

        if (mu > 0.0) {
                pFinal = C0 + C1 * mu + C2 * mu * mu + C3 * mu * mu * mu; // + (C4 + C5 * mu + C6 * mu * mu) * e;
        } else {
                pFinal = C0 + C1 * mu + C3 * mu * mu * mu; //  + (C4 + C5 * mu) * e;
        }
        pFinal = -pFinal; // we want the mean stress, not the pressure.


        //printf("pFinal = %f\n", pFinal);
        p_rate = (pFinal - pInitial) / dt;

}

/* ----------------------------------------------------------------------
 Tait EOS based on current density vs. reference density.

 input: (1) reference sound speed
 (2) equilibrium mass density
 (3) current mass density

 output:(1) pressure
 (2) current speed of sound
 ------------------------------------------------------------------------- */
void TaitEOS_density(const double exponent, const double c0_reference, const double rho_reference, const double rho_current,
                double &pressure, double &sound_speed) {

        double B = rho_reference * c0_reference * c0_reference / exponent;
        double tmp = pow(rho_current / rho_reference, exponent);
        pressure = B * (tmp - 1.0);
        double bulk_modulus = B * tmp * exponent; // computed as rho * d(pressure)/d(rho)
        sound_speed = sqrt(bulk_modulus / rho_current);

//      if (fabs(pressure) > 0.01) {
//              printf("tmp = %f, press=%f, K=%f\n", tmp, pressure, bulk_modulus);
//      }

}

/* ----------------------------------------------------------------------
 perfect gas EOS
 input: gamma -- adiabatic index (ratio of specific heats)
 J -- determinant of deformation gradient
 volume0 -- reference configuration volume of particle
 energy -- energy of particle
 pInitial -- initial pressure of the particle
 d -- isotropic part of the strain rate tensor,
 dt -- time-step size

 output: final pressure pFinal, pressure rate p_rate
 ------------------------------------------------------------------------- */
void PerfectGasEOS(const double gamma, const double vol, const double mass, const double energy, double &pFinal, double &c0) {

        /*
         * perfect gas EOS is p = (gamma - 1) rho e
         */

        if (energy > 0.0) {

                pFinal = (1.0 - gamma) * energy / vol;
//printf("gamma = %f, vol%f, e=%g ==> p=%g\n", gamma, vol, energy, *pFinal__/1.0e-9);

                c0 = sqrt((gamma - 1.0) * energy / mass);

        } else {
                pFinal = c0 = 0.0;
        }

}

/* ----------------------------------------------------------------------
 linear strength model for use with linear elasticity
 input: lambda, mu : Lame parameters
 input: sigmaInitial_dev, d_dev: initial stress deviator, deviatoric part of the strain rate tensor
 input: dt: time-step
 output:  sigmaFinal_dev, sigmaFinal_dev_rate__: final stress deviator and its rate.
 ------------------------------------------------------------------------- */
void LinearStrength(const double mu, const Matrix3d sigmaInitial_dev, const Matrix3d d_dev, const double dt,
                Matrix3d &sigmaFinal_dev__, Matrix3d &sigma_dev_rate__) {

        /*
         * deviatoric rate of unrotated stress
         */
        sigma_dev_rate__ = 2.0 * mu * d_dev;

        /*
         * elastic update to the deviatoric stress
         */
        sigmaFinal_dev__ = sigmaInitial_dev + dt * sigma_dev_rate__;
}

/* ----------------------------------------------------------------------
 linear strength model for use with linear elasticity
 input: lambda, mu : Lame parameters
 input: F: deformation gradient
 output:  total stress tensor, deviator + pressure
 ------------------------------------------------------------------------- */
//void PairTlsph::LinearStrengthDefgrad(double lambda, double mu, Matrix3d F, Matrix3d *T) {
//      Matrix3d E, PK2, eye, sigma, S, tau;
//
//      eye.setIdentity();
//
//      E = 0.5 * (F * F.transpose() - eye); // strain measure E = 0.5 * (B - I) = 0.5 * (F * F^T - I)
//      tau = lambda * E.trace() * eye + 2.0 * mu * E; // Kirchhoff stress, work conjugate to above strain
//      sigma = tau / F.determinant(); // convert Kirchhoff stress to Cauchy stress
//
////printf("l=%f, mu=%f, sigma xy = %f\n", lambda, mu, sigma(0,1));
//
////    E = 0.5 * (F.transpose() * F - eye); // Green-Lagrange Strain E = 0.5 * (C - I)
////    S = lambda * E.trace() * eye + 2.0 * mu * Deviator(E); // PK2 stress
////    tau = F * S * F.transpose(); // convert PK2 to Kirchhoff stress
////    sigma = tau / F.determinant();
//
//      //*T = sigma;
//
//      /*
//       * neo-hookean model due to Bonet
//       */
////    lambda = mu = 100.0;
////    // left Cauchy-Green Tensor, b = F.F^T
//      double J = F.determinant();
//      double logJ = log(J);
//      Matrix3d b;
//      b = F * F.transpose();
//
//      sigma = (mu / J) * (b - eye) + (lambda / J) * logJ * eye;
//      *T = sigma;
//}
/* ----------------------------------------------------------------------
 linear strength model for use with linear elasticity
 input: lambda, mu : Lame parameters
 input: sigmaInitial_dev, d_dev: initial stress deviator, deviatoric part of the strain rate tensor
 input: dt: time-step
 output:  sigmaFinal_dev, sigmaFinal_dev_rate__: final stress deviator and its rate.
 ------------------------------------------------------------------------- */
void LinearPlasticStrength(const double G, const double yieldStress, const Matrix3d sigmaInitial_dev, const Matrix3d d_dev,
                const double dt, Matrix3d &sigmaFinal_dev__, Matrix3d &sigma_dev_rate__, double &plastic_strain_increment) {

        Matrix3d sigmaTrial_dev, dev_rate;
        double J2;

        /*
         * deviatoric rate of unrotated stress
         */
        dev_rate = 2.0 * G * d_dev;

        /*
         * perform a trial elastic update to the deviatoric stress
         */
        sigmaTrial_dev = sigmaInitial_dev + dt * dev_rate; // increment stress deviator using deviatoric rate

        /*
         * check yield condition
         */
        J2 = sqrt(3. / 2.) * sigmaTrial_dev.norm();

        if (J2 < yieldStress) {
                /*
                 * no yielding has occurred.
                 * final deviatoric stress is trial deviatoric stress
                 */
                sigma_dev_rate__ = dev_rate;
                sigmaFinal_dev__ = sigmaTrial_dev;
                plastic_strain_increment = 0.0;
                //printf("no yield\n");

        } else {
                //printf("yiedl\n");
                /*
                 * yielding has occurred
                 */
                plastic_strain_increment = (J2 - yieldStress) / (3.0 * G);

                /*
                 * new deviatoric stress:
                 * obtain by scaling the trial stress deviator
                 */
                sigmaFinal_dev__ = (yieldStress / J2) * sigmaTrial_dev;

                /*
                 * new deviatoric stress rate
                 */
                sigma_dev_rate__ = sigmaFinal_dev__ - sigmaInitial_dev;
                //printf("yielding has occurred.\n");
        }
}

/* ----------------------------------------------------------------------
 Johnson Cook Material Strength model
 input:
 G : shear modulus
 cp : heat capacity
 espec : energy / mass
 A : initial yield stress under quasi-static / room temperature conditions
 B : proportionality factor for plastic strain dependency
 a : exponent for plastic strain dpendency
 C : proportionality factor for logarithmic plastic strain rate dependency
 epdot0 : dimensionality factor for plastic strain rate dependency
 T : current temperature
 T0 : reference (room) temperature
 Tmelt : melting temperature
 input: sigmaInitial_dev, d_dev: initial stress deviator, deviatoric part of the strain rate tensor
 input: dt: time-step
 output:  sigmaFinal_dev, sigmaFinal_dev_rate__: final stress deviator and its rate.
 ------------------------------------------------------------------------- */
void JohnsonCookStrength(const double G, const double cp, const double espec, const double A, const double B, const double a,
                const double C, const double epdot0, const double T0, const double Tmelt, const double /*M*/, const double dt, const double ep,
                const double epdot, const Matrix3d sigmaInitial_dev, const Matrix3d d_dev, Matrix3d &sigmaFinal_dev__,
                Matrix3d &sigma_dev_rate__, double &plastic_strain_increment) {

        Matrix3d sigmaTrial_dev, dev_rate;
        double J2, yieldStress;

        double deltaT = espec / cp;
        double TH = deltaT / (Tmelt - T0);
        TH = MAX(TH, 0.0);
        double epdot_ratio = epdot / epdot0;
        epdot_ratio = MAX(epdot_ratio, 1.0);
        //printf("current temperature delta is %f, TH=%f\n", deltaT, TH);

        yieldStress = (A + B * pow(ep, a)) * (1.0 + C * log(epdot_ratio)); // * (1.0 - pow(TH, M));

        /*
         * deviatoric rate of unrotated stress
         */
        dev_rate = 2.0 * G * d_dev;

        /*
         * perform a trial elastic update to the deviatoric stress
         */
        sigmaTrial_dev = sigmaInitial_dev + dt * dev_rate; // increment stress deviator using deviatoric rate

        /*
         * check yield condition
         */
        J2 = sqrt(3. / 2.) * sigmaTrial_dev.norm();

        if (J2 < yieldStress) {
                /*
                 * no yielding has occurred.
                 * final deviatoric stress is trial deviatoric stress
                 */
                sigma_dev_rate__ = dev_rate;
                sigmaFinal_dev__ = sigmaTrial_dev;
                plastic_strain_increment = 0.0;
                //printf("no yield\n");

        } else {
                //printf("yiedl\n");
                /*
                 * yielding has occurred
                 */
                plastic_strain_increment = (J2 - yieldStress) / (3.0 * G);

                /*
                 * new deviatoric stress:
                 * obtain by scaling the trial stress deviator
                 */
                sigmaFinal_dev__ = (yieldStress / J2) * sigmaTrial_dev;

                /*
                 * new deviatoric stress rate
                 */
                sigma_dev_rate__ = sigmaFinal_dev__ - sigmaInitial_dev;
                //printf("yielding has occurred.\n");
        }
}

/* ----------------------------------------------------------------------
 isotropic maximum strain damage model
 input:
 current strain
 maximum value of allowed principal strain

 output:
 return value is true if any eigenvalue of the current strain exceeds the allowed principal strain

 ------------------------------------------------------------------------- */

bool IsotropicMaxStrainDamage(const Matrix3d E, const double maxStrain) {

        /*
         * compute Eigenvalues of strain matrix
         */
        SelfAdjointEigenSolver < Matrix3d > es;
        es.compute(E); // compute eigenvalue and eigenvectors of strain

        double max_eigenvalue = es.eigenvalues().maxCoeff();

        if (max_eigenvalue > maxStrain) {
                return true;
        } else {
                return false;
        }
}

/* ----------------------------------------------------------------------
 isotropic maximum stress damage model
 input:
 current stress
 maximum value of allowed principal stress

 output:
 return value is true if any eigenvalue of the current stress exceeds the allowed principal stress

 ------------------------------------------------------------------------- */

bool IsotropicMaxStressDamage(const Matrix3d S, const double maxStress) {

        /*
         * compute Eigenvalues of strain matrix
         */
        SelfAdjointEigenSolver < Matrix3d > es;
        es.compute(S); // compute eigenvalue and eigenvectors of strain

        double max_eigenvalue = es.eigenvalues().maxCoeff();

        if (max_eigenvalue > maxStress) {
                return true;
        } else {
                return false;
        }
}

/* ----------------------------------------------------------------------
 Johnson-Cook failure model
 input:


 output:


 ------------------------------------------------------------------------- */

double JohnsonCookFailureStrain(const double p, const Matrix3d Sdev, const double d1, const double d2, const double d3,
                const double d4, const double epdot0, const double epdot) {



        double vm = sqrt(3. / 2.) * Sdev.norm(); // von-Mises equivalent stress
        if (vm < 0.0) {
                cout << "this is sdev " << endl << Sdev << endl;
                printf("vm=%f < 0.0, surely must be an error\n", vm);
                exit(1);
        }

        // determine stress triaxiality
        double triax = p / (vm + 0.01 * fabs(p)); // have softening in denominator to avoid division by zero
        if (triax < 0.0) {
                triax = 0.0;
        } else if (triax > 3.0) {
                triax = 3.0;
        }

        // Johnson-Cook failure strain, dependence on stress triaxiality
        double jc_failure_strain = d1 + d2 * exp(d3 * triax);

        // include strain rate dependency if parameter d4 is defined and current plastic strain rate exceeds reference strain rate
        if (d4 > 0.0) { //
                if (epdot > epdot0) {
                        double epdot_ratio = epdot / epdot0;
                        jc_failure_strain *= (1.0 + d4 * log(epdot_ratio));
                        //printf("epsdot=%f, epsdot0=%f, factor = %f\n", epdot, epdot0, (1.0 + d4 * log(epdot_ratio)));
                        //exit(1);

                }
        }

        return jc_failure_strain;

}
