/* -*- c++ -*- ----------------------------------------------------------
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

#ifndef SMD_MATERIAL_MODELS_H
#define SMD_MATERIAL_MODELS_H

#include <Eigen/Eigen>

/*
 * EOS models
 */
void LinearEOS(double lambda, double pInitial, double d, double dt, double &pFinal, double &p_rate);
void ShockEOS(double rho, double rho0, double e, double e0, double c0, double S, double Gamma,
              double pInitial, double dt, double &pFinal, double &p_rate);
void polynomialEOS(double rho, double rho0, double e, double C0, double C1, double C2, double C3,
                   double C4, double C5, double C6, double pInitial, double dt, double &pFinal,
                   double &p_rate);
void TaitEOS_density(const double exponent, const double c0_reference, const double rho_reference,
                     const double rho_current, double &pressure, double &sound_speed);
void PerfectGasEOS(const double gamma, const double vol, const double mass, const double energy,
                   double &pFinal__, double &c0);

/*
 * Material strength models
 */
void LinearStrength(const double mu, const Eigen::Matrix3d &sigmaInitial_dev,
                    const Eigen::Matrix3d &d_dev, const double dt,
                    Eigen::Matrix3d &sigmaFinal_dev__, Eigen::Matrix3d &sigma_dev_rate__);
void LinearPlasticStrength(const double G, const double yieldStress,
                           const Eigen::Matrix3d &sigmaInitial_dev, const Eigen::Matrix3d &d_dev,
                           const double dt, Eigen::Matrix3d &sigmaFinal_dev__,
                           Eigen::Matrix3d &sigma_dev_rate__, double &plastic_strain_increment);
void JohnsonCookStrength(const double G, const double cp, const double espec, const double A,
                         const double B, const double a, const double C, const double epdot0,
                         const double T0, const double Tmelt, const double M, const double dt,
                         const double ep, const double epdot,
                         const Eigen::Matrix3d &sigmaInitial_dev, const Eigen::Matrix3d &d_dev,
                         Eigen::Matrix3d &sigmaFinal_dev__, Eigen::Matrix3d &sigma_dev_rate__,
                         double &plastic_strain_increment);

/*
 * Damage models
 */

bool IsotropicMaxStrainDamage(const Eigen::Matrix3d &E, const double maxStrain);
bool IsotropicMaxStressDamage(const Eigen::Matrix3d &E, const double maxStrain);
double JohnsonCookFailureStrain(const double p, const Eigen::Matrix3d &Sdev, const double d1,
                                const double d2, const double d3, const double d4,
                                const double epdot0, const double epdot);

#endif /* SMD_MATERIAL_MODELS_H_ */
