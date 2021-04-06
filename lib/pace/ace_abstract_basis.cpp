/*
 * Performant implementation of atomic cluster expansion and interface to LAMMPS
 *
 * Copyright 2021  (c) Yury Lysogorskiy^1, Cas van der Oord^2, Anton Bochkarev^1,
 * Sarath Menon^1, Matteo Rinaldi^1, Thomas Hammerschmidt^1, Matous Mrovec^1,
 * Aidan Thompson^3, Gabor Csanyi^2, Christoph Ortner^4, Ralf Drautz^1
 *
 * ^1: Ruhr-University Bochum, Bochum, Germany
 * ^2: University of Cambridge, Cambridge, United Kingdom
 * ^3: Sandia National Laboratories, Albuquerque, New Mexico, USA
 * ^4: University of British Columbia, Vancouver, BC, Canada
 *
 *
 * See the LICENSE file.
 * This FILENAME is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// Created by Lysogorskiy Yury on 28.04.2020.

#include "ace_abstract_basis.h"

////embedding function
////case nemb = 1 only implementation
////F = sign(x)*(  ( 1 - exp(-(w*x)^3) )*abs(x)^m + ((1/w)^(m-1))*exp(-(w*x)^3)*abs(x) )
//// !! no prefactor wpre
void Fexp(DOUBLE_TYPE x, DOUBLE_TYPE m, DOUBLE_TYPE &F, DOUBLE_TYPE &DF) {
    DOUBLE_TYPE w = 1.e6;
    DOUBLE_TYPE eps = 1e-10;

    DOUBLE_TYPE lambda = pow(1.0 / w, m - 1.0);
    if (abs(x) > eps) {
        DOUBLE_TYPE g;
        DOUBLE_TYPE a = abs(x);
        DOUBLE_TYPE am = pow(a, m);
        DOUBLE_TYPE w3x3 = pow(w * a, 3);
        DOUBLE_TYPE sign_factor = (signbit(x) ? -1 : 1);
        if (w3x3 > 30.0)
            g = 0.0;
        else
            g = exp(-w3x3);

        DOUBLE_TYPE omg = 1.0 - g;
        F = sign_factor * (omg * am + lambda * g * a);
        DOUBLE_TYPE dg = -3.0 * w * w * w * a * a * g;
        DF = m * pow(a, m - 1.0) * omg - am * dg + lambda * dg * a + lambda * g;
    } else {
        F = lambda * x;
        DF = lambda;
    }
}


//Scaled-shifted embedding function
//F = sign(x)*(  ( 1 - exp(-(w*x)^3) )*abs(x)^m + ((1/w)^(m-1))*exp(-(w*x)^3)*abs(x) )
// !! no prefactor wpre
void FexpShiftedScaled(DOUBLE_TYPE rho, DOUBLE_TYPE mexp, DOUBLE_TYPE &F, DOUBLE_TYPE &DF) {
    DOUBLE_TYPE eps = 1e-10;
    DOUBLE_TYPE a, xoff, yoff, nx, exprho;

    if (abs(mexp - 1.0) < eps) {
        F = rho;
        DF = 1;
    } else {
        a = abs(rho);
        exprho = exp(-a);
        nx = 1. / mexp;
        xoff = pow(nx, (nx / (1.0 - nx))) * exprho;
        yoff = pow(nx, (1 / (1.0 - nx))) * exprho;
        DOUBLE_TYPE sign_factor = (signbit(rho) ? -1 : 1);
        F = sign_factor * (pow(xoff + a, mexp) - yoff);
        DF = yoff + mexp * (-xoff + 1.0) * pow(xoff + a, mexp - 1.);
    }
}

void ACEAbstractBasisSet::inner_cutoff(DOUBLE_TYPE rho_core, DOUBLE_TYPE rho_cut, DOUBLE_TYPE drho_cut,
                                       DOUBLE_TYPE &fcut, DOUBLE_TYPE &dfcut) {

    DOUBLE_TYPE rho_low = rho_cut - drho_cut;
    if (rho_core >= rho_cut) {
        fcut = 0;
        dfcut = 0;
    } else if (rho_core <= rho_low) {
        fcut = 1;
        dfcut = 0;
    } else {
        fcut = 0.5 * (1 + cos(M_PI * (rho_core - rho_low) / drho_cut));
        dfcut = -0.5 * sin(M_PI * (rho_core - rho_low) / drho_cut) * M_PI / drho_cut;
    }
}

void ACEAbstractBasisSet::FS_values_and_derivatives(Array1D<DOUBLE_TYPE> &rhos, DOUBLE_TYPE &value,
                                                    Array1D<DOUBLE_TYPE> &derivatives, DENSITY_TYPE ndensity) {
    DOUBLE_TYPE F, DF = 0, wpre, mexp;
    for (int p = 0; p < ndensity; p++) {
        wpre = FS_parameters.at(p * ndensity + 0);
        mexp = FS_parameters.at(p * ndensity + 1);
        if (this->npoti == "FinnisSinclair")
            Fexp(rhos(p), mexp, F, DF);
        else if (this->npoti == "FinnisSinclairShiftedScaled")
            FexpShiftedScaled(rhos(p), mexp, F, DF);
        value += F * wpre; // * weight (wpre)
        derivatives(p) = DF * wpre;// * weight (wpre)
    }
}

void ACEAbstractBasisSet::_clean() {

    delete[] elements_name;
    elements_name = nullptr;
    delete radial_functions;
    radial_functions = nullptr;
}

ACEAbstractBasisSet::ACEAbstractBasisSet(const ACEAbstractBasisSet &other) {
    ACEAbstractBasisSet::_copy_scalar_memory(other);
    ACEAbstractBasisSet::_copy_dynamic_memory(other);
}

ACEAbstractBasisSet &ACEAbstractBasisSet::operator=(const ACEAbstractBasisSet &other) {
    if (this != &other) {
        // deallocate old memory
        ACEAbstractBasisSet::_clean();
        //copy scalar values
        ACEAbstractBasisSet::_copy_scalar_memory(other);
        //copy dynamic memory
        ACEAbstractBasisSet::_copy_dynamic_memory(other);
    }
    return *this;
}

ACEAbstractBasisSet::~ACEAbstractBasisSet() {
    ACEAbstractBasisSet::_clean();
}

void ACEAbstractBasisSet::_copy_scalar_memory(const ACEAbstractBasisSet &src) {
    deltaSplineBins = src.deltaSplineBins;
    FS_parameters = src.FS_parameters;
    npoti = src.npoti;

    nelements = src.nelements;
    rankmax = src.rankmax;
    ndensitymax = src.ndensitymax;
    nradbase = src.nradbase;
    lmax = src.lmax;
    nradmax = src.nradmax;
    cutoffmax = src.cutoffmax;

    spherical_harmonics = src.spherical_harmonics;

    rho_core_cutoffs = src.rho_core_cutoffs;
    drho_core_cutoffs = src.drho_core_cutoffs;


    E0vals = src.E0vals;
}

void ACEAbstractBasisSet::_copy_dynamic_memory(const ACEAbstractBasisSet &src) {//allocate new memory
    if (src.elements_name == nullptr)
        throw runtime_error("Could not copy ACEAbstractBasisSet::elements_name - array not initialized");
    elements_name = new string[nelements];
    //copy
    for (SPECIES_TYPE mu = 0; mu < nelements; ++mu) {
        elements_name[mu] = src.elements_name[mu];
    }
    radial_functions = src.radial_functions->clone();
}

SPECIES_TYPE ACEAbstractBasisSet::get_species_index_by_name(const string &elemname) {
    for (SPECIES_TYPE t = 0; t < nelements; t++) {
        if (this->elements_name[t] == elemname)
            return t;
    }
    return -1;
}