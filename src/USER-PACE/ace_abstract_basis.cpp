//
// Created by Lysogorskiy Yury on 28.04.2020.
//


#include "ace_abstract_basis.h"

//embedding function
//case nemb = 1 only implementation
//F = sign(x)*(  ( ( 1 - exp(-w*x**2) )*abs(x) )^m +  m*exp(-w*x**2)*abs(x) )
// !! no prefactor wpre
void Fexp(DOUBLE_TYPE rho, DOUBLE_TYPE mexp, DOUBLE_TYPE &F, DOUBLE_TYPE &DF) {
    DOUBLE_TYPE w = 10.0;
    DOUBLE_TYPE eps = 1e-10;

    if (abs(rho) > eps) {
        DOUBLE_TYPE g, a, omg, y2, y1 = w * rho * rho;
        DOUBLE_TYPE sign_factor = (signbit(rho) ? -1 : 1);
        if (y1 > 30.0) g = 0;
        else g = exp(-y1);

        omg = 1. - g;
        a = abs(rho);
        y1 = pow(omg * a, mexp);
        y2 = mexp * g * a;
        F = sign_factor * (y1 + y2);

        DOUBLE_TYPE dg, da, dy, dy1, dy2;
        dg = -2.0 * w * rho * g;
        da = sign_factor;
        if (abs(y1) < eps) dy = 0.;
        else dy = mexp * y1 / (omg * a);


        dy1 = dy * (-dg * a + omg * da);
        dy2 = mexp * (dg * a + g * da);
        DF = sign_factor * (dy1 + dy2);

    } else {
        F = mexp * rho;
        DF = mexp;
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
    DOUBLE_TYPE F, DF = 0;
    for (int p = 0; p < ndensity; p++) {
        Fexp(rhos(p), FS_parameters.at(p * ndensity + 1), F, DF);
        value += F * FS_parameters.at(p * ndensity + 0); // * weight
        derivatives(p) = DF * FS_parameters.at(p * ndensity + 0);// * weight
    }
}

void ACEAbstractBasisSet::_clean() {
    if (elements_name != nullptr) {
        delete[] elements_name;
        elements_name = nullptr;
    }
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

void ACEAbstractBasisSet::_copy_scalar_memory(const ACEAbstractBasisSet &other) {
    ntot = other.ntot;
    FS_parameters = other.FS_parameters;

    nelements = other.nelements;
    rankmax = other.rankmax;
    ndensitymax = other.ndensitymax;
    nradbase = other.nradbase;
    lmax = other.lmax;
    nradmax = other.nradmax;
    cutoffmax = other.cutoffmax;

    radial_functions = other.radial_functions;
    spherical_harmonics = other.spherical_harmonics;

    rho_core_cutoffs = other.rho_core_cutoffs;
    drho_core_cutoffs = other.drho_core_cutoffs;


}

void ACEAbstractBasisSet::_copy_dynamic_memory(const ACEAbstractBasisSet &other) {//allocate new memory

    elements_name = new string[nelements];
    //copy
    for (SPECIES_TYPE mu = 0; mu < nelements; ++mu) {
        elements_name[mu] = other.elements_name[mu];
    }
}

SPECIES_TYPE ACEAbstractBasisSet::get_species_index_by_name(const string &elemname) {
    for (SPECIES_TYPE t = 0; t < nelements; t++) {
        if (this->elements_name[t] == elemname)
            return t;
    }
    return -1;
}