//
// Created by yury on 28.04.2020.
//

#ifndef ACE_EVALUATOR_ACE_ABSTRACT_BASIS_H
#define ACE_EVALUATOR_ACE_ABSTRACT_BASIS_H

#include <vector>
#include <string>
#include "ace_types.h"
#include "ace_radial.h"
#include "ace_contigous_array.h"
#include "ace_spherical_cart.h"

using namespace std;

class ACEAbstractBasisSet {
protected:
    vector<DOUBLE_TYPE> FS_parameters;// for Eq.(53)

    int ntot = 0;

    // routines for copying and cleaning dynamic memory of the class (see. Rule of Three)
    virtual void _clean(); //must be idempotent for safety
    virtual void _copy_dynamic_memory(const ACEAbstractBasisSet &other);

    virtual void _copy_scalar_memory(const ACEAbstractBasisSet &other);

public:
    SPECIES_TYPE nelements = 0; //number of elements
    RANK_TYPE rankmax = 0;
    DENSITY_TYPE ndensitymax = 0;
    NS_TYPE nradbase = 0;
    LS_TYPE lmax = 0;
    NS_TYPE nradmax = 0;
    DOUBLE_TYPE cutoffmax = 0;

    string *elements_name = nullptr; //mapping index(0..nelements) -> element_symbol

    ACERadialFunctions radial_functions;
    ACECartesianSphericalHarmonics spherical_harmonics;

    //energy-based inner cuttoff
    Array1D<DOUBLE_TYPE> rho_core_cutoffs;
    Array1D<DOUBLE_TYPE> drho_core_cutoffs;

    ACEAbstractBasisSet() = default;

    // copy constructor, operator= and destructor (see. Rule of Three)
    ACEAbstractBasisSet(const ACEAbstractBasisSet &other);

    ACEAbstractBasisSet &operator=(const ACEAbstractBasisSet &other);

    //virtual destructor
    virtual ~ACEAbstractBasisSet();


    void FS_values_and_derivatives(Array1D<DOUBLE_TYPE> &rhos, DOUBLE_TYPE &value, Array1D<DOUBLE_TYPE> &derivatives,
                                   DENSITY_TYPE ndensity);

    static void inner_cutoff(DOUBLE_TYPE rho_core, DOUBLE_TYPE rho_cut, DOUBLE_TYPE drho_cut, DOUBLE_TYPE &fcut,
                             DOUBLE_TYPE &dfcut);


    virtual void save(const string &filename) = 0;

    virtual void load(string filename) = 0;

    SPECIES_TYPE get_species_index_by_name(const string &elemname);
};

#include <vector>
#include "ace_c_basisfunction.h"
#include "ace_radial.h"
#include "ace_spherical_cart.h"
#include "ace_types.h"

#endif //ACE_EVALUATOR_ACE_ABSTRACT_BASIS_H
