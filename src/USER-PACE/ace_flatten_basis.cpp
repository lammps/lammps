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
 */


// Created by yury on 28.04.2020.

#include "ace_flatten_basis.h"

ACEFlattenBasisSet::ACEFlattenBasisSet(const ACEFlattenBasisSet &other) {
    _copy_scalar_memory(other);
    _copy_dynamic_memory(other);
}

ACEFlattenBasisSet &ACEFlattenBasisSet::operator=(const ACEFlattenBasisSet &other) {
    if (this != &other) {
        _clean();
        _copy_scalar_memory(other);
        _copy_dynamic_memory(other);
    }
    return *this;
}


ACEFlattenBasisSet::~ACEFlattenBasisSet() {
    ACEFlattenBasisSet::_clean();
}

void ACEFlattenBasisSet::_clean() {
    //chained call of base class method
    ACEAbstractBasisSet::_clean();
    _clean_contiguous_arrays();
    _clean_basissize_arrays();

}

void ACEFlattenBasisSet::_clean_basissize_arrays() {
    delete[] total_basis_size;
    total_basis_size = nullptr;

    delete[] total_basis_size_rank1;
    total_basis_size_rank1 = nullptr;
}

void ACEFlattenBasisSet::_clean_contiguous_arrays() {
    delete[] full_ns_rank1;
    full_ns_rank1 = nullptr;

    delete[] full_ls_rank1;
    full_ls_rank1 = nullptr;

    delete[] full_mus_rank1;
    full_mus_rank1 = nullptr;

    delete[] full_ms_rank1;
    full_ms_rank1 = nullptr;

    //////

    delete[] full_ns;
    full_ns = nullptr;

    delete[] full_ls;
    full_ls = nullptr;

    delete[] full_mus;
    full_mus = nullptr;

    delete[] full_ms;
    full_ms = nullptr;
}

void ACEFlattenBasisSet::_copy_scalar_memory(const ACEFlattenBasisSet &src) {
    ACEAbstractBasisSet::_copy_scalar_memory(src);

    rank_array_total_size_rank1 = src.rank_array_total_size_rank1;
    coeff_array_total_size_rank1 = src.coeff_array_total_size_rank1;
    rank_array_total_size = src.rank_array_total_size;
    ms_array_total_size = src.ms_array_total_size;
    coeff_array_total_size = src.coeff_array_total_size;

    max_B_array_size = src.max_B_array_size;
    max_dB_array_size = src.max_dB_array_size;
    num_ms_combinations_max = src.num_ms_combinations_max;
}

void ACEFlattenBasisSet::_copy_dynamic_memory(const ACEFlattenBasisSet &src) { //allocate new memory

    ACEAbstractBasisSet::_copy_dynamic_memory(src);

    if (src.total_basis_size_rank1 == nullptr)
        throw runtime_error("Could not copy ACEFlattenBasisSet::total_basis_size_rank1 - array not initialized");
    if (src.total_basis_size == nullptr)
        throw runtime_error("Could not copy ACEFlattenBasisSet::total_basis_size - array not initialized");

    delete[] total_basis_size_rank1;
    total_basis_size_rank1 = new SHORT_INT_TYPE[nelements];
    delete[] total_basis_size;
    total_basis_size = new SHORT_INT_TYPE[nelements];

    //copy
    for (SPECIES_TYPE mu = 0; mu < nelements; ++mu) {
        total_basis_size_rank1[mu] = src.total_basis_size_rank1[mu];
        total_basis_size[mu] = src.total_basis_size[mu];
    }
}

