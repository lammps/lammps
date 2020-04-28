//
// Created by yury on 28.04.2020.
//

#include "ace_c_basis.h"
#include "ace_flatten_basis.h"

ACEFlattenBasisSet::~ACEFlattenBasisSet() {
    ACEFlattenBasisSet::_clean();
}

void ACEFlattenBasisSet::_clean() {
    //chained call of base class method
    ACEAbstractBasisSet::_clean();

    if (full_ns != nullptr) {
        delete[] full_ns;
        full_ns = nullptr;
    }
    if (full_ls != nullptr) {
        delete[] full_ls;
        full_ls = nullptr;
    }
    if (full_Xs != nullptr) {
        delete[] full_Xs;
        full_Xs = nullptr;
    }
    if (full_ms != nullptr) {
        delete[] full_ms;
        full_ms = nullptr;
    }

    if (full_ns_rank1 != nullptr) {
        delete[] full_ns_rank1;
        full_ns_rank1 = nullptr;
    }

    if (full_ls_rank1 != nullptr) {
        delete[] full_ls_rank1;
        full_ls_rank1 = nullptr;
    }

    if (full_Xs_rank1 != nullptr) {
        delete[] full_Xs_rank1;
        full_Xs_rank1 = nullptr;
    }

    if (full_ms_rank1 != nullptr) {
        delete[] full_ms_rank1;
        full_ms_rank1 = nullptr;
    }

    if (total_basis_size != nullptr) {
        delete[] total_basis_size;
        total_basis_size = nullptr;
    }

    if (total_basis_size_rank1 != nullptr) {
        delete[] total_basis_size_rank1;
        total_basis_size_rank1 = nullptr;
    }


}

void ACEFlattenBasisSet::_copy_scalar_memory(const ACEFlattenBasisSet &other) {
    ACEAbstractBasisSet::_copy_scalar_memory(other);

    rank_array_total_size_rank1 = other.rank_array_total_size_rank1;
    coeff_array_total_size_rank1 = other.coeff_array_total_size_rank1;
    rank_array_total_size = other.rank_array_total_size;
    ms_array_total_size = other.ms_array_total_size;
    coeff_array_total_size = other.coeff_array_total_size;

    max_B_array_size = other.max_B_array_size;
    max_dB_array_size = other.max_dB_array_size;
    num_ms_combinations_max = other.num_ms_combinations_max;
}

void ACEFlattenBasisSet::_copy_dynamic_memory(const ACEFlattenBasisSet &other) {//allocate new memory
    ACEAbstractBasisSet::_copy_dynamic_memory(other);

    __copy_packed_arrays(other);

    total_basis_size = new SHORT_INT_TYPE[nelements];
    total_basis_size_rank1 = new SHORT_INT_TYPE[nelements];

    //copy
    for (SPECIES_TYPE mu = 0; mu < nelements; ++mu) {
        total_basis_size_rank1[mu] = other.total_basis_size_rank1[mu];
        total_basis_size[mu] = other.total_basis_size[mu];
    }
}

void ACEFlattenBasisSet::__copy_packed_arrays(const ACEFlattenBasisSet &other) {
    full_ns_rank1 = new NS_TYPE[other.rank_array_total_size_rank1];
    full_ls_rank1 = new LS_TYPE[other.rank_array_total_size_rank1];
    full_Xs_rank1 = new SPECIES_TYPE[other.rank_array_total_size_rank1];

    full_ms_rank1 = new MS_TYPE[other.rank_array_total_size_rank1];

    for (size_t i = 0; i < other.rank_array_total_size_rank1; ++i) {
        full_ns_rank1[i] = other.full_ns_rank1[i];
        full_Xs_rank1[i] = other.full_Xs_rank1[i];
        full_ls_rank1[i] = other.full_ls_rank1[i];
        full_ms_rank1[i] = other.full_ms_rank1[i];
    }

    //arrays and its sizes for rank > 1 basis functions for packed basis

    full_ns = new NS_TYPE[other.rank_array_total_size];
    full_ls = new LS_TYPE[other.rank_array_total_size];
    full_Xs = new SPECIES_TYPE[other.rank_array_total_size];

    for (size_t i = 0; i < other.rank_array_total_size; ++i) {
        full_ns[i] = other.full_ns[i];
        full_ls[i] = other.full_ls[i];
        full_Xs[i] = other.full_Xs[i];
    }

    full_ms = new MS_TYPE[other.ms_array_total_size];

    for (size_t i = 0; i < other.ms_array_total_size; ++i) {
        full_ms[i] = other.full_ms[i];
    }

}

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