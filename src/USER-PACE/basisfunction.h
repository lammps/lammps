//
// Created by lysogy36 on 26.02.20.
//

#ifndef ACE_BASISFUNCTION_H
#define ACE_BASISFUNCTION_H

#include <cstring>

#include "ace_types.h"

//macros for copying the member-array from "other" object
#define basis_mem_copy(other, array, size, type) if(other.array) { \
    if(!is_proxy) delete[] array;\
    array = new type[(size)];\
    is_proxy = false;\
    memcpy(array, other.array, (size) * sizeof(type));\
}


struct C_tilde_B_basis_function {
    // computed combinations of (m1, m2, ...) which have non-zero general Clebsch-Gordan coefficient
    MS_TYPE *ms = nullptr; // flattened array of (m1, m2, ..., m_rank)_1 ... (m1, m2, ..., m_rank)_num_of_ms_combinations of size rank*num_of_ms_combinations
    DOUBLE_TYPE *ctildes = nullptr; // C_tilde coefficients [num_of_ms_combinations][1..ndensity] - expansion used for the coeff_1, .. coeff_num_of_ms_combinations
    SPECIES_TYPE *mus = nullptr; //elements of neighbours atoms, lexicographically sorted, size (rank) index: 0..rank-1
    NS_TYPE *ns = nullptr; // indexes for radial part: n_1..n_rank
    LS_TYPE *ls = nullptr; // indexes for l-part: l1, l2, ..., l_rank
    SHORT_INT_TYPE num_of_ms_combinations; // number of combinations length of  ctildes array

    RANK_TYPE rank; // number of A's in multiplication, i.e. rank = 3 corresponds to (l1,l2,l3), four-body interaction
    DENSITY_TYPE ndensity; // density index
    SPECIES_TYPE mu0; // element of central atom

    bool is_half_ms_basis;
    bool is_proxy = false; //whether object is owner of the ms, ns, ls, mus, ctildes memory or just proxy;

    C_tilde_B_basis_function() = default;

    void _copy_from(const C_tilde_B_basis_function &other) {
        rank = other.rank;
        ndensity = other.ndensity;
        mu0 = other.mu0;
        num_of_ms_combinations = other.num_of_ms_combinations;
        is_half_ms_basis = other.is_half_ms_basis;
        is_proxy = other.is_proxy;

        if (!is_proxy) {
            basis_mem_copy(other, mus, rank, SPECIES_TYPE)
            basis_mem_copy(other, ns, rank, NS_TYPE)
            basis_mem_copy(other, ls, rank, LS_TYPE)
            basis_mem_copy(other, ms, num_of_ms_combinations * rank, MS_TYPE)
            basis_mem_copy(other, ctildes, num_of_ms_combinations * ndensity, DOUBLE_TYPE)
        } else { // if src object is proxy, then just copy pointers
            mus = other.mus;
            ns = other.ns;
            ls = other.ls;
            ms = other.ms;
            ctildes = other.ctildes;
        }
    }

    C_tilde_B_basis_function(const C_tilde_B_basis_function &other) {
#ifdef BASIS_FUNCTION_LIFE_CYCLE
        cout << "C_tilde_B_basis_function::copy constructor " << endl;
#endif
        _copy_from(other);
    }

    void _clean() {
        if (!is_proxy) {
            delete[] mus;
            mus = nullptr;
            delete[] ns;
            ns = nullptr;
            delete[] ls;
            ls = nullptr;
            delete[] ms;
            ms = nullptr;
            delete[] ctildes;
            ctildes = nullptr;
        }
    }

    ~C_tilde_B_basis_function() {
#ifdef BASIS_FUNCTION_LIFE_CYCLE
        cout << "C_tilde_B_basis_function::destructor " << endl;
#endif
        _clean();
    }

    C_tilde_B_basis_function &operator=(const C_tilde_B_basis_function &other) {
#ifdef BASIS_FUNCTION_LIFE_CYCLE
        cout << "C_tilde_B_basis_function::operator = " << endl;
#endif
        _clean();
        _copy_from(other);
        return *this;
    }

};

void print_C_tilde_B_basis_function(const C_tilde_B_basis_function &func);

#endif //ACE_BASISFUNCTION_H
