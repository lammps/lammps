//
// Created by lysogy36 on 26.02.20.
//

#ifndef ACE_ACE_C_BASISFUNCTION_H
#define ACE_ACE_C_BASISFUNCTION_H

#include <cstring>

#include "ace_types.h"

//macros for copying the member-array from "other" object for C-tilde and B-basis
#define basis_mem_copy(other, array, size, type) if(other.array) { \
    if(!is_proxy) delete[] array;\
    array = new type[(size)];\
    is_proxy = false;\
    memcpy(array, other.array, (size) * sizeof(type));\
}


struct ACECTildeBasisFunction {

    // flattened array of  computed combinations of (m1, m2, ..., m_rank)
    // which have non-zero general Clebsch-Gordan coefficient
    // (m1, m2, ..., m_rank)_1 ... (m1, m2, ..., m_rank)_num_ms_combs
    // size =  num_ms_combs * rank
    // effective shape [num_ms_combs][rank]
    MS_TYPE *ms_combs = nullptr;

    // c_tilde coefficients for all densities
    // size = num_ms_combs*ndensity
    // effective shape [num_ms_combs][ndensity]
    DOUBLE_TYPE *ctildes = nullptr;

    // species types of neighbours atoms, should be lexicographically sorted
    // size  = rank
    // effective shape [rank]
    SPECIES_TYPE *mus = nullptr;

    // "ns" - indexes for radial part
    // size  = rank
    // effective shape [rank]
    NS_TYPE *ns = nullptr;

    // "ls" - indexes for l-part: l1, l2, ..., l_rank
    // size  = rank
    // effective shape [rank]
    LS_TYPE *ls = nullptr;

    // number of ms-combinations, zero-th effective dimension of ms_combs arrays
    SHORT_INT_TYPE num_ms_combs;

    // rank, number of atomic base functions "A"s in basis function B
    RANK_TYPE rank;

    // number of densities
    DENSITY_TYPE ndensity;

    SPECIES_TYPE mu0; // element type of central atom

    // wheter ms array contains only "non-negative" ms-combinations.
    // positive ms-combination is when the first non-zero m is positive (0 1 -1)
    // negative ms-combination is when the first non-zero m is negative (0 -1 1)
    bool is_half_ms_basis;

    //flag, whether object is "owner" (i.e. reponsible for memory cleaning) of
    // the ms, ns, ls, mus, ctildes memory or just a proxy
    bool is_proxy = false;

    //default contructor
    ACECTildeBasisFunction() = default;

    // Because the ACECTildeBasisFunction contains dynamically allocated arrays, the Rule of Three should be
    // fullfilled, i.e. copy constructor (copy the dynamic arrays), operator= (release previous arrays and
    // copy the new dynamic arrays) and destructor (release all dynamically allocated memory)

    //copy constructor
    ACECTildeBasisFunction(const ACECTildeBasisFunction &other) {
        _copy_from(other);
    }

    //operator=
    ACECTildeBasisFunction &operator=(const ACECTildeBasisFunction &other) {
        _clean();
        _copy_from(other);
        return *this;
    }

    //destructor
    ~ACECTildeBasisFunction() {
        _clean();
    }

    void _copy_from(const ACECTildeBasisFunction &other) {
        rank = other.rank;
        ndensity = other.ndensity;
        mu0 = other.mu0;
        num_ms_combs = other.num_ms_combs;
        is_half_ms_basis = other.is_half_ms_basis;
        is_proxy = other.is_proxy;

        if (!is_proxy) { // if target not proxy, then copy the whole arrays
            basis_mem_copy(other, mus, rank, SPECIES_TYPE)
            basis_mem_copy(other, ns, rank, NS_TYPE)
            basis_mem_copy(other, ls, rank, LS_TYPE)
            basis_mem_copy(other, ms_combs, num_ms_combs * rank, MS_TYPE)
            basis_mem_copy(other, ctildes, num_ms_combs * ndensity, DOUBLE_TYPE)
        } else { // if target object is proxy, then just copy pointers
            mus = other.mus;
            ns = other.ns;
            ls = other.ls;
            ms_combs = other.ms_combs;
            ctildes = other.ctildes;
        }
    }

    void _clean() {
        //release memory if the structure is not proxy
        if (!is_proxy) {
            delete[] mus;
            delete[] ns;
            delete[] ls;
            delete[] ms_combs;
            delete[] ctildes;
        }

        mus = nullptr;
        ns = nullptr;
        ls = nullptr;
        ms_combs = nullptr;
        ctildes = nullptr;
    }

};

void print_C_tilde_B_basis_function(const ACECTildeBasisFunction &func);

#endif //ACE_ACE_C_BASISFUNCTION_H
