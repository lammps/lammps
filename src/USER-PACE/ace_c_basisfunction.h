//
// Created by Yury Lysogorskiy on 26.02.20.
//

#ifndef ACE_C_BASISFUNCTION_H
#define ACE_C_BASISFUNCTION_H

#include <cstring>
#include <sstream>

#include "ace_types.h"

//macros for copying the member-array from "other" object for C-tilde and B-basis
#define basis_mem_copy(other, array, size, type) if(other.array) { \
    if(!is_proxy) delete[] array;\
    array = new type[(size)];\
    is_proxy = false;\
    memcpy(array, other.array, (size) * sizeof(type));\
}


struct ACEAbstractBasisFunction {
    // flattened array of  computed combinations of (m1, m2, ..., m_rank)
    // which have non-zero general Clebsch-Gordan coefficient
    // (m1, m2, ..., m_rank)_1 ... (m1, m2, ..., m_rank)_num_ms_combs
    // size =  num_ms_combs * rank
    // effective shape [num_ms_combs][rank]
    MS_TYPE *ms_combs = nullptr;


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


    // number of combinations length of  ms_cg array
    SHORT_INT_TYPE num_ms_combs = 0;

    // rank, number of atomic base functions "A"s in basis function B
    RANK_TYPE rank = 0;

    // rankL=rank-2
    RANK_TYPE rankL = 0;

    // number of densities
    DENSITY_TYPE ndensity = 0;

    // species type of central atom
    SPECIES_TYPE mu0 = 0;

    // whether ms array contains only "non-negative" ms-combinations.
    // positive ms-combination is when the first non-zero m is positive (0 1 -1)
    // negative ms-combination is when the first non-zero m is negative (0 -1 1)
    bool is_half_ms_basis = false;

    //flag, whether object is "owner" (i.e. responsible for memory cleaning) of
    // the ms, ns, ls, mus, ctildes memory or just a proxy
    bool is_proxy = false;

    virtual void _copy_from(const ACEAbstractBasisFunction &other) {

        rank = other.rank;
        ndensity = other.ndensity;
        mu0 = other.mu0;
        num_ms_combs = other.num_ms_combs;
        is_half_ms_basis = other.is_half_ms_basis;
        is_proxy = false;

        basis_mem_copy(other, mus, rank, SPECIES_TYPE)
        basis_mem_copy(other, ns, rank, NS_TYPE)
        basis_mem_copy(other, ls, rank, LS_TYPE)
        basis_mem_copy(other, ms_combs, num_ms_combs * rank, MS_TYPE)
    }

    virtual void _clean() {
        //release memory if the structure is not proxy
        if (!is_proxy) {
            delete[] mus;
            delete[] ns;
            delete[] ls;
            delete[] ms_combs;
        }

        mus = nullptr;
        ns = nullptr;
        ls = nullptr;
        ms_combs = nullptr;
    }

};

struct ACECTildeBasisFunction : public ACEAbstractBasisFunction {

    // c_tilde coefficients for all densities
    // size = num_ms_combs*ndensity
    // effective shape [num_ms_combs][ndensity]
    DOUBLE_TYPE *ctildes = nullptr;

    //default constructor
    ACECTildeBasisFunction() = default;

    // Because the ACECTildeBasisFunction contains dynamically allocated arrays, the Rule of Three should be
    // fulfilled, i.e. copy constructor (copy the dynamic arrays), operator= (release previous arrays and
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
        ACEAbstractBasisFunction::_copy_from(other);
        is_proxy = false;
        basis_mem_copy(other, ctildes, num_ms_combs * ndensity, DOUBLE_TYPE)
    }

    void _clean() override {
        ACEAbstractBasisFunction::_clean();
        //release memory if the structure is not proxy
        if (!is_proxy) {
            delete[] ctildes;
        }
        ctildes = nullptr;
    }

    void print() const {
        printf("C_tilde_B_basis_function: ndensity= %d, mu0 = %d ", this->ndensity, this->mu0);
        printf(" XS=(");
        for (RANK_TYPE r = 0; r < this->rank; r++)
            printf("%d ", int(this->mus[r]));
        printf("), ns=(");
        for (RANK_TYPE r = 0; r < this->rank; r++)
            printf("%d ", this->ns[r]);
        printf("), ls=(");
        for (RANK_TYPE r = 0; r < this->rank; r++)
            printf("%d ", this->ls[r]);

        printf("), %d m_s combinations: {\n", this->num_ms_combs);
        for (int i = 0; i < this->num_ms_combs; i++) {
            printf("\t< ");
            for (RANK_TYPE r = 0; r < this->rank; r++)
                printf("%d ", this->ms_combs[i * this->rank + r]);
            printf(" >: c_tilde=");
            for (DENSITY_TYPE p = 0; p < this->ndensity; ++p)
                printf(" %f ", this->ctildes[i * this->ndensity + p]);
            printf("\n");
        }
        if (this->is_proxy)
            printf("proxy ");
        if (this->is_half_ms_basis)
            printf("half_ms_basis ");
        printf("}\n");
    }
};


#endif //ACE_C_BASISFUNCTION_H
