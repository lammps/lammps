//
// Created by yury on 28.04.2020.
//

#ifndef ACE_EVALUATOR_ACE_FLATTEN_BASIS_H
#define ACE_EVALUATOR_ACE_FLATTEN_BASIS_H

#include <vector>

#include "ace_c_basisfunction.h"
#include "ace_radial.h"
#include "ace_spherical_cart.h"
#include "ace_types.h"
#include "ace_abstract_basis.h"

class ACEFlattenBasisSet : public ACEAbstractBasisSet {
protected:
    //arrays and its sizes for rank = 1 basis functions for packed basis
    size_t rank_array_total_size_rank1 = 0;
    size_t coeff_array_total_size_rank1 = 0;
    //ns contiguous package
    NS_TYPE *full_ns_rank1 = nullptr;
    //ls contiguous package
    LS_TYPE *full_ls_rank1 = nullptr;
    //mus contiguous package
    SPECIES_TYPE *full_Xs_rank1 = nullptr;
    //m_s contiguous package
    MS_TYPE *full_ms_rank1 = nullptr;


    //arrays and its sizes for rank > 1 basis functions for packed basis
    size_t rank_array_total_size = 0;
    size_t ms_array_total_size = 0;
    size_t coeff_array_total_size = 0;
    //ns contiguous package
    NS_TYPE *full_ns = nullptr;
    //ls contiguous package
    LS_TYPE *full_ls = nullptr;
    //mus contiguous package
    SPECIES_TYPE *full_Xs = nullptr;
    //m_s contiguous package
    MS_TYPE *full_ms = nullptr;

    // rearrange basis functions in contiguous memory to optimize cache access
    virtual void pack_flatten_basis() = 0;

    // routines for copying and cleaning dynamic memory of the class (see. Rule of Three)
    void _clean() override; //must be idempotent for safety

    void __copy_packed_arrays(const ACEFlattenBasisSet &other);

    void _copy_dynamic_memory(const ACEFlattenBasisSet &other);

    void _copy_scalar_memory(const ACEFlattenBasisSet &other);

    virtual void flatten_basis() = 0;

public:

    //1D flat array basis representation: [mu]
    SHORT_INT_TYPE *total_basis_size = nullptr;
    SHORT_INT_TYPE *total_basis_size_rank1 = nullptr;

    // maximum over elements array size for B[func_ind][ms_ind] and dB[func_ind][ms_ind][r]
    size_t max_B_array_size = 0;
    size_t max_dB_array_size = 0;

    SHORT_INT_TYPE num_ms_combinations_max = 0;

    ACEFlattenBasisSet() = default;

    // copy constructor, operator= and destructor (see. Rule of Three)
    ACEFlattenBasisSet(const ACEFlattenBasisSet &other);

    ACEFlattenBasisSet &operator=(const ACEFlattenBasisSet &other);

    ~ACEFlattenBasisSet() override;
};

#endif //ACE_EVALUATOR_ACE_FLATTEN_BASIS_H
