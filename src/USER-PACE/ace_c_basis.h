#ifndef ACE_C_BASIS_H
#define ACE_C_BASIS_H


#include "ace_flatten_basis.h"

class ACECTildeBasisSet : public ACEFlattenBasisSet {
protected:
    //C_tilde coefficients contiguous package
    //[mu]
    DOUBLE_TYPE *full_c_tildes_rank1 = nullptr;
    DOUBLE_TYPE *full_c_tildes = nullptr;

    void pack_flatten_basis() override;

    // routines for copying and cleaning dynamic memory of the class (see. Rule of Three)
    void _clean() override;

    void _copy_dynamic_memory(const ACECTildeBasisSet &other);

    void _copy_scalar_memory(const ACECTildeBasisSet &other);

    void flatten_basis() override {};
public:
    //[mu][func_ind]
    ACECTildeBasisFunction **basis_rank1 = nullptr;
    ACECTildeBasisFunction **basis = nullptr;

    SHORT_INT_TYPE num_ctilde_max = 0;

    ACECTildeBasisSet() = default;

    // copy constructor, operator= and destructor (see. Rule of Three)
    ACECTildeBasisSet(const ACECTildeBasisSet &other);

    ACECTildeBasisSet &operator=(const ACECTildeBasisSet &other);

    ~ACECTildeBasisSet() override;

    void save(const string &filename) override;
    void load(string filename) override;
};

#endif //ACE_C_BASIS_H
