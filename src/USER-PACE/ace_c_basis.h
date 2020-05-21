#ifndef ACE_C_BASIS_H
#define ACE_C_BASIS_H


#include "ace_flatten_basis.h"

typedef vector<vector<ACECTildeBasisFunction>> C_tilde_full_basis_vector2d;

class ACECTildeBasisSet : public ACEFlattenBasisSet {
protected:
    //2D-vector-like basis representation:[mu][func_ind]
    C_tilde_full_basis_vector2d mu0_ctilde_basis_vector;
public:
    //[mu][func_ind]
    ACECTildeBasisFunction **basis_rank1 = nullptr;
    ACECTildeBasisFunction **basis = nullptr;

    //C_tilde coefficients contiguous package, size: coeff_array_total_size_rank1
    DOUBLE_TYPE *full_c_tildes_rank1 = nullptr;
    //C_tilde coefficients contiguous package, size: coeff_array_total_size
    DOUBLE_TYPE *full_c_tildes = nullptr;

    SHORT_INT_TYPE num_ctilde_max = 0;

    ACECTildeBasisSet() = default;

    // copy constructor, operator= and destructor (see. Rule of Three)
    ACECTildeBasisSet(const ACECTildeBasisSet &other);

    ACECTildeBasisSet &operator=(const ACECTildeBasisSet &other);

    ~ACECTildeBasisSet() override;

    // routines for copying and cleaning dynamic memory of the class (see. Rule of Three)
    void _clean() override;

    void _copy_dynamic_memory(const ACECTildeBasisSet &src);

    void _copy_scalar_memory(const ACECTildeBasisSet &src);

    void _clean_contiguous_arrays();

    void save(const string &filename) override;

    void load(string filename) override;

    void pack_flatten_basis() override;

    void compute_array_sizes(ACECTildeBasisFunction **basis_rank1, ACECTildeBasisFunction **basis);

    void _clean_basis_arrays();

    void flatten_basis(C_tilde_full_basis_vector2d &mu0_ctilde_basis_vector);

    void flatten_basis() override {};
};

#endif //ACE_C_BASIS_H
