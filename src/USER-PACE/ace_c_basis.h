#ifndef ACE_C_BASIS_H
#define ACE_C_BASIS_H

#include <vector>

#include "ace_c_basisfunction.h"
#include "ace_radial.h"
#include "ace_spherical_cart.h"
#include "ace_types.h"

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

    //TODO: implement write basis to file
    virtual void save(const string &filename) = 0;

    virtual void load(string filename) = 0;

    SPECIES_TYPE get_species_index_by_name(const string &elemname);
};

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

    //TODO: implement write basis to file
    void save(const string &filename) override;
    void load(string filename) override;
};

#endif //ACE_C_BASIS_H
