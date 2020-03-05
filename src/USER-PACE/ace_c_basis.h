#ifndef ACE_ACE_C_BASIS_H
#define ACE_ACE_C_BASIS_H

#include <vector>

#include "ace_types.h"
#include "ace_radial.h"
#include "ace_spherical_cart.h"
#include "ace_c_basisfunction.h"


class ACEBasisSet {
protected:
    //arrays and its sizes for rank = 1 basis functions
    size_t rank_array_total_size_rank1 = 0;
    size_t coeff_array_total_size_rank1 = 0;
    NS_TYPE *full_ns_rank1 = nullptr;
    LS_TYPE *full_ls_rank1 = nullptr;
    SPECIES_TYPE *full_Xs_rank1 = nullptr;
    MS_TYPE *full_ms_rank1 = nullptr;
    DOUBLE_TYPE *full_c_tildes_rank1 = nullptr;

    //arrays and its sizes for rank > 1 basis functions
    size_t rank_array_total_size = 0;
    size_t ms_array_total_size = 0;
    size_t coeff_array_total_size = 0;
    NS_TYPE *full_ns = nullptr;
    LS_TYPE *full_ls = nullptr;
    SPECIES_TYPE *full_Xs = nullptr;
    MS_TYPE *full_ms = nullptr;
    DOUBLE_TYPE *full_c_tildes = nullptr;

    vector<DOUBLE_TYPE> parameters;// for Eq.(53)

    DOUBLE_TYPE cutoff = 10.;
    int ntot;

    // rearrange basis functions in contiguous memory for optimize cache access
    void pack_flatten_basis();

    // routines for copying and cleaning dynamic memory of the class (see. Rule of Three)
    void _clean();

    void _copy_dynamic_memory(const ACEBasisSet &other);

    void _copy_scalar_memory(const ACEBasisSet &other);

public:
    SPECIES_TYPE nelements; //number of elements
    string *elements_name; //mapping index(0..nelements) -> element_symbol

    RANK_TYPE rankmax;
    DENSITY_TYPE ndensitymax;

    NS_TYPE nradbase;
    LS_TYPE lmax;
    NS_TYPE nradmax;
    DOUBLE_TYPE cutoffmax;


    //1D flat array basis representation: [mu]
    SHORT_INT_TYPE *total_basis_size = nullptr, *total_basis_size_rank1 = nullptr;

    //[mu][func_ind]
    ACECTildeBasisFunction **basis_rank1 = nullptr;
    ACECTildeBasisFunction **basis = nullptr;

    // maximum over elements array size for B[func_ind][ms_ind] and dB[func_ind][ms_ind][r]
    size_t max_B_array_size = 0;
    size_t max_dB_array_size = 0;

    ACERadialFunctions radial_functions;
    ACECartesianSphericalHarmonics spherical_harmonics;

    SHORT_INT_TYPE num_ctilde_max;
    SHORT_INT_TYPE num_ms_combinations_max;

    ACEBasisSet() = default;

    // copy constructor, operator= and destructor (see. Rule of Three)
    ACEBasisSet(const ACEBasisSet &other);

    ACEBasisSet &operator=(const ACEBasisSet &other);

    ~ACEBasisSet();


    void FS_values_and_derivatives(Array1D<DOUBLE_TYPE> &rhos, DOUBLE_TYPE &value, Array1D<DOUBLE_TYPE> &derivatives,
                                   DENSITY_TYPE ndensity);

    //TODO: implement write basis to file
    void save(string filename);

    void load(string filename);

};

#endif //ACE_ACE_C_BASIS_H
