#ifndef ACE_C_BASIS_H
#define ACE_C_BASIS_H

#include "ace_flatten_basis.h"

typedef vector<vector<ACECTildeBasisFunction>> C_tilde_full_basis_vector2d;

/**
 * ACE basis set of C-tilde basis functions
 */
class ACECTildeBasisSet : public ACEFlattenBasisSet {
public:

    ACECTildeBasisFunction **basis_rank1 = nullptr; ///< two-dimensional array of first-rank basis function with indices: [species index][func index]
    ACECTildeBasisFunction **basis = nullptr;  ///< two-dimensional array of higher rank basis function with indices: [species index][func index]

    DOUBLE_TYPE *full_c_tildes_rank1 = nullptr; ///< C_tilde coefficients contiguous package, size: coeff_array_total_size_rank1
    DOUBLE_TYPE *full_c_tildes = nullptr; ///< C_tilde coefficients contiguous package, size: coeff_array_total_size

    //TODO: remove?
    SHORT_INT_TYPE num_ctilde_max = 0;

    /**
     * Default constructor
     */
    ACECTildeBasisSet() = default;

    /**
    * Constructor from .ace file
    */
    ACECTildeBasisSet(const string filename);

    /**
     * Copy constructor (see. Rule of Three)
     * @param other
     */
    ACECTildeBasisSet(const ACECTildeBasisSet &other);

    /**
     *  operator= (see. Rule of Three)
     * @param other
     * @return
     */
    ACECTildeBasisSet &operator=(const ACECTildeBasisSet &other);

    /**
     * Destructor  (see. Rule of Three)
     */
    ~ACECTildeBasisSet() override;

    /**
     * Cleaning dynamic memory of the class (see. Rule of Three)
     */
    void _clean() override;

    /**
     * Copying and cleaning dynamic memory of the class (see. Rule of Three)
     * @param src
     */
    void _copy_dynamic_memory(const ACECTildeBasisSet &src);

    /**
     * Copying scalar variables
     * @param src
     */
    void _copy_scalar_memory(const ACECTildeBasisSet &src);

    /**
     * Clean contiguous arrays (full_c_tildes_rank1, full_c_tildes) and those of base class
     */
    void _clean_contiguous_arrays() override ;

    /**
     * Save potential to .ace file
     * @param filename .ace file name
     */
    void save(const string &filename) override;

    /**
     * Load potential from .ace
     * @param filename .ace file name
     */
    void load(const string filename) override;

    /**
     * Load the ACE type radial basis
     */
    void _load_radial_ChebExpCos(FILE *fptr,
                                 const string filename,
                                 const string radbasename);

    void _load_radial_SHIPsBasic(FILE *fptr,
                                 const string filename,
                                 const string radbasename);

    /**
     * Re-pack the constituent dynamic arrays of all basis functions in contiguous arrays
     */
    void pack_flatten_basis() override;

    /**
     * Computes flatten array sizes
     * @param basis_rank1 two-dimensional array of first-rank ACECTildeBasisFunctions
     * @param basis two-dimensional array of higher-rank ACECTildeBasisFunctions
     */
    void compute_array_sizes(ACECTildeBasisFunction** basis_rank1, ACECTildeBasisFunction** basis);

    /**
     * Clean basis arrays  'basis_rank1' and  'basis'
     */
    void _clean_basis_arrays();

    /**
     * Pack two-dimensional vector of ACECTildeBasisFunction into 1D dynami array with all basis functions
     * @param mu0_ctilde_basis_vector vector<vector<ACECTildeBasisFunction>>
     */
    void flatten_basis(C_tilde_full_basis_vector2d& mu0_ctilde_basis_vector);

    /**
     * Empty stub implementation
     */
    void flatten_basis() override{};
};

#endif //ACE_C_BASIS_H
