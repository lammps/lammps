/*
 * Performant implementation of atomic cluster expansion and interface to LAMMPS
 *
 * Copyright 2021  (c) Yury Lysogorskiy^1, Cas van der Oord^2, Anton Bochkarev^1,
 * Sarath Menon^1, Matteo Rinaldi^1, Thomas Hammerschmidt^1, Matous Mrovec^1,
 * Aidan Thompson^3, Gabor Csanyi^2, Christoph Ortner^4, Ralf Drautz^1
 *
 * ^1: Ruhr-University Bochum, Bochum, Germany
 * ^2: University of Cambridge, Cambridge, United Kingdom
 * ^3: Sandia National Laboratories, Albuquerque, New Mexico, USA
 * ^4: University of British Columbia, Vancouver, BC, Canada
 *
 *
 * See the LICENSE file.
 * This FILENAME is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// Created by Yury Lysogorskiy on 1.04.20.

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
    void _load_radial_ACERadial(FILE *fptr,
                                const string filename,
                                const string radbasename);

    void _load_radial_SHIPsBasic(FILE * fptr, 
                                 const string filename, 
                                 const string radbasename ); 

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
