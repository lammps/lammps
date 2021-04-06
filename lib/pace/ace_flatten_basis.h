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

// Created by Lysogroskiy Yury on 28.04.2020.

#ifndef ACE_EVALUATOR_ACE_FLATTEN_BASIS_H
#define ACE_EVALUATOR_ACE_FLATTEN_BASIS_H


#include "ace_abstract_basis.h"
#include "ace_c_basisfunction.h"
#include "ace_radial.h"
#include "ace_spherical_cart.h"
#include "ace_types.h"

/**
 * Basis set with basis function attributes, i.e. \f$ \mathbf{n}, \mathbf{l}, \mathbf{m}\f$, etc.
 * packed into contiguous arrays for the cache-friendly memory layout.
 */
class ACEFlattenBasisSet : public ACEAbstractBasisSet {
public:
    //arrays and its sizes for rank = 1 basis functions for packed basis

    size_t rank_array_total_size_rank1 = 0; ///< size for full_ns_rank1, full_ls_rank1, full_Xs_rank1
    size_t coeff_array_total_size_rank1 = 0; ///< size for full coefficients array (depends on B or C-Tilde basis)

    NS_TYPE *full_ns_rank1 = nullptr; ///<  ns contiguous package [rank_array_total_size_rank1]
    LS_TYPE *full_ls_rank1 = nullptr; ///< ls contiguous package [rank_array_total_size_rank1]
    SPECIES_TYPE *full_mus_rank1 = nullptr; ///< mus contiguous package [rank_array_total_size_rank1]
    MS_TYPE *full_ms_rank1 = nullptr; ///< m_s contiguous package[rank_array_total_size_rank1]

    //arrays and its sizes for rank > 1 basis functions for packed basis
    size_t rank_array_total_size = 0; ///< size for full_ns, full_ls, full_Xs
    size_t ms_array_total_size = 0; ///< size for full_ms array
    size_t coeff_array_total_size = 0;///< size for full coefficients arrays (depends on B- or C- basis)

    NS_TYPE *full_ns = nullptr; ///< ns contiguous package [rank_array_total_size]
    LS_TYPE *full_ls = nullptr; ///<  ls contiguous package [rank_array_total_size]
    SPECIES_TYPE *full_mus = nullptr; ///<  mus contiguous package [rank_array_total_size]
    MS_TYPE *full_ms = nullptr; ///< //m_s contiguous package [ms_array_total_size]

    /**
     * Rearrange basis functions in contiguous memory to optimize cache access
     */
    virtual void pack_flatten_basis() = 0;

    virtual void flatten_basis() = 0;

    //1D flat array basis representation: [mu]
    SHORT_INT_TYPE *total_basis_size_rank1 = nullptr; ///< per-species type array of total_basis_rank1[mu] sizes
    SHORT_INT_TYPE *total_basis_size = nullptr; ///< per-species type array of total_basis[mu] sizes

    size_t max_B_array_size = 0; ///< maximum over elements array size for B[func_ind][ms_ind]
    size_t max_dB_array_size = 0; ///<  maximum over elements array size for dB[func_ind][ms_ind][r]

    SHORT_INT_TYPE num_ms_combinations_max = 0; ///< maximum number of ms combinations  among all basis functions

    /**
     * Default constructor
     */
    ACEFlattenBasisSet() = default;

    /**
     * Copy constructor (see Rule of Three)
     * @param other
     */
    ACEFlattenBasisSet(const ACEFlattenBasisSet &other);

    /**
     * operator= (see Rule of Three)
     * @param other
     * @return
     */
    ACEFlattenBasisSet &operator=(const ACEFlattenBasisSet &other);

    /**
     * Destructor  (see Rule of Three)
     */
    ~ACEFlattenBasisSet() override;

    // routines for copying and cleaning dynamic memory of the class (see Rule of Three)
    /**
     * Cleaning dynamic memory of the class (see. Rule of Three),
     * must be idempotent for safety
     */
    void _clean() override;

    /**
    * Copying scalar variables
    * @param src
    */
    void _copy_scalar_memory(const ACEFlattenBasisSet &src);

    /**
    * Copying and cleaning dynamic memory of the class (see. Rule of Three)
    * @param src
    */
    void _copy_dynamic_memory(const ACEFlattenBasisSet &src);

    /**
    * Clean contiguous arrays
    */
    virtual void _clean_contiguous_arrays();

    /**
     * Release dynamic arrays with basis set sizes
     */
    void _clean_basissize_arrays();
};

#endif //ACE_EVALUATOR_ACE_FLATTEN_BASIS_H
