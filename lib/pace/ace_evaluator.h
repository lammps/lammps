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


// Created by Yury Lysogorskiy on 31.01.20.

#ifndef ACE_EVALUATOR_H
#define ACE_EVALUATOR_H

#include "ace_abstract_basis.h"
#include "ace_arraynd.h"
#include "ace_array2dlm.h"
#include "ace_c_basis.h"
#include "ace_complex.h"
#include "ace_timing.h"
#include "ace_types.h"

/**
 * Basic evaluator class, that should accept the basis set and implement the "compute_atom" method using given basis set.
 */
class ACEEvaluator {
protected:

    Array2D<DOUBLE_TYPE> A_rank1 = Array2D<DOUBLE_TYPE>("A_rank1"); ///< 2D-array for storing A's for rank=1, shape: A(mu_j,n)
    Array4DLM<ACEComplex> A = Array4DLM<ACEComplex>("A"); ///< 4D array with (l,m) last indices  for storing A's for rank>1: A(mu_j, n, l, m)

    Array1D<DOUBLE_TYPE> rhos = Array1D<DOUBLE_TYPE>("rhos"); ///< densities \f$ \rho^{(p)} \f$(ndensity), p  = 0 .. ndensity-1
    Array1D<DOUBLE_TYPE> dF_drho = Array1D<DOUBLE_TYPE>("dF_drho"); ///< derivatives of cluster functional wrt. densities, index = 0 .. ndensity-1

    /**
     * Initialize internal arrays according to basis set sizes
     * @param basis_set
     */
    void init(ACEAbstractBasisSet *basis_set);

public:
    // set of timers for code profiling

    ACETimer loop_over_neighbour_timer; ///< timer for loop over neighbours when constructing A's for single central atom
    ACETimer per_atom_calc_timer; ///< timer for single compute_atom call


    ACETimer forces_calc_loop_timer; ///< timer for forces calculations for single central atom
    ACETimer forces_calc_neighbour_timer; ///< timer for loop over neighbour atoms for force calculations

    ACETimer energy_calc_timer; ///< timer for energy calculation
    ACETimer total_time_calc_timer; ///< timer for total calculations of all atoms within given atomic environment system

    /**
     * Initialize all timers
     */
    void init_timers();

    /**
     * Mapping from external atoms types, i.e. LAMMPS, to internal SPECIES_TYPE, used in basis functions
     */
    Array1D<int> element_type_mapping = Array1D<int>("element_type_mapping");


    DOUBLE_TYPE e_atom = 0; ///< energy of current atom, including core-repulsion

    /**
     * temporary array for the pair forces between current atom_i and its neighbours atom_k
     * neighbours_forces(k,3),  k = 0..num_of_neighbours(atom_i)-1
     */
    Array2D<DOUBLE_TYPE> neighbours_forces = Array2D<DOUBLE_TYPE>("neighbours_forces");

    ACEEvaluator() = default;

    virtual ~ACEEvaluator() = default;

     /**
      * The key method to compute energy and forces for atom 'i'.
      * Method will update the  "e_atom" variable and "neighbours_forces(jj, alpha)" array
      *
      * @param i atom index
      * @param x atomic positions array of the real and ghost atoms, shape: [atom_ind][3]
      * @param type  atomic types array of the real and ghost atoms, shape: [atom_ind]
      * @param jnum  number of neighbours of atom_i
      * @param jlist array of neighbour indices, shape: [jnum]
      */
    virtual void compute_atom(int i, DOUBLE_TYPE **x, const SPECIES_TYPE *type, const int jnum, const int *jlist) = 0;

    /**
     * Resize all caches over neighbours atoms
     * @param max_jnum  maximum number of neighbours
     */
    virtual void resize_neighbours_cache(int max_jnum) = 0;

#ifdef EXTRA_C_PROJECTIONS
    /**
     * 2D array to store projections of basis function for rank = 1, shape: [func_ind][ndensity]
     */
    Array2D<DOUBLE_TYPE> basis_projections_rank1 = Array2D<DOUBLE_TYPE>("basis_projections_rank1");

    /**
     * 2D array to store projections of basis function for rank > 1, shape: [func_ind][ndensity]
     */
    Array2D<DOUBLE_TYPE> basis_projections = Array2D<DOUBLE_TYPE>("basis_projections");
#endif
};

//TODO: split into separate file
/**
 * Evaluator for C-tilde basis set, that should accept the basis set and implement the "compute_atom" method using given basis set.
 */
class ACECTildeEvaluator : public ACEEvaluator {

    /**
     * Weights \f$ \omega_{i \mu n 0 0} \f$ for rank = 1, see Eq.(10) from implementation notes,
     * 'i' is fixed for the current atom, shape: [nelements][nradbase]
     */
    Array2D<DOUBLE_TYPE> weights_rank1 = Array2D<DOUBLE_TYPE>("weights_rank1");

    /**
     * Weights \f$ \omega_{i \mu n l m} \f$ for rank > 1, see Eq.(10) from implementation notes,
     * 'i' is fixed for the current atom, shape: [nelements][nradbase][l=0..lmax, m]
     */
    Array4DLM<ACEComplex> weights = Array4DLM<ACEComplex>("weights");

    /**
     * cache for gradients of \f$ g(r)\f$: grad_phi(jj,n)=A2DLM(l,m)
     * shape:[max_jnum][nradbase]
     */
    Array2D<DOUBLE_TYPE> DG_cache = Array2D<DOUBLE_TYPE>("DG_cache");


    /**
     * cache for \f$ R_{nl}(r)\f$
     * shape:[max_jnum][nradbase][0..lmax]
     */
    Array3D<DOUBLE_TYPE> R_cache = Array3D<DOUBLE_TYPE>("R_cache");
    /**
     * cache for derivatives of \f$ R_{nl}(r)\f$
     * shape:[max_jnum][nradbase][0..lmax]
     */
    Array3D<DOUBLE_TYPE> DR_cache = Array3D<DOUBLE_TYPE>("DR_cache");
    /**
     * cache for \f$ Y_{lm}(\hat{r})\f$
     * shape:[max_jnum][0..lmax][m]
     */
    Array3DLM<ACEComplex> Y_cache = Array3DLM<ACEComplex>("Y_cache");
    /**
     * cache for \f$ \nabla Y_{lm}(\hat{r})\f$
     * shape:[max_jnum][0..lmax][m]
     */
    Array3DLM<ACEDYcomponent> DY_cache = Array3DLM<ACEDYcomponent>("dY_dense_cache");

    /**
     * cache for derivatives of hard-core repulsion
     * shape:[max_jnum]
     */
    Array1D<DOUBLE_TYPE> DCR_cache = Array1D<DOUBLE_TYPE>("DCR_cache");

    /**
    * Partial derivatives \f$ dB_{i \mu n l m t}^{(r)} \f$  with sequential numbering over [func_ind][ms_ind][r],
    * shape:[func_ms_r_ind]
    */
    Array1D<ACEComplex> dB_flatten = Array1D<ACEComplex>("dB_flatten");

    /**
     * pointer to the ACEBasisSet object
     */
    ACECTildeBasisSet *basis_set = nullptr;

    /**
     * Initialize internal arrays according to basis set sizes
     * @param basis_set
     */
    void init(ACECTildeBasisSet *basis_set);

public:

    ACECTildeEvaluator() = default;

    explicit ACECTildeEvaluator(ACECTildeBasisSet &bas) {
        set_basis(bas);
    }

    /**
     * set the basis function to the ACE evaluator
     * @param bas
     */
    void set_basis(ACECTildeBasisSet &bas);

    /**
     * The key method to compute energy and forces for atom 'i'.
     * Method will update the  "e_atom" variable and "neighbours_forces(jj, alpha)" array
     *
     * @param i atom index
     * @param x atomic positions array of the real and ghost atoms, shape: [atom_ind][3]
     * @param type  atomic types array of the real and ghost atoms, shape: [atom_ind]
     * @param jnum  number of neighbours of atom_i
     * @param jlist array of neighbour indices, shape: [jnum]
     */
    void compute_atom(int i, DOUBLE_TYPE **x, const SPECIES_TYPE *type, const int jnum, const int *jlist) override;

    /**
     * Resize all caches over neighbours atoms
     * @param max_jnum  maximum number of neighbours
     */
    void resize_neighbours_cache(int max_jnum) override;
};


#endif //ACE_EVALUATOR_H