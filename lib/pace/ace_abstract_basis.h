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


// Created by Lysogorskiy Yury on 28.04.2020.

#ifndef ACE_EVALUATOR_ACE_ABSTRACT_BASIS_H
#define ACE_EVALUATOR_ACE_ABSTRACT_BASIS_H

#include <vector>
#include <string>

#include "ace_c_basisfunction.h"
#include "ace_contigous_array.h"
#include "ace_radial.h"
#include "ace_spherical_cart.h"
#include "ace_types.h"

using namespace std;

/**
 * Abstract basis set class
 */
class ACEAbstractBasisSet {
public:
    SPECIES_TYPE nelements = 0;        ///< number of elements in basis set
    RANK_TYPE rankmax = 0;             ///< maximum value of rank
    DENSITY_TYPE ndensitymax = 0;      ///< maximum number of densities \f$ \rho^{(p)} \f$
    NS_TYPE nradbase = 0; ///< maximum number of radial \f$\textbf{basis}\f$ function \f$ g_{k}(r) \f$
    LS_TYPE lmax = 0;  ///< \f$ l_\textrm{max} \f$ - maximum value of orbital moment \f$ l \f$
    NS_TYPE nradmax = 0;  ///< maximum number \f$ n \f$ of radial function \f$ R_{nl}(r) \f$
    DOUBLE_TYPE cutoffmax = 0;  ///< maximum value of cutoff distance among all species in basis set
    DOUBLE_TYPE deltaSplineBins = 0;  ///< Spline interpolation density

    string npoti = "FinnisSinclair"; ///< FS and embedding function combination

    string *elements_name = nullptr; ///< Array of elements name for mapping from index (0..nelements-1) to element symbol (string)

    AbstractRadialBasis *radial_functions = nullptr; ///< object to work with radial functions
    ACECartesianSphericalHarmonics spherical_harmonics; ///< object to work with spherical harmonics in Cartesian representation


    Array1D<DOUBLE_TYPE> rho_core_cutoffs; ///< energy-based inner cut-off
    Array1D<DOUBLE_TYPE> drho_core_cutoffs; ///< decay of energy-based inner cut-off

    vector<DOUBLE_TYPE> FS_parameters;  ///< parameters for cluster functional, see Eq.(3) in implementation notes or Eq.(53) in <A HREF="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.99.014104">  PRB 99, 014104 (2019) </A>

    // E0 values
    Array1D<DOUBLE_TYPE> E0vals;

    /**
     * Default empty constructor
     */
    ACEAbstractBasisSet() = default;

    // copy constructor, operator= and destructor (see. Rule of Three)

    /**
     * Copy constructor (see. Rule of Three)
     * @param other
     */
    ACEAbstractBasisSet(const ACEAbstractBasisSet &other);

    /**
     * operator=  (see. Rule of Three)
     * @param other
     * @return
     */
    ACEAbstractBasisSet &operator=(const ACEAbstractBasisSet &other);

    /**
     * virtual destructor (see. Rule of Three)
     */
    virtual ~ACEAbstractBasisSet();

    /**
     * Computing cluster functional \f$ F(\rho_i^{(1)}, \dots, \rho_i^{(P)})  \f$
     * and its derivatives  \f$ (\partial F/\partial\rho_i^{(1)}, \dots, \partial F/\partial \rho_i^{(P)} ) \f$
     * @param rhos array with densities \f$ \rho^{(p)} \f$
     * @param value (out) return value of cluster functional
     * @param derivatives (out) array of derivatives  \f$ (\partial F/\partial\rho_i^{(1)}, \dots, \partial F/\partial \rho_i^{(P)} )  \f$
     * @param ndensity  number \f$ P \f$ of densities to use
     */
    void FS_values_and_derivatives(Array1D<DOUBLE_TYPE> &rhos, DOUBLE_TYPE &value, Array1D<DOUBLE_TYPE> &derivatives,
                                   DENSITY_TYPE ndensity);

    /**
     * Computing hard core pairwise repulsive potential \f$ f_{cut}(\rho_i^{(\textrm{core})})\f$ and its derivative,
     * see Eq.(29) of implementation notes
     * @param rho_core value of \f$ \rho_i^{(\textrm{core})} \f$
     * @param rho_cut  \f$ \rho_{cut}^{\mu_i} \f$ value
     * @param drho_cut \f$ \Delta_{cut}^{\mu_i} \f$ value
     * @param fcut (out) return inner cutoff function
     * @param dfcut (out) return derivative of inner cutoff function
     */
    static void inner_cutoff(DOUBLE_TYPE rho_core, DOUBLE_TYPE rho_cut, DOUBLE_TYPE drho_cut, DOUBLE_TYPE &fcut,
                             DOUBLE_TYPE &dfcut);


    /**
     * Virtual method to save potential to file
     * @param filename file name
     */
    virtual void save(const string &filename) = 0;

    /**
     * Virtual method to load potential from file
     * @param filename file name
     */
    virtual void load(const string filename) = 0;

    /**
     * Get the species index by its element name
     * @param elemname element name
     * @return species index
     */
    SPECIES_TYPE get_species_index_by_name(const string &elemname);


    // routines for copying and cleaning dynamic memory of the class (see. Rule of Three)

    /**
     * Routine for clean the dynamically allocated memory\n
     * IMPORTANT! It must be idempotent for safety.
     */
    virtual void _clean();

    /**
     * Copy dynamic memory from src. Must be override and extended in derived classes!
     * @param src source object to copy from
     */
    virtual void _copy_dynamic_memory(const ACEAbstractBasisSet &src);

    /**
     * Copy scalar values from src. Must be override and extended in derived classes!
     * @param src source object to copy from
     */
    virtual void _copy_scalar_memory(const ACEAbstractBasisSet &src);
};

void Fexp(DOUBLE_TYPE rho, DOUBLE_TYPE mexp, DOUBLE_TYPE &F, DOUBLE_TYPE &DF);

void FexpShiftedScaled(DOUBLE_TYPE rho, DOUBLE_TYPE mexp, DOUBLE_TYPE &F, DOUBLE_TYPE &DF);

#endif //ACE_EVALUATOR_ACE_ABSTRACT_BASIS_H
