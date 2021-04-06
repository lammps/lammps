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

// Created by Yury Lysogorskiy on 26.02.20.

#ifndef ACE_C_BASISFUNCTION_H
#define ACE_C_BASISFUNCTION_H

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "ace_types.h"

//macros for copying the member-array from "other" object for C-tilde and B-basis
#define basis_mem_copy(other, array, size, type) if(other.array) { \
    if(!is_proxy) delete[] array;\
    array = new type[(size)];\
    is_proxy = false;\
    memcpy(array, other.array, (size) * sizeof(type));\
}

/**
 * Abstract basis function, that stores general quantities
 */
struct ACEAbstractBasisFunction {
    /**
    * flattened array of  computed combinations of (m1, m2, ..., m_rank)
    * which have non-zero general Clebsch-Gordan coefficient:
    *  \f$ \mathbf{m}_1, \dots, \mathbf{m}_\mathrm{num ms combs}\f$ =
    * \f$ (m_1, m_2, \dots, m_{rank})_1,  \dots, (m_1, m_2, \dots, m_{rank})_{\mathrm{num ms combs}} \f$,
    * size =  num_ms_combs * rank,
    * effective shape: [num_ms_combs][rank]
    */
    MS_TYPE *ms_combs = nullptr;

    /**
    * species types of neighbours atoms \f$ \mathbf{\mu} = (\mu_1, \mu_2, \dots, \mu_{rank}) \f$,
    * should be lexicographically sorted,
    * size  = rank,
    * effective shape: [rank]
    */
    SPECIES_TYPE *mus = nullptr;

    /**
    * indices for radial part \f$ \mathbf{n} = (n_1, n_2, \dots, n_{rank}) \f$,
    * should be lexicographically sorted,
    * size  = rank,
    * effective shape: [rank]
    */
    NS_TYPE *ns = nullptr;


    /**
    * indices for radial part \f$ \mathbf{l} = (l_1, l_2, \dots, l_{rank}) \f$,
    * should be lexicographically sorted,
    * size  = rank,
    * effective shape: [rank]
    */
    LS_TYPE *ls = nullptr;

    SHORT_INT_TYPE num_ms_combs = 0; ///< number of different ms-combinations

    RANK_TYPE rank = 0; ///< number of atomic base functions "A"s in basis function product B

    DENSITY_TYPE ndensity = 0; ///< number of densities

    SPECIES_TYPE mu0 = 0; ///< species type of central atom

    /**
    * whether ms array contains only "non-negative" ms-combinations.
    * positive ms-combination is when the first non-zero m is positive (0 1 -1)
    * negative ms-combination is when the first non-zero m is negative (0 -1 1)
    */
    bool is_half_ms_basis = false;

    /*
     * flag, whether object is "owner" (i.e. responsible for memory cleaning) of
    * the ms, ns, ls, mus and other dynamically allocated arrases or just a proxy for them
     */
    bool is_proxy = false;

    /**
     * Copying static and dynamic memory variables from another ACEAbstractBasisFunction.
     * Always copy the dynamic memory, even if the source is a proxy object
     * @param other
     */
    virtual void _copy_from(const ACEAbstractBasisFunction &other) {
        rank = other.rank;
        ndensity = other.ndensity;
        mu0 = other.mu0;
        num_ms_combs = other.num_ms_combs;
        is_half_ms_basis = other.is_half_ms_basis;
        is_proxy = false;

        basis_mem_copy(other, mus, rank, SPECIES_TYPE)
        basis_mem_copy(other, ns, rank, NS_TYPE)
        basis_mem_copy(other, ls, rank, LS_TYPE)
        basis_mem_copy(other, ms_combs, num_ms_combs * rank, MS_TYPE)
    }

    /**
     * Clean the dynamically allocated memory if object is responsible for it
     */
    virtual void _clean() {
        //release memory if the structure is not proxy
        if (!is_proxy) {
            delete[] mus;
            delete[] ns;
            delete[] ls;
            delete[] ms_combs;
        }

        mus = nullptr;
        ns = nullptr;
        ls = nullptr;
        ms_combs = nullptr;
    }

};

/**
 * Representation of C-tilde basis function, i.e. the function that is summed up over a group of B-functions
 * that differs only by intermediate coupling orbital moment \f$ L \f$ and coefficients.
 */
struct ACECTildeBasisFunction : public ACEAbstractBasisFunction {

    /**
     * c_tilde coefficients for all densities,
     * size = num_ms_combs*ndensity,
     * effective shape [num_ms_combs][ndensity]
     */
    DOUBLE_TYPE *ctildes = nullptr;

    /*
     * Default constructor
     */
    ACECTildeBasisFunction() = default;

    // Because the ACECTildeBasisFunction contains dynamically allocated arrays, the Rule of Three should be
    // fulfilled, i.e. copy constructor (copy the dynamic arrays), operator= (release previous arrays and
    // copy the new dynamic arrays) and destructor (release all dynamically allocated memory)

    /**
     * Copy constructor, to fulfill the Rule of Three.
     * Always copy the dynamic memory, even if the source is a proxy object.
     */
    ACECTildeBasisFunction(const ACECTildeBasisFunction &other) {
        _copy_from(other);
    }

    /*
     * operator=, to fulfill the Rule of Three.
     * Always copy the dynamic memory, even if the source is a proxy object
     */
    ACECTildeBasisFunction &operator=(const ACECTildeBasisFunction &other) {
        _clean();
        _copy_from(other);
        return *this;
    }

    /*
     * Destructor
     */
    ~ACECTildeBasisFunction() {
        _clean();
    }

    /**
     * Copy from another object, always copy the memory, even if the source is a proxy object
     * @param other
     */
    void _copy_from(const ACECTildeBasisFunction &other) {
        ACEAbstractBasisFunction::_copy_from(other);
        is_proxy = false;
        basis_mem_copy(other, ctildes, num_ms_combs * ndensity, DOUBLE_TYPE)
    }

    /**
     * Clean the dynamically allocated memory if object is responsible for it
     */
    void _clean() override {
        ACEAbstractBasisFunction::_clean();
        //release memory if the structure is not proxy
        if (!is_proxy) {
            delete[] ctildes;
        }
        ctildes = nullptr;
    }

    /**
     * Print the information about basis function to cout, in order to ease the output redirection.
     */
    void print() const {
        using namespace std;
        cout << "ACECTildeBasisFunction: ndensity= " << (int) this->ndensity << ", mu0 = " << (int) this->mu0 << " ";
        cout << " mus=(";
        for (RANK_TYPE r = 0; r < this->rank; r++)
            cout << (int) this->mus[r] << " ";
        cout << "), ns=(";
        for (RANK_TYPE r = 0; r < this->rank; r++)
            cout << this->ns[r] << " ";
        cout << "), ls=(";
        for (RANK_TYPE r = 0; r < this->rank; r++)
            cout << this->ls[r] << " ";

        cout << "), " << this->num_ms_combs << " m_s combinations: {" << endl;

        for (int i = 0; i < this->num_ms_combs; i++) {
            cout << "\t< ";
            for (RANK_TYPE r = 0; r < this->rank; r++)
                cout << this->ms_combs[i * this->rank + r] << " ";
            cout << " >: c_tilde=";
            for (DENSITY_TYPE p = 0; p < this->ndensity; ++p)
                cout << " " << this->ctildes[i * this->ndensity + p] << " ";
            cout << endl;
        }
        if (this->is_proxy)
            cout << "proxy ";
        if (this->is_half_ms_basis)
            cout << "half_ms_basis";
        cout << "}" << endl;
    }
};

#endif //ACE_C_BASISFUNCTION_H
