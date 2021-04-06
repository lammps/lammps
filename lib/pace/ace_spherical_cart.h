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

// Created by  Ralf Drautz, Yury Lysogorskiy

#ifndef ACE_SPHERICAL_CART_H
#define ACE_SPHERICAL_CART_H

#include <cmath>

#include "ace_arraynd.h"
#include "ace_array2dlm.h"
#include "ace_complex.h"
#include "ace_types.h"


using namespace std;

const DOUBLE_TYPE sq1o4pi = 0.28209479177387814347; // sqrt(1/(4*pi))
const DOUBLE_TYPE sq4pi = 3.54490770181103176384; // sqrt(4*pi)
const DOUBLE_TYPE sq3 = 1.73205080756887719318;//sqrt(3), numpy
const DOUBLE_TYPE sq3o2 = 1.22474487139158894067;//sqrt(3/2), numpy

//definition of common factor for spherical harmonics = Y00
//const DOUBLE_TYPE Y00 = sq1o4pi;
const DOUBLE_TYPE Y00 = 1;

/**
Class to store spherical harmonics and their associated functions. \n
All the associated members such as \f$ P_{lm}, Y_{lm}\f$ etc are one dimensional arrays of length (L+1)*(L+2)/2. \n
The value that corresponds to a particular l, m configuration can be accessed through a \code ylm(l,m) \endcode \n
*/
class ACECartesianSphericalHarmonics {
public:

    /**
    int, the number of spherical harmonics to be found
    */
    LS_TYPE lmax;

    /**
     * Default constructor
     */
    ACECartesianSphericalHarmonics() = default;

    /**
     * Parametrized constructor.  Dynamically initialises all the arrays.
     * @param lmax maximum orbital moment
     */
    explicit ACECartesianSphericalHarmonics(LS_TYPE lmax);

    /**
     * Initialize internal arrays and precompute necessary coefficients
     * @param lm maximum orbital moment
     */
    void init(LS_TYPE lm);

    /**
     * Destructor
     */
    ~ACECartesianSphericalHarmonics();

    /**
     * Precompute necessaary helper arrays Precomputes the value of \f$ a_{lm}, b_{lm}, c_l, d_l \f$
     */
    void pre_compute();

    /**
    Function that computes \f$ \bar{P}_{lm} \f$ for the corresponding lmax value
    Input is \f$ \hat{r}_z \f$ which is the $z$-component of the bond direction.

    For each \f$ \hat{r}_z \f$, this computes the whole range of \f$ \bar{P}_{lm} \f$ values
    and its derivatives upto the lmax specified, which is a member of the class.

    @param rz, DOUBLE_TYPE

    @returns None
    */
    void compute_barplm(DOUBLE_TYPE rz, LS_TYPE lmaxi);

    /**
    Function that computes \f$ Y_{lm} \f$ for the corresponding lmax value
    Input is the bond-directon vector \f$ \hat{r}_x, \hat{r}_y, \hat{r}_z \f$

    Each \f$ Y_{lm} \f$ value is a ACEComplex object with real and imaginary parts. This function also
    finds the derivatives, which are stored in the Dycomponent class, with each component being a
    ACEComplex object.

    @param rx, DOUBLE_TYPE
    @param ry, DOUBLE_TYPE
    @param rz, DOUBLE_TYPE
    @param lmaxi, int
    */
    void compute_ylm(DOUBLE_TYPE rx, DOUBLE_TYPE ry, DOUBLE_TYPE rz, LS_TYPE lmaxi);

    Array2DLM<DOUBLE_TYPE> alm;
    Array2DLM<DOUBLE_TYPE> blm;
    Array1D<DOUBLE_TYPE> cl;
    Array1D<DOUBLE_TYPE> dl;

    Array2DLM<DOUBLE_TYPE> plm;
    Array2DLM<DOUBLE_TYPE> dplm;

    Array2DLM<ACEComplex> ylm; ///< Values of all spherical harmonics after \code compute_ylm(rx,ry,rz, lmaxi) \endcode call
    Array2DLM<ACEDYcomponent> dylm;///< Values of gradients of all spherical harmonics after \code compute_ylm(rx,ry,rz, lmaxi) \endcode call

};


#endif
