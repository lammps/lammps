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


// Created by Christoph Ortner on 03.06.2020

#ifndef SHIPs_RADIAL_FUNCTIONS_H
#define SHIPs_RADIAL_FUNCTIONS_H

#include "ace_arraynd.h"
#include "ace_types.h"
#include "ace_radial.h"

class SHIPsRadPolyBasis {

public:

    // transform parameters
    int p = 0;
    DOUBLE_TYPE r0 = 0.0;

    // cutoff parameters
    DOUBLE_TYPE rcut = 0.0;
    DOUBLE_TYPE xl = 0.0;
    DOUBLE_TYPE xr = 0.0;
    int pl = 0;
    int pr = 0;

    // basis size
    size_t maxn = 0;

    // recursion parameters
    Array1D<DOUBLE_TYPE> A = Array1D<DOUBLE_TYPE>("SHIPs radial basis: A");
    Array1D<DOUBLE_TYPE> B = Array1D<DOUBLE_TYPE>("SHIPs radial basis: B");
    Array1D<DOUBLE_TYPE> C = Array1D<DOUBLE_TYPE>("SHIPs radial basis: C");

    // temporary storage for evaluating the basis
    Array1D<DOUBLE_TYPE> P = Array1D<DOUBLE_TYPE>("SHIPs radial basis: P");
    Array1D<DOUBLE_TYPE> dP_dr = Array1D<DOUBLE_TYPE>("SHIPs radial basis: dP");

//////////////////////////////////

    SHIPsRadPolyBasis() = default;

    ~SHIPsRadPolyBasis() = default;

    // distance transform
    void transform(const DOUBLE_TYPE r, DOUBLE_TYPE &x_out, DOUBLE_TYPE &dx_out) const;

    // cutoff function
    void fcut(const DOUBLE_TYPE x, DOUBLE_TYPE &f_out, DOUBLE_TYPE &df_out) const;

    void fread(FILE *fptr);

    void _init(DOUBLE_TYPE r0, int p, DOUBLE_TYPE rcut,
               DOUBLE_TYPE xl, DOUBLE_TYPE xr,
               int pl, int pr, size_t maxn);

    void calcP(DOUBLE_TYPE r, size_t maxn, SPECIES_TYPE z1, SPECIES_TYPE z2);

    size_t get_maxn();

};




class SHIPsRadialFunctions : public AbstractRadialBasis {
public:

    // radial basis 
    SHIPsRadPolyBasis radbasis; 

    // pair potential basis 
    bool haspair = false; 
    SHIPsRadPolyBasis pairbasis; 

    // pair potential coefficients 
    Array1D<DOUBLE_TYPE> paircoeffs = Array1D<DOUBLE_TYPE>("SHIPs pairpot coeffs: paircoeffs");

    // spline parameters for repulsive core
    DOUBLE_TYPE ri = 0.0;
    DOUBLE_TYPE e0 = 0.0;
    DOUBLE_TYPE A = 0.0;
    DOUBLE_TYPE B = 0.0;

//////////////////////////////////

    SHIPsRadialFunctions() = default;

    ~SHIPsRadialFunctions() override = default;


    void fread(FILE *fptr);

    void load(string fname);

    size_t get_maxn();
    DOUBLE_TYPE get_rcut(); 

    bool has_pair(); 

    void init(NS_TYPE nradb, LS_TYPE lmax, NS_TYPE nradial, DOUBLE_TYPE deltaSplineBins, SPECIES_TYPE nelements,
              DOUBLE_TYPE cutoff,
              string radbasename) override;

    void
    evaluate(DOUBLE_TYPE r, NS_TYPE nradbase_c, NS_TYPE nradial_c, SPECIES_TYPE mu_i, SPECIES_TYPE mu_j,
             bool calc_second_derivatives = false) override;

    void
    evaluate_pair(DOUBLE_TYPE r, SPECIES_TYPE mu_i, SPECIES_TYPE mu_j,
                  bool calc_second_derivatives = false);

    void setuplookupRadspline() override;

    SHIPsRadialFunctions *clone() const override {
        return new SHIPsRadialFunctions(*this);
    };

    /**
     * Helper method, that populate `fr` and `dfr` 2D-arrays (n,l) with P(n), dP_dr  for given coordinate r
     * @param r
     * @param maxn
     * @param z1
     * @param z2
     */
    void fill_Rnl(DOUBLE_TYPE r, NS_TYPE maxn, SPECIES_TYPE z1, SPECIES_TYPE z2);

    void fill_gk(DOUBLE_TYPE r, NS_TYPE maxn, SPECIES_TYPE z1, SPECIES_TYPE z2);
};


#endif
