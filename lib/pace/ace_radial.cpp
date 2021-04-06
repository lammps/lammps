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


#include <cmath>
#include <functional>
#include <stdexcept>

#include "ace_radial.h"

const DOUBLE_TYPE pi = 3.14159265358979323846264338327950288419; // pi

ACERadialFunctions::ACERadialFunctions(NS_TYPE nradb, LS_TYPE lmax, NS_TYPE nradial, DOUBLE_TYPE deltaSplineBins,
                                       SPECIES_TYPE nelements,
                                       DOUBLE_TYPE cutoff, string radbasename) {
    init(nradb, lmax, nradial, deltaSplineBins, nelements, cutoff, radbasename);
}

void ACERadialFunctions::init(NS_TYPE nradb, LS_TYPE lmax, NS_TYPE nradial, DOUBLE_TYPE deltaSplineBins,
                              SPECIES_TYPE nelements,
                              DOUBLE_TYPE cutoff, string radbasename) {
    this->nradbase = nradb;
    this->lmax = lmax;
    this->nradial = nradial;
    this->deltaSplineBins = deltaSplineBins;
    this->nelements = nelements;
    this->cutoff = cutoff;
    this->radbasename = radbasename;

    gr.init(nradbase, "gr");
    dgr.init(nradbase, "dgr");
    d2gr.init(nradbase, "d2gr");


    fr.init(nradial, lmax + 1, "fr");
    dfr.init(nradial, lmax + 1, "dfr");
    d2fr.init(nradial, lmax + 1, "d2fr");


    cheb.init(nradbase + 1, "cheb");
    dcheb.init(nradbase + 1, "dcheb");
    cheb2.init(nradbase + 1, "cheb2");


    splines_gk.init(nelements, nelements, "splines_gk");
    splines_rnl.init(nelements, nelements, "splines_rnl");
    splines_hc.init(nelements, nelements, "splines_hc");

    lambda.init(nelements, nelements, "lambda");
    lambda.fill(1.);

    cut.init(nelements, nelements, "cut");
    cut.fill(1.);

    dcut.init(nelements, nelements, "dcut");
    dcut.fill(1.);

    crad.init(nelements, nelements, nradial, (lmax + 1), nradbase, "crad");
    crad.fill(0.);

    //hard-core repulsion
    prehc.init(nelements, nelements, "prehc");
    prehc.fill(0.);

    lambdahc.init(nelements, nelements, "lambdahc");
    lambdahc.fill(1.);

}


/**
Function that computes Chebyshev polynomials of first and second kind
 to setup the radial functions and the derivatives

@param n, x

@returns cheb1, dcheb1
*/
void ACERadialFunctions::calcCheb(NS_TYPE n, DOUBLE_TYPE x) {
    if (n < 0) {
        char s[1024];
        sprintf(s, "The order n of the polynomials should be positive %d\n", n);
        throw std::invalid_argument(s);
    }
    DOUBLE_TYPE twox = 2.0 * x;
    cheb(0) = 1.;
    dcheb(0) = 0.;
    cheb2(0) = 1.;

    if (nradbase > 1) {
        cheb(1) = x;
        cheb2(1) = twox;
    }
    for (NS_TYPE m = 1; m <= n - 1; m++) {
        cheb(m + 1) = twox * cheb(m) - cheb(m - 1);
        cheb2(m + 1) = twox * cheb2(m) - cheb2(m - 1);
    }
    for (NS_TYPE m = 1; m <= n; m++) {
        dcheb(m) = m * cheb2(m - 1);
    }
#ifdef DEBUG_RADIAL
    for ( NS_TYPE  m=0; m<=n; m++ ) {
        printf(" m %d cheb %f dcheb %f \n", m, cheb(m), dcheb(m));
    }
#endif
}

/**
Function that computes radial basis.

@param lam, nradbase, cut, dcut, r

@returns gr, dgr
*/
void ACERadialFunctions::radbase(DOUBLE_TYPE lam, DOUBLE_TYPE cut, DOUBLE_TYPE dcut, DOUBLE_TYPE r) {
    /*lam is given by the formula (24), that contains cut */

    if (r < cut) {
        if (radbasename == "ChebExpCos") {
            chebExpCos(lam, cut, dcut, r);
        } else if (radbasename == "ChebPow") {
            chebPow(lam, cut, dcut, r);
        } else if (radbasename == "ChebLinear") {
            chebLinear(lam, cut, dcut, r);
        } else {
            throw invalid_argument("Unknown radial basis function name: " + radbasename);
        }
    } else {
        gr.fill(0);
        dgr.fill(0);
    }
}

/***
 *  Radial function: ChebExpCos, cheb exp scaling including cos envelope
 * @param lam function parameter
 * @param cut cutoff distance
 * @param r function input argument
 * @return fills in gr and dgr arrays
 */
void
ACERadialFunctions::chebExpCos(DOUBLE_TYPE lam, DOUBLE_TYPE cut, DOUBLE_TYPE dcut, DOUBLE_TYPE r) {
    DOUBLE_TYPE y2, y1, x, dx;
    DOUBLE_TYPE env, denv, fcut, dfcut;
    /* scaled distance x and derivative*/
    y1 = exp(-lam * r / cut);
    y2 = exp(-lam);
    x = 1.0 - 2.0 * ((y1 - y2) / (1 - y2));
    dx = 2 * (lam / cut) * (y1 / (1 - y2));
    /* calculation of Chebyshev polynomials from the recursion */
    calcCheb(nradbase - 1, x);
    gr(0) = cheb(0);
    dgr(0) = dcheb(0) * dx;
    for (NS_TYPE n = 2; n <= nradbase; n++) {
        gr(n - 1) = 0.5 - 0.5 * cheb(n - 1);
        dgr(n - 1) = -0.5 * dcheb(n - 1) * dx;
    }
    env = 0.5 * (1.0 + cos(M_PI * r / cut));
    denv = -0.5 * sin(M_PI * r / cut) * M_PI / cut;
    for (NS_TYPE n = 0; n < nradbase; n++) {
        dgr(n) = gr(n) * denv + dgr(n) * env;
        gr(n) = gr(n) * env;
    }
    // for radtype = 3 a smooth cut is already included in the basis function
    dx = cut - dcut;
    if (r > dx) {
        fcut = 0.5 * (1.0 + cos(M_PI * (r - dx) / dcut));
        dfcut = -0.5 * sin(M_PI * (r - dx) / dcut) * M_PI / dcut;
        for (NS_TYPE n = 0; n < nradbase; n++) {
            dgr(n) = gr(n) * dfcut + dgr(n) * fcut;
            gr(n) = gr(n) * fcut;
        }
    }
}

/***
*  Radial function: ChebPow, Radial function: ChebPow
* - argument of Chebyshev polynomials
* x = 2.0*( 1.0 - (1.0 - r/rcut)^lam ) - 1.0
* - radial function
* gr(n) = ( 1.0 - Cheb(n) )/2.0, n = 1,...,nradbase
* - the function fulfills:
* gr(n) = 0 at rcut
* dgr(n) = 0 at rcut for lam >= 1
* second derivative zero at rcut for lam >= 2
* -> the radial function does not require a separate cutoff function
* - corresponds to radial basis radtype=5 in Fortran code
*
* @param lam function parameter
* @param cut cutoff distance
* @param r function input argument
* @return fills in gr and dgr arrays
*/
void
ACERadialFunctions::chebPow(DOUBLE_TYPE lam, DOUBLE_TYPE cut, DOUBLE_TYPE dcut, DOUBLE_TYPE r) {
    DOUBLE_TYPE y, dy, x, dx;
    /* scaled distance x and derivative*/
    y = (1.0 - r / cut);
    dy = pow(y, (lam - 1.0));
    y = dy * y;
    dy = -lam / cut * dy;

    x = 2.0 * (1.0 - y) - 1.0;
    dx = -2.0 * dy;
    calcCheb(nradbase, x);
    for (NS_TYPE n = 1; n <= nradbase; n++) {
        gr(n - 1) = 0.5 - 0.5 * cheb(n);
        dgr(n - 1) = -0.5 * dcheb(n) * dx;
    }
}


void
ACERadialFunctions::chebLinear(DOUBLE_TYPE lam, DOUBLE_TYPE cut, DOUBLE_TYPE dcut, DOUBLE_TYPE r) {
    DOUBLE_TYPE x, dx;
    /* scaled distance x and derivative*/
    x = (1.0 - r / cut);
    dx = -1 / cut;
    calcCheb(nradbase, x);
    for (NS_TYPE n = 1; n <= nradbase; n++) {
        gr(n - 1) = 0.5 - 0.5 * cheb(n);
        dgr(n - 1) = -0.5 * dcheb(n) * dx;
    }
}

/**
Function that computes radial functions.

@param nradbase, nelements, elei, elej

@returns fr, dfr
*/
void ACERadialFunctions::radfunc(SPECIES_TYPE elei, SPECIES_TYPE elej) {
    DOUBLE_TYPE frval, dfrval;
    for (NS_TYPE n = 0; n < nradial; n++) {
        for (LS_TYPE l = 0; l <= lmax; l++) {
            frval = 0.0;
            dfrval = 0.0;
            for (NS_TYPE k = 0; k < nradbase; k++) {
                frval += crad(elei, elej, n, l, k) * gr(k);
                dfrval += crad(elei, elej, n, l, k) * dgr(k);
            }
            fr(n, l) = frval;
            dfr(n, l) = dfrval;
        }
    }
}


void ACERadialFunctions::all_radfunc(SPECIES_TYPE mu_i, SPECIES_TYPE mu_j, DOUBLE_TYPE r) {
    DOUBLE_TYPE lam = lambda(mu_i, mu_j);
    DOUBLE_TYPE r_cut = cut(mu_i, mu_j);
    DOUBLE_TYPE dr_cut = dcut(mu_i, mu_j);
    // set up radial functions
    radbase(lam, r_cut, dr_cut, r); //update gr, dgr
    radfunc(mu_i, mu_j); // update fr(nr, l),  dfr(nr, l)
}


void ACERadialFunctions::setuplookupRadspline() {
    using namespace std::placeholders;
    DOUBLE_TYPE lam, r_cut, dr_cut;
    DOUBLE_TYPE cr_c, dcr_c, pre, lamhc;

    // at r = rcut + eps the function and its derivatives is zero
    for (SPECIES_TYPE elei = 0; elei < nelements; elei++) {
        for (SPECIES_TYPE elej = 0; elej < nelements; elej++) {

            lam = lambda(elei, elej);
            r_cut = cut(elei, elej);
            dr_cut = dcut(elei, elej);

            splines_gk(elei, elej).setupSplines(gr.get_size(),
                                                std::bind(&ACERadialFunctions::radbase, this, lam, r_cut, dr_cut,
                                                          _1),//update gr, dgr
                                                gr.get_data(),
                                                dgr.get_data(), deltaSplineBins, cutoff);

            splines_rnl(elei, elej).setupSplines(fr.get_size(),
                                                 std::bind(&ACERadialFunctions::all_radfunc, this, elei, elej,
                                                           _1), // update fr(nr, l),  dfr(nr, l)
                                                 fr.get_data(),
                                                 dfr.get_data(), deltaSplineBins, cutoff);


            pre = prehc(elei, elej);
            lamhc = lambdahc(elei, elej);
//            radcore(r, pre, lamhc, cutoff, cr_c, dcr_c);
            splines_hc(elei, elej).setupSplines(1,
                                                std::bind(&ACERadialFunctions::radcore, _1, pre, lamhc, cutoff,
                                                          std::ref(cr_c), std::ref(dcr_c)),
                                                &cr_c,
                                                &dcr_c, deltaSplineBins, cutoff);
        }
    }

}

/**
Function that gets radial function from look-up table using splines.

@param r, nradbase_c, nradial_c, lmax, mu_i, mu_j

@returns fr, dfr, gr, dgr, cr, dcr
*/
void
ACERadialFunctions::evaluate(DOUBLE_TYPE r, NS_TYPE nradbase_c, NS_TYPE nradial_c, SPECIES_TYPE mu_i,
                             SPECIES_TYPE mu_j, bool calc_second_derivatives) {
    auto &spline_gk = splines_gk(mu_i, mu_j);
    auto &spline_rnl = splines_rnl(mu_i, mu_j);
    auto &spline_hc = splines_hc(mu_i, mu_j);

    spline_gk.calcSplines(r, calc_second_derivatives); // populate  splines_gk.values, splines_gk.derivatives;
    for (NS_TYPE nr = 0; nr < nradbase_c; nr++) {
        gr(nr) = spline_gk.values(nr);
        dgr(nr) = spline_gk.derivatives(nr);
        if (calc_second_derivatives)
            d2gr(nr) = spline_gk.second_derivatives(nr);
    }

    spline_rnl.calcSplines(r, calc_second_derivatives);
    for (size_t ind = 0; ind < fr.get_size(); ind++) {
        fr.get_data(ind) = spline_rnl.values.get_data(ind);
        dfr.get_data(ind) = spline_rnl.derivatives.get_data(ind);
        if (calc_second_derivatives)
            d2fr.get_data(ind) = spline_rnl.second_derivatives.get_data(ind);
    }

    spline_hc.calcSplines(r, calc_second_derivatives);
    cr = spline_hc.values(0);
    dcr = spline_hc.derivatives(0);
    if (calc_second_derivatives)
        d2cr = spline_hc.second_derivatives(0);
}


void
ACERadialFunctions::radcore(DOUBLE_TYPE r, DOUBLE_TYPE pre, DOUBLE_TYPE lambda, DOUBLE_TYPE cutoff, DOUBLE_TYPE &cr,
                            DOUBLE_TYPE &dcr) {
/* pseudocode for hard core repulsion
in:
 r: distance
 pre: prefactor: read from input, depends on pair of atoms mu_i mu_j
 lambda: exponent: read from input, depends on pair of atoms mu_i mu_j
 cutoff: cutoff distance: read from input, depends on pair of atoms mu_i mu_j
out:
cr: hard core repulsion
dcr: derivative of hard core repulsion

 function
 \$f f_{core} = pre \exp( - \lambda r^2 ) / r   \$f

*/

    DOUBLE_TYPE r2, lr2, y, x0, env, denv;

//   repulsion strictly positive and decaying
    pre = abs(pre);
    lambda = abs(lambda);

    r2 = r * r;
    lr2 = lambda * r2;
    if (lr2 < 50.0) {
        y = exp(-lr2);
        cr = pre * y / r;
        dcr = -pre * y * (2.0 * lr2 + 1.0) / r2;

        x0 = r / cutoff;
        env = 0.5 * (1.0 + cos(pi * x0));
        denv = -0.5 * sin(pi * x0) * pi / cutoff;
        dcr = cr * denv + dcr * env;
        cr = cr * env;
    } else {
        cr = 0.0;
        dcr = 0.0;
    }

}

void
ACERadialFunctions::evaluate_range(vector<DOUBLE_TYPE> r_vec, NS_TYPE nradbase_c, NS_TYPE nradial_c, SPECIES_TYPE mu_i,
                                   SPECIES_TYPE mu_j) {
    if (nradbase_c > nradbase)
        throw invalid_argument("nradbase_c couldn't be larger than nradbase");
    if (nradial_c > nradial)
        throw invalid_argument("nradial_c couldn't be larger than nradial");
    if (mu_i > nelements)
        throw invalid_argument("mu_i couldn't be larger than nelements");
    if (mu_j > nelements)
        throw invalid_argument("mu_j couldn't be larger than nelements");

    gr_vec.resize(r_vec.size(), nradbase_c);
    dgr_vec.resize(r_vec.size(), nradbase_c);
    d2gr_vec.resize(r_vec.size(), nradbase_c);

    fr_vec.resize(r_vec.size(), fr.get_dim(0), fr.get_dim(1));
    dfr_vec.resize(r_vec.size(), fr.get_dim(0), fr.get_dim(1));
    d2fr_vec.resize(r_vec.size(), fr.get_dim(0), fr.get_dim(1));

    for (size_t i = 0; i < r_vec.size(); i++) {
        DOUBLE_TYPE r = r_vec[i];
        this->evaluate(r, nradbase_c, nradial_c, mu_i, mu_j, true);
        for (NS_TYPE nr = 0; nr < nradbase_c; nr++) {
            gr_vec(i, nr) = gr(nr);
            dgr_vec(i, nr) = dgr(nr);
            d2gr_vec(i, nr) = d2gr(nr);
        }

        for (NS_TYPE nr = 0; nr < nradial_c; nr++) {
            for (LS_TYPE l = 0; l <= lmax; l++) {
                fr_vec(i, nr, l) = fr(nr, l);
                dfr_vec(i, nr, l) = dfr(nr, l);
                d2fr_vec(i, nr, l) = d2fr(nr, l);

            }
        }
    }

}

void SplineInterpolator::setupSplines(int num_of_functions, RadialFunctions func,
                                      DOUBLE_TYPE *values,
                                      DOUBLE_TYPE *dvalues, DOUBLE_TYPE deltaSplineBins, DOUBLE_TYPE cutoff) {

    this->deltaSplineBins = deltaSplineBins;
    this->cutoff = cutoff;
    this->ntot = static_cast<int>(cutoff / deltaSplineBins);

    DOUBLE_TYPE r, c[4];
    this->num_of_functions = num_of_functions;
    this->values.resize(num_of_functions);
    this->derivatives.resize(num_of_functions);
    this->second_derivatives.resize(num_of_functions);

    Array1D<DOUBLE_TYPE> f1g(num_of_functions);
    Array1D<DOUBLE_TYPE> f1gd1(num_of_functions);
    f1g.fill(0);
    f1gd1.fill(0);

    nlut = ntot;
    DOUBLE_TYPE f0, f1, f0d1, f1d1;
    int idx;

    // cutoff is global cutoff
    rscalelookup = (DOUBLE_TYPE) nlut / cutoff;
    invrscalelookup = 1.0 / rscalelookup;

    lookupTable.init(ntot + 1, num_of_functions, 4);
    if (values == nullptr & num_of_functions > 0)
        throw invalid_argument("SplineInterpolator::setupSplines: values could not be null");
    if (dvalues == nullptr & num_of_functions > 0)
        throw invalid_argument("SplineInterpolator::setupSplines: dvalues could not be null");

    for (int n = nlut; n >= 1; n--) {
        r = invrscalelookup * DOUBLE_TYPE(n);
        func(r); //populate values and dvalues arrays
        for (int func_id = 0; func_id < num_of_functions; func_id++) {
            f0 = values[func_id];
            f1 = f1g(func_id);
            f0d1 = dvalues[func_id] * invrscalelookup;
            f1d1 = f1gd1(func_id);
            // evaluate coefficients
            c[0] = f0;
            c[1] = f0d1;
            c[2] = 3.0 * (f1 - f0) - f1d1 - 2.0 * f0d1;
            c[3] = -2.0 * (f1 - f0) + f1d1 + f0d1;
            // store coefficients
            for (idx = 0; idx <= 3; idx++)
                lookupTable(n, func_id, idx) = c[idx];

            // evaluate function values and derivatives at current position
            f1g(func_id) = c[0];
            f1gd1(func_id) = c[1];
        }
    }


}


void SplineInterpolator::calcSplines(DOUBLE_TYPE r, bool calc_second_derivatives) {
    DOUBLE_TYPE wl, wl2, wl3, w2l1, w3l2, w4l2;
    DOUBLE_TYPE c[4];
    int func_id, idx;
    DOUBLE_TYPE x = r * rscalelookup;
    int nl = static_cast<int>(floor(x));

    if (nl <= 0)
        throw std::invalid_argument("Encountered very small distance. Stopping.");

    if (nl < nlut) {
        wl = x - DOUBLE_TYPE(nl);
        wl2 = wl * wl;
        wl3 = wl2 * wl;
        w2l1 = 2.0 * wl;
        w3l2 = 3.0 * wl2;
        w4l2 = 6.0 * wl;
        for (func_id = 0; func_id < num_of_functions; func_id++) {
            for (idx = 0; idx <= 3; idx++) {
                c[idx] = lookupTable(nl, func_id, idx);
            }
            values(func_id) = c[0] + c[1] * wl + c[2] * wl2 + c[3] * wl3;
            derivatives(func_id) = (c[1] + c[2] * w2l1 + c[3] * w3l2) * rscalelookup;
            if (calc_second_derivatives)
                second_derivatives(func_id) = (c[2] + c[3] * w4l2) * rscalelookup * rscalelookup * 2;
        }
    } else { // fill with zeroes
        values.fill(0);
        derivatives.fill(0);
        if (calc_second_derivatives)
            second_derivatives.fill(0);
    }
}

void SplineInterpolator::calcSplines(DOUBLE_TYPE r, SHORT_INT_TYPE func_ind) {
    DOUBLE_TYPE wl, wl2, wl3, w2l1, w3l2;
    DOUBLE_TYPE c[4];
    int idx;
    DOUBLE_TYPE x = r * rscalelookup;
    int nl = static_cast<int>(floor(x));

    if (nl <= 0)
        throw std::invalid_argument("Encountered very small distance. Stopping.");

    if (nl < nlut) {
        wl = x - DOUBLE_TYPE(nl);
        wl2 = wl * wl;
        wl3 = wl2 * wl;
        w2l1 = 2.0 * wl;
        w3l2 = 3.0 * wl2;

        for (idx = 0; idx <= 3; idx++) {
            c[idx] = lookupTable(nl, func_ind, idx);
        }
        values(func_ind) = c[0] + c[1] * wl + c[2] * wl2 + c[3] * wl3;
        derivatives(func_ind) = (c[1] + c[2] * w2l1 + c[3] * w3l2) * rscalelookup;

    } else { // fill with zeroes
        values(func_ind) = 0;
        derivatives(func_ind) = 0;
    }
}
