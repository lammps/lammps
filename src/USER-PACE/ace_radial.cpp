#include <cmath>

#include "ace_radial.h"
#include "ace_utils.h"

/**
Constructor for ACERadialFunctions.

@param None

@returns None
*/
ACERadialFunctions::ACERadialFunctions(NS_TYPE nradb, LS_TYPE lmax, NS_TYPE nradial, int ntot, SPECIES_TYPE nelements,
                                       DOUBLE_TYPE cutoff) {
    init(nradb, lmax, nradial, ntot, nelements, cutoff);
}

void ACERadialFunctions::init(NS_TYPE nradb, LS_TYPE lmax, NS_TYPE nradial, int ntot, SPECIES_TYPE nelements,
                              DOUBLE_TYPE cutoff) {
    int s;
    this->nradbase = nradb;
    this->lmax = lmax;
    this->nradial = nradial;
    this->ntot = ntot;
    this->nelements = nelements;
    this->cutoff = cutoff;

    gr.init(nradbase + 1, "gr");
    dgr.init(nradbase + 1, "dgr");

    f1g.init(nradbase + 1, "f1g");
    f1gd1.init(nradbase + 1, "f1gd1");

    fr.init(nradbase, lmax + 1, "fr");
    dfr.init(nradbase, lmax + 1, "dfr");

    f1f.init(lmax + 1, nradbase, "f1f");
    f1fd1.init(lmax + 1, nradbase, "f1fd1");


    cheb.init(nradbase + 1, "cheb");
    dcheb.init(nradbase + 1, "dcheb");
    cheb2.init(nradbase + 1, "cheb2");


    //lutfrs.init(4, nradmax, 1 + lmax, ntot, nelements, nelements, "lutfrs");
    lutfrs.init(nelements, nelements, ntot, 1 + lmax, nradial, 4, "lutfrs");
    //lutgrs.init(4, nradbase, ntot, nelements, nelements, "lutgrs");
    lutgrs.init(nelements, nelements, ntot, nradbase, 4, "lutgrs");
    lambda.init(nelements, nelements, "lambda");
    lambda.fill(1.);

    cut.init(nelements, nelements, "cut");
    cut.fill(1.);

    dcut.init(nelements, nelements, "dcut");
    dcut.fill(1.);

    crad.init(nelements, nelements, (lmax + 1), nradial, nradbase, "crad");
    crad.fill(0.);

}

/**
Destructor for ACERadialFunctions.

@param None

@returns None
*/
ACERadialFunctions::~ACERadialFunctions() {
}

/**
Function that computes Chebyshev polynomials of first and second kind
 to setup the radial functions and the derivatives

@param n, x

@returns cheb1, dcheb1
*/
void ACERadialFunctions::calcCheb(NS_TYPE n, DOUBLE_TYPE x) {
    if (n < 0) {
        fprintf(stderr, "The order n of the polynomials should be positive %d\n", n);
        exit(EXIT_FAILURE);
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

@param lam, nradbase, r_c, d_r_c, r

@returns gr, dgr
*/
void ACERadialFunctions::radbase(DOUBLE_TYPE lam, DOUBLE_TYPE r_c, DOUBLE_TYPE d_r_c, DOUBLE_TYPE r) {
    /*lam is given by the formula (24), that contains r_c */
    DOUBLE_TYPE x1, x2, x, dx;
    DOUBLE_TYPE env, denv, fcut, dfcut;
    if (r >= r_c) {
        gr.fill(0);
        dgr.fill(0);
    }
    /* scaled distance x and derivative*/
    x1 = exp(-lam);
    x2 = exp(-lam * r / r_c);
    x = 1.0 - 2 * ((x2 - x1) / (1 - x1));
    dx = 2 * (lam / r_c) * (x2 / (1 - x1));
    /* calculation of Chebyshev polynomials from the recursion */
    calcCheb(nradbase - 1, x);
#ifdef DEBUG_RADIAL
    for(int ii = 0; ii<nradbase; ii++)
        printf("cheb(%d) = %f, dcheb(%d) = %f\n",ii, cheb(ii), ii, dcheb(ii));
#endif
    gr(0) = cheb(0);
    dgr(0) = dcheb(0) * dx;
    for (NS_TYPE n = 2; n <= nradbase; n++) {
        gr(n - 1) = 0.5 - 0.5 * cheb(n - 1);
        dgr(n - 1) = -0.5 * dcheb(n - 1) * dx;
#ifdef DEBUG_RADIAL
        printf("1: n %d  gr[n] %f  dgr[n] %f\n", n-1,gr(n-1), dgr(n-1));
#endif
    }
    env = 0.5 * (1.0 + cos(pi * r / r_c));
    denv = -0.5 * sin(pi * r / r_c) * pi / r_c;
    for (NS_TYPE n = 0; n < nradbase; n++) {
        dgr(n) = gr(n) * denv + dgr(n) * env;
        gr(n) = gr(n) * env;
#ifdef DEBUG_RADIAL
        printf("2: n %d  gr[n] %f  dgr[n] %f\n", n,gr(n), dgr(n));
#endif
    }
    // for radtype = 3 a smooth cutoff is already included in the basis function
    dx = r_c - d_r_c;
    if (r > dx) {
        fcut = 0.5 * (1.0 + cos(pi * (r - dx) / d_r_c));
        dfcut = -0.5 * sin(pi * (r - dx) / d_r_c) * pi / d_r_c;
        for (NS_TYPE n = 0; n < nradbase; n++) {
            dgr(n) = gr(n) * dfcut + dgr(n) * fcut;
            gr(n) = gr(n) * fcut;
#ifdef DEBUG_RADIAL
            printf("3: n %d  gr[n] %f  dgr[n] %f\n", n,gr(n), dgr(n));
#endif
        }
    }
}

/**
Function that computes radial functions.

@param nradbase, nelements, elei, elej

@returns fr, dfr
*/
void ACERadialFunctions::radfunc(SPECIES_TYPE elei, SPECIES_TYPE elej) {
    DOUBLE_TYPE frval, dfrval;
    for (NS_TYPE nr = 0; nr < nradial; nr++) {
        for (LS_TYPE l = 0; l <= lmax; l++) {
            frval = 0.0;
            dfrval = 0.0;
            for (NS_TYPE idx = 0; idx < nradbase; idx++) {
                frval += crad(elei, elej, l, nr, idx) * gr(idx);
                dfrval += crad(elei, elej, l, nr, idx) * dgr(idx);
            }
            fr(nr, l) = frval;
            dfr(nr, l) = dfrval;
        }
    }
}

/**
Function that sets up the look-up tables for spline-representation of radial functions.

@param ntot

@returns None
*/
void ACERadialFunctions::setuplookupRadspline() {
    SPECIES_TYPE elei, elej;
    int s, n, l, idx;
    NS_TYPE nr;
    DOUBLE_TYPE x, lam, cu, dc, f0, f1, f1d1, f0d1;
    DOUBLE_TYPE c[4];
    nlut = ntot;
    // cutoff is global cutoff
    rscalelookup = (DOUBLE_TYPE) nlut / cutoff;
#ifdef DEBUG_RADIAL
    printf("rscalelookup=%f\n", rscalelookup);
#endif
    invrscalelookup = 1.0 / rscalelookup;
    lutfrs.fill(0.0);
    lutgrs.fill(0.0);
    // at r = rcut + eps the function and its derivatives is zero
    s = nradbase;
    f1g.fill(0);
    f1gd1.fill(0);
    for (elei = 0; elei < nelements; elei++) {
        for (elej = 0; elej < nelements; elej++) {
            for (n = nlut - 1; n >= 0; n--) {
                x = invrscalelookup * DOUBLE_TYPE(n);
                lam = lambda(elei, elej);
                cu = cut(elei, elej);
                dc = dcut(elei, elej);
                // set up radial functions
                radbase(lam, cu, dc, x);
                radfunc(elei, elej);
                for (nr = 0; nr < nradbase; nr++) {
                    f0 = gr(nr);
                    f1 = f1g(nr);
                    f0d1 = dgr(nr) * invrscalelookup;
                    f1d1 = f1gd1(nr);
                    // evaluate coefficients
                    c[0] = f0;
                    c[1] = f0d1;
                    c[2] = 3.0 * (f1 - f0) - f1d1 - 2.0 * f0d1;
                    c[3] = -2.0 * (f1 - f0) + f1d1 + f0d1;
                    // store coefficients
                    for (idx = 0; idx <= 3; idx++) {
                        lutgrs(elei, elej, n, nr, idx) = c[idx];
                    }
                    // evalute function values and derivatives at current position
                    f1g(nr) = c[0];
                    f1gd1(nr) = c[1];
                }
                for (nr = 0; nr < nradial; nr++) {
                    for (l = 0; l <= lmax; l++) {
                        f0 = fr(nr, l);
                        f1 = f1f(l, nr);
                        f0d1 = dfr(nr, l) * invrscalelookup;
                        f1d1 = f1fd1(l, nr);
                        // evaluate coefficients
                        c[0] = f0;
                        c[1] = f0d1;
                        c[2] = 3.0 * (f1 - f0) - f1d1 - 2.0 * f0d1;
                        c[3] = -2.0 * (f1 - f0) + f1d1 + f0d1;
                        // store coefficients
                        for (idx = 0; idx <= 3; idx++) {
                            lutfrs(elei, elej, n, l, nr, idx) = c[idx];
                        }
                        // evalute and store function values and derivatives at current position
                        f1f(l, nr) = c[0];
                        f1fd1(l, nr) = c[1];
                    }
                }
            }
        }
    }
}

/**
Function that gets radial function from look-up table using splines.

@param r, nradbase, nradial, lmax, elei, elej

@returns fr, dfr, gr, dgr
*/
void
ACERadialFunctions::lookupRadspline(DOUBLE_TYPE r, NS_TYPE nradbase, NS_TYPE nradial, SPECIES_TYPE elei,
                                    SPECIES_TYPE elej) {
//  when to declare variables in header vs on-the-fly?
    DOUBLE_TYPE x;
    int s, nr, nl, l, idx;
    DOUBLE_TYPE wl, wl2, wl3, w2l1, w3l2;
    DOUBLE_TYPE c[4];
    x = r * rscalelookup;
    nl = static_cast<int>(floor(x));
    if (nl <= 0) {
        fprintf(stderr, "Encountered very small distance.\n Stopping.");
        exit(EXIT_FAILURE);
    }
    if (nl < nlut) {
        wl = x - DOUBLE_TYPE(nl);
        wl2 = wl * wl;
        wl3 = wl2 * wl;
        w2l1 = 2.0 * wl;
        w3l2 = 3.0 * wl2;
        for (nr = 0; nr < nradbase; nr++) {
            for (idx = 0; idx <= 3; idx++) {
                c[idx] = lutgrs(elei, elej, nl, nr, idx);
            }
            gr(nr) = c[0] + c[1] * wl + c[2] * wl2 + c[3] * wl3;
            dgr(nr) = (c[1] + c[2] * w2l1 + c[3] * w3l2) * rscalelookup;
        }
        for (nr = 0; nr < nradial; nr++) {
            for (l = 0; l <= lmax; l++) {
                for (idx = 0; idx <= 3; idx++) {
                    c[idx] = lutfrs(elei, elej, nl, l, nr, idx);
                }
                fr(nr, l) = c[0] + c[1] * wl + c[2] * wl2 + c[3] * wl3;
                dfr(nr, l) = (c[1] + c[2] * w2l1 + c[3] * w3l2) * rscalelookup;
            }
        }
    } else {
        gr.fill(0);
        dgr.fill(0);
        fr.fill(0);
        dfr.fill(0);
    }
}