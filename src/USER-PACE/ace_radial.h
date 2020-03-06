#ifndef ACE_RADIAL_FUNCTIONS_H
#define ACE_RADIAL_FUNCTIONS_H

#include "ace_types.h"
#include "multiarray/ace_arraynd.h"

using namespace std;

/**
Class to store radial functions and their associated functions. \n
*/
class ACERadialFunctions {
public:
    ACERadialFunctions() = default;

    ACERadialFunctions(NS_TYPE nradb, LS_TYPE lmax, NS_TYPE nradial, int ntot, SPECIES_TYPE nelements,
                       DOUBLE_TYPE cutoff);

    void init(NS_TYPE nradb, LS_TYPE lmax, NS_TYPE nradial, int ntot, SPECIES_TYPE nelements, DOUBLE_TYPE cutoff);

    ~ACERadialFunctions();

    void calcCheb(NS_TYPE n, DOUBLE_TYPE x);

    void radbase(DOUBLE_TYPE lam, DOUBLE_TYPE r_c, DOUBLE_TYPE d_r_c, DOUBLE_TYPE r);

    static void
    radcore(DOUBLE_TYPE r, DOUBLE_TYPE pre, DOUBLE_TYPE lambda, DOUBLE_TYPE cutoff, DOUBLE_TYPE &cr, DOUBLE_TYPE &dcr);

    void setuplookupRadspline();

    void radfunc(SPECIES_TYPE elei, SPECIES_TYPE elej);

    void lookupRadspline(DOUBLE_TYPE r, NS_TYPE nradbase, NS_TYPE nradial, SPECIES_TYPE elei, SPECIES_TYPE elej);

    /**
    Arrays to store radial functions.
    */
    Array1D<DOUBLE_TYPE> gr; // g_k(r) functions
    Array1D<DOUBLE_TYPE> dgr; //derivatives of g_k(r) functions
    Array2D<DOUBLE_TYPE> fr;  //R_nl(r) functions
    Array2D<DOUBLE_TYPE> dfr; //derivatives of R_nl(r) functions

    //hard-core repulsion
    DOUBLE_TYPE cr, dcr;

    /**
    Arrays to store Chebyshev polynomials.
    */
    Array1D<DOUBLE_TYPE> cheb; // Chebyshev polynomials of the first kind
    Array1D<DOUBLE_TYPE> dcheb; // derivatives Chebyshev polynomials of the first kind
    Array1D<DOUBLE_TYPE> cheb2; // Chebyshev polynomials of the second kind
//TODO make look-up tables an independent class
    /**
    Variables for look-up tables.
    */
    DOUBLE_TYPE rscalelookup;
    DOUBLE_TYPE invrscalelookup;
    int nlut;
    /**
    Arrays for look-up tables.
    */
    Array6D<DOUBLE_TYPE> lutfrs;
    Array5D<DOUBLE_TYPE> lutgrs;
    Array4D<DOUBLE_TYPE> luthcs; // hard-core repulsion



    Array5D<DOUBLE_TYPE> crad;
    Array2D<DOUBLE_TYPE> lambda;
    Array2D<DOUBLE_TYPE> cut;
    Array2D<DOUBLE_TYPE> dcut;
    Array2D<DOUBLE_TYPE> f1f;
    Array2D<DOUBLE_TYPE> f1fd1;

    Array1D<DOUBLE_TYPE> f1g;
    Array1D<DOUBLE_TYPE> f1gd1;

    //hard-core repulsion coefficients
    Array2D<DOUBLE_TYPE> prehc;
    Array2D<DOUBLE_TYPE> lambdahc;


//--------------------------------------------------------------------------

    SPECIES_TYPE nelements;
    LS_TYPE lmax;
    NS_TYPE nradial, nradbase;
    DOUBLE_TYPE cutoff;


    /**
    Number of bins for look-up tables.
    */
    int ntot = 10000;
//--------------------------------------------------------------------------
};

#endif