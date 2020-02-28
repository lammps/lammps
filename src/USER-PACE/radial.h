#ifndef ACE_RADIAL_FUNCTIONS_H
#define ACE_RADIAL_FUNCTIONS_H

#include "ace_types.h"
#include "multiarray_auto.h"

using namespace std;

/**
Class to store radial functions and their associated functions. \n
*/
class RFunctions {
public:
    RFunctions() = default;

    RFunctions(NS_TYPE nradb, LS_TYPE lmax, NS_TYPE nradial, int ntot, SPECIES_TYPE nelements, DOUBLE_TYPE cutoff);

    void init(NS_TYPE nradb, LS_TYPE lmax, NS_TYPE nradial, int ntot, SPECIES_TYPE nelements, DOUBLE_TYPE cutoff);

    ~RFunctions();

    void calcCheb(NS_TYPE n, DOUBLE_TYPE x);

    void radbase(DOUBLE_TYPE lam, DOUBLE_TYPE r_c, DOUBLE_TYPE d_r_c, DOUBLE_TYPE r);

    void setuplookupRadspline();

    void radfunc(SPECIES_TYPE elei, SPECIES_TYPE elej);

    void lookupRadspline(DOUBLE_TYPE r, NS_TYPE nradbase, NS_TYPE nradial, SPECIES_TYPE elei, SPECIES_TYPE elej);

    /**
    Arrays to store radial functions.
    */
    Array1D<DOUBLE_TYPE> gr;
    Array1D<DOUBLE_TYPE> dgr;
    Array2D<DOUBLE_TYPE> fr;//done
    Array2D<DOUBLE_TYPE> dfr;//done

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
    Array5D<DOUBLE_TYPE> crad;
    Array3D<DOUBLE_TYPE> cradl0; //crad for l = 0; then nradial2 = nradbase

    Array2D<DOUBLE_TYPE> lambda;
    Array2D<DOUBLE_TYPE> cut;
    Array2D<DOUBLE_TYPE> dcut;
    Array2D<DOUBLE_TYPE> f1f; //done
    Array2D<DOUBLE_TYPE> f1fd1;

    Array1D<DOUBLE_TYPE> f1g;
    Array1D<DOUBLE_TYPE> f1gd1;


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