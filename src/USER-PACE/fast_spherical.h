#ifndef ACE_FAST_SPHERICAL_H
#define ACE_FAST_SPHERICAL_H

#include <cmath>

#include "ace_types.h"
#include "complex.h"
#include "ace_utils.h"

#include "USER-PACE/include/multiarray/multiarray2dlm.h"
#include "USER-PACE/include/multiarray/multiarray_auto.h"

using namespace std;


const DOUBLE_TYPE sq3o4pi = 0.48860251190291992158; // sqrt(3/(4*pi))
const DOUBLE_TYPE sq1o4pi = 0.28209479177387814347; // sqrt(1/(4*pi))
const DOUBLE_TYPE sq3o8pi = 0.34549414947133547927; // sqrt(3/(8*pi))

//definition of common factor for spherical harmonics = Y00
const DOUBLE_TYPE Y00 = sq1o4pi;

/**
Class to store spherical harmonics and their associated functions. \n
All the associated members such as \f$ P_{lm}, Y_{lm}\f$ etc are one dimensional arrays of length (L+1)*(L+2)/2. \n
The value that corresponds to a particular l, m configuration can be accessed through a preprocessor directive as \n
\code ylm[at(l,m)] \endcode \n
which can access the (m+(l*(l+1))/2) value from the one dimensional array.
*/
class FastSHarmonics {
public:

    /**
    int, the number of spherical harmonics to be found
    */
    LS_TYPE lmax;

    FastSHarmonics() = default;

    explicit FastSHarmonics(LS_TYPE lmax);

    void init(LS_TYPE lm);

    ~FastSHarmonics();

    void pre_compute();

    void compute_barplm(DOUBLE_TYPE rz, LS_TYPE lmaxi);

    void compute_ylm(DOUBLE_TYPE rx, DOUBLE_TYPE ry, DOUBLE_TYPE rz, LS_TYPE lmaxi);

    Array2DLM <DOUBLE_TYPE> alm;
    Array2DLM <DOUBLE_TYPE> blm;
    Array1D <DOUBLE_TYPE> cl;
    Array1D <DOUBLE_TYPE> dl;

    Array2DLM <DOUBLE_TYPE> plm;
    Array2DLM <DOUBLE_TYPE> dplm;

    Array2DLM <Complex> ylm;
    Array2DLM <Dycomponent> dylm;

};


#endif
