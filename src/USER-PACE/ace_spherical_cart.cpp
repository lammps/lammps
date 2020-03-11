#include <cmath>

#include "ace_spherical_cart.h"


/**
Constructor for SHarmonics. Dynamically initialises all the arrays.

@param lmax, int

The value of lmax

@returns None
*/
ACECartesianSphericalHarmonics::ACECartesianSphericalHarmonics(LS_TYPE lm) {
    init(lm);
}

void ACECartesianSphericalHarmonics::init(LS_TYPE lm) {
    lmax = lm;

    alm.init(lmax, "alm");
    blm.init(lmax, "blm");
    cl.init(lmax + 1);
    dl.init(lmax + 1);

    plm.init(lmax, "plm");
    dplm.init(lmax, "dplm");

    ylm.init(lmax, "ylm");
    dylm.init(lmax, "dylm");

    pre_compute();
}

/**
Destructor for SHarmonics. Frees the memory of all the arrays.

@param None

@returns None
*/
ACECartesianSphericalHarmonics::~ACECartesianSphericalHarmonics() {}


/**
Precomputes the value of \f$ a_{lm}, b_{lm} \f$ values.

@param None

@returns None
*/
void ACECartesianSphericalHarmonics::pre_compute() {

    DOUBLE_TYPE a, b;
    DOUBLE_TYPE lsq, ld, l1, l2;
    DOUBLE_TYPE msq;

    for (LS_TYPE l = 1; l <= lmax; l++) {
        lsq = l * l;
        ld = 2 * l;
        l1 = (4 * lsq - 1);
        l2 = lsq - ld + 1;
        for (MS_TYPE m = 0; m < l - 1; m++) {
            msq = m * m;
            a = sqrt((DOUBLE_TYPE(l1)) / (DOUBLE_TYPE(lsq - msq)));
            b = -sqrt((DOUBLE_TYPE(l2 - msq)) / (DOUBLE_TYPE(4 * l2 - 1)));
            alm(l, m) = a;
            blm(l, m) = b;
        }
    }

    for (LS_TYPE l = 1; l <= lmax; l++) {
        cl(l) = -sqrt(1.0 + 0.5 / (DOUBLE_TYPE(l)));
        dl(l) = sqrt(DOUBLE_TYPE(2 * (l - 1) + 3));
    }
}

/**
Function that computes \f$ \bar{P}_{lm} \f$ for the corresponding lmax value
Input is \f$ \hat{r}_z \f$ which is the $z$-component of the bond direction.

For each \f$ \hat{r}_z \f$, this computes the whole range of \f$ \bar{P}_{lm} \f$ values
and its derivatives upto the lmax specified, which is a member of the class.

@param rz, DOUBLE_TYPE

@returns None
*/
void ACECartesianSphericalHarmonics::compute_barplm(DOUBLE_TYPE rz, LS_TYPE lmaxi) {

    // requires -1 <= rz <= 1 , NO CHECKING IS PERFORMED !!!!!!!!!

//    DOUBLE_TYPE sq1o2pi  = 0.39894228040143267794; // sqrt(1/(2*pi))
//    DOUBLE_TYPE sq3o2pi  = 0.69098829894267095853; // sqrt(3/(2*pi))
// prefactors include 1/sqrt(2) factor compared to reference


    DOUBLE_TYPE t;

    // l=0, m=0
    //plm(0, 0) = Y00/sq1o4pi; //= sq1o4pi;
    plm(0, 0) = Y00; //= 1;
    dplm(0, 0) = 0.0;

    if (lmaxi > 0) {

        // l=1, m=0
        plm(1, 0) = Y00 * sq3 * rz;
        dplm(1, 0) = Y00 * sq3;

        // l=1, m=1
        plm(1, 1) = -sq3o2 * Y00;
        dplm(1, 1) = 0.0;

        // loop l = 2, lmax
        for (LS_TYPE l = 2; l <= lmaxi; l++) {
            for (MS_TYPE m = 0; m < l - 1; m++) {
                plm(l, m) = alm(l, m) * (rz * plm(l - 1, m) + blm(l, m) * plm(l - 2, m));
                dplm(l, m) = alm(l, m) * (plm(l - 1, m) + rz * dplm(l - 1, m) + blm(l, m) * dplm(l - 2, m));
            }
            t = dl(l) * plm(l - 1, l - 1);
            plm(l, l - 1) = t * rz;
            dplm(l, l - 1) = t;
            plm(l, l) = cl(l) * plm(l - 1, l - 1);
            dplm(l, l) = 0.0;
        }
    }
}  //end compute_barplm


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

@returns None
*/
void ACECartesianSphericalHarmonics::compute_ylm(DOUBLE_TYPE rx, DOUBLE_TYPE ry, DOUBLE_TYPE rz, LS_TYPE lmaxi) {

    // requires rx^2 + ry^2 + rz^2 = 1 , NO CHECKING IS PERFORMED !!!!!!!!!

    DOUBLE_TYPE real;
    DOUBLE_TYPE img;
    MS_TYPE m;
    ACEComplex phase;
    ACEComplex phasem, mphasem1;
    ACEComplex dyx, dyy, dyz;
    ACEComplex rdy;

    phase.real = rx;
    phase.img = ry;
    //compute barplm
    compute_barplm(rz, lmaxi);

    //m = 0
    m = 0;
    for (LS_TYPE l = 0; l <= lmaxi; l++) {

        ylm(l, m).real = plm(l, m);
        ylm(l, m).img = 0.0;

        dyz.real = dplm(l, m);
        rdy.real = dyz.real * rz;

        dylm(l, m).a[0].real = -rdy.real * rx;
        dylm(l, m).a[0].img = 0.0;
        dylm(l, m).a[1].real = -rdy.real * ry;
        dylm(l, m).a[1].img = 0.0;
        dylm(l, m).a[2].real = dyz.real - rdy.real * rz;
        dylm(l, m).a[2].img = 0;
    }
    //m = 0
    m = 1;
    for (LS_TYPE l = 1; l <= lmaxi; l++) {

        ylm(l, m) = phase * plm(l, m);

//        std::cout << "Re ylm(" << l << "," << m <<")= " << ylm(l, m).real << std::endl;
//        std::cout << "Im ylm(" << l << "," << m <<")= " << ylm(l, m).img << std::endl;

        dyx.real = plm(l, m);
        dyx.img = 0.0;
        dyy.real = 0.0;
        dyy.img = plm(l, m);
        dyz.real = phase.real * dplm(l, m);
        dyz.img = phase.img * dplm(l, m);

        rdy.real = rx * dyx.real + +rz * dyz.real;
        rdy.img = ry * dyy.img + rz * dyz.img;

        dylm(l, m).a[0].real = dyx.real - rdy.real * rx;
        dylm(l, m).a[0].img = -rdy.img * rx;
        dylm(l, m).a[1].real = -rdy.real * ry;
        dylm(l, m).a[1].img = dyy.img - rdy.img * ry;
        dylm(l, m).a[2].real = dyz.real - rdy.real * rz;
        dylm(l, m).a[2].img = dyz.img - rdy.img * rz;
    }

    // m > 1
    phasem = phase;
    for (MS_TYPE m = 2; m <= lmaxi; m++) {

        mphasem1.real = phasem.real * DOUBLE_TYPE(m);
        mphasem1.img = phasem.img * DOUBLE_TYPE(m);
        phasem = phasem * phase;

        for (LS_TYPE l = m; l <= lmaxi; l++) {

            ylm(l, m).real = phasem.real * plm(l, m);
            ylm(l, m).img = phasem.img * plm(l, m);

            dyx = mphasem1 * plm(l, m);
            dyy.real = -dyx.img;
            dyy.img = dyx.real;
            dyz = phasem * dplm(l, m);

            rdy.real = rx * dyx.real + ry * dyy.real + rz * dyz.real;
            rdy.img = rx * dyx.img + ry * dyy.img + rz * dyz.img;

            dylm(l, m).a[0].real = dyx.real - rdy.real * rx;
            dylm(l, m).a[0].img = dyx.img - rdy.img * rx;
            dylm(l, m).a[1].real = dyy.real - rdy.real * ry;
            dylm(l, m).a[1].img = dyy.img - rdy.img * ry;
            dylm(l, m).a[2].real = dyz.real - rdy.real * rz;
            dylm(l, m).a[2].img = dyz.img - rdy.img * rz;
        }
    }

    //fill-in m<0
//    for (LS_TYPE l = 1; l <= lmaxi; l++) {
//        for (MS_TYPE m = 1; m <= l; m++) {
//            auto phase = DOUBLE_TYPE(abs(m) % 2 == 0 ? 1 : -1);
//
//            //Y_l,{-m} = (-1)^m (Y_l,{m})c.c.
//            ylm(l, -m) = ylm(l, m) * phase;
//            ylm(l, -m).img *= -1;
//
//            //DY_l,{-m} = (-1)^m (DY_l,{m})c.c.
//            dylm(l, -m) = dylm(l, m);
//            dylm(l, -m) *= phase;
//            //complex conjugate
//            dylm(l, -m).a[0].img *= -1;
//            dylm(l, -m).a[1].img *= -1;
//            dylm(l, -m).a[2].img *= -1;
//
//
//        }
//    }
}

