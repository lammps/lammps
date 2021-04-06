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

#include <cmath>

#include "ace_spherical_cart.h"

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
Destructor for ACECartesianSphericalHarmonics.

@param None

@returns None
*/
ACECartesianSphericalHarmonics::~ACECartesianSphericalHarmonics() {}


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


void ACECartesianSphericalHarmonics::compute_barplm(DOUBLE_TYPE rz, LS_TYPE lmaxi) {

    // requires -1 <= rz <= 1 , NO CHECKING IS PERFORMED !!!!!!!!!
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

}

