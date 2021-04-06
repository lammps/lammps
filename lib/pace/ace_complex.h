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

#ifndef ACE_COMPLEX_H
#define ACE_COMPLEX_H


/**
A custom data structure for complex numbers and overloaded operations with them.
*/
struct ACEComplex {
public:
    /**
    Double, real part of the complex number
    */
    DOUBLE_TYPE real;
    /**
    Double, imaginary part of the complex number
    */
    DOUBLE_TYPE img;

    ACEComplex &operator=(const ACEComplex &rhs) = default;

    ACEComplex &operator=(const DOUBLE_TYPE &rhs) {
        this->real = rhs;
        this->img = 0.;
        return *this;
    }

    /**
    Overloading of arithmetical operator += ACEComplex
    */
    ACEComplex &operator+=(const ACEComplex &rhs) {
        this->real += rhs.real;
        this->img += rhs.img;
        return *this; // return the result by reference
    }

    /**
    Overloading of arithmetical operator += DOUBLE_TYPE
    */
    ACEComplex &operator+=(const DOUBLE_TYPE &rhs) {
        this->real += rhs;
        return *this; // return the result by reference
    }

    /**
    Overloading of arithmetical operator *= DOUBLE_TYPE
    */
    ACEComplex &operator*=(const DOUBLE_TYPE &rhs) {
        this->real *= rhs;
        this->img *= rhs;
        return *this; // return the result by reference
    }

    /**
    Overloading of arithmetical operator *= ACEComplex
    */
    ACEComplex &operator*=(const ACEComplex &rhs) {
        DOUBLE_TYPE old_real = this->real;
        this->real = old_real * rhs.real - this->img * rhs.img;
        this->img = old_real * rhs.img + this->img * rhs.real;
        return *this; // return the result by reference
    }

    /**
    Overloading of arithmetical operator *  ACEComplex
    */
    ACEComplex operator*(const ACEComplex &cm2) const {
        ACEComplex res{real * cm2.real - img * cm2.img,
                       real * cm2.img + img * cm2.real};
        return res;
    }

    /*
     * Return complex conjugated copy of itself
     */
    ACEComplex conjugated() const {
        ACEComplex res{real, -img};
        return res;
    }

    /*
     * Complex conjugate itself inplace
     */
    void conjugate() {
        img = -img;
    }

    /*
     *  Multiplication by ACEComplex and return real-part only
     */
    DOUBLE_TYPE real_part_product(const ACEComplex &cm2) const {
        return real * cm2.real - img * cm2.img;
    }

    /*
     *  Multiplication by DOUBLE_TYPE and return real-part only
     */
    DOUBLE_TYPE real_part_product(const DOUBLE_TYPE &cm2) const {
        return real * cm2;
    }

    /*
     *  Overloading of arithmetical operator * by  DOUBLE_TYPE
     */
    ACEComplex operator*(const DOUBLE_TYPE &cm2) const {
        ACEComplex res{real * cm2,
                       img * cm2};
        return res;
    }

    /*
     *  Overloading of arithmetical operator +  ACEComplex
     */
    ACEComplex operator+(const ACEComplex &cm2) const {
        ACEComplex res{real + cm2.real,
                       img + cm2.img};
        return res;
    }

    /*
     *  Overloading of arithmetical operator + with DOUBLE_TYPE
     */
    ACEComplex operator+(const DOUBLE_TYPE &cm2) const {
        ACEComplex res{real + cm2, img};
        return res;
    }

    /*
     *  Overloading of arithmetical operator == ACEComplex
     */
    bool operator==(const ACEComplex &c2) const {
        return (real == c2.real) && (img == c2.img);
    }

    /*
     *  Overloading of arithmetical operator == DOUBLE_TYPE
     */
    bool operator==(const DOUBLE_TYPE &d2) const {
        return (real == d2) && (img == 0.);
    }

    /*
     *  Overloading of arithmetical operator != ACEComplex
     */
    bool operator!=(const ACEComplex &c2) const {
        return (real != c2.real) || (img != c2.img);
    }

    /*
     *  Overloading of arithmetical operator != DOUBLE_TYPE
     */
    bool operator!=(const DOUBLE_TYPE &d2) const {
        return (real != d2) || (img != 0.);
    }

};

/*
 * double * complex for commutativity with complex * double
 */
inline ACEComplex operator*(const DOUBLE_TYPE &real, const ACEComplex &cm) {
    return cm * real;
}

/*
 * double + complex for commutativity with complex + double
 */
inline ACEComplex operator+(const DOUBLE_TYPE &real, const ACEComplex &cm) {
    return cm + real;
}

/**
A structure to store the derivative of \f$ Y_{lm} \f$
*/
struct ACEDYcomponent {
public:
    /**
    complex, contains the three components of derivative of Ylm,
    \f$ \frac{dY_{lm}}{dx}, \frac{dY_{lm}}{dy} and \frac{dY_{lm}}{dz}\f$
    */
    ACEComplex a[3];

    /*
     *  Overloading of arithmetical operator*= DOUBLE_TYPE
     */
    ACEDYcomponent &operator*=(const DOUBLE_TYPE &rhs) {
        this->a[0] *= rhs;
        this->a[1] *= rhs;
        this->a[2] *= rhs;
        return *this;
    }

    /*
     *  Overloading of arithmetical operator * ACEComplex
     */
    ACEDYcomponent operator*(const ACEComplex &rhs) const {
        ACEDYcomponent res;
        res.a[0] = this->a[0] * rhs;
        res.a[1] = this->a[1] * rhs;
        res.a[2] = this->a[2] * rhs;
        return res;
    }

    /*
     *  Overloading of arithmetical operator * DOUBLE_TYPE
     */
    ACEDYcomponent operator*(const DOUBLE_TYPE &rhs) const {
        ACEDYcomponent res;
        res.a[0] = this->a[0] * rhs;
        res.a[1] = this->a[1] * rhs;
        res.a[2] = this->a[2] * rhs;
        return res;
    }

    /*
     *  Return conjugated copy of itself
     */
    ACEDYcomponent conjugated() const {
        ACEDYcomponent res;
        res.a[0] = this->a[0].conjugated();
        res.a[1] = this->a[1].conjugated();
        res.a[2] = this->a[2].conjugated();
        return res;
    }

    /*
     *  Conjugated itself in-place
     */
    void conjugate() {
        this->a[0].conjugate();
        this->a[1].conjugate();
        this->a[2].conjugate();
    }

};


#endif //ACE_COMPLEX_H
