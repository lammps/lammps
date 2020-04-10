//
// Created by Yury Lysogorskiy on 26.02.20.
//
#ifndef ACE_COMPLEX_H
#define ACE_COMPLEX_H


/**
A custom data structure for complex numbers.
Associated functions can be added later if required.
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
Overloading of arithmetical operators
*/

    ACEComplex &operator+=(const ACEComplex &rhs) {
        this->real += rhs.real;
        this->img += rhs.img;
        return *this; // return the result by reference
    }

    ACEComplex &operator+=(const DOUBLE_TYPE &rhs) {
        this->real += rhs;
        return *this; // return the result by reference
    }


    ACEComplex &operator*=(const DOUBLE_TYPE &rhs) {
        this->real *= rhs;
        this->img *= rhs;
        return *this; // return the result by reference
    }

    ACEComplex &operator*=(const ACEComplex &rhs) {
        DOUBLE_TYPE old_real = this->real;
        this->real = old_real * rhs.real - this->img * rhs.img;
        this->img = old_real * rhs.img + this->img * rhs.real;
        return *this; // return the result by reference
    }

    ACEComplex operator*(const ACEComplex &cm2) const {
        ACEComplex res{real * cm2.real - img * cm2.img,
                       real * cm2.img + img * cm2.real};
        return res;
    }

    //complex conjugation
    ACEComplex conjugated() const {
        ACEComplex res{real, -img};
        return res;
    }

    void conjugate() {
        img = -img;
    }

    // real_part_product for real-part only multiplication
    DOUBLE_TYPE real_part_product(const ACEComplex &cm2) const {
        return real * cm2.real - img * cm2.img;
    }


    DOUBLE_TYPE operator%(const DOUBLE_TYPE &cm2) const {
        return real * cm2;
    }

    DOUBLE_TYPE real_part_product(const DOUBLE_TYPE &cm2) const {
        return real * cm2;
    }

    ACEComplex operator*(const DOUBLE_TYPE &cm2) const {
        ACEComplex res{real * cm2,
                       img * cm2};
        return res;
    }

    ACEComplex operator+(const ACEComplex &cm2) const {
        ACEComplex res{real + cm2.real,
                       img + cm2.img};
        return res;
    }

    ACEComplex operator+(const DOUBLE_TYPE &cm2) const {
        ACEComplex res{real + cm2, img};
        return res;
    }

    bool operator==(const ACEComplex &c2) const {
        return (real == c2.real) && (img == c2.img);
    }

    bool operator==(const DOUBLE_TYPE &d2) const {
        return (real == d2) && (img == 0.);
    }

    bool operator!=(const ACEComplex &c2) const {
        return (real != c2.real) || (img != c2.img);
    }

    bool operator!=(const DOUBLE_TYPE &d2) const {
        return (real != d2) || (img != 0.);
    }

};

// double * complex for commutativity with complex * double
inline ACEComplex operator*(const DOUBLE_TYPE &real, const ACEComplex &cm) {
    return cm * real;
}

// double + complex for commutativity with complex + double
inline ACEComplex operator+(const DOUBLE_TYPE &real, const ACEComplex &cm) {
    return cm + real;
}

/**
A structure to store the derivative of \f$ Y_{lm} \f$
*/
struct Dycomponent {
public:
    /**
    complex, contains the three components of derivative of Ylm,
    \f$ \frac{dY_{lm}}{dx}, \frac{dY_{lm}}{dy} and \frac{dY_{lm}}{dz}\f$
    */
    ACEComplex a[3];

    Dycomponent &operator*=(const DOUBLE_TYPE &rhs) {
        this->a[0] *= rhs;
        this->a[1] *= rhs;
        this->a[2] *= rhs;
        return *this;
    }

    Dycomponent operator*(const ACEComplex &rhs) const {
        Dycomponent res;
        res.a[0] = this->a[0] * rhs;
        res.a[1] = this->a[1] * rhs;
        res.a[2] = this->a[2] * rhs;
        return res;
    }

    Dycomponent operator*(const DOUBLE_TYPE &rhs) const {
        Dycomponent res;
        res.a[0] = this->a[0] * rhs;
        res.a[1] = this->a[1] * rhs;
        res.a[2] = this->a[2] * rhs;
        return res;
    }

    Dycomponent conjugated() const {
        Dycomponent res;
        res.a[0] = this->a[0].conjugated();
        res.a[1] = this->a[1].conjugated();
        res.a[2] = this->a[2].conjugated();
        return res;
    }

    void conjugate() {
        this->a[0].conjugate();
        this->a[1].conjugate();
        this->a[2].conjugate();
    }

};


#endif //ACE_COMPLEX_H
