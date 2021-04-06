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


// Created by Yury Lysogorskiy on 11.01.20.


#ifndef ACE_ARRAY2DLM_H
#define ACE_ARRAY2DLM_H

#include <stdexcept>
#include <string>

#include "ace_arraynd.h"
#include "ace_contigous_array.h"
#include "ace_types.h"

using namespace std;

/**
 * Contiguous array to organize values by \f$ (l,m) \f$ indiced (orbital moment and its projection).
 * Only \f$ l_\textrm{max}\f$ should be provided, \f$ m = -l, \dots,l \f$
 * for \f$ l = 0, \dots, l_\textrm{max}\f$
 * @tparam T type of values to store
 */
template<typename T>
class Array2DLM : public ContiguousArrayND<T> {
    using ContiguousArrayND<T>::array_name;
    using ContiguousArrayND<T>::data;
    using ContiguousArrayND<T>::size;

    LS_TYPE lmax = 0; ///< orbital dimension \f$ l_{max} \f$

    bool is_proxy = false; ///< flag to show, if object is owning the memory or just represent it (proxying)dimensions

public:
    /**
     *  Default empty constructor
     */
    Array2DLM() = default;

    /**
     *  Parametrized constructor
     * @param lmax maximum value of \f$ l \f$
     * @param array_name name of the array
     */
    explicit Array2DLM(LS_TYPE lmax, string array_name = "Array2DLM") {
        init(lmax, array_name);
    }


    /**
     * Constructor to create slices-proxy array, i.e. to provide access to the memory, but not to own it.
     * @param lmax  maximum value of \f$ l \f$
     * @param data_ptr pointer to original data
     * @param array_name name of the array
     */
    Array2DLM(LS_TYPE lmax, T *data_ptr, string array_name = "Array2DLM") {
        this->lmax = lmax;
        this->size = (lmax + 1) * (lmax + 1);
        this->data = data_ptr;
        this->array_name = array_name;
        is_proxy = true;
    };

    /**
     * Destructor
     */
    ~Array2DLM() {
        if (!is_proxy) {
            if (data != nullptr) delete[] data;
        }
        data = nullptr;
    }

    /**
     * Initialize array, allocate memory
     * @param lmax maximum value of l
     * @param array_name name of the array
     */
    void init(LS_TYPE lmax, string array_name = "Array2DLM") {
        if (is_proxy) {
            char s[1024];
            sprintf(s, "Could not re-initialize proxy-array %s\n", this->array_name.c_str());
            throw logic_error(s);
        }
        this->lmax = lmax;
        this->array_name = array_name;
        //for m = -l .. l
        if (size != (lmax + 1) * (lmax + 1)) {
            size = (lmax + 1) * (lmax + 1);
            if (data) delete[] data;
            data = new T[size]{};
            memset(data, 0.0, size * sizeof(T));
        } else {
            memset(data, 0, size * sizeof(T));
        }
    }

#ifdef MULTIARRAY_INDICES_CHECK
/**
 * Check if indices (l,m) are within array
 */
    void check_indices(LS_TYPE l, MS_TYPE m) const {

        if ((l < 0) | (l > lmax)) {
            fprintf(stderr, "%s: Index l = %d out of range (0, %d)\n", array_name.c_str(), l, lmax);
            exit(EXIT_FAILURE);
        }

        if ((m < -l) | (m > l)) {
            fprintf(stderr, "%s: Index m = %d out of range (%d, %d)\n", array_name.c_str(), m, -l, l);
            exit(EXIT_FAILURE);
        }
        size_t ii = l * (l + 1) + m;
        if (ii >= size) {
            fprintf(stderr, "%s: index = %ld out of range %ld\n", array_name.c_str(), ii, size);
            exit(EXIT_FAILURE);
        }
    }
#endif

    /**
     * Accessing the array value by index (l,m) for reading
     * @param l
     * @param m
     * @return array value
     */
    inline const T &operator()(LS_TYPE l, MS_TYPE m) const {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(l, m);
#endif
        //l^2 + l + m
        return data[l * (l + 1) + m];
    }

    /**
     * Accessing the array value by index (l,m) for writing
     * @param l
     * @param m
     * @return array value
     */
    inline T &operator()(LS_TYPE l, MS_TYPE m) {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(l, m);
#endif
        //l^2 + l + m
        return data[l * (l + 1) + m];
    }

    /**
     * Convert array to STL vector<vector<T>> container
     * @return  vector<vector<T>> container
     */
    vector<vector<T>> to_vector() const {
        vector<vector<T>> res;
        res.resize(lmax + 1);

        for (int i = 0; i < lmax + 1; i++) {
            res[i].resize(i + 1);
            for (int j = 0; j < i + 1; j++) {
                res[i][j] = operator()(i, j);
            }
        }
        return res;
    }
};

/**
 * Contiguous array to organize values by \f$ (i_0, l , m) \f$ indices.
 * Only \f$ d_{0}, l_\textrm{max}\f$ should be provided: \f$ m = -l, \dots,l \f$
 * for \f$ l = 0, \dots, l_\textrm{max}\f$
 * @tparam T type of values to store
 */
template<typename T>
class Array3DLM : public ContiguousArrayND<T> {
    using ContiguousArrayND<T>::array_name;
    using ContiguousArrayND<T>::data;
    using ContiguousArrayND<T>::size;

    LS_TYPE lmax = 0; ///< orbital dimension \f$ l_{max} \f$


    size_t dim[1] = {0}; ///< linear dimension \f$ d_{0} \f$

    size_t s[1] = {0}; ///< strides for linear dimensions

    Array1D<Array2DLM<T> *> _proxy_slices; ///< slices representation
public:
    /**
     *  Default empty constructor
     */
    Array3DLM() = default;

    /**
     *  Parametrized constructor
     * @param array_name name of the array
     */
    Array3DLM(string array_name) {
        this->array_name = array_name;
    };

    /**
     *  Parametrized constructor
     * @param d0 maximum value of \f$ i_0 \f$
     * @param lmax maximum value of \f$ l \f$
     * @param array_name name of the array
     */
    explicit Array3DLM(size_t d0, LS_TYPE lmax, string array_name = "Array3DLM") {
        init(d0, lmax, array_name);
    }

    /**
     * Initialize array and its slices
     * @param d0 maximum value of \f$ i_0 \f$
     * @param lmax  maximum value of \f$ l \f$
     * @param array_name name of the array
     */
    void init(size_t d0, LS_TYPE lmax, string array_name = "Array3DLM") {
        this->array_name = array_name;
        this->lmax = lmax;
        dim[0] = d0;
        s[0] = lmax * lmax;
        if (size != s[0] * dim[0]) {
            size = s[0] * dim[0];
            if (data) delete[] data;
            data = new T[size]{};
            memset(data, 0, size * sizeof(T));
        } else {
            memset(data, 0, size * sizeof(T));
        }

        _proxy_slices.set_array_name(array_name + "_proxy");
        //arrange proxy-slices
        _clear_proxies();
        _proxy_slices.resize(dim[0]);
        for (size_t i0 = 0; i0 < dim[0]; ++i0) {
            _proxy_slices(i0) = new Array2DLM<T>(this->lmax, &this->data[i0 * s[0]],
                                                 array_name + "_slice");
        }
    }

    /**
     * Release pointers to slices
     */
    void _clear_proxies() {
        for (size_t i0 = 0; i0 < _proxy_slices.get_dim(0); ++i0) {
            delete _proxy_slices(i0);
            _proxy_slices(i0) = nullptr;
        }
    }

    /**
     * Destructor, clear proxies
     */
    ~Array3DLM() {
        _clear_proxies();
    }

    /**
     * Resize array to new dimensions
     * @param d0
     * @param lmax
     */
    void resize(size_t d0, LS_TYPE lmax) {
        _clear_proxies();
        init(d0, lmax, this->array_name);
    }

    /**
     * Get array dimensions
     * @param d dimension index
     * @return  dimension along axis 'd'
     */
    size_t get_dim(int d) const {
        return dim[d];
    }

#ifdef MULTIARRAY_INDICES_CHECK
    /**
     * Check if indices (i0, l,m) are within array
     */
    void check_indices(size_t i0, LS_TYPE l, MS_TYPE m) const {
        if ((l < 0) | (l > lmax)) {
            fprintf(stderr, "%s: Index l = %d out of range (0, %d)\n", array_name.c_str(), l, lmax);
            exit(EXIT_FAILURE);
        }

        if ((m < -l) | (m > l)) {
            fprintf(stderr, "%s: Index m = %d out of range (%d, %d)\n", array_name.c_str(), m, -l, l);
            exit(EXIT_FAILURE);
        }

        if ((i0 < 0) | (i0 >= dim[0])) {
            fprintf(stderr, "%s: index i0 = %ld out of range (0, %ld)\n", array_name.c_str(), i0, dim[0] - 1);
            exit(EXIT_FAILURE);
        }

        size_t ii = i0 * s[0] + l * (l + 1) + m;
        if (ii >= size) {
            fprintf(stderr, "%s: index = %ld out of range %ld\n", array_name.c_str(), ii, size);
            exit(EXIT_FAILURE);
        }
    }
#endif

    /**
     * Accessing the array value by index (i0,l,m) for reading
     * @param i0
     * @param l
     * @param m
     * @return array value
     */
    inline const T &operator()(size_t i0, LS_TYPE l, MS_TYPE m) const {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, l, m);
#endif
        return data[i0 * s[0] + l * (l + 1) + m];
    }

    /**
     * Accessing the array value by index (i0,l,m) for writing
     * @param i0
     * @param l
     * @param m
     * @return array value
     */
    inline T &operator()(size_t i0, LS_TYPE l, MS_TYPE m) {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, l, m);
#endif
        return data[i0 * s[0] + l * (l + 1) + m];
    }

    /**
     * Return proxy Array2DLM pointing to i0, l=0, m=0 to read
     * @param i0
     * @return proxy Array2DLM pointing to i0, l=0, m=0
     */
    inline const Array2DLM<T> &operator()(size_t i0) const {
        return *_proxy_slices(i0);
    }

    /**
     * Return proxy Array2DLM pointing to i0, l=0, m=0 to write
     * @param i0
     * @return proxy Array2DLM pointing to i0, l=0, m=0
     */
    inline Array2DLM<T> &operator()(size_t i0) {
        return *_proxy_slices(i0);
    }
};


/**
 * Contiguous array to organize values by \f$ (i_0, i_1, l , m) \f$ indices.
 * Only \f$ d_{0}, d_{1}, l_\textrm{max}\f$ should be provided: \f$ m = -l, \dots,l \f$
 * for \f$ l = 0, \dots, l_\textrm{max}\f$
 * @tparam T type of values to store
 */
template<typename T>
class Array4DLM : public ContiguousArrayND<T> {
    using ContiguousArrayND<T>::array_name;
    using ContiguousArrayND<T>::data;
    using ContiguousArrayND<T>::size;

    LS_TYPE lmax = 0; ///< orbital dimension \f$ l_{max} \f$
    size_t dim[2] = {0, 0}; ///< linear dimension \f$ d_{0}, d_{1} \f$
    size_t s[2] = {0, 0}; ///< strides for linear dimensions

    Array2D<Array2DLM<T> *> _proxy_slices; ///< slices representation
public:
    /**
     *  Default empty constructor
     */
    Array4DLM() = default;

    /**
     *  Parametrized constructor
     * @param array_name name of the array
     */
    Array4DLM(string array_name) {
        this->array_name = array_name;
    };

    /**
     *  Parametrized constructor
     * @param d0 maximum value of \f$ i_0 \f$
     * @param d1 maximum value of \f$ i_1 \f$
     * @param lmax maximum value of \f$ l \f$
     * @param array_name name of the array
     */
    explicit Array4DLM(size_t d0, size_t d1, LS_TYPE lmax, string array_name = "Array4DLM") {
        init(d0, d1, lmax, array_name);
    }

    /**
     * Initialize array, reallocate memory and its slices
     * @param d0 maximum value of \f$ i_0 \f$
     * @param d1 maximum value of \f$ i_1 \f$
     * @param lmax  maximum value of \f$ l \f$
     * @param array_name name of the array
     */
    void init(size_t d0, size_t d1, LS_TYPE lmax, string array_name = "Array4DLM") {
        this->array_name = array_name;
        this->lmax = lmax;
        dim[1] = d1;
        dim[0] = d0;
        s[1] = lmax * lmax;
        s[0] = s[1] * dim[1];
        if (size != s[0] * dim[0]) {
            size = s[0] * dim[0];
            if (data) delete[] data;
            data = new T[size]{};
            memset(data, 0, size * sizeof(T));
        } else {
            memset(data, 0, size * sizeof(T));
        }

        _proxy_slices.set_array_name(array_name + "_proxy");
        //release old memory if there is any
        _clear_proxies();
        //arrange proxy-slices
        _proxy_slices.resize(dim[0], dim[1]);
        for (size_t i0 = 0; i0 < dim[0]; ++i0)
            for (size_t i1 = 0; i1 < dim[1]; ++i1) {
                _proxy_slices(i0, i1) = new Array2DLM<T>(this->lmax, &this->data[i0 * s[0] + i1 * s[1]],
                                                         array_name + "_slice");
            }
    }

    /**
     * Release pointers to slices
     */
    void _clear_proxies() {

        for (size_t i0 = 0; i0 < _proxy_slices.get_dim(0); ++i0)
            for (size_t i1 = 0; i1 < _proxy_slices.get_dim(1); ++i1) {
                delete _proxy_slices(i0, i1);
                _proxy_slices(i0, i1) = nullptr;
            }
    }

    /**
     * Destructor, clear proxies
     */
    ~Array4DLM() {
        _clear_proxies();
    }

    /**
     * Deallocate memory, reallocate with the new dimensions
     * @param d0
     * @param lmax
     */
    void resize(size_t d0, size_t d1, LS_TYPE lmax) {
        _clear_proxies();
        init(d0, d1, lmax, this->array_name);
    }

    /**
      * Get array dimensions
      * @param d dimension index
      * @return  dimension along axis 'd'
      */
    size_t get_dim(int d) const {
        return dim[d];
    }

#ifdef MULTIARRAY_INDICES_CHECK
    /**
     * Check if indices (i0, l,m) are within array
     */
    void check_indices(size_t i0, size_t i1, LS_TYPE l, MS_TYPE m) const {
        if ((l < 0) | (l > lmax)) {
            fprintf(stderr, "%s: Index l = %d out of range (0, %d)\n", array_name.c_str(), l, lmax);
            exit(EXIT_FAILURE);
        }

        if ((m < -l) | (m > l)) {
            fprintf(stderr, "%s: Index m = %d out of range (%d, %d)\n", array_name.c_str(), m, -l, l);
            exit(EXIT_FAILURE);
        }

        if ((i0 < 0) | (i0 >= dim[0])) {
            fprintf(stderr, "%s: index i0 = %ld out of range (0, %ld)\n", array_name.c_str(), i0, dim[0] - 1);
            exit(EXIT_FAILURE);
        }


        if ((i1 < 0) | (i1 >= dim[1])) {
            fprintf(stderr, "%s: index i1 = %ld out of range (0, %ld)\n", array_name.c_str(), i1, dim[1] - 1);
            exit(EXIT_FAILURE);
        }

        size_t ii = i0 * s[0] + i1 * s[1] + l * (l + 1) + m;
        if (ii >= size) {
            fprintf(stderr, "%s: index = %ld out of range %ld\n", array_name.c_str(), ii, size);
            exit(EXIT_FAILURE);
        }
    }
#endif

    /**
     * Accessing the array value by index (i0,l,m) for reading
     * @param i0
     * @param i1
     * @param l
     * @param m
     * @return array value
     */
    inline const T &operator()(size_t i0, size_t i1, LS_TYPE l, MS_TYPE m) const {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, i1, l, m);
#endif
        return data[i0 * s[0] + i1 * s[1] + l * (l + 1) + m];
    }

    /**
     * Accessing the array value by index (i0,l,m) for writing
     * @param i0
     * @param i1
     * @param l
     * @param m
     * @return array value
     */
    inline T &operator()(size_t i0, size_t i1, LS_TYPE l, MS_TYPE m) {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, i1, l, m);
#endif
        return data[i0 * s[0] + i1 * s[1] + l * (l + 1) + m];
    }

    /**
     * Return proxy Array2DLM pointing to i0, i1, l=0, m=0 to read
     * @param i0
     * @param i1
     * @return proxy Array2DLM pointing to i0, l=0, m=0
     */
    inline const Array2DLM<T> &operator()(size_t i0, size_t i1) const {
        return *_proxy_slices(i0, i1);
    }

    /**
     * Return proxy Array2DLM pointing to i0, i1, l=0, m=0 to write
     * @param i0
     * @param i1
     * @return proxy Array2DLM pointing to i0, l=0, m=0
     */
    inline Array2DLM<T> &operator()(size_t i0, size_t i1) {
        return *_proxy_slices(i0, i1);
    }
};

#endif //ACE_ARRAY2DLM_H