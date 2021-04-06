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


//automatically generate source code

#ifndef ACE_MULTIARRAY_H
#define ACE_MULTIARRAY_H

#include <cstring>
#include <vector>
#include <stdexcept>

#include "ace_contigous_array.h"

using namespace std;


/**
 * Multidimensional (1 - dimensional) array of type T with contiguous memory layout.
 * If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
 * be performed before accessing memory. By default this is turned off.
 * @tparam T data type
 */
template<class T>
class Array1D : public ContiguousArrayND<T> {
    using ContiguousArrayND<T>::array_name;
    using ContiguousArrayND<T>::data;
    using ContiguousArrayND<T>::size;
    using ContiguousArrayND<T>::is_proxy_;

    size_t dim[1] = {0}; ///< dimensions
    size_t s[1] = {0}; ///< strides
    int ndim = 1; ///< number of dimensions
public:

    /**
     * Default empty constructor
     */
    Array1D() = default;

    /**
     * Parametrized constructor with array name
     * @param array_name name of array (for error logging)
     */
    Array1D(const string &array_name) { this->array_name = array_name; }

    /**
    * Parametrized constructor
    * @param d0,... array sizes for different dimensions
    * @param array_name string name of the array
    */
    Array1D(size_t d0, const string &array_name = "Array1D", T *new_data = nullptr) {
        init(d0, array_name, new_data);
    }

    /**
    * Setup the dimensions and strides of array
    * @param d0,... array sizes for different dimensions
    */
    void set_dimensions(size_t d0) {

        dim[0] = d0;

        s[0] = 1;

        size = s[0] * dim[0];
    };

    /**
    * Initialize array storage, dimensions and strides
    * @param d0,... array sizes for different dimensions
    * @param array_name string name of the array
    */
    void init(size_t d0, const string &array_name = "Array1D", T *new_data = nullptr) {
        this->array_name = array_name;

        size_t old_size = size;
        set_dimensions(d0);

        bool new_is_proxy = (new_data != nullptr);

        if (new_is_proxy) {
            if (!is_proxy_) delete[] data;
            data = new_data;
        } else {
            if (size != old_size) {
                T *old_data = data; //preserve the pointer to the old data
                data = new T[size]; // allocate new data
                if (old_data != nullptr) { //
                    size_t min_size = old_size < size ? old_size : size;
                    memcpy(data, old_data, min_size * sizeof(T));
                    if (!is_proxy_) delete[] old_data;
                }
                //memset(data, 0, size * sizeof(T));
            } else {
                //memset(data, 0, size * sizeof(T));
            }
        }

        is_proxy_ = new_is_proxy;
    }

    /**
    * Resize array
    * @param d0,... array sizes for different dimensions
    */
    void resize(size_t d0) {
        init(d0, this->array_name);
    }

    /**
    * Reshape array without reset the data
    * @param d0,... array sizes for different dimensions
    */
    void reshape(size_t d0) {
        //check data size consistency
        size_t new_size = d0;
        if (new_size != size)
            throw invalid_argument("Couldn't reshape array when the size is not conserved");
        set_dimensions(d0);
    }

    /**
    * Get array size in dimension "d"
    * @param d dimension index
    */
    size_t get_dim(int d) const {
        return dim[d];
    }

#ifdef MULTIARRAY_INDICES_CHECK

    /**
    * Check indices validity. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    void check_indices(size_t i0) const {

        if ((i0 < 0) | (i0 >= dim[0])) {
            char buf[1024];
            sprintf(buf, "%s: index i0=%ld out of range (0, %ld)\n", array_name.c_str(), i0, dim[0] - 1);
            throw std::out_of_range(buf);
        }

    }

#endif

    /**
    * Array access operator() for reading. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    inline const T &operator()(size_t i0) const {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0);
#endif

        return data[i0];

    }

    /**
    * Array access operator() for writing. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    inline T &operator()(size_t i0) {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0);
#endif

        return data[i0];

    }

    /**
    * Array comparison operator
    * @param other
    */
    bool operator==(const Array1D &other) const {
        //compare dimensions
        for (int d = 0; d < ndim; d++) {
            if (this->dim[d] != other.dim[d])
                return false;
        }
        return ContiguousArrayND<T>::operator==(other);
    }

    /**
    * Convert to nested vector<vector<...<T>> container
    * @return vector container
    */
    vector<T> to_vector() const {
        vector<T> res;

        res.resize(dim[0]);

        for (int i0 = 0; i0 < dim[0]; i0++) {
            res[i0] = operator()(i0);

        }


        return res;
    } // end to_vector()


    /**
    * Set values to vector<vector<...<T>> container
    * @param vec container
    */
    void set_vector(const vector<T> &vec) {
        size_t d0 = 0;
        d0 = vec.size();


        init(d0, array_name);
        for (int i0 = 0; i0 < dim[0]; i0++) {
            operator()(i0) = vec.at(i0);

        }

    }

    /**
    * Parametrized constructor from  vector<vector<...<T>> container
    * @param vec container
    * @param array_name array name
    */
    Array1D(const vector<T> &vec, const string &array_name = "Array1D") {
        this->set_vector(vec);
        this->array_name = array_name;
    }

    /**
    * operator= to vector<vector<...<T>> container
    * @param vec container
    */
    Array1D &operator=(const vector<T> &vec) {
        this->set_vector(vec);
        return *this;
    }


    vector<size_t> get_shape() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = dim[d];
        return sh;
    }

    vector<size_t> get_strides() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = s[d];
        return sh;
    }

    vector<size_t> get_memory_strides() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = s[d] * sizeof(T);
        return sh;
    }
};


/**
 * Multidimensional (2 - dimensional) array of type T with contiguous memory layout.
 * If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
 * be performed before accessing memory. By default this is turned off.
 * @tparam T data type
 */
template<class T>
class Array2D : public ContiguousArrayND<T> {
    using ContiguousArrayND<T>::array_name;
    using ContiguousArrayND<T>::data;
    using ContiguousArrayND<T>::size;
    using ContiguousArrayND<T>::is_proxy_;

    size_t dim[2] = {0}; ///< dimensions
    size_t s[2] = {0}; ///< strides
    int ndim = 2; ///< number of dimensions
public:

    /**
     * Default empty constructor
     */
    Array2D() = default;

    /**
     * Parametrized constructor with array name
     * @param array_name name of array (for error logging)
     */
    Array2D(const string &array_name) { this->array_name = array_name; }

    /**
    * Parametrized constructor
    * @param d0,... array sizes for different dimensions
    * @param array_name string name of the array
    */
    Array2D(size_t d0, size_t d1, const string &array_name = "Array2D", T *new_data = nullptr) {
        init(d0, d1, array_name, new_data);
    }

    /**
    * Setup the dimensions and strides of array
    * @param d0,... array sizes for different dimensions
    */
    void set_dimensions(size_t d0, size_t d1) {

        dim[0] = d0;
        dim[1] = d1;

        s[1] = 1;
        s[0] = s[1] * dim[1];

        size = s[0] * dim[0];
    };

    /**
    * Initialize array storage, dimensions and strides
    * @param d0,... array sizes for different dimensions
    * @param array_name string name of the array
    */
    void init(size_t d0, size_t d1, const string &array_name = "Array2D", T *new_data = nullptr) {
        this->array_name = array_name;

        size_t old_size = size;
        set_dimensions(d0, d1);

        bool new_is_proxy = (new_data != nullptr);

        if (new_is_proxy) {
            if (!is_proxy_) delete[] data;
            data = new_data;
        } else {
            if (size != old_size) {
                T *old_data = data; //preserve the pointer to the old data
                data = new T[size]; // allocate new data
                if (old_data != nullptr) { //
                    size_t min_size = old_size < size ? old_size : size;
                    memcpy(data, old_data, min_size * sizeof(T));
                    if (!is_proxy_) delete[] old_data;
                }
                //memset(data, 0, size * sizeof(T));
            } else {
                //memset(data, 0, size * sizeof(T));
            }
        }

        is_proxy_ = new_is_proxy;
    }

    /**
    * Resize array
    * @param d0,... array sizes for different dimensions
    */
    void resize(size_t d0, size_t d1) {
        init(d0, d1, this->array_name);
    }

    /**
    * Reshape array without reset the data
    * @param d0,... array sizes for different dimensions
    */
    void reshape(size_t d0, size_t d1) {
        //check data size consistency
        size_t new_size = d0 * d1;
        if (new_size != size)
            throw invalid_argument("Couldn't reshape array when the size is not conserved");
        set_dimensions(d0, d1);
    }

    /**
    * Get array size in dimension "d"
    * @param d dimension index
    */
    size_t get_dim(int d) const {
        return dim[d];
    }

#ifdef MULTIARRAY_INDICES_CHECK

    /**
    * Check indices validity. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    void check_indices(size_t i0, size_t i1) const {

        if ((i0 < 0) | (i0 >= dim[0])) {
            char buf[1024];
            sprintf(buf, "%s: index i0=%ld out of range (0, %ld)\n", array_name.c_str(), i0, dim[0] - 1);
            throw std::out_of_range(buf);
        }

        if ((i1 < 0) | (i1 >= dim[1])) {
            char buf[1024];
            sprintf(buf, "%s: index i1=%ld out of range (0, %ld)\n", array_name.c_str(), i1, dim[1] - 1);
            throw std::out_of_range(buf);
        }

    }

#endif

    /**
    * Array access operator() for reading. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    inline const T &operator()(size_t i0, size_t i1) const {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, i1);
#endif

        return data[i0 * s[0] + i1];

    }

    /**
    * Array access operator() for writing. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    inline T &operator()(size_t i0, size_t i1) {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, i1);
#endif

        return data[i0 * s[0] + i1];

    }

    /**
    * Array comparison operator
    * @param other
    */
    bool operator==(const Array2D &other) const {
        //compare dimensions
        for (int d = 0; d < ndim; d++) {
            if (this->dim[d] != other.dim[d])
                return false;
        }
        return ContiguousArrayND<T>::operator==(other);
    }

    /**
    * Convert to nested vector<vector<...<T>> container
    * @return vector container
    */
    vector<vector<T>> to_vector() const {
        vector<vector<T>> res;

        res.resize(dim[0]);

        for (int i0 = 0; i0 < dim[0]; i0++) {
            res[i0].resize(dim[1]);


            for (int i1 = 0; i1 < dim[1]; i1++) {
                res[i0][i1] = operator()(i0, i1);

            }
        }


        return res;
    } // end to_vector()


    /**
    * Set values to vector<vector<...<T>> container
    * @param vec container
    */
    void set_vector(const vector<vector<T>> &vec) {
        size_t d0 = 0;
        size_t d1 = 0;
        d0 = vec.size();

        if (d0 > 0) {
            d1 = vec.at(0).size();


        }

        init(d0, d1, array_name);
        for (int i0 = 0; i0 < dim[0]; i0++) {
            if (vec.at(i0).size() != d1)
                throw std::invalid_argument("Vector size is not constant at dimension 1");

            for (int i1 = 0; i1 < dim[1]; i1++) {
                operator()(i0, i1) = vec.at(i0).at(i1);

            }
        }

    }

    /**
    * Parametrized constructor from  vector<vector<...<T>> container
    * @param vec container
    * @param array_name array name
    */
    Array2D(const vector<vector<T>> &vec, const string &array_name = "Array2D") {
        this->set_vector(vec);
        this->array_name = array_name;
    }

    /**
    * operator= to vector<vector<...<T>> container
    * @param vec container
    */
    Array2D &operator=(const vector<vector<T>> &vec) {
        this->set_vector(vec);
        return *this;
    }


    /**
    * operator= to flatten vector<T> container
    * @param vec container
    */
    Array2D &operator=(const vector<T> &vec) {
        this->set_flatten_vector(vec);
        return *this;
    }


    vector<size_t> get_shape() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = dim[d];
        return sh;
    }

    vector<size_t> get_strides() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = s[d];
        return sh;
    }

    vector<size_t> get_memory_strides() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = s[d] * sizeof(T);
        return sh;
    }
};


/**
 * Multidimensional (3 - dimensional) array of type T with contiguous memory layout.
 * If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
 * be performed before accessing memory. By default this is turned off.
 * @tparam T data type
 */
template<class T>
class Array3D : public ContiguousArrayND<T> {
    using ContiguousArrayND<T>::array_name;
    using ContiguousArrayND<T>::data;
    using ContiguousArrayND<T>::size;
    using ContiguousArrayND<T>::is_proxy_;

    size_t dim[3] = {0}; ///< dimensions
    size_t s[3] = {0}; ///< strides
    int ndim = 3; ///< number of dimensions
public:

    /**
     * Default empty constructor
     */
    Array3D() = default;

    /**
     * Parametrized constructor with array name
     * @param array_name name of array (for error logging)
     */
    Array3D(const string &array_name) { this->array_name = array_name; }

    /**
    * Parametrized constructor
    * @param d0,... array sizes for different dimensions
    * @param array_name string name of the array
    */
    Array3D(size_t d0, size_t d1, size_t d2, const string &array_name = "Array3D", T *new_data = nullptr) {
        init(d0, d1, d2, array_name, new_data);
    }

    /**
    * Setup the dimensions and strides of array
    * @param d0,... array sizes for different dimensions
    */
    void set_dimensions(size_t d0, size_t d1, size_t d2) {

        dim[0] = d0;
        dim[1] = d1;
        dim[2] = d2;

        s[2] = 1;
        s[1] = s[2] * dim[2];
        s[0] = s[1] * dim[1];

        size = s[0] * dim[0];
    };

    /**
    * Initialize array storage, dimensions and strides
    * @param d0,... array sizes for different dimensions
    * @param array_name string name of the array
    */
    void init(size_t d0, size_t d1, size_t d2, const string &array_name = "Array3D", T *new_data = nullptr) {
        this->array_name = array_name;

        size_t old_size = size;
        set_dimensions(d0, d1, d2);

        bool new_is_proxy = (new_data != nullptr);

        if (new_is_proxy) {
            if (!is_proxy_) delete[] data;
            data = new_data;
        } else {
            if (size != old_size) {
                T *old_data = data; //preserve the pointer to the old data
                data = new T[size]; // allocate new data
                if (old_data != nullptr) { //
                    size_t min_size = old_size < size ? old_size : size;
                    memcpy(data, old_data, min_size * sizeof(T));
                    if (!is_proxy_) delete[] old_data;
                }
                //memset(data, 0, size * sizeof(T));
            } else {
                //memset(data, 0, size * sizeof(T));
            }
        }

        is_proxy_ = new_is_proxy;
    }

    /**
    * Resize array
    * @param d0,... array sizes for different dimensions
    */
    void resize(size_t d0, size_t d1, size_t d2) {
        init(d0, d1, d2, this->array_name);
    }

    /**
    * Reshape array without reset the data
    * @param d0,... array sizes for different dimensions
    */
    void reshape(size_t d0, size_t d1, size_t d2) {
        //check data size consistency
        size_t new_size = d0 * d1 * d2;
        if (new_size != size)
            throw invalid_argument("Couldn't reshape array when the size is not conserved");
        set_dimensions(d0, d1, d2);
    }

    /**
    * Get array size in dimension "d"
    * @param d dimension index
    */
    size_t get_dim(int d) const {
        return dim[d];
    }

#ifdef MULTIARRAY_INDICES_CHECK

    /**
    * Check indices validity. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    void check_indices(size_t i0, size_t i1, size_t i2) const {

        if ((i0 < 0) | (i0 >= dim[0])) {
            char buf[1024];
            sprintf(buf, "%s: index i0=%ld out of range (0, %ld)\n", array_name.c_str(), i0, dim[0] - 1);
            throw std::out_of_range(buf);
        }

        if ((i1 < 0) | (i1 >= dim[1])) {
            char buf[1024];
            sprintf(buf, "%s: index i1=%ld out of range (0, %ld)\n", array_name.c_str(), i1, dim[1] - 1);
            throw std::out_of_range(buf);
        }

        if ((i2 < 0) | (i2 >= dim[2])) {
            char buf[1024];
            sprintf(buf, "%s: index i2=%ld out of range (0, %ld)\n", array_name.c_str(), i2, dim[2] - 1);
            throw std::out_of_range(buf);
        }

    }

#endif

    /**
    * Array access operator() for reading. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    inline const T &operator()(size_t i0, size_t i1, size_t i2) const {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, i1, i2);
#endif

        return data[i0 * s[0] + i1 * s[1] + i2];

    }

    /**
    * Array access operator() for writing. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    inline T &operator()(size_t i0, size_t i1, size_t i2) {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, i1, i2);
#endif

        return data[i0 * s[0] + i1 * s[1] + i2];

    }

    /**
    * Array comparison operator
    * @param other
    */
    bool operator==(const Array3D &other) const {
        //compare dimensions
        for (int d = 0; d < ndim; d++) {
            if (this->dim[d] != other.dim[d])
                return false;
        }
        return ContiguousArrayND<T>::operator==(other);
    }

    /**
    * Convert to nested vector<vector<...<T>> container
    * @return vector container
    */
    vector<vector<vector<T>>> to_vector() const {
        vector<vector<vector<T>>> res;

        res.resize(dim[0]);

        for (int i0 = 0; i0 < dim[0]; i0++) {
            res[i0].resize(dim[1]);


            for (int i1 = 0; i1 < dim[1]; i1++) {
                res[i0][i1].resize(dim[2]);


                for (int i2 = 0; i2 < dim[2]; i2++) {
                    res[i0][i1][i2] = operator()(i0, i1, i2);

                }
            }
        }


        return res;
    } // end to_vector()


    /**
    * Set values to vector<vector<...<T>> container
    * @param vec container
    */
    void set_vector(const vector<vector<vector<T>>> &vec) {
        size_t d0 = 0;
        size_t d1 = 0;
        size_t d2 = 0;
        d0 = vec.size();

        if (d0 > 0) {
            d1 = vec.at(0).size();
            if (d1 > 0) {
                d2 = vec.at(0).at(0).size();


            }
        }

        init(d0, d1, d2, array_name);
        for (int i0 = 0; i0 < dim[0]; i0++) {
            if (vec.at(i0).size() != d1)
                throw std::invalid_argument("Vector size is not constant at dimension 1");

            for (int i1 = 0; i1 < dim[1]; i1++) {
                if (vec.at(i0).at(i1).size() != d2)
                    throw std::invalid_argument("Vector size is not constant at dimension 2");

                for (int i2 = 0; i2 < dim[2]; i2++) {
                    operator()(i0, i1, i2) = vec.at(i0).at(i1).at(i2);

                }
            }
        }

    }

    /**
    * Parametrized constructor from  vector<vector<...<T>> container
    * @param vec container
    * @param array_name array name
    */
    Array3D(const vector<vector<vector<T>>> &vec, const string &array_name = "Array3D") {
        this->set_vector(vec);
        this->array_name = array_name;
    }

    /**
    * operator= to vector<vector<...<T>> container
    * @param vec container
    */
    Array3D &operator=(const vector<vector<vector<T>>> &vec) {
        this->set_vector(vec);
        return *this;
    }


    /**
    * operator= to flatten vector<T> container
    * @param vec container
    */
    Array3D &operator=(const vector<T> &vec) {
        this->set_flatten_vector(vec);
        return *this;
    }


    vector<size_t> get_shape() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = dim[d];
        return sh;
    }

    vector<size_t> get_strides() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = s[d];
        return sh;
    }

    vector<size_t> get_memory_strides() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = s[d] * sizeof(T);
        return sh;
    }
};


/**
 * Multidimensional (4 - dimensional) array of type T with contiguous memory layout.
 * If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
 * be performed before accessing memory. By default this is turned off.
 * @tparam T data type
 */
template<class T>
class Array4D : public ContiguousArrayND<T> {
    using ContiguousArrayND<T>::array_name;
    using ContiguousArrayND<T>::data;
    using ContiguousArrayND<T>::size;
    using ContiguousArrayND<T>::is_proxy_;

    size_t dim[4] = {0}; ///< dimensions
    size_t s[4] = {0}; ///< strides
    int ndim = 4; ///< number of dimensions
public:

    /**
     * Default empty constructor
     */
    Array4D() = default;

    /**
     * Parametrized constructor with array name
     * @param array_name name of array (for error logging)
     */
    Array4D(const string &array_name) { this->array_name = array_name; }

    /**
    * Parametrized constructor
    * @param d0,... array sizes for different dimensions
    * @param array_name string name of the array
    */
    Array4D(size_t d0, size_t d1, size_t d2, size_t d3, const string &array_name = "Array4D", T *new_data = nullptr) {
        init(d0, d1, d2, d3, array_name, new_data);
    }

    /**
    * Setup the dimensions and strides of array
    * @param d0,... array sizes for different dimensions
    */
    void set_dimensions(size_t d0, size_t d1, size_t d2, size_t d3) {

        dim[0] = d0;
        dim[1] = d1;
        dim[2] = d2;
        dim[3] = d3;

        s[3] = 1;
        s[2] = s[3] * dim[3];
        s[1] = s[2] * dim[2];
        s[0] = s[1] * dim[1];

        size = s[0] * dim[0];
    };

    /**
    * Initialize array storage, dimensions and strides
    * @param d0,... array sizes for different dimensions
    * @param array_name string name of the array
    */
    void init(size_t d0, size_t d1, size_t d2, size_t d3, const string &array_name = "Array4D", T *new_data = nullptr) {
        this->array_name = array_name;

        size_t old_size = size;
        set_dimensions(d0, d1, d2, d3);

        bool new_is_proxy = (new_data != nullptr);

        if (new_is_proxy) {
            if (!is_proxy_) delete[] data;
            data = new_data;
        } else {
            if (size != old_size) {
                T *old_data = data; //preserve the pointer to the old data
                data = new T[size]; // allocate new data
                if (old_data != nullptr) { //
                    size_t min_size = old_size < size ? old_size : size;
                    memcpy(data, old_data, min_size * sizeof(T));
                    if (!is_proxy_) delete[] old_data;
                }
                //memset(data, 0, size * sizeof(T));
            } else {
                //memset(data, 0, size * sizeof(T));
            }
        }

        is_proxy_ = new_is_proxy;
    }

    /**
    * Resize array
    * @param d0,... array sizes for different dimensions
    */
    void resize(size_t d0, size_t d1, size_t d2, size_t d3) {
        init(d0, d1, d2, d3, this->array_name);
    }

    /**
    * Reshape array without reset the data
    * @param d0,... array sizes for different dimensions
    */
    void reshape(size_t d0, size_t d1, size_t d2, size_t d3) {
        //check data size consistency
        size_t new_size = d0 * d1 * d2 * d3;
        if (new_size != size)
            throw invalid_argument("Couldn't reshape array when the size is not conserved");
        set_dimensions(d0, d1, d2, d3);
    }

    /**
    * Get array size in dimension "d"
    * @param d dimension index
    */
    size_t get_dim(int d) const {
        return dim[d];
    }

#ifdef MULTIARRAY_INDICES_CHECK

    /**
    * Check indices validity. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    void check_indices(size_t i0, size_t i1, size_t i2, size_t i3) const {

        if ((i0 < 0) | (i0 >= dim[0])) {
            char buf[1024];
            sprintf(buf, "%s: index i0=%ld out of range (0, %ld)\n", array_name.c_str(), i0, dim[0] - 1);
            throw std::out_of_range(buf);
        }

        if ((i1 < 0) | (i1 >= dim[1])) {
            char buf[1024];
            sprintf(buf, "%s: index i1=%ld out of range (0, %ld)\n", array_name.c_str(), i1, dim[1] - 1);
            throw std::out_of_range(buf);
        }

        if ((i2 < 0) | (i2 >= dim[2])) {
            char buf[1024];
            sprintf(buf, "%s: index i2=%ld out of range (0, %ld)\n", array_name.c_str(), i2, dim[2] - 1);
            throw std::out_of_range(buf);
        }

        if ((i3 < 0) | (i3 >= dim[3])) {
            char buf[1024];
            sprintf(buf, "%s: index i3=%ld out of range (0, %ld)\n", array_name.c_str(), i3, dim[3] - 1);
            throw std::out_of_range(buf);
        }

    }

#endif

    /**
    * Array access operator() for reading. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    inline const T &operator()(size_t i0, size_t i1, size_t i2, size_t i3) const {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, i1, i2, i3);
#endif

        return data[i0 * s[0] + i1 * s[1] + i2 * s[2] + i3];

    }

    /**
    * Array access operator() for writing. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    inline T &operator()(size_t i0, size_t i1, size_t i2, size_t i3) {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, i1, i2, i3);
#endif

        return data[i0 * s[0] + i1 * s[1] + i2 * s[2] + i3];

    }

    /**
    * Array comparison operator
    * @param other
    */
    bool operator==(const Array4D &other) const {
        //compare dimensions
        for (int d = 0; d < ndim; d++) {
            if (this->dim[d] != other.dim[d])
                return false;
        }
        return ContiguousArrayND<T>::operator==(other);
    }

    /**
    * Convert to nested vector<vector<...<T>> container
    * @return vector container
    */
    vector<vector<vector<vector<T>>>> to_vector() const {
        vector<vector<vector<vector<T>>>> res;

        res.resize(dim[0]);

        for (int i0 = 0; i0 < dim[0]; i0++) {
            res[i0].resize(dim[1]);


            for (int i1 = 0; i1 < dim[1]; i1++) {
                res[i0][i1].resize(dim[2]);


                for (int i2 = 0; i2 < dim[2]; i2++) {
                    res[i0][i1][i2].resize(dim[3]);


                    for (int i3 = 0; i3 < dim[3]; i3++) {
                        res[i0][i1][i2][i3] = operator()(i0, i1, i2, i3);

                    }
                }
            }
        }


        return res;
    } // end to_vector()


    /**
    * Set values to vector<vector<...<T>> container
    * @param vec container
    */
    void set_vector(const vector<vector<vector<vector<T>>>> &vec) {
        size_t d0 = 0;
        size_t d1 = 0;
        size_t d2 = 0;
        size_t d3 = 0;
        d0 = vec.size();

        if (d0 > 0) {
            d1 = vec.at(0).size();
            if (d1 > 0) {
                d2 = vec.at(0).at(0).size();
                if (d2 > 0) {
                    d3 = vec.at(0).at(0).at(0).size();


                }
            }
        }

        init(d0, d1, d2, d3, array_name);
        for (int i0 = 0; i0 < dim[0]; i0++) {
            if (vec.at(i0).size() != d1)
                throw std::invalid_argument("Vector size is not constant at dimension 1");

            for (int i1 = 0; i1 < dim[1]; i1++) {
                if (vec.at(i0).at(i1).size() != d2)
                    throw std::invalid_argument("Vector size is not constant at dimension 2");

                for (int i2 = 0; i2 < dim[2]; i2++) {
                    if (vec.at(i0).at(i1).at(i2).size() != d3)
                        throw std::invalid_argument("Vector size is not constant at dimension 3");

                    for (int i3 = 0; i3 < dim[3]; i3++) {
                        operator()(i0, i1, i2, i3) = vec.at(i0).at(i1).at(i2).at(i3);

                    }
                }
            }
        }

    }

    /**
    * Parametrized constructor from  vector<vector<...<T>> container
    * @param vec container
    * @param array_name array name
    */
    Array4D(const vector<vector<vector<vector<T>>>> &vec, const string &array_name = "Array4D") {
        this->set_vector(vec);
        this->array_name = array_name;
    }

    /**
    * operator= to vector<vector<...<T>> container
    * @param vec container
    */
    Array4D &operator=(const vector<vector<vector<vector<T>>>> &vec) {
        this->set_vector(vec);
        return *this;
    }


    /**
    * operator= to flatten vector<T> container
    * @param vec container
    */
    Array4D &operator=(const vector<T> &vec) {
        this->set_flatten_vector(vec);
        return *this;
    }


    vector<size_t> get_shape() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = dim[d];
        return sh;
    }

    vector<size_t> get_strides() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = s[d];
        return sh;
    }

    vector<size_t> get_memory_strides() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = s[d] * sizeof(T);
        return sh;
    }
};


/**
 * Multidimensional (5 - dimensional) array of type T with contiguous memory layout.
 * If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
 * be performed before accessing memory. By default this is turned off.
 * @tparam T data type
 */
template<class T>
class Array5D : public ContiguousArrayND<T> {
    using ContiguousArrayND<T>::array_name;
    using ContiguousArrayND<T>::data;
    using ContiguousArrayND<T>::size;
    using ContiguousArrayND<T>::is_proxy_;

    size_t dim[5] = {0}; ///< dimensions
    size_t s[5] = {0}; ///< strides
    int ndim = 5; ///< number of dimensions
public:

    /**
     * Default empty constructor
     */
    Array5D() = default;

    /**
     * Parametrized constructor with array name
     * @param array_name name of array (for error logging)
     */
    Array5D(const string &array_name) { this->array_name = array_name; }

    /**
    * Parametrized constructor
    * @param d0,... array sizes for different dimensions
    * @param array_name string name of the array
    */
    Array5D(size_t d0, size_t d1, size_t d2, size_t d3, size_t d4, const string &array_name = "Array5D",
            T *new_data = nullptr) {
        init(d0, d1, d2, d3, d4, array_name, new_data);
    }

    /**
    * Setup the dimensions and strides of array
    * @param d0,... array sizes for different dimensions
    */
    void set_dimensions(size_t d0, size_t d1, size_t d2, size_t d3, size_t d4) {

        dim[0] = d0;
        dim[1] = d1;
        dim[2] = d2;
        dim[3] = d3;
        dim[4] = d4;

        s[4] = 1;
        s[3] = s[4] * dim[4];
        s[2] = s[3] * dim[3];
        s[1] = s[2] * dim[2];
        s[0] = s[1] * dim[1];

        size = s[0] * dim[0];
    };

    /**
    * Initialize array storage, dimensions and strides
    * @param d0,... array sizes for different dimensions
    * @param array_name string name of the array
    */
    void init(size_t d0, size_t d1, size_t d2, size_t d3, size_t d4, const string &array_name = "Array5D",
              T *new_data = nullptr) {
        this->array_name = array_name;

        size_t old_size = size;
        set_dimensions(d0, d1, d2, d3, d4);

        bool new_is_proxy = (new_data != nullptr);

        if (new_is_proxy) {
            if (!is_proxy_) delete[] data;
            data = new_data;
        } else {
            if (size != old_size) {
                T *old_data = data; //preserve the pointer to the old data
                data = new T[size]; // allocate new data
                if (old_data != nullptr) { //
                    size_t min_size = old_size < size ? old_size : size;
                    memcpy(data, old_data, min_size * sizeof(T));
                    if (!is_proxy_) delete[] old_data;
                }
                //memset(data, 0, size * sizeof(T));
            } else {
                //memset(data, 0, size * sizeof(T));
            }
        }

        is_proxy_ = new_is_proxy;
    }

    /**
    * Resize array
    * @param d0,... array sizes for different dimensions
    */
    void resize(size_t d0, size_t d1, size_t d2, size_t d3, size_t d4) {
        init(d0, d1, d2, d3, d4, this->array_name);
    }

    /**
    * Reshape array without reset the data
    * @param d0,... array sizes for different dimensions
    */
    void reshape(size_t d0, size_t d1, size_t d2, size_t d3, size_t d4) {
        //check data size consistency
        size_t new_size = d0 * d1 * d2 * d3 * d4;
        if (new_size != size)
            throw invalid_argument("Couldn't reshape array when the size is not conserved");
        set_dimensions(d0, d1, d2, d3, d4);
    }

    /**
    * Get array size in dimension "d"
    * @param d dimension index
    */
    size_t get_dim(int d) const {
        return dim[d];
    }

#ifdef MULTIARRAY_INDICES_CHECK

    /**
    * Check indices validity. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    void check_indices(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4) const {

        if ((i0 < 0) | (i0 >= dim[0])) {
            char buf[1024];
            sprintf(buf, "%s: index i0=%ld out of range (0, %ld)\n", array_name.c_str(), i0, dim[0] - 1);
            throw std::out_of_range(buf);
        }

        if ((i1 < 0) | (i1 >= dim[1])) {
            char buf[1024];
            sprintf(buf, "%s: index i1=%ld out of range (0, %ld)\n", array_name.c_str(), i1, dim[1] - 1);
            throw std::out_of_range(buf);
        }

        if ((i2 < 0) | (i2 >= dim[2])) {
            char buf[1024];
            sprintf(buf, "%s: index i2=%ld out of range (0, %ld)\n", array_name.c_str(), i2, dim[2] - 1);
            throw std::out_of_range(buf);
        }

        if ((i3 < 0) | (i3 >= dim[3])) {
            char buf[1024];
            sprintf(buf, "%s: index i3=%ld out of range (0, %ld)\n", array_name.c_str(), i3, dim[3] - 1);
            throw std::out_of_range(buf);
        }

        if ((i4 < 0) | (i4 >= dim[4])) {
            char buf[1024];
            sprintf(buf, "%s: index i4=%ld out of range (0, %ld)\n", array_name.c_str(), i4, dim[4] - 1);
            throw std::out_of_range(buf);
        }

    }

#endif

    /**
    * Array access operator() for reading. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    inline const T &operator()(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4) const {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, i1, i2, i3, i4);
#endif

        return data[i0 * s[0] + i1 * s[1] + i2 * s[2] + i3 * s[3] + i4];

    }

    /**
    * Array access operator() for writing. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    inline T &operator()(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4) {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, i1, i2, i3, i4);
#endif

        return data[i0 * s[0] + i1 * s[1] + i2 * s[2] + i3 * s[3] + i4];

    }

    /**
    * Array comparison operator
    * @param other
    */
    bool operator==(const Array5D &other) const {
        //compare dimensions
        for (int d = 0; d < ndim; d++) {
            if (this->dim[d] != other.dim[d])
                return false;
        }
        return ContiguousArrayND<T>::operator==(other);
    }

    /**
    * Convert to nested vector<vector<...<T>> container
    * @return vector container
    */
    vector<vector<vector<vector<vector<T>>>>> to_vector() const {
        vector<vector<vector<vector<vector<T>>>>> res;

        res.resize(dim[0]);

        for (int i0 = 0; i0 < dim[0]; i0++) {
            res[i0].resize(dim[1]);


            for (int i1 = 0; i1 < dim[1]; i1++) {
                res[i0][i1].resize(dim[2]);


                for (int i2 = 0; i2 < dim[2]; i2++) {
                    res[i0][i1][i2].resize(dim[3]);


                    for (int i3 = 0; i3 < dim[3]; i3++) {
                        res[i0][i1][i2][i3].resize(dim[4]);


                        for (int i4 = 0; i4 < dim[4]; i4++) {
                            res[i0][i1][i2][i3][i4] = operator()(i0, i1, i2, i3, i4);

                        }
                    }
                }
            }
        }


        return res;
    } // end to_vector()


    /**
    * Set values to vector<vector<...<T>> container
    * @param vec container
    */
    void set_vector(const vector<vector<vector<vector<vector<T>>>>> &vec) {
        size_t d0 = 0;
        size_t d1 = 0;
        size_t d2 = 0;
        size_t d3 = 0;
        size_t d4 = 0;
        d0 = vec.size();

        if (d0 > 0) {
            d1 = vec.at(0).size();
            if (d1 > 0) {
                d2 = vec.at(0).at(0).size();
                if (d2 > 0) {
                    d3 = vec.at(0).at(0).at(0).size();
                    if (d3 > 0) {
                        d4 = vec.at(0).at(0).at(0).at(0).size();


                    }
                }
            }
        }

        init(d0, d1, d2, d3, d4, array_name);
        for (int i0 = 0; i0 < dim[0]; i0++) {
            if (vec.at(i0).size() != d1)
                throw std::invalid_argument("Vector size is not constant at dimension 1");

            for (int i1 = 0; i1 < dim[1]; i1++) {
                if (vec.at(i0).at(i1).size() != d2)
                    throw std::invalid_argument("Vector size is not constant at dimension 2");

                for (int i2 = 0; i2 < dim[2]; i2++) {
                    if (vec.at(i0).at(i1).at(i2).size() != d3)
                        throw std::invalid_argument("Vector size is not constant at dimension 3");

                    for (int i3 = 0; i3 < dim[3]; i3++) {
                        if (vec.at(i0).at(i1).at(i2).at(i3).size() != d4)
                            throw std::invalid_argument("Vector size is not constant at dimension 4");

                        for (int i4 = 0; i4 < dim[4]; i4++) {
                            operator()(i0, i1, i2, i3, i4) = vec.at(i0).at(i1).at(i2).at(i3).at(i4);

                        }
                    }
                }
            }
        }

    }

    /**
    * Parametrized constructor from  vector<vector<...<T>> container
    * @param vec container
    * @param array_name array name
    */
    Array5D(const vector<vector<vector<vector<vector<T>>>>> &vec, const string &array_name = "Array5D") {
        this->set_vector(vec);
        this->array_name = array_name;
    }

    /**
    * operator= to vector<vector<...<T>> container
    * @param vec container
    */
    Array5D &operator=(const vector<vector<vector<vector<vector<T>>>>> &vec) {
        this->set_vector(vec);
        return *this;
    }


    /**
    * operator= to flatten vector<T> container
    * @param vec container
    */
    Array5D &operator=(const vector<T> &vec) {
        this->set_flatten_vector(vec);
        return *this;
    }


    vector<size_t> get_shape() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = dim[d];
        return sh;
    }

    vector<size_t> get_strides() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = s[d];
        return sh;
    }

    vector<size_t> get_memory_strides() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = s[d] * sizeof(T);
        return sh;
    }
};


/**
 * Multidimensional (6 - dimensional) array of type T with contiguous memory layout.
 * If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
 * be performed before accessing memory. By default this is turned off.
 * @tparam T data type
 */
template<class T>
class Array6D : public ContiguousArrayND<T> {
    using ContiguousArrayND<T>::array_name;
    using ContiguousArrayND<T>::data;
    using ContiguousArrayND<T>::size;
    using ContiguousArrayND<T>::is_proxy_;

    size_t dim[6] = {0}; ///< dimensions
    size_t s[6] = {0}; ///< strides
    int ndim = 6; ///< number of dimensions
public:

    /**
     * Default empty constructor
     */
    Array6D() = default;

    /**
     * Parametrized constructor with array name
     * @param array_name name of array (for error logging)
     */
    Array6D(const string &array_name) { this->array_name = array_name; }

    /**
    * Parametrized constructor
    * @param d0,... array sizes for different dimensions
    * @param array_name string name of the array
    */
    Array6D(size_t d0, size_t d1, size_t d2, size_t d3, size_t d4, size_t d5, const string &array_name = "Array6D",
            T *new_data = nullptr) {
        init(d0, d1, d2, d3, d4, d5, array_name, new_data);
    }

    /**
    * Setup the dimensions and strides of array
    * @param d0,... array sizes for different dimensions
    */
    void set_dimensions(size_t d0, size_t d1, size_t d2, size_t d3, size_t d4, size_t d5) {

        dim[0] = d0;
        dim[1] = d1;
        dim[2] = d2;
        dim[3] = d3;
        dim[4] = d4;
        dim[5] = d5;

        s[5] = 1;
        s[4] = s[5] * dim[5];
        s[3] = s[4] * dim[4];
        s[2] = s[3] * dim[3];
        s[1] = s[2] * dim[2];
        s[0] = s[1] * dim[1];

        size = s[0] * dim[0];
    };

    /**
    * Initialize array storage, dimensions and strides
    * @param d0,... array sizes for different dimensions
    * @param array_name string name of the array
    */
    void init(size_t d0, size_t d1, size_t d2, size_t d3, size_t d4, size_t d5, const string &array_name = "Array6D",
              T *new_data = nullptr) {
        this->array_name = array_name;

        size_t old_size = size;
        set_dimensions(d0, d1, d2, d3, d4, d5);

        bool new_is_proxy = (new_data != nullptr);

        if (new_is_proxy) {
            if (!is_proxy_) delete[] data;
            data = new_data;
        } else {
            if (size != old_size) {
                T *old_data = data; //preserve the pointer to the old data
                data = new T[size]; // allocate new data
                if (old_data != nullptr) { //
                    size_t min_size = old_size < size ? old_size : size;
                    memcpy(data, old_data, min_size * sizeof(T));
                    if (!is_proxy_) delete[] old_data;
                }
                //memset(data, 0, size * sizeof(T));
            } else {
                //memset(data, 0, size * sizeof(T));
            }
        }

        is_proxy_ = new_is_proxy;
    }

    /**
    * Resize array
    * @param d0,... array sizes for different dimensions
    */
    void resize(size_t d0, size_t d1, size_t d2, size_t d3, size_t d4, size_t d5) {
        init(d0, d1, d2, d3, d4, d5, this->array_name);
    }

    /**
    * Reshape array without reset the data
    * @param d0,... array sizes for different dimensions
    */
    void reshape(size_t d0, size_t d1, size_t d2, size_t d3, size_t d4, size_t d5) {
        //check data size consistency
        size_t new_size = d0 * d1 * d2 * d3 * d4 * d5;
        if (new_size != size)
            throw invalid_argument("Couldn't reshape array when the size is not conserved");
        set_dimensions(d0, d1, d2, d3, d4, d5);
    }

    /**
    * Get array size in dimension "d"
    * @param d dimension index
    */
    size_t get_dim(int d) const {
        return dim[d];
    }

#ifdef MULTIARRAY_INDICES_CHECK

    /**
    * Check indices validity. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    void check_indices(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4, size_t i5) const {

        if ((i0 < 0) | (i0 >= dim[0])) {
            char buf[1024];
            sprintf(buf, "%s: index i0=%ld out of range (0, %ld)\n", array_name.c_str(), i0, dim[0] - 1);
            throw std::out_of_range(buf);
        }

        if ((i1 < 0) | (i1 >= dim[1])) {
            char buf[1024];
            sprintf(buf, "%s: index i1=%ld out of range (0, %ld)\n", array_name.c_str(), i1, dim[1] - 1);
            throw std::out_of_range(buf);
        }

        if ((i2 < 0) | (i2 >= dim[2])) {
            char buf[1024];
            sprintf(buf, "%s: index i2=%ld out of range (0, %ld)\n", array_name.c_str(), i2, dim[2] - 1);
            throw std::out_of_range(buf);
        }

        if ((i3 < 0) | (i3 >= dim[3])) {
            char buf[1024];
            sprintf(buf, "%s: index i3=%ld out of range (0, %ld)\n", array_name.c_str(), i3, dim[3] - 1);
            throw std::out_of_range(buf);
        }

        if ((i4 < 0) | (i4 >= dim[4])) {
            char buf[1024];
            sprintf(buf, "%s: index i4=%ld out of range (0, %ld)\n", array_name.c_str(), i4, dim[4] - 1);
            throw std::out_of_range(buf);
        }

        if ((i5 < 0) | (i5 >= dim[5])) {
            char buf[1024];
            sprintf(buf, "%s: index i5=%ld out of range (0, %ld)\n", array_name.c_str(), i5, dim[5] - 1);
            throw std::out_of_range(buf);
        }

    }

#endif

    /**
    * Array access operator() for reading. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    inline const T &operator()(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4, size_t i5) const {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, i1, i2, i3, i4, i5);
#endif

        return data[i0 * s[0] + i1 * s[1] + i2 * s[2] + i3 * s[3] + i4 * s[4] + i5];

    }

    /**
    * Array access operator() for writing. If preprocessor macro MULTIARRAY_INDICES_CHECK is defined, then the check of index will
    * be performed before accessing memory. By default this is turned off.
    * @param i0,... indices
    */
    inline T &operator()(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4, size_t i5) {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, i1, i2, i3, i4, i5);
#endif

        return data[i0 * s[0] + i1 * s[1] + i2 * s[2] + i3 * s[3] + i4 * s[4] + i5];

    }

    /**
    * Array comparison operator
    * @param other
    */
    bool operator==(const Array6D &other) const {
        //compare dimensions
        for (int d = 0; d < ndim; d++) {
            if (this->dim[d] != other.dim[d])
                return false;
        }
        return ContiguousArrayND<T>::operator==(other);
    }

    /**
    * Convert to nested vector<vector<...<T>> container
    * @return vector container
    */
    vector<vector<vector<vector<vector<vector<T>>>>>> to_vector() const {
        vector<vector<vector<vector<vector<vector<T>>>>>> res;

        res.resize(dim[0]);

        for (int i0 = 0; i0 < dim[0]; i0++) {
            res[i0].resize(dim[1]);


            for (int i1 = 0; i1 < dim[1]; i1++) {
                res[i0][i1].resize(dim[2]);


                for (int i2 = 0; i2 < dim[2]; i2++) {
                    res[i0][i1][i2].resize(dim[3]);


                    for (int i3 = 0; i3 < dim[3]; i3++) {
                        res[i0][i1][i2][i3].resize(dim[4]);


                        for (int i4 = 0; i4 < dim[4]; i4++) {
                            res[i0][i1][i2][i3][i4].resize(dim[5]);


                            for (int i5 = 0; i5 < dim[5]; i5++) {
                                res[i0][i1][i2][i3][i4][i5] = operator()(i0, i1, i2, i3, i4, i5);

                            }
                        }
                    }
                }
            }
        }


        return res;
    } // end to_vector()


    /**
    * Set values to vector<vector<...<T>> container
    * @param vec container
    */
    void set_vector(const vector<vector<vector<vector<vector<vector<T>>>>>> &vec) {
        size_t d0 = 0;
        size_t d1 = 0;
        size_t d2 = 0;
        size_t d3 = 0;
        size_t d4 = 0;
        size_t d5 = 0;
        d0 = vec.size();

        if (d0 > 0) {
            d1 = vec.at(0).size();
            if (d1 > 0) {
                d2 = vec.at(0).at(0).size();
                if (d2 > 0) {
                    d3 = vec.at(0).at(0).at(0).size();
                    if (d3 > 0) {
                        d4 = vec.at(0).at(0).at(0).at(0).size();
                        if (d4 > 0) {
                            d5 = vec.at(0).at(0).at(0).at(0).at(0).size();


                        }
                    }
                }
            }
        }

        init(d0, d1, d2, d3, d4, d5, array_name);
        for (int i0 = 0; i0 < dim[0]; i0++) {
            if (vec.at(i0).size() != d1)
                throw std::invalid_argument("Vector size is not constant at dimension 1");

            for (int i1 = 0; i1 < dim[1]; i1++) {
                if (vec.at(i0).at(i1).size() != d2)
                    throw std::invalid_argument("Vector size is not constant at dimension 2");

                for (int i2 = 0; i2 < dim[2]; i2++) {
                    if (vec.at(i0).at(i1).at(i2).size() != d3)
                        throw std::invalid_argument("Vector size is not constant at dimension 3");

                    for (int i3 = 0; i3 < dim[3]; i3++) {
                        if (vec.at(i0).at(i1).at(i2).at(i3).size() != d4)
                            throw std::invalid_argument("Vector size is not constant at dimension 4");

                        for (int i4 = 0; i4 < dim[4]; i4++) {
                            if (vec.at(i0).at(i1).at(i2).at(i3).at(i4).size() != d5)
                                throw std::invalid_argument("Vector size is not constant at dimension 5");

                            for (int i5 = 0; i5 < dim[5]; i5++) {
                                operator()(i0, i1, i2, i3, i4, i5) = vec.at(i0).at(i1).at(i2).at(i3).at(i4).at(i5);

                            }
                        }
                    }
                }
            }
        }

    }

    /**
    * Parametrized constructor from  vector<vector<...<T>> container
    * @param vec container
    * @param array_name array name
    */
    Array6D(const vector<vector<vector<vector<vector<vector<T>>>>>> &vec, const string &array_name = "Array6D") {
        this->set_vector(vec);
        this->array_name = array_name;
    }

    /**
    * operator= to vector<vector<...<T>> container
    * @param vec container
    */
    Array6D &operator=(const vector<vector<vector<vector<vector<vector<T>>>>>> &vec) {
        this->set_vector(vec);
        return *this;
    }


    /**
    * operator= to flatten vector<T> container
    * @param vec container
    */
    Array6D &operator=(const vector<T> &vec) {
        this->set_flatten_vector(vec);
        return *this;
    }


    vector<size_t> get_shape() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = dim[d];
        return sh;
    }

    vector<size_t> get_strides() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = s[d];
        return sh;
    }

    vector<size_t> get_memory_strides() const {
        vector<size_t> sh(ndim);
        for (int d = 0; d < ndim; d++)
            sh[d] = s[d] * sizeof(T);
        return sh;
    }
};


#endif //ACE_MULTIARRAY_H
