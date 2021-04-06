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
#ifndef ACE_CONTIGUOUSARRAYND_H
#define ACE_CONTIGUOUSARRAYND_H

#include <string>

#include "ace_types.h"

using namespace std;

/**
 * Common predecessor class to represent multidimensional array of type T
 * and store it in memory contiguous form
 *
 * @tparam T data type
 */
template<typename T>
class ContiguousArrayND {
protected:
    T *data = nullptr; ///< pointer to contiguous data
    size_t size = 0; ///< total array size
    string array_name = "Array"; ///<array name
    bool is_proxy_ = false; ///< array is proxy (wrapper) and not owner of the memory
public:

    /**
     * Default empty constructor
     */
    ContiguousArrayND() = default;


    /**
     *  Constructor with array name
     * @param array_name name of array (for error logging)
     */
    ContiguousArrayND(string array_name) : array_name(array_name) {};

    /**
     * Copy constructor
     * @param other other ContiguousArrayND
     */
    ContiguousArrayND(const ContiguousArrayND &other) : array_name(other.array_name), size(other.size), is_proxy_(other.is_proxy_) {
#ifdef MULTIARRAY_LIFE_CYCLE
        cout<<array_name<<"::copy constructor"<<endl;
#endif
        if(!is_proxy_) { //if not the proxy, then copy the values
            if (size > 0) {
                data = new T[size];
                for (size_t ind = 0; ind < size; ind++)
                    data[ind] = other.data[ind];
            }
        } else { //is proxy, then copy the pointer
            data = other.data;
        }
    }

    /**
     * Overload operator=
     * @param other another  ContiguousArrayND
     * @return itself
     */

    ContiguousArrayND &operator=(const ContiguousArrayND &other) {
#ifdef MULTIARRAY_LIFE_CYCLE
        cout<<array_name<<"::operator="<<endl;
#endif
        if (this != &other) {
            array_name = other.array_name;
            size = other.size;
            is_proxy_ = other.is_proxy_;
            if(!is_proxy_) { //if not the proxy, then copy the values
                if (size > 0) {

                    if(data!=nullptr) delete[] data;
                    data = new T[size];

                    for (size_t ind = 0; ind < size; ind++)
                        data[ind] = other.data[ind];
                }
            } else { //is proxy, then copy the pointer
                data = other.data;
            }
        }
        return *this;
    }


    //TODO: make destructor virtual, check the destructors in inherited classes

    /**
     * Destructor
     */
    ~ContiguousArrayND() {
#ifdef MULTIARRAY_LIFE_CYCLE
        cout<<array_name<<"::~destructor"<<endl;
#endif
        if(! is_proxy_) {
            delete[] data;
        }
        data = nullptr;
    }

    /**
     * Set array name
     * @param name array name
     */
    void set_array_name(const string &name) {
        this->array_name = name;
    }

    /**
     * Get total number of elements in array (its size)
     * @return array size
     */
    size_t get_size() const {
        return size;
    }

    /**
     * Fill array with value
     * @param value value to fill
     */
    void fill(T value) {
        for (size_t ind = 0; ind < size; ind++)
            data[ind] = value;
    }

    /**
     * Get array data at absolute index ind for reading
     * @param ind absolute index
     * @return array value
     */
    inline const T &get_data(size_t ind) const {
#ifdef MULTIARRAY_INDICES_CHECK
        if ((ind < 0) | (ind >= size)) {
            printf("%s: get_data ind=%d out of range (0, %d)\n", array_name, ind, size);
            exit(EXIT_FAILURE);
        }
#endif
        return data[ind];
    }

    /**
     * Get array data at absolute index ind for writing
     * @param ind absolute index
     * @return array value
     */
    inline T &get_data(size_t ind) {
#ifdef MULTIARRAY_INDICES_CHECK
        if ((ind < 0) | (ind >= size)) {
            printf("%s: get_data ind=%ld out of range (0, %ld)\n", array_name.c_str(), ind, size);
            exit(EXIT_FAILURE);
        }
#endif
        return data[ind];
    }

    /**
     * Get array data pointer
     * @return data array pointer
     */
    inline T* get_data() const {
        return data;
    }

    /**
     * Overload comparison operator==
     * Compare the total size and array values elementwise.
     *
     * @param other another array
     * @return
     */
    bool operator==(const ContiguousArrayND &other) const {
        if (this->size != other.size)
            return false;

        for (size_t i = 0; i < this->size; ++i)
            if (this->data[i] != other.data[i])
                return false;

        return true;
    }


    /**
    * Convert to flatten vector<T> container
    * @return vector container
    */
    vector<T> to_flatten_vector() const {
        vector<T> res;

        res.resize(size);
        size_t vec_ind = 0;

        for (int vec_ind = 0; vec_ind < size; vec_ind++)
                res.at(vec_ind) = data[vec_ind];

        return res;
    } // end to_flatten_vector()


    /**
    * Set values from flatten vector<T> container
    * @param vec container
    */
    void set_flatten_vector(const vector<T> &vec) {
        if (vec.size() != size)
            throw std::invalid_argument("Flatten vector size is not consistent with expected size");
        for (size_t i = 0; i < size; i++) {
            data[i] = vec[i];
        }
    }

    bool is_proxy(){
        return is_proxy_;
    }

};


#endif //ACE_CONTIGUOUSARRAYND_H