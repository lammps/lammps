//
// Created by Yury Lysogorskiy on 11.01.20.
//
#ifndef ACE_ARRAY2DLM_H
#define ACE_ARRAY2DLM_H

#include <string>

#include "ace_types.h"
#include "ace_contigous_array.h"
#include "ace_arraynd.h"
#include <stdexcept>
using namespace std;

template<typename T>
class Array2DLM : public ContiguousArrayND<T> {
    using ContiguousArrayND<T>::array_name;
    using ContiguousArrayND<T>::data;
    using ContiguousArrayND<T>::size;
    // dimensions
    LS_TYPE lmax = 0;

    bool is_proxy = false;

public:
    // default empty constructor
    Array2DLM() = default;

    // parametrized constructor
    explicit Array2DLM(LS_TYPE lmax, string array_name = "Array2DLM") {
        init(lmax, array_name);
    }

    //constructor for create slices-proxy
    Array2DLM(LS_TYPE lmax, T *data_ptr, string array_name = "Array2DLM") {
        this->lmax = lmax;
        this->size = (lmax + 1) * (lmax + 1);
        this->data = data_ptr;
        this->array_name = array_name;
        is_proxy = true;
    };

    //destructor
    ~Array2DLM() {
        if (!is_proxy) {
            if (data != nullptr) delete[] data;
        }
        data = nullptr;
    }

    //initialize array
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

    inline const T &operator()(LS_TYPE l, MS_TYPE m) const {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(l, m);
#endif
        //l^2 + l + m
        return data[l * (l + 1) + m];
    }

    inline T &operator()(LS_TYPE l, MS_TYPE m) {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(l, m);
#endif
        //l^2 + l + m
        return data[l * (l + 1) + m];
    }
};


template<typename T>
class Array3DLM : public ContiguousArrayND<T> {
    using ContiguousArrayND<T>::array_name;
    using ContiguousArrayND<T>::data;
    using ContiguousArrayND<T>::size;
    // dimensions
    LS_TYPE lmax = 0;
    // dimensions
    size_t dim[1] = {0};
    // strides
    size_t s[1] = {0};

    Array1D<Array2DLM<T> *> _proxy_slices;
public:
    // default empty constructor
    Array3DLM() = default;

    Array3DLM(string array_name) {
        this->array_name = array_name;
    };

    // parametrized constructor
    explicit Array3DLM(size_t d0, LS_TYPE lmax, string array_name = "Array3DLM") {
        init(d0, lmax, array_name);
    }

    //initialize array
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

    void _clear_proxies() {
        for (size_t i0 = 0; i0 < _proxy_slices.get_dim(0); ++i0) {
            delete _proxy_slices(i0);
            _proxy_slices(i0) = nullptr;
        }
    }

    ~Array3DLM() {
        _clear_proxies();
    }

    void resize(size_t d0, LS_TYPE lmax) {
        _clear_proxies();
        init(d0, lmax, this->array_name);
    }

    size_t get_dim(int d) const {
        return dim[d];
    }

#ifdef MULTIARRAY_INDICES_CHECK

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

    inline const T &operator()(size_t i0, LS_TYPE l, MS_TYPE m) const {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, l, m);
#endif
        return data[i0 * s[0] + l * (l + 1) + m];
    }

    inline T &operator()(size_t i0, LS_TYPE l, MS_TYPE m) {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, l, m);
#endif
        return data[i0 * s[0] + l * (l + 1) + m];
    }

    inline const Array2DLM<T> &operator()(size_t i0) const {
        //return proxy Array2DLM pointing to i0,i1,l=0,m=0
        return *_proxy_slices(i0);
    }

    inline Array2DLM<T> &operator()(size_t i0) {
        return *_proxy_slices(i0);
    }
};


template<typename T>
class Array4DLM : public ContiguousArrayND<T> {
    using ContiguousArrayND<T>::array_name;
    using ContiguousArrayND<T>::data;
    using ContiguousArrayND<T>::size;
    // dimensions
    LS_TYPE lmax = 0;
    // dimensions
    size_t dim[2] = {0, 0};
    // strides
    size_t s[2] = {0, 0};

    Array2D<Array2DLM<T> *> _proxy_slices;
public:
    // default empty constructor
    Array4DLM() = default;

    Array4DLM(string array_name) {
        this->array_name = array_name;
    };

    // parametrized constructor
    explicit Array4DLM(size_t d0, size_t d1, LS_TYPE lmax, string array_name = "Array4DLM") {
        init(d0, d1, lmax, array_name);
    }

    //initialize array
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

    void _clear_proxies() {

        for (size_t i0 = 0; i0 < _proxy_slices.get_dim(0); ++i0)
            for (size_t i1 = 0; i1 < _proxy_slices.get_dim(1); ++i1) {
                delete _proxy_slices(i0, i1);
                _proxy_slices(i0, i1) = nullptr;
            }
    }

    ~Array4DLM() {
        _clear_proxies();
    }

    void resize(size_t d0, size_t d1, LS_TYPE lmax) {
        _clear_proxies();
        init(d0, d1, lmax, this->array_name);
    }

    size_t get_dim(int d) const {
        return dim[d];
    }

#ifdef MULTIARRAY_INDICES_CHECK

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

    inline const T &operator()(size_t i0, size_t i1, LS_TYPE l, MS_TYPE m) const {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, i1, l, m);
#endif
        return data[i0 * s[0] + i1 * s[1] + l * (l + 1) + m];
    }

    inline T &operator()(size_t i0, size_t i1, LS_TYPE l, MS_TYPE m) {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices(i0, i1, l, m);
#endif
        return data[i0 * s[0] + i1 * s[1] + l * (l + 1) + m];
    }

    inline const Array2DLM<T> &operator()(size_t i0, size_t i1) const {
        //return proxy Array2DLM pointing to i0,i1,l=0,m=0
        return *_proxy_slices(i0, i1);
    }

    inline Array2DLM<T> &operator()(size_t i0, size_t i1) {
        return *_proxy_slices(i0, i1);
    }
};

#endif //ACE_ARRAY2DLM_H