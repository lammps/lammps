//
// Created by Yury Lysogorskiy on 11.01.20.
//
#ifndef ACE_CONTIGOUSARRAYND_H
#define ACE_CONTIGOUSARRAYND_H

#include <string>
#include "ace_types.h"

using namespace std;


template<typename T>
class ContiguousArrayND {
protected:
    string array_name = "Array";
    T *data = nullptr;
    size_t size = 0;
public:

    //default empty constructor
    ContiguousArrayND() = default;

    //default empty constructor
    ContiguousArrayND(string array_name) : array_name(array_name) {};

    //copy constructor
    ContiguousArrayND(const ContiguousArrayND &other) : array_name(other.array_name), size(other.size) {
#ifdef MULTIARRAY_LIFE_CYCLE
        cout<<array_name<<"::copy constructor"<<endl;
#endif
        if (size > 0) {
            data = new T[size];
            for (size_t ind = 0; ind < size; ind++)
                data[ind] = other.data[ind];
        }
    }

    //operator=
    ContiguousArrayND &operator=(const ContiguousArrayND &other) {
#ifdef MULTIARRAY_LIFE_CYCLE
        cout<<array_name<<"::operator="<<endl;
#endif
        array_name = other.array_name;
        size = other.size;
        if (size > 0) {
            data = new T[size];
            for (size_t ind = 0; ind < size; ind++)
                data[ind] = other.data[ind];
        }
        return *this;
    }


    //destructor
    ~ContiguousArrayND() {
#ifdef MULTIARRAY_LIFE_CYCLE
        cout<<array_name<<"::~destructor"<<endl;
#endif
        if (data != nullptr) delete[] data;
        data = nullptr;
    }

    size_t get_size() const {
        return size;
    }

    void fill(T value) {
        for (size_t ind = 0; ind < size; ind++)
            data[ind] = value;
    }

    inline const T &get_data(size_t ind) const {
#ifdef MULTIARRAY_INDICES_CHECK
        if ((ind < 0) | (ind >= size)) {
            printf("%s: get_data ind=%d out of range (0, %d)\n", array_name, ind, size);
            exit(EXIT_FAILURE);
        }
#endif
        return data[ind];
    }

    inline T &get_data(size_t ind) {
#ifdef MULTIARRAY_INDICES_CHECK
        if ((ind < 0) | (ind >= size)) {
            printf("%s: get_data ind=%d out of range (0, %d)\n", array_name, ind, size);
            exit(EXIT_FAILURE);
        }
#endif
        return data[ind];
    }

};


#endif //ACE_CONTIGOUSARRAYND_H