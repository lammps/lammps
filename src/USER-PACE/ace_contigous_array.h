//
// Created by Yury Lysogorskiy on 11.01.20.
//
#ifndef ACE_CONTIGUOUSARRAYND_H
#define ACE_CONTIGUOUSARRAYND_H

#include <string>
#include "ace_types.h"

using namespace std;


template<typename T>
class ContiguousArrayND {
protected:
    T *data = nullptr;
    size_t size = 0;
    string array_name = "Array";
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
        if (this != &other) {
            array_name = other.array_name;
            size = other.size;
            if (size > 0) {
                data = new T[size];
                for (size_t ind = 0; ind < size; ind++)
                    data[ind] = other.data[ind];
            }
        }
        return *this;
    }


    //TODO: make destructor virtual, check the destructors in inherited classes
    //destructor
    ~ContiguousArrayND() {
#ifdef MULTIARRAY_LIFE_CYCLE
        cout<<array_name<<"::~destructor"<<endl;
#endif
        if (data != nullptr) delete[] data;
        data = nullptr;
    }

    void set_array_name(const string &name) {
        this->array_name = name;
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

    bool operator==(const ContiguousArrayND &other) const {
        if (this->size != other.size)
            return false;


        for (size_t i = 0; i < this->size; ++i) {
            if (this->data[i] != other.data[i])
                return false;
        }

        return true;
    }
};


#endif //ACE_CONTIGUOUSARRAYND_H