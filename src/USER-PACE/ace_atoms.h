//
// Created by lysogy36 on 31.01.20.
//

#ifndef ACE_ACE_ATOMS_H
#define ACE_ACE_ATOMS_H

#include <cmath>

#include "ace_types.h"
#include "multiarray/ace_arraynd.h"

#define sqr(x) (x)*(x)

struct ACEAtomicEnvironment {
    int n_atoms_tot;
    DOUBLE_TYPE **x = nullptr;
    SPECIES_TYPE *species_type = nullptr;
    int *num_neighbours = nullptr;
    int **neighbour_list = nullptr;

    ACEAtomicEnvironment() = default;

    ACEAtomicEnvironment(int n_atoms) {
        n_atoms_tot = n_atoms;
        x = new DOUBLE_TYPE *[n_atoms_tot];
        for (int i = 0; i < n_atoms_tot; i++) {
            x[i] = new DOUBLE_TYPE[3];
        }

        species_type = new SPECIES_TYPE[n_atoms_tot];
        for (int i = 0; i < n_atoms_tot; i++) {
            species_type[i] = 0;
        }
    }

    void compute_neighbour_list(DOUBLE_TYPE cutoff = 100.) {
        if (num_neighbours != nullptr) delete[] num_neighbours;

        if (neighbour_list != nullptr) {
            for (int i = 0; i < n_atoms_tot; i++) {
                delete[] neighbour_list[i];
            }
            delete[] neighbour_list;
        }

        num_neighbours = new int[n_atoms_tot];
        neighbour_list = new int *[n_atoms_tot];

        for (int i = 0; i < n_atoms_tot; i++) {

            int num_neigh = 0;
            for (int j = 0; j < n_atoms_tot; j++) {
                if (i != j) {
                    if (sqrt(sqr(x[i][0] - x[j][0]) + sqr(x[i][1] - x[j][1]) + sqr(x[i][2] - x[j][2])) <= cutoff)
                        num_neigh++;
                }
            }

            num_neighbours[i] = num_neigh;
            neighbour_list[i] = new int[num_neigh];
            num_neigh = 0;
            for (int j = 0; j < n_atoms_tot; j++) {
                if (i != j) {
                    if (sqrt(sqr(x[i][0] - x[j][0]) + sqr(x[i][1] - x[j][1]) + sqr(x[i][2] - x[j][2])) <= cutoff) {
                        neighbour_list[i][num_neigh] = j;
                        num_neigh++;
                    }
                }
            }

        }

    }

    void _clean() {

        for (int i = 0; i < n_atoms_tot; i++) {
            delete[] x[i];
        }

        delete[] x;

        for (int i = 0; i < n_atoms_tot; i++) {
            delete[] neighbour_list[i];

        }
        delete[] neighbour_list;
        neighbour_list = nullptr;

        delete[] species_type;
        species_type = nullptr;

        delete[] num_neighbours;
        num_neighbours = nullptr;

    }

    void _copy_from(const ACEAtomicEnvironment &other) {
        n_atoms_tot = other.n_atoms_tot;

        x = new DOUBLE_TYPE *[n_atoms_tot];
        for (int i = 0; i < n_atoms_tot; i++) {
            x[i] = new DOUBLE_TYPE[3];
        }

        neighbour_list = new int *[n_atoms_tot];

        num_neighbours = new int[n_atoms_tot];
        species_type = new SPECIES_TYPE[n_atoms_tot];

        for (int i = 0; i < n_atoms_tot; i++) {
            x[i][0] = other.x[i][0];
            x[i][1] = other.x[i][1];
            x[i][2] = other.x[i][2];

            num_neighbours[i] = other.num_neighbours[i];
            species_type[i] = other.species_type[i];

            neighbour_list[i] = new int[num_neighbours[i]];
            for (int j = 0; j < num_neighbours[i]; j++) {
                neighbour_list[i][j] = other.neighbour_list[i][j];
            }
        }
    }

    ACEAtomicEnvironment(const ACEAtomicEnvironment &other) {
        _copy_from(other);
    }

    ACEAtomicEnvironment &operator=(const ACEAtomicEnvironment &other) {
        _clean();
        _copy_from(other);
        return *this;
    }

    ~ACEAtomicEnvironment() {
        _clean();
    }

};

ACEAtomicEnvironment create_linear_chain(int n, int axis = 2);

ACEAtomicEnvironment create_cube(const DOUBLE_TYPE dr, const DOUBLE_TYPE cube_side_length);

ACEAtomicEnvironment create_bcc(const DOUBLE_TYPE lat);

ACEAtomicEnvironment
create_supercell(ACEAtomicEnvironment &simple_cell, DOUBLE_TYPE lx, DOUBLE_TYPE ly, DOUBLE_TYPE lz, int nx, int ny,
                 int nz);

typedef Array2D<DOUBLE_TYPE> Matrix;

Matrix rotation_matrix(DOUBLE_TYPE theta, DOUBLE_TYPE theta1, DOUBLE_TYPE theta2);

void rotate_structure(ACEAtomicEnvironment &env, Matrix &rotation_matrix);

#endif //ACE_ACE_ATOMS_H
