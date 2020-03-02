//
// Created by lysogy36 on 31.01.20.
//
#include "ace_types.h"
#include "ace_atoms.h"
#include "multiarray/ace_arraynd.h"


ACEAtomicEnvironment create_cube(const DOUBLE_TYPE dr, const DOUBLE_TYPE cube_side_length) {
    int n_atoms = 0;
    for (DOUBLE_TYPE x = -cube_side_length / 2; x <= cube_side_length / 2 + dr / 2; x += dr)
        for (DOUBLE_TYPE y = -cube_side_length / 2; y <= cube_side_length / 2 + dr / 2; y += dr)
            for (DOUBLE_TYPE z = -cube_side_length / 2; z <= cube_side_length / 2 + dr / 2; z += dr)
                n_atoms++;

    ACEAtomicEnvironment a(n_atoms);
    int i = 0;
    for (DOUBLE_TYPE x = -cube_side_length / 2; x <= cube_side_length / 2 + dr / 2; x += dr)
        for (DOUBLE_TYPE y = -cube_side_length / 2; y <= cube_side_length / 2 + dr / 2; y += dr)
            for (DOUBLE_TYPE z = -cube_side_length / 2; z <= cube_side_length / 2 + dr / 2; z += dr) {
                a.x[i][0] = x;
                a.x[i][1] = y;
                a.x[i][2] = z;
                i++;
            }

    a.compute_neighbour_list();

    return a;

}

ACEAtomicEnvironment create_linear_chain(const int n, const int axis) {
    ACEAtomicEnvironment a(n);
    for (int i = 0; i < a.n_atoms_tot; i++) {
        //a.x[i] = new DOUBLE_TYPE[3];
        a.x[i][0] = 0;
        a.x[i][1] = 0;
        a.x[i][2] = 0;
        a.x[i][axis] = (1. * i - (double) n / 2.);
    }

    a.compute_neighbour_list();

    return a;
}

ACEAtomicEnvironment
create_supercell(ACEAtomicEnvironment &simple_cell, DOUBLE_TYPE lx, DOUBLE_TYPE ly, DOUBLE_TYPE lz, int nx, int ny,
                 int nz) {
    int number_of_cells = nx * ny * nz;
    ACEAtomicEnvironment a(simple_cell.n_atoms_tot * number_of_cells);

    int at_i = 0;
    for (int ix = 0; ix < nx; ix++)
        for (int iy = 0; iy < ny; iy++)
            for (int iz = 0; iz < nz; iz++) {
                for (int simple_at_i = 0; simple_at_i < simple_cell.n_atoms_tot; simple_at_i++) {
                    a.x[at_i][0] = simple_cell.x[simple_at_i][0] + lx * ix;
                    a.x[at_i][1] = simple_cell.x[simple_at_i][1] + ly * iy;
                    a.x[at_i][2] = simple_cell.x[simple_at_i][2] + lz * iz;
                    at_i++;
                }
            }
    a.compute_neighbour_list();
    return a;
}

ACEAtomicEnvironment create_bcc(const DOUBLE_TYPE lat) {
    ACEAtomicEnvironment a(9);

    a.x[0][0] = -lat / 2.;
    a.x[0][1] = -lat / 2.;
    a.x[0][2] = -lat / 2.;


    a.x[1][0] = lat / 2.;
    a.x[1][1] = -lat / 2.;
    a.x[1][2] = -lat / 2.;

    a.x[2][0] = -lat / 2.;
    a.x[2][1] = lat / 2.;
    a.x[2][2] = -lat / 2.;

    a.x[3][0] = lat / 2.;
    a.x[3][1] = lat / 2.;
    a.x[3][2] = -lat / 2.;


    a.x[4][0] = -lat / 2.;
    a.x[4][1] = -lat / 2.;
    a.x[4][2] = lat / 2.;

    a.x[5][0] = lat / 2.;
    a.x[5][1] = -lat / 2.;
    a.x[5][2] = lat / 2.;

    a.x[6][0] = -lat / 2.;
    a.x[6][1] = lat / 2.;
    a.x[6][2] = lat / 2.;

    a.x[7][0] = lat / 2.;
    a.x[7][1] = lat / 2.;
    a.x[7][2] = lat / 2.;

    a.x[8][0] = 0;
    a.x[8][1] = 0;
    a.x[8][2] = 0;

    a.compute_neighbour_list();

    return a;
}

typedef Array2D<DOUBLE_TYPE> Matrix;

Matrix rotation_matrix(DOUBLE_TYPE theta, DOUBLE_TYPE theta1, DOUBLE_TYPE theta2) {
    DOUBLE_TYPE Rx[3][3] = {{1, 0,          0},
                            {0, cos(theta), -sin(theta)},
                            {0, sin(theta), cos(theta)},
    };

    DOUBLE_TYPE Ry[3][3] = {{cos(theta1),  0, sin(theta1)},
                            {0,            1, 0},
                            {-sin(theta1), 0, cos(theta1)}
    };

    DOUBLE_TYPE Rz[3][3] = {{cos(theta2), -sin(theta2), 0},
                            {sin(theta2), cos(theta2),  0},
                            {0,           0,            1}
    };

    DOUBLE_TYPE R[3][3] = {0};
    int i, j, k;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++)
                R[i][j] += Rx[i][k] * Ry[k][j];
        }

    DOUBLE_TYPE R2[3][3] = {0};
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++)
                R2[i][j] += R[i][k] * Rz[k][j];
        }

    Matrix res(3, 3);
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            res(i, j) = R2[i][j];
        }
    }

    return res;
}


void rotate_structure(ACEAtomicEnvironment &env, Matrix &rotation_matrix) {
    int nat, i, j, k;

    for (nat = 0; nat < env.n_atoms_tot; nat++) {
        DOUBLE_TYPE r[3] = {0};
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                r[i] += rotation_matrix(i, j) * env.x[nat][j];
            }
        }
        for (i = 0; i < 3; i++)
            env.x[nat][i] = r[i];
    }
}