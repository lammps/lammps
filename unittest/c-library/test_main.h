/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef TEST_MAIN_H
#define TEST_MAIN_H

extern bool verbose;

#define EXPECT_FP_LE_WITH_EPS(val1, val2, eps)                \
    do {                                                      \
        const double diff = fabs(val1 - val2);                \
        const double div  = std::min(fabs(val1), fabs(val2)); \
        const double err  = (div == 0.0) ? diff : diff / div; \
        EXPECT_PRED_FORMAT2(::testing::DoubleLE, err, eps);   \
    } while (0);

#endif
