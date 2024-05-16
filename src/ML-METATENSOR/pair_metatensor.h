/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#ifdef PAIR_CLASS
// clang-format off
PairStyle(metatensor, PairMetatensor);
// clang-format on
#else

#ifndef LMP_PAIR_METATENSOR_H
#define LMP_PAIR_METATENSOR_H

#include "pair.h"

namespace LAMMPS_NS {
class MetatensorSystemAdaptor;
struct PairMetatensorData;

class PairMetatensor : public Pair {
public:
    PairMetatensor(class LAMMPS *);
    ~PairMetatensor();

    void compute(int, int) override;
    void settings(int, char **) override;
    void coeff(int, char **) override;
    void init_style() override;
    double init_one(int, int) override;
    void init_list(int id, NeighList *ptr) override;

    void allocate();
private:
    PairMetatensorData* mts_data;

    // mapping from LAMMPS types to metatensor types
    int32_t* type_mapping;
};

}    // namespace LAMMPS_NS

#endif
#endif
