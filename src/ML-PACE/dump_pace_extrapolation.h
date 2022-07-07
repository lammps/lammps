/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(pace/extrapolation,DumpPACEExtrapolation);
// clang-format on
#else

#ifndef LMP_DUMP_PACE_AL_H
#define LMP_DUMP_PACE_AL_H

#include "pair.h"
#include "dump_custom.h"
#include "pair_pace_extrapolation.h"

namespace LAMMPS_NS {
    // forward declaration
    class PairPACEExtrapolation;

    class DumpPACEExtrapolation : public DumpCustom {
        const std::string SPECIES_TYPE_FNAME = "species_types.dat";
        PairPACEExtrapolation *pairPaceExtrapolation;
    public:
        DumpPACEExtrapolation(class LAMMPS *lmp, int nargs, char **argv);

        void write() override;
    };

}    // namespace LAMMPS_NS

#endif
#endif
