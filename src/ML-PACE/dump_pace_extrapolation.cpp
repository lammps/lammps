// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "dump_pace_extrapolation.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "modify.h"


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpPACEExtrapolation::DumpPACEExtrapolation(struct LAMMPS *lmp, int nargs, char **argv) : DumpCustom(lmp, nargs,
                                                                                                      argv) {
    pairPaceExtrapolation = (PairPACEExtrapolation *) force->pair_match("pace/extrapolation", 1);
    if (!pairPaceExtrapolation)
        error->all(FLERR, "Dump pace/extrapolation requires a `pace/extrapolation` pair style");

    // Save atomtypes to elements mapping into SPECIES_TYPE_FNAME to be used later in Python for loading extrapolative structures
    FILE *species_type_file = fopen(SPECIES_TYPE_FNAME.c_str(), "w");
    const int n = atom->ntypes;
    for (int i = 0; i < n; i++) {
        auto elemname = pairPaceExtrapolation->element_names[i].c_str();
        fprintf(species_type_file, "%s ", elemname);
    }
    fclose(species_type_file);
}

void DumpPACEExtrapolation::write() {
    int current_time_step = update->ntimestep;
    // dump only if
    // 1) extrapolation grades were computed on current timestep AND
    // 2) max extrapolation grade > gamma_lower_bound
    if (current_time_step == pairPaceExtrapolation->bevaluator_timestep) {
        if (pairPaceExtrapolation->max_gamma_grade_per_structure > pairPaceExtrapolation->gamma_lower_bound) {
            DumpCustom::write();
            MPI_Barrier(world);
        }
    }
}
