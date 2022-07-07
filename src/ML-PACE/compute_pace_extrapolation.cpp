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

//
// Created by Yury Lysogorskiy on 23.06.22.
//

#include "compute_pace_extrapolation.h"
#include "pair_pace_extrapolation.h"


#include "comm.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;

ComputePACEExtrapolation::ComputePACEExtrapolation(class LAMMPS *lmp, int narg, char **arg) :
        Compute(lmp, narg, arg) {
    if (narg < 3) error->all(FLERR, "Illegal compute pace/extrapolation command");
    peratom_flag = 1;
    size_peratom_cols = 0;
    scalar_flag = 1; // compute max of gamma
}

ComputePACEExtrapolation::~ComputePACEExtrapolation() {
}

void ComputePACEExtrapolation::init() {
    if (force->pair == nullptr)
        error->all(FLERR, "Compute pace/extrapolation requires a pair style pace/extrapolation be defined");

    if ((modify->get_compute_by_style("pace/extrapolation").size() > 1) && (comm->me == 0))
        error->warning(FLERR, "More than one instance of compute pace/atom");

    pair_pace_extrapolation = (PairPACEExtrapolation *) force->pair_match("pace/extrapolation", 1);
    if (!pair_pace_extrapolation)
        error->all(FLERR, "Compute pace/extrapolation requires a `pace/extrapolation` pair style");
}

void ComputePACEExtrapolation::invoke_compute_extrapolation_grades() {
    bigint current_timestep = update->ntimestep;
    pair_pace_extrapolation->bevaluator_timestep_shift = current_timestep;
    int old_vflag_fdotr = pair_pace_extrapolation->vflag_fdotr;
    pair_pace_extrapolation->vflag_fdotr = 0;
    pair_pace_extrapolation->is_set_energies_forces = false;

    pair_pace_extrapolation->compute(0, 0);

    pair_pace_extrapolation->is_set_energies_forces = true;
    pair_pace_extrapolation->vflag_fdotr = old_vflag_fdotr;
}

double ComputePACEExtrapolation::compute_scalar() {
    invoked_scalar = update->ntimestep;

    // check the coherence of bevaluator_timestep (when extrapolation grades are computed) and actual timestep
    // if not coherent, change pair->bevaluator_timestep_shift to current timestep
    // and call invoke_compute_extrapolation_grades without updating energies and forces
    if (invoked_scalar != pair_pace_extrapolation->bevaluator_timestep) {
        //utils::logmesg(lmp,"[ComputePaceAtom::compute_scalar] Reseting timestep shift to {} (pace timestep={}) and recomputing\n",invoked_scalar,pair->bevaluator_timestep);
        invoke_compute_extrapolation_grades();
    }

    scalar = pair_pace_extrapolation->max_gamma_grade_per_structure;
    return scalar;
}

void ComputePACEExtrapolation::compute_peratom() {
    invoked_peratom = update->ntimestep;

    // check the coherence of bevaluator_timestep (when extrapolation grades are computed) and actual timestep
    // if not coherent, change pair->bevaluator_timestep_shift to current timestep
    // and call invoke_compute_extrapolation_grades without updating energies and forces
    if (invoked_peratom != pair_pace_extrapolation->bevaluator_timestep) {
        //utils::logmesg(lmp,"[ComputePaceAtom::compute_peratom] Reseting timestep shift to {} (pace timestep={}) and recomputing\n",invoked_peratom,pair->bevaluator_timestep);
        invoke_compute_extrapolation_grades();
    }

    vector_atom = pair_pace_extrapolation->extrapolation_grade_gamma;
}