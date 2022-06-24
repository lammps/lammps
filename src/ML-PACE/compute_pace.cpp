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

#include "compute_pace.h"
#include "pair_pace_al.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;

ComputePaceAtom::ComputePaceAtom(class LAMMPS *lmp, int narg, char **arg) :
        Compute(lmp, narg, arg) {
    if (narg < 3) error->all(FLERR, "Illegal compute pace/atom command");
    peratom_flag = 1;
}

ComputePaceAtom::~ComputePaceAtom() {
}

void ComputePaceAtom::init() {
    if (force->pair == nullptr)
        error->all(FLERR, "Compute pace/atom requires a pair style be defined");

    int count = 0;
    for (int i = 0; i < modify->ncompute; i++)
        if (strcmp(modify->compute[i]->style, "pace/atom") == 0) count++;
    if (count > 1 && comm->me == 0) error->warning(FLERR, "More than one compute pace/atom");

    pair_pace_al = force->pair_match("pace/al", 1);
    if (!pair_pace_al)
        error->all(FLERR, "Compute pace/atom requires a `pace/al` pair style");

}

void ComputePaceAtom::compute_peratom() {
    invoked_peratom = update->ntimestep;
    auto pair = (PairPACEActiveLearning *) pair_pace_al;
    if (invoked_peratom != pair->bevaluator_timestep)
        error->all(FLERR, "PACE/gamma was not computed on needed timestep");
    vector_atom = pair->extrapolation_grade_gamma;
    scalar = pair->max_gamma_grade_per_structure;
}