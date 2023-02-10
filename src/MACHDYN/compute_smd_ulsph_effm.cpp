// clang-format off
/* ----------------------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the MACHDYN package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 https://www.lammps.org/, Sandia National Laboratories
 LAMMPS development team: developers@lammps.org

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#include "compute_smd_ulsph_effm.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSMD_Ulsph_Effm::ComputeSMD_Ulsph_Effm(LAMMPS *lmp, int narg, char **arg) :
                Compute(lmp, narg, arg) {
        if (narg != 3)
                error->all(FLERR, "Illegal compute smd/ulsph_effm command");
        if (atom->contact_radius_flag != 1)
                error->all(FLERR,
                                "compute smd/ulsph_effm command requires atom_style with contact_radius (e.g. smd)");

        peratom_flag = 1;
        size_peratom_cols = 0;

        nmax = 0;
        dt_vector = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSMD_Ulsph_Effm::~ComputeSMD_Ulsph_Effm() {
        memory->sfree(dt_vector);
}

/* ---------------------------------------------------------------------- */

void ComputeSMD_Ulsph_Effm::init() {

        int count = 0;
        for (int i = 0; i < modify->ncompute; i++)
                if (strcmp(modify->compute[i]->style, "smd/ulsph_effm") == 0)
                        count++;
        if (count > 1 && comm->me == 0)
                error->warning(FLERR, "More than one compute smd/ulsph_effm");
}

/* ---------------------------------------------------------------------- */

void ComputeSMD_Ulsph_Effm::compute_peratom() {
        invoked_peratom = update->ntimestep;

        // grow rhoVector array if necessary

        if (atom->nmax > nmax) {
                memory->sfree(dt_vector);
                nmax = atom->nmax;
                dt_vector = (double *) memory->smalloc(nmax * sizeof(double),
                                "atom:tlsph_dt_vector");
                vector_atom = dt_vector;
        }

        int itmp = 0;
        auto particle_dt = (double *) force->pair->extract("smd/ulsph/effective_modulus_ptr", itmp);
        if (particle_dt == nullptr)
          error->all(FLERR, "compute smd/ulsph_effm failed to access particle_dt array");

        int *mask = atom->mask;
        int nlocal = atom->nlocal;

        for (int i = 0; i < nlocal; i++) {
                if (mask[i] & groupbit) {
                        dt_vector[i] = particle_dt[i];
                } else {
                        dt_vector[i] = 0.0;
                }
        }
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based array
 ------------------------------------------------------------------------- */

double ComputeSMD_Ulsph_Effm::memory_usage() {
        double bytes = (double)nmax * sizeof(double);
        return bytes;
}
