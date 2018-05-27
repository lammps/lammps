/* ----------------------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the USER-SMD package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#include <cstring>
#include "compute_smd_tlsph_dt.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "pair.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSMDTlsphDt::ComputeSMDTlsphDt(LAMMPS *lmp, int narg, char **arg) :
                Compute(lmp, narg, arg) {
        if (narg != 3)
                error->all(FLERR, "Illegal compute smd/tlsph_dt command");
        if (atom->contact_radius_flag != 1)
                error->all(FLERR,
                                "compute smd/tlsph_dt command requires atom_style with contact_radius (e.g. smd)");

        peratom_flag = 1;
        size_peratom_cols = 0;

        nmax = 0;
        dt_vector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeSMDTlsphDt::~ComputeSMDTlsphDt() {
        memory->sfree(dt_vector);
}

/* ---------------------------------------------------------------------- */

void ComputeSMDTlsphDt::init() {

        int count = 0;
        for (int i = 0; i < modify->ncompute; i++)
                if (strcmp(modify->compute[i]->style, "smd/tlsph_dt") == 0)
                        count++;
        if (count > 1 && comm->me == 0)
                error->warning(FLERR, "More than one compute smd/tlsph_dt");
}

/* ---------------------------------------------------------------------- */

void ComputeSMDTlsphDt::compute_peratom() {
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
        double *particle_dt = (double *) force->pair->extract("smd/tlsph/particle_dt_ptr",
                        itmp);
        if (particle_dt == NULL) {
                error->all(FLERR,
                                "compute smd/tlsph_dt failed to access particle_dt array");
        }

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

double ComputeSMDTlsphDt::memory_usage() {
        double bytes = nmax * sizeof(double);
        return bytes;
}
