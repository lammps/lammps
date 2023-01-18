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

#include <cstring>
#include <Eigen/Eigen>
#include "compute_smd_ulsph_strain_rate.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "pair.h"

using namespace Eigen;
using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

ComputeSMDULSPHStrainRate::ComputeSMDULSPHStrainRate(LAMMPS *lmp, int narg, char **arg) :
                Compute(lmp, narg, arg) {
        if (narg != 3)
                error->all(FLERR, "Illegal compute smd/ulsph_strain_rate command");

        peratom_flag = 1;
        size_peratom_cols = 6;

        nmax = 0;
        strain_rate_array = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSMDULSPHStrainRate::~ComputeSMDULSPHStrainRate() {
        memory->sfree(strain_rate_array);
}

/* ---------------------------------------------------------------------- */

void ComputeSMDULSPHStrainRate::init() {

        int count = 0;
        for (int i = 0; i < modify->ncompute; i++)
                if (strcmp(modify->compute[i]->style, "smd/ulsph_strain_rate") == 0)
                        count++;
        if (count > 1 && comm->me == 0)
                error->warning(FLERR, "More than one compute smd/ulsph_strain_rate");
}

/* ---------------------------------------------------------------------- */

void ComputeSMDULSPHStrainRate::compute_peratom() {
        invoked_peratom = update->ntimestep;
        int *mask = atom->mask;

        // grow vector array if necessary

        if (atom->nmax > nmax) {
                memory->destroy(strain_rate_array);
                nmax = atom->nmax;
                memory->create(strain_rate_array, nmax, size_peratom_cols, "stresstensorVector");
                array_atom = strain_rate_array;
        }

        int itmp = 0;
        auto L = (Matrix3d *) force->pair->extract("smd/ulsph/velocityGradient_ptr", itmp);
        if (L == nullptr) {
                error->all(FLERR,
                                "compute smd/ulsph_strain_rate could not access any velocity gradients. Are the matching pair styles present?");
        }
        int nlocal = atom->nlocal;
        Matrix3d D;

        for (int i = 0; i < nlocal; i++) {
                if (mask[i] & groupbit) {
                        D = 0.5 * (L[i] + L[i].transpose());
                        strain_rate_array[i][0] = D(0, 0); // xx
                        strain_rate_array[i][1] = D(1, 1); // yy
                        strain_rate_array[i][2] = D(2, 2); // zz
                        strain_rate_array[i][3] = D(0, 1); // xy
                        strain_rate_array[i][4] = D(0, 2); // xz
                        strain_rate_array[i][5] = D(1, 2); // yz
                } else {
                        strain_rate_array[i][0] = 0.0;
                        strain_rate_array[i][1] = 0.0;
                        strain_rate_array[i][2] = 0.0;
                        strain_rate_array[i][3] = 0.0;
                        strain_rate_array[i][4] = 0.0;
                        strain_rate_array[i][5] = 0.0;
                }
        }
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based array
 ------------------------------------------------------------------------- */

double ComputeSMDULSPHStrainRate::memory_usage() {
        double bytes = (double)size_peratom_cols * nmax * sizeof(double);
        return bytes;
}
