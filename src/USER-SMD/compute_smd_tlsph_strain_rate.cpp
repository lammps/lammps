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
#include <Eigen/Eigen>
#include "compute_smd_tlsph_strain_rate.h"
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

ComputeSMDTLSPHStrainRate::ComputeSMDTLSPHStrainRate(LAMMPS *lmp, int narg, char **arg) :
                Compute(lmp, narg, arg) {
        if (narg != 3)
                error->all(FLERR, "Illegal compute smd/ulsph_strain_rate command");

        peratom_flag = 1;
        size_peratom_cols = 6;

        nmax = 0;
        strain_rate_array = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeSMDTLSPHStrainRate::~ComputeSMDTLSPHStrainRate() {
        memory->sfree(strain_rate_array);
}

/* ---------------------------------------------------------------------- */

void ComputeSMDTLSPHStrainRate::init() {

        int count = 0;
        for (int i = 0; i < modify->ncompute; i++)
                if (strcmp(modify->compute[i]->style, "smd/ulsph_strain_rate") == 0)
                        count++;
        if (count > 1 && comm->me == 0)
                error->warning(FLERR, "More than one compute smd/ulsph_strain_rate");
}

/* ---------------------------------------------------------------------- */

void ComputeSMDTLSPHStrainRate::compute_peratom() {
        invoked_peratom = update->ntimestep;

        // grow vector array if necessary

        if (atom->nmax > nmax) {
                memory->destroy(strain_rate_array);
                nmax = atom->nmax;
                memory->create(strain_rate_array, nmax, size_peratom_cols, "stresstensorVector");
                array_atom = strain_rate_array;
        }

        int itmp = 0;
        Matrix3d *D = (Matrix3d *) force->pair->extract("smd/tlsph/strain_rate_ptr", itmp);
        if (D == NULL) {
                error->all(FLERR,
                                "compute smd/tlsph_strain_rate could not access strain rate. Are the matching pair styles present?");
        }

        int nlocal = atom->nlocal;

        for (int i = 0; i < nlocal; i++) {

                strain_rate_array[i][0] = D[i](0, 0); // xx
                strain_rate_array[i][1] = D[i](1, 1); // yy
                strain_rate_array[i][2] = D[i](2, 2); // zz
                strain_rate_array[i][3] = D[i](0, 1); // xy
                strain_rate_array[i][4] = D[i](0, 2); // xz
                strain_rate_array[i][5] = D[i](1, 2); // yz
        }
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based array
 ------------------------------------------------------------------------- */

double ComputeSMDTLSPHStrainRate::memory_usage() {
        double bytes = size_peratom_cols * nmax * sizeof(double);
        return bytes;
}
