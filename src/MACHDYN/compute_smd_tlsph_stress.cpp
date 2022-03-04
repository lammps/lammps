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
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#include "compute_smd_tlsph_stress.h"
#include <cmath>
#include <cstring>
#include <Eigen/Eigen>
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


/*
 * deviator of a tensor
 */
static Matrix3d Deviator(Matrix3d M) {
        Matrix3d eye;
        eye.setIdentity();
        eye *= M.trace() / 3.0;
        return M - eye;
}

/* ---------------------------------------------------------------------- */

ComputeSMDTLSPHStress::ComputeSMDTLSPHStress(LAMMPS *lmp, int narg, char **arg) :
                Compute(lmp, narg, arg) {
        if (narg != 3)
                error->all(FLERR, "Illegal compute smd/tlsph_stress command");

        peratom_flag = 1;
        size_peratom_cols = 7;

        nmax = 0;
        stress_array = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSMDTLSPHStress::~ComputeSMDTLSPHStress() {
        memory->sfree(stress_array);
}

/* ---------------------------------------------------------------------- */

void ComputeSMDTLSPHStress::init() {

        int count = 0;
        for (int i = 0; i < modify->ncompute; i++)
                if (strcmp(modify->compute[i]->style, "smd/tlsph_stress") == 0)
                        count++;
        if (count > 1 && comm->me == 0)
                error->warning(FLERR, "More than one compute smd/tlsph_stress");
}

/* ---------------------------------------------------------------------- */

void ComputeSMDTLSPHStress::compute_peratom() {
        invoked_peratom = update->ntimestep;
        Matrix3d stress_deviator;
        double von_mises_stress;

        // grow vector array if necessary

        if (atom->nmax > nmax) {
                memory->destroy(stress_array);
                nmax = atom->nmax;
                memory->create(stress_array, nmax, size_peratom_cols, "stresstensorVector");
                array_atom = stress_array;
        }

        int itmp = 0;
        Matrix3d *T = (Matrix3d *) force->pair->extract("smd/tlsph/stressTensor_ptr", itmp);
        if (T == nullptr) {
                error->all(FLERR, "compute smd/tlsph_stress could not access stress tensors. Are the matching pair styles present?");
        }
        int nlocal = atom->nlocal;
        int *mask = atom->mask;

        for (int i = 0; i < nlocal; i++) {
                if (mask[i] & groupbit) {
                        stress_deviator = Deviator(T[i]);
                        von_mises_stress = sqrt(3. / 2.) * stress_deviator.norm();
                        stress_array[i][0] = T[i](0, 0); // xx
                        stress_array[i][1] = T[i](1, 1); // yy
                        stress_array[i][2] = T[i](2, 2); // zz
                        stress_array[i][3] = T[i](0, 1); // xy
                        stress_array[i][4] = T[i](0, 2); // xz
                        stress_array[i][5] = T[i](1, 2); // yz
                        stress_array[i][6] = von_mises_stress;
                } else {
                        for (int j = 0; j < size_peratom_cols; j++) {
                                stress_array[i][j] = 0.0;
                        }
                }
        }
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based array
 ------------------------------------------------------------------------- */

double ComputeSMDTLSPHStress::memory_usage() {
        double bytes = (double)size_peratom_cols * nmax * sizeof(double);
        return bytes;
}
