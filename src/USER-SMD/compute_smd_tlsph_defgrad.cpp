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

#include "compute_smd_tlsph_defgrad.h"
#include <cstring>
#include <Eigen/Eigen>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace Eigen;
using namespace std;
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSMDTLSPHDefgrad::ComputeSMDTLSPHDefgrad(LAMMPS *lmp, int narg, char **arg) :
                Compute(lmp, narg, arg) {
        if (narg != 3)
                error->all(FLERR, "Illegal compute smd/tlsph_defgrad command");

        peratom_flag = 1;
        size_peratom_cols = 10;

        nmax = 0;
        defgradVector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeSMDTLSPHDefgrad::~ComputeSMDTLSPHDefgrad() {
        memory->sfree(defgradVector);
}

/* ---------------------------------------------------------------------- */

void ComputeSMDTLSPHDefgrad::init() {

        int count = 0;
        for (int i = 0; i < modify->ncompute; i++)
                if (strcmp(modify->compute[i]->style, "smd/tlsph_defgrad") == 0)
                        count++;
        if (count > 1 && comm->me == 0)
                error->warning(FLERR, "More than one compute smd/tlsph_defgrad");
}

/* ---------------------------------------------------------------------- */

void ComputeSMDTLSPHDefgrad::compute_peratom() {
        double **defgrad = atom->smd_data_9;
        Matrix3d F;
        invoked_peratom = update->ntimestep;

        // grow vector array if necessary
        if (atom->nmax > nmax) {
                memory->destroy(defgradVector);
                nmax = atom->nmax;
                memory->create(defgradVector, nmax, size_peratom_cols, "defgradVector");
                array_atom = defgradVector;
        }

        int *mask = atom->mask;
        int nlocal = atom->nlocal;

        for (int i = 0; i < nlocal; i++) {
                if (mask[i] & groupbit) {
                        F(0, 0) = defgrad[i][0];
                        F(0, 1) = defgrad[i][1];
                        F(0, 2) = defgrad[i][2];
                        F(1, 0) = defgrad[i][3];
                        F(1, 1) = defgrad[i][4];
                        F(1, 2) = defgrad[i][5];
                        F(2, 0) = defgrad[i][6];
                        F(2, 1) = defgrad[i][7];
                        F(2, 2) = defgrad[i][8];

                        defgradVector[i][0] = F(0, 0);
                        defgradVector[i][1] = F(0, 1);
                        defgradVector[i][2] = F(0, 2);
                        defgradVector[i][3] = F(1, 0);
                        defgradVector[i][4] = F(1, 1);
                        defgradVector[i][5] = F(1, 2);
                        defgradVector[i][6] = F(2, 0);
                        defgradVector[i][7] = F(2, 1);
                        defgradVector[i][8] = F(2, 2);
                        defgradVector[i][9] = F.determinant();
                } else {
                        defgradVector[i][0] = 1.0;
                        defgradVector[i][1] = 0.0;
                        defgradVector[i][2] = 0.0;
                        defgradVector[i][3] = 0.0;
                        defgradVector[i][4] = 1.0;
                        defgradVector[i][5] = 0.0;
                        defgradVector[i][6] = 0.0;
                        defgradVector[i][7] = 0.0;
                        defgradVector[i][8] = 1.0;
                        defgradVector[i][9] = 1.0;
                }
        }
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based array
 ------------------------------------------------------------------------- */

double ComputeSMDTLSPHDefgrad::memory_usage() {
        double bytes = size_peratom_cols * nmax * sizeof(double);
        return bytes;
}
