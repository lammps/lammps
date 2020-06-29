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

#include "compute_smd_triangle_vertices.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace std;
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSMDTriangleVertices::ComputeSMDTriangleVertices(LAMMPS *lmp, int narg, char **arg) :
        Compute(lmp, narg, arg) {
    if (narg != 3)
        error->all(FLERR, "Illegal compute smd/triangle_vertices command");

    peratom_flag = 1;
    size_peratom_cols = 9;

    nmax = 0;
    outputVector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeSMDTriangleVertices::~ComputeSMDTriangleVertices() {
    memory->sfree(outputVector);
}

/* ---------------------------------------------------------------------- */

void ComputeSMDTriangleVertices::init() {

    int count = 0;
    for (int i = 0; i < modify->ncompute; i++)
        if (strcmp(modify->compute[i]->style, "smd/triangle_vertices") == 0)
            count++;
    if (count > 1 && comm->me == 0)
        error->warning(FLERR, "More than one compute smd/triangle_vertices");
}

/* ---------------------------------------------------------------------- */

void ComputeSMDTriangleVertices::compute_peratom() {

        double **smd_data_9 = atom->smd_data_9;
        tagint *mol = atom->molecule;

    invoked_peratom = update->ntimestep;

    // grow vector array if necessary

    if (atom->nmax > nmax) {
        memory->destroy(outputVector);
        nmax = atom->nmax;
        memory->create(outputVector, nmax, size_peratom_cols, "defgradVector");
        array_atom = outputVector;
    }

    /*
     * triangle vertices are stored using the smd_data_9 array ...
     * this is a hack but ok for now as I do not have to create additional storage space
     * all triangle particles have molecule id >= 65535
     */

    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
        if ((mask[i] & groupbit) && (mol[i] >= 65535) ){
            outputVector[i][0] = smd_data_9[i][0];
            outputVector[i][1] = smd_data_9[i][1];
            outputVector[i][2] = smd_data_9[i][2];
            outputVector[i][3] = smd_data_9[i][3];
            outputVector[i][4] = smd_data_9[i][4];
            outputVector[i][5] = smd_data_9[i][5];
            outputVector[i][6] = smd_data_9[i][6];
            outputVector[i][7] = smd_data_9[i][7];
            outputVector[i][8] = smd_data_9[i][8];
        } else {
            for (int j = 0; j < size_peratom_cols; j++) {
                outputVector[i][j] = 0.0;
            }
        }
    }
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based array
 ------------------------------------------------------------------------- */

double ComputeSMDTriangleVertices::memory_usage() {
    double bytes = size_peratom_cols * nmax * sizeof(double);
    return bytes;
}
