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

#include "compute_smd_tlsph_shape.h"
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
using namespace std;
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSmdTlsphShape::ComputeSmdTlsphShape(LAMMPS *lmp, int narg, char **arg) :
                Compute(lmp, narg, arg) {
        if (narg != 3)
                error->all(FLERR, "Illegal compute smd/tlsph_strain command");

        peratom_flag = 1;
        size_peratom_cols = 7;

        nmax = 0;
        strainVector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeSmdTlsphShape::~ComputeSmdTlsphShape() {
        memory->sfree(strainVector);
}

/* ---------------------------------------------------------------------- */

void ComputeSmdTlsphShape::init() {

        int count = 0;
        for (int i = 0; i < modify->ncompute; i++)
                if (strcmp(modify->compute[i]->style, "smd/tlsph_strain") == 0)
                        count++;
        if (count > 1 && comm->me == 0)
                error->warning(FLERR, "More than one compute smd/tlsph_strain");
}

/* ---------------------------------------------------------------------- */

void ComputeSmdTlsphShape::compute_peratom() {
        double *contact_radius = atom->contact_radius;
        invoked_peratom = update->ntimestep;

        // grow vector array if necessary

        if (atom->nmax > nmax) {
                memory->destroy(strainVector);
                nmax = atom->nmax;
                memory->create(strainVector, nmax, size_peratom_cols, "strainVector");
                array_atom = strainVector;
        }

        int itmp = 0;
        Matrix3d *R = (Matrix3d *) force->pair->extract("smd/tlsph/rotation_ptr", itmp);
        if (R == NULL) {
                error->all(FLERR, "compute smd/tlsph_shape failed to access rotation array");
        }

        Matrix3d *F = (Matrix3d *) force->pair->extract("smd/tlsph/Fincr_ptr", itmp);
        if (F == NULL) {
                error->all(FLERR, "compute smd/tlsph_shape failed to access deformation gradient array");
        }

        int *mask = atom->mask;
        int nlocal = atom->nlocal;
        Matrix3d E, eye;
        eye.setIdentity();
        Quaterniond q;

        for (int i = 0; i < nlocal; i++) {
                if (mask[i] & groupbit) {

                        E = 0.5 * (F[i].transpose() * F[i] - eye); // Green-Lagrange strain
                        strainVector[i][0] = contact_radius[i] * (1.0 + E(0, 0));
                        strainVector[i][1] = contact_radius[i] * (1.0 + E(1, 1));
                        strainVector[i][2] = contact_radius[i] * (1.0 + E(2, 2));

                        q = R[i]; // convert pure rotation matrix to quaternion
                        strainVector[i][3] = q.w();
                        strainVector[i][4] = q.x();
                        strainVector[i][5] = q.y();
                        strainVector[i][6] = q.z();
                } else {
                        for (int j = 0; j < size_peratom_cols; j++) {
                                strainVector[i][j] = 0.0;
                        }
                }
        }
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based array
 ------------------------------------------------------------------------- */

double ComputeSmdTlsphShape::memory_usage() {
        double bytes = size_peratom_cols * nmax * sizeof(double);
        return bytes;
}
