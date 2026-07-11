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

#include "fix_smd_move_triangulated_surface.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "math_const.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <Eigen/Eigen>

using namespace Eigen;
using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace std;

/* ---------------------------------------------------------------------- */

FixSMDMoveTriSurf::FixSMDMoveTriSurf(LAMMPS *lmp, int narg, char **arg) :
                Fix(lmp, narg, arg) {

        if (atom->smd_flag != 1) {
                error->all(FLERR, "fix fix smd/move_tri_surf command requires atom_style smd");
        }

        if (narg < 3)
                error->all(FLERR, "Illegal number of arguments for fix fix smd/move_tri_surf command");

        rotateFlag = linearFlag = wiggleFlag = false;
        wiggle_direction  = 1.0;
        wiggle_max_travel = 0.0;

        int iarg = 3;

        if (comm->me == 0) {
                printf("\n>>========>>========>>========>>========>>========>>========>>========>>========\n");
                printf("fix fix smd/move_tri_surf is active for group: %s \n", arg[1]);
        }
        while (true) {

                if (iarg >= narg) {
                        break;
                }

                if (strcmp(arg[iarg], "*LINEAR") == 0) {
                        linearFlag = true;
                        if (comm->me == 0) {
                                printf("... will move surface in a linear fashion\n");
                        }

                        iarg++;
                        if (iarg == narg) {
                                error->all(FLERR, "expected three floats for velocity following *LINEAR");
                        }
                        vx = utils::numeric(FLERR, arg[iarg],false,lmp);

                        iarg++;
                        if (iarg == narg) {
                                error->all(FLERR, "expected three floats for velocity following *LINEAR");
                        }
                        vy = utils::numeric(FLERR, arg[iarg],false,lmp);

                        iarg++;
                        if (iarg == narg) {
                                error->all(FLERR, "expected three floats for velocity following *LINEAR");
                        }
                        vz = utils::numeric(FLERR, arg[iarg],false,lmp);

                } else if (strcmp(arg[iarg], "*WIGGLE") == 0) {
                        wiggleFlag = true;
                        if (comm->me == 0) {
                                printf("... will move surface in wiggle fashion\n");
                        }

                        iarg++;
                        if (iarg == narg) {
                                error->all(FLERR, "expected 4 floats following *WIGGLE : vx vy vz max_travel");
                        }
                        vx = utils::numeric(FLERR, arg[iarg],false,lmp);

                        iarg++;
                        if (iarg == narg) {
                                error->all(FLERR, "expected 4 floats following *WIGGLE : vx vy vz max_travel");
                        }
                        vy = utils::numeric(FLERR, arg[iarg],false,lmp);

                        iarg++;
                        if (iarg == narg) {
                                error->all(FLERR, "expected 4 floats following *WIGGLE : vx vy vz max_travel");
                        }
                        vz = utils::numeric(FLERR, arg[iarg],false,lmp);

                        iarg++;
                        if (iarg == narg) {
                                error->all(FLERR, "expected 4 floats following *WIGGLE : vx vy vz max_travel");
                        }
                        wiggle_max_travel = utils::numeric(FLERR, arg[iarg],false,lmp);

                } else if (strcmp(arg[iarg], "*ROTATE") == 0) {
                        rotateFlag = true;

                        iarg++;
                        if (iarg == narg) {
                                error->all(FLERR, "expected 7 floats following *ROTATE: origin, rotation axis, and rotation period");
                        }
                        origin(0) = utils::numeric(FLERR, arg[iarg],false,lmp);

                        iarg++;
                        if (iarg == narg) {
                                error->all(FLERR, "expected 7 floats following *ROTATE: origin, rotation axis, and rotation period");
                        }
                        origin(1) = utils::numeric(FLERR, arg[iarg],false,lmp);

                        iarg++;
                        if (iarg == narg) {
                                error->all(FLERR, "expected 7 floats following *ROTATE: origin, rotation axis, and rotation period");
                        }
                        origin(2) = utils::numeric(FLERR, arg[iarg],false,lmp);

                        iarg++;
                        if (iarg == narg) {
                                error->all(FLERR, "expected 7 floats following *ROTATE: origin, rotation axis, and rotation period");
                        }
                        rotation_axis(0) = utils::numeric(FLERR, arg[iarg],false,lmp);

                        iarg++;
                        if (iarg == narg) {
                                error->all(FLERR, "expected 7 floats following *ROTATE: origin, rotation axis, and rotation period");
                        }
                        rotation_axis(1) = utils::numeric(FLERR, arg[iarg],false,lmp);

                        iarg++;
                        if (iarg == narg) {
                                error->all(FLERR, "expected 7 floats following *ROTATE: origin, rotation axis, and rotation period");
                        }
                        rotation_axis(2) = utils::numeric(FLERR, arg[iarg],false,lmp);

                        iarg++;
                        if (iarg == narg) {
                                error->all(FLERR, "expected 7 floats following *ROTATE: origin, rotation axis, and rotation period");
                        }
                        rotation_period = utils::numeric(FLERR, arg[iarg],false,lmp);

                        /*
                         * construct rotation matrix
                         */

                        u_cross(0, 0) = 0.0;
                        u_cross(0, 1) = -rotation_axis(2);
                        u_cross(0, 2) = rotation_axis(1);

                        u_cross(1, 0) = rotation_axis(2);
                        u_cross(1, 1) = 0.0;
                        u_cross(1, 2) = -rotation_axis(0);

                        u_cross(2, 0) = -rotation_axis(1);
                        u_cross(2, 1) = rotation_axis(0);
                        u_cross(2, 2) = 0.0;

                        uxu = rotation_axis * rotation_axis.transpose();

                        if (comm->me == 0) {
                                printf("will rotate with period %f\n", rotation_period);
                        }

                } else {
                        char msg[128];
                        snprintf(msg,128, "Illegal keyword for fix smd/move_tri_surf: %s\n", arg[iarg]);
                        error->all(FLERR, msg);
                }

                iarg++;

        }

        if (comm->me == 0) {
                printf(">>========>>========>>========>>========>>========>>========>>========>>========\n\n");
        }

        // set comm sizes needed by this fix
        comm_forward = 12;

        //atom->add_callback(Atom::GROW);
        //atom->add_callback(Atom::RESTART);

        time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixSMDMoveTriSurf::setmask() {
        int mask = 0;
        mask |= INITIAL_INTEGRATE;
        //mask |= PRE_EXCHANGE;
        return mask;
}

/* ---------------------------------------------------------------------- */

void FixSMDMoveTriSurf::init() {
        dtv = update->dt;
}

/* ----------------------------------------------------------------------
 ------------------------------------------------------------------------- */

void FixSMDMoveTriSurf::initial_integrate(int /*vflag*/) {
        double **x = atom->x;
        double **x0 = atom->x0;
        double **v = atom->v;
        double **vest = atom->vest;
        double **smd_data_9 = atom->smd_data_9;
        tagint *mol = atom->molecule;

        int *mask = atom->mask;
        int nlocal = atom->nlocal;
        double phi;
        int i;
        Matrix3d eye, Rot;
        eye.setIdentity();

        Vector3d v1, v2, v3, n, point, rotated_point, R, vel;

        if (igroup == atom->firstgroup)
                nlocal = atom->nfirst;

        if (linearFlag) { // translate particles
                for (i = 0; i < nlocal; i++) {
                        if (mask[i] & groupbit) {

                                v[i][0] = vx;
                                v[i][1] = vy;
                                v[i][2] = vz;

                                vest[i][0] = vx;
                                vest[i][1] = vy;
                                vest[i][2] = vz;

                                x[i][0] += dtv * vx;
                                x[i][1] += dtv * vy;
                                x[i][2] += dtv * vz;

                                /*
                                 * if this is a triangle, move the vertices as well
                                 */

                                if (mol[i] >= 65535) {
                                        smd_data_9[i][0] += dtv * vx;
                                        smd_data_9[i][1] += dtv * vy;
                                        smd_data_9[i][2] += dtv * vz;

                                        smd_data_9[i][3] += dtv * vx;
                                        smd_data_9[i][4] += dtv * vy;
                                        smd_data_9[i][5] += dtv * vz;

                                        smd_data_9[i][6] += dtv * vx;
                                        smd_data_9[i][7] += dtv * vy;
                                        smd_data_9[i][8] += dtv * vz;
                                }

                        }
                }
        }

        if (wiggleFlag) { // wiggle particles forward and backward

                wiggle_travel += sqrt(vx * vx + vy * vy + vz * vz ) * dtv;
                double wiggle_vx = wiggle_direction * vx;
                double wiggle_vy = wiggle_direction * vy;
                double wiggle_vz = wiggle_direction * vz;

                //printf("wiggle vz is %f, wiggle_max_travel is %f, dir=%f\n", wiggle_vz, wiggle_max_travel, wiggle_direction);

                for (i = 0; i < nlocal; i++) {
                        if (mask[i] & groupbit) {

                                v[i][0] = wiggle_vx;
                                v[i][1] = wiggle_vy;
                                v[i][2] = wiggle_vz;

                                vest[i][0] = wiggle_vx;
                                vest[i][1] = wiggle_vy;
                                vest[i][2] = wiggle_vz;

                                x[i][0] += dtv * wiggle_vx;
                                x[i][1] += dtv * wiggle_vy;
                                x[i][2] += dtv * wiggle_vz;

                                /*
                                 * if this is a triangle, move the vertices as well
                                 */

                                if (mol[i] >= 65535) {
                                        smd_data_9[i][0] += dtv * wiggle_vx;
                                        smd_data_9[i][1] += dtv * wiggle_vy;
                                        smd_data_9[i][2] += dtv * wiggle_vz;

                                        smd_data_9[i][3] += dtv * wiggle_vx;
                                        smd_data_9[i][4] += dtv * wiggle_vy;
                                        smd_data_9[i][5] += dtv * wiggle_vz;

                                        smd_data_9[i][6] += dtv * wiggle_vx;
                                        smd_data_9[i][7] += dtv * wiggle_vy;
                                        smd_data_9[i][8] += dtv * wiggle_vz;
                                }

                        }
                }

                if (wiggle_travel >= wiggle_max_travel) {
                        wiggle_direction *= -1.0;
                        wiggle_travel = 0.0;
                }
        }

        if (rotateFlag) { // rotate particles
                Vector3d xnew, R_new, x_correct;

                /*
                 * rotation angle and matrix form of Rodrigues' rotation formula
                 */

                phi = MY_2PI * dtv / rotation_period;
                //printf("dt=%f, phi =%f, T=%f\n", dtv, phi, rotation_period);
                Rot = cos(phi) * eye + sin(phi) * u_cross + (1.0 - cos(phi)) * uxu;

                for (i = 0; i < nlocal; i++) {
                        if (mask[i] & groupbit) {

                                /*
                                 * generate vector R from origin to point which is to be rotated
                                 */
                                point << x[i][0], x[i][1], x[i][2];
                                R = point - origin;

                                /*
                                 * rotate vector and shift away from origin
                                 */
                                rotated_point = Rot * R + origin;

                                /*
                                 * determine velocity
                                 */
                                vel = (rotated_point - point) / update->dt;

                                /*
                                 * assign new velocities and coordinates
                                 */
                                v[i][0] = vel(0);
                                v[i][1] = vel(1);
                                v[i][2] = vel(2);

                                vest[i][0] = vel(0);
                                vest[i][1] = vel(1);
                                vest[i][2] = vel(2);

                                x[i][0] = rotated_point(0);
                                x[i][1] = rotated_point(1);
                                x[i][2] = rotated_point(2);

                                /*
                                 * if this is a triangle, rotate the vertices as well
                                 */

                                if (mol[i] >= 65535) {

                                        v1 << smd_data_9[i][0], smd_data_9[i][1], smd_data_9[i][2];
                                        R = v1 - origin;
                                        rotated_point = Rot * R + origin;
                                        vel = (rotated_point - v1) / update->dt;
                                        smd_data_9[i][0] = rotated_point(0);
                                        smd_data_9[i][1] = rotated_point(1);
                                        smd_data_9[i][2] = rotated_point(2);
                                        v1 = rotated_point;

                                        v2 << smd_data_9[i][3], smd_data_9[i][4], smd_data_9[i][5];
                                        R = v2 - origin;
                                        rotated_point = Rot * R + origin;
                                        vel = (rotated_point - v2) / update->dt;
                                        smd_data_9[i][3] = rotated_point(0);
                                        smd_data_9[i][4] = rotated_point(1);
                                        smd_data_9[i][5] = rotated_point(2);
                                        v2 = rotated_point;

                                        v3 << smd_data_9[i][6], smd_data_9[i][7], smd_data_9[i][8];
                                        R = v3 - origin;
                                        rotated_point = Rot * R + origin;
                                        vel = (rotated_point - v3) / update->dt;
                                        smd_data_9[i][6] = rotated_point(0);
                                        smd_data_9[i][7] = rotated_point(1);
                                        smd_data_9[i][8] = rotated_point(2);
                                        v3 = rotated_point;

                                        // recalculate triangle normal
                                        n = (v2 - v1).cross(v2 - v3);
                                        n /= n.norm();
                                        x0[i][0] = n(0);
                                        x0[i][1] = n(1);
                                        x0[i][2] = n(2);

                                }

                        }
                }
        }

        // we changed smd_data_9, x0. perform communication to ghosts
        comm->forward_comm(this);

}

/* ---------------------------------------------------------------------- */

void FixSMDMoveTriSurf::reset_dt() {
        dtv = update->dt;
}

/* ---------------------------------------------------------------------- */

int FixSMDMoveTriSurf::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/) {
        int i, j, m;
        double **x0 = atom->x0;
        double **smd_data_9 = atom->smd_data_9;

        //printf("in FixSMDIntegrateTlsph::pack_forward_comm\n");
        m = 0;
        for (i = 0; i < n; i++) {
                j = list[i];
                buf[m++] = x0[j][0];
                buf[m++] = x0[j][1];
                buf[m++] = x0[j][2];

                buf[m++] = smd_data_9[j][0];
                buf[m++] = smd_data_9[j][1];
                buf[m++] = smd_data_9[j][2];
                buf[m++] = smd_data_9[j][3];
                buf[m++] = smd_data_9[j][4];
                buf[m++] = smd_data_9[j][5];
                buf[m++] = smd_data_9[j][6];
                buf[m++] = smd_data_9[j][7];
                buf[m++] = smd_data_9[j][8];

        }
        return m;
}

/* ---------------------------------------------------------------------- */

void FixSMDMoveTriSurf::unpack_forward_comm(int n, int first, double *buf) {
        int i, m, last;
        double **x0 = atom->x0;
        double **smd_data_9 = atom->smd_data_9;

        //printf("in FixSMDMoveTriSurf::unpack_forward_comm\n");
        m = 0;
        last = first + n;
        for (i = first; i < last; i++) {
                x0[i][0] = buf[m++];
                x0[i][1] = buf[m++];
                x0[i][2] = buf[m++];

                smd_data_9[i][0] = buf[m++];
                smd_data_9[i][1] = buf[m++];
                smd_data_9[i][2] = buf[m++];
                smd_data_9[i][3] = buf[m++];
                smd_data_9[i][4] = buf[m++];
                smd_data_9[i][5] = buf[m++];
                smd_data_9[i][6] = buf[m++];
                smd_data_9[i][7] = buf[m++];
                smd_data_9[i][8] = buf[m++];
        }
}
