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

#include "fix_smd_adjust_dt.h"


#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "fix.h"
#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double BIG = 1.0e20;

/* ---------------------------------------------------------------------- */

FixSMDTlsphDtReset::FixSMDTlsphDtReset(LAMMPS *lmp, int narg, char **arg) :
                Fix(lmp, narg, arg) {
        if (narg != 4) error->all(FLERR, "Illegal fix smd/adjust_dt command");

        // set time_depend, else elapsed time accumulation can be messed up

        time_depend = 1;
        scalar_flag = 1;
        vector_flag = 1;
        size_vector = 2;
        global_freq = 1;
        extscalar = 0;
        extvector = 0;
        restart_global = 1; // this fix stores global (i.e., not per-atom) info: elaspsed time

        safety_factor = utils::numeric(FLERR,arg[3],false,lmp);

        // initializations
        t_elapsed = 0.0;
}

/* ---------------------------------------------------------------------- */

int FixSMDTlsphDtReset::setmask() {
        int mask = 0;
        mask |= INITIAL_INTEGRATE;
        mask |= END_OF_STEP;
        return mask;
}

/* ---------------------------------------------------------------------- */

void FixSMDTlsphDtReset::init() {
        dt = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixSMDTlsphDtReset::setup(int /*vflag*/) {
        end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixSMDTlsphDtReset::initial_integrate(int /*vflag*/) {

        //printf("in adjust_dt: dt = %20.10f\n", update->dt);

        t_elapsed += update->dt;
}

/* ---------------------------------------------------------------------- */

void FixSMDTlsphDtReset::end_of_step() {
        double dtmin = BIG;
        int itmp = 0;

        /*
         * extract minimum CFL timestep from TLSPH and ULSPH pair styles
         */

        auto dtCFL_TLSPH = (double *) force->pair->extract("smd/tlsph/dtCFL_ptr", itmp);
        auto dtCFL_ULSPH = (double *) force->pair->extract("smd/ulsph/dtCFL_ptr", itmp);
        auto dt_TRI = (double *) force->pair->extract("smd/tri_surface/stable_time_increment_ptr", itmp);
        auto dt_HERTZ = (double *) force->pair->extract("smd/hertz/stable_time_increment_ptr", itmp);
        auto dt_PERI_IPMB = (double *) force->pair->extract("smd/peri_ipmb/stable_time_increment_ptr", itmp);

        if ((dtCFL_TLSPH == nullptr) && (dtCFL_ULSPH == nullptr) && (dt_TRI == nullptr) && (dt_HERTZ == nullptr)
                        && (dt_PERI_IPMB == nullptr)) {
                error->all(FLERR, "fix smd/adjust_dt failed to access a valid dtCFL");
        }

        if (dtCFL_TLSPH != nullptr) {
                dtmin = MIN(dtmin, *dtCFL_TLSPH);
        }

        if (dtCFL_ULSPH != nullptr) {
                dtmin = MIN(dtmin, *dtCFL_ULSPH);
        }

        if (dt_TRI != nullptr) {
                dtmin = MIN(dtmin, *dt_TRI);
        }

        if (dt_HERTZ != nullptr) {
                dtmin = MIN(dtmin, *dt_HERTZ);
        }

        if (dt_PERI_IPMB != nullptr) {
                dtmin = MIN(dtmin, *dt_PERI_IPMB);
        }

//      double **v = atom->v;
//      double **f = atom->f;
//      double *rmass = atom->rmass;
//      double *radius = atom->radius;
//      int *mask = atom->mask;
//      int nlocal = atom->nlocal;
//      double dtv, dtf, dtsq;
//      double vsq, fsq, massinv, xmax;
//      double delx, dely, delz, delr;

//      for (int i = 0; i < nlocal; i++) {
//              if (mask[i] & groupbit) {
//                      xmax = 0.005 * radius[i];
//                      massinv = 1.0 / rmass[i];
//                      vsq = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
//                      fsq = f[i][0] * f[i][0] + f[i][1] * f[i][1] + f[i][2] * f[i][2];
//                      dtv = dtf = BIG;
//                      if (vsq > 0.0)
//                              dtv = xmax / sqrt(vsq);
//                      if (fsq > 0.0)
//                              dtf = sqrt(2.0 * xmax / (sqrt(fsq) * massinv));
//                      dt = MIN(dtv, dtf);
//                      dtmin = MIN(dtmin, dt);
//                      dtsq = dt * dt;
//                      delx = dt * v[i][0] + 0.5 * dtsq * massinv * f[i][0];
//                      dely = dt * v[i][1] + 0.5 * dtsq * massinv * f[i][1];
//                      delz = dt * v[i][2] + 0.5 * dtsq * massinv * f[i][2];
//                      delr = sqrt(delx * delx + dely * dely + delz * delz);
//                      if (delr > xmax)
//                              dt *= xmax / delr;
//                      dtmin = MIN(dtmin, dt);
//
////                    xmax = 0.05 * radius[i];
////                    massinv = 1.0 / rmass[i];
////                    fsq = f[i][0] * f[i][0] + f[i][1] * f[i][1] + f[i][2] * f[i][2];
////                    dtf = BIG;
////                    if (fsq > 0.0)
////                            dtf = sqrt(2.0 * xmax / (sqrt(fsq) * massinv));
////                    dtmin = MIN(dtmin, dtf);
//              }
//      }

        dtmin *= safety_factor; // apply safety factor
        MPI_Allreduce(&dtmin, &dt, 1, MPI_DOUBLE, MPI_MIN, world);

        if (update->ntimestep == 0) {
                dt = 1.0e-16;
        }

        //printf("dtmin is now: %f, dt is now%f\n", dtmin, dt);


        update->dt = dt;
        update->dt_default = 0;
        if (force->pair)
                force->pair->reset_dt();
        for (int i = 0; i < modify->nfix; i++)
                modify->fix[i]->reset_dt();
}

/* ---------------------------------------------------------------------- */

double FixSMDTlsphDtReset::compute_scalar() {
        return t_elapsed;
}

/* ----------------------------------------------------------------------
 pack entire state of Fix into one write
 ------------------------------------------------------------------------- */

void FixSMDTlsphDtReset::write_restart(FILE *fp) {
        int n = 0;
        double list[1];
        list[n++] = t_elapsed;

        if (comm->me == 0) {
                int size = n * sizeof(double);
                fwrite(&size, sizeof(int), 1, fp);
                fwrite(list, sizeof(double), n, fp);
        }
}

/* ----------------------------------------------------------------------
 use state info from restart file to restart the Fix
 ------------------------------------------------------------------------- */

void FixSMDTlsphDtReset::restart(char *buf) {
        int n = 0;
        auto list = (double *) buf;
        t_elapsed = list[n++];
}

