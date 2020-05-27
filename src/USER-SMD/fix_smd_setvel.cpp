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

#include "fix_smd_setvel.h"
#include <mpi.h>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum {
        NONE, CONSTANT, EQUAL, ATOM
};

/* ---------------------------------------------------------------------- */

FixSMDSetVel::FixSMDSetVel(LAMMPS *lmp, int narg, char **arg) :
                Fix(lmp, narg, arg) {
        if (narg < 6)
                error->all(FLERR, "Illegal fix setvelocity command");

        dynamic_group_allow = 1;
        vector_flag = 1;
        size_vector = 3;
        global_freq = 1;
        extvector = 1;

        xstr = ystr = zstr = NULL;

        if (strstr(arg[3], "v_") == arg[3]) {
                int n = strlen(&arg[3][2]) + 1;
                xstr = new char[n];
                strcpy(xstr, &arg[3][2]);
        } else if (strcmp(arg[3], "NULL") == 0) {
                xstyle = NONE;
        } else {
                xvalue = force->numeric(FLERR, arg[3]);
                xstyle = CONSTANT;
        }
        if (strstr(arg[4], "v_") == arg[4]) {
                int n = strlen(&arg[4][2]) + 1;
                ystr = new char[n];
                strcpy(ystr, &arg[4][2]);
        } else if (strcmp(arg[4], "NULL") == 0) {
                ystyle = NONE;
        } else {
                yvalue = force->numeric(FLERR, arg[4]);
                ystyle = CONSTANT;
        }
        if (strstr(arg[5], "v_") == arg[5]) {
                int n = strlen(&arg[5][2]) + 1;
                zstr = new char[n];
                strcpy(zstr, &arg[5][2]);
        } else if (strcmp(arg[5], "NULL") == 0) {
                zstyle = NONE;
        } else {
                zvalue = force->numeric(FLERR, arg[5]);
                zstyle = CONSTANT;
        }

        // optional args

        iregion = -1;
        idregion = NULL;

        int iarg = 6;
        while (iarg < narg) {
                if (strcmp(arg[iarg], "region") == 0) {
                        if (iarg + 2 > narg)
                                error->all(FLERR, "Illegal fix setvelocity command");
                        iregion = domain->find_region(arg[iarg + 1]);
                        if (iregion == -1)
                                error->all(FLERR, "Region ID for fix setvelocity does not exist");
                        int n = strlen(arg[iarg + 1]) + 1;
                        idregion = new char[n];
                        strcpy(idregion, arg[iarg + 1]);
                        iarg += 2;
                } else
                        error->all(FLERR, "Illegal fix setvelocity command");
        }

        force_flag = 0;
        foriginal[0] = foriginal[1] = foriginal[2] = 0.0;

        maxatom = atom->nmax;
        memory->create(sforce, maxatom, 3, "setvelocity:sforce");
}

/* ---------------------------------------------------------------------- */

FixSMDSetVel::~FixSMDSetVel() {
        delete[] xstr;
        delete[] ystr;
        delete[] zstr;
        delete[] idregion;
        memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixSMDSetVel::setmask() {
        int mask = 0;
        //mask |= INITIAL_INTEGRATE;
        mask |= POST_FORCE;
        return mask;
}

/* ---------------------------------------------------------------------- */

void FixSMDSetVel::init() {
        // check variables

        if (xstr) {
                xvar = input->variable->find(xstr);
                if (xvar < 0)
                        error->all(FLERR, "Variable name for fix setvelocity does not exist");
                if (input->variable->equalstyle(xvar))
                        xstyle = EQUAL;
                else if (input->variable->atomstyle(xvar))
                        xstyle = ATOM;
                else
                        error->all(FLERR, "Variable for fix setvelocity is invalid style");
        }
        if (ystr) {
                yvar = input->variable->find(ystr);
                if (yvar < 0)
                        error->all(FLERR, "Variable name for fix setvelocity does not exist");
                if (input->variable->equalstyle(yvar))
                        ystyle = EQUAL;
                else if (input->variable->atomstyle(yvar))
                        ystyle = ATOM;
                else
                        error->all(FLERR, "Variable for fix setvelocity is invalid style");
        }
        if (zstr) {
                zvar = input->variable->find(zstr);
                if (zvar < 0)
                        error->all(FLERR, "Variable name for fix setvelocity does not exist");
                if (input->variable->equalstyle(zvar))
                        zstyle = EQUAL;
                else if (input->variable->atomstyle(zvar))
                        zstyle = ATOM;
                else
                        error->all(FLERR, "Variable for fix setvelocity is invalid style");
        }

        // set index and check validity of region

        if (iregion >= 0) {
                iregion = domain->find_region(idregion);
                if (iregion == -1)
                        error->all(FLERR, "Region ID for fix setvelocity does not exist");
        }

        if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
                varflag = ATOM;
        else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
                varflag = EQUAL;
        else
                varflag = CONSTANT;

        // cannot use non-zero forces for a minimization since no energy is integrated
        // use fix addforce instead

        int flag = 0;
        if (update->whichflag == 2) {
                if (xstyle == EQUAL || xstyle == ATOM)
                        flag = 1;
                if (ystyle == EQUAL || ystyle == ATOM)
                        flag = 1;
                if (zstyle == EQUAL || zstyle == ATOM)
                        flag = 1;
                if (xstyle == CONSTANT && xvalue != 0.0)
                        flag = 1;
                if (ystyle == CONSTANT && yvalue != 0.0)
                        flag = 1;
                if (zstyle == CONSTANT && zvalue != 0.0)
                        flag = 1;
        }
        if (flag)
                error->all(FLERR, "Cannot use non-zero forces in an energy minimization");
}

/* ---------------------------------------------------------------------- */

void FixSMDSetVel::setup(int vflag) {
        if (strstr(update->integrate_style, "verlet"))
                post_force(vflag);
        else
      error->all(FLERR,"Fix smd/setvel does not support RESPA");
}

/* ---------------------------------------------------------------------- */

void FixSMDSetVel::min_setup(int vflag) {
        post_force(vflag);
}

/* ---------------------------------------------------------------------- */

//void FixSMDSetVel::initial_integrate(int vflag) {
void FixSMDSetVel::post_force(int /*vflag*/) {
        double **x = atom->x;
        double **f = atom->f;
        double **v = atom->v;
        double **vest = atom->vest;
        int *mask = atom->mask;
        int nlocal = atom->nlocal;

        // update region if necessary

        Region *region = NULL;
        if (iregion >= 0) {
                region = domain->regions[iregion];
                region->prematch();
        }

        // reallocate sforce array if necessary

        if (varflag == ATOM && atom->nmax > maxatom) {
                maxatom = atom->nmax;
                memory->destroy(sforce);
                memory->create(sforce, maxatom, 3, "setvelocity:sforce");
        }

        foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
        force_flag = 0;

        if (varflag == CONSTANT) {
                for (int i = 0; i < nlocal; i++)
                        if (mask[i] & groupbit) {
                                if (region && !region->match(x[i][0], x[i][1], x[i][2]))
                                        continue;
                                foriginal[0] += f[i][0];
                                foriginal[1] += f[i][1];
                                foriginal[2] += f[i][2];
                                if (xstyle) {
                                        v[i][0] = xvalue;
                                        vest[i][0] = xvalue;
                                        f[i][0] = 0.0;
                                }
                                if (ystyle) {
                                        v[i][1] = yvalue;
                                        vest[i][1] = yvalue;
                                        f[i][1] = 0.0;
                                }
                                if (zstyle) {
                                        v[i][2] = zvalue;
                                        vest[i][2] = zvalue;
                                        f[i][2] = 0.0;
                                }
                        }

                // variable force, wrap with clear/add

        } else {

                modify->clearstep_compute();

                if (xstyle == EQUAL)
                        xvalue = input->variable->compute_equal(xvar);
                else if (xstyle == ATOM)
                        input->variable->compute_atom(xvar, igroup, &sforce[0][0], 3, 0);
                if (ystyle == EQUAL)
                        yvalue = input->variable->compute_equal(yvar);
                else if (ystyle == ATOM)
                        input->variable->compute_atom(yvar, igroup, &sforce[0][1], 3, 0);
                if (zstyle == EQUAL)
                        zvalue = input->variable->compute_equal(zvar);
                else if (zstyle == ATOM)
                        input->variable->compute_atom(zvar, igroup, &sforce[0][2], 3, 0);

                modify->addstep_compute(update->ntimestep + 1);

                //printf("setting velocity at timestep %d\n", update->ntimestep);

                for (int i = 0; i < nlocal; i++)
                        if (mask[i] & groupbit) {
                                if (region && !region->match(x[i][0], x[i][1], x[i][2]))
                                        continue;
                                foriginal[0] += f[i][0];
                                foriginal[1] += f[i][1];
                                foriginal[2] += f[i][2];
                                if (xstyle == ATOM) {
                                        vest[i][0] = v[i][0] = sforce[i][0];
                                        f[i][0] = 0.0;
                                } else if (xstyle) {
                                        vest[i][0] = v[i][0] = xvalue;
                                        f[i][0] = 0.0;
                                }

                                if (ystyle == ATOM) {
                                        vest[i][1] = v[i][1] = sforce[i][1];
                                        f[i][1] = 0.0;
                                } else if (ystyle) {
                                        vest[i][1] = v[i][1] = yvalue;
                                        f[i][1] = 0.0;
                                }

                                if (zstyle == ATOM) {
                                        vest[i][2] = v[i][2] = sforce[i][2];
                                        f[i][2] = 0.0;
                                } else if (zstyle) {
                                        vest[i][2] = v[i][2] = zvalue;
                                        f[i][2] = 0.0;
                                }

                        }
        }
}

/* ----------------------------------------------------------------------
 return components of total force on fix group before force was changed
 ------------------------------------------------------------------------- */

double FixSMDSetVel::compute_vector(int n) {
// only sum across procs one time

        if (force_flag == 0) {
                MPI_Allreduce(foriginal, foriginal_all, 3, MPI_DOUBLE, MPI_SUM, world);
                force_flag = 1;
        }
        return foriginal_all[n];
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based array
 ------------------------------------------------------------------------- */

double FixSMDSetVel::memory_usage() {
        double bytes = 0.0;
        if (varflag == ATOM)
                bytes = atom->nmax * 3 * sizeof(double);
        return bytes;
}
