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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "atom_vec_smd.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "force.h"
#include "fix.h"
#include "fix_adapt.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;
#define NMAT_FULL 9
#define NMAT_SYMM 6

/* ---------------------------------------------------------------------- */

AtomVecSMD::AtomVecSMD(LAMMPS *lmp) :
                AtomVec(lmp) {
        molecular = 0;

        comm_x_only = 0;
        comm_f_only = 0;
        size_forward = 6; // variables that are changed by time integration
        size_reverse = 4; // f[3] + de
        size_border = 31;
        size_velocity = 6; // v + vest
        size_data_atom = 13; // 7 + 3 x0 + 3 x
        size_data_vel = 4;
        xcol_data = 11;

        atom->radius_flag = 1;
        atom->rmass_flag = 1;
        atom->vfrac_flag = 1;
        atom->contact_radius_flag = 1;
        atom->molecule_flag = 1;
        atom->smd_data_9_flag = 1;
        atom->e_flag = 1;
        atom->vest_flag = 1;
        atom->smd_stress_flag = 1;
        atom->eff_plastic_strain_flag = 1;
        atom->x0_flag = 1;
        atom->damage_flag = 1;
        atom->eff_plastic_strain_rate_flag = 1;

        forceclearflag = 1;

        atom->smd_flag = 1;
}

/* ---------------------------------------------------------------------- */

void AtomVecSMD::init() {
        AtomVec::init();

        // do nothing here
}

/* ----------------------------------------------------------------------
 grow atom arrays
 n = 0 grows arrays by a chunk
 n > 0 allocates arrays to size n
 ------------------------------------------------------------------------- */

void AtomVecSMD::grow(int n) {
        if (n == 0)
                grow_nmax();
        else
                nmax = n;
        atom->nmax = nmax;
        if (nmax < 0 || nmax > MAXSMALLINT)
                error->one(FLERR, "Per-processor system is too big");

        //printf("in grow, nmax is now %d\n", nmax);

        tag = memory->grow(atom->tag, nmax, "atom:tag");
        type = memory->grow(atom->type, nmax, "atom:type");
        mask = memory->grow(atom->mask, nmax, "atom:mask");
        image = memory->grow(atom->image, nmax, "atom:image");
        x = memory->grow(atom->x, nmax, 3, "atom:x");
        v = memory->grow(atom->v, nmax, 3, "atom:v");

        f = memory->grow(atom->f, nmax * comm->nthreads, 3, "atom:f");
        de = memory->grow(atom->de, nmax * comm->nthreads, "atom:de");

        vfrac = memory->grow(atom->vfrac, nmax, "atom:vfrac");
        rmass = memory->grow(atom->rmass, nmax, "atom:rmass");
        x0 = memory->grow(atom->x0, nmax, 3, "atom:x0");
        radius = memory->grow(atom->radius, nmax, "atom:radius");
        contact_radius = memory->grow(atom->contact_radius, nmax, "atom:contact_radius");
        molecule = memory->grow(atom->molecule, nmax, "atom:molecule");
        smd_data_9 = memory->grow(atom->smd_data_9, nmax, NMAT_FULL, "atom:defgrad_old");
        e = memory->grow(atom->e, nmax, "atom:e");
        vest = memory->grow(atom->vest, nmax, 3, "atom:vest");
        tlsph_stress = memory->grow(atom->smd_stress, nmax, NMAT_SYMM, "atom:tlsph_stress");
        eff_plastic_strain = memory->grow(atom->eff_plastic_strain, nmax, "atom:eff_plastic_strain");
        eff_plastic_strain_rate = memory->grow(atom->eff_plastic_strain_rate, nmax, "atom:eff_plastic_strain_rate");
        damage = memory->grow(atom->damage, nmax, "atom:damage");

        if (atom->nextra_grow)
                for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
                        modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
 reset local array ptrs
 ------------------------------------------------------------------------- */

void AtomVecSMD::grow_reset() {
        tag = atom->tag;
        type = atom->type;
        mask = atom->mask;
        image = atom->image;
        x = atom->x;
        v = atom->v;
        f = atom->f;
        radius = atom->radius;
        rmass = atom->rmass;

        vfrac = atom->vfrac;
        x0 = atom->x0;
        contact_radius = atom->contact_radius;
        molecule = atom->molecule;
        smd_data_9 = atom->smd_data_9;
        e = atom->e;
        de = atom->de;
        tlsph_stress = atom->smd_stress;
        eff_plastic_strain = atom->eff_plastic_strain;
        eff_plastic_strain_rate = atom->eff_plastic_strain_rate;
        damage = atom->damage;
        vest = atom->vest;
}

/* ----------------------------------------------------------------------
 copy atom I info to atom J
 ------------------------------------------------------------------------- */

void AtomVecSMD::copy(int i, int j, int delflag) {
        tag[j] = tag[i];
        type[j] = type[i];
        mask[j] = mask[i];
        image[j] = image[i];
        x[j][0] = x[i][0];
        x[j][1] = x[i][1];
        x[j][2] = x[i][2];
        v[j][0] = v[i][0];
        v[j][1] = v[i][1];
        v[j][2] = v[i][2];

        vfrac[j] = vfrac[i];
        rmass[j] = rmass[i];
        x0[j][0] = x0[i][0];
        x0[j][1] = x0[i][1];
        x0[j][2] = x0[i][2];
        radius[j] = radius[i];
        contact_radius[j] = contact_radius[i];
        molecule[j] = molecule[i];
        e[j] = e[i];
        eff_plastic_strain[j] = eff_plastic_strain[i];
        eff_plastic_strain_rate[j] = eff_plastic_strain_rate[i];
        vest[j][0] = vest[i][0];
        vest[j][1] = vest[i][1];
        vest[j][2] = vest[i][2];

        for (int k = 0; k < NMAT_FULL; k++) {
                smd_data_9[j][k] = smd_data_9[i][k];
        }

        for (int k = 0; k < NMAT_SYMM; k++) {
                tlsph_stress[j][k] = tlsph_stress[i][k];
        }

        damage[j] = damage[i];

        if (atom->nextra_grow)
                for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
                        modify->fix[atom->extra_grow[iextra]]->copy_arrays(i, j, delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVecSMD::pack_comm(int /*n*/, int * /*list*/, double * /*buf*/, int /*pbc_flag*/, int * /*pbc*/) {
        error->one(FLERR, "atom vec tlsph can only be used with ghost velocities turned on");
        return -1;
}

/* ---------------------------------------------------------------------- */

int AtomVecSMD::pack_comm_vel(int n, int *list, double *buf, int pbc_flag, int *pbc) {
        // communicate quantities to ghosts, which are changed by time-integration AND are required on ghost atoms.

        //no need to pack stress or defgrad information here, as these quantities are not required for ghost atoms.
        // Inside pair_style tlsph, these quantities are computed and communicated to ghosts.

        // no need to communicate x0 here, as it is not changed by time integration
        // if x0 is changed when the ref config is updated, this communication is performed in the fix_integrate/tlsph
        // similarily, rmass could be removed here.
        // radius should be communicated here for future time-integration of the radius with ulsph (not implemented yet)
        int i, j, m;
        double dx, dy, dz, dvx, dvy, dvz;

        m = 0;
        if (pbc_flag == 0) {
                for (i = 0; i < n; i++) {
                        j = list[i];
                        buf[m++] = x[j][0];
                        buf[m++] = x[j][1];
                        buf[m++] = x[j][2]; //3
                        buf[m++] = radius[j];
                        buf[m++] = vfrac[j]; // 5
                        buf[m++] = v[j][0];
                        buf[m++] = v[j][1];
                        buf[m++] = v[j][2]; // 8

                        buf[m++] = vest[j][0];
                        buf[m++] = vest[j][1];
                        buf[m++] = vest[j][2]; // 11
                        buf[m++] = e[j]; // 12

                }
        } else {
                if (domain->triclinic == 0) {
                        dx = pbc[0] * domain->xprd;
                        dy = pbc[1] * domain->yprd;
                        dz = pbc[2] * domain->zprd;
                } else {
                        dx = pbc[0] * domain->xprd + pbc[5] * domain->xy + pbc[4] * domain->xz;
                        dy = pbc[1] * domain->yprd + pbc[3] * domain->yz;
                        dz = pbc[2] * domain->zprd;
                }
                if (!deform_vremap) {
                        for (i = 0; i < n; i++) {
                                j = list[i];
                                buf[m++] = x[j][0] + dx;
                                buf[m++] = x[j][1] + dy;
                                buf[m++] = x[j][2] + dz;
                                buf[m++] = radius[j];
                                buf[m++] = vfrac[j];
                                buf[m++] = v[j][0];
                                buf[m++] = v[j][1];
                                buf[m++] = v[j][2]; // 8

                                buf[m++] = vest[j][0];
                                buf[m++] = vest[j][1];
                                buf[m++] = vest[j][2]; // 11
                                buf[m++] = e[j]; // 12

                        }
                } else {
                        dvx = pbc[0] * h_rate[0] + pbc[5] * h_rate[5] + pbc[4] * h_rate[4];
                        dvy = pbc[1] * h_rate[1] + pbc[3] * h_rate[3];
                        dvz = pbc[2] * h_rate[2];
//                      printf("\ndvx = %f, dvy=%f, dvz=%f\n", dvx, dvy, dvz);
//                      printf("dx = %f, dy=%f, dz=%f\n", dx, dy, dz);
                        for (i = 0; i < n; i++) {
                                j = list[i];
                                buf[m++] = x[j][0] + dx;
                                buf[m++] = x[j][1] + dy;
                                buf[m++] = x[j][2] + dz;
                                buf[m++] = radius[j];
                                buf[m++] = vfrac[j];
                                if (mask[i] & deform_groupbit) {
                                        buf[m++] = v[j][0] + dvx;
                                        buf[m++] = v[j][1] + dvy;
                                        buf[m++] = v[j][2] + dvz;
                                        buf[m++] = vest[j][0] + dvx;
                                        buf[m++] = vest[j][1] + dvy;
                                        buf[m++] = vest[j][2] + dvz;
                                } else {
                                        buf[m++] = v[j][0];
                                        buf[m++] = v[j][1];
                                        buf[m++] = v[j][2]; // 8
                                        buf[m++] = vest[j][0];
                                        buf[m++] = vest[j][1];
                                        buf[m++] = vest[j][2]; // 11
                                }

                                buf[m++] = e[j]; // 12

                        }
                }
        }

        return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSMD::pack_comm_hybrid(int n, int *list, double *buf) {
        int i, j, m;

        m = 0;
        for (i = 0; i < n; i++) {
                j = list[i];
                buf[m++] = radius[j];
                buf[m++] = vfrac[j];
                buf[m++] = vest[j][0];
                buf[m++] = vest[j][1];
                buf[m++] = vest[j][2];
                buf[m++] = e[j];
        }
        return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecSMD::unpack_comm(int /*n*/, int /*first*/, double * /*buf*/) {
        error->one(FLERR, "atom vec tlsph can only be used with ghost velocities turned on");
}

/* ---------------------------------------------------------------------- */

void AtomVecSMD::unpack_comm_vel(int n, int first, double *buf) {
        int i, m, last;

        m = 0;
        last = first + n;
        for (i = first; i < last; i++) {
                x[i][0] = buf[m++];
                x[i][1] = buf[m++];
                x[i][2] = buf[m++]; //3
                radius[i] = buf[m++];
                vfrac[i] = buf[m++]; // 5
                v[i][0] = buf[m++];
                v[i][1] = buf[m++];
                v[i][2] = buf[m++]; // 8

                vest[i][0] = buf[m++];
                vest[i][1] = buf[m++];
                vest[i][2] = buf[m++]; // 11
                e[i] = buf[m++];

        }
}

/* ---------------------------------------------------------------------- */

int AtomVecSMD::unpack_comm_hybrid(int n, int first, double *buf) {
        int i, m, last;

        m = 0;
        last = first + n;
        for (i = first; i < last; i++) {
                radius[i] = buf[m++];
                vfrac[i] = buf[m++];
                vest[i][0] = buf[m++];
                vest[i][1] = buf[m++];
                vest[i][2] = buf[m++];
                e[i] = buf[m++];
        }
        return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSMD::pack_reverse(int n, int first, double *buf) {
        int i, m, last;

        printf("in pack_reverse\n");

        m = 0;
        last = first + n;
        for (i = first; i < last; i++) {
                buf[m++] = f[i][0];
                buf[m++] = f[i][1];
                buf[m++] = f[i][2];
                buf[m++] = de[i];
        }
        return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSMD::pack_reverse_hybrid(int n, int first, double *buf) {
        int i, m, last;

        m = 0;
        last = first + n;
        for (i = first; i < last; i++) {
                buf[m++] = de[i];
        }
        return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecSMD::unpack_reverse(int n, int *list, double *buf) {
        int i, j, m;

        m = 0;
        for (i = 0; i < n; i++) {
                j = list[i];
                f[j][0] += buf[m++];
                f[j][1] += buf[m++];
                f[j][2] += buf[m++];
                de[j] += buf[m++];
        }
}

/* ---------------------------------------------------------------------- */

int AtomVecSMD::unpack_reverse_hybrid(int n, int *list, double *buf) {
        int i, j, m;

        m = 0;
        for (i = 0; i < n; i++) {
                j = list[i];
                de[j] += buf[m++];
        }
        return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSMD::pack_border(int /*n*/, int * /*list*/, double * /*buf*/, int /*pbc_flag*/, int * /*pbc*/) {
        error->one(FLERR, "atom vec tlsph can only be used with ghost velocities turned on");
        return -1;
}

/* ---------------------------------------------------------------------- */

int AtomVecSMD::pack_border_vel(int n, int *list, double *buf, int pbc_flag, int *pbc) {
        int i, j, m;
        double dx, dy, dz, dvx, dvy, dvz;

        //printf("AtomVecSMD::pack_border_vel\n");

        m = 0;
        if (pbc_flag == 0) {
                for (i = 0; i < n; i++) {
                        j = list[i];
                        buf[m++] = x[j][0];
                        buf[m++] = x[j][1];
                        buf[m++] = x[j][2]; // 3
                        buf[m++] = x0[j][0];
                        buf[m++] = x0[j][1];
                        buf[m++] = x0[j][2]; // 6
                        buf[m++] = ubuf(tag[j]).d;
                        buf[m++] = ubuf(type[j]).d;
                        buf[m++] = ubuf(mask[j]).d;
                        buf[m++] = ubuf(molecule[j]).d; // 10
                        buf[m++] = radius[j];
                        buf[m++] = rmass[j];
                        buf[m++] = vfrac[j];
                        buf[m++] = contact_radius[j];
                        buf[m++] = e[j];
                        buf[m++] = eff_plastic_strain[j]; // 16

                        for (int k = 0; k < NMAT_FULL; k++) {
                                buf[m++] = smd_data_9[j][k];
                        } // 25

                        for (int k = 0; k < NMAT_SYMM; k++) {
                                buf[m++] = tlsph_stress[j][k];
                        } // 31

                        buf[m++] = v[j][0];
                        buf[m++] = v[j][1];
                        buf[m++] = v[j][2]; // 34
                        buf[m++] = vest[j][0];
                        buf[m++] = vest[j][1];
                        buf[m++] = vest[j][2]; // 37
                }
        } else {

                if (domain->triclinic == 0) {
                        dx = pbc[0] * domain->xprd;
                        dy = pbc[1] * domain->yprd;
                        dz = pbc[2] * domain->zprd;
                } else {
                        dx = pbc[0];
                        dy = pbc[1];
                        dz = pbc[2];
                }
                if (!deform_vremap) {
                        //printf("dx = %f\n", dx);
                        for (i = 0; i < n; i++) {
                                j = list[i];
                                buf[m++] = x[j][0] + dx;
                                buf[m++] = x[j][1] + dy;
                                buf[m++] = x[j][2] + dz; // 3
                                buf[m++] = x0[j][0]; // this is correct
                                buf[m++] = x0[j][1];
                                buf[m++] = x0[j][2]; // 6
                                buf[m++] = ubuf(tag[j]).d;
                                buf[m++] = ubuf(type[j]).d;
                                buf[m++] = ubuf(mask[j]).d;
                                buf[m++] = ubuf(molecule[j]).d; // 10
                                buf[m++] = radius[j];
                                buf[m++] = rmass[j];
                                buf[m++] = vfrac[j];
                                buf[m++] = contact_radius[j];
                                buf[m++] = e[j];
                                buf[m++] = eff_plastic_strain[j]; // 17

                                for (int k = 0; k < NMAT_FULL; k++) {
                                        buf[m++] = smd_data_9[j][k];
                                } // 26

                                for (int k = 0; k < NMAT_SYMM; k++) {
                                        buf[m++] = tlsph_stress[j][k];
                                } // 32

                                buf[m++] = v[j][0];
                                buf[m++] = v[j][1];
                                buf[m++] = v[j][2]; // 35
                                buf[m++] = vest[j][0];
                                buf[m++] = vest[j][1];
                                buf[m++] = vest[j][2]; // 38

                        }
                } else {
                        dvx = pbc[0] * h_rate[0] + pbc[5] * h_rate[5] + pbc[4] * h_rate[4];
                        dvy = pbc[1] * h_rate[1] + pbc[3] * h_rate[3];
                        dvz = pbc[2] * h_rate[2];
//                      printf("\ndvx = %f, dvy=%f, dvz=%f\n", dvx, dvy, dvz);
//                      printf("dx = %f, dy=%f, dz=%f\n", dx, dy, dz);
                        for (i = 0; i < n; i++) {
                                j = list[i];
                                buf[m++] = x[j][0] + dx;
                                buf[m++] = x[j][1] + dy;
                                buf[m++] = x[j][2] + dz; // 3
                                buf[m++] = x0[j][0];
                                buf[m++] = x0[j][1];
                                buf[m++] = x0[j][2]; // 6
                                buf[m++] = ubuf(tag[j]).d;
                                buf[m++] = ubuf(type[j]).d;
                                buf[m++] = ubuf(mask[j]).d;
                                buf[m++] = ubuf(molecule[j]).d; // 10
                                buf[m++] = radius[j];
                                buf[m++] = rmass[j];
                                buf[m++] = vfrac[j];
                                buf[m++] = contact_radius[j];
                                buf[m++] = e[j];
                                buf[m++] = eff_plastic_strain[j]; // 16

                                for (int k = 0; k < NMAT_FULL; k++) {
                                        buf[m++] = smd_data_9[j][k];
                                } // 25

                                for (int k = 0; k < NMAT_SYMM; k++) {
                                        buf[m++] = tlsph_stress[j][k];
                                } // 31

                                if (mask[i] & deform_groupbit) {
                                        buf[m++] = v[j][0] + dvx;
                                        buf[m++] = v[j][1] + dvy;
                                        buf[m++] = v[j][2] + dvz; // 34
                                        buf[m++] = vest[j][0] + dvx;
                                        buf[m++] = vest[j][1] + dvy;
                                        buf[m++] = vest[j][2] + dvz; // 37

                                } else {
                                        buf[m++] = v[j][0];
                                        buf[m++] = v[j][1];
                                        buf[m++] = v[j][2]; // 34
                                        buf[m++] = vest[j][0];
                                        buf[m++] = vest[j][1];
                                        buf[m++] = vest[j][2]; // 37
                                }

                        }
                }
        }

        if (atom->nextra_border)
                for (int iextra = 0; iextra < atom->nextra_border; iextra++)
                        m += modify->fix[atom->extra_border[iextra]]->pack_border(n, list, &buf[m]);

        return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSMD::pack_border_hybrid(int n, int *list, double *buf) {
        int i, j, m;

        m = 0;
        for (i = 0; i < n; i++) {
                j = list[i];

                buf[m++] = x0[j][0];
                buf[m++] = x0[j][1];
                buf[m++] = x0[j][2]; // 3
                buf[m++] = ubuf(molecule[j]).d; // 4
                buf[m++] = radius[j];
                buf[m++] = rmass[j];
                buf[m++] = vfrac[j];
                buf[m++] = contact_radius[j];
                buf[m++] = e[j];
                buf[m++] = eff_plastic_strain[j]; // 11

                for (int k = 0; k < NMAT_FULL; k++) {
                        buf[m++] = smd_data_9[j][k];
                } // 20

                for (int k = 0; k < NMAT_SYMM; k++) {
                        buf[m++] = tlsph_stress[j][k];
                } // 26

        }
        return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecSMD::unpack_border(int /*n*/, int /*first*/, double * /*buf*/) {
        error->one(FLERR, "atom vec tlsph can only be used with ghost velocities turned on");
}

/* ---------------------------------------------------------------------- */

void AtomVecSMD::unpack_border_vel(int n, int first, double *buf) {
        int i, m, last;

        m = 0;
        last = first + n;
        for (i = first; i < last; i++) {
                if (i == nmax)
                        grow(0);
                x[i][0] = buf[m++];
                x[i][1] = buf[m++];
                x[i][2] = buf[m++]; // 3
                x0[i][0] = buf[m++];
                x0[i][1] = buf[m++];
                x0[i][2] = buf[m++]; // 6
                tag[i] = (tagint) ubuf(buf[m++]).i;
                type[i] = (int) ubuf(buf[m++]).i;
                mask[i] = (int) ubuf(buf[m++]).i;
                molecule[i] = (tagint) ubuf(buf[m++]).i; // 10

                radius[i] = buf[m++];
                rmass[i] = buf[m++];
                vfrac[i] = buf[m++];
                contact_radius[i] = buf[m++];
                e[i] = buf[m++];
                eff_plastic_strain[i] = buf[m++]; // 16

                for (int k = 0; k < NMAT_FULL; k++) {
                        smd_data_9[i][k] = buf[m++];
                } // 25

                for (int k = 0; k < NMAT_SYMM; k++) {
                        tlsph_stress[i][k] = buf[m++];
                } // 31

                v[i][0] = buf[m++];
                v[i][1] = buf[m++];
                v[i][2] = buf[m++]; // 34
                vest[i][0] = buf[m++];
                vest[i][1] = buf[m++];
                vest[i][2] = buf[m++]; // 37
        }

        if (atom->nextra_border)
                for (int iextra = 0; iextra < atom->nextra_border; iextra++)
                        m += modify->fix[atom->extra_border[iextra]]->unpack_border(n, first, &buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecSMD::unpack_border_hybrid(int n, int first, double *buf) {
        int i, m, last;

        m = 0;
        last = first + n;
        for (i = first; i < last; i++) {
                x0[i][0] = buf[m++];
                x0[i][1] = buf[m++];
                x0[i][2] = buf[m++]; // 3
                molecule[i] = (tagint) ubuf(buf[m++]).i; // 4
                radius[i] = buf[m++];
                rmass[i] = buf[m++];
                vfrac[i] = buf[m++];
                contact_radius[i] = buf[m++];
                e[i] = buf[m++];
                eff_plastic_strain[i] = buf[m++]; // 11

                for (int k = 0; k < NMAT_FULL; k++) {
                        smd_data_9[i][k] = buf[m++];
                } // 20

                for (int k = 0; k < NMAT_SYMM; k++) {
                        tlsph_stress[i][k] = buf[m++];
                } // 26
        }
        return m;
}

/* ----------------------------------------------------------------------
 pack data for atom I for sending to another proc
 xyz must be 1st 3 values, so comm::exchange() can test on them
 ------------------------------------------------------------------------- */

int AtomVecSMD::pack_exchange(int i, double *buf) {
        int m = 1;

        //printf("in AtomVecSMD::pack_exchange tag %d\n", tag[i]);

        buf[m++] = x[i][0];
        buf[m++] = x[i][1];
        buf[m++] = x[i][2]; // 3
        buf[m++] = x0[i][0];
        buf[m++] = x0[i][1];
        buf[m++] = x0[i][2]; // 6
        buf[m++] = ubuf(tag[i]).d;
        buf[m++] = ubuf(type[i]).d;
        buf[m++] = ubuf(mask[i]).d;
        buf[m++] = ubuf(image[i]).d;
        buf[m++] = ubuf(molecule[i]).d; // 11
        buf[m++] = radius[i];
        buf[m++] = rmass[i];
        buf[m++] = vfrac[i];
        buf[m++] = contact_radius[i];
        buf[m++] = e[i];
        buf[m++] = eff_plastic_strain[i]; // 18
        buf[m++] = eff_plastic_strain_rate[i]; // 19

        for (int k = 0; k < NMAT_FULL; k++) {
                buf[m++] = smd_data_9[i][k];
        } // 27

        for (int k = 0; k < NMAT_SYMM; k++) {
                buf[m++] = tlsph_stress[i][k];
        } // 33

        buf[m++] = v[i][0];
        buf[m++] = v[i][1];
        buf[m++] = v[i][2]; // 36
        buf[m++] = vest[i][0];
        buf[m++] = vest[i][1];
        buf[m++] = vest[i][2]; // 39

        buf[m++] = damage[i];

        if (atom->nextra_grow)
                for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
                        m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i, &buf[m]);

        buf[0] = m;
        return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSMD::unpack_exchange(double *buf) {
        int nlocal = atom->nlocal;
        if (nlocal == nmax)
                grow(0);

        int m = 1;

        x[nlocal][0] = buf[m++];
        x[nlocal][1] = buf[m++];
        x[nlocal][2] = buf[m++]; // 3
        x0[nlocal][0] = buf[m++];
        x0[nlocal][1] = buf[m++];
        x0[nlocal][2] = buf[m++]; // 6
        tag[nlocal] = (tagint) ubuf(buf[m++]).i;
        type[nlocal] = (int) ubuf(buf[m++]).i;
        mask[nlocal] = (int) ubuf(buf[m++]).i;
        image[nlocal] = (imageint) ubuf(buf[m++]).i;
        molecule[nlocal] = (tagint) ubuf(buf[m++]).i; // 11

        radius[nlocal] = buf[m++];
        rmass[nlocal] = buf[m++];
        vfrac[nlocal] = buf[m++];
        contact_radius[nlocal] = buf[m++];
        e[nlocal] = buf[m++];
        eff_plastic_strain[nlocal] = buf[m++]; // 18
        eff_plastic_strain_rate[nlocal] = buf[m++]; // 19

        for (int k = 0; k < NMAT_FULL; k++) {
                smd_data_9[nlocal][k] = buf[m++];
        } // 27

        for (int k = 0; k < NMAT_SYMM; k++) {
                tlsph_stress[nlocal][k] = buf[m++];
        } // 33

        v[nlocal][0] = buf[m++];
        v[nlocal][1] = buf[m++];
        v[nlocal][2] = buf[m++]; // 36
        vest[nlocal][0] = buf[m++];
        vest[nlocal][1] = buf[m++];
        vest[nlocal][2] = buf[m++]; // 39

        damage[nlocal] = buf[m++]; //40

        if (atom->nextra_grow)
                for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
                        m += modify->fix[atom->extra_grow[iextra]]->unpack_exchange(nlocal, &buf[m]);

        atom->nlocal++;
        return m;
}

/* ----------------------------------------------------------------------
 size of restart data for all atoms owned by this proc
 include extra data stored by fixes
 ------------------------------------------------------------------------- */

int AtomVecSMD::size_restart() {
        int i;

        int nlocal = atom->nlocal;
        int n = 43 * nlocal; // count pack_restart + 1 (size of buffer)

        if (atom->nextra_restart)
                for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
                        for (i = 0; i < nlocal; i++)
                                n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

        return n;
}

/* ----------------------------------------------------------------------
 pack atom I's data for restart file including extra quantities
 xyz must be 1st 3 values, so that read_restart can test on them
 molecular types may be negative, but write as positive
 ------------------------------------------------------------------------- */
int AtomVecSMD::pack_restart(int i, double *buf) {
        int m = 1; // 1

        buf[m++] = x[i][0];
        buf[m++] = x[i][1];
        buf[m++] = x[i][2]; // 4
        buf[m++] = x0[i][0];
        buf[m++] = x0[i][1];
        buf[m++] = x0[i][2]; // 7
        buf[m++] = ubuf(tag[i]).d;
        buf[m++] = ubuf(type[i]).d;
        buf[m++] = ubuf(mask[i]).d; // 10
        buf[m++] = ubuf(image[i]).d;
        buf[m++] = ubuf(molecule[i]).d;
        buf[m++] = radius[i];
        buf[m++] = rmass[i];
        buf[m++] = vfrac[i]; // 15
        buf[m++] = contact_radius[i];
        buf[m++] = e[i];
        buf[m++] = eff_plastic_strain[i];
        buf[m++] = eff_plastic_strain_rate[i]; // 19

        for (int k = 0; k < NMAT_FULL; k++) {
                buf[m++] = smd_data_9[i][k];
        } // 28

        for (int k = 0; k < NMAT_SYMM; k++) {
                buf[m++] = tlsph_stress[i][k];
        } // 34

        buf[m++] = v[i][0];
        buf[m++] = v[i][1];
        buf[m++] = v[i][2]; // 37
        buf[m++] = vest[i][0];
        buf[m++] = vest[i][1];
        buf[m++] = vest[i][2]; // 40

        buf[m++] = damage[i]; // 41

        if (atom->nextra_restart)
                for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
                        m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i, &buf[m]);

        buf[0] = m;
        return m;
}

/* ----------------------------------------------------------------------
 unpack data for one atom from restart file including extra quantities
 ------------------------------------------------------------------------- */

int AtomVecSMD::unpack_restart(double *buf) {
        int nlocal = atom->nlocal;
        if (nlocal == nmax) {
                grow(0);
                if (atom->nextra_store)
                        memory->grow(atom->extra, nmax, atom->nextra_store, "atom:extra");
        }

        int m = 1;

        x[nlocal][0] = buf[m++];
        x[nlocal][1] = buf[m++];
        x[nlocal][2] = buf[m++]; // 3
        x0[nlocal][0] = buf[m++];
        x0[nlocal][1] = buf[m++];
        x0[nlocal][2] = buf[m++]; // 6
        tag[nlocal] = (tagint) ubuf(buf[m++]).i;
        type[nlocal] = (int) ubuf(buf[m++]).i;
        mask[nlocal] = (int) ubuf(buf[m++]).i;
        image[nlocal] = (imageint) ubuf(buf[m++]).i;
        molecule[nlocal] = (tagint) ubuf(buf[m++]).i; // 11

        radius[nlocal] = buf[m++];
        rmass[nlocal] = buf[m++];
        vfrac[nlocal] = buf[m++]; //14
        contact_radius[nlocal] = buf[m++]; //15
        e[nlocal] = buf[m++];
        eff_plastic_strain[nlocal] = buf[m++]; // 18
        eff_plastic_strain_rate[nlocal] = buf[m++]; // 29

        for (int k = 0; k < NMAT_FULL; k++) {
                smd_data_9[nlocal][k] = buf[m++];
        } // 28

        for (int k = 0; k < NMAT_SYMM; k++) {
                tlsph_stress[nlocal][k] = buf[m++];
        } // 34

        v[nlocal][0] = buf[m++];
        v[nlocal][1] = buf[m++];
        v[nlocal][2] = buf[m++]; // 37
        vest[nlocal][0] = buf[m++];
        vest[nlocal][1] = buf[m++];
        vest[nlocal][2] = buf[m++]; // 40

        damage[nlocal] = buf[m++]; //41

        //printf("nlocal in restart is %d\n", nlocal);

        double **extra = atom->extra;
        if (atom->nextra_store) {
                int size = static_cast<int>(buf[0]) - m;
                for (int i = 0; i < size; i++)
                        extra[nlocal][i] = buf[m++];
        }

        atom->nlocal++;

        //printf("returning m=%d in unpack_restart\n", m);

        return m;
}

/* ----------------------------------------------------------------------
 create one atom of itype at coord
 set other values to defaults
 ------------------------------------------------------------------------- */

void AtomVecSMD::create_atom(int itype, double *coord) {
        int nlocal = atom->nlocal;
        if (nlocal == nmax) {
                printf("nlocal = %d, nmax = %d, calling grow\n", nlocal, nmax);
                grow(0);
                printf("... finished growing\n");
        }

        tag[nlocal] = 0;
        type[nlocal] = itype;
        x[nlocal][0] = coord[0];
        x[nlocal][1] = coord[1];
        x[nlocal][2] = coord[2];
        x0[nlocal][0] = coord[0];
        x0[nlocal][1] = coord[1];
        x0[nlocal][2] = coord[2];
        mask[nlocal] = 1;
        image[nlocal] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
        v[nlocal][0] = 0.0;
        v[nlocal][1] = 0.0;
        v[nlocal][2] = 0.0;
        vest[nlocal][0] = 0.0;
        vest[nlocal][1] = 0.0;
        vest[nlocal][2] = 0.0;

        vfrac[nlocal] = 1.0;
        rmass[nlocal] = 1.0;
        radius[nlocal] = 0.5;
        contact_radius[nlocal] = 0.5;
        molecule[nlocal] = 1;
        e[nlocal] = 0.0;
        eff_plastic_strain[nlocal] = 0.0;
        eff_plastic_strain_rate[nlocal] = 0.0;

        for (int k = 0; k < NMAT_FULL; k++) {
                smd_data_9[nlocal][k] = 0.0;
        }
        smd_data_9[nlocal][0] = 1.0; // xx
        smd_data_9[nlocal][4] = 1.0; // yy
        smd_data_9[nlocal][8] = 1.0; // zz

        for (int k = 0; k < NMAT_SYMM; k++) {
                tlsph_stress[nlocal][k] = 0.0;
        }

        damage[nlocal] = 0.0;

        atom->nlocal++;
}

/* ----------------------------------------------------------------------
 unpack one line from Atoms section of data file
 initialize other atom quantities
 ------------------------------------------------------------------------- */

void AtomVecSMD::data_atom(double *coord, imageint imagetmp, char **values) {
        int nlocal = atom->nlocal;
        if (nlocal == nmax)
                grow(0);

        tag[nlocal] = ATOTAGINT(values[0]);

        type[nlocal] = force->inumeric(FLERR,values[1]);
        if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
                error->one(FLERR, "Invalid atom type in Atoms section of data file");

        molecule[nlocal] = ATOTAGINT(values[2]);
        if (molecule[nlocal] <= 0)
                error->one(FLERR, "Invalid molecule in Atoms section of data file");

        vfrac[nlocal] = force->numeric(FLERR,values[3]);
        if (vfrac[nlocal] < 0.0)
                error->one(FLERR, "Invalid volume in Atoms section of data file");

        rmass[nlocal] = force->numeric(FLERR,values[4]);
        if (rmass[nlocal] == 0.0)
                error->one(FLERR, "Invalid mass in Atoms section of data file");

        radius[nlocal] = force->numeric(FLERR,values[5]);
        if (radius[nlocal] < 0.0)
                error->one(FLERR, "Invalid radius in Atoms section of data file");

        contact_radius[nlocal] = force->numeric(FLERR,values[6]);
        if (contact_radius[nlocal] < 0.0)
                error->one(FLERR, "Invalid contact radius in Atoms section of data file");

        e[nlocal] = 0.0;

        x0[nlocal][0] = force->numeric(FLERR,values[7]);
        x0[nlocal][1] = force->numeric(FLERR,values[8]);
        x0[nlocal][2] = force->numeric(FLERR,values[9]);

        x[nlocal][0] = coord[0];
        x[nlocal][1] = coord[1];
        x[nlocal][2] = coord[2];

        image[nlocal] = imagetmp;

        mask[nlocal] = 1;
        v[nlocal][0] = 0.0;
        v[nlocal][1] = 0.0;
        v[nlocal][2] = 0.0;
        vest[nlocal][0] = 0.0;
        vest[nlocal][1] = 0.0;
        vest[nlocal][2] = 0.0;

        damage[nlocal] = 0.0;

        eff_plastic_strain[nlocal] = 0.0;
        eff_plastic_strain_rate[nlocal] = 0.0;

        for (int k = 0; k < NMAT_FULL; k++) {
                smd_data_9[nlocal][k] = 0.0;
        }

        for (int k = 0; k < NMAT_SYMM; k++) {
                tlsph_stress[nlocal][k] = 0.0;
        }

        smd_data_9[nlocal][0] = 1.0; // xx
        smd_data_9[nlocal][4] = 1.0; // yy
        smd_data_9[nlocal][8] = 1.0; // zz

        atom->nlocal++;
}

/* ----------------------------------------------------------------------
 unpack hybrid quantities from one line in Atoms section of data file
 initialize other atom quantities for this sub-style
 ------------------------------------------------------------------------- */

int AtomVecSMD::data_atom_hybrid(int /*nlocal*/, char **/*values*/) {
        error->one(FLERR, "hybrid atom style functionality not yet implemented for atom style smd");
        return -1;
}

/* ----------------------------------------------------------------------
 unpack one line from Velocities section of data file
 ------------------------------------------------------------------------- */

void AtomVecSMD::data_vel(int m, char **values) {
        v[m][0] = force->numeric(FLERR,values[0]);
        v[m][1] = force->numeric(FLERR,values[1]);
        v[m][2] = force->numeric(FLERR,values[2]);
        vest[m][0] = force->numeric(FLERR,values[0]);
        vest[m][1] = force->numeric(FLERR,values[1]);
        vest[m][2] = force->numeric(FLERR,values[2]);
}

/* ----------------------------------------------------------------------
 unpack hybrid quantities from one line in Velocities section of data file
 ------------------------------------------------------------------------- */

int AtomVecSMD::data_vel_hybrid(int /*m*/, char **/*values*/) {
        error->one(FLERR, "hybrid atom style functionality not yet implemented for atom style smd");
        return 0;
}

/* ----------------------------------------------------------------------
 pack atom info for data file including 3 image flags
 ------------------------------------------------------------------------- */

void AtomVecSMD::pack_data(double **buf) {
        int nlocal = atom->nlocal;
        for (int i = 0; i < nlocal; i++) {
                buf[i][0] = ubuf(tag[i]).d;
                buf[i][1] = ubuf(type[i]).d;
                buf[i][2] = ubuf(molecule[i]).d;
                buf[i][3] = vfrac[i];
                buf[i][4] = rmass[i];
                buf[i][5] = radius[i];
                buf[i][6] = contact_radius[i];

                buf[i][7] = x[i][0];
                buf[i][8] = x[i][1];
                buf[i][9] = x[i][2];

                buf[i][10] = x0[i][0];
                buf[i][11] = x0[i][1];
                buf[i][12] = x0[i][2];

                buf[i][13] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
                buf[i][14] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
                buf[i][15] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
        }
}

/* ----------------------------------------------------------------------
 pack hybrid atom info for data file
 ------------------------------------------------------------------------- */

int AtomVecSMD::pack_data_hybrid(int /*i*/, double * /*buf*/) {
        error->one(FLERR, "hybrid atom style functionality not yet implemented for atom style smd");
        return -1;
}

/* ----------------------------------------------------------------------
 write atom info to data file including 3 image flags
 ------------------------------------------------------------------------- */

void AtomVecSMD::write_data(FILE *fp, int n, double **buf) {
        for (int i = 0; i < n; i++)
                fprintf(fp,
                TAGINT_FORMAT
                " %d %d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e\n", (tagint) ubuf(buf[i][0]).i,
                                (int) ubuf(buf[i][1]).i, (int) ubuf(buf[i][2]).i, buf[i][3], buf[i][4], buf[i][5], buf[i][6], buf[i][7], buf[i][8],
                                buf[i][9], buf[i][10], buf[i][11], buf[i][12]);
}

/* ----------------------------------------------------------------------
 write hybrid atom info to data file
 ------------------------------------------------------------------------- */

int AtomVecSMD::write_data_hybrid(FILE * /*fp*/, double * /*buf*/) {
        error->one(FLERR, "hybrid atom style functionality not yet implemented for atom style smd");
        return -1;
}

/* ----------------------------------------------------------------------
 pack velocity info for data file
 ------------------------------------------------------------------------- */

void AtomVecSMD::pack_vel(double **buf) {
        int nlocal = atom->nlocal;
        for (int i = 0; i < nlocal; i++) {
                buf[i][0] = ubuf(tag[i]).d;
                buf[i][1] = v[i][0];
                buf[i][2] = v[i][1];
                buf[i][3] = v[i][2];
        }
}

/* ----------------------------------------------------------------------
 pack hybrid velocity info for data file
 ------------------------------------------------------------------------- */

int AtomVecSMD::pack_vel_hybrid(int /*i*/, double * /*buf*/) {
        error->one(FLERR, "hybrid atom style functionality not yet implemented for atom style smd");
        return 0;
}

/* ----------------------------------------------------------------------
 write velocity info to data file
 ------------------------------------------------------------------------- */

void AtomVecSMD::write_vel(FILE *fp, int n, double **buf) {
        for (int i = 0; i < n; i++)
                fprintf(fp, TAGINT_FORMAT
                " %-1.16e %-1.16e %-1.16e\n", (tagint) ubuf(buf[i][0]).i, buf[i][1], buf[i][2], buf[i][3]);
}

/* ----------------------------------------------------------------------
 write hybrid velocity info to data file
 ------------------------------------------------------------------------- */

int AtomVecSMD::write_vel_hybrid(FILE * /*fp*/, double * /*buf*/) {
        error->one(FLERR, "hybrid atom style functionality not yet implemented for atom style smd");
        return 3;
}

/* ----------------------------------------------------------------------
 return # of bytes of allocated memory
 ------------------------------------------------------------------------- */

bigint AtomVecSMD::memory_usage() {
        bigint bytes = 0;

        if (atom->memcheck("tag"))
                bytes += memory->usage(tag, nmax);
        if (atom->memcheck("type"))
                bytes += memory->usage(type, nmax);
        if (atom->memcheck("molecule"))
                bytes += memory->usage(molecule, nmax);
        if (atom->memcheck("mask"))
                bytes += memory->usage(mask, nmax);
        if (atom->memcheck("image"))
                bytes += memory->usage(image, nmax);
        if (atom->memcheck("x"))
                bytes += memory->usage(x, nmax, 3);
        if (atom->memcheck("v"))
                bytes += memory->usage(v, nmax, 3);
        if (atom->memcheck("vest"))
                bytes += memory->usage(vest, nmax, 3);
        if (atom->memcheck("f"))
                bytes += memory->usage(f, nmax * comm->nthreads, 3);

        if (atom->memcheck("radius"))
                bytes += memory->usage(radius, nmax);
        if (atom->memcheck("contact_radius"))
                bytes += memory->usage(contact_radius, nmax);
        if (atom->memcheck("vfrac"))
                bytes += memory->usage(vfrac, nmax);
        if (atom->memcheck("rmass"))
                bytes += memory->usage(rmass, nmax);
        if (atom->memcheck("eff_plastic_strain"))
                bytes += memory->usage(eff_plastic_strain, nmax);
        if (atom->memcheck("eff_plastic_strain_rate"))
                bytes += memory->usage(eff_plastic_strain_rate, nmax);
        if (atom->memcheck("e"))
                bytes += memory->usage(e, nmax);
        if (atom->memcheck("de"))
                bytes += memory->usage(de, nmax);

        if (atom->memcheck("smd_data_9"))
                bytes += memory->usage(smd_data_9, nmax, NMAT_FULL);
        if (atom->memcheck("tlsph_stress"))
                bytes += memory->usage(tlsph_stress, nmax, NMAT_SYMM);

        if (atom->memcheck("damage"))
                bytes += memory->usage(damage, nmax);

        return bytes;
}

/* ---------------------------------------------------------------------- */

void AtomVecSMD::force_clear(int n, size_t nbytes) {
        //printf("clearing force on atom %d", n);
        memset(&de[n], 0, nbytes);
        memset(&f[0][0], 0, 3 * nbytes);
}
