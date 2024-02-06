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

/* ----------------------------------------------------------------------
   Contributing author: Mike Parks (SNL)
------------------------------------------------------------------------- */

#include "pair_smd_triangulated_surface.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>
#include <Eigen/Eigen>

using namespace std;
using namespace LAMMPS_NS;
using namespace Eigen;

/* ---------------------------------------------------------------------- */

PairTriSurf::PairTriSurf(LAMMPS *lmp) :
                Pair(lmp) {

        onerad_dynamic = onerad_frozen = maxrad_dynamic = maxrad_frozen = nullptr;
        bulkmodulus = nullptr;
        kn = nullptr;
        scale = 1.0;
}

/* ---------------------------------------------------------------------- */

PairTriSurf::~PairTriSurf() {

        if (allocated) {
                memory->destroy(setflag);
                memory->destroy(cutsq);
                memory->destroy(bulkmodulus);
                memory->destroy(kn);

                delete[] onerad_dynamic;
                delete[] onerad_frozen;
                delete[] maxrad_dynamic;
                delete[] maxrad_frozen;
        }
}

/* ---------------------------------------------------------------------- */

void PairTriSurf::compute(int eflag, int vflag) {
        int i, j, ii, jj, inum, jnum, itype, jtype;
        double rsq, r, evdwl, fpair;
        int *ilist, *jlist, *numneigh, **firstneigh;
        double rcut, r_geom, delta, r_tri, r_particle, touch_distance, dt_crit;
        int tri, particle;
        Vector3d normal, x1, x2, x3, x4, x13, x23, x43, w, cp, x4cp, vnew, v_old;
        ;
        Vector3d xi, x_center, dx;
        Matrix2d C;
        Vector2d w2d, rhs;

        evdwl = 0.0;
        ev_init(eflag, vflag);

        tagint *mol = atom->molecule;
        double **f = atom->f;
        double **smd_data_9 = atom->smd_data_9;
        double **x = atom->x;
        double **x0 = atom->x0;
        double **v = atom->v;
        double *rmass = atom->rmass;
        int *type = atom->type;
        int nlocal = atom->nlocal;
        double *radius = atom->contact_radius;
        double rcutSq;
        Vector3d offset;

        int newton_pair = force->newton_pair;
        int periodic = (domain->xperiodic || domain->yperiodic || domain->zperiodic);

        inum = list->inum;
        ilist = list->ilist;
        numneigh = list->numneigh;
        firstneigh = list->firstneigh;

        int max_neighs = 0;
        stable_time_increment = 1.0e22;

        // loop over neighbors of my atoms using a half neighbor list
        for (ii = 0; ii < inum; ii++) {
                i = ilist[ii];
                itype = type[i];
                jlist = firstneigh[i];
                jnum = numneigh[i];
                max_neighs = MAX(max_neighs, jnum);

                for (jj = 0; jj < jnum; jj++) {
                        j = jlist[jj];

                        j &= NEIGHMASK;

                        jtype = type[j];

                        /*
                         * decide which one of i, j is triangle and which is particle
                         */
                        if ((mol[i] < 65535) && (mol[j] >= 65535)) {
                                particle = i;
                                tri = j;
                        } else if ((mol[j] < 65535) && (mol[i] >= 65535)) {
                                particle = j;
                                tri = i;
                        } else {
                                error->one(FLERR, "unknown case");
                        }

                        //x_center << x[tri][0], x[tri][1], x[tri][2]; // center of triangle
                        x_center(0) = x[tri][0];
                        x_center(1) = x[tri][1];
                        x_center(2) = x[tri][2];
                        //x4 << x[particle][0], x[particle][1], x[particle][2];
                        x4(0) = x[particle][0];
                        x4(1) = x[particle][1];
                        x4(2) = x[particle][2];
                        dx = x_center - x4; //
                        if (periodic) {
                                domain->minimum_image(dx(0), dx(1), dx(2));
                        }
                        rsq = dx.squaredNorm();

                        r_tri = scale * radius[tri];
                        r_particle = scale * radius[particle];
                        rcut = r_tri + r_particle;
                        rcutSq = rcut * rcut;

                        //printf("type i=%d, type j=%d, r=%f, ri=%f, rj=%f\n", itype, jtype, sqrt(rsq), ri, rj);

                        if (rsq < rcutSq) {

                                /*
                                 * gather triangle information
                                 */
                                normal(0) = x0[tri][0];
                                normal(1) = x0[tri][1];
                                normal(2) = x0[tri][2];

                                /*
                                 * distance check: is particle closer than its radius to the triangle plane?
                                 */
                                if (fabs(dx.dot(normal)) < radius[particle]) {
                                        /*
                                         * get other two triangle vertices
                                         */
                                        x1(0) = smd_data_9[tri][0];
                                        x1(1) = smd_data_9[tri][1];
                                        x1(2) = smd_data_9[tri][2];
                                        x2(0) = smd_data_9[tri][3];
                                        x2(1) = smd_data_9[tri][4];
                                        x2(2) = smd_data_9[tri][5];
                                        x3(0) = smd_data_9[tri][6];
                                        x3(1) = smd_data_9[tri][7];
                                        x3(2) = smd_data_9[tri][8];

                                        PointTriangleDistance(x4, x1, x2, x3, cp, r);

                                        /*
                                         * distance to closest point
                                         */
                                        x4cp = x4 - cp;

                                        /*
                                         * flip normal to point in direction of x4cp
                                         */

                                        if (x4cp.dot(normal) < 0.0) {
                                                normal *= -1.0;
                                        }

                                        /*
                                         * penalty force pushes particle away from triangle
                                         */
                                        if (r < 1.0 * radius[particle]) {

                                                delta = radius[particle] - r; // overlap distance
                                                r_geom = radius[particle];
                                                fpair = 1.066666667e0 * bulkmodulus[itype][jtype] * delta * sqrt(delta * r_geom);
                                                dt_crit = 3.14 * sqrt(rmass[particle] / (fpair / delta));
                                                stable_time_increment = MIN(stable_time_increment, dt_crit);

                                                evdwl = r * fpair * 0.4e0 * delta; // GCG 25 April: this expression conserves total energy

                                                fpair /= (r + 1.0e-2 * radius[particle]); // divide by r + softening and multiply with non-normalized distance vector

                                                if (particle < nlocal) {
                                                        f[particle][0] += x4cp(0) * fpair;
                                                        f[particle][1] += x4cp(1) * fpair;
                                                        f[particle][2] += x4cp(2) * fpair;
                                                }

                                                if (tri < nlocal) {
                                                        f[tri][0] -= x4cp(0) * fpair;
                                                        f[tri][1] -= x4cp(1) * fpair;
                                                        f[tri][2] -= x4cp(2) * fpair;
                                                }

                                                if (evflag) {
                                                        ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, x4cp(0), x4cp(1), x4cp(2));
                                                }

                                        }

                                        /*
                                         * if particle comes too close to triangle, reflect its velocity and explicitly move it away
                                         */

                                        touch_distance = 1.0 * radius[particle];
                                        if (r < touch_distance) {

                                                /*
                                                 * reflect velocity if it points toward triangle
                                                 */

                                                normal = x4cp / r;

                                                //v_old << v[particle][0], v[particle][1], v[particle][2];
                                                v_old(0) = v[particle][0];
                                                v_old(1) = v[particle][1];
                                                v_old(2) = v[particle][2];
                                                if (v_old.dot(normal) < 0.0) {
                                                        //printf("flipping velocity\n");
                                                        vnew = 1.0 * (-2.0 * v_old.dot(normal) * normal + v_old);
                                                        v[particle][0] = vnew(0);
                                                        v[particle][1] = vnew(1);
                                                        v[particle][2] = vnew(2);
                                                }

                                                //printf("moving particle on top of triangle\n");
                                                x[particle][0] = cp(0) + touch_distance * normal(0);
                                                x[particle][1] = cp(1) + touch_distance * normal(1);
                                                x[particle][2] = cp(2) + touch_distance * normal(2);
                                        }

                                }
                        }
                }
        }

//      int max_neighs_all = 0;
//      MPI_Allreduce(&max_neighs, &max_neighs_all, 1, MPI_INT, MPI_MAX, world);
//      if (comm->me == 0) {
//              printf("max. neighs in tri pair is %d\n", max_neighs_all);
//      }
//
//              double stable_time_increment_all = 0.0;
//              MPI_Allreduce(&stable_time_increment, &stable_time_increment_all, 1, MPI_DOUBLE, MPI_MIN, world);
//              if (comm->me == 0) {
//                      printf("stable time step tri pair is %f\n", stable_time_increment_all);
//              }
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairTriSurf::allocate() {
        allocated = 1;
        int n = atom->ntypes;

        memory->create(setflag, n + 1, n + 1, "pair:setflag");
        for (int i = 1; i <= n; i++)
                for (int j = i; j <= n; j++)
                        setflag[i][j] = 0;

        memory->create(bulkmodulus, n + 1, n + 1, "pair:kspring");
        memory->create(kn, n + 1, n + 1, "pair:kn");

        memory->create(cutsq, n + 1, n + 1, "pair:cutsq"); // always needs to be allocated, even with granular neighborlist

        onerad_dynamic = new double[n + 1];
        onerad_frozen = new double[n + 1];
        maxrad_dynamic = new double[n + 1];
        maxrad_frozen = new double[n + 1];
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairTriSurf::settings(int narg, char **arg) {
        if (narg != 1)
                error->all(FLERR, "Illegal number of args for pair_style smd/tri_surface");

        scale = utils::numeric(FLERR, arg[0],false,lmp);
        if (comm->me == 0) {
                printf("\n>>========>>========>>========>>========>>========>>========>>========>>========\n");
                printf("SMD/TRI_SURFACE CONTACT SETTINGS:\n");
                printf("... effective contact radius is scaled by %f\n", scale);
                printf(">>========>>========>>========>>========>>========>>========>>========>>========\n");
        }

}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairTriSurf::coeff(int narg, char **arg) {
        if (narg != 3)
                error->all(FLERR, "Incorrect args for pair coefficients");
        if (!allocated)
                allocate();

        int ilo, ihi, jlo, jhi;
        utils::bounds(FLERR,arg[0], 1,atom->ntypes, ilo, ihi, error);
        utils::bounds(FLERR,arg[1], 1,atom->ntypes, jlo, jhi, error);

        double bulkmodulus_one = utils::numeric(FLERR,arg[2],false,lmp);

        // set short-range force constant
        double kn_one = 0.0;
        if (domain->dimension == 3) {
                kn_one = (16. / 15.) * bulkmodulus_one; //assuming poisson ratio = 1/4 for 3d
        } else {
                kn_one = 0.251856195 * (2. / 3.) * bulkmodulus_one; //assuming poisson ratio = 1/3 for 2d
        }

        int count = 0;
        for (int i = ilo; i <= ihi; i++) {
                for (int j = MAX(jlo, i); j <= jhi; j++) {
                        bulkmodulus[i][j] = bulkmodulus_one;
                        kn[i][j] = kn_one;
                        setflag[i][j] = 1;
                        count++;
                }
        }

        if (count == 0)
                error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairTriSurf::init_one(int i, int j) {

        if (!allocated)
                allocate();

        if (setflag[i][j] == 0)
                error->all(FLERR, "All pair coeffs are not set");

        bulkmodulus[j][i] = bulkmodulus[i][j];
        kn[j][i] = kn[i][j];

        // cutoff = sum of max I,J radii for
        // dynamic/dynamic & dynamic/frozen interactions, but not frozen/frozen

        double cutoff = maxrad_dynamic[i] + maxrad_dynamic[j];
        cutoff = MAX(cutoff, maxrad_frozen[i] + maxrad_dynamic[j]);
        cutoff = MAX(cutoff, maxrad_dynamic[i] + maxrad_frozen[j]);

        if (comm->me == 0) {
                printf("cutoff for pair smd/smd/tri_surface = %f\n", cutoff);
        }
        return cutoff;
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairTriSurf::init_style() {
        int i;

        // error checks

        if (!atom->contact_radius_flag)
                error->all(FLERR, "Pair style smd/smd/tri_surface requires atom style with contact_radius");

        neighbor->add_request(this, NeighConst::REQ_SIZE);

        // set maxrad_dynamic and maxrad_frozen for each type
        // include future Fix pour particles as dynamic

        for (i = 1; i <= atom->ntypes; i++)
                onerad_dynamic[i] = onerad_frozen[i] = 0.0;

        double *radius = atom->radius;
        int *type = atom->type;
        int nlocal = atom->nlocal;

        for (i = 0; i < nlocal; i++) {
                onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]], radius[i]);
        }

        MPI_Allreduce(&onerad_dynamic[1], &maxrad_dynamic[1], atom->ntypes, MPI_DOUBLE, MPI_MAX, world);
        MPI_Allreduce(&onerad_frozen[1], &maxrad_frozen[1], atom->ntypes, MPI_DOUBLE, MPI_MAX, world);
}

/* ----------------------------------------------------------------------
 neighbor callback to inform pair style of neighbor list to use
 optional granular history list
 ------------------------------------------------------------------------- */

void PairTriSurf::init_list(int id, NeighList *ptr) {
        if (id == 0)
                list = ptr;
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based arrays
 ------------------------------------------------------------------------- */

double PairTriSurf::memory_usage() {

        return 0.0;
}

/*
 * distance between triangle and point
 */
/*
 function [dist,PP0] = pointTriangleDistance(TRI,P)
 % calculate distance between a point and a triangle in 3D
 % SYNTAX
 %   dist = pointTriangleDistance(TRI,P)
 %   [dist,PP0] = pointTriangleDistance(TRI,P)
 %
 % DESCRIPTION
 %   Calculate the distance of a given point P from a triangle TRI.
 %   Point P is a row vector of the form 1x3. The triangle is a matrix
 %   formed by three rows of points TRI = [P1;P2;P3] each of size 1x3.
 %   dist = pointTriangleDistance(TRI,P) returns the distance of the point P
 %   to the triangle TRI.
 %   [dist,PP0] = pointTriangleDistance(TRI,P) additionally returns the
 %   closest point PP0 to P on the triangle TRI.
 %
 % Author: Gwendolyn Fischer
 % Release: 1.0
 % Release date: 09/02/02
 % Release: 1.1 Fixed Bug because of normalization
 % Release: 1.2 Fixed Bug because of typo in region 5 20101013
 % Release: 1.3 Fixed Bug because of typo in region 2 20101014

 % Possible extension could be a version tailored not to return the distance
 % and additionally the closest point, but instead return only the closest
 % point. Could lead to a small speed gain.

 % Example:
 % %% The Problem
 % P0 = [0.5 -0.3 0.5];
 %
 % P1 = [0 -1 0];
 % P2 = [1  0 0];
 % P3 = [0  0 0];
 %
 % vertices = [P1; P2; P3];
 % faces = [1 2 3];
 %
 % %% The Engine
 % [dist,PP0] = pointTriangleDistance([P1;P2;P3],P0);
 %
 % %% Visualization
 % [x,y,z] = sphere(20);
 % x = dist*x+P0(1);
 % y = dist*y+P0(2);
 % z = dist*z+P0(3);
 %
 % figure
 % hold all
 % patch('Vertices',vertices,'Faces',faces,'FaceColor','r','FaceAlpha',0.8);
 % plot3(P0(1),P0(2),P0(3),'b*');
 % plot3(PP0(1),PP0(2),PP0(3),'*g')
 % surf(x,y,z,'FaceColor','b','FaceAlpha',0.3)
 % view(3)

 % The algorithm is based on
 % "David Eberly, 'Distance Between Point and Triangle in 3D',
 % Geometric Tools, LLC, (1999)"
 % https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
 %
 %        ^t
 %  \     |
 %   \reg2|
 %    \   |
 %     \  |
 %      \ |
 %       \|
 %        *P2
 %        |\
%        | \
%  reg3  |  \ reg1
 %        |   \
%        |reg0\
%        |     \
%        |      \ P1
 % -------*-------*------->s
 %        |P0      \
%  reg4  | reg5    \ reg6
 */

//void PairTriSurf::PointTriangleDistance(const Vector3d P, const Vector3d TRI1, const Vector3d TRI2, const Vector3d TRI3,
//              Vector3d &CP, double &dist) {
//
//      Vector3d B, E0, E1, D;
//      double a, b, c, d, e, f;
//      double det, s, t, sqrDistance, tmp0, tmp1, numer, denom, invDet;
//
//      // rewrite triangle in normal form
//      B = TRI1;
//      E0 = TRI2 - B;
//      E1 = TRI3 - B;
//
//      D = B - P;
//      a = E0.dot(E0);
//      b = E0.dot(E1);
//      c = E1.dot(E1);
//      d = E0.dot(D);
//      e = E1.dot(D);
//      f = D.dot(D);
//
//      det = a * c - b * b;
//      //% do we have to use abs here?
//      s = b * e - c * d;
//      t = b * d - a * e;
//
//      //% Terible tree of conditionals to determine in which region of the diagram
//      //% shown above the projection of the point into the triangle-plane lies.
//      if ((s + t) <= det) {
//              if (s < 0) {
//                      if (t < 0) {
//                              // %region4
//                              if (d < 0) {
//                                      t = 0;
//                                      if (-d >= a) {
//                                              s = 1;
//                                              sqrDistance = a + 2 * d + f;
//                                      } else {
//                                              s = -d / a;
//                                              sqrDistance = d * s + f;
//                                      }
//                              } else {
//                                      s = 0;
//                                      if (e >= 0) {
//                                              t = 0;
//                                              sqrDistance = f;
//                                      } else {
//                                              if (-e >= c) {
//                                                      t = 1;
//                                                      sqrDistance = c + 2 * e + f;
//                                              } else {
//                                                      t = -e / c;
//                                                      sqrDistance = e * t + f;
//                                              }
//                                      }
//                              }
//                              // end % of region 4
//                      } else {
//                              // % region 3
//                              s = 0;
//                              if (e >= 0) {
//                                      t = 0;
//                                      sqrDistance = f;
//                              } else {
//                                      if (-e >= c) {
//                                              t = 1;
//                                              sqrDistance = c + 2 * e + f;
//                                      } else {
//                                              t = -e / c;
//                                              sqrDistance = e * t + f;
//                                      }
//                              }
//                      }
//                      // end of region 3
//              } else {
//                      if (t < 0) {
//                              //% region 5
//                              t = 0;
//                              if (d >= 0) {
//                                      s = 0;
//                                      sqrDistance = f;
//                              } else {
//                                      if (-d >= a) {
//                                              s = 1;
//                                              sqrDistance = a + 2 * d + f;
//                                      } else {
//                                              s = -d / a;
//                                              sqrDistance = d * s + f;
//                                      }
//                              }
//                      } else {
//                              // region 0
//                              invDet = 1 / det;
//                              s = s * invDet;
//                              t = t * invDet;
//                              sqrDistance = s * (a * s + b * t + 2 * d) + t * (b * s + c * t + 2 * e) + f;
//                      }
//              }
//      } else {
//              if (s < 0) {
//                      // % region 2
//                      tmp0 = b + d;
//                      tmp1 = c + e;
//                      if (tmp1 > tmp0) { //% minimum on edge s+t=1
//                              numer = tmp1 - tmp0;
//                              denom = a - 2 * b + c;
//                              if (numer >= denom) {
//                                      s = 1;
//                                      t = 0;
//                                      sqrDistance = a + 2 * d + f;
//                              } else {
//                                      s = numer / denom;
//                                      t = 1 - s;
//                                      sqrDistance = s * (a * s + b * t + 2 * d) + t * (b * s + c * t + 2 * e) + f;
//                              }
//                      } else
//                              // % minimum on edge s=0
//                              s = 0;
//                      if (tmp1 <= 0) {
//                              t = 1;
//                              sqrDistance = c + 2 * e + f;
//                      } else {
//                              if (e >= 0) {
//                                      t = 0;
//                                      sqrDistance = f;
//                              } else {
//                                      t = -e / c;
//                                      sqrDistance = e * t + f;
//                              }
//                      }
//              } //end % of region     2
//              else {
//                      if (t < 0) {
//                              // %region6
//                              tmp0 = b + e;
//                              tmp1 = a + d;
//                              if (tmp1 > tmp0) {
//                                      numer = tmp1 - tmp0;
//                                      denom = a - 2 * b + c;
//                                      if (numer >= denom) {
//                                              t = 1;
//                                              s = 0;
//                                              sqrDistance = c + 2 * e + f;
//                                      } else {
//                                              t = numer / denom;
//                                              s = 1 - t;
//                                              sqrDistance = s * (a * s + b * t + 2 * d) + t * (b * s + c * t + 2 * e) + f;
//                                      }
//                              } else {
//                                      t = 0;
//                                      if (tmp1 <= 0) {
//                                              s = 1;
//                                              sqrDistance = a + 2 * d + f;
//                                      } else {
//                                              if (d >= 0) {
//                                                      s = 0;
//                                                      sqrDistance = f;
//                                              } else {
//                                                      s = -d / a;
//                                                      sqrDistance = d * s + f;
//                                              }
//                                      }
//                              } // % end region 6
//                      } else {
//                              //% region 1
//                              numer = c + e - b - d;
//                              if (numer <= 0) {
//                                      s = 0;
//                                      t = 1;
//                                      sqrDistance = c + 2 * e + f;
//                              } else {
//                                      denom = a - 2 * b + c;
//                                      if (numer >= denom) {
//                                              s = 1;
//                                              t = 0;
//                                              sqrDistance = a + 2 * d + f;
//                                      } else {
//                                              s = numer / denom;
//                                              t = 1 - s;
//                                              sqrDistance = s * (a * s + b * t + 2 * d) + t * (b * s + c * t + 2 * e) + f;
//                                      }
//                              } //% end of region 1
//                      }
//              }
//      }
//
//      // % account for numerical round-off error
//      if (sqrDistance < 0) {
//              sqrDistance = 0;
//      }
//
//      dist = sqrt(sqrDistance);
//
//      // closest point
//      CP = B + s * E0 + t * E1;
//
//}
/*
 * % The algorithm is based on
 % "David Eberly, 'Distance Between Point and Triangle in 3D',
 % Geometric Tools, LLC, (1999)"
 % https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
 */

void PairTriSurf::PointTriangleDistance(const Vector3d& sourcePosition, const Vector3d& TRI0, const Vector3d& TRI1,
                const Vector3d& TRI2, Vector3d &CP, double &dist) {

        Vector3d edge0 = TRI1 - TRI0;
        Vector3d edge1 = TRI2 - TRI0;
        Vector3d v0 = TRI0 - sourcePosition;

        double a = edge0.dot(edge0);
        double b = edge0.dot(edge1);
        double c = edge1.dot(edge1);
        double d = edge0.dot(v0);
        double e = edge1.dot(v0);

        double det = a * c - b * b;
        double s = b * e - c * d;
        double t = b * d - a * e;

        if (s + t < det) {
                if (s < 0.f) {
                        if (t < 0.f) {
                                if (d < 0.f) {
                                        s = clamp(-d / a, 0.f, 1.f);
                                        t = 0.f;
                                } else {
                                        s = 0.f;
                                        t = clamp(-e / c, 0.f, 1.f);
                                }
                        } else {
                                s = 0.f;
                                t = clamp(-e / c, 0.f, 1.f);
                        }
                } else if (t < 0.f) {
                        s = clamp(-d / a, 0.f, 1.f);
                        t = 0.f;
                } else {
                        float invDet = 1.f / det;
                        s *= invDet;
                        t *= invDet;
                }
        } else {
                if (s < 0.f) {
                        float tmp0 = b + d;
                        float tmp1 = c + e;
                        if (tmp1 > tmp0) {
                                float numer = tmp1 - tmp0;
                                float denom = a - 2 * b + c;
                                s = clamp(numer / denom, 0.f, 1.f);
                                t = 1 - s;
                        } else {
                                t = clamp(-e / c, 0.f, 1.f);
                                s = 0.f;
                        }
                } else if (t < 0.f) {
                        if (a + d > b + e) {
                                float numer = c + e - b - d;
                                float denom = a - 2 * b + c;
                                s = clamp(numer / denom, 0.f, 1.f);
                                t = 1 - s;
                        } else {
                                s = clamp(-e / c, 0.f, 1.f);
                                t = 0.f;
                        }
                } else {
                        float numer = c + e - b - d;
                        float denom = a - 2 * b + c;
                        s = clamp(numer / denom, 0.f, 1.f);
                        t = 1.f - s;
                }
        }

        CP = TRI0 + s * edge0 + t * edge1;
        dist = (CP - sourcePosition).norm();

}

double PairTriSurf::clamp(const double a, const double min, const double max) {
        if (a < min) {
                return min;
        } else if (a > max) {
                return max;
        } else {
                return a;
        }
}

void *PairTriSurf::extract(const char *str, int &/*i*/) {
        //printf("in PairTriSurf::extract\n");
        if (strcmp(str, "smd/tri_surface/stable_time_increment_ptr") == 0) {
                return (void *) &stable_time_increment;
        }

        return nullptr;

}
