// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Dac Nguyen (ndactrung@gmail.com)
   Ref: Wang, Yu, Langston, Fraige, Particle shape effects in discrete
   element modelling of cohesive angular particles, Granular Matter 2011,
   13:1-12.
   Note: The current implementation has not taken into account
         the contact history for friction forces.
------------------------------------------------------------------------- */

#include "pair_body_rounded_polyhedron_omp.h"

#include "atom.h"
#include "atom_vec_body.h"
#include "body_rounded_polyhedron.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "suffix.h"

#include "omp_compat.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

namespace {
constexpr double EPSILON = 1.0e-3;
constexpr int MAX_CONTACTS = 32;

enum {EE_INVALID=0,EE_NONE,EE_INTERACT};
enum {EF_INVALID=0,EF_NONE,EF_PARALLEL,EF_SAME_SIDE_OF_FACE,
      EF_INTERSECT_INSIDE,EF_INTERSECT_OUTSIDE};
}

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolyhedronOMP::PairBodyRoundedPolyhedronOMP(LAMMPS *lmp) :
  PairBodyRoundedPolyhedron(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  // grow per-atom arrays if needed and initialize

  if (nall > nmax) {
    memory->destroy(dnum);
    memory->destroy(dfirst);
    memory->destroy(ednum);
    memory->destroy(edfirst);
    memory->destroy(facnum);
    memory->destroy(facfirst);
    memory->destroy(enclosing_radius);
    memory->destroy(rounded_radius);
    nmax = nall;
    memory->create(dnum, nmax, "pair:dnum");
    memory->create(dfirst, nmax, "pair:dfirst");
    memory->create(ednum, nmax, "pair:ednum");
    memory->create(edfirst, nmax, "pair:edfirst");
    memory->create(facnum, nmax, "pair:facnum");
    memory->create(facfirst, nmax, "pair:facfirst");
    memory->create(enclosing_radius, nmax, "pair:enclosing_radius");
    memory->create(rounded_radius, nmax, "pair:rounded_radius");
  }

  ndiscrete = nedge = nface = 0;
  for (int i = 0; i < nall; i++) dnum[i] = ednum[i] = facnum[i] = 0;

  // pre-initialize body2space for all body atoms serially before the parallel region
  // body2space modifies shared arrays (discrete, edge, face, dnum, dfirst, ...)
  // and must NOT be called from within the parallel region

  int *body = atom->body;
  for (int i = 0; i < nall; i++) {
    if (body[i] >= 0) body2space(i);
  }

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    if (evflag) {
      if (eflag) {
        eval<1,1>(ifrom, ito, thr);
      } else {
        eval<1,0>(ifrom, ito, thr);
      }
    } else {
      eval<0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG>
void PairBodyRoundedPolyhedronOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  int i, j, ii, jj, jnum, itype, jtype;
  int ni, nj, npi, npj, ifirst, jfirst;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl;
  double facc[3];
  double rsq, eradi, eradj, r;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;

  double **x = atom->x;
  double **v = atom->v;
  auto *const f = (dbl3_t *) thr->get_f()[0];
  auto *const torque = (dbl3_t *) thr->get_torque()[0];
  double **angmom = atom->angmom;
  int *const body = atom->body;
  int *const type = atom->type;
  const int nlocal = atom->nlocal;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // thread-local visited-vertex flag array replaces discrete[*][6]

  double *vflag_thr = nullptr;
  if (dmax > 0) memory->create(vflag_thr, dmax, "pair:vflag_thr");

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (body[i] >= 0) {
      npi = dnum[i];
      ifirst = dfirst[i];
      eradi = enclosing_radius[i];
    }

    for (jj = 0; jj < jnum; ++jj) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      evdwl = 0.0;
      facc[0] = facc[1] = facc[2] = 0;

      if (body[i] < 0 || body[j] < 0) continue;

      npj = dnum[j];
      jfirst = dfirst[j];
      eradj = enclosing_radius[j];

      r = sqrt(rsq);
      if (r > eradi + eradj + cut_inner) continue;

      // sphere-sphere interaction

      if (npi == 1 && npj == 1) {
        sphere_against_sphere_thr(i, j, itype, jtype, delx, dely, delz, rsq,
                                  v, f, EVFLAG, thr);
        continue;
      }

      // reset thread-local visited-vertex flags for both bodies

      if (vflag_thr) {
        for (ni = 0; ni < npi; ni++) vflag_thr[ifirst+ni] = 0.0;
        for (nj = 0; nj < npj; nj++) vflag_thr[jfirst+nj] = 0.0;
      }

      // one body is a sphere

      if (npj == 1) {
        sphere_against_face_thr(i, j, itype, jtype, x, v, f, torque,
                                angmom, EVFLAG, thr);
        sphere_against_edge_thr(i, j, itype, jtype, x, v, f, torque,
                                angmom, vflag_thr, EVFLAG, thr);
        continue;
      } else if (npi == 1) {
        sphere_against_face_thr(j, i, jtype, itype, x, v, f, torque,
                                angmom, EVFLAG, thr);
        sphere_against_edge_thr(j, i, jtype, itype, x, v, f, torque,
                                angmom, vflag_thr, EVFLAG, thr);
        continue;
      }

      // polyhedron-polyhedron interaction

      int num_contacts = 0;
      Contact contact_list[MAX_CONTACTS];

      edge_against_face_thr(i, j, itype, jtype, x, contact_list, num_contacts,
                            evdwl, facc, f, torque, v, angmom, vflag_thr);
      edge_against_face_thr(j, i, jtype, itype, x, contact_list, num_contacts,
                            evdwl, facc, f, torque, v, angmom, vflag_thr);
      edge_against_edge_thr(i, j, itype, jtype, x, contact_list, num_contacts,
                            evdwl, facc, f, torque, v, angmom);

      if (num_contacts > 0)
        rescale_cohesive_forces_thr(x, f, torque, contact_list, num_contacts,
                                    itype, jtype, facc);

      if (EVFLAG)
        ev_tally_xyz_thr(this, i, j, nlocal, force->newton_pair,
                         evdwl, 0.0, facc[0], facc[1], facc[2],
                         delx, dely, delz, thr);

    } // end for jj
  } // end for ii

  if (vflag_thr) memory->destroy(vflag_thr);
}

/* ----------------------------------------------------------------------
   Sphere-sphere interaction using per-thread force array
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronOMP::sphere_against_sphere_thr(
    int ibody, int jbody, int itype, int jtype,
    double delx, double dely, double delz, double rsq,
    double **v, dbl3_t *f, int evflag, ThrData *thr)
{
  double rradi, rradj, contact_dist;
  double vr1, vr2, vr3, vnnr, vn1, vn2, vn3, vt1, vt2, vt3;
  double rij, rsqinv, R, fx, fy, fz, fn[3], ft[3], fpair, energy;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  rradi = rounded_radius[ibody];
  rradj = rounded_radius[jbody];
  contact_dist = rradi + rradj;

  rij = sqrt(rsq);
  R = rij - contact_dist;

  energy = 0;
  kernel_force(R, itype, jtype, energy, fpair);

  fx = delx*fpair/rij;
  fy = dely*fpair/rij;
  fz = delz*fpair/rij;

  if (R <= 0) {
    rsqinv = 1.0/rsq;
    vr1 = v[ibody][0] - v[jbody][0];
    vr2 = v[ibody][1] - v[jbody][1];
    vr3 = v[ibody][2] - v[jbody][2];

    vnnr = vr1*delx + vr2*dely + vr3*delz;
    vn1 = delx*vnnr * rsqinv;
    vn2 = dely*vnnr * rsqinv;
    vn3 = delz*vnnr * rsqinv;

    vt1 = vr1 - vn1;
    vt2 = vr2 - vn2;
    vt3 = vr3 - vn3;

    fn[0] = -c_n * vn1;
    fn[1] = -c_n * vn2;
    fn[2] = -c_n * vn3;

    ft[0] = -c_t * vt1;
    ft[1] = -c_t * vt2;
    ft[2] = -c_t * vt3;

    fx += fn[0] + ft[0];
    fy += fn[1] + ft[1];
    fz += fn[2] + ft[2];
  }

  f[ibody].x += fx;
  f[ibody].y += fy;
  f[ibody].z += fz;

  if (newton_pair || jbody < nlocal) {
    f[jbody].x -= fx;
    f[jbody].y -= fy;
    f[jbody].z -= fz;
  }

  if (evflag)
    ev_tally_xyz_thr(this, ibody, jbody, nlocal, newton_pair,
                     energy, 0.0, fx, fy, fz, delx, dely, delz, thr);
}

/* ----------------------------------------------------------------------
   Sphere-edge interaction using per-thread force/torque and vflag_thr
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronOMP::sphere_against_edge_thr(
    int ibody, int jbody, int itype, int jtype,
    double **x, double **v, dbl3_t *f, dbl3_t *torque,
    double **angmom, double *vflag_thr, int evflag, ThrData *thr)
{
  int ni, nei, ifirst, iefirst, npi1, npi2, ibonus;
  double xi1[3], xi2[3], vti[3], h[3], fn[3], ft[3], d, t;
  double delx, dely, delz, rsq, rij, rsqinv, R, fx, fy, fz, fpair, energy;
  double rradi, rradj, contact_dist;
  double vr1, vr2, vr3, vnnr, vn1, vn2, vn3, vt1, vt2, vt3;
  double *quat, *inertia;
  AtomVecBody::Bonus *bonus;

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  ifirst = dfirst[ibody];
  iefirst = edfirst[ibody];
  nei = ednum[ibody];

  rradi = rounded_radius[ibody];
  rradj = rounded_radius[jbody];
  contact_dist = rradi + rradj;

  for (ni = 0; ni < nei; ni++) {
    npi1 = static_cast<int>(edge[iefirst+ni][0]);
    npi2 = static_cast<int>(edge[iefirst+ni][1]);

    xi1[0] = x[ibody][0] + discrete[ifirst+npi1][0];
    xi1[1] = x[ibody][1] + discrete[ifirst+npi1][1];
    xi1[2] = x[ibody][2] + discrete[ifirst+npi1][2];

    xi2[0] = x[ibody][0] + discrete[ifirst+npi2][0];
    xi2[1] = x[ibody][1] + discrete[ifirst+npi2][1];
    xi2[2] = x[ibody][2] + discrete[ifirst+npi2][2];

    project_pt_line(x[jbody], xi1, xi2, h, d, t);

    if (d > contact_dist + cut_inner) continue;
    if (t < 0 || t > 1) continue;

    if (fabs(t) < EPSILON) {
      if (vflag_thr && static_cast<int>(vflag_thr[ifirst+npi1]) == 1)
        continue;
      else {
        h[0] = xi1[0]; h[1] = xi1[1]; h[2] = xi1[2];
        if (vflag_thr) vflag_thr[ifirst+npi1] = 1.0;
      }
    }

    if (fabs(t-1) < EPSILON) {
      if (vflag_thr && static_cast<int>(vflag_thr[ifirst+npi2]) == 1)
        continue;
      else {
        h[0] = xi2[0]; h[1] = xi2[1]; h[2] = xi2[2];
        if (vflag_thr) vflag_thr[ifirst+npi2] = 1.0;
      }
    }

    delx = h[0] - x[jbody][0];
    dely = h[1] - x[jbody][1];
    delz = h[2] - x[jbody][2];
    rsq = delx*delx + dely*dely + delz*delz;
    rsqinv = (rsq == 0.0) ? 0.0 : 1.0/rsq;
    rij = sqrt(rsq);
    R = rij - contact_dist;

    energy = 0;
    kernel_force(R, itype, jtype, energy, fpair);

    fx = delx*fpair/rij;
    fy = dely*fpair/rij;
    fz = delz*fpair/rij;

    if (R <= 0) {
      ibonus = atom->body[ibody];
      bonus = &avec->bonus[ibonus];
      quat = bonus->quat;
      inertia = bonus->inertia;
      total_velocity(h, x[ibody], v[ibody], angmom[ibody], inertia, quat, vti);

      vr1 = vti[0] - v[jbody][0];
      vr2 = vti[1] - v[jbody][1];
      vr3 = vti[2] - v[jbody][2];

      vnnr = vr1*delx + vr2*dely + vr3*delz;
      vn1 = delx*vnnr * rsqinv;
      vn2 = dely*vnnr * rsqinv;
      vn3 = delz*vnnr * rsqinv;

      vt1 = vr1 - vn1;
      vt2 = vr2 - vn2;
      vt3 = vr3 - vn3;

      fn[0] = -c_n * vn1;
      fn[1] = -c_n * vn2;
      fn[2] = -c_n * vn3;

      ft[0] = -c_t * vt1;
      ft[1] = -c_t * vt2;
      ft[2] = -c_t * vt3;

      fx += fn[0] + ft[0];
      fy += fn[1] + ft[1];
      fz += fn[2] + ft[2];
    }

    f[ibody].x += fx;
    f[ibody].y += fy;
    f[ibody].z += fz;
    sum_torque(x[ibody], h, fx, fy, fz, &torque[ibody].x);

    if (newton_pair || jbody < nlocal) {
      f[jbody].x -= fx;
      f[jbody].y -= fy;
      f[jbody].z -= fz;
    }

    if (evflag)
      ev_tally_xyz_thr(this, ibody, jbody, nlocal, newton_pair,
                       energy, 0.0, fx, fy, fz, delx, dely, delz, thr);
  }
}

/* ----------------------------------------------------------------------
   Sphere-face interaction using per-thread force/torque
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronOMP::sphere_against_face_thr(
    int ibody, int jbody, int itype, int jtype,
    double **x, double **v, dbl3_t *f, dbl3_t *torque,
    double **angmom, int evflag, ThrData *thr)
{
  int ni, nfi, inside, ifirst, iffirst, npi1, npi2, npi3, ibonus, tmp;
  double xi1[3], xi2[3], xi3[3], ui[3], vi[3], vti[3], n[3], h[3], fn[3], ft[3], d;
  double delx, dely, delz, rsq, rij, rsqinv, R, fx, fy, fz, fpair, energy;
  double rradi, rradj, contact_dist;
  double vr1, vr2, vr3, vnnr, vn1, vn2, vn3, vt1, vt2, vt3;
  double *quat, *inertia;
  AtomVecBody::Bonus *bonus;

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  ifirst = dfirst[ibody];
  iffirst = facfirst[ibody];
  nfi = facnum[ibody];

  rradi = rounded_radius[ibody];
  rradj = rounded_radius[jbody];
  contact_dist = rradi + rradj;

  for (ni = 0; ni < nfi; ni++) {
    npi1 = static_cast<int>(face[iffirst+ni][0]);
    npi2 = static_cast<int>(face[iffirst+ni][1]);
    npi3 = static_cast<int>(face[iffirst+ni][2]);

    xi1[0] = x[ibody][0] + discrete[ifirst+npi1][0];
    xi1[1] = x[ibody][1] + discrete[ifirst+npi1][1];
    xi1[2] = x[ibody][2] + discrete[ifirst+npi1][2];

    xi2[0] = x[ibody][0] + discrete[ifirst+npi2][0];
    xi2[1] = x[ibody][1] + discrete[ifirst+npi2][1];
    xi2[2] = x[ibody][2] + discrete[ifirst+npi2][2];

    xi3[0] = x[ibody][0] + discrete[ifirst+npi3][0];
    xi3[1] = x[ibody][1] + discrete[ifirst+npi3][1];
    xi3[2] = x[ibody][2] + discrete[ifirst+npi3][2];

    MathExtra::sub3(xi2, xi1, ui);
    MathExtra::sub3(xi3, xi1, vi);
    MathExtra::cross3(ui, vi, n);
    MathExtra::norm3(n);

    if (opposite_sides(n, xi1, x[ibody], x[jbody]) == 0) continue;

    project_pt_plane(x[jbody], xi1, xi2, xi3, h, d, inside);

    inside_polygon(ibody, ni, x[ibody], h, nullptr, inside, tmp);
    if (inside == 0) continue;

    delx = h[0] - x[jbody][0];
    dely = h[1] - x[jbody][1];
    delz = h[2] - x[jbody][2];
    rsq = delx*delx + dely*dely + delz*delz;
    rij = sqrt(rsq);
    R = rij - contact_dist;

    energy = 0;
    kernel_force(R, itype, jtype, energy, fpair);

    fx = delx*fpair/rij;
    fy = dely*fpair/rij;
    fz = delz*fpair/rij;

    if (R <= 0) {
      rsqinv = 1.0/rsq;
      ibonus = atom->body[ibody];
      bonus = &avec->bonus[ibonus];
      quat = bonus->quat;
      inertia = bonus->inertia;
      total_velocity(h, x[ibody], v[ibody], angmom[ibody], inertia, quat, vti);

      vr1 = vti[0] - v[jbody][0];
      vr2 = vti[1] - v[jbody][1];
      vr3 = vti[2] - v[jbody][2];

      vnnr = vr1*delx + vr2*dely + vr3*delz;
      vn1 = delx*vnnr * rsqinv;
      vn2 = dely*vnnr * rsqinv;
      vn3 = delz*vnnr * rsqinv;

      vt1 = vr1 - vn1;
      vt2 = vr2 - vn2;
      vt3 = vr3 - vn3;

      fn[0] = -c_n * vn1;
      fn[1] = -c_n * vn2;
      fn[2] = -c_n * vn3;

      ft[0] = -c_t * vt1;
      ft[1] = -c_t * vt2;
      ft[2] = -c_t * vt3;

      fx += fn[0] + ft[0];
      fy += fn[1] + ft[1];
      fz += fn[2] + ft[2];
    }

    f[ibody].x += fx;
    f[ibody].y += fy;
    f[ibody].z += fz;
    sum_torque(x[ibody], h, fx, fy, fz, &torque[ibody].x);

    if (newton_pair || jbody < nlocal) {
      f[jbody].x -= fx;
      f[jbody].y -= fy;
      f[jbody].z -= fz;
    }

    if (evflag)
      ev_tally_xyz_thr(this, ibody, jbody, nlocal, newton_pair,
                       energy, 0.0, fx, fy, fz, delx, dely, delz, thr);
  }
}

/* ----------------------------------------------------------------------
   Edge-edge interactions using per-thread arrays
------------------------------------------------------------------------- */

int PairBodyRoundedPolyhedronOMP::edge_against_edge_thr(
    int ibody, int jbody, int itype, int jtype,
    double **x, Contact *contact_list, int &num_contacts,
    double &evdwl, double *facc,
    dbl3_t *f, dbl3_t *torque, double **v, double **angmom)
{
  int ni, nei, nj, nej, interact;
  double rradi, rradj, energy;

  nei = ednum[ibody];
  rradi = rounded_radius[ibody];
  nej = ednum[jbody];
  rradj = rounded_radius[jbody];

  energy = 0;
  interact = EE_NONE;

  for (ni = 0; ni < nei; ni++) {
    for (nj = 0; nj < nej; nj++) {
      interact = interaction_edge_to_edge_thr(
          ibody, ni, x[ibody], rradi,
          jbody, nj, x[jbody], rradj,
          itype, jtype, cut_inner,
          contact_list, num_contacts, energy, facc,
          f, torque, x, v, angmom);
    }
  }

  evdwl += energy;
  return interact;
}

/* ----------------------------------------------------------------------
   Edge-face interactions using per-thread arrays
------------------------------------------------------------------------- */

int PairBodyRoundedPolyhedronOMP::edge_against_face_thr(
    int ibody, int jbody, int itype, int jtype,
    double **x, Contact *contact_list, int &num_contacts,
    double &evdwl, double *facc,
    dbl3_t *f, dbl3_t *torque, double **v, double **angmom,
    double *vflag_thr)
{
  int ni, nei, nj, nfj, interact;
  double rradi, rradj, energy;

  nei = ednum[ibody];
  rradi = rounded_radius[ibody];
  nfj = facnum[jbody];
  rradj = rounded_radius[jbody];

  energy = 0;
  interact = EF_NONE;

  for (ni = 0; ni < nei; ni++) {
    for (nj = 0; nj < nfj; nj++) {
      interact = interaction_face_to_edge_thr(
          jbody, nj, x[jbody], rradj,
          ibody, ni, x[ibody], rradi,
          itype, jtype, cut_inner,
          contact_list, num_contacts, energy, facc,
          f, torque, x, v, angmom, vflag_thr);
    }
  }

  evdwl += energy;
  return interact;
}

/* ----------------------------------------------------------------------
   Interaction between an edge from jbody and a face from ibody
   Uses vflag_thr instead of discrete[*][6] for vertex visited flag
------------------------------------------------------------------------- */

int PairBodyRoundedPolyhedronOMP::interaction_face_to_edge_thr(
    int ibody, int face_index, double *xmi, double rounded_radius_i,
    int jbody, int edge_index, double *xmj, double rounded_radius_j,
    int itype, int jtype, double cut_inner,
    Contact *contact_list, int &num_contacts, double &energy, double *facc,
    dbl3_t *f, dbl3_t *torque, double **x, double **v, double **angmom,
    double *vflag_thr)
{
  if (face_index >= facnum[ibody]) return EF_INVALID;

  int ifirst, iffirst, jfirst, npi1, npi2, npi3;
  int jefirst, npj1, npj2;
  double xi1[3], xi2[3], xi3[3], xpj1[3], xpj2[3], ui[3], vi[3], n[3];

  ifirst = dfirst[ibody];
  iffirst = facfirst[ibody];
  npi1 = static_cast<int>(face[iffirst+face_index][0]);
  npi2 = static_cast<int>(face[iffirst+face_index][1]);
  npi3 = static_cast<int>(face[iffirst+face_index][2]);

  xi1[0] = xmi[0] + discrete[ifirst+npi1][0];
  xi1[1] = xmi[1] + discrete[ifirst+npi1][1];
  xi1[2] = xmi[2] + discrete[ifirst+npi1][2];

  xi2[0] = xmi[0] + discrete[ifirst+npi2][0];
  xi2[1] = xmi[1] + discrete[ifirst+npi2][1];
  xi2[2] = xmi[2] + discrete[ifirst+npi2][2];

  xi3[0] = xmi[0] + discrete[ifirst+npi3][0];
  xi3[1] = xmi[1] + discrete[ifirst+npi3][1];
  xi3[2] = xmi[2] + discrete[ifirst+npi3][2];

  MathExtra::sub3(xi2, xi1, ui);
  MathExtra::sub3(xi3, xi1, vi);
  MathExtra::cross3(ui, vi, n);
  MathExtra::norm3(n);

  double xc[3], dot, ans[3];
  xc[0] = (xi1[0] + xi2[0] + xi3[0]) / 3.0;
  xc[1] = (xi1[1] + xi2[1] + xi3[1]) / 3.0;
  xc[2] = (xi1[2] + xi2[2] + xi3[2]) / 3.0;
  MathExtra::sub3(xc, xmi, ans);
  dot = MathExtra::dot3(ans, n);
  if (dot < 0) MathExtra::negate3(n);

  jfirst = dfirst[jbody];
  jefirst = edfirst[jbody];
  npj1 = static_cast<int>(edge[jefirst+edge_index][0]);
  npj2 = static_cast<int>(edge[jefirst+edge_index][1]);

  xpj1[0] = xmj[0] + discrete[jfirst+npj1][0];
  xpj1[1] = xmj[1] + discrete[jfirst+npj1][1];
  xpj1[2] = xmj[2] + discrete[jfirst+npj1][2];

  xpj2[0] = xmj[0] + discrete[jfirst+npj2][0];
  xpj2[1] = xmj[1] + discrete[jfirst+npj2][1];
  xpj2[2] = xmj[2] + discrete[jfirst+npj2][2];

  if (opposite_sides(n, xi1, xmi, xpj1) == 0 &&
      opposite_sides(n, xi1, xmi, xpj2) == 0)
    return EF_NONE;

  double hi1[3], hi2[3], d1, d2, contact_dist;
  int inside1 = 0;
  int inside2 = 0;

  int interact = edge_face_intersect(xi1, xi2, xi3, xpj1, xpj2,
                                     hi1, hi2, d1, d2, inside1, inside2);

  inside_polygon(ibody, face_index, xmi, hi1, hi2, inside1, inside2);

  contact_dist = rounded_radius_i + rounded_radius_j;

  if (interact == EF_SAME_SIDE_OF_FACE || interact == EF_PARALLEL) {

    if (d1 > contact_dist + cut_inner && d2 > contact_dist + cut_inner)
      return EF_NONE;

    int num_outside = 0;
    int jflag = 1;

    if (d1 <= contact_dist + cut_inner) {
      if (inside1) {
        if (!vflag_thr || static_cast<int>(vflag_thr[jfirst+npj1]) == 0) {
          pair_force_and_torque_thr(jbody, ibody, xpj1, hi1, d1, contact_dist,
                                    jtype, itype, x, v, f, torque, angmom,
                                    jflag, energy, facc);
          if (d1 <= contact_dist) {
            contact_list[num_contacts].ibody = ibody;
            contact_list[num_contacts].jbody = jbody;
            contact_list[num_contacts].xi[0] = hi1[0];
            contact_list[num_contacts].xi[1] = hi1[1];
            contact_list[num_contacts].xi[2] = hi1[2];
            contact_list[num_contacts].xj[0] = xpj1[0];
            contact_list[num_contacts].xj[1] = xpj1[1];
            contact_list[num_contacts].xj[2] = xpj1[2];
            contact_list[num_contacts].type = 0;
            contact_list[num_contacts].separation = d1 - contact_dist;
            contact_list[num_contacts].unique = 1;
            num_contacts++;
          }
          if (vflag_thr) vflag_thr[jfirst+npj1] = 1.0;
        }
      } else {
        num_outside++;
      }
    }

    if (d2 <= contact_dist + cut_inner) {
      if (inside2) {
        if (!vflag_thr || static_cast<int>(vflag_thr[jfirst+npj2]) == 0) {
          pair_force_and_torque_thr(jbody, ibody, xpj2, hi2, d2, contact_dist,
                                    jtype, itype, x, v, f, torque, angmom,
                                    jflag, energy, facc);
          if (d2 <= contact_dist) {
            contact_list[num_contacts].ibody = ibody;
            contact_list[num_contacts].jbody = jbody;
            contact_list[num_contacts].xi[0] = hi2[0];
            contact_list[num_contacts].xi[1] = hi2[1];
            contact_list[num_contacts].xi[2] = hi2[2];
            contact_list[num_contacts].xj[0] = xpj2[0];
            contact_list[num_contacts].xj[1] = xpj2[1];
            contact_list[num_contacts].xj[2] = xpj2[2];
            contact_list[num_contacts].type = 0;
            contact_list[num_contacts].separation = d2 - contact_dist;
            contact_list[num_contacts].unique = 1;
            num_contacts++;
          }
          if (vflag_thr) vflag_thr[jfirst+npj2] = 1.0;
        }
      } else {
        num_outside++;
      }
    }

    if (num_outside == 2) interact = EF_INTERSECT_OUTSIDE;

  } else if (interact == EF_INTERSECT_OUTSIDE) {
    // edge endpoints both projected outside the face triangle;
    // edge-vs-face-edge interactions are handled by edge_against_edge

  } else if (interact == EF_INTERSECT_INSIDE) {
    int jflag = 1;
    if (d1 < d2)
      pair_force_and_torque_thr(jbody, ibody, xpj1, hi1, d1, contact_dist,
                                jtype, itype, x, v, f, torque, angmom,
                                jflag, energy, facc);
    else
      pair_force_and_torque_thr(jbody, ibody, xpj2, hi2, d2, contact_dist,
                                jtype, itype, x, v, f, torque, angmom,
                                jflag, energy, facc);
  }

  return interact;
}

/* ----------------------------------------------------------------------
   Interaction between two edges using per-thread arrays
------------------------------------------------------------------------- */

int PairBodyRoundedPolyhedronOMP::interaction_edge_to_edge_thr(
    int ibody, int edge_index_i, double *xmi, double rounded_radius_i,
    int jbody, int edge_index_j, double *xmj, double rounded_radius_j,
    int itype, int jtype, double cut_inner,
    Contact *contact_list, int &num_contacts, double &energy, double *facc,
    dbl3_t *f, dbl3_t *torque, double **x, double **v, double **angmom)
{
  int ifirst, iefirst, jfirst, jefirst, npi1, npi2, npj1, npj2, interact;
  double xi1[3], xi2[3], xpj1[3], xpj2[3];
  double r, t1, t2, h1[3], h2[3];
  double contact_dist;

  ifirst = dfirst[ibody];
  iefirst = edfirst[ibody];
  npi1 = static_cast<int>(edge[iefirst+edge_index_i][0]);
  npi2 = static_cast<int>(edge[iefirst+edge_index_i][1]);

  xi1[0] = xmi[0] + discrete[ifirst+npi1][0];
  xi1[1] = xmi[1] + discrete[ifirst+npi1][1];
  xi1[2] = xmi[2] + discrete[ifirst+npi1][2];

  xi2[0] = xmi[0] + discrete[ifirst+npi2][0];
  xi2[1] = xmi[1] + discrete[ifirst+npi2][1];
  xi2[2] = xmi[2] + discrete[ifirst+npi2][2];

  jfirst = dfirst[jbody];
  jefirst = edfirst[jbody];
  npj1 = static_cast<int>(edge[jefirst+edge_index_j][0]);
  npj2 = static_cast<int>(edge[jefirst+edge_index_j][1]);

  xpj1[0] = xmj[0] + discrete[jfirst+npj1][0];
  xpj1[1] = xmj[1] + discrete[jfirst+npj1][1];
  xpj1[2] = xmj[2] + discrete[jfirst+npj1][2];

  xpj2[0] = xmj[0] + discrete[jfirst+npj2][0];
  xpj2[1] = xmj[1] + discrete[jfirst+npj2][1];
  xpj2[2] = xmj[2] + discrete[jfirst+npj2][2];

  contact_dist = rounded_radius_i + rounded_radius_j;

  int jflag = 1;
  distance_bt_edges(xpj1, xpj2, xi1, xi2, h1, h2, t1, t2, r);

  interact = EE_NONE;

  double rmin = MIN(rounded_radius_i, rounded_radius_j);
  if (r < EPSILON*rmin) return interact;

  if (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1 &&
      r < contact_dist + cut_inner) {
    pair_force_and_torque_thr(jbody, ibody, h1, h2, r, contact_dist,
                              jtype, itype, x, v, f, torque, angmom,
                              jflag, energy, facc);
    interact = EE_INTERACT;
    if (r <= contact_dist) {
      contact_list[num_contacts].ibody = ibody;
      contact_list[num_contacts].jbody = jbody;
      contact_list[num_contacts].xi[0] = h2[0];
      contact_list[num_contacts].xi[1] = h2[1];
      contact_list[num_contacts].xi[2] = h2[2];
      contact_list[num_contacts].xj[0] = h1[0];
      contact_list[num_contacts].xj[1] = h1[1];
      contact_list[num_contacts].xj[2] = h1[2];
      contact_list[num_contacts].type = 1;
      contact_list[num_contacts].separation = r - contact_dist;
      contact_list[num_contacts].unique = 1;
      num_contacts++;
    }
  }
  return interact;
}

/* ----------------------------------------------------------------------
   Compute force and torque between two bodies given a pair of contact
   points; uses per-thread f and torque arrays
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronOMP::pair_force_and_torque_thr(
    int ibody, int jbody, double *pi, double *pj, double r,
    double contact_dist, int itype, int jtype,
    double **x, double **v, dbl3_t *f, dbl3_t *torque, double **angmom,
    int jflag, double &energy, double *facc)
{
  double delx, dely, delz, R, fx, fy, fz, fpair;

  delx = pi[0] - pj[0];
  dely = pi[1] - pj[1];
  delz = pi[2] - pj[2];
  R = r - contact_dist;

  kernel_force(R, itype, jtype, energy, fpair);

  fx = delx*fpair/r;
  fy = dely*fpair/r;
  fz = delz*fpair/r;

  if (R <= 0) {
    contact_forces_thr(ibody, jbody, pi, pj, delx, dely, delz, fx, fy, fz,
                       x, v, angmom, f, torque, facc);
  } else {
    f[ibody].x += fx;
    f[ibody].y += fy;
    f[ibody].z += fz;
    sum_torque(x[ibody], pi, fx, fy, fz, &torque[ibody].x);

    facc[0] += fx; facc[1] += fy; facc[2] += fz;

    if (jflag) {
      f[jbody].x -= fx;
      f[jbody].y -= fy;
      f[jbody].z -= fz;
      sum_torque(x[jbody], pj, -fx, -fy, -fz, &torque[jbody].x);
    }
  }
}

/* ----------------------------------------------------------------------
   Compute contact forces; uses per-thread f and torque arrays
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronOMP::contact_forces_thr(
    int ibody, int jbody, double *xi, double *xj,
    double delx, double dely, double delz,
    double fx, double fy, double fz,
    double **x, double **v, double **angmom,
    dbl3_t *f, dbl3_t *torque, double *facc)
{
  int ibonus, jbonus;
  double rsq, rsqinv, vr1, vr2, vr3, vnnr, vn1, vn2, vn3, vt1, vt2, vt3;
  double fn[3], ft[3], vi[3], vj[3];
  double *quat, *inertia;
  AtomVecBody::Bonus *bonus;

  ibonus = atom->body[ibody];
  bonus = &avec->bonus[ibonus];
  quat = bonus->quat;
  inertia = bonus->inertia;
  total_velocity(xi, x[ibody], v[ibody], angmom[ibody], inertia, quat, vi);

  jbonus = atom->body[jbody];
  bonus = &avec->bonus[jbonus];
  quat = bonus->quat;
  inertia = bonus->inertia;
  total_velocity(xj, x[jbody], v[jbody], angmom[jbody], inertia, quat, vj);

  rsq = delx*delx + dely*dely + delz*delz;
  rsqinv = 1.0/rsq;

  vr1 = vi[0] - vj[0];
  vr2 = vi[1] - vj[1];
  vr3 = vi[2] - vj[2];

  vnnr = vr1*delx + vr2*dely + vr3*delz;
  vn1 = delx*vnnr * rsqinv;
  vn2 = dely*vnnr * rsqinv;
  vn3 = delz*vnnr * rsqinv;

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  fn[0] = -c_n * vn1;
  fn[1] = -c_n * vn2;
  fn[2] = -c_n * vn3;

  ft[0] = -c_t * vt1;
  ft[1] = -c_t * vt2;
  ft[2] = -c_t * vt3;

  fx = fn[0] + ft[0] + mu * fx;
  fy = fn[1] + ft[1] + mu * fy;
  fz = fn[2] + ft[2] + mu * fz;

  f[ibody].x += fx;
  f[ibody].y += fy;
  f[ibody].z += fz;
  sum_torque(x[ibody], xi, fx, fy, fz, &torque[ibody].x);

  f[jbody].x -= fx;
  f[jbody].y -= fy;
  f[jbody].z -= fz;
  sum_torque(x[jbody], xj, -fx, -fy, -fz, &torque[jbody].x);

  facc[0] += fx; facc[1] += fy; facc[2] += fz;
}

/* ----------------------------------------------------------------------
   Rescale cohesive forces for all contacts; uses per-thread f and torque
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronOMP::rescale_cohesive_forces_thr(
    double **x, dbl3_t *f, dbl3_t *torque,
    Contact *contact_list, int &num_contacts,
    int itype, int jtype, double *facc)
{
  int m, ibody, jbody;
  double delx, dely, delz, fx, fy, fz, R, fpair, r, contact_area;

  int num_unique_contacts = 0;
  if (num_contacts == 1) {
    num_unique_contacts = 1;
    contact_area = 0;
  } else if (num_contacts == 2) {
    num_unique_contacts = 2;
    contact_area = num_contacts * A_ua;
  } else {
    find_unique_contacts(contact_list, num_contacts);

    double xc[3], dx, dy, dz;
    xc[0] = xc[1] = xc[2] = 0;
    num_unique_contacts = 0;
    for (m = 0; m < num_contacts; m++) {
      if (contact_list[m].unique == 0) continue;
      xc[0] += contact_list[m].xi[0];
      xc[1] += contact_list[m].xi[1];
      xc[2] += contact_list[m].xi[2];
      num_unique_contacts++;
    }

    const double dble_unique =
      (num_unique_contacts > 0) ? (double) num_unique_contacts : 1.0;
    xc[0] /= dble_unique;
    xc[1] /= dble_unique;
    xc[2] /= dble_unique;

    contact_area = 0.0;
    for (m = 0; m < num_contacts; m++) {
      if (contact_list[m].unique == 0) continue;
      dx = contact_list[m].xi[0] - xc[0];
      dy = contact_list[m].xi[1] - xc[1];
      dz = contact_list[m].xi[2] - xc[2];
      contact_area += (dx*dx + dy*dy + dz*dz);
    }
    contact_area *= (MY_PI / dble_unique);
  }

  double j_a = contact_area / (num_unique_contacts * A_ua);
  if (j_a < 1.0) j_a = 1.0;

  for (m = 0; m < num_contacts; m++) {
    if (contact_list[m].unique == 0) continue;

    ibody = contact_list[m].ibody;
    jbody = contact_list[m].jbody;

    delx = contact_list[m].xi[0] - contact_list[m].xj[0];
    dely = contact_list[m].xi[1] - contact_list[m].xj[1];
    delz = contact_list[m].xi[2] - contact_list[m].xj[2];
    r = sqrt(delx*delx + dely*dely + delz*delz);
    R = contact_list[m].separation;

    double energy = 0;
    kernel_force(R, itype, jtype, energy, fpair);

    fpair *= j_a;
    fx = delx*fpair/r;
    fy = dely*fpair/r;
    fz = delz*fpair/r;

    f[ibody].x += fx;
    f[ibody].y += fy;
    f[ibody].z += fz;
    sum_torque(x[ibody], contact_list[m].xi, fx, fy, fz, &torque[ibody].x);

    f[jbody].x -= fx;
    f[jbody].y -= fy;
    f[jbody].z -= fz;
    sum_torque(x[jbody], contact_list[m].xj, -fx, -fy, -fz, &torque[jbody].x);

    facc[0] += fx; facc[1] += fy; facc[2] += fz;
  }
}

/* ---------------------------------------------------------------------- */

double PairBodyRoundedPolyhedronOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairBodyRoundedPolyhedron::memory_usage();
  return bytes;
}
