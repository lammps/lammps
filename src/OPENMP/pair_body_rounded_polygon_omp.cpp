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
------------------------------------------------------------------------- */

#include "pair_body_rounded_polygon_omp.h"

#include "atom.h"
#include "atom_vec_body.h"
#include "body_rounded_polygon.h"
#include "comm.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "suffix.h"

#include "omp_compat.h"

#include <cmath>

using namespace LAMMPS_NS;

namespace {
constexpr double EPSILON = 1.0e-3;
constexpr int MAX_CONTACTS = 4;
constexpr int EFF_CONTACTS = 2;

enum { INVALID=0, NONE=1, VERTEXI=2, VERTEXJ=3, EDGE=4 };
}

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolygonOMP::PairBodyRoundedPolygonOMP(LAMMPS *lmp) :
  PairBodyRoundedPolygon(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolygonOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  // grow per-atom body data structures if needed and initialize

  if (nall > nmax) {
    memory->destroy(dnum);
    memory->destroy(dfirst);
    memory->destroy(ednum);
    memory->destroy(edfirst);
    memory->destroy(enclosing_radius);
    memory->destroy(rounded_radius);
    nmax = nall;
    memory->create(dnum, nmax, "pair:dnum");
    memory->create(dfirst, nmax, "pair:dfirst");
    memory->create(ednum, nmax, "pair:ednum");
    memory->create(edfirst, nmax, "pair:edfirst");
    memory->create(enclosing_radius, nmax, "pair:enclosing_radius");
    memory->create(rounded_radius, nmax, "pair:rounded_radius");
  }

  ndiscrete = nedge = 0;
  for (int i = 0; i < nall; i++) dnum[i] = ednum[i] = 0;

  // pre-initialize body2space for all body atoms serially
  // body2space modifies shared arrays (discrete, dnum, dfirst, edge, ednum, edfirst)
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
        eval<1,1,1>(ifrom, ito, thr);
      } else {
        eval<1,0,1>(ifrom, ito, thr);
      }
    } else {
      eval<0,0,1>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairBodyRoundedPolygonOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i, j, ii, jj, jnum, itype, jtype;
  int ni, nj, npi, npj, ifirst, jfirst;
  int nei, nej, iefirst, jefirst;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl;
  double rsq, r, radi, radj, k_nij, k_naij;
  double facc[3];
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;

  const double * const * _noalias const x = atom->x;
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  auto * _noalias const torque = (dbl3_t *) thr->get_torque()[0];
  const double * const * _noalias const v = atom->v;
  const double * const * _noalias const angmom = atom->angmom;
  const double * _noalias const radius = atom->radius;
  const tagint * _noalias const tag = atom->tag;
  const int * _noalias const body = atom->body;
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // allocate thread-local scratch arrays for vertex and edge forces
  // these replace the shared discrete[*][3..5] and edge[*][2..4] members
  // to avoid data races between threads

  double *vforce_thr = nullptr;
  double *eforce_thr = nullptr;
  if (dmax > 0) memory->create(vforce_thr, dmax * 3, "pair:vforce_thr");
  if (edmax > 0) memory->create(eforce_thr, edmax * 3, "pair:eforce_thr");

  // loop over neighbors of my atoms
  // body sub-particle data has already been initialized serially

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (body[i] >= 0) {
      npi = dnum[i];
      ifirst = dfirst[i];
      nei = ednum[i];
      iefirst = edfirst[i];
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      radj = radius[j];

      evdwl = 0.0;
      facc[0] = facc[1] = facc[2] = 0;

      if (body[i] < 0 || body[j] < 0) continue;

      npj = dnum[j];
      jfirst = dfirst[j];
      nej = ednum[j];
      jefirst = edfirst[j];

      k_nij = k_n[itype][jtype];
      k_naij = k_na[itype][jtype];

      r = sqrt(rsq);
      if (r > radi + radj + cut_inner) continue;

      if (npi == 1 && npj == 1) {
        sphere_against_sphere_thr(i, j, delx, dely, delz, rsq, k_nij, k_naij, v, f, thr);
        continue;
      }

      // reset thread-local vertex and edge force scratch arrays

      if (vforce_thr) {
        for (ni = 0; ni < npi; ni++) {
          vforce_thr[(ifirst+ni)*3+0] = 0;
          vforce_thr[(ifirst+ni)*3+1] = 0;
          vforce_thr[(ifirst+ni)*3+2] = 0;
        }
        for (nj = 0; nj < npj; nj++) {
          vforce_thr[(jfirst+nj)*3+0] = 0;
          vforce_thr[(jfirst+nj)*3+1] = 0;
          vforce_thr[(jfirst+nj)*3+2] = 0;
        }
      }

      if (eforce_thr) {
        for (ni = 0; ni < nei; ni++) {
          eforce_thr[(iefirst+ni)*3+0] = 0;
          eforce_thr[(iefirst+ni)*3+1] = 0;
          eforce_thr[(iefirst+ni)*3+2] = 0;
        }
        for (nj = 0; nj < nej; nj++) {
          eforce_thr[(jefirst+nj)*3+0] = 0;
          eforce_thr[(jefirst+nj)*3+1] = 0;
          eforce_thr[(jefirst+nj)*3+2] = 0;
        }
      }

      int num_contacts, done;
      double delta_a, j_a;
      Contact contact_list[MAX_CONTACTS];

      num_contacts = 0;

      // check interaction between i's vertices and j's edges

      vertex_against_edge_thr(i, j, k_nij, k_naij, x, f, torque, tag,
                               contact_list, num_contacts, evdwl, facc,
                               vforce_thr, eforce_thr, thr);

      // check interaction between j's vertices and i's edges

      vertex_against_edge_thr(j, i, k_nij, k_naij, x, f, torque, tag,
                               contact_list, num_contacts, evdwl, facc,
                               vforce_thr, eforce_thr, thr);

      if (num_contacts >= 2) {

        done = 0;
        for (int m = 0; m < num_contacts-1; m++) {
          for (int n = m+1; n < num_contacts; n++) {
            delta_a = contact_separation(contact_list[m], contact_list[n]);
            if (delta_a > 0) {
              j_a = delta_a / (EFF_CONTACTS * delta_ua);
              if (j_a < 1.0) j_a = 1.0;

              contact_forces_thr(contact_list[m], j_a, x, v, angmom, f, torque,
                                 evdwl, facc, vforce_thr, eforce_thr, thr);
              contact_forces_thr(contact_list[n], j_a, x, v, angmom, f, torque,
                                 evdwl, facc, vforce_thr, eforce_thr, thr);
              done = 1;
              break;
            }
          }
          if (done == 1) break;
        }

      } else if (num_contacts == 1) {

        contact_forces_thr(contact_list[0], 1.0, x, v, angmom, f, torque,
                           evdwl, facc, vforce_thr, eforce_thr, thr);
      }

      if (EVFLAG) ev_tally_xyz_thr(this, i, j, nlocal, NEWTON_PAIR, evdwl, 0.0,
                                   facc[0], facc[1], facc[2], delx, dely, delz, thr);

    } // end for jj
  } // end for ii

  if (vforce_thr) memory->destroy(vforce_thr);
  if (eforce_thr) memory->destroy(eforce_thr);
}

/* ----------------------------------------------------------------------
   Sphere-sphere interaction using per-thread force array
------------------------------------------------------------------------- */

void PairBodyRoundedPolygonOMP::sphere_against_sphere_thr(int i, int j,
                                double delx, double dely, double delz,
                                double rsq, double k_n, double k_na,
                                const double * const * v, dbl3_t *f, ThrData *thr)
{
  double rradi, rradj;
  double vr1, vr2, vr3, vnnr, vn1, vn2, vn3, vt1, vt2, vt3;
  double rij, rsqinv, R, fx, fy, fz, fn[3], ft[3], fpair, shift, energy;
  int nlocal = atom->nlocal;

  rradi = rounded_radius[i];
  rradj = rounded_radius[j];

  rsqinv = 1.0/rsq;
  rij = sqrt(rsq);
  R = rij - (rradi + rradj);
  shift = k_na * cut_inner;

  energy = 0;

  if (R <= 0) {
    fpair = -k_n * R - shift;
    energy = (0.5 * k_n * R + shift) * R;
  } else if (R <= cut_inner) {
    fpair = k_na * R - shift;
    energy = (-0.5 * k_na * R + shift) * R;
  } else fpair = 0.0;

  fx = delx*fpair/rij;
  fy = dely*fpair/rij;
  fz = delz*fpair/rij;

  double rmin = MIN(rradi, rradj);
  if (R <= EPSILON*rmin) {

    vr1 = v[i][0] - v[j][0];
    vr2 = v[i][1] - v[j][1];
    vr3 = v[i][2] - v[j][2];

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

  // newton_pair is always 1 for this style (enforced by init_style)
  f[i].x += fx;
  f[i].y += fy;
  f[i].z += fz;

  if (j < nlocal) {
    f[j].x -= fx;
    f[j].y -= fy;
    f[j].z -= fz;
  }

  ev_tally_xyz_thr(this, i, j, nlocal, 1,
                   energy, 0.0, fx, fy, fz, delx, dely, delz, thr);
}

/* ----------------------------------------------------------------------
   Vertex-edge interaction using thread-local scratch arrays
------------------------------------------------------------------------- */

int PairBodyRoundedPolygonOMP::vertex_against_edge_thr(int i, int j,
                                  double k_n, double k_na,
                                  const double * const *x,
                                  dbl3_t *f, dbl3_t *torque,
                                  const tagint *tag,
                                  Contact *contact_list, int &num_contacts,
                                  double &evdwl, double *facc,
                                  double *vforce_thr, double *eforce_thr,
                                  ThrData * /*thr*/)
{
  int ni, npi, ifirst;
  int nj, jfirst, nej, jefirst;
  double xpi[3], xpj[3], dist, eradj, rradi, rradj;
  double fx, fy, fz, energy;
  int interact;

  npi = dnum[i];
  ifirst = dfirst[i];
  rradi = rounded_radius[i];

  jfirst = dfirst[j];
  nej = ednum[j];
  jefirst = edfirst[j];
  eradj = enclosing_radius[j];
  rradj = rounded_radius[j];

  energy = 0;
  interact = 0;

  for (ni = 0; ni < npi; ni++) {

    xpi[0] = x[i][0] + discrete[ifirst+ni][0];
    xpi[1] = x[i][1] + discrete[ifirst+ni][1];
    xpi[2] = x[i][2] + discrete[ifirst+ni][2];

    {
      double dx = xpi[0] - x[j][0];
      double dy = xpi[1] - x[j][1];
      double dz = xpi[2] - x[j][2];
      dist = sqrt(dx*dx + dy*dy + dz*dz);
    }

    if (dist > eradj + rradj + rradi + cut_inner) continue;

    int mode, contact, p2vertex;
    double d, R, hi[3], t, delx, dely, delz, fpair, shift;
    double rij;

    for (nj = 0; nj < nej; nj++) {

      mode = compute_distance_to_vertex(j, nj, const_cast<double *>(x[j]), rradj,
                                        xpi, rradi, cut_inner,
                                        d, hi, t, contact);

      if (mode == INVALID || mode == NONE) continue;

      if (mode == VERTEXI || mode == VERTEXJ) {

        interact = 1;

        if (mode == VERTEXI) p2vertex = (int)edge[jefirst+nj][0];
        else p2vertex = (int)edge[jefirst+nj][1];

        xpj[0] = x[j][0] + discrete[jfirst+p2vertex][0];
        xpj[1] = x[j][1] + discrete[jfirst+p2vertex][1];
        xpj[2] = x[j][2] + discrete[jfirst+p2vertex][2];

        delx = xpi[0] - xpj[0];
        dely = xpi[1] - xpj[1];
        delz = xpi[2] - xpj[2];

        rij = sqrt(delx*delx + dely*dely + delz*delz);
        R = rij - (rradi + rradj);
        shift = k_na * cut_inner;

        if (R <= 0) {
          fpair = -k_n * R - shift;
          energy += (0.5 * k_n * R + shift) * R;
        } else if (R <= cut_inner) {
          fpair = k_na * R - shift;
          energy += (-0.5 * k_na * R + shift) * R;
        } else fpair = 0.0;

        fx = delx*fpair/rij;
        fy = dely*fpair/rij;
        fz = delz*fpair/rij;

        if (tag[i] < tag[j] || npi == 1) {

          f[i].x += fx;
          f[i].y += fy;
          f[i].z += fz;
          sum_torque(const_cast<double *>(x[i]), xpi, fx, fy, fz, &torque[i].x);

          f[j].x -= fx;
          f[j].y -= fy;
          f[j].z -= fz;
          sum_torque(const_cast<double *>(x[j]), xpj, -fx, -fy, -fz, &torque[j].x);

          facc[0] += fx; facc[1] += fy; facc[2] += fz;
        }

      } else if (mode == EDGE) {

        interact = 1;

        delx = xpi[0] - hi[0];
        dely = xpi[1] - hi[1];
        delz = xpi[2] - hi[2];

        R = d - (rradi + rradj);
        shift = k_na * cut_inner;

        if (R <= 0) {
          fpair = -k_n * R - shift;
          energy += (0.5 * k_n * R + shift) * R;
        } else if (R <= cut_inner) {
          fpair = k_na * R - shift;
          energy += (-0.5 * k_na * R + shift) * R;
        } else fpair = 0.0;

        fx = delx*fpair/d;
        fy = dely*fpair/d;
        fz = delz*fpair/d;

        if (contact == 1) {

          contact_list[num_contacts].ibody = i;
          contact_list[num_contacts].jbody = j;
          contact_list[num_contacts].vertex = ni;
          contact_list[num_contacts].edge = nj;
          contact_list[num_contacts].xv[0] = xpi[0];
          contact_list[num_contacts].xv[1] = xpi[1];
          contact_list[num_contacts].xv[2] = xpi[2];
          contact_list[num_contacts].xe[0] = hi[0];
          contact_list[num_contacts].xe[1] = hi[1];
          contact_list[num_contacts].xe[2] = hi[2];
          contact_list[num_contacts].separation = R;
          num_contacts++;

          // store forces to thread-local scratch arrays for later accumulation

          if (vforce_thr) {
            vforce_thr[(ifirst+ni)*3+0] = fx;
            vforce_thr[(ifirst+ni)*3+1] = fy;
            vforce_thr[(ifirst+ni)*3+2] = fz;
          }

          if (eforce_thr) {
            eforce_thr[(jefirst+nj)*3+0] = -fx;
            eforce_thr[(jefirst+nj)*3+1] = -fy;
            eforce_thr[(jefirst+nj)*3+2] = -fz;
          }

        } else {

          // accumulate force and torque to both bodies directly

          f[i].x += fx;
          f[i].y += fy;
          f[i].z += fz;
          sum_torque(const_cast<double *>(x[i]), xpi, fx, fy, fz, &torque[i].x);

          f[j].x -= fx;
          f[j].y -= fy;
          f[j].z -= fz;
          sum_torque(const_cast<double *>(x[j]), hi, -fx, -fy, -fz, &torque[j].x);

          facc[0] += fx; facc[1] += fy; facc[2] += fz;
        }
      }
    } // end for edges of body j
  } // end for vertices of body i

  evdwl += energy;

  return interact;
}

/* ----------------------------------------------------------------------
   Contact forces using thread-local scratch arrays
------------------------------------------------------------------------- */

void PairBodyRoundedPolygonOMP::contact_forces_thr(Contact &contact, double j_a,
                                  const double * const *x,
                                  const double * const *v,
                                  const double * const *angmom,
                                  dbl3_t *f, dbl3_t *torque,
                                  double & /*evdwl*/, double *facc,
                                  double *vforce_thr, double *eforce_thr,
                                  ThrData * /*thr*/)
/* const_cast is safe: total_velocity and sum_torque do not modify their
   pointer arguments but lack const in their signatures */
{
  int ibody, jbody, ibonus, jbonus, ifirst, jefirst, ni, nj;
  double fx, fy, fz, delx, dely, delz, rsq, rsqinv;
  double vr1, vr2, vr3, vnnr, vn1, vn2, vn3, vt1, vt2, vt3;
  double fn[3], ft[3], vi[3], vj[3];
  double *quat, *inertia;
  AtomVecBody::Bonus *bonus;

  ibody = contact.ibody;
  jbody = contact.jbody;

  ibonus = atom->body[ibody];
  bonus = &avec->bonus[ibonus];
  quat = bonus->quat;
  inertia = bonus->inertia;
  total_velocity(contact.xv, const_cast<double *>(x[ibody]),
                 const_cast<double *>(v[ibody]),
                 const_cast<double *>(angmom[ibody]), inertia, quat, vi);

  jbonus = atom->body[jbody];
  bonus = &avec->bonus[jbonus];
  quat = bonus->quat;
  inertia = bonus->inertia;
  total_velocity(contact.xe, const_cast<double *>(x[jbody]),
                 const_cast<double *>(v[jbody]),
                 const_cast<double *>(angmom[jbody]), inertia, quat, vj);

  delx = contact.xv[0] - contact.xe[0];
  dely = contact.xv[1] - contact.xe[1];
  delz = contact.xv[2] - contact.xe[2];
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

  ifirst = dfirst[ibody];
  ni = contact.vertex;

  double vfx = 0, vfy = 0, vfz = 0;
  if (vforce_thr) {
    vfx = vforce_thr[(ifirst+ni)*3+0];
    vfy = vforce_thr[(ifirst+ni)*3+1];
    vfz = vforce_thr[(ifirst+ni)*3+2];
  }

  fx = vfx * j_a + fn[0] + ft[0] + mu * vfx;
  fy = vfy * j_a + fn[1] + ft[1] + mu * vfy;
  fz = vfz * j_a + fn[2] + ft[2] + mu * vfz;

  f[ibody].x += fx;
  f[ibody].y += fy;
  f[ibody].z += fz;
  sum_torque(const_cast<double *>(x[ibody]), contact.xv, fx, fy, fz, &torque[ibody].x);

  facc[0] += fx; facc[1] += fy; facc[2] += fz;

  jefirst = edfirst[jbody];
  nj = contact.edge;

  double efx = 0, efy = 0, efz = 0;
  if (eforce_thr) {
    efx = eforce_thr[(jefirst+nj)*3+0];
    efy = eforce_thr[(jefirst+nj)*3+1];
    efz = eforce_thr[(jefirst+nj)*3+2];
  }

  fx = efx * j_a - fn[0] - ft[0] + mu * efx;
  fy = efy * j_a - fn[1] - ft[1] + mu * efy;
  fz = efz * j_a - fn[2] - ft[2] + mu * efz;

  f[jbody].x += fx;
  f[jbody].y += fy;
  f[jbody].z += fz;
  sum_torque(const_cast<double *>(x[jbody]), contact.xe, fx, fy, fz, &torque[jbody].x);
}

/* ---------------------------------------------------------------------- */

double PairBodyRoundedPolygonOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairBodyRoundedPolygon::memory_usage();
  return bytes;
}
