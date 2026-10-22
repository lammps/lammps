// clang-format off
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
   Contributing author: Trung Dac Nguyen (ndactrung@gmail.com)
   Ref: Wang, Yu, Langston, Fraige, Particle shape effects in discrete
   element modelling of cohesive angular particles, Granular Matter 2011,
   13:1-12.
   Note: The current implementation has not taken into account
         the contact history for friction forces.
------------------------------------------------------------------------- */

#include "pair_body_rounded_polyhedron.h"

#include "atom.h"
#include "atom_vec_body.h"
#include "body_rounded_polyhedron.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "fix_neigh_history.h"
#include "fix_store_atom.h"
#include "force.h"
#include "group.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <utility>

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr int DELTA = 10000;
static constexpr double BIG = 1.0e20;
static constexpr double EPSILON = 1.0e-3; // dimensionless threshold (dot products, end point checks, contact checks)
static constexpr double PARALLEL_TOL = 1.0e-12;   // 1 - |cos| below which two edges are parallel
static constexpr double CONE_TOL = 1.0e-10;       // tolerance of the normal cone tests
static constexpr double PATCH_ANGLE = 5.0 * DEG2RAD;   // see patch_factor()
static constexpr double LINE_ANGLE = 5.0 * DEG2RAD;    // see edge_edge_weight()
static const double COS_PATCH_ANGLE = cos(PATCH_ANGLE);
static const double COS_HALF_PATCH_ANGLE = cos(0.5 * PATCH_ANGLE);
static const double COS_LINE_ANGLE = cos(LINE_ANGLE);
static constexpr int MAX_FACE_SIZE = 4;   // maximum number of vertices per face (same as BodyRoundedPolyhedron)
static constexpr int NFNC = 12;           // per-body force and torque of the j_a scaling, and of damping
static constexpr char id_fix_store_prefix[] = "BODY_ROUNDED_POLYHEDRON_WORK_";

//#define _POLYHEDRON_DEBUG

enum {EE_INVALID=0,EE_NONE,EE_INTERACT};
enum {EF_INVALID=0,EF_NONE,EF_INTERACT};

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolyhedron::PairBodyRoundedPolyhedron(LAMMPS *lmp) :
    Pair(lmp), avec(nullptr), bptr(nullptr)
{
  dmax = nmax = 0;
  discrete = nullptr;
  dnum = dfirst = nullptr;

  edmax = ednummax = 0;
  edge = nullptr;
  ednum = edfirst = nullptr;

  facmax = facnummax = 0;
  face = nullptr;
  facnorm = facplane = nullptr;
  facsize = nullptr;
  facnum = facfirst = nullptr;

  enclosing_radius = nullptr;
  rounded_radius = nullptr;
  maxrad = nullptr;

  single_enable = 0;
  restartinfo = 0;

  // work done by the forces that do not derive from the energy,
  // accessible via compute pair

  nextra = 2;
  pvector = new double[nextra];
  w_ja = w_diss = 0.0;

  fnc = nullptr;
  nmax_fnc = 0;
  id_fix_store = nullptr;
  fix_store = nullptr;
  comm_reverse = NFNC;

  c_n = 0.1;
  c_t = 0.2;
  mu = 0.0;
  A_ua = 1.0;

  k_n = nullptr;
  k_na = nullptr;
  k_t = nullptr;

  // the tangential displacement is stored only if requested, see settings()
  // the surfaces of two bodies interact beyond contact up to cut_inner,
  // so the history is kept for all pairs in the neighbor list

  history = 0;
  dt = 0.0;
  id_history = nullptr;
  fix_history = nullptr;
  beyond_contact = 1;
}

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolyhedron::~PairBodyRoundedPolyhedron()
{
  memory->destroy(discrete);
  memory->destroy(dnum);
  memory->destroy(dfirst);

  memory->destroy(edge);
  memory->destroy(ednum);
  memory->destroy(edfirst);

  memory->destroy(face);
  memory->destroy(facnorm);
  memory->destroy(facplane);
  memory->destroy(facsize);
  memory->destroy(facnum);
  memory->destroy(facfirst);

  memory->destroy(enclosing_radius);
  memory->destroy(rounded_radius);
  memory->destroy(maxrad);

  delete[] pvector;
  memory->destroy(fnc);
  if (id_fix_store && modify) modify->delete_fix(id_fix_store);
  delete[] id_fix_store;

  if (modify) {
    if (fix_history) modify->delete_fix(id_history);
    else history_dummy_fix(0);
  }
  delete[] id_history;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(k_n);
    memory->destroy(k_na);
    memory->destroy(k_t);
  }
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,facc[3];
  double rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **torque = atom->torque;
  double **angmom = atom->angmom;
  double *radius = atom->radius;
  int *body = atom->body;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // grow the per-atom lists if necessary and initialize

  if (atom->nmax > nmax) {
    memory->destroy(dnum);
    memory->destroy(dfirst);
    memory->destroy(ednum);
    memory->destroy(edfirst);
    memory->destroy(facnum);
    memory->destroy(facfirst);
    memory->destroy(enclosing_radius);
    memory->destroy(rounded_radius);
    nmax = atom->nmax;
    memory->create(dnum,nmax,"pair:dnum");
    memory->create(dfirst,nmax,"pair:dfirst");
    memory->create(ednum,nmax,"pair:ednum");
    memory->create(edfirst,nmax,"pair:edfirst");
    memory->create(facnum,nmax,"pair:facnum");
    memory->create(facfirst,nmax,"pair:facfirst");
    memory->create(enclosing_radius,nmax,"pair:enclosing_radius");
    memory->create(rounded_radius,nmax,"pair:rounded_radius");
  }

  ndiscrete = nedge = nface = 0;
  for (i = 0; i < nall; i++)
    dnum[i] = ednum[i] = facnum[i] = 0;

  // per-body forces and torques that do not derive from the energy

  if (atom->nmax > nmax_fnc) {
    memory->destroy(fnc);
    nmax_fnc = atom->nmax;
    memory->create(fnc,nmax_fnc,NFNC,"pair:fnc");
  }
  for (i = 0; i < nall; i++)
    for (int k = 0; k < NFNC; k++) fnc[i][k] = 0.0;

  // tangential displacements of the pairs in the neighbor list

  int *touch = nullptr, **firsttouch = nullptr;
  double *allshear = nullptr, **firstshear = nullptr;
  if (history) {
    firsttouch = fix_history->firstflag;
    firstshear = fix_history->firstvalue;
  }
  scratch.shear = nullptr;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    if (history) {
      touch = firsttouch[i];
      allshear = firstshear[i];
    }

    if ((body[i] >= 0) && (dnum[i] == 0)) body2space(i);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      // body/body interactions

      evdwl = 0.0;
      facc[0] = facc[1] = facc[2] = 0;

      // the tangential displacement is reset unless the pair is in contact

      if (history) {
        scratch.shear = &allshear[3*jj];
        scratch.shear_i = i;
        scratch.touched = 0;
        scratch.fnsum = 0.0;
        for (int k = 0; k < 3; k++) scratch.pcsum[k] = scratch.vtsum[k] = scratch.nsum[k] = 0.0;
      }

      // no interaction, radius = enclosing + rounded radius

      double r = sqrt(rsq);
      if ((body[i] >= 0) && (body[j] >= 0) && (r <= radius[i] + radius[j] + cut_inner)) {
        if (dnum[j] == 0) body2space(j);
        pair_interaction(i, j, delx, dely, delz, rsq, x, v, angmom, f, torque, fnc,
                         scratch, evdwl, facc);
        if (history) history_friction(i, j, x, f, torque, fnc, scratch, facc);

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,evdwl,0.0,
                                 facc[0],facc[1],facc[2],delx,dely,delz);
      }

      if (history) {
        touch[jj] = scratch.touched;
        if (!scratch.touched) scratch.shear[0] = scratch.shear[1] = scratch.shear[2] = 0.0;
      }

    } // end for jj
  }

  if (vflag_fdotr) virial_fdotr_compute();

  work_nonconservative();
}

/* ----------------------------------------------------------------------
   Compute the interaction between bodies i and j:
   accumulate the forces and torques to f and torque, the forces and torques
   that do not derive from the energy to fnc, and return the energy in evdwl
   and the total force on body i in facc
   f, torque, fnc and the scratch space s may be per-thread storage
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::pair_interaction(int i, int j, double delx, double dely,
                                                 double delz, double rsq, double **x,
                                                 double **v, double **angmom, double **f,
                                                 double **torque, double **fnc, Scratch &s,
                                                 double &evdwl, double *facc)
{
  int itype = atom->type[i];
  int jtype = atom->type[j];
  int npi = dnum[i];
  int ifirst = dfirst[i];
  int npj = dnum[j];
  int jfirst = dfirst[j];
  std::vector<Contact> &contacts = s.contacts;

  // sphere-sphere interaction

  if (npi == 1 && npj == 1) {
    sphere_against_sphere(i, j, itype, jtype, delx, dely, delz, rsq, x, v, angmom, f, torque,
                          fnc, s, evdwl, facc);
    return;
  }

  // reset the flags of the vertices already interacted with

  if ((int) s.vertex_done.size() < ndiscrete) s.vertex_done.resize(ndiscrete);
  for (int ni = 0; ni < npi; ni++) s.vertex_done[ifirst+ni] = 0;
  for (int nj = 0; nj < npj; nj++) s.vertex_done[jfirst+nj] = 0;

  // one of the two bodies is a sphere
  // one friction force per pair of bodies, see friction_forces()

  if ((npj == 1) || (npi == 1)) {
    contacts.clear();
    if (npj == 1) {
      sphere_against_polyhedron(i, j, itype, jtype, x, v, f, torque, angmom, fnc, s,
                                evdwl, facc);
    } else {

      // the force on body j is returned, facc is the force on body i

      double fj[3] = {0.0, 0.0, 0.0};
      sphere_against_polyhedron(j, i, jtype, itype, x, v, f, torque, angmom, fnc, s,
                                evdwl, fj);
      facc[0] -= fj[0];
      facc[1] -= fj[1];
      facc[2] -= fj[2];
    }
    if (!contacts.empty()) {
      friction_forces(itype, jtype, x, v, angmom, f, torque, fnc, i, facc, s);
    }
    return;
  }

  // reach of the vertices of each body towards the other body along the
  // line between the centers: all points of body j are at least
  // r - max reach of j away from the center plane of body i, so a feature of
  // body i whose vertices do not reach r - max reach of j - cut, cannot
  // interact with body j and vice versa, and the bodies cannot interact at all
  // if the gap between their extents exceeds cut
  // cut includes a margin for rounding, beyond which all interactions vanish

  double r = sqrt(rsq);
  double cut = (rounded_radius[i] + rounded_radius[j] + cut_inner) * (1.0 + EPSILON);
  if ((int) s.reach.size() < ndiscrete) s.reach.resize(ndiscrete);
  s.ibody = i;
  s.jbody = j;
  if (r > 0.0) {
    double u[3] = {-delx/r, -dely/r, -delz/r};
    double reachi = MathExtra::dot3(discrete[ifirst], u);
    for (int ni = 0; ni < npi; ni++) {
      s.reach[ifirst+ni] = MathExtra::dot3(discrete[ifirst+ni], u);
      reachi = MAX(reachi, s.reach[ifirst+ni]);
    }
    double reachj = -MathExtra::dot3(discrete[jfirst], u);
    for (int nj = 0; nj < npj; nj++) {
      s.reach[jfirst+nj] = -MathExtra::dot3(discrete[jfirst+nj], u);
      reachj = MAX(reachj, s.reach[jfirst+nj]);
    }
    if (r - reachi - reachj > cut) return;
    s.reach_min_i = r - reachj - cut;
    s.reach_min_j = r - reachi - cut;
  } else {
    for (int ni = 0; ni < npi; ni++) s.reach[ifirst+ni] = 0.0;
    for (int nj = 0; nj < npj; nj++) s.reach[jfirst+nj] = 0.0;
    s.reach_min_i = s.reach_min_j = -cut;
  }

  contacts.clear();

  // facc is the force on body i, fj collects the forces on body j
  // returned by the routines below, which is subtracted at the end

  double fj[3] = {0.0, 0.0, 0.0};

  // corners of the regions of overlap of nearly parallel faces, first, since
  // the weights of the other contacts depend on them, see vertex_face_weight()

  face_face_patches(i, j, itype, jtype, x, v, f, torque, angmom, fnc, s, evdwl, facc);

  // vertices inside the core of the other body, for large overlaps

  vertex_against_core(i, j, itype, jtype, x, v, f, torque, angmom, fnc, s, evdwl, facc);
  vertex_against_core(j, i, jtype, itype, x, v, f, torque, angmom, fnc, s, evdwl, fj);

  // check interaction between i's edges and j' faces
  #ifdef _POLYHEDRON_DEBUG
  printf("INTERACTION between edges of %d vs. faces of %d:\n", i, j);
  #endif
  edge_against_face(i, j, itype, jtype, x, v, f, torque, angmom,
                    fnc, s, evdwl, facc);

  // check interaction between j's edges and i' faces
  #ifdef _POLYHEDRON_DEBUG
  printf("\nINTERACTION between edges of %d vs. faces of %d:\n", j, i);
  #endif
  edge_against_face(j, i, jtype, itype, x, v, f, torque, angmom,
                    fnc, s, evdwl, fj);

  // vertices without a vertex-face interaction interact with the edges,
  // and else with the vertices of the other body, see Fig. 7, Wang et al.

  vertex_against_edge(i, j, itype, jtype, x, v, f, torque, angmom, fnc, s, evdwl, facc);
  vertex_against_edge(j, i, jtype, itype, x, v, f, torque, angmom, fnc, s, evdwl, fj);
  vertex_against_vertex(i, j, itype, jtype, x, v, f, torque, angmom, fnc, s, evdwl, facc);

  // check interaction between i's edges and j' edges, which returns
  // the force on body j
  #ifdef _POLYHEDRON_DEBUG
  printf("INTERACTION between edges of %d vs. edges of %d:\n", i, j);
  #endif
  edge_against_edge(i, j, itype, jtype, x, v, f, torque, angmom,
                    fnc, s, evdwl, fj);


  // estimate the contact area
  // also consider point contacts and line contacts

  if (!contacts.empty()) {
    rescale_cohesive_forces(x, v, angmom, f, torque, fnc, contacts, itype, jtype, i, evdwl,
                            facc);

    // one friction force per pair of bodies, see friction_forces()

    friction_forces(itype, jtype, x, v, angmom, f, torque, fnc, i, facc, s);
  }

  facc[0] -= fj[0];
  facc[1] -= fj[1];
  facc[2] -= fj[2];
}

/* ----------------------------------------------------------------------
   accumulate the work done by the forces that do not derive from the
   energy: the part of the contact forces added by the j_a scaling (w_ja),
   and damping and friction (w_diss), so that the total energy minus
   w_ja and w_diss is conserved
   with velocity Verlet, the work over a time step is
     0.5 * dt * (F(n) + F(n+1)) . v(n+1/2)
   with F(n) of the previous step stored per atom with fix STORE/ATOM
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::work_nonconservative()
{
  if (force->newton_pair) comm->reverse_comm(this);

  double **fprev = fix_store->astore;
  double **v = atom->v;
  double **angmom = atom->angmom;
  int *body = atom->body;
  int nlocal = atom->nlocal;
  double halfdt = 0.5 * dt;
  double omega[3],ex[3],ey[3],ez[3];

  for (int i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;

    // no time step has been taken during setup, nor during minimization

    if (!update->setupflag && (update->whichflag == 1)) {
      AtomVecBody::Bonus *bonus = &avec->bonus[body[i]];
      MathExtra::q_to_exyz(bonus->quat,ex,ey,ez);
      MathExtra::angmom_to_omega(angmom[i],ex,ey,ez,bonus->inertia,omega);
      for (int k = 0; k < 3; k++) {
        w_ja += halfdt * ((fprev[i][k] + fnc[i][k]) * v[i][k] +
                          (fprev[i][3+k] + fnc[i][3+k]) * omega[k]);
        w_diss += halfdt * ((fprev[i][6+k] + fnc[i][6+k]) * v[i][k] +
                            (fprev[i][9+k] + fnc[i][9+k]) * omega[k]);
      }
    }
    for (int k = 0; k < NFNC; k++) fprev[i][k] = fnc[i][k];
  }

  pvector[0] = w_ja;
  pvector[1] = w_diss;
}

/* ---------------------------------------------------------------------- */

int PairBodyRoundedPolyhedron::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++)
    for (int k = 0; k < NFNC; k++) buf[m++] = fnc[i][k];
  return m;
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    for (int k = 0; k < NFNC; k++) fnc[j][k] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(k_n,n+1,n+1,"pair:k_n");
  memory->create(k_na,n+1,n+1,"pair:k_na");
  memory->create(k_t,n+1,n+1,"pair:k_t");
  memory->create(maxrad,n+1,"pair:maxrad");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::settings(int narg, char **arg)
{
  if (narg < 5) error->all(FLERR,"Illegal pair_style command");

  c_n = utils::numeric(FLERR,arg[0],false,lmp);
  c_t = utils::numeric(FLERR,arg[1],false,lmp);
  mu = utils::numeric(FLERR,arg[2],false,lmp);
  A_ua = utils::numeric(FLERR,arg[3],false,lmp);
  cut_inner = utils::numeric(FLERR,arg[4],false,lmp);

  if (A_ua < 0) A_ua = 1;

  int history_one = 0;
  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"history") == 0) {
      history_one = 1;
      iarg++;
    } else error->all(FLERR, iarg, "Unknown pair_style body/rounded/polyhedron keyword {}",
                      arg[iarg]);
  }

  // create a placeholder for fix NEIGH_HISTORY, which replaces it in init_style(),
  // so that the order of the fixes follows the input, or remove the history

  if (history_one && !history) {
    history_dummy_fix(1);
  } else if (!history_one && history) {
    if (fix_history) modify->delete_fix(id_history);
    else history_dummy_fix(0);
    fix_history = nullptr;
  }
  history = history_one;
}

/* ----------------------------------------------------------------------
   create (flag = 1) or delete (flag = 0) the placeholder of fix NEIGH_HISTORY
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::history_dummy_fix(int flag)
{
  std::string id_dummy = "NEIGH_HISTORY_BODY_DUMMY" + std::to_string(instance_me);
  if (flag) modify->add_fix(id_dummy + " all DUMMY");
  else if (modify->get_fix_by_id(id_dummy)) modify->delete_fix(id_dummy);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double k_n_one = utils::numeric(FLERR,arg[2],false,lmp);
  double k_na_one = utils::numeric(FLERR,arg[3],false,lmp);

  // the tangential stiffness defaults to 2/7 of the normal stiffness,
  // as in pair style gran/hooke/history

  double k_t_one = 2.0/7.0 * k_n_one;
  if (narg == 5) k_t_one = utils::numeric(FLERR,arg[4],false,lmp);
  if (k_t_one < 0.0)
    error->all(FLERR, 4, "Tangential stiffness of pair style body/rounded/polyhedron "
               "must not be negative");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      k_n[i][j] = k_n_one;
      k_na[i][j] = k_na_one;
      k_t[i][j] = k_t_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::init_style()
{
  avec = dynamic_cast<AtomVecBody *>(atom->style_match("body"));
  if (!avec) error->all(FLERR,"Pair body/rounded/polyhedron requires atom style body");
  if (strcmp(avec->bptr->style,"rounded/polyhedron") != 0)
    error->all(FLERR,"Pair body/rounded/polyhedron requires body style rounded/polyhedron");
  bptr = dynamic_cast<BodyRoundedPolyhedron *>(avec->bptr);

  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style body/rounded/polyhedron requires newton pair on");

  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair body/rounded/polyhedron requires ghost atoms store velocity");

  // the neighbor list keeps the regular cutoff, which includes cut_inner,
  // also with the contact history

  if (history) neighbor->add_request(this, NeighConst::REQ_HISTORY);
  else neighbor->add_request(this);

  reset_dt();

  // on the first init, fix NEIGH_HISTORY replaces the placeholder created in
  // settings(), so that its position in the list of fixes is preserved

  if (history) {
    delete[] id_history;
    id_history = utils::strdup(fmt::format("NEIGH_HISTORY_BODY{}", instance_index()));
    if (!fix_history) {
      fix_history = dynamic_cast<FixNeighHistory *>(
        modify->replace_fix("NEIGH_HISTORY_BODY_DUMMY" + std::to_string(instance_me),
                            fmt::format("{} all NEIGH_HISTORY 3", id_history), 1));
      fix_history->pair = this;
    } else {
      fix_history = dynamic_cast<FixNeighHistory *>(modify->get_fix_by_id(id_history));
      if (!fix_history) error->all(FLERR, "Could not find pair fix neigh history ID");
    }
  }

  // per-atom storage of the forces and torques that do not derive
  // from the energy at the previous step, see work_nonconservative()

  if (!id_fix_store) {
    id_fix_store = utils::strdup(std::string(id_fix_store_prefix) + std::to_string(instance_me));
    fix_store = dynamic_cast<FixStoreAtom *>(
      modify->add_fix(fmt::format("{} {} STORE/ATOM {} 0 0 0", id_fix_store,
                                  group->names[0], NFNC)));

    // zero the stored values of the bodies inserted during a run, e.g. by
    // fix pour or fix deposit, instead of keeping those of a former atom

    fix_store->create_attribute = 1;
  } else {
    fix_store = dynamic_cast<FixStoreAtom *>(modify->get_fix_by_id(id_fix_store));
    if (!fix_store)
      error->all(FLERR, "Could not find internal fix STORE/ATOM id {}", id_fix_store);
  }

  // find the maximum radius (enclosing + rounded) for each atom type

  int i, itype;
  int* body = atom->body;
  double *radius = atom->radius;
  int* type = atom->type;
  int ntypes = atom->ntypes;
  int nlocal = atom->nlocal;

  if (atom->nmax > nmax) {
    memory->destroy(dnum);
    memory->destroy(dfirst);
    memory->destroy(ednum);
    memory->destroy(edfirst);
    memory->destroy(facnum);
    memory->destroy(facfirst);
    memory->destroy(enclosing_radius);
    memory->destroy(rounded_radius);
    nmax = atom->nmax;
    memory->create(dnum,nmax,"pair:dnum");
    memory->create(dfirst,nmax,"pair:dfirst");
    memory->create(ednum,nmax,"pair:ednum");
    memory->create(edfirst,nmax,"pair:edfirst");
    memory->create(facnum,nmax,"pair:facnum");
    memory->create(facfirst,nmax,"pair:facfirst");
    memory->create(enclosing_radius,nmax,"pair:enclosing_radius");
    memory->create(rounded_radius,nmax,"pair:rounded_radius");
  }

  ndiscrete = nedge = nface = 0;
  for (i = 0; i < nlocal; i++)
    dnum[i] = ednum[i] = facnum[i] = 0;

  double *mrad = nullptr;
  memory->create(mrad,ntypes+1,"pair:mrad");
  for (i = 1; i <= ntypes; i++)
    maxrad[i] = mrad[i] = 0;

  Fix *fixpour = nullptr;
  auto pours = modify->get_fix_by_style("^pour");
  if (!pours.empty()) fixpour = pours[0];

  Fix *fixdep = nullptr;
  auto deps = modify->get_fix_by_style("^deposit");
  if (!deps.empty()) fixdep = deps[0];

  for (i = 1; i <= ntypes; i++) {
    mrad[i] = 0.0;
    if (fixpour) {
      itype = i;
      mrad[i] = *((double *) fixpour->extract("radius",itype));
    }
    if (fixdep) {
      itype = i;
      mrad[i] = *((double *) fixdep->extract("radius",itype));
    }
  }

  for (i = 0; i < nlocal; i++) {
    itype = type[i];
    if ((body[i] >= 0) && (radius[i] > mrad[itype])) mrad[itype] = radius[i];
  }

  MPI_Allreduce(&mrad[1],&maxrad[1],ntypes,MPI_DOUBLE,MPI_MAX,world);

  memory->destroy(mrad);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBodyRoundedPolyhedron::init_one(int i, int j)
{
  k_n[j][i] = k_n[i][j];
  k_na[j][i] = k_na[i][j];
  k_t[j][i] = k_t[i][j];

  // the surfaces of two bodies interact up to cut_inner

  return (maxrad[i]+maxrad[j]+cut_inner);
}

/* ----------------------------------------------------------------------
   convert N sub-particles in body I to space frame using current quaternion
   store sub-particle space-frame displacements from COM in discrete list
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::body2space(int i)
{
  int ibonus = atom->body[i];
  AtomVecBody::Bonus *bonus = &avec->bonus[ibonus];
  int nsub = bptr->nsub(bonus);
  double *coords = bptr->coords(bonus);
  int body_num_edges = bptr->nedges(bonus);
  double* edge_ends = bptr->edges(bonus);
  int body_num_faces = bptr->nfaces(bonus);
  double* face_pts = bptr->faces(bonus);
  double eradius = bptr->enclosing_radius(bonus);
  double rradius = bptr->rounded_radius(bonus);

  // get the number of sub-particles (vertices)
  // and the index of the first vertex of my body in the list

  dnum[i] = nsub;
  dfirst[i] = ndiscrete;

  // grow the vertex list if necessary

  if (ndiscrete + nsub > dmax) {
    dmax += DELTA;
    memory->grow(discrete,dmax,3,"pair:discrete");
  }

  double p[3][3];
  MathExtra::quat_to_mat(bonus->quat,p);

  for (int m = 0; m < nsub; m++) {
    MathExtra::matvec(p,&coords[3*m],discrete[ndiscrete]);
    ndiscrete++;
  }

  // get the number of edges (vertices)
  // and the index of the first edge of my body in the list

  ednum[i] = body_num_edges;
  edfirst[i] = nedge;

  // grow the edge list if necessary
  // the 2 columns are for vertex indices within body

  if (nedge + body_num_edges > edmax) {
    edmax += DELTA;
    memory->grow(edge,edmax,2,"pair:edge");
  }

  if ((body_num_edges > 0) && (edge_ends == nullptr))
    error->one(FLERR,"Inconsistent edge data for body of atom {}", atom->tag[i]);

  for (int m = 0; m < body_num_edges; m++) {
    edge[nedge][0] = static_cast<int>(edge_ends[2*m+0]);
    edge[nedge][1] = static_cast<int>(edge_ends[2*m+1]);
    nedge++;
  }

  // get the number of faces and the index of the first face

  facnum[i] = body_num_faces;
  facfirst[i] = nface;

  // grow the face list if necessary
  // the first 3 columns are for vertex indices within body, the last 3 for forces

  if (nface + body_num_faces > facmax) {
    facmax += DELTA;
    memory->grow(face,facmax,MAX_FACE_SIZE,"pair:face");
    memory->grow(facnorm,facmax,3,"pair:facnorm");
    memory->grow(facplane,facmax,6,"pair:facplane");
    memory->grow(facsize,facmax,"pair:facsize");
  }

  if ((body_num_faces > 0) && (face_pts == nullptr))
    error->one(FLERR,"Inconsistent face data for body of atom {}", atom->tag[i]);

  for (int m = 0; m < body_num_faces; m++) {
    for (int k = 0; k < MAX_FACE_SIZE; k++)
      face[nface][k] = static_cast<int>(face_pts[MAX_FACE_SIZE*m+k]);
    nface++;
  }

  // geometry of the faces, which is used many times for each pair of bodies,
  // see face_size(), face_normal(), and nearest_face()

  double **x = atom->x;
  for (int m = facfirst[i]; m < nface; m++) {
    int n = 0;
    while ((n < MAX_FACE_SIZE) && (static_cast<int>(face[m][n]) >= 0)) n++;
    facsize[m] = n;

    // outward unit normal from the vertices relative to the center

    double *x1 = discrete[dfirst[i]+static_cast<int>(face[m][0])];
    double *x2 = discrete[dfirst[i]+static_cast<int>(face[m][1])];
    double *x3 = discrete[dfirst[i]+static_cast<int>(face[m][2])];
    double u[3], v[3];
    MathExtra::sub3(x2, x1, u);
    MathExtra::sub3(x3, x1, v);
    MathExtra::cross3(u, v, facnorm[m]);
    MathExtra::norm3(facnorm[m]);
    if (MathExtra::dot3(x1, facnorm[m]) < 0.0) MathExtra::negate3(facnorm[m]);

    // first vertex and outward unit normal from the vertices in space

    double xi1[3], xi2[3], xi3[3], xc[3], ans[3];
    double *nf = &facplane[m][3];
    MathExtra::add3(x[i], x1, xi1);
    MathExtra::add3(x[i], x2, xi2);
    MathExtra::add3(x[i], x3, xi3);
    MathExtra::sub3(xi2, xi1, u);
    MathExtra::sub3(xi3, xi1, v);
    MathExtra::cross3(u, v, nf);
    MathExtra::norm3(nf);
    xc[0] = (xi1[0] + xi2[0] + xi3[0])/3.0;
    xc[1] = (xi1[1] + xi2[1] + xi3[1])/3.0;
    xc[2] = (xi1[2] + xi2[2] + xi3[2])/3.0;
    MathExtra::sub3(xc, x[i], ans);
    if (MathExtra::dot3(ans, nf) < 0) MathExtra::negate3(nf);
    MathExtra::copy3(xi1, facplane[m]);
  }

  enclosing_radius[i] = eradius;
  rounded_radius[i] = rradius;
}

/* ----------------------------------------------------------------------
   Interaction between two spheres with different radii
   according to the 2D model from Fraige et al.
---------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::sphere_against_sphere(int ibody, int jbody,
  int itype, int jtype, double delx, double dely, double delz, double rsq,
  double** x, double** v, double** angmom, double** f, double** torque, double** fnc,
  Scratch &s, double &evdwl, double* facc)
{
  double rradi,rradj,contact_dist;
  double rij,R,fx,fy,fz,fpair,energy;
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

  f[ibody][0] += fx;
  f[ibody][1] += fy;
  f[ibody][2] += fz;

  if (newton_pair || jbody < nlocal) {
    f[jbody][0] -= fx;
    f[jbody][1] -= fy;
    f[jbody][2] -= fz;
  }

  evdwl += energy;
  facc[0] += fx; facc[1] += fy; facc[2] += fz;

  // damping and friction at the contact point between the surfaces

  if (R <= 0) {
    double n[3] = {delx/rij, dely/rij, delz/rij};
    double pc[3];
    contact_point(x[ibody], x[jbody], n, rradi, rradj, pc);
    double fne = -k_n[itype][jtype] * R;
    damping_friction(ibody, jbody, pc, n, fne, 1, 1, x, v, angmom, f, torque, fnc, ibody, facc,
                     &s);
  }
}

/* ----------------------------------------------------------------------
   Interaction bt a polyhedron (ibody) and a sphere (jbody): the rounded
   polyhedron and the sphere touch at a single point, the point of the
   polyhedron nearest to the center of the sphere, which is the nearest
   of the projections of the center onto the faces, if inside the face
   and in front of it, and onto the edges, limited to the ends of the edge.
   A contact per face, edge, or vertex instead would miss the vertices,
   or count the same point twice at the boundary of a face and an edge
---------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::sphere_against_polyhedron(int ibody, int jbody,
  int itype, int jtype, double** x, double** v, double** f, double** torque,
  double** angmom, double** fnc, Scratch &s, double &evdwl, double* facc)
{
  int ifirst = dfirst[ibody];
  int iefirst = edfirst[ibody];
  int iffirst = facfirst[ibody];
  double contact_dist = rounded_radius[ibody] + rounded_radius[jbody];
  double xi1[3], xi2[3], xi3[3], ui[3], vi[3], n[3], h[3], hmin[3], d, t;
  double dmin = -1.0;
  int inside, tmp;

  // the center of the sphere is inside the core of the polyhedron, for large
  // overlaps: push it out through the nearest face

  if (facnum[ibody] > 0) {
    double sd = nearest_face(ibody, x[ibody], x[jbody], n);
    if (sd <= 0.0) {
      sd = MIN(sd, -EPSILON*EPSILON*contact_dist);
      for (int k = 0; k < 3; k++) h[k] = x[jbody][k] - sd*n[k];
      sphere_point_contact(ibody, jbody, itype, jtype, h, x, v, f, torque, angmom, fnc, s,
                           evdwl, facc, 1);
      return;
    }
  }

  // faces in front of the center of the sphere with the projection inside

  for (int nf = 0; nf < facnum[ibody]; nf++) {
    MathExtra::add3(x[ibody], discrete[ifirst+static_cast<int>(face[iffirst+nf][0])], xi1);
    MathExtra::add3(x[ibody], discrete[ifirst+static_cast<int>(face[iffirst+nf][1])], xi2);
    MathExtra::add3(x[ibody], discrete[ifirst+static_cast<int>(face[iffirst+nf][2])], xi3);
    MathExtra::sub3(xi2, xi1, ui);
    MathExtra::sub3(xi3, xi1, vi);
    MathExtra::cross3(ui, vi, n);
    MathExtra::norm3(n);
    if (opposite_sides(n, xi1, x[ibody], x[jbody]) == 0) continue;
    project_pt_plane(x[jbody], xi1, xi2, xi3, h, d, inside);
    inside_polygon(ibody, nf, x[ibody], h, nullptr, inside, tmp);
    if (!inside) continue;
    if ((dmin < 0.0) || (d < dmin)) {
      dmin = d;
      MathExtra::copy3(h, hmin);
    }
  }

  // edges, the nearest point is an end of the edge if the projection is not inside

  for (int ne = 0; ne < ednum[ibody]; ne++) {
    MathExtra::add3(x[ibody], discrete[ifirst+static_cast<int>(edge[iefirst+ne][0])], xi1);
    MathExtra::add3(x[ibody], discrete[ifirst+static_cast<int>(edge[iefirst+ne][1])], xi2);
    project_pt_line(x[jbody], xi1, xi2, h, d, t);
    if (t < 0.0) {
      MathExtra::copy3(xi1, h);
      d = sqrt(MathExtra::distsq3(x[jbody], xi1));
    } else if (t > 1.0) {
      MathExtra::copy3(xi2, h);
      d = sqrt(MathExtra::distsq3(x[jbody], xi2));
    }
    if ((dmin < 0.0) || (d < dmin)) {
      dmin = d;
      MathExtra::copy3(h, hmin);
    }
  }

  if ((dmin < 0.0) || (dmin > contact_dist + cut_inner)) return;
  sphere_point_contact(ibody, jbody, itype, jtype, hmin, x, v, f, torque, angmom, fnc, s, evdwl,
                       facc);
}

/* ----------------------------------------------------------------------
   Force between the point h of a polyhedron (ibody) and a sphere (jbody)
   inside = 1 if the center of the sphere is inside the core of the
   polyhedron, then the distance to h is negative
---------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::sphere_point_contact(int ibody, int jbody,
  int itype, int jtype, double *h, double** x, double** v, double** f, double** torque,
  double** angmom, double** fnc, Scratch &s, double &evdwl, double* facc, int inside)
{
  double delx,dely,delz,rsq,rij,R,fx,fy,fz,fpair,energy;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double rradi = rounded_radius[ibody];
  double rradj = rounded_radius[jbody];
  double contact_dist = rradi + rradj;

  delx = h[0] - x[jbody][0];
  dely = h[1] - x[jbody][1];
  delz = h[2] - x[jbody][2];
  rsq = delx*delx + dely*dely + delz*delz;
  if (rsq == 0.0) return;
  rij = inside ? -sqrt(rsq) : sqrt(rsq);
  R = rij - contact_dist;

  energy = 0;
  kernel_force(R, itype, jtype, energy, fpair);

  fx = delx*fpair/rij;
  fy = dely*fpair/rij;
  fz = delz*fpair/rij;

  if (R <= 0) { // in contact

    // damping at the contact point between the surfaces, the friction
    // force is computed once per pair of bodies from the contacts

    double nrm[3] = {delx/rij, dely/rij, delz/rij};
    double pc[3];
    contact_point(h, x[jbody], nrm, rradi, rradj, pc);
    damping_friction(ibody, jbody, pc, nrm, 0.0, 1, 0, x, v, angmom, f, torque, fnc, ibody,
                     facc);

    Contact c;
    c.ibody = ibody;
    c.jbody = jbody;
    MathExtra::copy3(h, c.xi);
    MathExtra::copy3(x[jbody], c.xj);
    c.type = 0;
    c.separation = R;
    c.r = rij;
    c.unique = 1;
    c.patch = 0;
    c.w = 1.0;
    s.contacts.push_back(c);
  }

  f[ibody][0] += fx;
  f[ibody][1] += fy;
  f[ibody][2] += fz;
  sum_torque(x[ibody], h, fx, fy, fz, torque[ibody]);

  if (newton_pair || jbody < nlocal) {
    f[jbody][0] -= fx;
    f[jbody][1] -= fy;
    f[jbody][2] -= fz;
  }

  evdwl += energy;
  facc[0] += fx; facc[1] += fy; facc[2] += fz;
}

/* ----------------------------------------------------------------------
   Return 1 if vertex ni of body ibody may interact with the other body
   of the pair, see pair_interaction(): a vertex, edge or face that does
   not reach far enough towards the other body is farther than the
   interaction range from all points of the other body
------------------------------------------------------------------------- */

int PairBodyRoundedPolyhedron::vertex_near(const Scratch &s, int ibody, int ni) const
{
  double reach_min = (ibody == s.ibody) ? s.reach_min_i : s.reach_min_j;
  return (s.reach[dfirst[ibody]+ni] >= reach_min) ? 1 : 0;
}

/* ----------------------------------------------------------------------
   Return 1 if edge ne of body ibody may interact with the other body
------------------------------------------------------------------------- */

int PairBodyRoundedPolyhedron::edge_near(const Scratch &s, int ibody, int ne) const
{
  int iefirst = edfirst[ibody];
  return (vertex_near(s, ibody, static_cast<int>(edge[iefirst+ne][0])) ||
          vertex_near(s, ibody, static_cast<int>(edge[iefirst+ne][1]))) ? 1 : 0;
}

/* ----------------------------------------------------------------------
   Return 1 if face nf of body ibody may interact with the other body
------------------------------------------------------------------------- */

int PairBodyRoundedPolyhedron::face_near(const Scratch &s, int ibody, int nf) const
{
  int iffirst = facfirst[ibody];
  for (int k = 0; k < MAX_FACE_SIZE; k++) {
    int np = static_cast<int>(face[iffirst+nf][k]);
    if (np < 0) break;
    if (vertex_near(s, ibody, np)) return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   Normal cones of the vertices and edges of body ibody: return 1 if the
   direction d from a vertex or edge towards the other body of a contact is
   in the normal cone of that feature, i.e. no incident edge of the vertex,
   or adjacent face of the edge, is closer to the other body in direction d.
   Else that edge or face is closer, and the contact is represented by its
   contacts instead: e.g. a vertex of a face in contact with a parallel face
   beyond an edge of that face would otherwise get a tilted normal towards
   the edge and push the bodies sideways.  Ties, e.g. for parallel faces,
   are accepted within CONE_TOL.
   only the point-like contacts of a vertex with an edge or a vertex are
   tested: the contacts of the vertices of a face with another face, and of
   crossing edges, are the corners of the contact region of two faces, which
   are in contact also if not nearest, when the faces are inclined
------------------------------------------------------------------------- */

int PairBodyRoundedPolyhedron::vertex_cone(int ibody, int nv, const double *d,
                                           const double *ul) const
{
  int ifirst = dfirst[ibody];
  int iefirst = edfirst[ibody];
  double dlen = MathExtra::len3(d);
  double e[3];

  for (int ne = 0; ne < ednum[ibody]; ne++) {
    int na = static_cast<int>(edge[iefirst+ne][0]);
    int nb = static_cast<int>(edge[iefirst+ne][1]);
    int nw;
    if (na == nv) nw = nb;
    else if (nb == nv) nw = na;
    else continue;
    MathExtra::sub3(discrete[ifirst+nw], discrete[ifirst+nv], e);
    double elen = MathExtra::len3(e);

    // an edge nearly parallel to the direction ul of an edge of the other body
    // touches it along a line, whose end is the vertex, see edge_edge_weight()

    if (ul && (fabs(MathExtra::dot3(e, ul)) > COS_LINE_ANGLE * elen * MathExtra::len3(ul)))
      continue;
    if (MathExtra::dot3(d, e) > CONE_TOL * dlen * elen) return 0;
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

int PairBodyRoundedPolyhedron::edge_cone(int ibody, int ne, const double *d) const
{
  int ifirst = dfirst[ibody];
  int iefirst = edfirst[ibody];
  int iffirst = facfirst[ibody];
  int na = static_cast<int>(edge[iefirst+ne][0]);
  int nb = static_cast<int>(edge[iefirst+ne][1]);
  double dlen = MathExtra::len3(d);
  double u[3], tf[3];
  MathExtra::sub3(discrete[ifirst+nb], discrete[ifirst+na], u);
  double uu = MathExtra::dot3(u, u);
  if (uu == 0.0) return 1;

  // the faces adjacent to the edge contain both of its end points,
  // tf is the direction from the edge into the face, perpendicular to the edge

  for (int nf = 0; nf < facnum[ibody]; nf++) {
    int hasa = 0, hasb = 0, nvf = 0;
    double xc[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < MAX_FACE_SIZE; k++) {
      int np = static_cast<int>(face[iffirst+nf][k]);
      if (np < 0) break;
      if (np == na) hasa = 1;
      if (np == nb) hasb = 1;
      MathExtra::add3(xc, discrete[ifirst+np], xc);
      nvf++;
    }
    if (!hasa || !hasb) continue;
    MathExtra::scale3(1.0/nvf, xc);
    MathExtra::sub3(xc, discrete[ifirst+na], tf);
    double s = MathExtra::dot3(tf, u) / uu;
    for (int k = 0; k < 3; k++) tf[k] -= s*u[k];
    if (MathExtra::dot3(d, tf) > CONE_TOL * dlen * MathExtra::len3(tf)) return 0;
  }
  return 1;
}

/* ----------------------------------------------------------------------
   Weights of the contacts, which scale their forces and energies.
   Two nearly parallel faces touch over the region where they overlap, which
   the contacts represent by its corners: the vertices of either face inside
   the other face, and the crossings of their edges.  The number of corners
   changes, e.g. when a vertex crosses an edge of the other face and two
   crossings replace it, or when one face is rotated against the other,
   which adds crossings near the middle of the edges, where the region
   hardly turns.  So these regions are found as polygons, and each corner is
   weighted by the exterior angle of the polygon there, divided by kappa =
   2 pi / (the larger number of vertices of the two faces), see
   face_face_patches(): since the exterior angles of a polygon add up to
   2 pi, the weights of the corners always add up to the same total and
   change continuously, a corner where the region hardly turns has a small
   weight, and the corners of a face inside a face of the same shape have
   a weight of one.  The regions apply fully up to half of PATCH_ANGLE
   between the faces, see patch_factor(), and fade out as they get narrower
   than the contact distance, or as the faces penetrate each other by more
   than the contact distance, see face_face_patches().  The other contacts
   are part of such a region to the same extent, and only the remaining
   weight applies to them, see feature_patch_factor().
   Similarly, two nearly parallel edges touch along the overlap of the edges,
   which the vertex contacts at its ends represent, while a crossing of the
   edges represents it once the edges are clearly inclined.  So the weight of
   a crossing grows from zero for parallel edges to one at an angle of
   LINE_ANGLE, see edge_edge_weight(), and the contacts of the vertices of
   the crossing edges keep the remaining weight, see vertex_edges_weight().
------------------------------------------------------------------------- */

static double smoothstep(double x)
{
  if (x <= 0.0) return 0.0;
  if (x >= 1.0) return 1.0;
  return x * x * (3.0 - 2.0 * x);
}

/* ----------------------------------------------------------------------
   factor for the weights of the corners of a region of overlap of two
   faces with the outward normals n1 and n2, from the angle between the
   faces: one up to half of PATCH_ANGLE, zero from PATCH_ANGLE on
------------------------------------------------------------------------- */

double PairBodyRoundedPolyhedron::patch_factor(const double *n1, const double *n2) const
{
  double c = -MathExtra::dot3(n1, n2);
  if (c <= COS_PATCH_ANGLE) return 0.0;
  if (c >= COS_HALF_PATCH_ANGLE) return 1.0;
  double phi = acos(MIN(1.0, c));
  return 1.0 - smoothstep(2.0 * phi / PATCH_ANGLE - 1.0);
}

/* ----------------------------------------------------------------------
   vertex of the polygon of a region of overlap of two faces i and j, see
   face_face_patches(), with coordinates u, v, its origin, and the line of
   the edge from it to the next vertex:
     type = 0: vertex a of face j, 1: vertex a of face i, 2: crossing of edge
            a of face j with edge b of face i, 3: on the boundary of the range
     src = 0: edge srci of face j, 1: edge srci of face i, 2: boundary of the range
   edge k of a face runs from its vertex k to vertex k+1
------------------------------------------------------------------------- */

namespace {
constexpr int MAXP = 2*MAX_FACE_SIZE + 2;   // vertices of a clipped polygon

struct PatchVertex {
  double u, v;
  int type, a, b;
  int src, srci;
};

/* ----------------------------------------------------------------------
   clip the convex polygon poly with np vertices by the half-plane
   hu*u + hv*v + h0 >= 0 along edge line of face i (line >= 0) or along the
   boundary of the range (line < 0), using the storage work, merge vertices
   closer than tol, and return the number of vertices of the clipped polygon
   nvi = number of vertices of face i
------------------------------------------------------------------------- */

double closest_segments(const double *p1, const double *q1, const double *p2, const double *q2,
                        double *c1, double *c2);

int clip_patch(PatchVertex *poly, int np, double hu, double hv, double h0, int line, int nvi,
               PatchVertex *work, double tol)
{
  int nw = 0;
  for (int m = 0; m < np; m++) {
    PatchVertex &pp = poly[m];
    PatchVertex &pq = poly[(m+1) % np];
    double sp = hu*pp.u + hv*pp.v + h0;
    double sq = hu*pq.u + hv*pq.v + h0;
    if (sp >= 0.0) {
      if (nw >= MAXP) break;
      work[nw++] = pp;
    }
    if ((sp >= 0.0) != (sq >= 0.0)) {

      // clipping a convex polygon by a half-plane adds at most one vertex,
      // guard against round-off with nearly collinear vertices

      if (nw >= MAXP) break;
      PatchVertex &pn = work[nw++];
      double t = sp / (sp - sq);
      pn.u = pp.u + t*(pq.u - pp.u);
      pn.v = pp.v + t*(pq.v - pp.v);

      // intersection of the line of the edge from pp with the clipping line

      pn.type = 3;
      pn.a = pn.b = -1;
      if (line >= 0) {
        if (pp.src == 0) {
          pn.type = 2;
          pn.a = pp.srci;
          pn.b = line;
        } else if (pp.src == 1) {
          if (line == (pp.srci+1) % nvi) {
            pn.type = 1;
            pn.a = line;
          } else if (pp.srci == (line+1) % nvi) {
            pn.type = 1;
            pn.a = pp.srci;
          }
        }
      }

      // the edge from the intersection runs along the clipping line when
      // leaving the half-plane, else along the edge from pp

      if (sp >= 0.0) {
        pn.src = (line >= 0) ? 1 : 2;
        pn.srci = line;
      } else {
        pn.src = pp.src;
        pn.srci = pp.srci;
      }
    }
  }

  // merge vertices that coincide, e.g. a vertex on an edge, keeping the edge
  // from the last of them

  int n = 0;
  for (int m = 0; m < nw; m++) {
    if ((n > 0) && (fabs(work[m].u - poly[n-1].u) < tol) && (fabs(work[m].v - poly[n-1].v) < tol)) {
      poly[n-1].src = work[m].src;
      poly[n-1].srci = work[m].srci;
      continue;
    }
    poly[n++] = work[m];
  }
  if ((n > 1) && (fabs(poly[n-1].u - poly[0].u) < tol) && (fabs(poly[n-1].v - poly[0].v) < tol))
    n--;
  return n;
}

/* ----------------------------------------------------------------------
   closest points c1 on the segment p1-q1 and c2 on the segment p2-q2,
   return their distance, see e.g. Ericson, Real-Time Collision Detection
------------------------------------------------------------------------- */

double closest_segments(const double *p1, const double *q1, const double *p2, const double *q2,
                        double *c1, double *c2)
{
  double d1[3], d2[3], r[3];
  MathExtra::sub3(q1, p1, d1);
  MathExtra::sub3(q2, p2, d2);
  MathExtra::sub3(p1, p2, r);
  double a = MathExtra::dot3(d1, d1);
  double e = MathExtra::dot3(d2, d2);
  double f = MathExtra::dot3(d2, r);
  double s, t;
  if ((a == 0.0) && (e == 0.0)) {
    s = t = 0.0;
  } else if (a == 0.0) {
    s = 0.0;
    t = MAX(0.0, MIN(1.0, f / e));
  } else {
    double cc = MathExtra::dot3(d1, r);
    if (e == 0.0) {
      t = 0.0;
      s = MAX(0.0, MIN(1.0, -cc / a));
    } else {
      double b = MathExtra::dot3(d1, d2);
      double denom = a*e - b*b;
      s = (denom > 0.0) ? MAX(0.0, MIN(1.0, (b*f - cc*e) / denom)) : 0.0;
      t = (b*s + f) / e;
      if (t < 0.0) {
        t = 0.0;
        s = MAX(0.0, MIN(1.0, -cc / a));
      } else if (t > 1.0) {
        t = 1.0;
        s = MAX(0.0, MIN(1.0, (b - cc) / a));
      }
    }
  }
  for (int k = 0; k < 3; k++) {
    c1[k] = p1[k] + s*d1[k];
    c2[k] = p2[k] + t*d2[k];
  }
  return sqrt(MathExtra::distsq3(c1, c2));
}
}    // namespace

/* ----------------------------------------------------------------------
   Corners of the regions of overlap of the nearly parallel faces of body
   ibody and body jbody, see the weights above: both faces are projected onto
   the plane halfway between them, with the unit normal nm from body i to
   body j, and the face of body j is clipped by the face of body i and by the
   range of interaction.  Each corner of the clipped polygon interacts with
   the weight patch_factor() times its exterior angle over kappa, like the
   vertex of one face at the corner with the other face, along the normal of
   that face, or like the two edges crossing at the corner, along their
   common normal, so that the forces derive from the energy of the corner.
   Since all corners follow from the same polygon, a vertex of one face that
   crosses an edge of the other face passes its exterior angle exactly on to
   the two crossings that replace it, also when the faces are not exactly
   parallel, and the corners on the boundary of the range do not interact.
   The factors of the pairs of faces are stored in s.patch for the weights
   of the other contacts, see feature_patch_factor(), which follow.
   the total force on body i is accumulated to facc
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::face_face_patches(int ibody, int jbody, int itype, int jtype,
  double** x, double** v, double** f, double** torque, double** angmom, double** fnc,
  Scratch &s, double &evdwl, double* facc)
{
  int ifirst = dfirst[ibody];
  int jfirst = dfirst[jbody];
  int iffirst = facfirst[ibody];
  int jffirst = facfirst[jbody];
  double contact_dist = rounded_radius[ibody] + rounded_radius[jbody];
  double tol = EPSILON*EPSILON*EPSILON*MAX(enclosing_radius[ibody], enclosing_radius[jbody]);
  double energy = 0.0;
  double ni[3], nj[3], nm[3], e1[3], e2[3], o[3], q[3], pa[3], pb[3], xa[3], xb[3];
  double vi[MAX_FACE_SIZE][3], vj[MAX_FACE_SIZE][3], pi2[MAX_FACE_SIZE][2];
  PatchVertex poly[MAXP], work[MAXP];

  // factors of the pairs of faces, see feature_patch_factor(), which are
  // stored once the first pair of faces is nearly parallel, since that is rare
  // the faces near the other body are determined only for such pairs, with
  // -1 for not yet determined

  int nfj = facnum[jbody];
  s.patch_any = 0;
  s.near_j.assign(nfj, -1);

  for (int fi = 0; fi < facnum[ibody]; fi++) {
    face_normal(ibody, fi, ni);
    int near_i = -1;
    for (int fj = 0; fj < facnum[jbody]; fj++) {
      face_normal(jbody, fj, nj);
      if (-MathExtra::dot3(ni, nj) <= COS_PATCH_ANGLE) continue;
      double sp = patch_factor(ni, nj);
      if (sp <= 0.0) continue;
      if (near_i < 0) near_i = face_near(s, ibody, fi);
      if (!near_i) break;
      if (s.near_j[fj] < 0) s.near_j[fj] = face_near(s, jbody, fj);
      if (!s.near_j[fj]) continue;

      // plane halfway between the faces with the basis e1, e2, and the
      // center o of face i as origin

      MathExtra::sub3(ni, nj, nm);
      MathExtra::norm3(nm);
      double ax[3] = {1.0, 0.0, 0.0};
      if (fabs(nm[0]) > 0.9) {
        ax[0] = 0.0;
        ax[1] = 1.0;
      }
      MathExtra::cross3(nm, ax, e1);
      MathExtra::norm3(e1);
      MathExtra::cross3(nm, e1, e2);

      // vertices of both faces in space, and their projections, oriented
      // counterclockwise

      int nvi = face_size(ibody, fi);
      int nvj = face_size(jbody, fj);
      o[0] = o[1] = o[2] = 0.0;
      for (int k = 0; k < nvi; k++) {
        MathExtra::add3(x[ibody], discrete[ifirst+static_cast<int>(face[iffirst+fi][k])], vi[k]);
        MathExtra::add3(o, vi[k], o);
      }
      MathExtra::scale3(1.0/nvi, o);
      for (int k = 0; k < nvj; k++)
        MathExtra::add3(x[jbody], discrete[jfirst+static_cast<int>(face[jffirst+fj][k])], vj[k]);
      MathExtra::copy3(vi[0], xa);
      MathExtra::copy3(vj[0], xb);

      double areai = 0.0, areaj = 0.0;
      for (int k = 0; k < nvi; k++) {
        MathExtra::sub3(vi[k], o, q);
        pi2[k][0] = MathExtra::dot3(q, e1);
        pi2[k][1] = MathExtra::dot3(q, e2);
      }
      for (int k = 0; k < nvj; k++) {
        MathExtra::sub3(vj[k], o, q);
        poly[k].u = MathExtra::dot3(q, e1);
        poly[k].v = MathExtra::dot3(q, e2);
      }
      for (int k = 0; k < nvi; k++) {
        int kn = (k+1) % nvi;
        areai += pi2[k][0]*pi2[kn][1] - pi2[kn][0]*pi2[k][1];
      }
      for (int k = 0; k < nvj; k++) {
        int kn = (k+1) % nvj;
        areaj += poly[k].u*poly[kn].v - poly[kn].u*poly[k].v;
      }
      if (areai < 0.0)
        for (int k = 0; k < nvi/2; k++) {
          std::swap(pi2[k][0], pi2[nvi-1-k][0]);
          std::swap(pi2[k][1], pi2[nvi-1-k][1]);
          for (int m = 0; m < 3; m++) std::swap(vi[k][m], vi[nvi-1-k][m]);
        }
      if (areaj < 0.0)
        for (int k = 0; k < nvj/2; k++) {
          std::swap(poly[k], poly[nvj-1-k]);
          for (int m = 0; m < 3; m++) std::swap(vj[k][m], vj[nvj-1-k][m]);
        }
      for (int k = 0; k < nvj; k++) {
        poly[k].type = 0;
        poly[k].a = k;
        poly[k].b = -1;
        poly[k].src = 0;
        poly[k].srci = k;
      }

      // clip the face of body j by the edges of the face of body i, and by the
      // range of interaction, where the gap between the planes of the faces
      // along nm is linear in the coordinates on the plane

      double cni = MathExtra::dot3(nm, ni);
      double cnj = MathExtra::dot3(nm, nj);
      double g[3];
      MathExtra::sub3(xb, o, q);
      g[0] = MathExtra::dot3(q, nj) / cnj;
      MathExtra::sub3(xa, o, q);
      g[0] -= MathExtra::dot3(q, ni) / cni;
      g[1] = MathExtra::dot3(e1, ni) / cni - MathExtra::dot3(e1, nj) / cnj;
      g[2] = MathExtra::dot3(e2, ni) / cni - MathExtra::dot3(e2, nj) / cnj;

      int np = nvj;
      for (int k = 0; (k <= nvi) && (np > 0); k++) {

        // half-plane hu*u + hv*v + h0 >= 0

        double hu, hv, h0;
        if (k < nvi) {
          double *pa2 = pi2[k];
          double *pb2 = pi2[(k+1) % nvi];
          hu = -(pb2[1] - pa2[1]);
          hv = pb2[0] - pa2[0];
          h0 = -(hu*pa2[0] + hv*pa2[1]);
        } else {
          hu = -g[1];
          hv = -g[2];
          h0 = contact_dist + cut_inner - g[0];
        }
        np = clip_patch(poly, np, hu, hv, h0, (k < nvi) ? k : -1, nvi, work, tol);
      }
      if (np < 3) continue;

      // the regions of overlap fade out as they get narrower than the
      // contact distance, e.g. for a body sliding off the face of another,
      // or two parallel faces that meet at their edges only, and as the faces
      // penetrate each other by more than the contact distance, which the
      // contacts of their vertices and edges then represent

      double gmin = BIG;
      for (int k = 0; k < np; k++) gmin = MIN(gmin, g[0] + g[1]*poly[k].u + g[2]*poly[k].v);
      sp *= smoothstep(1.0 + gmin / contact_dist);

      double width = BIG;
      for (int k = 0; k < np; k++) {
        PatchVertex &pk = poly[k];
        PatchVertex &pn = poly[(k+1) % np];
        double ex = pn.u - pk.u, ey = pn.v - pk.v;
        double elen = sqrt(ex*ex + ey*ey);
        if (elen == 0.0) continue;
        double dmax = 0.0;
        for (int m = 0; m < np; m++)
          dmax = MAX(dmax, (ex*(poly[m].v - pk.v) - ey*(poly[m].u - pk.u)) / elen);
        width = MIN(width, dmax);
      }
      sp *= smoothstep(width / contact_dist);
      if (sp <= 0.0) continue;
      if (!s.patch_any) {
        s.patch.assign(static_cast<std::size_t>(facnum[ibody])*nfj, 0.0);
        s.patch_any = 1;
      }
      s.patch[fi*nfj+fj] = sp;

      // corners of the region of overlap, with the geometry of the contact of
      // their vertex with the other face, or of their two crossing edges, whose
      // forces derive from their energy

      double kappa = MY_2PI / MAX(nvi, nvj);
      for (int k = 0; k < np; k++) {
        PatchVertex &pp = poly[(k+np-1) % np];
        PatchVertex &pc = poly[k];
        PatchVertex &pn = poly[(k+1) % np];
        if (pc.type == 3) continue;
        double d1x = pc.u - pp.u, d1y = pc.v - pp.v;
        double d2x = pn.u - pc.u, d2y = pn.v - pc.v;
        double ext = atan2(d1x*d2y - d1y*d2x, d1x*d2x + d1y*d2y);
        if (ext <= 0.0) continue;
        double w = sp * ext / kappa;

        double r;
        if (pc.type == 0) {
          MathExtra::sub3(vj[pc.a], xa, q);
          r = MathExtra::dot3(q, ni);
          MathExtra::copy3(vj[pc.a], pb);
          for (int m = 0; m < 3; m++) pa[m] = pb[m] - r*ni[m];
        } else if (pc.type == 1) {
          MathExtra::sub3(vi[pc.a], xb, q);
          r = MathExtra::dot3(q, nj);
          MathExtra::copy3(vi[pc.a], pa);
          for (int m = 0; m < 3; m++) pb[m] = pa[m] - r*nj[m];
        } else {

          // the nearest points of the edges themselves, since those of nearly
          // parallel lines through the edges may be far away from the edges

          r = closest_segments(vj[pc.a], vj[(pc.a+1) % nvj], vi[pc.b], vi[(pc.b+1) % nvi],
                               pb, pa);
          MathExtra::sub3(pb, pa, q);
          if (MathExtra::dot3(q, nm) < 0.0) r = -r;
        }
        if (r > contact_dist + cut_inner) continue;

        // the corner lies on the other face, or the edges touch: the force
        // is along the normal of that face, or of the plane halfway between
        // the faces, also for a vanishing distance

        if (fabs(r) < tol) {
          if (r == 0.0) r = EPSILON*tol;
          if (pc.type == 0) {
            for (int m = 0; m < 3; m++) pa[m] = pb[m] - r*ni[m];
          } else if (pc.type == 1) {
            for (int m = 0; m < 3; m++) pb[m] = pa[m] - r*nj[m];
          } else {
            for (int m = 0; m < 3; m++) pb[m] = pa[m] + r*nm[m];
          }
        }

        pair_force_and_torque(ibody, jbody, pa, pb, r, contact_dist, itype, jtype,
                              x, v, f, torque, angmom, fnc, 1, energy, facc, w);

        if (r <= contact_dist) {
          Contact c;
          c.ibody = ibody;
          c.jbody = jbody;
          MathExtra::copy3(pa, c.xi);
          MathExtra::copy3(pb, c.xj);
          c.type = 2;
          c.separation = r - contact_dist;
          c.r = r;
          c.unique = 1;
          c.w = w;
          c.patch = 1;
          s.contacts.push_back(c);
        }
      }
    }
  }

  evdwl += energy;
}

/* ----------------------------------------------------------------------
   Return 1 if face nf of body ibody contains vertex nv (if nv >= 0) or
   edge ne (if ne >= 0)
------------------------------------------------------------------------- */

int PairBodyRoundedPolyhedron::face_has(int ibody, int nf, int nv, int ne) const
{
  int iefirst = edfirst[ibody];
  int iffirst = facfirst[ibody];
  int na = nv, nb = nv;
  if (ne >= 0) {
    na = static_cast<int>(edge[iefirst+ne][0]);
    nb = static_cast<int>(edge[iefirst+ne][1]);
  }
  int hasa = 0, hasb = 0;
  for (int k = 0; k < face_size(ibody, nf); k++) {
    int np = static_cast<int>(face[iffirst+nf][k]);
    if (np == na) hasa = 1;
    if (np == nb) hasb = 1;
  }
  return (hasa && hasb) ? 1 : 0;
}

/* ----------------------------------------------------------------------
   largest factor of a region of overlap, see face_face_patches(), of a face
   of body ibody at vertex iv or edge ie (the other one is -1) and a face of
   body jbody at vertex jv, edge je, or the face jf (the others are -1): the
   contact of these two features is part of that region to this extent
------------------------------------------------------------------------- */

double PairBodyRoundedPolyhedron::feature_patch_factor(const Scratch &s, int ibody, int iv, int ie,
                                                       int jbody, int jv, int je, int jf) const
{
  if (!s.patch_any || s.patch.empty()) return 0.0;
  int nfj = facnum[s.jbody];
  double smax = 0.0;
  for (int mf = 0; mf < facnum[ibody]; mf++) {
    if (!face_has(ibody, mf, iv, ie)) continue;
    for (int nf = 0; nf < facnum[jbody]; nf++) {
      if (jf >= 0) {
        if (nf != jf) continue;
      } else if (!face_has(jbody, nf, jv, je)) continue;
      double sp = (ibody == s.ibody) ? s.patch[mf*nfj+nf] : s.patch[nf*nfj+mf];
      smax = MAX(smax, sp);
    }
  }
  return smax;
}

/* ----------------------------------------------------------------------
   weight of the contact of vertex nv of body ibody with face nf of body
   jbody, which is part of a region of overlap of two nearly parallel faces
   to the extent of feature_patch_factor(), see face_face_patches()
------------------------------------------------------------------------- */

double PairBodyRoundedPolyhedron::vertex_face_weight(const Scratch &s, int ibody, int nv, int jbody,
                                                     int nf) const
{
  return 1.0 - feature_patch_factor(s, ibody, nv, -1, jbody, -1, -1, nf);
}

/* ----------------------------------------------------------------------
   weight of the contact at the crossing of edge ei of body ibody and edge
   ej of body jbody: the crossing of two edges, whose weight vanishes for
   parallel edges, see above, unless it is part of a region of overlap of
   two nearly parallel faces next to the edges, see face_face_patches()
------------------------------------------------------------------------- */

double PairBodyRoundedPolyhedron::edge_edge_weight(const Scratch &s, int ibody, int ei, int jbody,
                                                   int ej) const
{
  int ifirst = dfirst[ibody];
  int jfirst = dfirst[jbody];
  int iefirst = edfirst[ibody];
  int jefirst = edfirst[jbody];
  double ui[3], uj[3];
  MathExtra::sub3(discrete[ifirst+static_cast<int>(edge[iefirst+ei][1])],
                  discrete[ifirst+static_cast<int>(edge[iefirst+ei][0])], ui);
  MathExtra::sub3(discrete[jfirst+static_cast<int>(edge[jefirst+ej][1])],
                  discrete[jfirst+static_cast<int>(edge[jefirst+ej][0])], uj);
  double c = fabs(MathExtra::dot3(ui, uj)) / (MathExtra::len3(ui) * MathExtra::len3(uj));
  double wline = (c <= COS_LINE_ANGLE) ? 1.0 : smoothstep(acos(MIN(1.0, c)) / LINE_ANGLE);
  return (1.0 - feature_patch_factor(s, ibody, -1, ei, jbody, -1, ej, -1)) * wline;
}

/* ----------------------------------------------------------------------
   Nearest points hi and hj, at a distance r, of edge ei of body ibody and
   edge ej of body jbody, return 1 if the edges interact as edges, i.e. the
   nearest points are inside the edges and within the interaction range.
   crossed = 1 if the edges have crossed each other, when the nearest point
   of the edge of body j is inside body i.
   the normal cones are not tested, see edge_cone(): two crossing edges of
   faces in contact are a corner of the contact region, also when the faces
   are inclined to each other, and a face is then closer than the edge on
   one side of the crossing
------------------------------------------------------------------------- */

int PairBodyRoundedPolyhedron::edge_edge_nearest(int ibody, int ei, int jbody, int ej,
                                                 double *hi, double *hj, double &r,
                                                 int &crossed) const
{
  double **x = atom->x;
  double xi1[3], xi2[3], xj1[3], xj2[3], ti, tj;
  int ifirst = dfirst[ibody];
  int jfirst = dfirst[jbody];
  int iefirst = edfirst[ibody];
  int jefirst = edfirst[jbody];
  MathExtra::add3(x[ibody], discrete[ifirst+static_cast<int>(edge[iefirst+ei][0])], xi1);
  MathExtra::add3(x[ibody], discrete[ifirst+static_cast<int>(edge[iefirst+ei][1])], xi2);
  MathExtra::add3(x[jbody], discrete[jfirst+static_cast<int>(edge[jefirst+ej][0])], xj1);
  MathExtra::add3(x[jbody], discrete[jfirst+static_cast<int>(edge[jefirst+ej][1])], xj2);
  distance_bt_edges(xj1, xj2, xi1, xi2, hj, hi, tj, ti, r);

  crossed = 0;
  double contact_dist = rounded_radius[ibody] + rounded_radius[jbody];
  double rmin = MIN(rounded_radius[ibody], rounded_radius[jbody]);
  if ((ti < 0) || (ti > 1) || (tj < 0) || (tj > 1) || (r >= contact_dist + cut_inner)) return 0;

  // unit common normal n of the two edges, pointing from body i to body j:
  // outward of the faces of body i at its edge and inward of those of body j,
  // or else along the line between the centers.  hj - hi is along n unless
  // the edges are parallel.  The sign of the separation along n tells
  // whether the edges have crossed, also when they touch or when the closest
  // points lie on the plane of a face: nearest_face() alone cannot decide
  // this, e.g. for an edge that pierces a face of body i next to its edge

  double ui[3], uj[3], n[3], d[3];
  MathExtra::sub3(xi2, xi1, ui);
  MathExtra::sub3(xj2, xj1, uj);
  MathExtra::cross3(ui, uj, n);
  MathExtra::sub3(hj, hi, d);
  double nlen = MathExtra::len3(n);

  // parallel edges: the direction between the closest points

  if (nlen <= PARALLEL_TOL * sqrt(MathExtra::lensq3(ui) * MathExtra::lensq3(uj))) {
    if (r < EPSILON*rmin) return 0;
    double nc[3];
    if (nearest_face(ibody, x[ibody], hj, nc) < 0.0) crossed = 1;
    return 1;
  }
  MathExtra::scale3(1.0/nlen, n);

  double score = 0.0, nf[3];
  for (int mf = 0; mf < facnum[ibody]; mf++) {
    if (!face_has(ibody, mf, -1, ei)) continue;
    face_normal(ibody, mf, nf);
    score += MathExtra::dot3(n, nf);
  }
  for (int mf = 0; mf < facnum[jbody]; mf++) {
    if (!face_has(jbody, mf, -1, ej)) continue;
    face_normal(jbody, mf, nf);
    score -= MathExtra::dot3(n, nf);
  }
  if (score == 0.0) {
    double c[3];
    MathExtra::sub3(x[jbody], x[ibody], c);
    score = MathExtra::dot3(n, c);
  }
  if (score < 0.0) MathExtra::negate3(n);

  // signed separation along n, which is not zero, so that the direction
  // of the force is n also when the edges touch

  double sep = MathExtra::dot3(d, n);
  if (sep == 0.0) sep = EPSILON*EPSILON*rmin;
  for (int k = 0; k < 3; k++) hj[k] = hi[k] + sep*n[k];
  r = fabs(sep);
  if (sep < 0.0) crossed = 1;
  return 1;
}

/* ----------------------------------------------------------------------
   Return 1 if edge ei of body ibody and edge ej of body jbody interact
   as edges, see edge_edge_nearest() and interaction_edge_to_edge()
------------------------------------------------------------------------- */

int PairBodyRoundedPolyhedron::edges_interact(int ibody, int ei, int jbody, int ej) const
{
  double hi[3], hj[3], r;
  int crossed;
  return edge_edge_nearest(ibody, ei, jbody, ej, hi, hj, r, crossed);
}

/* ----------------------------------------------------------------------
   Return the largest weight of the edge-edge interactions of an edge of body
   ibody that ends at vertex ni with edge ej of body jbody, or with an edge
   of body jbody that ends at vertex nj if ej < 0, zero if there is none:
   the contact of the vertex is represented by these edge-edge interactions
   to the extent of their weight, see edge_edge_weight(), e.g. fully when an
   edge crosses the other edge near the vertex, and not at all for a vertex
   at the end of two parallel edges, whose crossing has a weight of zero
------------------------------------------------------------------------- */

double PairBodyRoundedPolyhedron::vertex_edges_weight(const Scratch &s, int ibody, int ni,
                                                      int jbody, int ej, int nj) const
{
  int iefirst = edfirst[ibody];
  int jefirst = edfirst[jbody];
  double wmax = 0.0;
  for (int ei = 0; ei < ednum[ibody]; ei++) {
    if ((static_cast<int>(edge[iefirst+ei][0]) != ni) &&
        (static_cast<int>(edge[iefirst+ei][1]) != ni)) continue;
    for (int e = 0; e < ednum[jbody]; e++) {
      if (ej >= 0) {
        if (e != ej) continue;
      } else if ((static_cast<int>(edge[jefirst+e][0]) != nj) &&
                 (static_cast<int>(edge[jefirst+e][1]) != nj)) continue;

      // the weight is cheaper than the test whether the edges interact, and
      // no weight is larger than one

      double w = edge_edge_weight(s, ibody, ei, jbody, e);
      if ((w > wmax) && edges_interact(ibody, ei, jbody, e)) {
        wmax = w;
        if (wmax >= 1.0) return 1.0;
      }
    }
  }
  return wmax;
}

/* ----------------------------------------------------------------------
   Interaction between the vertices of body i and the edges of body j,
   see step 4 in Fig. 7, Wang et al.: a vertex of body i that has no
   interaction with a face of body j interacts with the nearest edge of
   body j whose nearest point to the vertex is inside the edge
   the vertices that interact are flagged in the scratch space
   the total force on body i is accumulated to facc
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::vertex_against_edge(int ibody, int jbody,
  int itype, int jtype, double** x, double** v, double** f, double** torque,
  double** angmom, double** fnc, Scratch &s, double &evdwl, double* facc)
{
  std::vector<int> &vertex_done = s.vertex_done;
  std::vector<Contact> &contacts = s.contacts;

  int ifirst = dfirst[ibody];
  int jfirst = dfirst[jbody];
  int jefirst = edfirst[jbody];
  double rradi = rounded_radius[ibody];
  double rradj = rounded_radius[jbody];
  double eradj = enclosing_radius[jbody];
  double contact_dist = rradi + rradj;
  double energy = 0.0;
  double xpi[3], xj1[3], xj2[3], u[3], w[3], h[3], hmin[3], n[3];

  for (int ni = 0; ni < dnum[ibody]; ni++) {
    if (vertex_done[ifirst+ni]) continue;
    if (!vertex_near(s, ibody, ni)) continue;
    MathExtra::add3(x[ibody], discrete[ifirst+ni], xpi);

    double dist = sqrt(MathExtra::distsq3(xpi, x[jbody]));
    if (dist > eradj + rradj + rradi + cut_inner) continue;

    // only an edge within the interaction range can be the nearest edge
    // that interacts, so test the more expensive conditions only for those

    double dmin = -1.0, wmin = 0.0;
    for (int ne = 0; ne < ednum[jbody]; ne++) {
      if (!edge_near(s, jbody, ne)) continue;
      MathExtra::add3(x[jbody], discrete[jfirst+static_cast<int>(edge[jefirst+ne][0])], xj1);
      MathExtra::add3(x[jbody], discrete[jfirst+static_cast<int>(edge[jefirst+ne][1])], xj2);
      MathExtra::sub3(xj2, xj1, u);
      MathExtra::sub3(xpi, xj1, w);
      double uu = MathExtra::dot3(u, u);
      if (uu == 0.0) continue;
      double t = MathExtra::dot3(w, u) / uu;
      if ((t <= 0.0) || (t >= 1.0)) continue;
      h[0] = xj1[0] + t*u[0];
      h[1] = xj1[1] + t*u[1];
      h[2] = xj1[2] + t*u[2];
      double d = sqrt(MathExtra::distsq3(xpi, h));
      if (d > contact_dist + cut_inner) continue;

      // the vertex and the edge must be the nearest features of their bodies
      // to each other, see vertex_cone(), a vertex inside body j is skipped below

      double dv[3], dm[3];
      MathExtra::sub3(h, xpi, dv);
      for (int k = 0; k < 3; k++) dm[k] = -dv[k];
      if (!vertex_cone(ibody, ni, dv, u) || !edge_cone(jbody, ne, dm)) continue;

      // an edge ending at the vertex interacts with this edge as an edge,
      // which represents the contact until its nearest point reaches the vertex

      if ((dmin < 0.0) || (d < dmin)) {
        double we = (1.0 - vertex_edges_weight(s, ibody, ni, jbody, ne, -1)) *
          (1.0 - feature_patch_factor(s, ibody, ni, -1, jbody, -1, ne, -1));
        if (we <= 0.0) continue;
        dmin = d;
        wmin = we;
        MathExtra::copy3(h, hmin);
      }
    }

    if (dmin <= 0.0) continue;

    // a vertex inside body j is handled by the vertex-face interactions

    if (nearest_face(jbody, x[jbody], xpi, n) < 0.0) continue;

    pair_force_and_torque(ibody, jbody, xpi, hmin, dmin, contact_dist, itype, jtype,
                          x, v, f, torque, angmom, fnc, 1, energy, facc, wmin);

    if (dmin <= contact_dist) {
      Contact c;
      c.ibody = ibody;
      c.jbody = jbody;
      MathExtra::copy3(xpi, c.xi);
      MathExtra::copy3(hmin, c.xj);
      c.type = 0;
      c.separation = dmin - contact_dist;
      c.r = dmin;
      c.unique = 1;
      c.patch = 0;
      c.w = wmin;
      contacts.push_back(c);
    }
    vertex_done[ifirst+ni] = 1;
  }

  evdwl += energy;
}

/* ----------------------------------------------------------------------
   Interaction between the vertices of body i and those of body j, see
   step 5 in Fig. 7, Wang et al.: a vertex of body i that has no interaction
   with a face or an edge of body j interacts with the nearest vertex of
   body j that has no such interaction with body i either
   each vertex interacts with at most one vertex of the other body
   the total force on body i is accumulated to facc
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::vertex_against_vertex(int ibody, int jbody,
  int itype, int jtype, double** x, double** v, double** f, double** torque,
  double** angmom, double** fnc, Scratch &s, double &evdwl, double* facc)
{
  std::vector<int> &vertex_done = s.vertex_done;
  std::vector<Contact> &contacts = s.contacts;

  int ifirst = dfirst[ibody];
  int jfirst = dfirst[jbody];
  double contact_dist = rounded_radius[ibody] + rounded_radius[jbody];
  double energy = 0.0;
  double xpi[3], xpj[3], xmin[3];

  for (int ni = 0; ni < dnum[ibody]; ni++) {
    if (vertex_done[ifirst+ni]) continue;
    if (!vertex_near(s, ibody, ni)) continue;
    MathExtra::add3(x[ibody], discrete[ifirst+ni], xpi);

    // only a vertex within the interaction range can be the nearest vertex
    // that interacts, so test the more expensive conditions only for those

    double dmin = -1.0, wmin = 0.0;
    int nmin = -1;
    for (int nj = 0; nj < dnum[jbody]; nj++) {
      if (vertex_done[jfirst+nj]) continue;
      if (!vertex_near(s, jbody, nj)) continue;
      MathExtra::add3(x[jbody], discrete[jfirst+nj], xpj);
      double d = sqrt(MathExtra::distsq3(xpi, xpj));
      if (d > contact_dist + cut_inner) continue;
      if ((dmin < 0.0) || (d < dmin)) {

        // both vertices must be the nearest features of their bodies to each
        // other, see vertex_cone(), unless one is inside the other body

        double dv[3], dm[3], n[3];
        MathExtra::sub3(xpj, xpi, dv);
        for (int k = 0; k < 3; k++) dm[k] = -dv[k];
        if ((!vertex_cone(ibody, ni, dv) || !vertex_cone(jbody, nj, dm)) &&
            (nearest_face(jbody, x[jbody], xpi, n) >= 0.0) &&
            (nearest_face(ibody, x[ibody], xpj, n) >= 0.0)) continue;
        double we = (1.0 - vertex_edges_weight(s, ibody, ni, jbody, -1, nj)) *
          (1.0 - feature_patch_factor(s, ibody, ni, -1, jbody, nj, -1, -1));
        if (we <= 0.0) continue;
        dmin = d;
        wmin = we;
        nmin = nj;
        MathExtra::copy3(xpj, xmin);
      }
    }

    if ((nmin < 0) || (dmin <= 0.0)) continue;

    pair_force_and_torque(ibody, jbody, xpi, xmin, dmin, contact_dist, itype, jtype,
                          x, v, f, torque, angmom, fnc, 1, energy, facc, wmin);

    if (dmin <= contact_dist) {
      Contact c;
      c.ibody = ibody;
      c.jbody = jbody;
      MathExtra::copy3(xpi, c.xi);
      MathExtra::copy3(xmin, c.xj);
      c.type = 0;
      c.separation = dmin - contact_dist;
      c.r = dmin;
      c.unique = 1;
      c.patch = 0;
      c.w = wmin;
      contacts.push_back(c);
    }
    vertex_done[ifirst+ni] = 1;
    vertex_done[jfirst+nmin] = 1;
  }

  evdwl += energy;
}

/* ----------------------------------------------------------------------
   Determine the interaction mode between i's edges against j's edges

   i = atom i (body i)
   j = atom j (body j)
   x      = atoms' coordinates
   f      = atoms' forces
   torque = atoms' torques
   tag    = atoms' tags
   contacts = list of contacts, the contacts between i's edges
              and j's edges are appended
   Return:

---------------------------------------------------------------------- */

int PairBodyRoundedPolyhedron::edge_against_edge(int ibody, int jbody,
  int itype, int jtype, double** x, double** v, double** f, double** torque,
  double** angmom, double** fnc, Scratch &s, double &evdwl, double* facc)
{
  int ni,nei,nj,nej,interact;
  double rradi,rradj,energy;

  nei = ednum[ibody];
  rradi = rounded_radius[ibody];
  nej = ednum[jbody];
  rradj = rounded_radius[jbody];

  energy = 0;
  interact = EE_NONE;

  // loop through body i's edges

  for (ni = 0; ni < nei; ni++) {
    if (!edge_near(s, ibody, ni)) continue;

    for (nj = 0; nj < nej; nj++) {
      if (!edge_near(s, jbody, nj)) continue;

      // compute the distance between the edge nj to the edge ni
      #ifdef _POLYHEDRON_DEBUG
      printf("Compute interaction between edge %d of body %d "
             "with edge %d of body %d:\n",
             nj, jbody, ni, ibody);
      #endif

      interact = interaction_edge_to_edge(ibody, ni, x[ibody], rradi,
                                          jbody, nj, x[jbody], rradj,
                                          itype, jtype, cut_inner,
                                          v, f, torque, angmom, fnc, s,
                                          energy, facc);
    }

  } // end for looping through the edges of body i

  evdwl += energy;

  return interact;
}

/* ----------------------------------------------------------------------
   Determine the interaction mode between i's edges against j's faces

   i = atom i (body i)
   j = atom j (body j)
   x      = atoms' coordinates
   f      = atoms' forces
   torque = atoms' torques
   tag    = atoms' tags
   contacts = list of contacts, the contacts between i's edges
              and j's faces are appended
   Return:

---------------------------------------------------------------------- */

int PairBodyRoundedPolyhedron::edge_against_face(int ibody, int jbody,
  int itype, int jtype, double** x, double** v, double** f, double** torque,
  double** angmom, double** fnc, Scratch &s, double &evdwl, double* facc)
{
  int ni,nei,nj,nfj,interact;
  double rradi,rradj,energy;

  nei = ednum[ibody];
  rradi = rounded_radius[ibody];
  nfj = facnum[jbody];
  rradj = rounded_radius[jbody];

  energy = 0;
  interact = EF_NONE;

  // loop through body i's edges

  for (ni = 0; ni < nei; ni++) {
    if (!edge_near(s, ibody, ni)) continue;

    // loop through body j's faces

    for (nj = 0; nj < nfj; nj++) {
      if (!face_near(s, jbody, nj)) continue;

      // compute the distance between the face nj to the edge ni
      #ifdef _POLYHEDRON_DEBUG
      printf("Compute interaction between face %d of body %d with "
             "edge %d of body %d:\n",
             nj, jbody, ni, ibody);
      #endif

      interact = interaction_face_to_edge(jbody, nj, x[jbody], rradj,
                                          ibody, ni, x[ibody], rradi,
                                          itype, jtype, cut_inner,
                                          v, f, torque, angmom, fnc, s,
                                          energy, facc);
    }

  } // end for looping through the edges of body i

  evdwl += energy;

  return interact;
}

/* -------------------------------------------------------------------------
  Compute the distance between an edge of body i and an edge from
  another body
  Input:
    ibody      = body i (i.e. atom i)
    face_index = face index of body i
    xmi        = atom i's coordinates (body i's center of mass)
    rounded_radius_i = rounded radius of the body i
    jbody      = body i (i.e. atom j)
    edge_index = coordinate of the tested edge from another body
    xmj        = atom j's coordinates (body j's center of mass)
    rounded_radius_j = rounded radius of the body j
    cut_inner  = cutoff for vertex-vertex and vertex-edge interaction
  Output:
    d          = Distance from a point x0 to an edge
    hi         = coordinates of the projection of x0 on the edge

  contact      = 0 no contact between the queried edge and the face
                 1 contact detected
  return
    INVALID if the face index is invalid
    NONE    if there is no interaction
------------------------------------------------------------------------- */

int PairBodyRoundedPolyhedron::interaction_edge_to_edge(int ibody,
  int edge_index_i,  double * /*xmi*/, double rounded_radius_i,
  int jbody, int edge_index_j, double * /*xmj*/, double rounded_radius_j,
  int itype, int jtype, double /*cut_inner*/, double** v, double** f,
  double** torque, double** angmom, double** fnc, Scratch &s,
  double &energy, double* facc)
{
  std::vector<Contact> &contacts = s.contacts;
  double** x = atom->x;

  // nearest points hi on the edge of body i and hj on the edge of body j

  double hi[3], hj[3], r;
  int crossed;
  if (!edge_edge_nearest(ibody, edge_index_i, jbody, edge_index_j, hi, hj, r, crossed))
    return EE_NONE;

  // the edges have crossed each other if the closest point on the edge
  // of body j is inside body i: use a negative distance so that the
  // overlap is deeper than the rounded radii and the force is repulsive

  if (crossed) r = -r;

  double w = edge_edge_weight(s, ibody, edge_index_i, jbody, edge_index_j);
  if (w <= 0.0) return EE_NONE;

  double contact_dist = rounded_radius_i + rounded_radius_j;
  int jflag = 1;
  pair_force_and_torque(jbody, ibody, hj, hi, r, contact_dist,
                        jtype, itype, x, v, f, torque, angmom,
                        fnc, jflag, energy, facc, w);

  if (r <= contact_dist) {
    // store the contact info
    Contact c;
    c.ibody = ibody;
    c.jbody = jbody;
    MathExtra::copy3(hi, c.xi);
    MathExtra::copy3(hj, c.xj);
    c.type = 1;
    c.separation = r - contact_dist;
    c.r = r;
    c.unique = 1;
    c.patch = 0;
    c.w = w;
    contacts.push_back(c);
  }
  return EE_INTERACT;
}

/* -------------------------------------------------------------------------
  Compute the interaction between a face of body i and an edge from
  another body
  Input:
    ibody      = body i (i.e. atom i)
    face_index = face index of body i
    xmi        = atom i's coordinates (body i's center of mass)
    rounded_radius_i = rounded radius of the body i
    jbody      = body i (i.e. atom j)
    edge_index = coordinate of the tested edge from another body
    xmj        = atom j's coordinates (body j's center of mass)
    rounded_radius_j = rounded radius of the body j
    cut_inner  = cutoff for vertex-vertex and vertex-edge interaction
  Output:
    d          = Distance from a point x0 to an edge
    hi         = coordinates of the projection of x0 on the edge

  contact      = 0 no contact between the queried edge and the face
                 1 contact detected
  return
    INVALID if the face index is invalid
    NONE    if there is no interaction
------------------------------------------------------------------------- */

int PairBodyRoundedPolyhedron::interaction_face_to_edge(int ibody,
  int face_index, double *xmi, double rounded_radius_i,
  int jbody, int edge_index, double *xmj, double rounded_radius_j,
  int itype, int jtype, double cut_inner, double** v, double** f,
  double** torque, double** angmom, double** fnc, Scratch &s,
  double &energy, double* facc)
{
  std::vector<Contact> &contacts = s.contacts;
  std::vector<int> &vertex_done = s.vertex_done;
  if (face_index >= facnum[ibody]) return EF_INVALID;

  int ifirst,iffirst,jfirst,npi1,npi2,npi3;
  int jefirst,npj1,npj2;
  double xi1[3],xi2[3],xi3[3],xpj1[3],xpj2[3],ui[3],vi[3],n[3];

  double** x = atom->x;

  ifirst = dfirst[ibody];
  iffirst = facfirst[ibody];
  npi1 = static_cast<int>(face[iffirst+face_index][0]);
  npi2 = static_cast<int>(face[iffirst+face_index][1]);
  npi3 = static_cast<int>(face[iffirst+face_index][2]);

  // compute the space-fixed coordinates for the vertices of the face

  xi1[0] = xmi[0] + discrete[ifirst+npi1][0];
  xi1[1] = xmi[1] + discrete[ifirst+npi1][1];
  xi1[2] = xmi[2] + discrete[ifirst+npi1][2];

  xi2[0] = xmi[0] + discrete[ifirst+npi2][0];
  xi2[1] = xmi[1] + discrete[ifirst+npi2][1];
  xi2[2] = xmi[2] + discrete[ifirst+npi2][2];

  xi3[0] = xmi[0] + discrete[ifirst+npi3][0];
  xi3[1] = xmi[1] + discrete[ifirst+npi3][1];
  xi3[2] = xmi[2] + discrete[ifirst+npi3][2];

  // find the normal unit vector of the face, ensure it point outward of the body

  MathExtra::sub3(xi2, xi1, ui);
  MathExtra::sub3(xi3, xi1, vi);
  MathExtra::cross3(ui, vi, n);
  MathExtra::norm3(n);

  double xc[3], dot, ans[3];
  xc[0] = (xi1[0] + xi2[0] + xi3[0])/3.0;
  xc[1] = (xi1[1] + xi2[1] + xi3[1])/3.0;
  xc[2] = (xi1[2] + xi2[2] + xi3[2])/3.0;
  MathExtra::sub3(xc, xmi, ans);
  dot = MathExtra::dot3(ans, n);
  if (dot < 0) MathExtra::negate3(n);

  // two ends of the edge from body j

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

  // signed distances of the two ends of the edge from the face plane
  // an end behind the face plane interacts with another face of body i,
  // or, if it is inside body i, see vertex_against_core()

  double hi1[3], hi2[3], contact_dist;
  int inside1 = 0;
  int inside2 = 0;

  contact_dist = rounded_radius_i + rounded_radius_j;
  double cut = contact_dist + cut_inner;
  double d1 = (xpj1[0]-xi1[0])*n[0] + (xpj1[1]-xi1[1])*n[1] + (xpj1[2]-xi1[2])*n[2];
  double d2 = (xpj2[0]-xi1[0])*n[0] + (xpj2[1]-xi1[1])*n[1] + (xpj2[2]-xi1[2])*n[2];
  int in1 = (d1 > 0.0) && (d1 <= cut) && (vertex_done[jfirst+npj1] == 0);
  int in2 = (d2 > 0.0) && (d2 <= cut) && (vertex_done[jfirst+npj2] == 0);
  if (!in1 && !in2) return EF_NONE;

  // projections of the two ends on the face plane, and whether they are
  // inside the face

  for (int k = 0; k < 3; k++) {
    hi1[k] = xpj1[k] - d1*n[k];
    hi2[k] = xpj2[k] - d2*n[k];
  }
  inside_polygon(ibody, face_index, xmi, hi1, hi2, inside1, inside2);

  int jflag = 1;

  // a vertex interacts only with the face of body i nearest to it, i.e. with
  // the largest signed distance, see nearest_face(): next to an acute edge
  // of body i, a vertex may project inside both faces of the edge

  double nc[3];
  double stol = EPSILON*EPSILON*contact_dist;
  if (in1 && inside1 && (d1 < nearest_face(ibody, xmi, xpj1, nc) - stol)) in1 = 0;
  if (in2 && inside2 && (d2 < nearest_face(ibody, xmi, xpj2, nc) - stol)) in2 = 0;

  // an end in the interaction zone whose projection is inside the face:
  // compute vertex-face interaction and accumulate force/torque to both bodies

  for (int k = 0; k < 2; k++) {
    if (!(k ? (in2 && inside2) : (in1 && inside1))) continue;
    int npj = k ? npj2 : npj1;
    double *xpj = k ? xpj2 : xpj1;
    double *hi = k ? hi2 : hi1;
    double d = k ? d2 : d1;

    double w = vertex_face_weight(s, jbody, npj, ibody, face_index);
    pair_force_and_torque(jbody, ibody, xpj, hi, d, contact_dist,
                          jtype, itype, x, v, f, torque, angmom,
                          fnc, jflag, energy, facc, w);

    if (d <= contact_dist) {
      // store the contact info
      Contact c;
      c.ibody = ibody;
      c.jbody = jbody;
      c.xi[0] = hi[0];
      c.xi[1] = hi[1];
      c.xi[2] = hi[2];
      c.xj[0] = xpj[0];
      c.xj[1] = xpj[1];
      c.xj[2] = xpj[2];
      c.type = 0;
      c.separation = d - contact_dist;
      c.r = d;
      c.unique = 1;
      c.patch = 0;
      c.w = w;
      contacts.push_back(c);
    }
    vertex_done[jfirst+npj] = 1;
  }

  return EF_INTERACT;
}

/* ----------------------------------------------------------------------
  Interaction of the vertices of body ibody that are inside the core of
  body jbody (the polyhedron without the rounded skin), which is the case
  when the overlap exceeds the contact distance: each such vertex is pushed
  out through the face of jbody nearest to it, with the (negative) signed
  distance to that face, so that the force keeps growing with the overlap
  the force on body ibody is accumulated to facc
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::vertex_against_core(int ibody, int jbody,
  int itype, int jtype, double** x, double** v, double** f, double** torque,
  double** angmom, double** fnc, Scratch &s, double &evdwl, double* facc)
{
  if (facnum[jbody] == 0) return;

  std::vector<Contact> &contacts = s.contacts;
  std::vector<int> &vertex_done = s.vertex_done;
  int ifirst = dfirst[ibody];
  double contact_dist = rounded_radius[ibody] + rounded_radius[jbody];
  double energy = 0.0;
  double xp[3], hp[3], nc[3];
  int nf;

  for (int ni = 0; ni < dnum[ibody]; ni++) {
    if (vertex_done[ifirst+ni] || !vertex_near(s, ibody, ni)) continue;

    MathExtra::add3(x[ibody], discrete[ifirst+ni], xp);
    double sd = nearest_face(jbody, x[jbody], xp, nc, &nf);
    if (sd > 0.0) continue;

    // avoid a zero distance when the vertex lies on the face plane

    sd = MIN(sd, -EPSILON*EPSILON*contact_dist);
    for (int k = 0; k < 3; k++) hp[k] = xp[k] - sd*nc[k];

    double w = vertex_face_weight(s, ibody, ni, jbody, nf);
    pair_force_and_torque(ibody, jbody, xp, hp, sd, contact_dist, itype, jtype, x, v, f,
                          torque, angmom, fnc, 1, energy, facc, w);

    Contact c;
    c.ibody = jbody;
    c.jbody = ibody;
    c.xi[0] = hp[0];
    c.xi[1] = hp[1];
    c.xi[2] = hp[2];
    c.xj[0] = xp[0];
    c.xj[1] = xp[1];
    c.xj[2] = xp[2];
    c.type = 0;
    c.separation = sd - contact_dist;
    c.r = sd;
    c.unique = 1;
    c.patch = 0;
    c.w = w;
    contacts.push_back(c);

    vertex_done[ifirst+ni] = 1;
  }

  evdwl += energy;
}

/* ----------------------------------------------------------------------
  Compute forces and torques between two bodies caused by the interaction
  between a pair of points on either bodies (similar to sphere-sphere)
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::pair_force_and_torque(int ibody, int jbody,
                 double* pi, double* pj, double r, double contact_dist,
                 int itype, int jtype, double** x,
                 double** /*v*/, double** f, double** torque, double** /*angmom*/,
                 double** /*fnc*/, int jflag, double& energy, double* facc, double w)
{
  double delx,dely,delz,R,fx,fy,fz,fpair;

  delx = pi[0] - pj[0];
  dely = pi[1] - pj[1];
  delz = pi[2] - pj[2];
  R = r - contact_dist;

  // the force and energy are scaled by the weight w of the interaction,
  // see vertex_face_weight()

  double e = 0.0;
  kernel_force(R, itype, jtype, e, fpair);

  fx = w*delx*fpair/r;
  fy = w*dely*fpair/r;
  fz = w*delz*fpair/r;

  #ifdef _POLYHEDRON_DEBUG
  printf("  - R = %f; r = %f; k_na = %f; shift = %f; fpair = %f;"
         " energy = %f; jflag = %d\n", R, r, k_na, shift, fpair,
         energy, jflag);
  #endif

  // contact: the elastic, cohesive and damping forces and the energy are
  // those of the unique contacts, see rescale_cohesive_forces()

  if (R > 0) {

    // accumulate force and torque to both bodies directly

    energy += w*e;

    f[ibody][0] += fx;
    f[ibody][1] += fy;
    f[ibody][2] += fz;
    sum_torque(x[ibody], pi, fx, fy, fz, torque[ibody]);

    facc[0] += fx; facc[1] += fy; facc[2] += fz;

    if (jflag) {
      f[jbody][0] -= fx;
      f[jbody][1] -= fy;
      f[jbody][2] -= fz;
      sum_torque(x[jbody], pj, -fx, -fy, -fz, torque[jbody]);
    }
  }
}

/* ----------------------------------------------------------------------
  Kernel force is model-dependent and can be derived for other styles
    here is the harmonic potential (linear piece-wise forces) in Wang et al.
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::kernel_force(double R, int itype, int jtype,
  double& energy, double& fpair)
{
  double fe, fc;
  energy += normal_force(R, itype, jtype, fe, fc);
  fpair = fe + fc;
}

/* ----------------------------------------------------------------------
   Normal force at the surface separation R, see Eq. 1, Wang et al.:
     fe = elastic (repulsive) force, k_n * delta_n, for R < 0
     fc = cohesive (attractive) force, -k_na * delta_na, for R <= cut_inner,
          where delta_na = cut_inner - R is the overlap of the cohesive regions,
          which keeps growing when the surfaces deform
   return the energy, which is zero at R = cut_inner
------------------------------------------------------------------------- */

double PairBodyRoundedPolyhedron::normal_force(double R, int itype, int jtype,
                                               double &fe, double &fc)
{
  fe = fc = 0.0;
  if (R > cut_inner) return 0.0;

  double kn = k_n[itype][jtype];
  double kna = k_na[itype][jtype];
  double dna = cut_inner - R;
  fc = -kna * dna;
  double energy = -0.5 * kna * dna * dna;
  if (R < 0.0) {
    fe = -kn * R;
    energy += 0.5 * kn * R * R;
  }
  return energy;
}

/* ----------------------------------------------------------------------
  Rescale the forces and torques for all the contacts
  the total force on body iref is accumulated to facc
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::rescale_cohesive_forces(double** x,
     double** v, double** angmom, double** f, double** torque, double** fnc,
     std::vector<Contact> &contacts, int itype, int jtype, int iref, double &evdwl,
     double* facc)
{
  int m,ibody,jbody;
  double delx,dely,delz,fx,fy,fz,R,fpair,r,contact_area,w;

  const int num_contacts = contacts.size();

  // contact area A_a = pi <r_k^2> from the distances r_k of the unique
  // contacts to their average position, see Eq. 6, Wang et al.
  // the averages and the number of contacts are weighted by the weights of
  // the contacts, see vertex_face_weight(), so that they change continuously

  if (num_contacts > 1) find_unique_contacts(contacts);

  double xc[3] = {0.0, 0.0, 0.0};
  double wsum = 0.0;
  for (m = 0; m < num_contacts; m++) {
    if (contacts[m].unique == 0) continue;
    w = contacts[m].w;
    xc[0] += w * contacts[m].xi[0];
    xc[1] += w * contacts[m].xi[1];
    xc[2] += w * contacts[m].xi[2];
    wsum += w;
  }
  if (wsum <= 0.0) return;
  MathExtra::scale3(1.0/wsum, xc);

  contact_area = 0.0;
  for (m = 0; m < num_contacts; m++) {
    if (contacts[m].unique == 0) continue;
    contact_area += contacts[m].w * MathExtra::distsq3(contacts[m].xi, xc);
  }
  contact_area *= MY_PI / wsum;

  double j_a = contact_area / (wsum * A_ua);
  if (j_a < 1.0) j_a = 1.0;
  for (m = 0; m < num_contacts; m++) {
    if (contacts[m].unique == 0) continue;

    ibody = contacts[m].ibody;
    jbody = contacts[m].jbody;
    w = contacts[m].w;

    delx = contacts[m].xi[0] - contacts[m].xj[0];
    dely = contacts[m].xi[1] - contacts[m].xj[1];
    delz = contacts[m].xi[2] - contacts[m].xj[2];
    r = contacts[m].r;
    R = contacts[m].separation;

    // only the cohesive force is scaled by j_a, see Eq. 6, Wang et al.

    double fe, fc;
    evdwl += w * normal_force(R, itype, jtype, fe, fc);
    fpair = w * (fe + j_a * fc);
    fx = delx*fpair/r;
    fy = dely*fpair/r;
    fz = delz*fpair/r;

    f[ibody][0] += fx;
    f[ibody][1] += fy;
    f[ibody][2] += fz;
    sum_torque(x[ibody], contacts[m].xi, fx, fy, fz, torque[ibody]);

    f[jbody][0] -= fx;
    f[jbody][1] -= fy;
    f[jbody][2] -= fz;
    sum_torque(x[jbody], contacts[m].xj, -fx, -fy, -fz, torque[jbody]);

    // facc is the force on body iref

    if (ibody == iref) {
      facc[0] += fx; facc[1] += fy; facc[2] += fz;
    } else {
      facc[0] -= fx; facc[1] -= fy; facc[2] -= fz;
    }

    // the part of the force added by the j_a scaling does not derive
    // from the energy

    double s = w * (j_a - 1.0) * fc / r;
    double fja[3] = {s*delx, s*dely, s*delz};
    fnc[ibody][0] += fja[0];
    fnc[ibody][1] += fja[1];
    fnc[ibody][2] += fja[2];
    fnc[jbody][0] -= fja[0];
    fnc[jbody][1] -= fja[1];
    fnc[jbody][2] -= fja[2];
    sum_torque(x[ibody], contacts[m].xi, fja[0], fja[1], fja[2], &fnc[ibody][3]);
    sum_torque(x[jbody], contacts[m].xj, -fja[0], -fja[1], -fja[2], &fnc[jbody][3]);

    // damping at the contact point between the surfaces, the friction force
    // is computed once per pair of bodies, see friction_forces()

    if (r != 0.0) {
      double n[3] = {delx/r, dely/r, delz/r};
      double pc[3];
      contact_point(contacts[m].xi, contacts[m].xj, n, rounded_radius[ibody],
                    rounded_radius[jbody], pc);
      damping_friction(ibody, jbody, pc, n, 0.0, 1, 0, x, v, angmom, f, torque, fnc, iref, facc,
                       nullptr, w);
    }
  }
}

/* ----------------------------------------------------------------------
  Contact point between the rounded surfaces of two bodies, halfway between
  the surface points on the line through the points pi on ibody and pj on
  jbody, where n is the unit normal pointing from jbody to ibody
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::contact_point(const double *pi, const double *pj,
                                              const double *n, double rradi, double rradj,
                                              double *pc)
{
  for (int k = 0; k < 3; k++)
    pc[k] = 0.5 * ((pi[k] - rradi * n[k]) + (pj[k] + rradj * n[k]));
}

/* ----------------------------------------------------------------------
  Damping and friction forces at the contact point pc between two bodies,
  from the relative velocity of the two bodies at that point:
    damping:  -c_n v_n - c_t v_t
    friction: magnitude mu * fne with the elastic normal force fne, opposite
              to v_t, capped by c_t * |v_t| so that it vanishes smoothly
              as sliding stops, see Eq. 4, Wang et al.
              with the contact history, the friction force is that of a
              tangential spring instead, see tangential_spring()
  n = unit normal pointing from jbody to ibody
  the forces act at pc on both bodies, so that they exert torques
  the total force on body iref is accumulated to facc
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::damping_friction(int ibody, int jbody, const double *pc,
  const double *n, double fne, int damping, int friction, double** x, double** v,
  double** angmom, double** f, double** torque, double** fnc, int iref, double* facc,
  Scratch *hs, double w)
{
  double vi[3], vj[3], vr[3], vn[3], vt[3], fdiss[3];
  AtomVecBody::Bonus *bonus;

  bonus = &avec->bonus[atom->body[ibody]];
  total_velocity(const_cast<double *>(pc), x[ibody], v[ibody], angmom[ibody], bonus->inertia,
                 bonus->quat, vi);
  bonus = &avec->bonus[atom->body[jbody]];
  total_velocity(const_cast<double *>(pc), x[jbody], v[jbody], angmom[jbody], bonus->inertia,
                 bonus->quat, vj);

  MathExtra::sub3(vi, vj, vr);
  double vnnr = MathExtra::dot3(vr, n);
  for (int k = 0; k < 3; k++) {
    vn[k] = vnnr * n[k];
    vt[k] = vr[k] - vn[k];
    fdiss[k] = 0.0;
  }

  // the damping is scaled by the weight w of the contact, see vertex_face_weight(),
  // the elastic normal force fne of the friction includes it

  if (damping)
    for (int k = 0; k < 3; k++) fdiss[k] = -w * (c_n * vn[k] + c_t * vt[k]);

  double vtmag = MathExtra::len3(vt);
  if (friction && hs && hs->shear) {

    // with the contact history, a single friction force of the pair of
    // bodies is computed from all its contacts in history_friction()

    if (fne > 0.0) {
      double sign = (ibody == hs->shear_i) ? 1.0 : -1.0;
      hs->fnsum += fne;
      for (int k = 0; k < 3; k++) {
        hs->pcsum[k] += fne * pc[k];
        hs->vtsum[k] += sign * fne * vt[k];
        hs->nsum[k] += sign * fne * n[k];
      }
    }
    if (!damping) return;
  } else if (friction && (fne > 0.0) && (vtmag > 0.0)) {
    double scale = MIN(mu * fne, c_t * vtmag) / vtmag;
    for (int k = 0; k < 3; k++) fdiss[k] -= scale * vt[k];
  }

  f[ibody][0] += fdiss[0];
  f[ibody][1] += fdiss[1];
  f[ibody][2] += fdiss[2];
  sum_torque(x[ibody], const_cast<double *>(pc), fdiss[0], fdiss[1], fdiss[2], torque[ibody]);

  f[jbody][0] -= fdiss[0];
  f[jbody][1] -= fdiss[1];
  f[jbody][2] -= fdiss[2];
  sum_torque(x[jbody], const_cast<double *>(pc), -fdiss[0], -fdiss[1], -fdiss[2], torque[jbody]);

  // damping and friction do not derive from the energy

  for (int k = 0; k < 3; k++) {
    fnc[ibody][6+k] += fdiss[k];
    fnc[jbody][6+k] -= fdiss[k];
  }
  sum_torque(x[ibody], const_cast<double *>(pc), fdiss[0], fdiss[1], fdiss[2], &fnc[ibody][9]);
  sum_torque(x[jbody], const_cast<double *>(pc), -fdiss[0], -fdiss[1], -fdiss[2],
             &fnc[jbody][9]);

  if (ibody == iref) {
    facc[0] += fdiss[0]; facc[1] += fdiss[1]; facc[2] += fdiss[2];
  } else {
    facc[0] -= fdiss[0]; facc[1] -= fdiss[1]; facc[2] -= fdiss[2];
  }
}

/* ----------------------------------------------------------------------
  Friction force at a contact point from a tangential spring, with the
  tangential displacement xi of the pair of bodies stored in s.shear,
  see e.g. Luding, Granular Matter 10, 235 (2008):
  xi is rotated into the tangent plane of the contact normal n, keeping its
  magnitude, since the normal changes, and xi is incremented by the
  tangential relative velocity vt times dt.  The spring force -k_t xi is
  limited to mu times the elastic normal force fne, and then xi is reduced
  accordingly (sliding), see history_friction().
  the force on body ibody is added to fs, xi refers to body s.shear_i
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::tangential_spring(int ibody, int jbody, const double *n,
                                                  const double *vt, double fne, Scratch &s,
                                                  double *fs)
{
  double sign = (ibody == s.shear_i) ? 1.0 : -1.0;
  double xi[3], ft[3];
  for (int k = 0; k < 3; k++) xi[k] = sign * s.shear[k];

  double xin = MathExtra::dot3(xi, n);
  double mag = MathExtra::len3(xi);
  for (int k = 0; k < 3; k++) xi[k] -= xin * n[k];
  double magt = MathExtra::len3(xi);
  if (magt > 0.0) MathExtra::scale3(mag/magt, xi);

  if (!update->setupflag)
    for (int k = 0; k < 3; k++) xi[k] += vt[k] * dt;

  double kt = k_t[atom->type[ibody]][atom->type[jbody]];
  for (int k = 0; k < 3; k++) ft[k] = -kt * xi[k];

  double ftmag = MathExtra::len3(ft);
  double ftmax = mu * MAX(fne, 0.0);
  if (ftmag > ftmax) {
    double scale = ftmax / ftmag;
    MathExtra::scale3(scale, ft);
    MathExtra::scale3(scale, xi);
  }

  for (int k = 0; k < 3; k++) {
    fs[k] += ft[k];
    s.shear[k] = sign * xi[k];
  }
  s.touched = 1;
}

/* ----------------------------------------------------------------------
  Friction force of the pair of bodies i and j with the contact history:
  a single tangential spring, see tangential_spring(), for all contacts of
  the pair, limited by mu times the sum of their elastic normal forces.
  It acts at the average of the contact points, with the average tangential
  relative velocity and normal, all weighted by the elastic normal forces of
  the contacts, as collected by damping_friction().  A single friction force
  at the contact with the largest overlap instead would jump between the
  contacts when that contact changes, e.g. for a body resting with a face on
  another body under a lateral load, so that the bodies do not come to rest.
  note: this deviates from Wang et al., where the single friction force of
  a pair acts at the contact point with the largest solid overlap, with the
  relative velocity and normal of that contact.  The limit mu F_ne uses the
  elastic normal force of the whole interaction.  Without the contact
  history, the friction force still acts at the contact with the largest
  overlap.
  the total force on body i is accumulated to facc
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::history_friction(int i, int j, double **x, double **f,
                                                 double **torque, double **fnc, Scratch &s,
                                                 double *facc)
{
  if (s.fnsum <= 0.0) return;

  double pc[3], vt[3], n[3], ft[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    pc[k] = s.pcsum[k] / s.fnsum;
    vt[k] = s.vtsum[k] / s.fnsum;
    n[k] = s.nsum[k];
  }

  // the tangential velocities of the contacts refer to their own normals

  double nmag = MathExtra::len3(n);
  if (nmag > 0.0) {
    MathExtra::scale3(1.0/nmag, n);
    double vtn = MathExtra::dot3(vt, n);
    for (int k = 0; k < 3; k++) vt[k] -= vtn * n[k];
  }

  tangential_spring(i, j, n, vt, s.fnsum, s, ft);

  f[i][0] += ft[0];
  f[i][1] += ft[1];
  f[i][2] += ft[2];
  sum_torque(x[i], pc, ft[0], ft[1], ft[2], torque[i]);

  f[j][0] -= ft[0];
  f[j][1] -= ft[1];
  f[j][2] -= ft[2];
  sum_torque(x[j], pc, -ft[0], -ft[1], -ft[2], torque[j]);

  // friction does not derive from the energy

  for (int k = 0; k < 3; k++) {
    fnc[i][6+k] += ft[k];
    fnc[j][6+k] -= ft[k];
  }
  sum_torque(x[i], pc, ft[0], ft[1], ft[2], &fnc[i][9]);
  sum_torque(x[j], pc, -ft[0], -ft[1], -ft[2], &fnc[j][9]);

  facc[0] += ft[0]; facc[1] += ft[1]; facc[2] += ft[2];
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::reset_dt()
{
  // with run_style respa, the pair style is computed with the time step
  // of its level

  dt = update->dt;
  if (utils::strmatch(update->integrate_style, "^respa")) {
    auto *respa = dynamic_cast<Respa *>(update->integrate);
    if (respa && (respa->level_pair >= 0)) dt = respa->step[respa->level_pair];
  }
}

/* ----------------------------------------------------------------------
  Friction force of a pair of bodies from its contacts in s.contacts:
  at the contact with the largest overlap, see Wang et al., or with the
  contact history from all unique contacts, see history_friction()
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::friction_forces(int itype, int jtype, double** x, double** v,
                                                double** angmom, double** f, double** torque,
                                                double** fnc, int iref, double* facc, Scratch &s)
{
  std::vector<Contact> &contacts = s.contacts;
  if (contacts.empty()) return;

  if (s.shear) {
    for (auto &contact : contacts)
      if (contact.unique)
        friction_force(contact, itype, jtype, x, v, angmom, f, torque, fnc, iref, facc, s);
    return;
  }

  // the contact with the largest elastic normal force, i.e. the largest
  // weighted overlap, so that the choice does not jump between coincident
  // contacts whose weights change with the orientation, see
  // vertex_face_weight()

  int mmax = -1;
  double fmax = 0.0;
  for (int m = 0; m < (int) contacts.size(); m++) {
    if (!contacts[m].unique) continue;
    double fm = -contacts[m].separation * contacts[m].w;
    if ((mmax < 0) || (fm > fmax)) {
      mmax = m;
      fmax = fm;
    }
  }
  if (mmax >= 0)
    friction_force(contacts[mmax], itype, jtype, x, v, angmom, f, torque, fnc, iref, facc, s);
}

/* ----------------------------------------------------------------------
  Friction force at a contact during gross sliding, see Eq. 4, Wang et al.:
  magnitude mu * F_ne with the elastic normal force F_ne, opposite to the
  tangential relative velocity, capped by c_t * |v_t| so that it vanishes
  smoothly as sliding stops instead of reversing the sliding direction
  within a time step, or with the contact history from a tangential spring,
  see tangential_spring()
  the total force on body iref is accumulated to facc
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::friction_force(Contact &contact, int itype, int jtype,
  double** x, double** v, double** angmom, double** f, double** torque, double** fnc,
  int iref, double* facc, Scratch &s)
{
  if ((contact.separation >= 0.0) || (contact.r == 0.0)) return;

  double n[3], pc[3];
  MathExtra::sub3(contact.xi, contact.xj, n);
  MathExtra::scale3(1.0/contact.r, n);
  contact_point(contact.xi, contact.xj, n, rounded_radius[contact.ibody],
                rounded_radius[contact.jbody], pc);

  double fne = -k_n[itype][jtype] * contact.separation * contact.w;
  damping_friction(contact.ibody, contact.jbody, pc, n, fne, 0, 1, x, v, angmom, f, torque, fnc,
                   iref, facc, &s);
}

/* ----------------------------------------------------------------------
  Accumulate torque to body from the force f=(fx,fy,fz) acting at point x
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::sum_torque(double* xm, double *x, double fx,
                                      double fy, double fz, double* torque)
{
  double rx = x[0] - xm[0];
  double ry = x[1] - xm[1];
  double rz = x[2] - xm[2];
  double tx = ry * fz - rz * fy;
  double ty = rz * fx - rx * fz;
  double tz = rx * fy - ry * fx;
  torque[0] += tx;
  torque[1] += ty;
  torque[2] += tz;
}

/* ----------------------------------------------------------------------
  Find the face of body ibody whose plane has the largest signed distance
  to the point q, the distance is positive when q is in front of the face
  (outside of the body) and negative when q is behind the face
  Input:
    ibody = body i (i.e. atom i)
    xmi   = atom i's coordinates (body i's center of mass)
    q     = tested point
  Output:
    n     = outward unit normal of the face plane
    nf    = index of the face, if not a null pointer
  return the signed distance from q to the face plane,
    which is negative only if q is inside the (convex) body
------------------------------------------------------------------------- */

double PairBodyRoundedPolyhedron::nearest_face(int ibody, double * /*xmi*/,
                                               const double *q, double *n, int *nf) const
{
  int iffirst = facfirst[ibody];
  double smax = 0.0;
  double ans[3];

  // planes of the faces with their outward unit normals, see body2space()

  for (int nf_index = 0; nf_index < facnum[ibody]; nf_index++) {
    double *plane = facplane[iffirst+nf_index];
    MathExtra::sub3(q, plane, ans);
    double s = MathExtra::dot3(ans, &plane[3]);
    if ((nf_index == 0) || (s > smax)) {
      smax = s;
      MathExtra::copy3(&plane[3], n);
      if (nf) *nf = nf_index;
    }
  }

  return smax;
}

/* ----------------------------------------------------------------------
  Test if two points a and b are in opposite sides of a plane defined by
  a normal vector n and a point x0
------------------------------------------------------------------------- */

int PairBodyRoundedPolyhedron::opposite_sides(double* n, double* x0,
                                           double* a, double* b)
{
  double m_a = n[0]*(a[0] - x0[0])+n[1]*(a[1] - x0[1])+n[2]*(a[2] - x0[2]);
  double m_b = n[0]*(b[0] - x0[0])+n[1]*(b[1] - x0[1])+n[2]*(b[2] - x0[2]);
  // equal to zero when either a or b is on the plane
  if (m_a * m_b <= 0)
    return 1;
  else
    return 0;
}

/* ----------------------------------------------------------------------
  Find the projection of q on the plane defined by point p and the normal
  unit vector n: q_proj = q - dot(q - p, n) * n
  and the distance d from q to the plane
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::project_pt_plane(const double* q,
                                        const double* p, const double* n,
                                        double* q_proj, double &d) const
{
  double dot, ans[3], n_p[3];
  n_p[0] = n[0]; n_p[1] = n[1]; n_p[2] = n[2];
  MathExtra::sub3(q, p, ans);
  dot = MathExtra::dot3(ans, n_p);
  MathExtra::scale3(dot, n_p);
  MathExtra::sub3(q, n_p, q_proj);
  MathExtra::sub3(q, q_proj, ans);
  d = MathExtra::len3(ans);
}

/* ----------------------------------------------------------------------
  Check if points q1 and q2 are inside a convex polygon, i.e. a face of
  a polyhedron
    ibody       = atom i's index
    face_index  = face index of the body
    xmi         = atom i's coordinates
    q1          = tested point on the face (e.g. the projection of a point)
    q2          = another point (can be a null pointer) for face-edge intersection
  Output:
    inside1     = 1 if q1 is inside the polygon, 0 otherwise
    inside2     = 1 if q2 is inside the polygon, 0 otherwise
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::inside_polygon(int ibody, int face_index,
                            double* xmi, const double* q1, const double* q2,
                            int& inside1, int& inside2) const

{
  int i,n,ifirst,iffirst,npi1,npi2;
  double xi1[3],xi2[3],u[3],v[3],costheta,anglesum1,anglesum2,magu,magv,rradi;

  ifirst = dfirst[ibody];
  iffirst = facfirst[ibody];
  rradi = rounded_radius[ibody];
  double rradsq = rradi*rradi;
  anglesum1 = anglesum2 = 0;
  int atvertex1 = 0, atvertex2 = 0;
  for (i = 0; i < MAX_FACE_SIZE; i++) {
    npi1 = static_cast<int>(face[iffirst+face_index][i]);
    if (npi1 < 0) break;
    n = i + 1;
    if (n <= MAX_FACE_SIZE - 1) {
      npi2 = static_cast<int>(face[iffirst+face_index][n]);
      if (npi2 < 0) npi2 = static_cast<int>(face[iffirst+face_index][0]);
    } else {
      npi2 = static_cast<int>(face[iffirst+face_index][0]);
    }

    xi1[0] = xmi[0] + discrete[ifirst+npi1][0];
    xi1[1] = xmi[1] + discrete[ifirst+npi1][1];
    xi1[2] = xmi[2] + discrete[ifirst+npi1][2];

    xi2[0] = xmi[0] + discrete[ifirst+npi2][0];
    xi2[1] = xmi[1] + discrete[ifirst+npi2][1];
    xi2[2] = xmi[2] + discrete[ifirst+npi2][2];

    MathExtra::sub3(xi1,q1,u);
    MathExtra::sub3(xi2,q1,v);
    magu = MathExtra::len3(u);
    magv = MathExtra::len3(v);

    // the point is at either vertices

    if (magu * magv < EPSILON*rradsq) atvertex1 = 1;
    else {
      costheta = MathExtra::dot3(u,v)/(magu*magv);
      anglesum1 += acos(MAX(-1.0, MIN(1.0, costheta)));
    }

    if (q2 != nullptr) {
      MathExtra::sub3(xi1,q2,u);
      MathExtra::sub3(xi2,q2,v);
      magu = MathExtra::len3(u);
      magv = MathExtra::len3(v);
      if (magu * magv < EPSILON*rradsq) atvertex2 = 1;
      else {
        costheta = MathExtra::dot3(u,v)/(magu*magv);
        anglesum2 += acos(MAX(-1.0, MIN(1.0, costheta)));
      }
    }
  }

  // a point at a vertex is inside

  if (atvertex1 || (fabs(anglesum1 - MY_2PI) < EPSILON)) inside1 = 1;
  else inside1 = 0;

  if (q2 != nullptr) {
    if (atvertex2 || (fabs(anglesum2 - MY_2PI) < EPSILON)) inside2 = 1;
    else inside2 = 0;
  }
}

/* ----------------------------------------------------------------------
  Find the projection of q on the plane defined by 3 points x1, x2 and x3
  returns the distance d from q to the plane and whether the projected
  point is inside the triangle defined by (x1, x2, x3)
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::project_pt_plane(const double* q,
      const double* x1, const double* x2, const double* x3, double* q_proj,
      double &d, int& inside) const
{
  double u[3],v[3],n[3];

  // plane normal vector

  MathExtra::sub3(x2, x1, u);
  MathExtra::sub3(x3, x1, v);
  MathExtra::cross3(u, v, n);
  MathExtra::norm3(n);

  // solve for the intersection between the line and the plane

  double m[3][3], invm[3][3], p[3], ans[3];
  m[0][0] = -n[0];
  m[0][1] = u[0];
  m[0][2] = v[0];

  m[1][0] = -n[1];
  m[1][1] = u[1];
  m[1][2] = v[1];

  m[2][0] = -n[2];
  m[2][1] = u[2];
  m[2][2] = v[2];

  MathExtra::sub3(q, x1, p);
  MathExtra::invert3(m, invm);
  MathExtra::matvec(invm, p, ans);

  double t = ans[0];
  q_proj[0] = q[0] + n[0] * t;
  q_proj[1] = q[1] + n[1] * t;
  q_proj[2] = q[2] + n[2] * t;

  // check if the projection point is inside the triangle
  // exclude the edges and vertices
  // edge-sphere and sphere-sphere interactions are handled separately

  inside = 0;
  if (ans[1] > 0 && ans[2] > 0 && ans[1] + ans[2] < 1) {
    inside = 1;
  }

  // distance from q to q_proj

  MathExtra::sub3(q, q_proj, ans);
  d = MathExtra::len3(ans);
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::project_pt_line(const double* q,
     const double* xi1, const double* xi2, double* h, double& d, double& t) const
{
  double u[3],v[3],r[3],s;

  MathExtra::sub3(xi2, xi1, u);
  MathExtra::norm3(u);
  MathExtra::sub3(q, xi1, v);

  s = MathExtra::dot3(u, v);
  h[0] = xi1[0] + s * u[0];
  h[1] = xi1[1] + s * u[1];
  h[2] = xi1[2] + s * u[2];

  MathExtra::sub3(q, h, r);
  d = MathExtra::len3(r);

  // fraction of h along the edge, u is the unit director of the edge

  t = s / sqrt(MathExtra::distsq3(xi1, xi2));
}

/* ----------------------------------------------------------------------
  compute the shortest distance between two edges (line segments)
  x1, x2: two endpoints of the first edge
  x3, x4: two endpoints of the second edge
  h1: the end point of the shortest segment perpendicular to both edges
      on the line (x1;x2)
  h2: the end point of the shortest segment perpendicular to both edges
      on the line (x3;x4)
  t1: fraction of h1 in the segment (x1,x2)
  t2: fraction of h2 in the segment (x3,x4)
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::distance_bt_edges(const double* x1,
                  const double* x2, const double* x3, const double* x4,
                  double* h1, double* h2, double& t1, double& t2, double& r) const
{
  double u[3],v[3],n[3],dot;

  // set the default returned values

  t1 = -2;
  t2 = 2;
  r = 0;

  // find the edge unit directors and their dot product

  MathExtra::sub3(x2, x1, u);
  MathExtra::norm3(u);
  MathExtra::sub3(x4, x3, v);
  MathExtra::norm3(v);
  dot = MathExtra::dot3(u,v);
  dot = fabs(dot);

  // parallel edges have no pair of nearest points inside both edges: the
  // ends of their overlap are vertices of either edge, and the contacts there
  // are those of the vertices, as for nearly parallel edges whose nearest
  // points are outside of the edges.  A single contact at an end or the middle
  // of the overlap instead would jump between these points as the edges move
  // along each other, and would get a tilted normal when they are offset
  // sideways, while the faces next to the edges are closer, see edge_cone().
  // only the distance r between the lines is returned

  if (1.0 - dot < PARALLEL_TOL) {
    double x13[3];
    MathExtra::sub3(x1, x3, x13);
    double s1 = MathExtra::dot3(x13, v);
    for (int k = 0; k < 3; k++) x13[k] -= s1*v[k];
    r = MathExtra::len3(x13);
    return;
  }

  // find the vector n perpendicular to both edges

  MathExtra::cross3(u, v, n);
  MathExtra::norm3(n);

  // find the intersection of the line (x3,x4) and the plane (x1,x2,n)
  // s = director of the line (x3,x4)
  // n_p = plane normal vector of the plane (x1,x2,n)

  double s[3], n_p[3];
  MathExtra::sub3(x4, x3, s);
  MathExtra::sub3(x2, x1, u);
  MathExtra::cross3(u, n, n_p);
  MathExtra::norm3(n_p);

  // solve for the intersection between the line and the plane

  double m[3][3], invm[3][3], p[3], ans[3];
  m[0][0] = -s[0];
  m[0][1] = u[0];
  m[0][2] = n[0];

  m[1][0] = -s[1];
  m[1][1] = u[1];
  m[1][2] = n[1];

  m[2][0] = -s[2];
  m[2][1] = u[2];
  m[2][2] = n[2];

  MathExtra::sub3(x3, x1, p);
  MathExtra::invert3(m, invm);
  MathExtra::matvec(invm, p, ans);

  t2 = ans[0];
  h2[0] = x3[0] + s[0] * t2;
  h2[1] = x3[1] + s[1] * t2;
  h2[2] = x3[2] + s[2] * t2;

  project_pt_plane(h2, x1, n, h1, r);

  // fraction of h1 along the edge from its projection, since a division by
  // a single coordinate difference is inaccurate when the edge is nearly
  // perpendicular to that axis

  double h1x1[3];
  MathExtra::sub3(h1, x1, h1x1);
  t1 = MathExtra::dot3(h1x1, u) / MathExtra::dot3(u, u);
}

/* ----------------------------------------------------------------------
  Calculate the total velocity of a point (vertex, a point on an edge):
    vi = vcm + omega ^ (p - xcm)
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::total_velocity(double* p, double *xcm,
  double* vcm, double *angmom, double *inertia, double *quat, double* vi)
{
  double r[3],omega[3],ex_space[3],ey_space[3],ez_space[3];
  r[0] = p[0] - xcm[0];
  r[1] = p[1] - xcm[1];
  r[2] = p[2] - xcm[2];
  MathExtra::q_to_exyz(quat,ex_space,ey_space,ez_space);
  MathExtra::angmom_to_omega(angmom,ex_space,ey_space,ez_space,
                             inertia,omega);
  vi[0] = omega[1]*r[2] - omega[2]*r[1] + vcm[0];
  vi[1] = omega[2]*r[0] - omega[0]*r[2] + vcm[1];
  vi[2] = omega[0]*r[1] - omega[1]*r[0] + vcm[2];
}

/* ----------------------------------------------------------------------
  Determine the length of the contact segment, i.e. the separation between
  2 contacts, should be extended for 3D models.
------------------------------------------------------------------------- */

double PairBodyRoundedPolyhedron::contact_separation(const Contact& c1,
                                                     const Contact& c2)
{
  double x1 = 0.5*(c1.xi[0] + c1.xj[0]);
  double y1 = 0.5*(c1.xi[1] + c1.xj[1]);
  double z1 = 0.5*(c1.xi[2] + c1.xj[2]);
  double x2 = 0.5*(c2.xi[0] + c2.xj[0]);
  double y2 = 0.5*(c2.xi[1] + c2.xj[1]);
  double z2 = 0.5*(c2.xi[2] + c2.xj[2]);
  double rsq = (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) + (z2 - z1)*(z2 - z1);
  return rsq;
}

/* ----------------------------------------------------------------------
   find the number of unique contacts
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::find_unique_contacts(std::vector<Contact> &contacts)
{
  int n = contacts.size();
  for (int i = 0; i < n - 1; i++) {

    for (int j = i + 1; j < n; j++) {
      if (contacts[i].unique == 0) continue;

      // the corners of a region of overlap are distinct, and nearby corners
      // have weights that add up, see face_face_patches()

      if (contacts[i].patch || contacts[j].patch) continue;
      double d = contact_separation(contacts[i], contacts[j]);
      int ibody = contacts[i].ibody;
      int jbody = contacts[i].jbody;
      double rradi = rounded_radius[ibody];
      double rradj = rounded_radius[jbody];
      double rmin = MIN(rradi, rradj);
      if (d < EPSILON*EPSILON*rmin*rmin) contacts[j].unique = 0;
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairBodyRoundedPolyhedron::memory_usage()
{
  double bytes = Pair::memory_usage();
  bytes += (double) nmax * 6 * sizeof(int);    // dnum+dfirst+ednum+edfirst+facnum+facfirst [nmax]
  return bytes;
}
