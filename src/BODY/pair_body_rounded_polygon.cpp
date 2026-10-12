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
   Ref: Fraige, Langston, Matchett and Dodds, Particuology 2008, 6:455-466
   Note: The current implementation has not taken into account
         the contact history for friction forces.
------------------------------------------------------------------------- */

#include "pair_body_rounded_polygon.h"

#include "atom.h"
#include "atom_vec_body.h"
#include "body_rounded_polygon.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "fix_neigh_history.h"
#include "fix_store_atom.h"
#include "force.h"
#include "group.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

static constexpr int DELTA = 10000;
static constexpr double EPSILON = 1.0e-3; // dimensionless threshold (dot products, end point checks, contact checks)
static constexpr int EFF_CONTACTS = 2;    // effective contacts for 2D models
static constexpr double BIG = 1.0e20;
static constexpr int NFNC = 12;           // per-body force and torque of the j_a scaling, and of damping
static constexpr char id_fix_store_prefix[] = "BODY_ROUNDED_POLYGON_WORK_";

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolygon::PairBodyRoundedPolygon(LAMMPS *lmp) :
    Pair(lmp), k_n(nullptr), k_na(nullptr), avec(nullptr), bptr(nullptr)
{
  dmax = nmax = 0;
  discrete = nullptr;
  dnum = dfirst = nullptr;

  edmax = ednummax = 0;
  edge = nullptr;
  ednum = edfirst = nullptr;

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
  delta_ua = 1.0;

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

PairBodyRoundedPolygon::~PairBodyRoundedPolygon()
{
  memory->destroy(discrete);
  memory->destroy(dnum);
  memory->destroy(dfirst);

  memory->destroy(edge);
  memory->destroy(ednum);
  memory->destroy(edfirst);

  memory->destroy(enclosing_radius);
  memory->destroy(rounded_radius);

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
    memory->destroy(maxrad);
  }
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolygon::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl;
  double rsq,r,radi,radj;
  double facc[3];
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
    memory->destroy(enclosing_radius);
    memory->destroy(rounded_radius);
    nmax = atom->nmax;
    memory->create(dnum,nmax,"pair:dnum");
    memory->create(dfirst,nmax,"pair:dfirst");
    memory->create(ednum,nmax,"pair:ednum");
    memory->create(edfirst,nmax,"pair:edfirst");
    memory->create(enclosing_radius,nmax,"pair:enclosing_radius");
    memory->create(rounded_radius,nmax,"pair:rounded_radius");
  }

  ndiscrete = nedge = 0;
  for (i = 0; i < nall; i++)
    dnum[i] = ednum[i] = 0;

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
    radi = radius[i];
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
      radj = radius[j];

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

      // no interaction

      // note: body/rounded/polyhedron additionally skips the pairs, and the
      // vertices, edges and faces, that do not reach far enough towards the
      // other body along the line between the centers to interact with it,
      // see PairBodyRoundedPolyhedron::pair_interaction().  The same test in
      // this pair style gave the same results, but no speedup for
      // examples/body/in.squares and in.wall2d: each vertex is already tested
      // cheaply against the enclosing circle, and in these dense systems only
      // 7-19% of the pairs passing the test below could be skipped as a whole.
      // It may still help for dilute systems.

      r = sqrt(rsq);
      if ((body[i] >= 0) && (body[j] >= 0) && (r <= radi + radj + cut_inner)) {
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

void PairBodyRoundedPolygon::pair_interaction(int i, int j, double delx, double dely,
                                              double delz, double rsq, double **x,
                                              double **v, double **angmom, double **f,
                                              double **torque, double **fnc, Scratch &s,
                                              double &evdwl, double *facc)
{
  int itype = atom->type[i];
  int jtype = atom->type[j];
  int npi = dnum[i];
  int npj = dnum[j];
  double k_nij = k_n[itype][jtype];
  double k_naij = k_na[itype][jtype];
  std::vector<Contact> &contacts = s.contacts;

  if (npi == 1 && npj == 1) {
    sphere_against_sphere(i, j, delx, dely, delz, rsq, k_nij, k_naij, x, v, angmom, f, torque,
                          fnc, s, evdwl, facc);
    return;
  }

  contacts.clear();

  // check interaction between i's vertices and j' edges

  vertex_against_edge(i, j, 1, k_nij, k_naij, x, f, torque, s, evdwl, facc);

  // check interaction between j's vertices and i' edges
  // this returns the force on body j in fj, facc is the force on body i

  double fj[3] = {0.0, 0.0, 0.0};
  vertex_against_edge(j, i, 0, k_nij, k_naij, x, f, torque, s, evdwl, fj);

  int num_contacts = contacts.size();

  // the contact forces are applied at up to two contacts, see Fraige et al.:
  // the two contacts farthest apart, whose separation is the contact length
  // that scales the cohesive forces, or else the contact with the largest
  // overlap.  contacts between two vertices are treated like vertex-edge
  // contacts.  there is one friction force per pair of bodies, at the applied
  // contact with the largest overlap, or with the contact history from both
  // applied contacts, see history_friction()

  if (num_contacts > 0) {
    int m0 = 0, n0 = -1;
    double j_a = 1.0;

    for (int m = 1; m < num_contacts; m++)
      if (contacts[m].separation < contacts[m0].separation) m0 = m;

    double rmin = MIN(rounded_radius[i], rounded_radius[j]);
    double delta_max = EPSILON*rmin;
    for (int m = 0; m < num_contacts-1; m++) {
      for (int n = m+1; n < num_contacts; n++) {
        double delta_a = contact_separation(contacts[m], contacts[n]);
        if (delta_a > delta_max) {
          delta_max = delta_a;
          m0 = m;
          n0 = n;
        }
      }
    }
    if (n0 >= 0) j_a = MAX(delta_max / (EFF_CONTACTS * delta_ua), 1.0);

    int friction_m = 1;
    if ((n0 >= 0) && (contacts[n0].separation < contacts[m0].separation)) friction_m = 0;
    int friction_n = 1 - friction_m;
    if (history) friction_m = friction_n = 1;

    contact_forces(contacts[m0], j_a, friction_m, x, v, angmom, f, torque, fnc, evdwl,
                   (contacts[m0].ibody == i) ? facc : fj, s);
    if (n0 >= 0)
      contact_forces(contacts[n0], j_a, friction_n, x, v, angmom, f, torque, fnc, evdwl,
                     (contacts[n0].ibody == i) ? facc : fj, s);
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

void PairBodyRoundedPolygon::work_nonconservative()
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

int PairBodyRoundedPolygon::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++)
    for (int k = 0; k < NFNC; k++) buf[m++] = fnc[i][k];
  return m;
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolygon::unpack_reverse_comm(int n, int *list, double *buf)
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

void PairBodyRoundedPolygon::allocate()
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

void PairBodyRoundedPolygon::settings(int narg, char **arg)
{
  if (narg < 5) error->all(FLERR,"Illegal pair_style command");

  c_n = utils::numeric(FLERR,arg[0],false,lmp);
  c_t = utils::numeric(FLERR,arg[1],false,lmp);
  mu = utils::numeric(FLERR,arg[2],false,lmp);
  delta_ua = utils::numeric(FLERR,arg[3],false,lmp);
  cut_inner = utils::numeric(FLERR,arg[4],false,lmp);

  if (delta_ua < 0) delta_ua = 1;

  int history_one = 0;
  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"history") == 0) {
      history_one = 1;
      iarg++;
    } else error->all(FLERR, iarg, "Unknown pair_style body/rounded/polygon keyword {}",
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

void PairBodyRoundedPolygon::history_dummy_fix(int flag)
{
  std::string id_dummy = "NEIGH_HISTORY_BODY_DUMMY" + std::to_string(instance_me);
  if (flag) modify->add_fix(id_dummy + " all DUMMY");
  else if (modify->get_fix_by_id(id_dummy)) modify->delete_fix(id_dummy);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBodyRoundedPolygon::coeff(int narg, char **arg)
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
    error->all(FLERR, 4, "Tangential stiffness of pair style body/rounded/polygon "
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

void PairBodyRoundedPolygon::init_style()
{
  avec = dynamic_cast<AtomVecBody *>(atom->style_match("body"));
  if (!avec)
    error->all(FLERR,"Pair body/rounded/polygon requires atom style body");
  if (strcmp(avec->bptr->style,"rounded/polygon") != 0)
    error->all(FLERR,"Pair body/rounded/polygon requires body style rounded/polygon");
  bptr = dynamic_cast<BodyRoundedPolygon *>(avec->bptr);

  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style body/rounded/polygon requires newton pair on");

  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair body/rounded/polygon requires ghost atoms store velocity");

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
    memory->destroy(enclosing_radius);
    memory->destroy(rounded_radius);
    nmax = atom->nmax;
    memory->create(dnum,nmax,"pair:dnum");
    memory->create(dfirst,nmax,"pair:dfirst");
    memory->create(ednum,nmax,"pair:ednum");
    memory->create(edfirst,nmax,"pair:edfirst");
    memory->create(enclosing_radius,nmax,"pair:enclosing_radius");
    memory->create(rounded_radius,nmax,"pair:rounded_radius");
  }

  ndiscrete = nedge = 0;
  for (i = 0; i < nlocal; i++)
    dnum[i] = ednum[i] = 0;

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

double PairBodyRoundedPolygon::init_one(int i, int j)
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

void PairBodyRoundedPolygon::body2space(int i)
{
  int ibonus = atom->body[i];
  AtomVecBody::Bonus *bonus = &avec->bonus[ibonus];
  int nsub = bptr->nsub(bonus);
  double *coords = bptr->coords(bonus);
  int body_num_edges = bptr->nedges(bonus);
  double* edge_ends = bptr->edges(bonus);
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
  // the first 2 columns are for vertex indices within body, then the
  // outward unit normal of the edge in the xy plane, zero for rods and disks
  // or for an edge of zero length, and the squared length of the edge,
  // see nearest_point()

  if (nedge + body_num_edges > edmax) {
    edmax += DELTA;
    memory->grow(edge,edmax,6,"pair:edge");
  }

  if ((body_num_edges > 0) && (edge_ends == nullptr))
    error->one(FLERR,"Inconsistent edge data for body of atom {}", atom->tag[i]);

  for (int m = 0; m < body_num_edges; m++) {
    int np1 = static_cast<int>(edge_ends[2*m+0]);
    int np2 = static_cast<int>(edge_ends[2*m+1]);
    edge[nedge][0] = np1;
    edge[nedge][1] = np2;
    const double *a = discrete[dfirst[i]+np1];
    const double *b = discrete[dfirst[i]+np2];
    double v[3], en[3] = {0.0, 0.0, 0.0};
    MathExtra::sub3(b, a, v);
    double vsq = MathExtra::lensq3(v);
    if ((nsub > 2) && (vsq > 0.0)) {
      en[0] = v[1];
      en[1] = -v[0];
      if (en[0]*a[0] + en[1]*a[1] < 0.0) MathExtra::negate3(en);
      MathExtra::norm3(en);
    }
    edge[nedge][2] = en[0];
    edge[nedge][3] = en[1];
    edge[nedge][4] = en[2];
    edge[nedge][5] = vsq;
    nedge++;
  }

  enclosing_radius[i] = eradius;
  rounded_radius[i] = rradius;
}

/* ----------------------------------------------------------------------
   Normal force at the surface separation R, see Eq. 1, Fraige et al.:
     fe = elastic (repulsive) force, k_n * delta_n, for R < 0
     fc = cohesive (attractive) force, -k_na * delta_na, for R <= cut_inner,
          where delta_na = cut_inner - R is the overlap of the cohesive regions,
          which keeps growing when the surfaces deform
   return the energy, which is zero at R = cut_inner
------------------------------------------------------------------------- */

double PairBodyRoundedPolygon::normal_force(double R, double kn, double kna,
                                            double &fe, double &fc)
{
  fe = fc = 0.0;
  if (R > cut_inner) return 0.0;

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
   Interaction between two spheres with different radii
   according to the 2D model from Fraige et al.
---------------------------------------------------------------------- */

void PairBodyRoundedPolygon::sphere_against_sphere(int i, int j,
                       double delx, double dely, double delz, double rsq,
                       double k_n, double k_na, double** x, double** v,
                       double** angmom, double** f, double** torque,
                       double** fnc, Scratch &s, double &evdwl, double* facc)
{
  double rradi,rradj;
  double rij,R,fx,fy,fz,fpair,fe,fc,energy;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  rradi = rounded_radius[i];
  rradj = rounded_radius[j];

  rij = sqrt(rsq);
  R = rij - (rradi + rradj);

  energy = normal_force(R, k_n, k_na, fe, fc);
  fpair = fe + fc;

  fx = delx*fpair/rij;
  fy = dely*fpair/rij;
  fz = delz*fpair/rij;

  f[i][0] += fx;
  f[i][1] += fy;
  f[i][2] += fz;

  if (newton_pair || j < nlocal) {
    f[j][0] -= fx;
    f[j][1] -= fy;
    f[j][2] -= fz;
  }

  evdwl += energy;
  facc[0] += fx; facc[1] += fy; facc[2] += fz;

  // damping and friction at the contact point between the surfaces

  double rmin = MIN(rradi, rradj);
  if (R <= EPSILON*rmin) {
    double n[3] = {delx/rij, dely/rij, delz/rij};
    double pc[3];
    contact_point(x[i], x[j], n, rradi, rradj, pc);
    damping_friction(i, j, pc, n, fe, 1, 1, x, v, angmom, f, torque, fnc, i, facc, &s);
  }
}

/* ----------------------------------------------------------------------
   Interaction between the vertices of body i and the nearest features
   (edges or vertices) of body j

   i = atom i (body i)
   j = atom j (body j)
   first  = 1 for the first of the two calls for a pair of bodies, where the
            contacts between two vertices are counted, 0 for the second call
   x      = atoms' coordinates
   f      = atoms' forces
   torque = atoms' torques
   s      = scratch space, the contacts between i's vertices
            and j's edges are appended to s.contacts
   each vertex interacts with the nearest point on the core of body j, on an
   edge or at a vertex, so that the result does not depend on the shape of
   the polygon.  A vertex-vertex interaction is counted only if each of the
   two vertices is the nearest feature of its body to the other vertex (the
   normal cones of the two vertices overlap), otherwise the vertex-edge
   interaction from the other side covers it.
   the energy of the contacts is added in contact_forces(), since only
   the applied contacts contribute
   Return:
     interact = 0 no interaction at all
                1 there's at least one case where i's vertices interacts
                  with j's edges
---------------------------------------------------------------------- */

int PairBodyRoundedPolygon::vertex_against_edge(int i, int j, int first,
                                                double k_n, double k_na,
                                                double** x, double** f,
                                                double** torque, Scratch &s,
                                                double &evdwl, double* facc)
{
  std::vector<Contact> &contacts = s.contacts;

  int npi = dnum[i];
  int ifirst = dfirst[i];
  double rradi = rounded_radius[i];
  double eradj = enclosing_radius[j];
  double rradj = rounded_radius[j];
  double rmin = MIN(rradi, rradj);

  double xpi[3], h[3], n[3], hh[3], nn[3], dist, d, R, fe, fc, fx, fy, fz;
  double energy = 0.0;
  int nv, nvv, interact = 0;

  // loop through body i's vertices

  for (int ni = 0; ni < npi; ni++) {

    // convert body-fixed coordinates to space-fixed, xi

    xpi[0] = x[i][0] + discrete[ifirst+ni][0];
    xpi[1] = x[i][1] + discrete[ifirst+ni][1];
    xpi[2] = x[i][2] + discrete[ifirst+ni][2];

    // the vertex is within the enclosing circle (sphere) of body j,
    // possibly interacting

    distance(xpi, x[j], dist);
    if (dist > eradj + rradj + rradi + cut_inner) continue;

    // nearest point h on the core of body j, at the signed distance d

    d = nearest_point(j, xpi, h, n, nv, rradi + rradj + cut_inner);
    R = d - (rradi + rradj);
    if (R > cut_inner) continue;

    // vertex-vertex: both vertices must be the nearest feature of their body
    // to the other one, and the interaction is counted in the first call only

    if (nv >= 0) {
      if (!first) continue;
      if (npi > 1) {
        nearest_point(i, h, hh, nn, nvv);
        if (nvv != ni) continue;
      }
    }

    interact = 1;

    // R > rc:     no interaction between the vertex and body j
    // 0 < R < rc: cohesion between the vertex and body j
    // R < 0:      deformation between the vertex and body j
    // the normal damping term -c_n * vn will be added later

    double evertex = normal_force(R, k_n, k_na, fe, fc);

    if (R < EPSILON*rmin) {

      // the vertex of body i contacts body j: store the forces with the
      // contact, to be rescaled and applied later, see pair_interaction()

      Contact c;
      c.ibody = i;
      c.jbody = j;
      c.vertex = ni;
      c.jvertex = nv;
      for (int k = 0; k < 3; k++) {
        c.xv[k] = xpi[k];
        c.xe[k] = h[k];
        c.n[k] = n[k];
        c.fe[k] = n[k]*fe;
        c.fc[k] = n[k]*fc;
      }
      c.separation = R;
      c.energy = evertex;
      contacts.push_back(c);

    } else {

      // cohesion without contact: accumulate the force and torque to both
      // bodies directly

      energy += evertex;
      fx = n[0]*fc;
      fy = n[1]*fc;
      fz = n[2]*fc;

      f[i][0] += fx;
      f[i][1] += fy;
      f[i][2] += fz;
      sum_torque(x[i], xpi, fx, fy, fz, torque[i]);

      f[j][0] -= fx;
      f[j][1] -= fy;
      f[j][2] -= fz;
      sum_torque(x[j], h, -fx, -fy, -fz, torque[j]);

      facc[0] += fx; facc[1] += fy; facc[2] += fz;
    }
  }

  evdwl += energy;

  return interact;
}

/* ----------------------------------------------------------------------
  Find the point h on the core (the polygon without the rounded skin) of body
  ibody nearest to the point xp
  return the signed distance from h to xp, negative if xp is inside the core,
  with the unit normal n pointing from h towards the outside at xp, and
  the index of the vertex of ibody at h, or -1 if h is inside an edge
  a distance larger than dcut may be returned without h and n
  inside the core, h is on the edge with the largest signed distance, which
  is the nearest one for a convex polygon, and n is its outward normal
------------------------------------------------------------------------- */

double PairBodyRoundedPolygon::nearest_point(int ibody, const double *xp, double *h,
                                             double *n, int &nv, double dcut)
{
  double **x = atom->x;
  const double *xm = x[ibody];
  int ifirst = dfirst[ibody];
  int iefirst = edfirst[ibody];
  int nedges = ednum[ibody];
  double a[3], b[3], v[3], u[3], p[3];

  // a disk, or a body without edges: the nearest vertex

  nv = 0;
  for (int k = 0; k < 3; k++) h[k] = xm[k] + discrete[ifirst][k];
  double dmin = MathExtra::distsq3(xp, h);
  if (nedges == 0) {
    for (int m = 1; m < dnum[ibody]; m++) {
      for (int k = 0; k < 3; k++) p[k] = xm[k] + discrete[ifirst+m][k];
      double dsq = MathExtra::distsq3(xp, p);
      if (dsq < dmin) {
        dmin = dsq;
        nv = m;
        MathExtra::copy3(p, h);
      }
    }
  }

  // signed distances of xp from the lines of the edges, with their outward
  // normals in the xy plane, see body2space(), which are zero for rods and
  // for edges of zero length: xp is inside the core of a polygon, whose
  // center of mass is inside the core, if all of them are negative

  double smax = -BIG;
  int emax = -1;
  for (int ne = 0; ne < nedges; ne++) {
    const double *eg = edge[iefirst+ne];
    if (eg[5] == 0.0) continue;
    const double *da = discrete[ifirst+static_cast<int>(eg[0])];
    for (int k = 0; k < 3; k++) u[k] = xp[k] - (xm[k] + da[k]);
    double sd = MathExtra::dot3(u, &eg[2]);
    if ((emax < 0) || (sd > smax)) {
      smax = sd;
      emax = ne;
    }
  }

  // the distance from a polygon is at least the largest signed distance:
  // return it if it exceeds dcut, without h and n

  if ((dnum[ibody] > 2) && (emax >= 0) && (smax > dcut)) {
    nv = -1;
    return smax;
  }

  // xp is inside the core: push it out through the nearest edge

  if ((dnum[ibody] > 2) && (emax >= 0) && (smax < 0.0)) {
    const double *en = &edge[iefirst+emax][2];
    for (int k = 0; k < 3; k++) {
      h[k] = xp[k] - smax * en[k];
      n[k] = en[k];
    }
    nv = -1;
    return smax;
  }

  // the nearest point over the edges xp is in front of, i.e. with a
  // non-negative signed distance, which include the nearest edge or both
  // edges at the nearest vertex of a convex polygon

  int first = 1;
  for (int ne = 0; ne < nedges; ne++) {
    const double *eg = edge[iefirst+ne];
    int np1 = static_cast<int>(eg[0]);
    int np2 = static_cast<int>(eg[1]);
    for (int k = 0; k < 3; k++) {
      a[k] = xm[k] + discrete[ifirst+np1][k];
      u[k] = xp[k] - a[k];
    }
    if (MathExtra::dot3(u, &eg[2]) < 0.0) continue;
    for (int k = 0; k < 3; k++) {
      b[k] = xm[k] + discrete[ifirst+np2][k];
      v[k] = b[k] - a[k];
    }
    double t = (eg[5] > 0.0) ? MathExtra::dot3(u, v) / eg[5] : 0.0;
    int nvp = -1;
    if (t <= 0.0) {
      MathExtra::copy3(a, p);
      nvp = np1;
    } else if (t >= 1.0) {
      MathExtra::copy3(b, p);
      nvp = np2;
    } else {
      for (int k = 0; k < 3; k++) p[k] = a[k] + t * v[k];
    }
    double dsq = MathExtra::distsq3(xp, p);
    if (first || (dsq < dmin)) {
      first = 0;
      dmin = dsq;
      nv = nvp;
      MathExtra::copy3(p, h);
    }
  }

  double d = sqrt(dmin);
  if (d == 0.0)
    error->one(FLERR, "A vertex of a body touches the core of body {}", atom->tag[ibody]);
  for (int k = 0; k < 3; k++) n[k] = (xp[k] - h[k]) / d;
  return d;
}

/* ----------------------------------------------------------------------
  Compute contact forces between two bodies
  modify the force stored at the vertex and edge in contact by j_a
  sum forces and torque to the corresponding bodies
  fn = normal friction component
  ft = tangential friction component (-c_t * v_t)
------------------------------------------------------------------------- */

void PairBodyRoundedPolygon::contact_forces(Contact& contact, double j_a,
                       int friction, double** x, double** v, double** angmom, double** f,
                       double** torque, double** fnc, double &evdwl,
                       double* facc, Scratch &s)
{
  int ibody = contact.ibody;
  int jbody = contact.jbody;

  // unit normal from the point on the edge to the vertex, and the contact
  // point between the rounded surfaces, where all contact forces act

  const double *n = contact.n;
  double pc[3];
  contact_point(contact.xv, contact.xe, n, rounded_radius[ibody], rounded_radius[jbody], pc);
  evdwl += contact.energy;

  // elastic force and cohesive force, only the latter is scaled by j_a,
  // see Eq. 5, Fraige et al.

  double fx = contact.fe[0] + contact.fc[0] * j_a;
  double fy = contact.fe[1] + contact.fc[1] * j_a;
  double fz = contact.fe[2] + contact.fc[2] * j_a;
  f[ibody][0] += fx;
  f[ibody][1] += fy;
  f[ibody][2] += fz;
  sum_torque(x[ibody], pc, fx, fy, fz, torque[ibody]);
  f[jbody][0] -= fx;
  f[jbody][1] -= fy;
  f[jbody][2] -= fz;
  sum_torque(x[jbody], pc, -fx, -fy, -fz, torque[jbody]);
  facc[0] += fx; facc[1] += fy; facc[2] += fz;

  // the part of the contact force added by the j_a scaling does not derive
  // from the energy: accumulate it per body to compute its work,
  // see work_nonconservative()

  double fja[3];
  for (int k = 0; k < 3; k++) {
    fja[k] = (j_a - 1.0) * contact.fc[k];
    fnc[ibody][k] += fja[k];
    fnc[jbody][k] -= fja[k];
  }
  sum_torque(x[ibody], pc, fja[0], fja[1], fja[2], &fnc[ibody][3]);
  sum_torque(x[jbody], pc, -fja[0], -fja[1], -fja[2], &fnc[jbody][3]);

  // damping and, if requested, friction at the contact point

  double fne = MathExtra::len3(contact.fe);
  damping_friction(ibody, jbody, pc, n, fne, 1, friction, x, v, angmom, f, torque, fnc,
                   ibody, facc, &s);
}

/* ----------------------------------------------------------------------
  Contact point between the rounded surfaces of two bodies, halfway between
  the surface points on the line through the points pi on ibody and pj on
  jbody, where n is the unit normal pointing from jbody to ibody
------------------------------------------------------------------------- */

void PairBodyRoundedPolygon::contact_point(const double *pi, const double *pj,
                                           const double *n, double rradi, double rradj,
                                           double *pc)
{
  for (int k = 0; k < 3; k++)
    pc[k] = 0.5 * ((pi[k] - rradi * n[k]) + (pj[k] + rradj * n[k]));
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

void PairBodyRoundedPolygon::tangential_spring(int ibody, int jbody, const double *n,
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
  note: this deviates from Fraige et al., where the single friction force of
  a pair acts at the contact point with the largest solid overlap, with the
  relative velocity and normal of that contact, as also in Wang et al.,
  Granular Matter 13, 1 (2011).  The limit mu F_ne uses the elastic normal
  force of the whole interaction.  Without the contact history, the friction
  force still acts at the contact with the largest overlap.
  the total force on body i is accumulated to facc
------------------------------------------------------------------------- */

void PairBodyRoundedPolygon::history_friction(int i, int j, double **x, double **f,
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

void PairBodyRoundedPolygon::reset_dt()
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
  Damping and friction forces at the contact point pc between two bodies,
  from the relative velocity of the two bodies at that point:
    damping:  -c_n v_n - c_t v_t
    friction: magnitude mu * fne with the elastic normal force fne, opposite
              to v_t, capped by c_t * |v_t| so that it vanishes smoothly
              as sliding stops, see Eq. 4, Fraige et al.
              with the contact history, the friction force is that of a
              tangential spring instead, see tangential_spring()
  n = unit normal pointing from jbody to ibody
  the forces act at pc on both bodies, so that they exert torques
  the total force on body iref is accumulated to facc
------------------------------------------------------------------------- */

void PairBodyRoundedPolygon::damping_friction(int ibody, int jbody, double *pc,
  const double *n, double fne, int damping, int friction, double** x, double** v,
  double** angmom, double** f, double** torque, double** fnc, int iref, double* facc,
  Scratch *hs)
{
  double vi[3], vj[3], vr[3], vn[3], vt[3], fdiss[3];
  AtomVecBody::Bonus *bonus;

  bonus = &avec->bonus[atom->body[ibody]];
  total_velocity(pc, x[ibody], v[ibody], angmom[ibody], bonus->inertia, bonus->quat, vi);
  bonus = &avec->bonus[atom->body[jbody]];
  total_velocity(pc, x[jbody], v[jbody], angmom[jbody], bonus->inertia, bonus->quat, vj);

  MathExtra::sub3(vi, vj, vr);
  double vnnr = MathExtra::dot3(vr, n);
  for (int k = 0; k < 3; k++) {
    vn[k] = vnnr * n[k];
    vt[k] = vr[k] - vn[k];
    fdiss[k] = 0.0;
  }

  if (damping)
    for (int k = 0; k < 3; k++) fdiss[k] = -c_n * vn[k] - c_t * vt[k];

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
  sum_torque(x[ibody], pc, fdiss[0], fdiss[1], fdiss[2], torque[ibody]);

  f[jbody][0] -= fdiss[0];
  f[jbody][1] -= fdiss[1];
  f[jbody][2] -= fdiss[2];
  sum_torque(x[jbody], pc, -fdiss[0], -fdiss[1], -fdiss[2], torque[jbody]);

  // damping and friction do not derive from the energy

  for (int k = 0; k < 3; k++) {
    fnc[ibody][6+k] += fdiss[k];
    fnc[jbody][6+k] -= fdiss[k];
  }
  sum_torque(x[ibody], pc, fdiss[0], fdiss[1], fdiss[2], &fnc[ibody][9]);
  sum_torque(x[jbody], pc, -fdiss[0], -fdiss[1], -fdiss[2], &fnc[jbody][9]);

  if (ibody == iref) {
    facc[0] += fdiss[0]; facc[1] += fdiss[1]; facc[2] += fdiss[2];
  } else {
    facc[0] -= fdiss[0]; facc[1] -= fdiss[1]; facc[2] -= fdiss[2];
  }
}

/* ----------------------------------------------------------------------
  Determine the length of the contact segment, i.e. the separation between
  2 contacts along the tangent of the first one
------------------------------------------------------------------------- */

double PairBodyRoundedPolygon::contact_separation(const Contact& c1,
                                                  const Contact& c2)
{
  return fabs(-c1.n[1] * (c2.xv[0] - c1.xv[0]) + c1.n[0] * (c2.xv[1] - c1.xv[1]));
}

/* ----------------------------------------------------------------------
  Accumulate torque to body from the force f=(fx,fy,fz) acting at point x
------------------------------------------------------------------------- */

void PairBodyRoundedPolygon::sum_torque(double* xm, double *x, double fx,
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
  Calculate the total velocity of a point (vertex, a point on an edge):
    vi = vcm + omega ^ (p - xcm)
------------------------------------------------------------------------- */

void PairBodyRoundedPolygon::total_velocity(double* p, double *xcm,
                              double* vcm, double *angmom, double *inertia,
                              double *quat, double* vi)
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

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolygon::distance(const double* x2, const double* x1,
                                      double& r)
{
  r = sqrt((x2[0] - x1[0]) * (x2[0] - x1[0])
    + (x2[1] - x1[1]) * (x2[1] - x1[1])
    + (x2[2] - x1[2]) * (x2[2] - x1[2]));
}

/* ---------------------------------------------------------------------- */

double PairBodyRoundedPolygon::memory_usage()
{
  double bytes = Pair::memory_usage();
  bytes += (double) nmax * 4 * sizeof(int);    // dnum + dfirst + ednum + edfirst [nmax]
  return bytes;
}
