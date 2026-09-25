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
#include "fix_store_atom.h"
#include "force.h"
#include "group.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

static constexpr int DELTA = 10000;
static constexpr double EPSILON = 1.0e-3; // dimensionless threshold (dot products, end point checks, contact checks)
static constexpr int EFF_CONTACTS = 2;    // effective contacts for 2D models
static constexpr int NFNC = 12;           // per-body force and torque of the j_a scaling, and of damping
static constexpr char id_fix_store_prefix[] = "BODY_ROUNDED_POLYGON_WORK_";

//#define _CONVEX_POLYGON
//#define _POLYGON_DEBUG

enum { INVALID=0, NONE=1, VERTEXI=2, VERTEXJ=3, EDGE=4 };

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

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(k_n);
    memory->destroy(k_na);
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

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

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

      if (body[i] < 0 || body[j] < 0) continue;

      if (dnum[j] == 0) body2space(j);

      // no interaction

      r = sqrt(rsq);
      if (r > radi + radj + cut_inner) continue;

      pair_interaction(i, j, delx, dely, delz, rsq, x, v, angmom, f, torque, fnc,
                       scratch, evdwl, facc);

      if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,evdwl,0.0,
                               facc[0],facc[1],facc[2],delx,dely,delz);

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
  tagint *tag = atom->tag;
  int itype = atom->type[i];
  int jtype = atom->type[j];
  int npi = dnum[i];
  int npj = dnum[j];
  double k_nij = k_n[itype][jtype];
  double k_naij = k_na[itype][jtype];
  std::vector<Contact> &contacts = s.contacts;

  if (npi == 1 && npj == 1) {
    sphere_against_sphere(i, j, delx, dely, delz, rsq, k_nij, k_naij, v, f, fnc,
                          evdwl, facc);
    return;
  }

  int num_contacts, done;
  double delta_a, j_a;

  contacts.clear();

  // check interaction between i's vertices and j' edges

  vertex_against_edge(i, j, k_nij, k_naij, x, f, torque, tag, s, evdwl, facc);

  // check interaction between j's vertices and i' edges
  // this returns the force on body j in fj, facc is the force on body i

  double fj[3] = {0.0, 0.0, 0.0};
  vertex_against_edge(j, i, k_nij, k_naij, x, f, torque, tag, s, evdwl, fj);

  num_contacts = contacts.size();

  // the contact forces are applied at up to two contacts, see Fraige et al.:
  // the first two contacts at different places, whose separation is the
  // contact length that scales the cohesive forces, or else a single contact.
  // contacts between two vertices are treated like vertex-edge contacts.
  // there is one friction force per pair of bodies, at the applied contact
  // with the largest overlap

  if (num_contacts > 0) {
    int m0 = 0, n0 = -1;
    j_a = 1.0;

    done = 0;
    for (int m = 0; (m < num_contacts-1) && !done; m++) {
      for (int n = m+1; n < num_contacts; n++) {
        delta_a = contact_separation(contacts[m], contacts[n]);
        if (delta_a > 0) {
          m0 = m;
          n0 = n;
          j_a = delta_a / (EFF_CONTACTS * delta_ua);
          if (j_a < 1.0) j_a = 1.0;
          done = 1;
          break;
        }
      }
    }

    int friction_m = 1;
    if ((n0 >= 0) && (contacts[n0].separation < contacts[m0].separation)) friction_m = 0;

    contact_forces(contacts[m0], j_a, friction_m, x, v, angmom, f, torque, fnc, evdwl,
                   (contacts[m0].ibody == i) ? facc : fj);
    if (n0 >= 0)
      contact_forces(contacts[n0], j_a, 1 - friction_m, x, v, angmom, f, torque, fnc, evdwl,
                     (contacts[n0].ibody == i) ? facc : fj);

    #ifdef _POLYGON_DEBUG
    printf("  Contacts %d and %d: j_a = %f\n", m0, n0, j_a);
    #endif
  }

  #ifdef _POLYGON_DEBUG
  int num_overlapping_contacts = 0;
  for (int m = 0; m < num_contacts-1; m++) {
    for (int n = m+1; n < num_contacts; n++) {
      double l = contact_separation(contacts[m], contacts[n]);
      if (l < EPSILON) num_overlapping_contacts++;
    }
  }
  printf("There are %d contacts detected, %d of which overlap.\n",
         num_contacts, num_overlapping_contacts);
  #endif

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
  double halfdt = 0.5 * update->dt;
  double omega[3],ex[3],ey[3],ez[3];

  for (int i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;

    // no time step has been taken during setup

    if (!update->setupflag) {
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

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      k_n[i][j] = k_n_one;
      k_na[i][j] = k_na_one;
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

  neighbor->add_request(this);

  // per-atom storage of the forces and torques that do not derive
  // from the energy at the previous step, see work_nonconservative()

  if (!id_fix_store) {
    id_fix_store = utils::strdup(std::string(id_fix_store_prefix) + std::to_string(instance_me));
    fix_store = dynamic_cast<FixStoreAtom *>(
      modify->add_fix(fmt::format("{} {} STORE/ATOM {} 0 0 0", id_fix_store,
                                  group->names[0], NFNC)));
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
                       double k_n, double k_na, double** v, double** f,
                       double** fnc, double &evdwl, double* facc)
{
  double rradi,rradj;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double rij,rsqinv,R,fx,fy,fz,fn[3],ft[3],fpair,fe,fc,energy;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  rradi = rounded_radius[i];
  rradj = rounded_radius[j];

  rsqinv = 1.0/rsq;
  rij = sqrt(rsq);
  R = rij - (rradi + rradj);

  energy = normal_force(R, k_n, k_na, fe, fc);
  fpair = fe + fc;

  fx = delx*fpair/rij;
  fy = dely*fpair/rij;
  fz = delz*fpair/rij;

  double rmin = MIN(rradi, rradj);
  if (R <= EPSILON*rmin) { // in contact

    // relative translational velocity

    vr1 = v[i][0] - v[j][0];
    vr2 = v[i][1] - v[j][1];
    vr3 = v[i][2] - v[j][2];

    // normal component

    vnnr = vr1*delx + vr2*dely + vr3*delz;
    vn1 = delx*vnnr * rsqinv;
    vn2 = dely*vnnr * rsqinv;
    vn3 = delz*vnnr * rsqinv;

    // tangential component

    vt1 = vr1 - vn1;
    vt2 = vr2 - vn2;
    vt3 = vr3 - vn3;

    // normal friction term at contact

    fn[0] = -c_n * vn1;
    fn[1] = -c_n * vn2;
    fn[2] = -c_n * vn3;

    // tangential friction term at contact,
    // excluding the tangential deformation term for now

    ft[0] = -c_t * vt1;
    ft[1] = -c_t * vt2;
    ft[2] = -c_t * vt3;

    fx += fn[0] + ft[0];
    fy += fn[1] + ft[1];
    fz += fn[2] + ft[2];

    // damping does not derive from the energy, the disks do not rotate

    for (int k = 0; k < 3; k++) {
      fnc[i][6+k] += fn[k] + ft[k];
      fnc[j][6+k] -= fn[k] + ft[k];
    }
  }

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
}

/* ----------------------------------------------------------------------
   Determine the interaction mode between i's vertices against j's edges

   i = atom i (body i)
   j = atom j (body j)
   x      = atoms' coordinates
   f      = atoms' forces
   torque = atoms' torques
   tag    = atoms' tags
   s      = scratch space, the contacts between i's vertices
            and j's edges are appended to s.contacts
   Return:
     interact = 0 no interaction at all
                1 there's at least one case where i's vertices interacts
                  with j's edges
---------------------------------------------------------------------- */

int PairBodyRoundedPolygon::vertex_against_edge(int i, int j,
                                                double k_n, double k_na,
                                                double** x, double** f,
                                                double** torque, tagint* tag,
                                                Scratch &s,
                                                double &evdwl, double* facc)
{
  std::vector<Contact> &contacts = s.contacts;
  std::vector<int> &vertex_done = s.vertex_done;

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

  if ((int) vertex_done.size() < dnum[j]) vertex_done.resize(dnum[j]);

  // loop through body i's vertices

  for (ni = 0; ni < npi; ni++) {

    // convert body-fixed coordinates to space-fixed, xi

    xpi[0] = x[i][0] + discrete[ifirst+ni][0];
    xpi[1] = x[i][1] + discrete[ifirst+ni][1];
    xpi[2] = x[i][2] + discrete[ifirst+ni][2];

    // compute the distance from the vertex to the COM of body j

    distance(xpi, x[j], dist);

    #ifdef _POLYGON_DEBUG
    printf("Distance between vertex %d of body %d (%0.1f %0.1f %0.1f) "
           "to body %d's COM: %f (cut = %0.1f)\n",
           ni, xpi[0], xpi[1], xpi[2], atom->tag[i], atom->tag[j], dist,
           eradj + rradi + rradj + cut_inner);
    #endif

    // the vertex is within the enclosing circle (sphere) of body j,
    // possibly interacting

    if (dist > eradj + rradj + rradi + cut_inner) continue;

    // a vertex of body j is shared by two edges and both can report it

    for (int m = 0; m < dnum[j]; m++) vertex_done[m] = 0;

    int mode, contact, p2vertex;
    double d, R, hi[3], t, delx, dely, delz, fe, fc;
    double rij;
    double rmin = MIN(rradi, rradj);

    // the vertex interacts with the edge of body j whose sector, as seen
    // from the center of body j, encloses the vertex, see Fig. 4b in
    // Fraige et al., or with all edges if there is no such sector

    int nsector = sector_edge(j, xpi);

    // loop through body j's edges

    for (nj = 0; nj < nej; nj++) {

      if ((nsector >= 0) && (nj != nsector)) continue;

      // compute the distance between the edge nj to the vertex xpi

      mode = compute_distance_to_vertex(j, nj, x[j], rradj,
                                        xpi, rradi, cut_inner,
                                        d, hi, t, contact);

      if (mode == INVALID || mode == NONE) continue;

      if (mode == VERTEXI || mode == VERTEXJ) {

        interact = 1;

        // vertex i interacts with a vertex of the edge, but does not contact

        if (mode == VERTEXI) p2vertex = (int)edge[jefirst+nj][0];
        else p2vertex = (int)edge[jefirst+nj][1];

        // count the interaction with this vertex of body j only once

        if (vertex_done[p2vertex]) continue;
        vertex_done[p2vertex] = 1;

        // double xj[3];
        // p2.body2space(p2vertex, xj);
        xpj[0] = x[j][0] + discrete[jfirst+p2vertex][0];
        xpj[1] = x[j][1] + discrete[jfirst+p2vertex][1];
        xpj[2] = x[j][2] + discrete[jfirst+p2vertex][2];

        delx = xpi[0] - xpj[0];
        dely = xpi[1] - xpj[1];
        delz = xpi[2] - xpj[2];

        // R = surface separation = rij shifted by the rounded radii
        // R = rij - (p1.rounded_radius + p2.rounded_radius);
        // note: the force is defined for R, not for rij
        // R > rc:     no interaction between vertex ni and p2vertex
        // 0 < R < rc: cohesion between vertex ni and p2vertex
        // R < 0:      deformation between vertex ni and p2vertex

        rij = sqrt(delx*delx + dely*dely + delz*delz);
        R = rij - (rradi + rradj);

        // the normal damping term -c_n * vn will be added later

        double evertex = normal_force(R, k_n, k_na, fe, fc);

        fx = delx*(fe + fc)/rij;
        fy = dely*(fe + fc)/rij;
        fz = delz*(fe + fc)/rij;

        #ifdef _POLYGON_DEBUG
        printf("  Interaction between vertex %d of %d and vertex %d of %d:",
               ni, tag[i], p2vertex, tag[j]);
        printf("    mode = %d; contact = %d; d = %f; rij = %f, t = %f\n",
               mode, contact, d, rij, t);
        printf("    R = %f; cut_inner = %f\n", R, cut_inner);
        printf("    fe = %f; fc = %f\n", fe, fc);
        #endif

        // add forces to body i and body j directly
        // avoid double counts this pair of vertices
        // i and j can be either local or ghost atoms (bodies)
        // probably need more work here when the vertices' interaction
        // are not symmetric, e.g. j interacts with the edge
        // consisting of i but in mode = EDGE instead of VERTEX*.
        // OR, for the time being assume that the edge length is
        // sufficiently greater than the rounded radius to distinguish
        // vertex-vertex from vertex-edge contact modes.
        // Special case: when i is a sphere, also accumulate

        if (tag[i] < tag[j] || npi == 1) {

          energy += evertex;

          if (R < EPSILON*rmin) {

            // vertex ni of body i contacts vertex p2vertex of body j:
            // store the forces with the contact like for a vertex-edge contact,
            // so that damping and friction apply

            Contact c;
            c.ibody = i;
            c.jbody = j;
            c.vertex = ni;
            c.edge = -1;
            c.xv[0] = xpi[0];
            c.xv[1] = xpi[1];
            c.xv[2] = xpi[2];
            c.xe[0] = xpj[0];
            c.xe[1] = xpj[1];
            c.xe[2] = xpj[2];
            c.separation = R;
            c.fe[0] = delx*fe/rij;
            c.fe[1] = dely*fe/rij;
            c.fe[2] = delz*fe/rij;
            c.fc[0] = delx*fc/rij;
            c.fc[1] = dely*fc/rij;
            c.fc[2] = delz*fc/rij;
            contacts.push_back(c);
            continue;
          }

          f[i][0] += fx;
          f[i][1] += fy;
          f[i][2] += fz;
          sum_torque(x[i], xpi, fx, fy, fz, torque[i]);

          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          sum_torque(x[j], xpj, -fx, -fy, -fz, torque[j]);

          facc[0] += fx; facc[1] += fy; facc[2] += fz;

          #ifdef _POLYGON_DEBUG
          printf("    from vertex-vertex: "
                 "force on vertex %d of body %d: fx %f fy %f fz %f\n"
                 "      torque body %d: %f %f %f\n"
                 "      torque body %d: %f %f %f\n", ni, tag[i], fx, fy, fz,
            tag[i],torque[i][0],torque[i][1],torque[i][2],
            tag[j],torque[j][0],torque[j][1],torque[j][2]);
        #endif
        }

        #ifdef _CONVEX_POLYGON
        // done with the edges from body j,
        // given that vertex ni interacts with only one vertex
        //   from one edge of body j
        break;
        #endif

      } else if (mode == EDGE) {

        interact = 1;

        // vertex i interacts with the edge

        delx = xpi[0] - hi[0];
        dely = xpi[1] - hi[1];
        delz = xpi[2] - hi[2];

        // R = surface separation = d shifted by the rounded radii
        // R = d - (p1.rounded_radius + p2.rounded_radius);
        // Note: the force is defined for R, not for d
        // R > rc:     no interaction between vertex i and edge j
        // 0 < R < rc: cohesion between vertex i and edge j
        // R < 0:      deformation between vertex i and edge j
        // rij = sqrt(delx*delx + dely*dely + delz*delz);

        R = d - (rradi + rradj);

        // the normal damping term -c_n * vn will be added later

        energy += normal_force(R, k_n, k_na, fe, fc);

        fx = delx*(fe + fc)/d;
        fy = dely*(fe + fc)/d;
        fz = delz*(fe + fc)/d;

        #ifdef _POLYGON_DEBUG
        printf("  Interaction between vertex %d of %d and edge %d of %d:",
               ni, tag[i], nj, tag[j]);
        printf("    mode = %d; contact = %d; d = %f; t = %f\n",
               mode, contact, d, t);
        printf("    R = %f; cut_inner = %f\n", R, cut_inner);
        printf("    fe = %f; fc = %f\n", fe, fc);
        #endif

        if (contact == 1) {

          // vertex ni of body i contacts with edge nj of body j

          // store the force with the contact to be rescaled later
          // the force must be stored per contact, not per vertex or edge,
          // since several vertices of body i can contact the same edge

          Contact c;
          c.ibody = i;
          c.jbody = j;
          c.vertex = ni;
          c.edge = nj;
          c.xv[0] = xpi[0];
          c.xv[1] = xpi[1];
          c.xv[2] = xpi[2];
          c.xe[0] = hi[0];
          c.xe[1] = hi[1];
          c.xe[2] = hi[2];
          c.separation = R;
          c.fe[0] = delx*fe/d;
          c.fe[1] = dely*fe/d;
          c.fe[2] = delz*fe/d;
          c.fc[0] = delx*fc/d;
          c.fc[1] = dely*fc/d;
          c.fc[2] = delz*fc/d;
          contacts.push_back(c);

        } else { // no contact

          // accumulate force and torque to both bodies directly

          f[i][0] += fx;
          f[i][1] += fy;
          f[i][2] += fz;
          sum_torque(x[i], xpi, fx, fy, fz, torque[i]);

          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          sum_torque(x[j], hi, -fx, -fy, -fz, torque[j]);

          facc[0] += fx; facc[1] += fy; facc[2] += fz;

          #ifdef _POLYGON_DEBUG
          printf("    from vertex-edge, no contact: "
                 "force on vertex %d of body %d: fx %f fy %f fz %f\n"
                 "      torque body %d: %f %f %f\n"
                 "      torque body %d: %f %f %f\n", ni, tag[i], fx, fy, fz,
                 tag[i],torque[i][0],torque[i][1],torque[i][2],
                 tag[j],torque[j][0],torque[j][1],torque[j][2]);
          #endif
        } // end if contact

        #ifdef _CONVEX_POLYGON
        // done with the edges from body j,
        // given that vertex ni interacts with only one edge from body j
        break;
        #endif
      } // end if mode

    } // end for looping through the edges of body j

  } // end for looping through the vertices of body i

  evdwl += energy;

  return interact;
}

/* ----------------------------------------------------------------------
  Find the edge of body ibody whose sector encloses the point xp, where the
  sector of an edge is bounded by the rays from the center of the body
  through the two vertices of the edge, see Fig. 4b in Fraige et al.
  return the edge index, or -1 if there is no such edge,
    e.g. for rods and disks, or for a non-convex polygon
------------------------------------------------------------------------- */

int PairBodyRoundedPolygon::sector_edge(int ibody, const double *xp)
{
  if (dnum[ibody] < 3) return -1;

  double **x = atom->x;
  int ifirst = dfirst[ibody];
  int iefirst = edfirst[ibody];
  double cx = xp[0] - x[ibody][0];
  double cy = xp[1] - x[ibody][1];

  for (int ne = 0; ne < ednum[ibody]; ne++) {
    int np1 = static_cast<int>(edge[iefirst+ne][0]);
    int np2 = static_cast<int>(edge[iefirst+ne][1]);
    double ax = discrete[ifirst+np1][0];
    double ay = discrete[ifirst+np1][1];
    double bx = discrete[ifirst+np2][0];
    double by = discrete[ifirst+np2][1];
    double s0 = ax*by - ay*bx;
    if (s0 == 0.0) continue;
    double s1 = ax*cy - ay*cx;
    double s2 = cx*by - cy*bx;
    if ((s1*s0 >= 0.0) && (s2*s0 >= 0.0)) return ne;
  }
  return -1;
}

/* -------------------------------------------------------------------------
  Compute the distance between an edge of body i and a vertex from
  another body
  Input:
    ibody      = body i (i.e. atom i)
    edge_index = edge index of body i
    xmi        = atom i's coordinates (body i's center of mass)
    x0         = coordinate of the tested vertex from another body
    x0_rounded_radius = rounded radius of the tested vertex
    cut_inner  = cutoff for vertex-vertex and vertex-edge interaction
  Output:
    d          = Distance from a point x0 to an edge
    hi         = coordinates of the projection of x0 on the edge
    t          = ratio to determine the relative position of hi
                 wrt xi and xj on the segment
  contact      = 0 no contact between the queried vertex and the edge
                 1 contact detected
  return
    INVALID if the edge index is invalid
    NONE    if there is no interaction
    VERTEXI if the tested vertex interacts with the first vertex of the edge
    VERTEXJ if the tested vertex interacts with the second vertex of the edge
    EDGE    if the tested vertex interacts with the edge
------------------------------------------------------------------------- */

int PairBodyRoundedPolygon::compute_distance_to_vertex(int ibody,
                                                int edge_index,
                                                double *xmi,
                                                double rounded_radius,
                                                double* x0,
                                                double x0_rounded_radius,
                                                double cut_inner,
                                                double &d,
                                                double hi[3],
                                                double &t,
                                                int &contact)
{
  if (edge_index >= ednum[ibody]) return INVALID;

  int mode,ifirst,iefirst,npi1,npi2;
  double xi1[3],xi2[3],u[3],v[3],uij[3];
  double udotv, magv, magucostheta;
  double delx,dely,delz;
  double rmin = MIN(rounded_radius, x0_rounded_radius);

  ifirst = dfirst[ibody];
  iefirst = edfirst[ibody];
  npi1 = static_cast<int>(edge[iefirst+edge_index][0]);
  npi2 = static_cast<int>(edge[iefirst+edge_index][1]);

  // compute the space-fixed coordinates for the vertices of the edge

  xi1[0] = xmi[0] + discrete[ifirst+npi1][0];
  xi1[1] = xmi[1] + discrete[ifirst+npi1][1];
  xi1[2] = xmi[2] + discrete[ifirst+npi1][2];

  xi2[0] = xmi[0] + discrete[ifirst+npi2][0];
  xi2[1] = xmi[1] + discrete[ifirst+npi2][1];
  xi2[2] = xmi[2] + discrete[ifirst+npi2][2];

  // u = x0 - xi1

  u[0] = x0[0] - xi1[0];
  u[1] = x0[1] - xi1[1];
  u[2] = x0[2] - xi1[2];

  // v = xi2 - xi1

  v[0] = xi2[0] - xi1[0];
  v[1] = xi2[1] - xi1[1];
  v[2] = xi2[2] - xi1[2];

  // dot product between u and v = magu * magv * costheta

  udotv = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  magv = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  magucostheta = udotv / magv;

  // uij is the unit vector pointing from xi to xj

  uij[0] = v[0] / magv;
  uij[1] = v[1] / magv;
  uij[2] = v[2] / magv;

  // position of the projection of x0 on the line (xi, xj)

  hi[0] = xi1[0] + magucostheta * uij[0];
  hi[1] = xi1[1] + magucostheta * uij[1];
  hi[2] = xi1[2] + magucostheta * uij[2];

  // distance from x0 to the line (xi, xj) = distance from x0 to hi

  distance(hi, x0, d);

  // determine the interaction mode
  // for 2D: a vertex can interact with one edge at most
  // for 3D: a vertex can interact with one face at most

  mode = NONE;
  contact = 0;

  if (d > rounded_radius + x0_rounded_radius + cut_inner) {

    // if the vertex is far away from the edge

    mode = NONE;

  } else {

    // check if x0 (the queried vertex) and xmi (the body's center of mass)
    // are on the different sides of the edge

    #ifdef _CONVEX_POLYGON
    int m = opposite_sides(xi1, xi2, x0, xmi);
    #else
    int m = 1;
    #endif

    if (m == 0) {

      // x0 and xmi are on not the opposite sides of the edge
      // leave xpi for another edge to detect

      mode = NONE;

    } else {

      // x0 and xmi are on the different sides
      // t is the ratio to detect if x0 is closer to the vertices xi or xj

      t = (magv > 0.0) ? magucostheta / magv : 0.0;

      double contact_dist = rounded_radius + x0_rounded_radius;
      if (t >= 0 && t <= 1) {
        mode = EDGE;
        if (d < contact_dist + EPSILON*rmin)
          contact = 1;

      } else { // t < 0 || t > 1: closer to either vertices of the edge

        if (t < 0) {
          // measure the distance from x0 to xi1
          delx = x0[0] - xi1[0];
          dely = x0[1] - xi1[1];
          delz = x0[2] - xi1[2];
          double dx0xi1 = sqrt(delx*delx + dely*dely + delz*delz);
          if (dx0xi1 > contact_dist + cut_inner)
            mode = NONE;
          else
            mode = VERTEXI;
        } else {
          // measure the distance from x0 to xi2
          delx = x0[0] - xi2[0];
          dely = x0[1] - xi2[1];
          delz = x0[2] - xi2[2];
          double dx0xi2 = sqrt(delx*delx + dely*dely + delz*delz);
          if (dx0xi2 > contact_dist + cut_inner)
            mode = NONE;
          else
            mode = VERTEXJ;
        }
      } // end if t >= 0 && t <= 1
    } // end if x0 and xmi are on the same side of the edge
  }

  return mode;
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
                       double** torque, double** fnc, double &/*evdwl*/,
                       double* facc)
{
  int ibody,jbody,ibonus,jbonus;
  double fx,fy,fz,delx,dely,delz,rsq,rsqinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double fn[3],ft[3],vi[3],vj[3];
  double *quat, *inertia;
  AtomVecBody::Bonus *bonus;

  ibody = contact.ibody;
  jbody = contact.jbody;

  // compute the velocity of the vertex in the space-fixed frame

  ibonus = atom->body[ibody];
  bonus = &avec->bonus[ibonus];
  quat = bonus->quat;
  inertia = bonus->inertia;
  total_velocity(contact.xv, x[ibody], v[ibody], angmom[ibody],
                 inertia, quat, vi);

  // compute the velocity of the point on the edge in the space-fixed frame

  jbonus = atom->body[jbody];
  bonus = &avec->bonus[jbonus];
  quat = bonus->quat;
  inertia = bonus->inertia;
  total_velocity(contact.xe, x[jbody], v[jbody], angmom[jbody],
                 inertia, quat, vj);

  // vector pointing from the vertex to the point on the edge

  delx = contact.xv[0] - contact.xe[0];
  dely = contact.xv[1] - contact.xe[1];
  delz = contact.xv[2] - contact.xe[2];
  rsq = delx*delx + dely*dely + delz*delz;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = vi[0] - vj[0];
  vr2 = vi[1] - vj[1];
  vr3 = vi[2] - vj[2];

  // normal component

  vnnr = vr1*delx + vr2*dely + vr3*delz;
  vn1 = delx*vnnr * rsqinv;
  vn2 = dely*vnnr * rsqinv;
  vn3 = delz*vnnr * rsqinv;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // normal friction term at contact

  fn[0] = -c_n * vn1;
  fn[1] = -c_n * vn2;
  fn[2] = -c_n * vn3;

  // tangential friction term at contact
  // excluding the tangential deformation term for now

  ft[0] = -c_t * vt1;
  ft[1] = -c_t * vt2;
  ft[2] = -c_t * vt3;

  // kinetic friction during gross sliding, see Eq. 4, Fraige et al.:
  // magnitude mu * |F_ne|, opposite to the tangential relative velocity,
  // capped by c_t * |v_t| so that it vanishes smoothly as sliding stops
  // instead of reversing the sliding direction within a time step

  // there is one friction force per pair of bodies, at the contact with
  // the largest overlap, and F_ne is the elastic normal force only

  double vtmag = sqrt(vt1*vt1 + vt2*vt2 + vt3*vt3);
  if (friction && (vtmag > 0.0)) {
    double fne = sqrt(contact.fe[0]*contact.fe[0] + contact.fe[1]*contact.fe[1] +
                      contact.fe[2]*contact.fe[2]);
    double scale = MIN(mu * fne, c_t * vtmag) / vtmag;
    ft[0] -= scale * vt1;
    ft[1] -= scale * vt2;
    ft[2] -= scale * vt3;
  }

  // the part of the contact force added by the j_a scaling, and damping and
  // friction do not derive from the energy: accumulate them per body
  // to compute their work, see work_nonconservative()

  double fja[3], fdiss[3];
  for (int k = 0; k < 3; k++) {
    fja[k] = (j_a - 1.0) * contact.fc[k];
    fdiss[k] = fn[k] + ft[k];
    fnc[ibody][k] += fja[k];
    fnc[jbody][k] -= fja[k];
    fnc[ibody][6+k] += fdiss[k];
    fnc[jbody][6+k] -= fdiss[k];
  }
  sum_torque(x[ibody], contact.xv, fja[0], fja[1], fja[2], &fnc[ibody][3]);
  sum_torque(x[jbody], contact.xe, -fja[0], -fja[1], -fja[2], &fnc[jbody][3]);
  sum_torque(x[ibody], contact.xv, fdiss[0], fdiss[1], fdiss[2], &fnc[ibody][9]);
  sum_torque(x[jbody], contact.xe, -fdiss[0], -fdiss[1], -fdiss[2], &fnc[jbody][9]);

  // only the cohesive force is scaled by j_a, see Eq. 5, Fraige et al.

  fx = contact.fe[0] + contact.fc[0] * j_a + fn[0] + ft[0];
  fy = contact.fe[1] + contact.fc[1] * j_a + fn[1] + ft[1];
  fz = contact.fe[2] + contact.fc[2] * j_a + fn[2] + ft[2];
  f[ibody][0] += fx;
  f[ibody][1] += fy;
  f[ibody][2] += fz;
  sum_torque(x[ibody], contact.xv, fx, fy, fz, torque[ibody]);

  // accumulate forces to the vertex only

  facc[0] += fx; facc[1] += fy; facc[2] += fz;

  fx = -contact.fe[0] - contact.fc[0] * j_a - fn[0] - ft[0];
  fy = -contact.fe[1] - contact.fc[1] * j_a - fn[1] - ft[1];
  fz = -contact.fe[2] - contact.fc[2] * j_a - fn[2] - ft[2];
  f[jbody][0] += fx;
  f[jbody][1] += fy;
  f[jbody][2] += fz;
  sum_torque(x[jbody], contact.xe, fx, fy, fz, torque[jbody]);

  #ifdef _POLYGON_DEBUG
  printf("From contact forces: vertex fx %f fy %f fz %f\n"
         "      torque body %d: %f %f %f\n"
         "      torque body %d: %f %f %f\n",
         contact.fe[0], contact.fe[1], contact.fe[2],
         atom->tag[ibody],torque[ibody][0],torque[ibody][1],torque[ibody][2],
         atom->tag[jbody],torque[jbody][0],torque[jbody][1],torque[jbody][2]);
  #endif
}

/* ----------------------------------------------------------------------
  Determine the length of the contact segment, i.e. the separation between
  2 contacts, should be extended for 3D models.
------------------------------------------------------------------------- */

double PairBodyRoundedPolygon::contact_separation(const Contact& c1,
                                                  const Contact& c2)
{
  double x1 = c1.xv[0];
  double y1 = c1.xv[1];
  double x2 = c1.xe[0];
  double y2 = c1.xe[1];
  double x3 = c2.xv[0];
  double y3 = c2.xv[1];

  int ibody = c1.ibody;
  int jbody = c1.ibody;
  double rradi = rounded_radius[ibody];
  double rradj = rounded_radius[jbody];
  double rmin = MIN(rradi, rradj);

  double delta_a = 0.0;
  if (fabs(x2 - x1) > EPSILON*rmin) {
    double A = (y2 - y1) / (x2 - x1);
    delta_a = fabs(y1 - A * x1 - y3 + A * x3) / sqrt(1 + A * A);
  } else {
    delta_a = fabs(x1 - x3);
  }

  return delta_a;
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
  Test if two points a and b are in opposite sides of the line that
  connects two points x1 and x2
------------------------------------------------------------------------- */

int PairBodyRoundedPolygon::opposite_sides(double* x1, double* x2,
                                           double* a, double* b)
{
  double m_a = (x1[1] - x2[1])*(a[0] - x1[0]) + (x2[0] - x1[0])*(a[1] - x1[1]);
  double m_b = (x1[1] - x2[1])*(b[0] - x1[0]) + (x2[0] - x1[0])*(b[1] - x1[1]);
  // equal to zero when either a or b is inline with the line x1-x2
  if (m_a * m_b <= 0)
    return 1;
  else
    return 0;
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
