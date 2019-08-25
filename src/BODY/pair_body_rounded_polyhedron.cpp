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

/* ----------------------------------------------------------------------
   Contributing author: Trung Dac Nguyen (ndactrung@gmail.com)
   Ref: Wang, Yu, Langston, Fraige, Particle shape effects in discrete
   element modelling of cohesive angular particles, Granular Matter 2011,
   13:1-12.
   Note: The current implementation has not taken into account
         the contact history for friction forces.
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_body_rounded_polyhedron.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_body.h"
#include "body_rounded_polyhedron.h"
#include "comm.h"
#include "force.h"
#include "fix.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "math_extra.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathExtra;
using namespace MathConst;

#define DELTA 10000
#define EPSILON 1e-3
#define MAX_FACE_SIZE 4  // maximum number of vertices per face (same as BodyRoundedPolyhedron)
#define MAX_CONTACTS 32  // for 3D models (including duplicated counts)

//#define _POLYHEDRON_DEBUG

enum {EE_INVALID=0,EE_NONE,EE_INTERACT};
enum {EF_INVALID=0,EF_NONE,EF_PARALLEL,EF_SAME_SIDE_OF_FACE,
      EF_INTERSECT_INSIDE,EF_INTERSECT_OUTSIDE};

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolyhedron::PairBodyRoundedPolyhedron(LAMMPS *lmp) : Pair(lmp)
{
  dmax = nmax = 0;
  discrete = NULL;
  dnum = dfirst = NULL;

  edmax = ednummax = 0;
  edge = NULL;
  ednum = edfirst = NULL;

  facmax = facnummax = 0;
  face = NULL;
  facnum = facfirst = NULL;

  enclosing_radius = NULL;
  rounded_radius = NULL;
  maxerad = NULL;

  single_enable = 0;
  restartinfo = 0;

  c_n = 0.1;
  c_t = 0.2;
  mu = 0.0;
  A_ua = 1.0;

  k_n = NULL;
  k_na = NULL;
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
  memory->destroy(facnum);
  memory->destroy(facfirst);

  memory->destroy(enclosing_radius);
  memory->destroy(rounded_radius);
  memory->destroy(maxerad);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(k_n);
    memory->destroy(k_na);
  }
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  int ni,nj,npi,npj,ifirst,jfirst,nei,nej,iefirst,jefirst;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,facc[3];
  double rsq,eradi,eradj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **torque = atom->torque;
  double **angmom = atom->angmom;
  int *body = atom->body;
  int *type = atom->type;
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

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (body[i] >= 0) {
      if (dnum[i] == 0) body2space(i);
      npi = dnum[i];
      ifirst = dfirst[i];
      nei = ednum[i];
      iefirst = edfirst[i];
      eradi = enclosing_radius[i];
     }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      // body/body interactions

      evdwl = 0.0;
      facc[0] = facc[1] = facc[2] = 0;

      if (body[i] < 0 || body[j] < 0) continue;

      if (dnum[j] == 0) body2space(j);
      npj = dnum[j];
      jfirst = dfirst[j];
      nej = ednum[j];
      jefirst = edfirst[j];
      eradj = enclosing_radius[j];

      // no interaction

      double r = sqrt(rsq);
      if (r > eradi + eradj + cut_inner) continue;

      // sphere-sphere interaction

      if (npi == 1 && npj == 1) {
        sphere_against_sphere(i, j, itype, jtype, delx, dely, delz,
                              rsq, v, f, evflag);
        continue;
      }

      // reset vertex and edge forces

      for (ni = 0; ni < npi; ni++) {
        discrete[ifirst+ni][3] = 0;
        discrete[ifirst+ni][4] = 0;
        discrete[ifirst+ni][5] = 0;
        discrete[ifirst+ni][6] = 0;
      }

      for (nj = 0; nj < npj; nj++) {
        discrete[jfirst+nj][3] = 0;
        discrete[jfirst+nj][4] = 0;
        discrete[jfirst+nj][5] = 0;
        discrete[jfirst+nj][6] = 0;
      }

      for (ni = 0; ni < nei; ni++) {
        edge[iefirst+ni][2] = 0;
        edge[iefirst+ni][3] = 0;
        edge[iefirst+ni][4] = 0;
        edge[iefirst+ni][5] = 0;
      }

      for (nj = 0; nj < nej; nj++) {
        edge[jefirst+nj][2] = 0;
        edge[jefirst+nj][3] = 0;
        edge[jefirst+nj][4] = 0;
        edge[jefirst+nj][5] = 0;
      }

      // one of the two bodies is a sphere

      if (npj == 1) {
        sphere_against_face(i, j, itype, jtype, x, v, f, torque,
                            angmom, evflag);
        sphere_against_edge(i, j, itype, jtype, x, v, f, torque,
                            angmom, evflag);
        continue;
      } else if (npi == 1) {
        sphere_against_face(j, i, jtype, itype, x, v, f, torque,
                            angmom, evflag);
        sphere_against_edge(j, i, jtype, itype, x, v, f, torque,
                            angmom, evflag);
        continue;
      }

      int interact, num_contacts;
      Contact contact_list[MAX_CONTACTS];

      num_contacts = 0;

      // check interaction between i's edges and j' faces
      #ifdef _POLYHEDRON_DEBUG
      printf("INTERACTION between edges of %d vs. faces of %d:\n", i, j);
      #endif
      interact = edge_against_face(i, j, itype, jtype, x, contact_list,
                                   num_contacts, evdwl, facc);

      // check interaction between j's edges and i' faces
      #ifdef _POLYHEDRON_DEBUG
      printf("\nINTERACTION between edges of %d vs. faces of %d:\n", j, i);
      #endif
      interact = edge_against_face(j, i, jtype, itype, x, contact_list,
                                   num_contacts, evdwl, facc);

      // check interaction between i's edges and j' edges
      #ifdef _POLYHEDRON_DEBUG
      printf("INTERACTION between edges of %d vs. edges of %d:\n", i, j);
      #endif
      interact = edge_against_edge(i, j, itype, jtype, x, contact_list,
                                   num_contacts, evdwl, facc);

      // estimate the contact area
      // also consider point contacts and line contacts

      if (num_contacts > 0) {
        rescale_cohesive_forces(x, f, torque, contact_list, num_contacts,
                                itype, jtype, facc);
      }

      if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,evdwl,0.0,
                               facc[0],facc[1],facc[2],delx,dely,delz);

    } // end for jj
  }

  if (vflag_fdotr) virial_fdotr_compute();
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
  memory->create(maxerad,n+1,"pair:maxerad");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::settings(int narg, char **arg)
{
  if (narg < 5) error->all(FLERR,"Illegal pair_style command");

  c_n = force->numeric(FLERR,arg[0]);
  c_t = force->numeric(FLERR,arg[1]);
  mu = force->numeric(FLERR,arg[2]);
  A_ua = force->numeric(FLERR,arg[3]);
  cut_inner = force->numeric(FLERR,arg[4]);

  if (A_ua < 0) A_ua = 1;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double k_n_one = force->numeric(FLERR,arg[2]);
  double k_na_one = force->numeric(FLERR,arg[3]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      k_n[i][j] = k_n_one;
      k_na[i][j] = k_na_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::init_style()
{
  avec = (AtomVecBody *) atom->style_match("body");
  if (!avec) error->all(FLERR,"Pair body/rounded/polyhedron requires "
                        "atom style body");
  if (strcmp(avec->bptr->style,"rounded/polyhedron") != 0)
    error->all(FLERR,"Pair body/rounded/polyhedron requires "
               "body style rounded/polyhedron");
  bptr = (BodyRoundedPolyhedron *) avec->bptr;

  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style body/rounded/polyhedron requires "
               "newton pair on");

  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair body/rounded/polyhedron requires "
               "ghost atoms store velocity");

  neighbor->request(this);

  // find the maximum enclosing radius for each atom type

  int i, itype;
  double eradi;
  int* body = atom->body;
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

  double *merad = NULL;
  memory->create(merad,ntypes+1,"pair:merad");
  for (i = 1; i <= ntypes; i++)
    maxerad[i] = merad[i] = 0;

  int ipour;
  for (ipour = 0; ipour < modify->nfix; ipour++)
    if (strcmp(modify->fix[ipour]->style,"pour") == 0) break;
  if (ipour == modify->nfix) ipour = -1;

  int idep;
  for (idep = 0; idep < modify->nfix; idep++)
    if (strcmp(modify->fix[idep]->style,"deposit") == 0) break;
  if (idep == modify->nfix) idep = -1;

  for (i = 1; i <= ntypes; i++) {
    merad[i] = 0.0;
    if (ipour >= 0) {
      itype = i;
      merad[i] =
        *((double *) modify->fix[ipour]->extract("radius",itype));
    }
    if (idep >= 0) {
      itype = i;
      merad[i] =
        *((double *) modify->fix[idep]->extract("radius",itype));
    }
  }

  for (i = 0; i < nlocal; i++) {
    itype = type[i];
    if (body[i] >= 0) {
      if (dnum[i] == 0) body2space(i);
      eradi = enclosing_radius[i];
      if (eradi > merad[itype]) merad[itype] = eradi;
    } else
      merad[itype] = 0;
  }

  MPI_Allreduce(&merad[1],&maxerad[1],ntypes,MPI_DOUBLE,MPI_MAX,world);

  memory->destroy(merad);

  sanity_check();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBodyRoundedPolyhedron::init_one(int i, int j)
{
  k_n[j][i] = k_n[i][j];
  k_na[j][i] = k_na[i][j];

  return (maxerad[i]+maxerad[j]);
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
  // the first 3 columns are for coords, the last 3 for forces

  if (ndiscrete + nsub > dmax) {
    dmax += DELTA;
    memory->grow(discrete,dmax,7,"pair:discrete");
  }

  double p[3][3];
  MathExtra::quat_to_mat(bonus->quat,p);

  for (int m = 0; m < nsub; m++) {
    MathExtra::matvec(p,&coords[3*m],discrete[ndiscrete]);
    discrete[ndiscrete][3] = 0;
    discrete[ndiscrete][4] = 0;
    discrete[ndiscrete][5] = 0;
    discrete[ndiscrete][6] = 0;
    ndiscrete++;
  }

  // get the number of edges (vertices)
  // and the index of the first edge of my body in the list

  ednum[i] = body_num_edges;
  edfirst[i] = nedge;

  // grow the edge list if necessary
  // the first 2 columns are for vertex indices within body, the last 3 for forces

  if (nedge + body_num_edges > edmax) {
    edmax += DELTA;
    memory->grow(edge,edmax,6,"pair:edge");
  }

  for (int m = 0; m < body_num_edges; m++) {
    edge[nedge][0] = static_cast<int>(edge_ends[2*m+0]);
    edge[nedge][1] = static_cast<int>(edge_ends[2*m+1]);
    edge[nedge][2] = 0;
    edge[nedge][3] = 0;
    edge[nedge][4] = 0;
    edge[nedge][5] = 0;
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
  }

  for (int m = 0; m < body_num_faces; m++) {
    for (int k = 0; k < MAX_FACE_SIZE; k++)
      face[nface][k] = static_cast<int>(face_pts[MAX_FACE_SIZE*m+k]);
    nface++;
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
  double** v, double** f, int evflag)
{
  double rradi,rradj,contact_dist;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double rij,rsqinv,R,fx,fy,fz,fn[3],ft[3],fpair,energy;
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

  if (R <= 0) { // in contact

    // relative translational velocity

    vr1 = v[ibody][0] - v[jbody][0];
    vr2 = v[ibody][1] - v[jbody][1];
    vr3 = v[ibody][2] - v[jbody][2];

    // normal component

    rsqinv = 1.0/rsq;
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
  }

  f[ibody][0] += fx;
  f[ibody][1] += fy;
  f[ibody][2] += fz;

  if (newton_pair || jbody < nlocal) {
    f[jbody][0] -= fx;
    f[jbody][1] -= fy;
    f[jbody][2] -= fz;
  }

  if (evflag) ev_tally_xyz(ibody,jbody,nlocal,newton_pair,
                           energy,0.0,fx,fy,fz,delx,dely,delz);
}

/* ----------------------------------------------------------------------
   Interaction bt the edges of a polyhedron (ibody) and a sphere (jbody)
---------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::sphere_against_edge(int ibody, int jbody,
  int itype, int jtype, double** x, double** v, double** f, double** torque,
  double** angmom, int evflag)
{
  int ni,nei,ifirst,iefirst,npi1,npi2,ibonus;
  double xi1[3],xi2[3],vti[3],h[3],fn[3],ft[3],d,t;
  double delx,dely,delz,rsq,rij,rsqinv,R,fx,fy,fz,fpair,energy;
  double rradi,rradj,contact_dist;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
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

    // compute the space-fixed coordinates for the vertices of the face

    xi1[0] = x[ibody][0] + discrete[ifirst+npi1][0];
    xi1[1] = x[ibody][1] + discrete[ifirst+npi1][1];
    xi1[2] = x[ibody][2] + discrete[ifirst+npi1][2];

    xi2[0] = x[ibody][0] + discrete[ifirst+npi2][0];
    xi2[1] = x[ibody][1] + discrete[ifirst+npi2][1];
    xi2[2] = x[ibody][2] + discrete[ifirst+npi2][2];

    // find the projection of the jbody's COM on the edge

    project_pt_line(x[jbody], xi1, xi2, h, d, t);

    if (d > contact_dist + cut_inner) continue;
    if (t < 0 || t > 1) continue;

    if (fabs(t) < EPSILON) {
      if (static_cast<int>(discrete[ifirst+npi1][6]) == 1)
        continue;
      else {
        h[0] = xi1[0];
        h[1] = xi1[1];
        h[2] = xi1[2];
        discrete[ifirst+npi1][6] = 1;
      }
    }

    if (fabs(t-1) < EPSILON) {
      if (static_cast<int>(discrete[ifirst+npi2][6]) == 1)
        continue;
      else {
        h[0] = xi2[0];
        h[1] = xi2[1];
        h[2] = xi2[2];
        discrete[ifirst+npi2][6] = 1;
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

    if (R <= 0) { // in contact

      // compute the velocity of the vertex in the space-fixed frame

      ibonus = atom->body[ibody];
      bonus = &avec->bonus[ibonus];
      quat = bonus->quat;
      inertia = bonus->inertia;
      total_velocity(h, x[ibody], v[ibody], angmom[ibody],
                     inertia, quat, vti);

      // relative translational velocity

      vr1 = vti[0] - v[jbody][0];
      vr2 = vti[1] - v[jbody][1];
      vr3 = vti[2] - v[jbody][2];

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
      // excluding the tangential deformation term

      ft[0] = -c_t * vt1;
      ft[1] = -c_t * vt2;
      ft[2] = -c_t * vt3;

      fx += fn[0] + ft[0];
      fy += fn[1] + ft[1];
      fz += fn[2] + ft[2];
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

    if (evflag) ev_tally_xyz(ibody,jbody,nlocal,newton_pair,
                           energy,0.0,fx,fy,fz,delx,dely,delz);
  }
}

/* ----------------------------------------------------------------------
   Interaction bt the faces of a polyhedron (ibody) and a sphere (jbody)
---------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::sphere_against_face(int ibody, int jbody,
 int itype, int jtype, double** x, double** v, double** f, double** torque,
 double** angmom, int evflag)
{
  int ni,nfi,inside,ifirst,iffirst,npi1,npi2,npi3,ibonus,tmp;
  double xi1[3],xi2[3],xi3[3],ui[3],vi[3],vti[3],n[3],h[3],fn[3],ft[3],d;
  double delx,dely,delz,rsq,rij,rsqinv,R,fx,fy,fz,fpair,energy;
  double rradi,rradj,contact_dist;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
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

    // compute the space-fixed coordinates for the vertices of the face

    xi1[0] = x[ibody][0] + discrete[ifirst+npi1][0];
    xi1[1] = x[ibody][1] + discrete[ifirst+npi1][1];
    xi1[2] = x[ibody][2] + discrete[ifirst+npi1][2];

    xi2[0] = x[ibody][0] + discrete[ifirst+npi2][0];
    xi2[1] = x[ibody][1] + discrete[ifirst+npi2][1];
    xi2[2] = x[ibody][2] + discrete[ifirst+npi2][2];

    xi3[0] = x[ibody][0] + discrete[ifirst+npi3][0];
    xi3[1] = x[ibody][1] + discrete[ifirst+npi3][1];
    xi3[2] = x[ibody][2] + discrete[ifirst+npi3][2];

    // find the normal unit vector of the face

    MathExtra::sub3(xi2, xi1, ui);
    MathExtra::sub3(xi3, xi1, vi);
    MathExtra::cross3(ui, vi, n);
    MathExtra::norm3(n);

    // skip if the COM of the two bodies are in the same side of the face

    if (opposite_sides(n, xi1, x[ibody], x[jbody]) == 0) continue;

    // find the projection of the sphere on the face

    project_pt_plane(x[jbody], xi1, xi2, xi3, h, d, inside);

    inside_polygon(ibody, ni, x[ibody], h, NULL, inside, tmp);
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

    if (R <= 0) { // in contact

      // compute the velocity of the vertex in the space-fixed frame

      ibonus = atom->body[ibody];
      bonus = &avec->bonus[ibonus];
      quat = bonus->quat;
      inertia = bonus->inertia;
      total_velocity(h, x[ibody], v[ibody], angmom[ibody],
                     inertia, quat, vti);

      // relative translational velocity

      vr1 = vti[0] - v[jbody][0];
      vr2 = vti[1] - v[jbody][1];
      vr3 = vti[2] - v[jbody][2];

      // normal component

      rsqinv = 1.0/rsq;
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

    if (evflag) ev_tally_xyz(ibody,jbody,nlocal,newton_pair,
                           energy,0.0,fx,fy,fz,delx,dely,delz);
  }
}

/* ----------------------------------------------------------------------
   Determine the interaction mode between i's edges against j's edges

   i = atom i (body i)
   j = atom j (body j)
   x      = atoms' coordinates
   f      = atoms' forces
   torque = atoms' torques
   tag    = atoms' tags
   contact_list = list of contacts
   num_contacts = number of contacts between i's edges and j's edges
   Return:

---------------------------------------------------------------------- */

int PairBodyRoundedPolyhedron::edge_against_edge(int ibody, int jbody,
  int itype, int jtype, double** x, Contact* contact_list, int &num_contacts,
  double &evdwl, double* facc)
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

    for (nj = 0; nj < nej; nj++) {

      // compute the distance between the edge nj to the edge ni
      #ifdef _POLYHEDRON_DEBUG
      printf("Compute interaction between edge %d of body %d "
             "with edge %d of body %d:\n",
             nj, jbody, ni, ibody);
      #endif

      interact = interaction_edge_to_edge(ibody, ni, x[ibody], rradi,
                                          jbody, nj, x[jbody], rradj,
                                          itype, jtype, cut_inner,
                                          contact_list, num_contacts,
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
   contact_list = list of contacts
   num_contacts = number of contacts between i's edges and j's faces
   Return:

---------------------------------------------------------------------- */

int PairBodyRoundedPolyhedron::edge_against_face(int ibody, int jbody,
  int itype, int jtype, double** x, Contact* contact_list, int &num_contacts,
  double &evdwl, double* facc)
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

    // loop through body j's faces

    for (nj = 0; nj < nfj; nj++) {

      // compute the distance between the face nj to the edge ni
      #ifdef _POLYHEDRON_DEBUG
      printf("Compute interaction between face %d of body %d with "
             "edge %d of body %d:\n",
             nj, jbody, ni, ibody);
      #endif

      interact = interaction_face_to_edge(jbody, nj, x[jbody], rradj,
                                          ibody, ni, x[ibody], rradi,
                                          itype, jtype, cut_inner,
                                          contact_list, num_contacts,
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
  int edge_index_i,  double *xmi, double rounded_radius_i,
  int jbody, int edge_index_j, double *xmj, double rounded_radius_j,
  int itype, int jtype, double cut_inner,
  Contact* contact_list, int &num_contacts, double &energy, double* facc)
{
  int ifirst,iefirst,jfirst,jefirst,npi1,npi2,npj1,npj2,interact;
  double xi1[3],xi2[3],xpj1[3],xpj2[3];
  double r,t1,t2,h1[3],h2[3];
  double contact_dist;

  double** x = atom->x;
  double** v = atom->v;
  double** f = atom->f;
  double** torque = atom->torque;
  double** angmom = atom->angmom;

  ifirst = dfirst[ibody];
  iefirst = edfirst[ibody];
  npi1 = static_cast<int>(edge[iefirst+edge_index_i][0]);
  npi2 = static_cast<int>(edge[iefirst+edge_index_i][1]);

  // compute the space-fixed coordinates for the edge ends

  xi1[0] = xmi[0] + discrete[ifirst+npi1][0];
  xi1[1] = xmi[1] + discrete[ifirst+npi1][1];
  xi1[2] = xmi[2] + discrete[ifirst+npi1][2];

  xi2[0] = xmi[0] + discrete[ifirst+npi2][0];
  xi2[1] = xmi[1] + discrete[ifirst+npi2][1];
  xi2[2] = xmi[2] + discrete[ifirst+npi2][2];

  // two ends of the edge from body j

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

  #ifdef _POLYHEDRON_DEBUG
  double ui[3],uj[3];
  MathExtra::sub3(xi1,xi2,ui);
  MathExtra::norm3(ui);
  MathExtra::sub3(xpj1,xpj2,uj);
  MathExtra::norm3(uj);
  double dot = MathExtra::dot3(ui, uj);
  printf("  edge npi1 = %d (%f %f %f); npi2 = %d (%f %f %f) vs."
         "  edge npj1 = %d (%f %f %f); npj2 = %d (%f %f %f): "
         "t1 = %f; t2 = %f; r = %f; dot = %f\n",
    npi1, xi1[0], xi1[1], xi1[2], npi2, xi2[0], xi2[1], xi2[2],
    npj1, xpj1[0], xpj1[1], xpj1[2], npj2, xpj2[0], xpj2[1], xpj2[2],
    t1, t2, r, dot);
  #endif

  interact = EE_NONE;

  // singularity case, ignore interactions

  if (r < EPSILON) return interact;

  // include the vertices for interactions

  if (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1 &&
      r < contact_dist + cut_inner) {
    pair_force_and_torque(jbody, ibody, h1, h2, r, contact_dist,
                          jtype, itype, x, v, f, torque, angmom,
                          jflag, energy, facc);

    interact = EE_INTERACT;
    if (r <= contact_dist) {
      // store the contact info
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
  } else {

  }

  return interact;
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
  int itype, int jtype, double cut_inner,
  Contact* contact_list, int &num_contacts, double &energy, double* facc)
{
  if (face_index >= facnum[ibody]) return EF_INVALID;

  int ifirst,iffirst,jfirst,npi1,npi2,npi3;
  int jefirst,npj1,npj2;
  double xi1[3],xi2[3],xi3[3],xpj1[3],xpj2[3],ui[3],vi[3],n[3];

  double** x = atom->x;
  double** v = atom->v;
  double** f = atom->f;
  double** torque = atom->torque;
  double** angmom = atom->angmom;

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

  // no interaction if two ends of the edge
  // are on the same side with the COM wrt the face

  if (opposite_sides(n, xi1, xmi, xpj1) == 0 &&
      opposite_sides(n, xi1, xmi, xpj2) == 0)
    return EF_NONE;

  // determine the intersection of the edge to the face

  double hi1[3], hi2[3], d1, d2, contact_dist;
  int inside1 = 0;
  int inside2 = 0;

  // enum {EF_PARALLEL=0,EF_SAME_SIDE_OF_FACE,
  //       EF_INTERSECT_INSIDE,EF_INTERSECT_OUTSIDE};

  int interact = edge_face_intersect(xi1, xi2, xi3, xpj1, xpj2,
                                     hi1, hi2, d1, d2, inside1, inside2);

  inside_polygon(ibody, face_index, xmi, hi1, hi2, inside1, inside2);

  contact_dist = rounded_radius_i + rounded_radius_j;

  // both endpoints are on the same side of, or parallel to, the face
  // and both are out of the interaction zone

  if (interact == EF_SAME_SIDE_OF_FACE || interact == EF_PARALLEL) {

    if (d1 > contact_dist + cut_inner && d2 > contact_dist + cut_inner)
      return EF_NONE;

    int num_outside = 0;
    int jflag = 1;

    #ifdef _POLYHEDRON_DEBUG
    if (interact == EF_SAME_SIDE_OF_FACE)
      printf(" - same side of face\n");
    else if (interact == EF_PARALLEL)
      printf(" - parallel\n");
    printf("     face: xi1 (%f %f %f) xi2 (%f %f %f) xi3 (%f %f %f)\n",
      xi1[0], xi1[1], xi1[2], xi2[0], xi2[1], xi2[2], xi3[0], xi3[1], xi3[2]);
    printf("     edge: xpj1 (%f %f %f) xpj2 (%f %f %f)\n",
      xpj1[0], xpj1[1], xpj1[2], xpj2[0], xpj2[1], xpj2[2]);
    #endif

    // xpj1 is in the interaction zone
    // and its projection on the face is inside the triangle
    // compute vertex-face interaction and accumulate force/torque to both bodies

    if (d1 <= contact_dist + cut_inner) {
      if (inside1) {
        if (static_cast<int>(discrete[jfirst+npj1][6]) == 0) {
          pair_force_and_torque(jbody, ibody, xpj1, hi1, d1, contact_dist,
                                jtype, itype, x, v, f, torque, angmom,
                                jflag, energy, facc);
          #ifdef _POLYHEDRON_DEBUG
          printf(" - compute pair force between vertex %d from edge %d of body %d "
                 "with face %d of body %d: d1 = %f\n",
            npj1, edge_index, jbody, face_index, ibody, d1);
          #endif

          if (d1 <= contact_dist) {
            // store the contact info
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

          discrete[jfirst+npj1][6] = 1;
        }
      } else {
        num_outside++;
      }
    }

    // xpj2 is in the interaction zone
    // and its projection on the face is inside the triangle
    // compute vertex-face interaction and accumulate force/torque to both bodies

    if (d2 <= contact_dist + cut_inner) {
      if (inside2) {
        if (static_cast<int>(discrete[jfirst+npj2][6]) == 0) {
          pair_force_and_torque(jbody, ibody, xpj2, hi2, d2, contact_dist,
                                jtype, itype, x, v, f, torque, angmom,
                                jflag, energy, facc);
          #ifdef _POLYHEDRON_DEBUG
          printf(" - compute pair force between vertex %d from edge %d of body %d "
                 "with face %d of body %d: d2 = %f\n",
                 npj2, edge_index, jbody, face_index, ibody, d2);
          #endif

          if (d2 <= contact_dist) {
            // store the contact info
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
          discrete[jfirst+npj2][6] = 1;
        }
      } else {
        num_outside++;
      }
    }

    // both ends have projection outside of the face
    // compute interaction between the edge with the three edges of the face

    if (num_outside == 2) {

      #ifdef _POLYHEDRON_DEBUG
      printf(" - outside = 2\n");
      printf(" - compute pair force between edge %d of body %d "
             "with 3 edges of face %d of body %d\n",
        edge_index, jbody, face_index, ibody);
      #endif

      interact = EF_INTERSECT_OUTSIDE;

    }

  } else if (interact == EF_INTERSECT_OUTSIDE) {

    // compute interaction between the edge with the three edges of the face

    #ifdef _POLYHEDRON_DEBUG
    printf(" - intersect outside triangle\n");
    printf(" - compute pair force between edge %d of body %d "
           "with face %d of body %d\n", edge_index, jbody, face_index, ibody);
    printf("     face: xi1 (%f %f %f) xi2 (%f %f %f) xi3 (%f %f %f)\n",
      xi1[0], xi1[1], xi1[2], xi2[0], xi2[1], xi2[2], xi3[0], xi3[1], xi3[2]);
    printf("     edge: xpj1 (%f %f %f) xpj2 (%f %f %f)\n",
      xpj1[0], xpj1[1], xpj1[2], xpj2[0], xpj2[1], xpj2[2]);

    #endif
  } else if (interact == EF_INTERSECT_INSIDE) {
    // need to do something here to resolve overlap!!
    // p is the intersection between the edge and the face
    int jflag = 1;
    if (d1 < d2)
      pair_force_and_torque(jbody, ibody, xpj1, hi1, d1, contact_dist,
                            jtype, itype, x, v, f, torque, angmom,
                            jflag, energy, facc);
    else
      pair_force_and_torque(jbody, ibody, xpj2, hi2, d2, contact_dist,
                            jtype, itype, x, v, f, torque, angmom,
                            jflag, energy, facc);
  }

  return interact;
}

/* ----------------------------------------------------------------------
  Compute forces and torques between two bodies caused by the interaction
  between a pair of points on either bodies (similar to sphere-sphere)
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::pair_force_and_torque(int ibody, int jbody,
                 double* pi, double* pj, double r, double contact_dist,
                 int itype, int jtype, double** x,
                 double** v, double** f, double** torque, double** angmom,
                 int jflag, double& energy, double* facc)
{
  double delx,dely,delz,R,fx,fy,fz,fpair;

  delx = pi[0] - pj[0];
  dely = pi[1] - pj[1];
  delz = pi[2] - pj[2];
  R = r - contact_dist;

  kernel_force(R, itype, jtype, energy, fpair);

  fx = delx*fpair/r;
  fy = dely*fpair/r;
  fz = delz*fpair/r;

  #ifdef _POLYHEDRON_DEBUG
  printf("  - R = %f; r = %f; k_na = %f; shift = %f; fpair = %f;"
         " energy = %f; jflag = %d\n", R, r, k_na, shift, fpair,
         energy, jflag);
  #endif

  if (R <= 0) {

    // contact: accumulate normal and tangential contact force components

    contact_forces(ibody, jbody, pi, pj, delx, dely, delz, fx, fy, fz,
                   x, v, angmom, f, torque, facc);
  } else {

    // accumulate force and torque to both bodies directly

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
  double kn = k_n[itype][jtype];
  double kna = k_na[itype][jtype];
  double shift = kna * cut_inner;
  double e = 0;
  if (R <= 0) {           // deformation occurs
    fpair = -kn * R - shift;
    e = (0.5 * kn * R + shift) * R;
  } else if (R <= cut_inner) {   // not deforming but cohesive ranges overlap
    fpair = kna * R - shift;
    e = (-0.5 * kna * R + shift) * R;
  } else fpair = 0.0;
  energy += e;
}

/* ----------------------------------------------------------------------
  Compute contact forces between two bodies
  modify the force stored at the vertex and edge in contact by j_a
  sum forces and torque to the corresponding bodies
  fx,fy,fz = unscaled cohesive forces
  fn = normal friction component
  ft = tangential friction component (-c_t * v_t)
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::contact_forces(int ibody, int jbody,
  double *xi, double *xj, double delx, double dely, double delz,
  double fx, double fy, double fz, double** x, double** v, double** angmom,
  double** f, double** torque, double* facc)
{
  int ibonus,jbonus;
  double rsq,rsqinv,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double fn[3],ft[3],vi[3],vj[3];
  double *quat, *inertia;
  AtomVecBody::Bonus *bonus;

  // compute the velocity of the vertex in the space-fixed frame

  ibonus = atom->body[ibody];
  bonus = &avec->bonus[ibonus];
  quat = bonus->quat;
  inertia = bonus->inertia;
  total_velocity(xi, x[ibody], v[ibody], angmom[ibody],
                 inertia, quat, vi);

  // compute the velocity of the point on the edge in the space-fixed frame

  jbonus = atom->body[jbody];
  bonus = &avec->bonus[jbonus];
  quat = bonus->quat;
  inertia = bonus->inertia;
  total_velocity(xj, x[jbody], v[jbody], angmom[jbody],
                 inertia, quat, vj);

  // vector pointing from the contact point on ibody to that on jbody

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

  // these are contact forces (F_n, F_t and F_ne) only
  // cohesive forces will be scaled by j_a after contact area is computed
  // mu * fne = tangential friction deformation during gross sliding
  // see Eq. 4, Fraige et al.

  fx = fn[0] + ft[0] + mu * fx;
  fy = fn[1] + ft[1] + mu * fy;
  fz = fn[2] + ft[2] + mu * fz;

  f[ibody][0] += fx;
  f[ibody][1] += fy;
  f[ibody][2] += fz;
  sum_torque(x[ibody], xi, fx, fy, fz, torque[ibody]);

  f[jbody][0] -= fx;
  f[jbody][1] -= fy;
  f[jbody][2] -= fz;
  sum_torque(x[jbody], xj, -fx, -fy, -fz, torque[jbody]);

  facc[0] += fx; facc[1] += fy; facc[2] += fz;

  #ifdef _POLYHEDRON_DEBUG
  printf("contact ibody = %d: f = %f %f %f; torque = %f %f %f\n", ibody,
     f[ibody][0], f[ibody][1], f[ibody][2],
     torque[ibody][0], torque[ibody][1], torque[ibody][2]);
  printf("contact jbody = %d: f = %f %f %f; torque = %f %f %f\n", jbody,
     f[jbody][0], f[jbody][1], f[jbody][2],
     torque[jbody][0], torque[jbody][1], torque[jbody][2]);
  #endif
}

/* ----------------------------------------------------------------------
  Rescale the forces and torques for all the contacts
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::rescale_cohesive_forces(double** x,
     double** f, double** torque, Contact* contact_list, int &num_contacts,
     int itype, int jtype, double* facc)
{
  int m,ibody,jbody;
  double delx,dely,delz,fx,fy,fz,R,fpair,r,contact_area;

  int num_unique_contacts = 0;
  if (num_contacts == 1) {
    num_unique_contacts = 1;
    contact_area = 0;
  } else if (num_contacts == 2) {
    num_unique_contacts = 2;
    contact_area = num_contacts * A_ua;
  } else {
    find_unique_contacts(contact_list, num_contacts);

    double xc[3],dx,dy,dz;
    xc[0] = xc[1] = xc[2] = 0;
    num_unique_contacts = 0;
    for (int m = 0; m < num_contacts; m++) {
      if (contact_list[m].unique == 0) continue;
      xc[0] += contact_list[m].xi[0];
      xc[1] += contact_list[m].xi[1];
      xc[2] += contact_list[m].xi[2];
      num_unique_contacts++;
    }

    xc[0] /= (double)num_unique_contacts;
    xc[1] /= (double)num_unique_contacts;
    xc[2] /= (double)num_unique_contacts;

    contact_area = 0.0;
    for (int m = 0; m < num_contacts; m++) {
      if (contact_list[m].unique == 0) continue;
      dx = contact_list[m].xi[0] - xc[0];
      dy = contact_list[m].xi[1] - xc[1];
      dz = contact_list[m].xi[2] - xc[2];
      contact_area += (dx*dx + dy*dy + dz*dz);
    }
    contact_area *= (MY_PI/(double)num_unique_contacts);
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

    f[ibody][0] += fx;
    f[ibody][1] += fy;
    f[ibody][2] += fz;
    sum_torque(x[ibody], contact_list[m].xi, fx, fy, fz, torque[ibody]);

    f[jbody][0] -= fx;
    f[jbody][1] -= fy;
    f[jbody][2] -= fz;
    sum_torque(x[jbody], contact_list[m].xj, -fx, -fy, -fz, torque[jbody]);

    facc[0] += fx; facc[1] += fy; facc[2] += fz;
  }
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
  Test if a line segment defined by two points a and b intersects with
  a triangle defined by three points x1, x2 and x3
------------------------------------------------------------------------- */

int PairBodyRoundedPolyhedron::edge_face_intersect(double* x1, double* x2,
               double* x3, double* a, double* b, double* h_a, double* h_b,
               double& d_a, double& d_b, int& inside_a, int& inside_b)
{
  double s[3], u[3], v[3], n[3];

  // line director

  MathExtra::sub3(b, a, s);

  // plane normal vector

  MathExtra::sub3(x2, x1, u);
  MathExtra::sub3(x3, x1, v);
  MathExtra::cross3(u, v, n);
  MathExtra::norm3(n);

  // find the projection of a and b to the plane and the corresponding distances

  project_pt_plane(a, x1, x2, x3, h_a, d_a, inside_a);

  project_pt_plane(b, x1, x2, x3, h_b, d_b, inside_b);

  // check if the line segment is parallel to the plane

  double dot = MathExtra::dot3(s, n);
  if (fabs(dot) < EPSILON) return EF_PARALLEL;

  // solve for the intersection between the line and the plane

  double m[3][3], invm[3][3], p[3], ans[3];
  m[0][0] = -s[0];
  m[0][1] = u[0];
  m[0][2] = v[0];

  m[1][0] = -s[1];
  m[1][1] = u[1];
  m[1][2] = v[1];

  m[2][0] = -s[2];
  m[2][1] = u[2];
  m[2][2] = v[2];

  MathExtra::sub3(a, x1, p);
  MathExtra::invert3(m, invm);
  MathExtra::matvec(invm, p, ans);

  // p is reused for the intersection point
  // s = b - a

  double t = ans[0];
  p[0] = a[0] + s[0] * t;
  p[1] = a[1] + s[1] * t;
  p[2] = a[2] + s[2] * t;

  // check if p is inside the triangle, excluding the edges and vertices
  // the edge-edge and edge-vertices are handled separately

  int inside = 0;
  if (ans[1] > 0 && ans[2] > 0 && ans[1] + ans[2] < 1)
    inside = 1;

  int interact;
  if (t < 0 || t > 1) {
    interact = EF_SAME_SIDE_OF_FACE;
  } else {
    if (inside == 1)
      interact = EF_INTERSECT_INSIDE;
    else
      interact = EF_INTERSECT_OUTSIDE;
  }

  return interact;
}

/* ----------------------------------------------------------------------
  Find the projection of q on the plane defined by point p and the normal
  unit vector n: q_proj = q - dot(q - p, n) * n
  and the distance d from q to the plane
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::project_pt_plane(const double* q,
                                        const double* p, const double* n,
                                        double* q_proj, double &d)
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
    q2          = another point (can be NULL) for face-edge intersection
  Output:
    inside1     = 1 if q1 is inside the polygon, 0 otherwise
    inside2     = 1 if q2 is inside the polygon, 0 otherwise
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::inside_polygon(int ibody, int face_index,
                            double* xmi, const double* q1, const double* q2,
                            int& inside1, int& inside2)

{
  int i,n,ifirst,iffirst,npi1,npi2;
  double xi1[3],xi2[3],u[3],v[3],costheta,anglesum1,anglesum2,magu,magv;

  ifirst = dfirst[ibody];
  iffirst = facfirst[ibody];
  anglesum1 = anglesum2 = 0;;
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

    if (magu * magv < EPSILON) inside1 = 1;
    else {
      costheta = MathExtra::dot3(u,v)/(magu*magv);
      anglesum1 += acos(costheta);
    }

    if (q2 != NULL) {
      MathExtra::sub3(xi1,q2,u);
      MathExtra::sub3(xi2,q2,v);
      magu = MathExtra::len3(u);
      magv = MathExtra::len3(v);
      if (magu * magv < EPSILON) inside2 = 1;
      else {
        costheta = MathExtra::dot3(u,v)/(magu*magv);
        anglesum2 += acos(costheta);
      }
    }
  }

  if (fabs(anglesum1 - MY_2PI) < EPSILON) inside1 = 1;
  else inside1 = 0;

  if (q2 != NULL) {
    if (fabs(anglesum2 - MY_2PI) < EPSILON) inside2 = 1;
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
      double &d, int& inside)
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
     const double* xi1, const double* xi2, double* h, double& d, double& t)
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

  if (fabs(xi2[0] - xi1[0]) > 0)
    t = (h[0] - xi1[0])/(xi2[0] - xi1[0]);
  else if (fabs(xi2[1] - xi1[1]) > 0)
    t = (h[1] - xi1[1])/(xi2[1] - xi1[1]);
  else if (fabs(xi2[2] - xi1[2]) > 0)
    t = (h[2] - xi1[2])/(xi2[2] - xi1[2]);
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
                  double* h1, double* h2, double& t1, double& t2, double& r)
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

  // check if two edges are parallel
  // find the two ends of the overlapping segment, if any

  if (fabs(dot - 1.0) < EPSILON) {

    double s1,s2,x13[3],x23[3],x13h[3];
    double t13,t23,t31,t41,x31[3],x41[3];
    t13=t23=t31=t41=0.0;

    MathExtra::sub3(x1,x3,x13); // x13 = x1 - x3
    MathExtra::sub3(x2,x3,x23); // x23 = x2 - x3

    s1 = MathExtra::dot3(x13,v);
    x13h[0] = x13[0] - s1*v[0];
    x13h[1] = x13[1] - s1*v[1];
    x13h[2] = x13[2] - s1*v[2];
    r = MathExtra::len3(x13h);

    // x13 is the projection of x1 on x3-x4

    x13[0] = x3[0] + s1*v[0];
    x13[1] = x3[1] + s1*v[1];
    x13[2] = x3[2] + s1*v[2];

    // x23 is the projection of x2 on x3-x4

    s2 = MathExtra::dot3(x23,v);
    x23[0] = x3[0] + s2*v[0];
    x23[1] = x3[1] + s2*v[1];
    x23[2] = x3[2] + s2*v[2];

    // find the fraction of the projection points on the edges

    if (fabs(x4[0] - x3[0]) > 0)
      t13 = (x13[0] - x3[0])/(x4[0] - x3[0]);
    else if (fabs(x4[1] - x3[1]) > 0)
      t13 = (x13[1] - x3[1])/(x4[1] - x3[1]);
    else if (fabs(x4[2] - x3[2]) > 0)
      t13 = (x13[2] - x3[2])/(x4[2] - x3[2]);

    if (fabs(x4[0] - x3[0]) > 0)
      t23 = (x23[0] - x3[0])/(x4[0] - x3[0]);
    else if (fabs(x4[1] - x3[1]) > 0)
      t23 = (x23[1] - x3[1])/(x4[1] - x3[1]);
    else if (fabs(x4[2] - x3[2]) > 0)
      t23 = (x23[2] - x3[2])/(x4[2] - x3[2]);

    if (fabs(x23[0] - x13[0]) > 0)
      t31 = (x3[0] - x13[0])/(x23[0] - x13[0]);
    else if (fabs(x23[1] - x13[1]) > 0)
      t31 = (x3[1] - x13[1])/(x23[1] - x13[1]);
    else if (fabs(x23[2] - x13[2]) > 0)
      t31 = (x3[2] - x13[2])/(x23[2] - x13[2]);

    // x31 is the projection of x3 on x1-x2

    x31[0] = x1[0] + t31*(x2[0] - x1[0]);
    x31[1] = x1[1] + t31*(x2[1] - x1[1]);
    x31[2] = x1[2] + t31*(x2[2] - x1[2]);

    if (fabs(x23[0] - x13[0]) > 0)
      t41 = (x4[0] - x13[0])/(x23[0] - x13[0]);
    else if (fabs(x23[1] - x13[1]) > 0)
      t41 = (x4[1] - x13[1])/(x23[1] - x13[1]);
    else if (fabs(x23[2] - x13[2]) > 0)
      t41 = (x4[2] - x13[2])/(x23[2] - x13[2]);

    // x41 is the projection of x4 on x1-x2

    x41[0] = x1[0] + t41*(x2[0] - x1[0]);
    x41[1] = x1[1] + t41*(x2[1] - x1[1]);
    x41[2] = x1[2] + t41*(x2[2] - x1[2]);

    // determine two ends from the overlapping segments

    int n1 = 0;
    int n2 = 0;
    if (t13 >= 0 && t13 <= 1) {
      h1[0] = x1[0];
      h1[1] = x1[1];
      h1[2] = x1[2];
      h2[0] = x13[0];
      h2[1] = x13[1];
      h2[2] = x13[2];
      t1 = 0;
      t2 = t13;
      n1++;
      n2++;
    }
    if (t23 >= 0 && t23 <= 1) {
      if (n1 == 0) {
        h1[0] = x2[0];
        h1[1] = x2[1];
        h1[2] = x2[2];
        h2[0] = x23[0];
        h2[1] = x23[1];
        h2[2] = x23[2];
        t1 = 1;
        t2 = t23;
        n1++;
        n2++;
      } else {
        h1[0] = (x1[0]+x2[0])/2;
        h1[1] = (x1[1]+x2[1])/2;
        h1[2] = (x1[2]+x2[2])/2;
        h2[0] = (x13[0]+x23[0])/2;
        h2[1] = (x13[1]+x23[1])/2;
        h2[2] = (x13[2]+x23[2])/2;
        t1 = 0.5;
        t2 = (t13+t23)/2;
        n1++;
        n2++;
      }
    }

    if (n1 == 0 && n2 == 0) {
      if (t31 >= 0 && t31 <= 1) {
        h1[0] = x31[0];
        h1[1] = x31[1];
        h1[2] = x31[2];
        h2[0] = x3[0];
        h2[1] = x3[1];
        h2[2] = x3[2];
        t1 = t31;
        t2 = 0;
        n1++;
        n2++;
      }
      if (t41 >= 0 && t41 <= 1) {
        if (n1 == 0) {
          h1[0] = x41[0];
          h1[1] = x41[1];
          h1[2] = x41[2];
          h2[0] = x4[0];
          h2[1] = x4[1];
          h2[2] = x4[2];
          t1 = t41;
          t2 = 1;
          n1++;
          n2++;
        } else {
          h1[0] = (x31[0]+x41[0])/2;
          h1[1] = (x31[1]+x41[1])/2;
          h1[2] = (x31[2]+x41[2])/2;
          h2[0] = (x3[0]+x4[0])/2;
          h2[1] = (x3[1]+x4[1])/2;
          h2[2] = (x3[2]+x4[2])/2;
          t1 = (t31+t41)/2;
          t2 = 0.5;
          n1++;
          n2++;
        }
      }
    }

    // if n1 == 0 and n2 == 0 at this point,
    // which means no overlapping segments bt two parallel edges,
    // return the default values of t1 and t2

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

  if (fabs(x2[0] - x1[0]) > 0)
    t1 = (h1[0] - x1[0])/(x2[0] - x1[0]);
  else if (fabs(x2[1] - x1[1]) > 0)
    t1 = (h1[1] - x1[1])/(x2[1] - x1[1]);
  else if (fabs(x2[2] - x1[2]) > 0)
    t1 = (h1[2] - x1[2])/(x2[2] - x1[2]);
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

void PairBodyRoundedPolyhedron::find_unique_contacts(Contact* contact_list,
                                                     int& num_contacts)
{
  int n = num_contacts;
  for (int i = 0; i < n - 1; i++) {

    for (int j = i + 1; j < n; j++) {
      if (contact_list[i].unique == 0) continue;
      double d = contact_separation(contact_list[i], contact_list[j]);
      if (d < EPSILON) contact_list[j].unique = 0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolyhedron::sanity_check()
{

  double x1[3],x2[3],x3[3],x4[3],h_a[3],h_b[3],d_a,d_b;
  double a[3],b[3],t_a,t_b;

  x1[0] = 0; x1[1] = 3; x1[2] = 0;
  x2[0] = 3; x2[1] = 0; x2[2] = 0;
  x3[0] = 4; x3[1] = 3; x3[2] = 0;
  x4[0] = 5; x4[1] = 3; x4[2] = 0;

  a[0] = 0; a[1] = 0; a[2] = 0;
  b[0] = 4; b[1] = 0; b[2] = 0;

  project_pt_line(a, x1, x2, h_a, d_a, t_a);
  project_pt_line(b, x1, x2, h_b, d_b, t_b);
/*
  printf("h_a: %f %f %f; h_b: %f %f %f; t_a = %f; t_b = %f; d = %f; d_b = %f\n",
    h_a[0], h_a[1], h_a[2], h_b[0], h_b[1], h_b[2], t_a, t_b, d_a, d_b);
*/
/*
  int inside_a, inside_b;
  int mode = edge_face_intersect(x1, x2, x3, a, b, h_a, h_b, d_a, d_b,
                                 inside_a, inside_b);

  double u[3],v[3],n[3];
  MathExtra::sub3(x2, x1, u);
  MathExtra::sub3(x3, x1, v);
  MathExtra::cross3(u, v, n);
  MathExtra::norm3(n);
*/
/*
  project_pt_plane(a, x1, x2, x3, h_a, d_a, inside_a);
  printf("h_a: %f %f %f; d = %f: inside %d\n",
    h_a[0], h_a[1], h_a[2], d_a, inside_a);
  project_pt_plane(b, x1, x2, x3, h_b, d_b, inside_b);
  printf("h_b: %f %f %f; d = %f: inside %d\n",
    h_b[0], h_b[1], h_b[2], d_b, inside_b);
*/
/*
  distance_bt_edges(x1, x2, x3, x4, h_a, h_b, t_a, t_b, d_a);
  printf("h_a: %f %f %f; h_b: %f %f %f; t_a = %f; t_b = %f; d = %f\n",
    h_a[0], h_a[1], h_a[2], h_b[0], h_b[1], h_b[2], t_a, t_b, d_a);
*/
}

