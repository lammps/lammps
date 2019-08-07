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
   Ref: Fraige, Langston, Matchett and Dodds, Particuology 2008, 6:455-466
   Note: The current implementation has not taken into account
         the contact history for friction forces.
------------------------------------------------------------------------- */

#include "pair_body_rounded_polygon.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_body.h"
#include "body_rounded_polygon.h"
#include "comm.h"
#include "force.h"
#include "fix.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000
#define EPSILON 1e-3
#define MAX_CONTACTS 4  // maximum number of contacts for 2D models
#define EFF_CONTACTS 2  // effective contacts for 2D models

//#define _CONVEX_POLYGON
//#define _POLYGON_DEBUG

enum {INVALID=0,NONE=1,VERTEXI=2,VERTEXJ=3,EDGE=4};

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolygon::PairBodyRoundedPolygon(LAMMPS *lmp) : Pair(lmp)
{
  dmax = nmax = 0;
  discrete = NULL;
  dnum = dfirst = NULL;

  edmax = ednummax = 0;
  edge = NULL;
  ednum = edfirst = NULL;

  enclosing_radius = NULL;
  rounded_radius = NULL;
  maxerad = NULL;

  single_enable = 0;
  restartinfo = 0;

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

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(k_n);
    memory->destroy(k_na);
    memory->destroy(maxerad);
  }
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolygon::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  int ni,nj,npi,npj,ifirst,jfirst;
  int nei,nej,iefirst,jefirst;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl;
  double rsq,rsqinv,r,radi,radj,eradi,eradj,rradi,rradj,k_nij,k_naij;
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
  tagint* tag = atom->tag;
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

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (body[i] >= 0) {
      if (dnum[i] == 0) body2space(i);
      npi = dnum[i];
      ifirst = dfirst[i];
      nei = ednum[i];
      iefirst = edfirst[i];
      eradi = enclosing_radius[i];
      rradi = rounded_radius[i];
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
      rradj = rounded_radius[j];

      k_nij = k_n[itype][jtype];
      k_naij = k_na[itype][jtype];

      // no interaction

      r = sqrt(rsq);
      if (r > radi + radj + cut_inner) continue;
      rsqinv = 1.0 / rsq;

      if (npi == 1 && npj == 1) {
        sphere_against_sphere(i, j, delx, dely, delz, rsq,
                            k_nij, k_naij, x, v, f, evflag);
        continue;
      }

      // reset vertex and edge forces

      for (ni = 0; ni < npi; ni++) {
        discrete[ifirst+ni][3] = 0;
        discrete[ifirst+ni][4] = 0;
        discrete[ifirst+ni][5] = 0;
      }

      for (nj = 0; nj < npj; nj++) {
        discrete[jfirst+nj][3] = 0;
        discrete[jfirst+nj][4] = 0;
        discrete[jfirst+nj][5] = 0;
      }

      for (ni = 0; ni < nei; ni++) {
        edge[iefirst+ni][2] = 0;
        edge[iefirst+ni][3] = 0;
        edge[iefirst+ni][4] = 0;
      }

      for (nj = 0; nj < nej; nj++) {
        edge[jefirst+nj][2] = 0;
        edge[jefirst+nj][3] = 0;
        edge[jefirst+nj][4] = 0;
      }

      int interact, num_contacts, done;
      double delta_a, j_a;
      Contact contact_list[MAX_CONTACTS];

      num_contacts = 0;

      // check interaction between i's vertices and j' edges

      interact = vertex_against_edge(i, j, k_nij, k_naij,
                                     x, f, torque, tag, contact_list,
                                     num_contacts, evdwl, facc);

      // check interaction between j's vertices and i' edges

      interact = vertex_against_edge(j, i, k_nij, k_naij,
                                     x, f, torque, tag, contact_list,
                                     num_contacts, evdwl, facc);

      if (num_contacts >= 2) {

        // find the first two distinct contacts

        done = 0;
        for (int m = 0; m < num_contacts-1; m++) {
          for (int n = m+1; n < num_contacts; n++) {
            delta_a = contact_separation(contact_list[m], contact_list[n]);
            if (delta_a > 0) {
              j_a = delta_a / (EFF_CONTACTS * delta_ua);
              if (j_a < 1.0) j_a = 1.0;

              // scale the force at both contacts

              contact_forces(contact_list[m], j_a, x, v, angmom, f, torque,
                             evdwl, facc);
              contact_forces(contact_list[n], j_a, x, v, angmom, f, torque,
                             evdwl, facc);
              done = 1;

              #ifdef _POLYGON_DEBUG
              printf("  Two separate contacts %d and %d: delta_a = %f; j_a = %f\n",
                m, n, delta_a, j_a);
              printf("    %d: vertex %d of body %d and edge %d of body %d; "
                     "xv = %f %f %f; xe = %f %f %f\n",
                     m, contact_list[m].vertex, contact_list[m].ibody,
                     contact_list[m].edge, contact_list[m].jbody,
                     contact_list[m].xv[0], contact_list[m].xv[1],
                     contact_list[m].xv[2], contact_list[m].xe[0],
                     contact_list[m].xe[1], contact_list[m].xe[2]);
              printf("    %d: vertex %d of body %d and edge %d of body %d; "
                     "xv = %f %f %f; xe = %f %f %f\n",
                     n, contact_list[n].vertex, contact_list[n].ibody,
                     contact_list[n].edge, contact_list[n].jbody,
                     contact_list[n].xv[0], contact_list[n].xv[1],
                     contact_list[n].xv[2], contact_list[n].xe[0],
                     contact_list[n].xe[1], contact_list[n].xe[2]);
              #endif

              break;
            }
          }
          if (done == 1) break;
        }


      } else if (num_contacts == 1) {

        // if there's only one contact, it should be handled here
        // since forces/torques have not been accumulated from vertex2edge()

        contact_forces(contact_list[0], 1.0, x, v, angmom, f, torque, evdwl, facc);

        #ifdef _POLYGON_DEBUG
        printf("One contact between vertex %d of body %d and edge %d of body %d:\n",
                contact_list[0].vertex, tag[contact_list[0].ibody],
                contact_list[0].edge, tag[contact_list[0].jbody]);
        printf("xv = %f %f %f; xe = %f %f %f\n",
               contact_list[0].xv[0], contact_list[0].xv[1], contact_list[0].xv[2],
               contact_list[0].xe[0], contact_list[0].xe[1], contact_list[0].xe[2]);
        #endif
      }

      #ifdef _POLYGON_DEBUG
      int num_overlapping_contacts = 0;
      for (int m = 0; m < num_contacts-1; m++) {
        for (int n = m+1; n < num_contacts; n++) {
          double l = contact_separation(contact_list[m], contact_list[n]);
          if (l < EPSILON) num_overlapping_contacts++;
        }
      }
      printf("There are %d contacts detected, %d of which overlap.\n",
             num_contacts, num_overlapping_contacts);
      #endif

      if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,evdwl,0.0,
                               facc[0],facc[1],facc[2],delx,dely,delz);

    } // end for jj
  }

  if (vflag_fdotr) virial_fdotr_compute();
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
  memory->create(maxerad,n+1,"pair:maxerad");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBodyRoundedPolygon::settings(int narg, char **arg)
{
  if (narg < 5) error->all(FLERR,"Illegal pair_style command");

  c_n = force->numeric(FLERR,arg[0]);
  c_t = force->numeric(FLERR,arg[1]);
  mu = force->numeric(FLERR,arg[2]);
  delta_ua = force->numeric(FLERR,arg[3]);
  cut_inner = force->numeric(FLERR,arg[4]);

  if (delta_ua < 0) delta_ua = 1;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBodyRoundedPolygon::coeff(int narg, char **arg)
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

void PairBodyRoundedPolygon::init_style()
{
  avec = (AtomVecBody *) atom->style_match("body");
  if (!avec)
    error->all(FLERR,"Pair body/rounded/polygon requires atom style body");
  if (strcmp(avec->bptr->style,"rounded/polygon") != 0)
    error->all(FLERR,"Pair body/rounded/polygon requires "
               "body style rounded/polygon");
  bptr = (BodyRoundedPolygon *) avec->bptr;

  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style body/rounded/polygon requires "
               "newton pair on");

  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair body/rounded/polygon requires "
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
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBodyRoundedPolygon::init_one(int i, int j)
{
  k_n[j][i] = k_n[i][j];
  k_na[j][i] = k_na[i][j];

  return (maxerad[i]+maxerad[j]);
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
  // the first 3 columns are for coords, the last 3 for forces

  if (ndiscrete + nsub > dmax) {
    dmax += DELTA;
    memory->grow(discrete,dmax,6,"pair:discrete");
  }

  double p[3][3];
  MathExtra::quat_to_mat(bonus->quat,p);

  for (int m = 0; m < nsub; m++) {
    MathExtra::matvec(p,&coords[3*m],discrete[ndiscrete]);
    discrete[ndiscrete][3] = 0;
    discrete[ndiscrete][4] = 0;
    discrete[ndiscrete][5] = 0;
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
    memory->grow(edge,edmax,5,"pair:edge");
  }

  for (int m = 0; m < body_num_edges; m++) {
    edge[nedge][0] = static_cast<int>(edge_ends[2*m+0]);
    edge[nedge][1] = static_cast<int>(edge_ends[2*m+1]);
    edge[nedge][2] = 0;
    edge[nedge][3] = 0;
    edge[nedge][4] = 0;
    nedge++;
  }

  enclosing_radius[i] = eradius;
  rounded_radius[i] = rradius;
}

/* ----------------------------------------------------------------------
   Interaction between two spheres with different radii
   according to the 2D model from Fraige et al.
---------------------------------------------------------------------- */

void PairBodyRoundedPolygon::sphere_against_sphere(int i, int j,
                       double delx, double dely, double delz, double rsq,
                       double k_n, double k_na, double** /*x*/, double** v,
                       double** f, int evflag)
{
  double eradi,eradj,rradi,rradj;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double rij,rsqinv,R,fx,fy,fz,fn[3],ft[3],fpair,shift,energy;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  eradi = enclosing_radius[i];
  rradi = rounded_radius[i];

  eradj = enclosing_radius[j];
  rradj = rounded_radius[j];

  rsqinv = 1.0/rsq;
  rij = sqrt(rsq);
  R = rij - (rradi + rradj);
  shift = k_na * cut_inner;

  energy = 0;

  if (R <= 0) {           // deformation occurs
    fpair = -k_n * R - shift;
    energy = (0.5 * k_n * R + shift) * R;
  } else if (R <= cut_inner) {   // not deforming but cohesive ranges overlap
    fpair = k_na * R - shift;
    energy = (-0.5 * k_na * R + shift) * R;
  } else fpair = 0.0;

  fx = delx*fpair/rij;
  fy = dely*fpair/rij;
  fz = delz*fpair/rij;

  if (R <= EPSILON) { // in contact

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

    // tangential friction term at contact
    // excluding the tangential deformation term

    ft[0] = -c_t * vt1;
    ft[1] = -c_t * vt2;
    ft[2] = -c_t * vt3;
  }

  f[i][0] += fx;
  f[i][1] += fy;
  f[i][2] += fz;

  if (newton_pair || j < nlocal) {
    f[j][0] -= fx;
    f[j][1] -= fy;
    f[j][2] -= fz;
  }

  if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                           energy,0.0,fx,fy,fz,delx,dely,delz);
}

/* ----------------------------------------------------------------------
   Determine the interaction mode between i's vertices against j's edges

   i = atom i (body i)
   j = atom j (body j)
   x      = atoms' coordinates
   f      = atoms' forces
   torque = atoms' torques
   tag    = atoms' tags
   contact_list = list of contacts
   num_contacts = number of contacts between i's vertices and j's edges
   Return:
     interact = 0 no interaction at all
                1 there's at least one case where i's vertices interacts
                  with j's edges
---------------------------------------------------------------------- */

int PairBodyRoundedPolygon::vertex_against_edge(int i, int j,
                                                double k_n, double k_na,
                                                double** x, double** f,
                                                double** torque, tagint* tag,
                                                Contact* contact_list,
                                                int &num_contacts,
                                                double &evdwl, double* facc)
{
  int ni, npi, ifirst, nei, iefirst;
  int nj, npj, jfirst, nej, jefirst;
  double xpi[3], xpj[3], dist, eradi, eradj, rradi, rradj;
  double fx, fy, fz, energy;
  int interact;

  npi = dnum[i];
  ifirst = dfirst[i];
  nei = ednum[i];
  iefirst = edfirst[i];
  eradi = enclosing_radius[i];
  rradi = rounded_radius[i];

  npj = dnum[j];
  jfirst = dfirst[j];
  nej = ednum[j];
  jefirst = edfirst[j];
  eradj = enclosing_radius[j];
  rradj = rounded_radius[j];

  energy = 0;
  interact = 0;

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

    int mode, contact, p2vertex;
    double d, R, hi[3], t, delx, dely, delz, fpair, shift;
    double rij;

    // loop through body j's edges

    for (nj = 0; nj < nej; nj++) {

      // compute the distance between the edge nj to the vertex xpi

      mode = compute_distance_to_vertex(j, nj, x[j], rradj,
                                        xpi, rradi, cut_inner,
                                        d, hi, t, contact);

      if (mode == INVALID || mode == NONE) continue;

      if (mode == VERTEXI || mode == VERTEXJ) {

        interact = 1;

        // vertex i interacts with a vertex of the edge, but does not contact

        if (mode == VERTEXI) p2vertex = edge[jefirst+nj][0];
        else if (mode == VERTEXJ) p2vertex = edge[jefirst+nj][1];

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
        shift = k_na * cut_inner;

        // the normal frictional term -c_n * vn will be added later

        if (R <= 0) {           // deformation occurs
          fpair = -k_n * R - shift;
          energy += (0.5 * k_n * R + shift) * R;
        } else if (R <= cut_inner) {   // not deforming but cohesive ranges overlap
          fpair = k_na * R - shift;
          energy += (-0.5 * k_na * R + shift) * R;
        } else fpair = 0.0;

        fx = delx*fpair/rij;
        fy = dely*fpair/rij;
        fz = delz*fpair/rij;

        #ifdef _POLYGON_DEBUG
        printf("  Interaction between vertex %d of %d and vertex %d of %d:",
               ni, tag[i], p2vertex, tag[j]);
        printf("    mode = %d; contact = %d; d = %f; rij = %f, t = %f\n",
               mode, contact, d, rij, t);
        printf("    R = %f; cut_inner = %f\n", R, cut_inner);
        printf("    fpair = %f\n", fpair);
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
        shift = k_na * cut_inner;

        // the normal frictional term -c_n * vn will be added later

        if (R <= 0) {           // deformation occurs
          fpair = -k_n * R - shift;
          energy += (0.5 * k_n * R + shift) * R;
        } else if (R <= cut_inner) {   // not deforming but cohesive ranges overlap
          fpair = k_na * R - shift;
          energy += (-0.5 * k_na * R + shift) * R;
        } else fpair = 0.0;

        fx = delx*fpair/d;
        fy = dely*fpair/d;
        fz = delz*fpair/d;

        #ifdef _POLYGON_DEBUG
        printf("  Interaction between vertex %d of %d and edge %d of %d:",
               ni, tag[i], nj, tag[j]);
        printf("    mode = %d; contact = %d; d = %f; t = %f\n",
               mode, contact, d, t);
        printf("    R = %f; cut_inner = %f\n", R, cut_inner);
        printf("    fpair = %f\n", fpair);
        #endif

        if (contact == 1) {

          // vertex ni of body i contacts with edge nj of body j

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

          // store forces to vertex ni and the edge nj
          // to be rescaled later

          discrete[ifirst+ni][3] = fx;
          discrete[ifirst+ni][4] = fy;
          discrete[ifirst+ni][5] = fz;

          edge[jefirst+nj][2] = -fx;
          edge[jefirst+nj][3] = -fy;
          edge[jefirst+nj][4] = -fz;

          #ifdef _POLYGON_DEBUG
          printf("  Stored forces at vertex and edge for accumulating later.\n");
          #endif

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

      if (fabs(xi2[0] - xi1[0]) > EPSILON)
        t = (hi[0] - xi1[0]) / (xi2[0] - xi1[0]);
      else if (fabs(xi2[1] - xi1[1]) > EPSILON)
        t = (hi[1] - xi1[1]) / (xi2[1] - xi1[1]);
      else if (fabs(xi2[2] - xi1[2]) > EPSILON)
        t = (hi[2] - xi1[2]) / (xi2[2] - xi1[2]);

      double contact_dist = rounded_radius + x0_rounded_radius;
      if (t >= 0 && t <= 1) {
        mode = EDGE;
        if (d < contact_dist + EPSILON)
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
                       double** x, double** v, double** angmom, double** f,
                       double** torque, double &/*evdwl*/, double* facc)
{
  int ibody,jbody,ibonus,jbonus,ifirst,jefirst,ni,nj;
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

  // only the cohesive force is scaled by j_a
  // mu * fne = tangential friction deformation during gross sliding
  // see Eq. 4, Fraige et al.

  ifirst = dfirst[ibody];
  ni = contact.vertex;

  fx = discrete[ifirst+ni][3] * j_a + fn[0] + ft[0] +
    mu * discrete[ifirst+ni][3];
  fy = discrete[ifirst+ni][4] * j_a + fn[1] + ft[1] +
    mu * discrete[ifirst+ni][4];
  fz = discrete[ifirst+ni][5] * j_a + fn[2] + ft[2] +
    mu * discrete[ifirst+ni][5];
  f[ibody][0] += fx;
  f[ibody][1] += fy;
  f[ibody][2] += fz;
  sum_torque(x[ibody], contact.xv, fx, fy, fz, torque[ibody]);

  // accumulate forces to the vertex only

  facc[0] += fx; facc[1] += fy; facc[2] += fz;

  // only the cohesive force is scaled by j_a
  // mu * fne = tangential friction deformation during gross sliding
  // Eq. 4, Fraige et al.

  jefirst = edfirst[jbody];
  nj = contact.edge;

  fx = edge[jefirst+nj][2] * j_a - fn[0] - ft[0] +
    mu * edge[jefirst+nj][2];
  fy = edge[jefirst+nj][3] * j_a - fn[1] - ft[1] +
    mu * edge[jefirst+nj][3];
  fz = edge[jefirst+nj][4] * j_a - fn[2] - ft[2] +
    mu * edge[jefirst+nj][4];
  f[jbody][0] += fx;
  f[jbody][1] += fy;
  f[jbody][2] += fz;
  sum_torque(x[jbody], contact.xe, fx, fy, fz, torque[jbody]);

  #ifdef _POLYGON_DEBUG
  printf("From contact forces: vertex fx %f fy %f fz %f\n"
         "      torque body %d: %f %f %f\n"
         "      torque body %d: %f %f %f\n",
         discrete[ifirst+ni][3], discrete[ifirst+ni][4], discrete[ifirst+ni][5],
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

  double delta_a = 0.0;
  if (fabs(x2 - x1) > EPSILON) {
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

