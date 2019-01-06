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
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "fix_wall_body_polygon.h"
#include "atom.h"
#include "atom_vec_body.h"
#include "body_rounded_polygon.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "respa.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{XPLANE=0,YPLANE=1,ZCYLINDER};    // XYZ PLANE need to be 0,1,2
enum{HOOKE,HOOKE_HISTORY};

enum {INVALID=0,NONE=1,VERTEX=2};
enum {FAR=0,XLO,XHI,YLO,YHI};

//#define _POLYGON_DEBUG
#define DELTA 10000
#define EPSILON 1e-2
#define BIG 1.0e20
#define MAX_CONTACTS 4  // maximum number of contacts for 2D models
#define EFF_CONTACTS 2  // effective contacts for 2D models

/* ---------------------------------------------------------------------- */

FixWallBodyPolygon::FixWallBodyPolygon(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix wall/body/polygon command");

  if (!atom->body_flag)
    error->all(FLERR,"Fix wall/body/polygon requires "
               "atom style body/rounded/polygon");

  restart_peratom = 1;
  create_attribute = 1;

  // wall/particle coefficients

  kn = force->numeric(FLERR,arg[3]);

  c_n = force->numeric(FLERR,arg[4]);
  if (strcmp(arg[5],"NULL") == 0) c_t = 0.5 * c_n;
  else c_t = force->numeric(FLERR,arg[5]);

  if (kn < 0.0 || c_n < 0.0 || c_t < 0.0)
    error->all(FLERR,"Illegal fix wall/body/polygon command");

  // wallstyle args

  int iarg = 6;
  if (strcmp(arg[iarg],"xplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix wall/body/polygon command");
    wallstyle = XPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = force->numeric(FLERR,arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = force->numeric(FLERR,arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"yplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix wall/body/polygon command");
    wallstyle = YPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = force->numeric(FLERR,arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = force->numeric(FLERR,arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"zcylinder") == 0) {
    if (narg < iarg+2) error->all(FLERR,"Illegal fix wall/body/polygon command");
    wallstyle = ZCYLINDER;
    lo = hi = 0.0;
    cylradius = force->numeric(FLERR,arg[iarg+1]);
    iarg += 2;
  }

  // check for trailing keyword/values

  wiggle = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"wiggle") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix wall/body/polygon command");
      if (strcmp(arg[iarg+1],"x") == 0) axis = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) axis = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) axis = 2;
      else error->all(FLERR,"Illegal fix wall/body/polygon command");
      amplitude = force->numeric(FLERR,arg[iarg+2]);
      period = force->numeric(FLERR,arg[iarg+3]);
      wiggle = 1;
      iarg += 4;
    } else error->all(FLERR,"Illegal fix wall/body/polygon command");
  }

  if (wallstyle == XPLANE && domain->xperiodic)
    error->all(FLERR,"Cannot use wall in periodic dimension");
  if (wallstyle == YPLANE && domain->yperiodic)
    error->all(FLERR,"Cannot use wall in periodic dimension");
  if (wallstyle == ZCYLINDER && (domain->xperiodic || domain->yperiodic))
    error->all(FLERR,"Cannot use wall in periodic dimension");

  if (wiggle && wallstyle == ZCYLINDER && axis != 2)
    error->all(FLERR,"Invalid wiggle direction for fix wall/body/polygon");

  // setup oscillations

  if (wiggle) omega = 2.0*MY_PI / period;

  time_origin = update->ntimestep;

  dmax = nmax = 0;
  discrete = NULL;
  dnum = dfirst = NULL;

  edmax = ednummax = 0;
  edge = NULL;
  ednum = edfirst = NULL;

  enclosing_radius = NULL;
  rounded_radius = NULL;
}

/* ---------------------------------------------------------------------- */

FixWallBodyPolygon::~FixWallBodyPolygon()
{
  memory->destroy(discrete);
  memory->destroy(dnum);
  memory->destroy(dfirst);

  memory->destroy(edge);
  memory->destroy(ednum);
  memory->destroy(edfirst);

  memory->destroy(enclosing_radius);
  memory->destroy(rounded_radius);
}

/* ---------------------------------------------------------------------- */

int FixWallBodyPolygon::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallBodyPolygon::init()
{
  dt = update->dt;

  avec = (AtomVecBody *) atom->style_match("body");
  if (!avec)
    error->all(FLERR,"Pair body/rounded/polygon requires atom style body");
  if (strcmp(avec->bptr->style,"rounded/polygon") != 0)
    error->all(FLERR,"Pair body/rounded/polygon requires "
               "body style rounded/polygon");
  bptr = (BodyRoundedPolygon *) avec->bptr;

  // set pairstyle from body/polygonular pair style

  if (force->pair_match("body/rounded/polygon",1))
    pairstyle = HOOKE;
  else error->all(FLERR,"Fix wall/body/polygon is incompatible with Pair style");
}

/* ---------------------------------------------------------------------- */

void FixWallBodyPolygon::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallBodyPolygon::post_force(int /*vflag*/)
{
  double vwall[3],dx,dy,dz,del1,del2,delxy,delr,rsq,eradi,rradi,wall_pos;
  int i,ni,npi,ifirst,nei,iefirst,side;
  double facc[3];

  // set position of wall to initial settings and velocity to 0.0
  // if wiggle, set wall position and velocity accordingly

  double wlo = lo;
  double whi = hi;
  vwall[0] = vwall[1] = vwall[2] = 0.0;
  if (wiggle) {
    double arg = omega * (update->ntimestep - time_origin) * dt;
    if (wallstyle == axis) {
      wlo = lo + amplitude - amplitude*cos(arg);
      whi = hi + amplitude - amplitude*cos(arg);
    }
    vwall[axis] = amplitude*omega*sin(arg);
  }

  // loop over all my atoms
  // rsq = distance from wall
  // dx,dy,dz = signed distance from wall
  // for rotating cylinder, reset vwall based on particle position
  // skip atom if not close enough to wall
  //   if wall was set to NULL, it's skipped since lo/hi are infinity
  // compute force and torque on atom if close enough to wall
  //   via wall potential matched to pair potential

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *body = atom->body;
  double *radius = atom->radius;
  double **torque = atom->torque;
  double **angmom = atom->angmom;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // grow the per-atom lists if necessary and initialize

  if (atom->nmax > nmax) {
    memory->destroy(dnum);
    memory->destroy(dfirst);
    memory->destroy(ednum);
    memory->destroy(edfirst);
    memory->destroy(enclosing_radius);
    memory->destroy(rounded_radius);
    nmax = atom->nmax;
    memory->create(dnum,nmax,"fix:dnum");
    memory->create(dfirst,nmax,"fix:dfirst");
    memory->create(ednum,nmax,"fix:ednum");
    memory->create(edfirst,nmax,"fix:edfirst");
    memory->create(enclosing_radius,nmax,"fix:enclosing_radius");
    memory->create(rounded_radius,nmax,"fix:rounded_radius");
  }

  ndiscrete = nedge = 0;
  for (i = 0; i < nlocal; i++)
    dnum[i] = ednum[i] = 0;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      if (body[i] < 0) continue;

      dx = dy = dz = 0.0;
      side = FAR;
      if (wallstyle == XPLANE) {
        del1 = x[i][0] - wlo;
        del2 = whi - x[i][0];
        if (del1 < del2) {
          dx = del1;
          wall_pos = wlo;
          side = XLO;
        } else {
          dx = -del2;
          wall_pos = whi;
          side = XHI;
        }
      } else if (wallstyle == YPLANE) {
        del1 = x[i][1] - wlo;
        del2 = whi - x[i][1];
        if (del1 < del2) {
          dy = del1;
          wall_pos = wlo;
          side = YLO;
        } else {
          dy = -del2;
          wall_pos = whi;
          side = YHI;
        }
      } else if (wallstyle == ZCYLINDER) {
        delxy = sqrt(x[i][0]*x[i][0] + x[i][1]*x[i][1]);
        delr = cylradius - delxy;
        if (delr > eradi) dz = cylradius;
        else {
          dx = -delr/delxy * x[i][0];
          dy = -delr/delxy * x[i][1];
        }
      }

      rsq = dx*dx + dy*dy + dz*dz;
      if (rsq > radius[i]*radius[i]) continue;

      if (dnum[i] == 0) body2space(i);
      npi = dnum[i];
      ifirst = dfirst[i];
      nei = ednum[i];
      iefirst = edfirst[i];
      eradi = enclosing_radius[i];
      rradi = rounded_radius[i];

      // reset vertex and edge forces

      for (ni = 0; ni < npi; ni++) {
        discrete[ifirst+ni][3] = 0;
        discrete[ifirst+ni][4] = 0;
        discrete[ifirst+ni][5] = 0;
      }

      for (ni = 0; ni < nei; ni++) {
        edge[iefirst+ni][2] = 0;
        edge[iefirst+ni][3] = 0;
        edge[iefirst+ni][4] = 0;
      }

      int interact, num_contacts, done;
      double delta_a, delta_ua, j_a;
      Contact contact_list[MAX_CONTACTS];

      num_contacts = 0;
      facc[0] = facc[1] = facc[2] = 0;
      interact = vertex_against_wall(i, wall_pos, x, f, torque, side,
                                     contact_list, num_contacts, facc);

      if (num_contacts >= 2) {

        // find the first two distinct contacts

        done = 0;
        for (int m = 0; m < num_contacts-1; m++) {
          for (int n = m+1; n < num_contacts; n++) {
            delta_a = contact_separation(contact_list[m], contact_list[n]);
            if (delta_a > 0) {
              delta_ua = 1.0;
              j_a = delta_a / (EFF_CONTACTS * delta_ua);
              if (j_a < 1.0) j_a = 1.0;

              // scale the force at both contacts

              contact_forces(contact_list[m], j_a, x, v, angmom, f, torque,
                             vwall, facc);
              contact_forces(contact_list[n], j_a, x, v, angmom, f, torque,
                             vwall, facc);
              done = 1;
              break;
            }
          }
          if (done == 1) break;
        }

      } else if (num_contacts == 1) {

        // if there's only one contact, it should be handled here
        // since forces/torques have not been accumulated from vertex2wall()

        contact_forces(contact_list[0], 1.0, x, v, angmom, f, torque,
                       vwall, facc);
      }
    } // group bit
  }
}

/* ---------------------------------------------------------------------- */

void FixWallBodyPolygon::reset_dt()
{
  dt = update->dt;
}

/* ----------------------------------------------------------------------
   convert N sub-particles in body I to space frame using current quaternion
   store sub-particle space-frame displacements from COM in discrete list
------------------------------------------------------------------------- */

void FixWallBodyPolygon::body2space(int i)
{
  int ibonus = atom->body[i];
  AtomVecBody::Bonus *bonus = &avec->bonus[ibonus];
  int nsub = bptr->nsub(bonus);
  double *coords = bptr->coords(bonus);
  int body_num_edges = bptr->nedges(bonus);
  double* vertices = bptr->edges(bonus);
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
    memory->grow(discrete,dmax,6,"fix:discrete");
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
  // the first 2 columns are for vertex indices within body,
  // the last 3 for forces

  if (nedge + body_num_edges > edmax) {
    edmax += DELTA;
    memory->grow(edge,edmax,5,"fix:edge");
  }

  for (int m = 0; m < body_num_edges; m++) {
    edge[nedge][0] = static_cast<int>(vertices[2*m+0]);
    edge[nedge][1] = static_cast<int>(vertices[2*m+1]);
    edge[nedge][2] = 0;
    edge[nedge][3] = 0;
    edge[nedge][4] = 0;
    nedge++;
  }

  enclosing_radius[i] = eradius;
  rounded_radius[i] = rradius;
}

/* ----------------------------------------------------------------------
   Determine the interaction mode between i's vertices against the wall

   i = atom i (body i)
   x      = atoms' coordinates
   f      = atoms' forces
   torque = atoms' torques
   Return:
     contact_list = list of contacts between i and the wall
     num_contacts = number of contacts between i's vertices and the wall
     interact = 0 no interaction with the wall
                1 there's at least one vertex of i interacts
                  with the wall
---------------------------------------------------------------------- */

int FixWallBodyPolygon::vertex_against_wall(int i, double wall_pos,
                double** x, double** f, double** torque, int side,
                Contact* contact_list, int &num_contacts, double* /*facc*/)
{
  int ni, npi, ifirst, interact;
  double xpi[3], eradi, rradi;
  double fx, fy, fz;

  npi = dnum[i];
  ifirst = dfirst[i];
  eradi = enclosing_radius[i];
  rradi = rounded_radius[i];

  interact = 0;

  // loop through body i's vertices

  for (ni = 0; ni < npi; ni++) {

    // convert body-fixed coordinates to space-fixed, xi

    xpi[0] = x[i][0] + discrete[ifirst+ni][0];
    xpi[1] = x[i][1] + discrete[ifirst+ni][1];
    xpi[2] = x[i][2] + discrete[ifirst+ni][2];

    int mode, contact;
    double d, R, hi[3], delx, dely, delz, fpair;
    double rij;

    // compute the distance from the vertex xpi to the wall

    mode = compute_distance_to_wall(xpi, rradi, wall_pos, side,
                                    d, hi, contact);

    if (mode == INVALID || mode == NONE) continue;

    if (mode == VERTEX) {

      interact = 1;

      // vertex i interacts with the wall

      delx = xpi[0] - hi[0];
      dely = xpi[1] - hi[1];
      delz = xpi[2] - hi[2];

      // R = surface separation = d shifted by the rounded radius
      // R = d - p1.rounded_radius;
      // note: the force is defined for R, not for d
      // R >  0: no interaction
      // R <= 0: deformation between vertex i and the wall

      rij = sqrt(delx*delx + dely*dely + delz*delz);
      R = rij - rradi;

      // the normal frictional term -c_n * vn will be added later

      if (R <= 0) {           // deformation occurs
        fpair = -kn * R;
      } else fpair = 0.0;

      fx = delx*fpair/rij;
      fy = dely*fpair/rij;
      fz = delz*fpair/rij;

      #ifdef _POLYGON_DEBUG
      printf("  Interaction between vertex %d of %d and wall:", ni);
      printf("    mode = %d; contact = %d; d = %f; rij = %f\n",
             mode, contact, d, rij);
      printf("    R = %f\n", R);
      printf("    fpair = %f\n", fpair);
      #endif

      if (contact == 1) {

        // vertex ni of body i contacts with edge nj of body j

        contact_list[num_contacts].ibody = i;
        contact_list[num_contacts].jbody = -1;
        contact_list[num_contacts].vertex = ni;
        contact_list[num_contacts].edge = -1;
        contact_list[num_contacts].xv[0] = xpi[0];
        contact_list[num_contacts].xv[1] = xpi[1];
        contact_list[num_contacts].xv[2] = xpi[2];
        contact_list[num_contacts].xe[0] = hi[0];
        contact_list[num_contacts].xe[1] = hi[1];
        contact_list[num_contacts].xe[2] = hi[2];
        contact_list[num_contacts].separation = R;
        num_contacts++;

        // store forces to vertex ni to be rescaled later,
        // if there are 2 contacts

        discrete[ifirst+ni][3] = fx;
        discrete[ifirst+ni][4] = fy;
        discrete[ifirst+ni][5] = fz;

        #ifdef _POLYGON_DEBUG
        printf("  Stored forces at vertex and edge for accumulating later.\n");
        #endif

      } else { // no contact

        // accumulate force and torque to the body directly

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        sum_torque(x[i], xpi, fx, fy, fz, torque[i]);

      } // end if contact

    } // end if mode

  } // end for looping through the vertices of body i

  return interact;
}

/* -------------------------------------------------------------------------
  Compute the distance between a vertex to the wall
  another body
  Input:
    x0         = coordinate of the tested vertex
    rradi      = rounded radius of the vertex
    wall_pos   = position of the wall
  Output:
    d          = Distance from a point x0 to an wall
    hi         = coordinates of the projection of x0 on the wall
  contact      = 0 no contact between the queried vertex and the wall
                 1 contact detected
  return NONE    if there is no interaction
         EDGE    if the tested vertex interacts with the wall
------------------------------------------------------------------------- */

int FixWallBodyPolygon::compute_distance_to_wall(double* x0, double rradi,
          double wall_pos, int side, double &d, double hi[3], int &contact)
{
  int mode;
  double delxy;

  // h0 = position of the projection of x0 on the wall
  if (wallstyle == XPLANE) {
    hi[0] = wall_pos;
    hi[1] = x0[1];
    hi[2] = x0[2];
  } else if (wallstyle == YPLANE) {
    hi[0] = x0[0];
    hi[1] = wall_pos;
    hi[2] = x0[2];
  } else if (wallstyle == ZCYLINDER) {
    delxy = sqrt(x0[0]*x0[0] + x0[1]*x0[1]);
    hi[0] = x0[0]*cylradius/delxy;
    hi[1] = x0[1]*cylradius/delxy;
    hi[2] = x0[2];
  }

  // distance from x0 to the wall = distance from x0 to hi

  distance(hi, x0, d);

  // determine the interaction mode

  if (d < rradi) {
    mode = VERTEX;
    contact = 1;
  } else {
    mode = NONE;
    if (side == XLO) {
      if (x0[0] < wall_pos) mode = VERTEX;
    } else if (side == XHI) {
      if (x0[0] > wall_pos) mode = VERTEX;
    } else if (side == YLO) {
      if (x0[1] < wall_pos) mode = VERTEX;
    } else if (side == YHI) {
      if (x0[1] > wall_pos) mode = VERTEX;
    }
  }

  if (mode == NONE) contact = 0;
  else contact = 1;

  return mode;
}

/* ----------------------------------------------------------------------
  Compute the contact forces between two bodies
  modify the force stored at the vertex and edge in contact by j_a
  sum forces and torque to the corresponding bodies
    fn = normal friction component
    ft = tangential friction component (-c_t * vrt)
------------------------------------------------------------------------- */

void FixWallBodyPolygon::contact_forces(Contact& contact, double j_a,
                      double** x, double** v, double** angmom, double** f,
                      double** torque, double* vwall, double* facc)
{
  int ibody,ibonus,ifirst, ni;
  double fx,fy,fz,delx,dely,delz,rsq,rsqinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double fn[3],ft[3],vi[3];
  double *quat, *inertia;
  AtomVecBody::Bonus *bonus;

  ibody = contact.ibody;

  // compute the velocity of the vertex in the space-fixed frame

  ibonus = atom->body[ibody];
  bonus = &avec->bonus[ibonus];
  quat = bonus->quat;
  inertia = bonus->inertia;
  total_velocity(contact.xv, x[ibody], v[ibody], angmom[ibody],
                 inertia, quat, vi);

  // vector pointing from the vertex to the point on the wall

  delx = contact.xv[0] - contact.xe[0];
  dely = contact.xv[1] - contact.xe[1];
  delz = contact.xv[2] - contact.xe[2];
  rsq = delx*delx + dely*dely + delz*delz;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = vi[0] - vwall[0];
  vr2 = vi[1] - vwall[1];
  vr3 = vi[2] - vwall[2];

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

  ifirst = dfirst[ibody];
  ni = contact.vertex;

  fx = discrete[ifirst+ni][3] * j_a + fn[0] + ft[0];
  fy = discrete[ifirst+ni][4] * j_a + fn[1] + ft[1];
  fz = discrete[ifirst+ni][5] * j_a + fn[2] + ft[2];
  f[ibody][0] += fx;
  f[ibody][1] += fy;
  f[ibody][2] += fz;
  sum_torque(x[ibody], contact.xv, fx, fy, fz, torque[ibody]);

  // accumulate forces to the vertex only

  facc[0] += fx; facc[1] += fy; facc[2] += fz;

  #ifdef _POLYGON_DEBUG
  printf("From contact forces: vertex fx %f fy %f fz %f\n"
         "      torque body %d: %f %f %f\n",
         discrete[ifirst+ni][3], discrete[ifirst+ni][4], discrete[ifirst+ni][5],
         atom->tag[ibody],torque[ibody][0],torque[ibody][1],torque[ibody][2]);
  #endif
}

/* ----------------------------------------------------------------------
  Determine the length of the contact segment, i.e. the separation between
  2 contacts
------------------------------------------------------------------------- */

double FixWallBodyPolygon::contact_separation(const Contact& c1,
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

void FixWallBodyPolygon::sum_torque(double* xm, double *x, double fx,
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

void FixWallBodyPolygon::total_velocity(double* p, double *xcm,
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

void FixWallBodyPolygon::distance(const double* x2, const double* x1,
                                  double& r) {
  r = sqrt((x2[0] - x1[0]) * (x2[0] - x1[0])
    + (x2[1] - x1[1]) * (x2[1] - x1[1])
    + (x2[2] - x1[2]) * (x2[2] - x1[2]));
}
