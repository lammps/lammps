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
------------------------------------------------------------------------- */

#include "fix_wall_body_polyhedron.h"
#include <cmath>
#include <cstring>
#include "atom.h"
#include "atom_vec_body.h"
#include "body_rounded_polyhedron.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{XPLANE=0,YPLANE=1,ZPLANE=2};    // XYZ PLANE need to be 0,1,2
enum{HOOKE,HOOKE_HISTORY};

enum {INVALID=0,NONE=1,VERTEX=2};
enum {FAR=0,XLO,XHI,YLO,YHI,ZLO,ZHI};

//#define _POLYHEDRON_DEBUG
#define DELTA 10000
#define EPSILON 1e-2
#define BIG 1.0e20
#define MAX_CONTACTS 4  // maximum number of contacts for 2D models
#define EFF_CONTACTS 2  // effective contacts for 2D models

/* ---------------------------------------------------------------------- */

FixWallBodyPolyhedron::FixWallBodyPolyhedron(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix wall/body/polyhedron command");

  if (!atom->body_flag)
    error->all(FLERR,"Fix wall/body/polyhedron requires "
               "atom style body/rounded/polyhedron");

  restart_peratom = 1;
  create_attribute = 1;
  wallstyle = -1;

  // wall/particle coefficients

  kn = utils::numeric(FLERR,arg[3],false,lmp);

  c_n = utils::numeric(FLERR,arg[4],false,lmp);
  if (strcmp(arg[5],"NULL") == 0) c_t = 0.5 * c_n;
  else c_t = utils::numeric(FLERR,arg[5],false,lmp);

  if (kn < 0.0 || c_n < 0.0 || c_t < 0.0)
    error->all(FLERR,"Illegal fix wall/body/polyhedron command");

  // wallstyle args

  int iarg = 6;
  if (strcmp(arg[iarg],"xplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix wall/body/polyhedron command");
    wallstyle = XPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = utils::numeric(FLERR,arg[iarg+2],false,lmp);
    iarg += 3;
  } else if (strcmp(arg[iarg],"yplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix wall/body/polyhedron command");
    wallstyle = YPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = utils::numeric(FLERR,arg[iarg+2],false,lmp);
    iarg += 3;
  } else if (strcmp(arg[iarg],"zplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix wall/body/polyhedron command");
    wallstyle = ZPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = utils::numeric(FLERR,arg[iarg+2],false,lmp);
    iarg += 3;
  } else error->all(FLERR,"Unknown wall style {}",arg[iarg]);

  // check for trailing keyword/values

  wiggle = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"wiggle") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix wall/body/polyhedron command");
      if (strcmp(arg[iarg+1],"x") == 0) axis = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) axis = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) axis = 2;
      else error->all(FLERR,"Illegal fix wall/body/polyhedron command");
      amplitude = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      period = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      wiggle = 1;
      iarg += 4;
    } else error->all(FLERR,"Illegal fix wall/body/polyhedron command");
  }

  if (wallstyle == XPLANE && domain->xperiodic)
    error->all(FLERR,"Cannot use wall in periodic dimension");
  if (wallstyle == YPLANE && domain->yperiodic)
    error->all(FLERR,"Cannot use wall in periodic dimension");
  if (wallstyle == ZPLANE && domain->zperiodic)
    error->all(FLERR,"Cannot use wall in periodic dimension");

  // setup oscillations

  if (wiggle) omega = 2.0*MY_PI / period;

  time_origin = update->ntimestep;

  dmax = nmax = 0;
  discrete = nullptr;
  dnum = dfirst = nullptr;

  edmax = ednummax = 0;
  edge = nullptr;
  ednum = edfirst = nullptr;

  facmax = facnummax = 0;
  face = nullptr;
  facnum = facfirst = nullptr;

  enclosing_radius = nullptr;
  rounded_radius = nullptr;
}

/* ---------------------------------------------------------------------- */

FixWallBodyPolyhedron::~FixWallBodyPolyhedron()
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
}

/* ---------------------------------------------------------------------- */

int FixWallBodyPolyhedron::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallBodyPolyhedron::init()
{
  dt = update->dt;

  avec = dynamic_cast<AtomVecBody *>(atom->style_match("body"));
  if (!avec)
    error->all(FLERR,"Pair body/rounded/polyhedron requires atom style body");
  if (strcmp(avec->bptr->style,"rounded/polyhedron") != 0)
    error->all(FLERR,"Pair body/rounded/polyhedron requires "
               "body style rounded/polyhedron");
  bptr = dynamic_cast<BodyRoundedPolyhedron *>(avec->bptr);

  // set pairstyle from body/polyhedronular pair style

  if (force->pair_match("body/rounded/polyhedron",1))
    pairstyle = HOOKE;
  else error->all(FLERR,"Fix wall/body/polyhedron is incompatible with Pair style");
}

/* ---------------------------------------------------------------------- */

void FixWallBodyPolyhedron::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallBodyPolyhedron::post_force(int /*vflag*/)
{
  double vwall[3],dx,dy,dz,del1,del2,rsq,wall_pos;
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
  //   if wall was set to a null pointer, it's skipped since lo/hi are infinity
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
    memory->destroy(facnum);
    memory->destroy(facfirst);
    memory->destroy(enclosing_radius);
    memory->destroy(rounded_radius);
    nmax = atom->nmax;
    memory->create(dnum,nmax,"fix:dnum");
    memory->create(dfirst,nmax,"fix:dfirst");
    memory->create(ednum,nmax,"fix:ednum");
    memory->create(edfirst,nmax,"fix:edfirst");
    memory->create(facnum,nmax,"fix:facnum");
    memory->create(facfirst,nmax,"fix:facfirst");
    memory->create(enclosing_radius,nmax,"fix:enclosing_radius");
    memory->create(rounded_radius,nmax,"fix:rounded_radius");
  }

  ndiscrete = nedge = nface = 0;
  for (i = 0; i < nlocal; i++)
    dnum[i] = ednum[i] = facnum[i] = 0;

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
      } else if (wallstyle == ZPLANE) {
        del1 = x[i][2] - wlo;
        del2 = whi - x[i][2];
        if (del1 < del2) {
          dy = del1;
          wall_pos = wlo;
          side = ZLO;
        } else {
          dy = -del2;
          wall_pos = whi;
          side = ZHI;
        }
      }

      rsq = dx*dx + dy*dy + dz*dz;
      if (rsq > radius[i]*radius[i]) continue;

      if (dnum[i] == 0) body2space(i);
      npi = dnum[i];
      ifirst = dfirst[i];
      nei = ednum[i];
      iefirst = edfirst[i];

      if (npi == 1) {
        sphere_against_wall(i, wall_pos, side, vwall, x, v, f, angmom, torque);
        continue;
      }

      // reset vertex and edge forces

      for (ni = 0; ni < npi; ni++) {
        discrete[ifirst+ni][3] = 0;
        discrete[ifirst+ni][4] = 0;
        discrete[ifirst+ni][5] = 0;
        discrete[ifirst+ni][6] = 0;
      }

      for (ni = 0; ni < nei; ni++) {
        edge[iefirst+ni][2] = 0;
        edge[iefirst+ni][3] = 0;
        edge[iefirst+ni][4] = 0;
        edge[iefirst+ni][5] = 0;
      }

      int num_contacts;
      Contact contact_list[MAX_CONTACTS];

      num_contacts = 0;
      facc[0] = facc[1] = facc[2] = 0;
      edge_against_wall(i, wall_pos, side, vwall, x, f, torque,
                        contact_list, num_contacts, facc);

    } // group bit
  }
}

/* ---------------------------------------------------------------------- */

void FixWallBodyPolyhedron::reset_dt()
{
  dt = update->dt;
}

/* ----------------------------------------------------------------------
   convert N sub-particles in body I to space frame using current quaternion
   store sub-particle space-frame displacements from COM in discrete list
------------------------------------------------------------------------- */

void FixWallBodyPolyhedron::body2space(int i)
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
    memory->grow(discrete,dmax,7,"fix:discrete");
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
  // the first 2 columns are for vertex indices within body,
  // the last 3 for forces

  if (nedge + body_num_edges > edmax) {
    edmax += DELTA;
    memory->grow(edge,edmax,6,"fix:edge");
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
    memory->grow(face,facmax,6,"pair:face");
  }

  for (int m = 0; m < body_num_faces; m++) {
    face[nface][0] = static_cast<int>(face_pts[3*m+0]);
    face[nface][1] = static_cast<int>(face_pts[3*m+1]);
    face[nface][2] = static_cast<int>(face_pts[3*m+2]);
    face[nface][3] = 0;
    face[nface][4] = 0;
    face[nface][5] = 0;
    nface++;
  }

  enclosing_radius[i] = eradius;
  rounded_radius[i] = rradius;
}

/* ----------------------------------------------------------------------
   Determine the interaction mode between a sphere against the wall

   i = atom i (body i)
   x      = atoms' coordinates
   f      = atoms' forces
   torque = atoms' torques
---------------------------------------------------------------------- */

int FixWallBodyPolyhedron::sphere_against_wall(int i, double wall_pos,
     int /*side*/, double* vwall, double** x, double** v, double** f,
     double** angmom, double** torque)
{
  int mode;
  double rradi,hi[3],d,delx,dely,delz,R,fpair,fx,fy,fz;

  rradi = rounded_radius[i];
  mode = NONE;

  if (wallstyle == XPLANE) {
    hi[0] = wall_pos;
    hi[1] = x[i][1];
    hi[2] = x[i][2];
  } else if (wallstyle == YPLANE) {
    hi[0] = x[i][0];
    hi[1] = wall_pos;
    hi[2] = x[i][2];
  } else if (wallstyle == ZPLANE) {
    hi[0] = x[i][0];
    hi[1] = x[i][1];
    hi[2] = wall_pos;
  }

  distance(hi, x[i], d);

  if (d <= rradi) {
    delx = x[i][0] - hi[0];
    dely = x[i][1] - hi[1];
    delz = x[i][2] - hi[2];
    R = d - rradi;

    fpair = -kn * R;

    fx = delx*fpair/d;
    fy = dely*fpair/d;
    fz = delz*fpair/d;

    contact_forces(i, 1.0, x[i], hi, delx, dely, delz,
                   fx, fy, fz, x, v, angmom, f, torque, vwall);
    mode = VERTEX;
  }

  return mode;
}

/* ----------------------------------------------------------------------
   Determine the interaction mode between i's vertices against the wall

   i = atom i (body i)
   x      = atoms' coordinates
   f      = atoms' forces
   torque = atoms' torques
   Output:
     contact_list = list of contacts between i and the wall
     num_contacts = number of contacts between i's vertices and the wall
   Return:
     number of contacts of the edge to the wall (0, 1 or 2)
---------------------------------------------------------------------- */

int FixWallBodyPolyhedron::edge_against_wall(int i, double wall_pos,
     int side, double* vwall, double** x, double** /*f*/, double** /*torque*/,
     Contact* /*contact_list*/, int &/*num_contacts*/, double* /*facc*/)
{
  int ni, nei, contact;
  double rradi;

  nei = ednum[i];
  rradi = rounded_radius[i];

  contact = 0;

  // loop through body i's edges

  for (ni = 0; ni < nei; ni++)
    compute_distance_to_wall(i, ni, x[i], rradi, wall_pos, side, vwall, contact);

  return contact;
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
         VERTEX  if the tested vertex interacts with the wall
------------------------------------------------------------------------- */

int FixWallBodyPolyhedron::compute_distance_to_wall(int ibody, int edge_index,
                        double *xmi, double rounded_radius_i, double wall_pos,
                        int /*side*/, double* vwall, int &contact)
{
  int mode,ifirst,iefirst,npi1,npi2;
  double d1,d2,xpi1[3],xpi2[3],hi[3];
  double fx,fy,fz,fpair,delx,dely,delz,R;

  double** x = atom->x;
  double** v = atom->v;
  double** f = atom->f;
  double** torque = atom->torque;
  double** angmom = atom->angmom;

  // two ends of the edge from body i

  ifirst = dfirst[ibody];
  iefirst = edfirst[ibody];
  npi1 = static_cast<int>(edge[iefirst+edge_index][0]);
  npi2 = static_cast<int>(edge[iefirst+edge_index][1]);

  xpi1[0] = xmi[0] + discrete[ifirst+npi1][0];
  xpi1[1] = xmi[1] + discrete[ifirst+npi1][1];
  xpi1[2] = xmi[2] + discrete[ifirst+npi1][2];

  xpi2[0] = xmi[0] + discrete[ifirst+npi2][0];
  xpi2[1] = xmi[1] + discrete[ifirst+npi2][1];
  xpi2[2] = xmi[2] + discrete[ifirst+npi2][2];

  // determine the intersection of the edge to the wall

  mode = NONE;
  double j_a = 1.0;

  if (wallstyle == XPLANE) {
    hi[0] = wall_pos;
    hi[1] = xpi1[1];
    hi[2] = xpi1[2];
  } else if (wallstyle == YPLANE) {
    hi[0] = xpi1[0];
    hi[1] = wall_pos;
    hi[2] = xpi1[2];
  } else if (wallstyle == ZPLANE) {
    hi[0] = xpi1[0];
    hi[1] = xpi1[1];
    hi[2] = wall_pos;
  }

  distance(hi, xpi1, d1);

  if (d1 <= rounded_radius_i && static_cast<int>(discrete[ifirst+npi1][6]) == 0) {
    delx = xpi1[0] - hi[0];
    dely = xpi1[1] - hi[1];
    delz = xpi1[2] - hi[2];
    R = d1 - rounded_radius_i;

    fpair = -kn * R;

    fx = delx*fpair/d1;
    fy = dely*fpair/d1;
    fz = delz*fpair/d1;

    contact_forces(ibody, j_a, xpi1, hi, delx, dely, delz,
                   fx, fy, fz, x, v, angmom, f, torque, vwall);
    discrete[ifirst+npi1][6] = 1;
    contact++;
    mode = VERTEX;
  }

  if (wallstyle == XPLANE) {
    hi[0] = wall_pos;
    hi[1] = xpi2[1];
    hi[2] = xpi2[2];
  } else if (wallstyle == YPLANE) {
    hi[0] = xpi2[0];
    hi[1] = wall_pos;
    hi[2] = xpi2[2];
  } else if (wallstyle == ZPLANE) {
    hi[0] = xpi2[0];
    hi[1] = xpi2[1];
    hi[2] = wall_pos;
  }

  distance(hi, xpi2, d2);

  if (d2 <= rounded_radius_i && static_cast<int>(discrete[ifirst+npi2][6]) == 0) {
    delx = xpi2[0] - hi[0];
    dely = xpi2[1] - hi[1];
    delz = xpi2[2] - hi[2];
    R = d2 - rounded_radius_i;

    fpair = -kn * R;

    fx = delx*fpair/d2;
    fy = dely*fpair/d2;
    fz = delz*fpair/d2;

    contact_forces(ibody, j_a, xpi2, hi, delx, dely, delz,
                   fx, fy, fz, x, v, angmom, f, torque, vwall);
    discrete[ifirst+npi2][6] = 1;
    contact++;
    mode = VERTEX;
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

void FixWallBodyPolyhedron::contact_forces(int ibody,
  double j_a, double *xi, double * /*xj*/, double delx, double dely, double delz,
  double fx, double fy, double fz, double** x, double** v, double** angmom,
  double** f, double** torque, double* vwall)
{
  int ibonus;
  double fxt,fyt,fzt,rsq,rsqinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double fn[3],ft[3],vi[3];
  double *quat, *inertia;
  AtomVecBody::Bonus *bonus;

  // compute the velocity of the vertex in the space-fixed frame

  ibonus = atom->body[ibody];
  bonus = &avec->bonus[ibonus];
  quat = bonus->quat;
  inertia = bonus->inertia;
  total_velocity(xi, x[ibody], v[ibody], angmom[ibody],
                 inertia, quat, vi);

  // vector pointing from the contact point on ibody to the wall

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

  fxt = fx; fyt = fy; fzt = fz;
  fx = fxt * j_a + fn[0] + ft[0];
  fy = fyt * j_a + fn[1] + ft[1];
  fz = fzt * j_a + fn[2] + ft[2];

  f[ibody][0] += fx;
  f[ibody][1] += fy;
  f[ibody][2] += fz;
  sum_torque(x[ibody], xi, fx, fy, fz, torque[ibody]);

  #ifdef _POLYHEDRON_DEBUG
  printf("From contact forces: vertex fx %f fy %f fz %f\n"
         "      torque body %d: %f %f %f\n"
         "      torque body %d: %f %f %f\n",
         fxt, fyt, fzt,
         atom->tag[ibody],torque[ibody][0],torque[ibody][1],torque[ibody][2],
         atom->tag[jbody],torque[jbody][0],torque[jbody][1],torque[jbody][2]);
  #endif
}

/* ----------------------------------------------------------------------
  Compute the contact forces between two bodies
  modify the force stored at the vertex and edge in contact by j_a
  sum forces and torque to the corresponding bodies
    fn = normal friction component
    ft = tangential friction component (-c_t * vrt)
------------------------------------------------------------------------- */

void FixWallBodyPolyhedron::contact_forces(Contact& contact, double j_a,
                      double** x, double** v, double** angmom, double** f,
                      double** torque, double* vwall, double* facc)
{
  int ibody,ibonus,ifirst,ni;
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

  #ifdef _POLYHEDRON_DEBUG
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

double FixWallBodyPolyhedron::contact_separation(const Contact& c1,
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

void FixWallBodyPolyhedron::sum_torque(double* xm, double *x, double fx,
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

void FixWallBodyPolyhedron::total_velocity(double* p, double *xcm,
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

void FixWallBodyPolyhedron::distance(const double* x2, const double* x1,
                                  double& r) {
  r = sqrt((x2[0] - x1[0]) * (x2[0] - x1[0])
    + (x2[1] - x1[1]) * (x2[1] - x1[1])
    + (x2[2] - x1[2]) * (x2[2] - x1[2]));
}
