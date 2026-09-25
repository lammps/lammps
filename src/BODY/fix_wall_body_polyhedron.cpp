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

#include "atom.h"
#include "atom_vec_body.h"
#include "body_rounded_polyhedron.h"
#include "domain.h"
#include "error.h"
#include "fix_wall.h"
#include "force.h"
#include "graphics.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{XPLANE=0,YPLANE=1,ZPLANE=2};    // XYZ PLANE need to be 0,1,2


static constexpr int DELTA = 10000;
static constexpr double BIG = 1.0e20;

/* ---------------------------------------------------------------------- */

FixWallBodyPolyhedron::FixWallBodyPolyhedron(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), avec(nullptr), bptr(nullptr), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 9) utils::missing_cmd_args(FLERR,"fix wall/body/polyhedron", error);

  if (!atom->body_flag)
    error->all(FLERR,"Fix wall/body/polyhedron requires atom style body/rounded/polyhedron");

  restart_peratom = 1;
  create_attribute = 1;
  wallstyle = -1;

  // wall/particle coefficients

  kn = utils::numeric(FLERR,arg[3],false,lmp);
  if (kn < 0.0) error->all(FLERR,3,"Illegal fix wall/body/polyhedron argument {}", kn);
  c_n = utils::numeric(FLERR,arg[4],false,lmp);
  if (c_n < 0.0)
    error->all(FLERR,4,"Illegal fix wall/body/polyhedron argument {}", c_n);
  if (strcmp(arg[5],"NULL") == 0) c_t = 0.5 * c_n;
  else c_t = utils::numeric(FLERR,arg[5],false,lmp);
  if (c_t < 0.0)
    error->all(FLERR,5,"Illegal fix wall/body/polyhedron argument {}", c_t);

  // wallstyle args

  numwalls = 0;
  if (strcmp(arg[6],"xplane") == 0) {
    wallstyle = XPLANE;
    if (strcmp(arg[7],"NULL") == 0) {
      lo = -BIG;
    } else {
      lo = utils::numeric(FLERR,arg[7],false,lmp);
      ++numwalls;
    }
    if (strcmp(arg[8],"NULL") == 0) {
      hi = BIG;
    } else {
      hi = utils::numeric(FLERR,arg[8],false,lmp);
      ++numwalls;
    }

  } else if (strcmp(arg[6],"yplane") == 0) {
    wallstyle = YPLANE;
    if (strcmp(arg[7],"NULL") == 0) {
      lo = -BIG;
    } else {
      lo = utils::numeric(FLERR,arg[7],false,lmp);
      ++numwalls;
    }
    if (strcmp(arg[8],"NULL") == 0) {
      hi = BIG;
    } else {
      hi = utils::numeric(FLERR,arg[8],false,lmp);
      ++numwalls;
    }

  } else if (strcmp(arg[6],"zplane") == 0) {
    wallstyle = ZPLANE;
    if (strcmp(arg[7],"NULL") == 0) {
      lo = -BIG;
    } else {
      lo = utils::numeric(FLERR,arg[7],false,lmp);
      ++numwalls;
    }
    if (strcmp(arg[8],"NULL") == 0) {
      hi = BIG;
    } else {
      hi = utils::numeric(FLERR,arg[8],false,lmp);
      ++numwalls;
    }
  } else error->all(FLERR, 6, "Unknown wall style {}",arg[6]);

  if ((wallstyle == XPLANE) && domain->xperiodic)
    error->all(FLERR, 6, "Cannot use wall in periodic dimension");
  if ((wallstyle == YPLANE) && domain->yperiodic)
    error->all(FLERR, 6, "Cannot use wall in periodic dimension");
  if ((wallstyle == ZPLANE) && domain->zperiodic)
    error->all(FLERR, 6, "Cannot use wall in periodic dimension");

  // check for trailing keyword/values

  wiggle = 0;
  int iarg = 9;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"wiggle") == 0) {
      if (iarg+4 > narg)
        utils::missing_cmd_args(FLERR,"fix wall/body/polyhedron wiggle", error);

      if (strcmp(arg[iarg+1],"x") == 0) axis = XPLANE;
      else if (strcmp(arg[iarg+1],"y") == 0) axis = YPLANE;
      else if (strcmp(arg[iarg+1],"z") == 0) axis = ZPLANE;
      else error->all(FLERR, iarg+1,
                      "Illegal fix wall/body/polyhedron wiggle direction {}", arg[iarg+1]);
      amplitude = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      period = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      wiggle = 1;
      iarg += 4;
    } else error->all(FLERR, iarg, "Unknown fix wall/body/polyhedron keyword {}", arg[iarg]);
  }

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

  // for rendering walls with dump image.
  if (numwalls > 0) {
    if (domain->dimension == 2) {
      // one cylinder object per wall to draw in 2d
      memory->create(imgobjs, numwalls, "fix_wall:imgobjs");
      memory->create(imgparms, numwalls, 8, "fix_wall:imgparms");
      for (int m = 0; m < numwalls; ++m) {
        imgobjs[m] = Graphics::CYLINDER;
        imgparms[m][0] = 1;    // use color of first atom type by default
      }
    } else {
      // two triangle objects per wall to draw in 3d
      memory->create(imgobjs, 2 * numwalls, "fix_wall:imgobjs");
      memory->create(imgparms, 2 * numwalls, 10, "fix_wall:imgparms");
      for (int m = 0; m < numwalls; ++m) {
        imgobjs[2 * m] = Graphics::TRIANGLE;
        imgobjs[2 * m + 1] = Graphics::TRIANGLE;
        imgparms[2 * m][0] = 1;        // use color of first atom type by default
        imgparms[2 * m + 1][0] = 1;    // use color of first atom type by default
      }
    }
  }
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

  memory->destroy(imgobjs);
  memory->destroy(imgparms);
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
    error->all(FLERR,Error::NOLASTLINE,"Pair body/rounded/polyhedron requires atom style body");
  if (strcmp(avec->bptr->style,"rounded/polyhedron") != 0)
    error->all(FLERR,Error::NOLASTLINE,
               "Pair body/rounded/polyhedron requires body style rounded/polyhedron");
  bptr = dynamic_cast<BodyRoundedPolyhedron *>(avec->bptr);

  if (!force->pair_match("^body/rounded/polyhedron",0))
    error->all(FLERR,Error::NOLASTLINE,"Fix wall/body/polyhedron is incompatible with Pair style");
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
  double vwall[3],del1,del2,wall_pos;
  int i,ni,npi,ifirst;

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

      // the nearer of the two walls, its inward unit normal, and the signed
      // distance of the center of the body from the wall, positive inside

      int dim = wallstyle;
      double nw[3] = {0.0, 0.0, 0.0};
      double scom;
      del1 = x[i][dim] - wlo;
      del2 = whi - x[i][dim];
      if (del1 < del2) {
        wall_pos = wlo;
        nw[dim] = 1.0;
        scom = del1;
      } else {
        wall_pos = whi;
        nw[dim] = -1.0;
        scom = del2;
      }
      if (scom > radius[i]) continue;

      if (dnum[i] == 0) body2space(i);
      npi = dnum[i];
      ifirst = dfirst[i];

      // every vertex of the body, or the center of a sphere,
      // interacts with the wall, using its signed distance from the wall

      for (ni = 0; ni < npi; ni++) {
        double xpi[3];
        MathExtra::add3(x[i], discrete[ifirst+ni], xpi);
        double sv = (xpi[dim] - wall_pos) * nw[dim];
        wall_force(i, xpi, nw, sv, vwall, x, v, angmom, f, torque);
      }
    } // group bit
  }

  // update wall image information
  int m = 0;
  if (wallstyle == XPLANE) {
    if (lo != -BIG) {
      FixWall::update_image_plane(m, FixWall::XLO, wlo, imgparms, domain);
      ++m;
    }
    if (hi != BIG) {
      FixWall::update_image_plane(m, FixWall::XHI, whi, imgparms, domain);
      ++m;
    }
  } else if (wallstyle == YPLANE) {
    if (lo != -BIG) {
      FixWall::update_image_plane(m, FixWall::YLO, wlo, imgparms, domain);
      ++m;
    }
    if (hi != BIG) {
      FixWall::update_image_plane(m, FixWall::YHI, whi, imgparms, domain);
      ++m;
    }
  } else if (wallstyle == ZPLANE) {
    if (lo != -BIG) {
      FixWall::update_image_plane(m, FixWall::ZLO, wlo, imgparms, domain);
      ++m;
    }
    if (hi != BIG) {
      FixWall::update_image_plane(m, FixWall::ZHI, whi, imgparms, domain);
      ++m;
    }
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
   Force between a vertex of body i, or the center of a sphere, and the wall
   xp = position of the vertex
   n  = inward unit normal of the wall
   sd = signed distance of the vertex from the wall, positive inside
   the rounded vertex overlaps the wall by -R = rounded radius - sd:
   elastic force k_n (-R) along n, and damping from the velocity of the
   body at the contact point relative to the wall, which both act at the
   contact point halfway between the rounded surface and the wall
------------------------------------------------------------------------- */

void FixWallBodyPolyhedron::wall_force(int i, const double *xp, const double *n, double sd,
                                       const double *vwall, double **x, double **v,
                                       double **angmom, double **f, double **torque)
{
  double rradi = rounded_radius[i];
  double R = sd - rradi;
  if (R >= 0.0) return;

  double pc[3], vi[3], vr[3], vn[3], vt[3], fw[3];
  for (int k = 0; k < 3; k++) pc[k] = xp[k] - 0.5 * (sd + rradi) * n[k];

  AtomVecBody::Bonus *bonus = &avec->bonus[atom->body[i]];
  total_velocity(pc, x[i], v[i], angmom[i], bonus->inertia, bonus->quat, vi);
  MathExtra::sub3(vi, vwall, vr);
  double vnnr = MathExtra::dot3(vr, n);
  for (int k = 0; k < 3; k++) {
    vn[k] = vnnr * n[k];
    vt[k] = vr[k] - vn[k];
    fw[k] = -kn * R * n[k] - c_n * vn[k] - c_t * vt[k];
  }

  f[i][0] += fw[0];
  f[i][1] += fw[1];
  f[i][2] += fw[2];
  sum_torque(x[i], pc, fw[0], fw[1], fw[2], torque[i]);
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

void FixWallBodyPolyhedron::total_velocity(const double* p, double *xcm, double* vcm, double *angmom,
                                           double *inertia, double *quat, double* vi)
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

void FixWallBodyPolyhedron::distance(const double* x2, const double* x1, double& r) {
  r = sqrt((x2[0] - x1[0]) * (x2[0] - x1[0])
    + (x2[1] - x1[1]) * (x2[1] - x1[1])
    + (x2[2] - x1[2]) * (x2[2] - x1[2]));
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image to render wall as plane
   data has been copied to dedicated storage during fix indent execution
------------------------------------------------------------------------- */

int FixWallBodyPolyhedron::image(int *&objs, double **&parms)
{
  objs = imgobjs;
  parms = imgparms;
  if (domain->dimension == 2) return numwalls;
  return 2*numwalls;
}

/* ---------------------------------------------------------------------- */

double FixWallBodyPolyhedron::memory_usage()
{
  return (double) nmax * 6 * sizeof(int);    // dnum+dfirst+ednum+edfirst+facnum+facfirst [nmax]
}
