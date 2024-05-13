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

#include "lattice.h"

#include "comm.h"
#include "domain.h"
#include "error.h"
#include "math_extra.h"
#include "memory.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathExtra;

static constexpr double BIG = 1.0e30;

/* ---------------------------------------------------------------------- */

Lattice::Lattice(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  nbasis = 0;
  basis = nullptr;

  // parse style arg

  if (narg < 1) utils::missing_cmd_args(FLERR, "lattice", error);

  if (strcmp(arg[0],"none") == 0) style = NONE;
  else if (strcmp(arg[0],"sc") == 0) style = SC;
  else if (strcmp(arg[0],"bcc") == 0) style = BCC;
  else if (strcmp(arg[0],"fcc") == 0) style = FCC;
  else if (strcmp(arg[0],"hcp") == 0) style = HCP;
  else if (strcmp(arg[0],"diamond") == 0) style = DIAMOND;
  else if (strcmp(arg[0],"sq") == 0) style = SQ;
  else if (strcmp(arg[0],"sq2") == 0) style = SQ2;
  else if (strcmp(arg[0],"hex") == 0) style = HEX;
  else if (strcmp(arg[0],"custom") == 0) style = CUSTOM;
  else error->all(FLERR,"Unknown lattice keyword: {}", arg[0]);

  if (style == NONE) {
    if (narg != 2) error->all(FLERR,"Illegal lattice command: expected 2 arguments but found {}", narg);

    xlattice = ylattice = zlattice = utils::numeric(FLERR,arg[1],false,lmp);
    if (xlattice <= 0.0) error->all(FLERR, "Invalid lattice none argument: {}", arg[1]);
    return;
  }

  // check that lattice matches dimension
  // style CUSTOM can be either 2d or 3d

  int dimension = domain->dimension;
  if (dimension == 2) {
    if (style == SC || style == BCC || style == FCC || style == HCP ||
        style == DIAMOND)
      error->all(FLERR,"Lattice style incompatible with simulation dimension");
  }
  if (dimension == 3) {
    if (style == SQ || style == SQ2 || style == HEX)
      error->all(FLERR,"Lattice style incompatible with simulation dimension");
  }

  // scale = conversion factor between lattice and box units

  if (narg < 2) utils::missing_cmd_args(FLERR, "lattice", error);
  scale = utils::numeric(FLERR,arg[1],false,lmp);
  if (scale <= 0.0) error->all(FLERR, "Invalid lattice {} argument: {}", arg[0], arg[1]);

  // set basis atoms for each style
  // x,y,z = fractional coords within unit cell
  // style CUSTOM will be defined by optional args

  if (style == SC) {
    add_basis(0.0,0.0,0.0);
  } else if (style == BCC) {
    add_basis(0.0,0.0,0.0);
    add_basis(0.5,0.5,0.5);
  } else if (style == FCC) {
    add_basis(0.0,0.0,0.0);
    add_basis(0.5,0.5,0.0);
    add_basis(0.5,0.0,0.5);
    add_basis(0.0,0.5,0.5);
  } else if (style == HCP) {
    add_basis(0.0,0.0,0.0);
    add_basis(0.5,0.5,0.0);
    add_basis(0.5,5.0/6.0,0.5);
    add_basis(0.0,1.0/3.0,0.5);
  } else if (style == SQ) {
    add_basis(0.0,0.0,0.0);
  } else if (style == SQ2) {
    add_basis(0.0,0.0,0.0);
    add_basis(0.5,0.5,0.0);
  } else if (style == HEX) {
    add_basis(0.0,0.0,0.0);
    add_basis(0.5,0.5,0.0);
  } else if (style == DIAMOND) {
    add_basis(0.0,0.0,0.0);
    add_basis(0.0,0.5,0.5);
    add_basis(0.5,0.0,0.5);
    add_basis(0.5,0.5,0.0);
    add_basis(0.25,0.25,0.25);
    add_basis(0.25,0.75,0.75);
    add_basis(0.75,0.25,0.75);
    add_basis(0.75,0.75,0.25);
  }

  // set defaults for optional args

  origin[0] = origin[1] = origin[2] = 0.0;

  orientx[0] = 1;  orientx[1] = 0;  orientx[2] = 0;
  orienty[0] = 0;  orienty[1] = 1;  orienty[2] = 0;
  orientz[0] = 0;  orientz[1] = 0;  orientz[2] = 1;

  int spaceflag = 0;

  a1[0] = 1.0;  a1[1] = 0.0;  a1[2] = 0.0;
  a2[0] = 0.0;  a2[1] = 1.0;  a2[2] = 0.0;
  a3[0] = 0.0;  a3[1] = 0.0;  a3[2] = 1.0;

  if (style == HEX) a2[1] = sqrt(3.0);
  if (style == HCP) {
    a2[1] = sqrt(3.0);
    a3[2] = sqrt(8.0/3.0);
  }

  // process optional args

  triclinic_general = 0;
  oriented = 0;

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"origin") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "lattice origin", error);
      origin[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      origin[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      origin[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if (origin[0] < 0.0 || origin[0] >= 1.0)
        error->all(FLERR, "Invalid lattice origin argument: {}", origin[0]);
      if (origin[1] < 0.0 || origin[1] >= 1.0)
        error->all(FLERR, "Invalid lattice origin argument: {}", origin[1]);
      if (origin[2] < 0.0 || origin[2] >= 1.0)
        error->all(FLERR, "Invalid lattice origin argument: {}", origin[2]);
      iarg += 4;

    } else if (strcmp(arg[iarg],"orient") == 0) {
      if (iarg+5 > narg) utils::missing_cmd_args(FLERR, "lattice orient", error);
      int dim = -1;
      if (strcmp(arg[iarg+1],"x") == 0) dim = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) dim = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) dim = 2;
      else error->all(FLERR,"Unknown lattice orient argument: {}", arg[iarg+1]);
      int *orient = nullptr;
      if (dim == 0) orient = orientx;
      else if (dim == 1) orient = orienty;
      else if (dim == 2) orient = orientz;
      orient[0] = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      orient[1] = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
      orient[2] = utils::inumeric(FLERR,arg[iarg+4],false,lmp);
      iarg += 5;

    } else if (strcmp(arg[iarg],"spacing") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "lattice spacing", error);
      spaceflag = 1;
      xlattice = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      ylattice = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      zlattice = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;

    } else if (strcmp(arg[iarg],"a1") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "lattice a1", error);
      if (style != CUSTOM)
        error->all(FLERR,
                   "Invalid a1 option in lattice command for non-custom style");
      a1[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      a1[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      a1[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;

    } else if (strcmp(arg[iarg],"a2") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "lattice a2", error);
      if (style != CUSTOM)
        error->all(FLERR,
                   "Invalid a2 option in lattice command for non-custom style");
      a2[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      a2[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      a2[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg],"a3") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "lattice a3", error);
      if (style != CUSTOM)
        error->all(FLERR,
                   "Invalid a3 option in lattice command for non-custom style");
      a3[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      a3[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      a3[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;

    } else if (strcmp(arg[iarg],"basis") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "lattice basis", error);
      if (style != CUSTOM)
        error->all(FLERR,
                   "Invalid basis option in lattice command for non-custom style");
      double x = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      double y = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      double z = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if (x < 0.0 || x >= 1.0)
        error->all(FLERR, "Invalid lattice basis argument: {}", x);
      if (y < 0.0 || y >= 1.0)
        error->all(FLERR, "Invalid lattice basis argument: {}", y);
      if (z < 0.0 || z >= 1.0)
        error->all(FLERR, "Invalid lattice basis argument: {}", z);
      add_basis(x,y,z);
      iarg += 4;

    } else if (strcmp(arg[iarg],"triclinic/general") == 0) {
      triclinic_general = 1;
      iarg++;

    } else error->all(FLERR,"Unknown lattice keyword: {}", arg[iarg]);
  }

  // check settings for errors

  if (nbasis == 0) error->all(FLERR,"No basis atoms in lattice");
  if (!orthogonal())
    error->all(FLERR,"Lattice orient vectors are not orthogonal");
  if (!right_handed_orientation())
    error->all(FLERR,"Lattice orient vectors are not right-handed");
  if (collinear())
    error->all(FLERR,"Lattice primitive vectors are collinear");

  // requirements for 2d system

  if (dimension == 2) {
    if (origin[2] != 0.0)
      error->all(FLERR,
                 "Lattice origin z coord must be 0.0 for 2d simulation");
    if (a1[2] != 0.0 || a2[2] != 0.0 || a3[0] != 0.0 || a3[1] != 0.0)
      error->all(FLERR,
                 "Lattice a1/a2/a3 vectors are not compatible with 2d simulation");
    if (orientx[2] != 0 || orienty[2] != 0 ||
        orientz[0] != 0 || orientz[1] != 0)
      error->all(FLERR,
                 "Lattice orient vectors are not compatible with 2d simulation");
    for (int i = 0; i < nbasis; i++)
      if (basis[i][2] != 0.0)
        error->all(FLERR,"Lattice basis atom z coords must be zero for 2d simulation");
  }

  // additional requirements for a general triclinic lattice
  // a123 prime are used to compute lattice spacings

  if (triclinic_general) {
    if (style != CUSTOM)
      error->all(FLERR,"Lattice triclinic/general must be style = CUSTOM");
    if (origin[0] != 0.0 || origin[1] != 0.0 || origin[2] != 0.0)
      error->all(FLERR,"Lattice triclinic/general must have default origin");
    int oriented = 0;
    if (orientx[0] != 1 || orientx[1] != 0 || orientx[2] != 0) oriented = 1;
    if (orienty[0] != 0 || orienty[1] != 1 || orienty[2] != 0) oriented = 1;
    if (orientz[0] != 0 || orientz[1] != 0 || orientz[2] != 1) oriented = 1;
    if (oriented)
      error->all(FLERR,"Lattice triclinic/general must have default orientation");
    if (dimension == 2 && (a3[0] != 0.0 || a3[1] != 0.0 || a3[2] != 1.0))
      error->all(FLERR,"Lattice triclinic/general a3 vector for a 2d simulation must be (0,0,1)");
    if (!right_handed_primitive())
      error->all(FLERR,"Lattice triclinic/general a1,a2,a3 must be right-handed");

    double rotmat[3][3];
    domain->general_to_restricted_rotation(a1,a2,a3,rotmat,a1_prime,a2_prime,a3_prime);
  }

  // user-defined lattice spacings must all be positive

  if (spaceflag) {
    if (xlattice <= 0.0 || ylattice <= 0.0 || zlattice <= 0.0)
      error->all(FLERR,"Lattice spacings are invalid");
  }

  // reset scale for LJ units (input scale is rho*)
  // scale = (Nbasis/(Vprimitive * rho*)) ^ (1/dim)
  // use fabs() in case a1,a2,a3 are not right-handed for general triclinic

  if (strcmp(update->unit_style,"lj") == 0) {
    double vec[3];
    MathExtra::cross3(a2,a3,vec);
    double volume = fabs(MathExtra::dot3(a1,vec));
    scale = pow(nbasis/volume/scale,1.0/dimension);
  }

  // initialize lattice <-> box transformation matrices

  setup_transform(a1,a2,a3);

  // automatic calculation of lattice spacings
  // convert 8 corners of primitive unit cell from lattice coords to box coords
  // min to max = bounding box around the pts in box space
  // xlattice,ylattice,zlattice = extent of bbox in box space
  // set xlattice,ylattice,zlattice to 0.0 initially
  //   since bbox uses them to shift origin (irrelevant for this computation)

  if (spaceflag == 0) {
    double xmin,ymin,zmin,xmax,ymax,zmax;
    xmin = ymin = zmin = BIG;
    xmax = ymax = zmax = -BIG;
    xlattice = ylattice = zlattice = 0.0;

    // for general triclinic, bounding box is around unit cell
    //   in restricted triclinic orientation, NOT general
    // this enables lattice spacings to be used for other commands (e.g. region)
    //   after create_box and create_atoms create the restricted triclnic system
    // reset transform used by bbox() to be based on rotated a123 prime vectors

    if (triclinic_general) setup_transform(a1_prime,a2_prime,a3_prime);

    bbox(0,0.0,0.0,0.0,xmin,ymin,zmin,xmax,ymax,zmax);
    bbox(0,1.0,0.0,0.0,xmin,ymin,zmin,xmax,ymax,zmax);
    bbox(0,0.0,1.0,0.0,xmin,ymin,zmin,xmax,ymax,zmax);
    bbox(0,1.0,1.0,0.0,xmin,ymin,zmin,xmax,ymax,zmax);
    bbox(0,0.0,0.0,1.0,xmin,ymin,zmin,xmax,ymax,zmax);
    bbox(0,1.0,0.0,1.0,xmin,ymin,zmin,xmax,ymax,zmax);
    bbox(0,0.0,1.0,1.0,xmin,ymin,zmin,xmax,ymax,zmax);
    bbox(0,1.0,1.0,1.0,xmin,ymin,zmin,xmax,ymax,zmax);

    // restore original general triclinic a123 transform

    if (triclinic_general) setup_transform(a1,a2,a3);

    xlattice = xmax - xmin;
    ylattice = ymax - ymin;
    zlattice = zmax - zmin;

  // user-defined lattice spacings

  } else {
    xlattice *= scale;
    ylattice *= scale;
    zlattice *= scale;
  }

  // print lattice spacings

  if (comm->me == 0)
    utils::logmesg(lmp,"Lattice spacing in x,y,z = {:.8} {:.8} {:.8}\n",
                   xlattice,ylattice,zlattice);
}

/* ---------------------------------------------------------------------- */

Lattice::~Lattice()
{
  memory->destroy(basis);
}

/* ----------------------------------------------------------------------
   return 1 if lattice is for a general triclinic simulation box
   queried by create_box and create_atoms
------------------------------------------------------------------------- */

int Lattice::is_general_triclinic()
{
  if (triclinic_general) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   check if 3 orientation vectors are mutually orthogonal
------------------------------------------------------------------------- */

int Lattice::orthogonal()
{
  if (orientx[0]*orienty[0] + orientx[1]*orienty[1] +
      orientx[2]*orienty[2]) return 0;
  if (orienty[0]*orientz[0] + orienty[1]*orientz[1] +
      orienty[2]*orientz[2]) return 0;
  if (orientx[0]*orientz[0] + orientx[1]*orientz[1] +
      orientx[2]*orientz[2]) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   check righthandedness of orientation vectors
   x cross y must be in same direction as z
------------------------------------------------------------------------- */

int Lattice::right_handed_orientation()
{
  int xy0 = orientx[1]*orienty[2] - orientx[2]*orienty[1];
  int xy1 = orientx[2]*orienty[0] - orientx[0]*orienty[2];
  int xy2 = orientx[0]*orienty[1] - orientx[1]*orienty[0];
  if (xy0*orientz[0] + xy1*orientz[1] + xy2*orientz[2] <= 0) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   check righthandedness of a1,a2,a3 primitive vectors
   x cross y must be in same direction as z
------------------------------------------------------------------------- */

int Lattice::right_handed_primitive()
{
  double vec[3];
  MathExtra::cross3(a1,a2,vec);
  if (MathExtra::dot3(vec,a3) <= 0.0) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   check collinearity of each pair of primitive vectors
   also checks if any primitive vector is zero-length
------------------------------------------------------------------------- */

int Lattice::collinear()
{
  double vec[3];
  MathExtra::cross3(a1,a2,vec);
  if (MathExtra::len3(vec) == 0.0) return 1;
  MathExtra::cross3(a2,a3,vec);
  if (MathExtra::len3(vec) == 0.0) return 1;
  MathExtra::cross3(a1,a3,vec);
  if (MathExtra::len3(vec) == 0.0) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   initialize lattice <-> box transformation matrices
------------------------------------------------------------------------- */

void Lattice::setup_transform(double *a, double *b, double *c)
{
  double length;

  // primitive = 3x3 matrix with primitive vectors as columns

  primitive[0][0] = a[0];
  primitive[1][0] = a[1];
  primitive[2][0] = a[2];
  primitive[0][1] = b[0];
  primitive[1][1] = b[1];
  primitive[2][1] = b[2];
  primitive[0][2] = c[0];
  primitive[1][2] = c[1];
  primitive[2][2] = c[2];

  // priminv = inverse of primitive

  double determinant = primitive[0][0]*primitive[1][1]*primitive[2][2] +
    primitive[0][1]*primitive[1][2]*primitive[2][0] +
    primitive[0][2]*primitive[1][0]*primitive[2][1] -
    primitive[0][0]*primitive[1][2]*primitive[2][1] -
    primitive[0][1]*primitive[1][0]*primitive[2][2] -
    primitive[0][2]*primitive[1][1]*primitive[2][0];

  if (determinant == 0.0)
    error->all(FLERR,"Degenerate lattice primitive vectors");

  priminv[0][0] = (primitive[1][1]*primitive[2][2] -
                   primitive[1][2]*primitive[2][1]) / determinant;
  priminv[1][0] = (primitive[1][2]*primitive[2][0] -
                   primitive[1][0]*primitive[2][2]) / determinant;
  priminv[2][0] = (primitive[1][0]*primitive[2][1] -
                   primitive[1][1]*primitive[2][0]) / determinant;

  priminv[0][1] = (primitive[0][2]*primitive[2][1] -
                   primitive[0][1]*primitive[2][2]) / determinant;
  priminv[1][1] = (primitive[0][0]*primitive[2][2] -
                   primitive[0][2]*primitive[2][0]) / determinant;
  priminv[2][1] = (primitive[0][1]*primitive[2][0] -
                   primitive[0][0]*primitive[2][1]) / determinant;

  priminv[0][2] = (primitive[0][1]*primitive[1][2] -
                   primitive[0][2]*primitive[1][1]) / determinant;
  priminv[1][2] = (primitive[0][2]*primitive[1][0] -
                   primitive[0][0]*primitive[1][2]) / determinant;
  priminv[2][2] = (primitive[0][0]*primitive[1][1] -
                   primitive[0][1]*primitive[1][0]) / determinant;

  // rotaterow = 3x3 matrix with normalized orient vectors as rows

  int lensq = orientx[0]*orientx[0] + orientx[1]*orientx[1] +
    orientx[2]*orientx[2];
  length = sqrt((double) lensq);
  if (length == 0.0) error->all(FLERR,"Zero-length lattice orient vector");

  rotaterow[0][0] = orientx[0] / length;
  rotaterow[0][1] = orientx[1] / length;
  rotaterow[0][2] = orientx[2] / length;

  lensq = orienty[0]*orienty[0] + orienty[1]*orienty[1] +
    orienty[2]*orienty[2];
  length = sqrt((double) lensq);
  if (length == 0.0) error->all(FLERR,"Zero-length lattice orient vector");

  rotaterow[1][0] = orienty[0] / length;
  rotaterow[1][1] = orienty[1] / length;
  rotaterow[1][2] = orienty[2] / length;

  lensq = orientz[0]*orientz[0] + orientz[1]*orientz[1] +
    orientz[2]*orientz[2];
  length = sqrt((double) lensq);
  if (length == 0.0) error->all(FLERR,"Zero-length lattice orient vector");

  rotaterow[2][0] = orientz[0] / length;
  rotaterow[2][1] = orientz[1] / length;
  rotaterow[2][2] = orientz[2] / length;

  // rotatecol = 3x3 matrix with normalized orient vectors as columns

  rotatecol[0][0] = rotaterow[0][0];
  rotatecol[1][0] = rotaterow[0][1];
  rotatecol[2][0] = rotaterow[0][2];

  rotatecol[0][1] = rotaterow[1][0];
  rotatecol[1][1] = rotaterow[1][1];
  rotatecol[2][1] = rotaterow[1][2];

  rotatecol[0][2] = rotaterow[2][0];
  rotatecol[1][2] = rotaterow[2][1];
  rotatecol[2][2] = rotaterow[2][2];
}

/* ----------------------------------------------------------------------
   convert lattice coords to box coords
   input x,y,z = point in lattice coords
   output x,y,z = point in box coords
   transformation: xyz_box = Rotate_row * scale * P * xyz_lattice + offset
     xyz_box = 3-vector of output box coords
     Rotate_row = 3x3 matrix = normalized orient vectors as rows
     scale = scale factor
     P = 3x3 matrix = primitive vectors as columns
     xyz_lattice = 3-vector of input lattice coords
     offset = 3-vector = (xlatt*origin[0], ylatt*origin[1], zlatt*origin[2])
------------------------------------------------------------------------- */

void Lattice::lattice2box(double &x, double &y, double &z)
{
  double x1 = primitive[0][0]*x + primitive[0][1]*y + primitive[0][2]*z;
  double y1 = primitive[1][0]*x + primitive[1][1]*y + primitive[1][2]*z;
  double z1 = primitive[2][0]*x + primitive[2][1]*y + primitive[2][2]*z;

  x1 *= scale;
  y1 *= scale;
  z1 *= scale;

  double xnew = rotaterow[0][0]*x1 + rotaterow[0][1]*y1 + rotaterow[0][2]*z1;
  double ynew = rotaterow[1][0]*x1 + rotaterow[1][1]*y1 + rotaterow[1][2]*z1;
  double znew = rotaterow[2][0]*x1 + rotaterow[2][1]*y1 + rotaterow[2][2]*z1;

  x = xnew + xlattice*origin[0];
  y = ynew + ylattice*origin[1];
  z = znew + zlattice*origin[2];
}

/* ----------------------------------------------------------------------
   convert box coords to lattice coords
   input x,y,z = point in box coords
   output x,y,z = point in lattice coords
   transformation: xyz_latt = P_inv * 1/scale * Rotate_col * (xyz_box - offset)
     xyz_lattice = 3-vector of output lattice coords
     P_inv = 3x3 matrix = inverse of primitive vectors as columns
     scale = scale factor
     Rotate_col = 3x3 matrix = normalized orient vectors as columns
     xyz_box = 3-vector of input box coords
     offset = 3-vector = (xlatt*origin[0], ylatt*origin[1], zlatt*origin[2])
------------------------------------------------------------------------- */

void Lattice::box2lattice(double &x, double &y, double &z)
{
  x -= xlattice*origin[0];
  y -= ylattice*origin[1];
  z -= zlattice*origin[2];

  double x1 = rotatecol[0][0]*x + rotatecol[0][1]*y + rotatecol[0][2]*z;
  double y1 = rotatecol[1][0]*x + rotatecol[1][1]*y + rotatecol[1][2]*z;
  double z1 = rotatecol[2][0]*x + rotatecol[2][1]*y + rotatecol[2][2]*z;

  x1 /= scale;
  y1 /= scale;
  z1 /= scale;

  x = priminv[0][0]*x1 + priminv[0][1]*y1 + priminv[0][2]*z1;
  y = priminv[1][0]*x1 + priminv[1][1]*y1 + priminv[1][2]*z1;
  z = priminv[2][0]*x1 + priminv[2][1]*y1 + priminv[2][2]*z1;
}

/* ----------------------------------------------------------------------
   add a basis atom to list
   x,y,z = fractional coords within unit cell
------------------------------------------------------------------------- */

void Lattice::add_basis(double x, double y, double z)
{
  memory->grow(basis,nbasis+1,3,"lattice:basis");
  basis[nbasis][0] = x;
  basis[nbasis][1] = y;
  basis[nbasis][2] = z;
  nbasis++;
}

/* ----------------------------------------------------------------------
   convert x,y,z from lattice coords to box coords (flag = 0)
   or from box coords to lattice coords (flag = 1)
   either way, use new point to expand bounding box (min to max)
------------------------------------------------------------------------- */

void Lattice::bbox(int flag, double x, double y, double z,
                   double &xmin, double &ymin, double &zmin,
                   double &xmax, double &ymax, double &zmax)
{
  if (flag == 0) lattice2box(x,y,z);
  else box2lattice(x,y,z);

  xmin = MIN(x,xmin);  ymin = MIN(y,ymin);  zmin = MIN(z,zmin);
  xmax = MAX(x,xmax);  ymax = MAX(y,ymax);  zmax = MAX(z,zmax);
}
