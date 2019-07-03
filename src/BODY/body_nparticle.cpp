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

#include "body_nparticle.h"
#include <cstring>
#include <cstdlib>
#include "my_pool_chunk.h"
#include "math_extra.h"
#include "atom_vec_body.h"
#include "atom.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define EPSILON 1.0e-7
enum{SPHERE,LINE,TRI};           // also in DumpImage

/* ---------------------------------------------------------------------- */

BodyNparticle::BodyNparticle(LAMMPS *lmp, int narg, char **arg) :
  Body(lmp, narg, arg), imflag(NULL), imdata(NULL)
{
  if (narg != 3) error->all(FLERR,"Invalid body nparticle command");

  int nmin = force->inumeric(FLERR,arg[1]);
  int nmax = force->inumeric(FLERR,arg[2]);
  if (nmin <= 0 || nmin > nmax)
    error->all(FLERR,"Invalid body nparticle command");

  size_forward = 0;
  size_border = 1 + 3*nmax;

  // NOTE: need to set appropriate nnbin param for dcp

  icp = new MyPoolChunk<int>(1,1);
  dcp = new MyPoolChunk<double>(3*nmin,3*nmax);

  memory->create(imflag,nmax,"body/nparticle:imflag");
  memory->create(imdata,nmax,4,"body/nparticle:imdata");
}

/* ---------------------------------------------------------------------- */

BodyNparticle::~BodyNparticle()
{
  delete icp;
  delete dcp;
  memory->destroy(imflag);
  memory->destroy(imdata);
}

/* ---------------------------------------------------------------------- */

int BodyNparticle::nsub(AtomVecBody::Bonus *bonus)
{
  return bonus->ivalue[0];
}

/* ---------------------------------------------------------------------- */

double *BodyNparticle::coords(AtomVecBody::Bonus *bonus)
{
  return bonus->dvalue;
}

/* ---------------------------------------------------------------------- */

int BodyNparticle::pack_border_body(AtomVecBody::Bonus *bonus, double *buf)
{
  int nsub = bonus->ivalue[0];
  buf[0] = nsub;
  memcpy(&buf[1],bonus->dvalue,3*nsub*sizeof(double));
  return 1+3*nsub;
}

/* ---------------------------------------------------------------------- */

int BodyNparticle::unpack_border_body(AtomVecBody::Bonus *bonus, double *buf)
{
  int nsub = static_cast<int> (buf[0]);
  bonus->ivalue[0] = nsub;
  memcpy(bonus->dvalue,&buf[1],3*nsub*sizeof(double));
  return 1+3*nsub;
}

/* ----------------------------------------------------------------------
   populate bonus data structure with data file values
------------------------------------------------------------------------- */

void BodyNparticle::data_body(int ibonus, int ninteger, int ndouble,
                              int *ifile, double *dfile)
{
  AtomVecBody::Bonus *bonus = &avec->bonus[ibonus];

  // set ninteger, ndouble in bonus and allocate 2 vectors of ints, doubles

  if (ninteger != 1)
    error->one(FLERR,"Incorrect # of integer values in "
               "Bodies section of data file");
  int nsub = ifile[0];
  if (nsub < 1)
    error->one(FLERR,"Incorrect integer value in "
               "Bodies section of data file");
  if (ndouble != 6 + 3*nsub)
    error->one(FLERR,"Incorrect # of floating-point values in "
               "Bodies section of data file");

  bonus->ninteger = 1;
  bonus->ivalue = icp->get(bonus->iindex);
  bonus->ivalue[0] = nsub;
  bonus->ndouble = 3*nsub;
  bonus->dvalue = dcp->get(bonus->ndouble,bonus->dindex);

  // diagonalize inertia tensor

  double tensor[3][3];
  tensor[0][0] = dfile[0];
  tensor[1][1] = dfile[1];
  tensor[2][2] = dfile[2];
  tensor[0][1] = tensor[1][0] = dfile[3];
  tensor[0][2] = tensor[2][0] = dfile[4];
  tensor[1][2] = tensor[2][1] = dfile[5];

  double *inertia = bonus->inertia;
  double evectors[3][3];
  int ierror = MathExtra::jacobi(tensor,inertia,evectors);
  if (ierror) error->one(FLERR,
                         "Insufficient Jacobi rotations for body nparticle");

  // if any principal moment < scaled EPSILON, set to 0.0

  double max;
  max = MAX(inertia[0],inertia[1]);
  max = MAX(max,inertia[2]);

  if (inertia[0] < EPSILON*max) inertia[0] = 0.0;
  if (inertia[1] < EPSILON*max) inertia[1] = 0.0;
  if (inertia[2] < EPSILON*max) inertia[2] = 0.0;

  // exyz_space = principal axes in space frame

  double ex_space[3],ey_space[3],ez_space[3];

  ex_space[0] = evectors[0][0];
  ex_space[1] = evectors[1][0];
  ex_space[2] = evectors[2][0];
  ey_space[0] = evectors[0][1];
  ey_space[1] = evectors[1][1];
  ey_space[2] = evectors[2][1];
  ez_space[0] = evectors[0][2];
  ez_space[1] = evectors[1][2];
  ez_space[2] = evectors[2][2];

  // enforce 3 evectors as a right-handed coordinate system
  // flip 3rd vector if needed

  double cross[3];
  MathExtra::cross3(ex_space,ey_space,cross);
  if (MathExtra::dot3(cross,ez_space) < 0.0) MathExtra::negate3(ez_space);

  // create initial quaternion

  MathExtra::exyz_to_q(ex_space,ey_space,ez_space,bonus->quat);

  // bonus->dvalue = sub-particle displacements in body frame

  double delta[3];

  int j = 6;
  int k = 0;
  for (int i = 0; i < nsub; i++) {
    delta[0] = dfile[j];
    delta[1] = dfile[j+1];
    delta[2] = dfile[j+2];
    MathExtra::transpose_matvec(ex_space,ey_space,ez_space,
                                delta,&bonus->dvalue[k]);
    j += 3;
    k += 3;
  }
}

/* ----------------------------------------------------------------------
   return radius of body particle defined by ifile/dfile params
   params are ordered as in data file
   called by Molecule class which needs single body size
------------------------------------------------------------------------- */

double BodyNparticle::radius_body(int /*ninteger*/, int ndouble,
                                  int *ifile, double *dfile)
{
  int nsub = ifile[0];
  if (nsub < 1)
    error->one(FLERR,"Incorrect integer value in "
               "Bodies section of data file");
  if (ndouble != 6 + 3*nsub)
    error->one(FLERR,"Incorrect # of floating-point values in "
               "Bodies section of data file");

  // sub-particle coords are relative to body center at (0,0,0)
  // offset = 6 for sub-particle coords

  double onerad;
  double maxrad = 0.0;
  double delta[3];

  int offset = 6;
  for (int i = 0; i < nsub; i++) {
    delta[0] = dfile[offset];
    delta[1] = dfile[offset+1];
    delta[2] = dfile[offset+2];
    offset += 3;
    onerad = MathExtra::len3(delta);
    maxrad = MAX(maxrad,onerad);
  }

  return maxrad;
}

/* ---------------------------------------------------------------------- */

int BodyNparticle::noutcol()
{
  return 3;
}

/* ---------------------------------------------------------------------- */

int BodyNparticle::noutrow(int ibonus)
{
  return avec->bonus[ibonus].ivalue[0];
}

/* ---------------------------------------------------------------------- */

void BodyNparticle::output(int ibonus, int m, double *values)
{
  AtomVecBody::Bonus *bonus = &avec->bonus[ibonus];

  double p[3][3];
  MathExtra::quat_to_mat(bonus->quat,p);
  MathExtra::matvec(p,&bonus->dvalue[3*m],values);

  double *x = atom->x[bonus->ilocal];
  values[0] += x[0];
  values[1] += x[1];
  values[2] += x[2];
}

/* ---------------------------------------------------------------------- */

int BodyNparticle::image(int ibonus, double flag1, double /*flag2*/,
                         int *&ivec, double **&darray)
{
  double p[3][3];
  double *x;

  AtomVecBody::Bonus *bonus = &avec->bonus[ibonus];
  int n = bonus->ivalue[0];

  for (int i = 0; i < n; i++) {
    imflag[i] = SPHERE;
    MathExtra::quat_to_mat(bonus->quat,p);
    MathExtra::matvec(p,&bonus->dvalue[3*i],imdata[i]);

    x = atom->x[bonus->ilocal];
    imdata[i][0] += x[0];
    imdata[i][1] += x[1];
    imdata[i][2] += x[2];
    imdata[i][3] = flag1;
  }

  ivec = imflag;
  darray = imdata;
  return n;
}
