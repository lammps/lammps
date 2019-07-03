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

#include "body_rounded_polyhedron.h"
#include <cstring>
#include <cstdlib>
#include "my_pool_chunk.h"
#include "atom_vec_body.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define EPSILON 1.0e-7
#define MAX_FACE_SIZE 4  // maximum number of vertices per face (for now)

enum{SPHERE,LINE};       // also in DumpImage

/* ---------------------------------------------------------------------- */

BodyRoundedPolyhedron::BodyRoundedPolyhedron(LAMMPS *lmp, int narg, char **arg) :
  Body(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Invalid body rounded/polygon command");

  // nmin and nmax are minimum and maximum number of vertices

  int nmin = force->inumeric(FLERR,arg[1]);
  int nmax = force->inumeric(FLERR,arg[2]);
  if (nmin <= 0 || nmin > nmax)
    error->all(FLERR,"Invalid body rounded/polyhedron command");

  size_forward = 0;

  // 3 integers: 1 for no. of vertices, 1 for no. of edges, 1 for no. of faces
  // 3*nmax doubles for vertex coordinates + 2*nmax doubles for edge ends +
  // (MAX_FACE_SIZE+1)*nmax for faces
  // 1 double for the enclosing radius
  // 1 double for the rounded radius

  size_border = 3 + 3*nmax + 2*nmax + MAX_FACE_SIZE*nmax + 1 + 1;

  // NOTE: need to set appropriate nnbin param for dcp

  icp = new MyPoolChunk<int>(1,3);
  dcp = new MyPoolChunk<double>(3*nmin+2+1+1,
                                3*nmax+2*nmax+MAX_FACE_SIZE*nmax+1+1);

  memory->create(imflag,2*nmax,"body/rounded/polyhedron:imflag");
  memory->create(imdata,2*nmax,7,"body/polyhedron:imdata");
}

/* ---------------------------------------------------------------------- */

BodyRoundedPolyhedron::~BodyRoundedPolyhedron()
{
  delete icp;
  delete dcp;
  memory->destroy(imflag);
  memory->destroy(imdata);
}

/* ---------------------------------------------------------------------- */

int BodyRoundedPolyhedron::nsub(AtomVecBody::Bonus *bonus)
{
  return bonus->ivalue[0];
}

/* ---------------------------------------------------------------------- */

double *BodyRoundedPolyhedron::coords(AtomVecBody::Bonus *bonus)
{
  return bonus->dvalue;
}

/* ---------------------------------------------------------------------- */

int BodyRoundedPolyhedron::nedges(AtomVecBody::Bonus *bonus)
{
  int nvertices = bonus->ivalue[0];
  int nedges = bonus->ivalue[1];
  //int nfaces = bonus->ivalue[2];
  if (nvertices == 1) return 0;
  else if (nvertices == 2) return 1;
  return nedges; //(nvertices+nfaces-2); // Euler's polyon formula: V-E+F=2
}

/* ---------------------------------------------------------------------- */

double *BodyRoundedPolyhedron::edges(AtomVecBody::Bonus *bonus)
{
  return bonus->dvalue+3*nsub(bonus);
}

/* ---------------------------------------------------------------------- */

int BodyRoundedPolyhedron::nfaces(AtomVecBody::Bonus *bonus)
{
  return bonus->ivalue[2];
}

/* ---------------------------------------------------------------------- */

double *BodyRoundedPolyhedron::faces(AtomVecBody::Bonus *bonus)
{
  int nvertices = bonus->ivalue[0];
  if (nvertices == 1 || nvertices == 2) return NULL;
  return bonus->dvalue+3*nsub(bonus)+2*nedges(bonus);
}

/* ---------------------------------------------------------------------- */

double BodyRoundedPolyhedron::enclosing_radius(struct AtomVecBody::Bonus *bonus)
{
  int nvertices = bonus->ivalue[0];
  if (nvertices == 1 || nvertices == 2)
  	return *(bonus->dvalue+3*nsub(bonus)+2);
  return *(bonus->dvalue+3*nsub(bonus) + 2*nedges(bonus) +
           MAX_FACE_SIZE*nfaces(bonus));
}

/* ---------------------------------------------------------------------- */

double BodyRoundedPolyhedron::rounded_radius(struct AtomVecBody::Bonus *bonus)
{
  int nvertices = bonus->ivalue[0];
  if (nvertices == 1 || nvertices == 2)
    return *(bonus->dvalue+3*nsub(bonus)+2+1);
  return *(bonus->dvalue+3*nsub(bonus) + 2*nedges(bonus) +
           MAX_FACE_SIZE*nfaces(bonus)+1);
}

/* ---------------------------------------------------------------------- */

int BodyRoundedPolyhedron::pack_border_body(AtomVecBody::Bonus *bonus, double *buf)
{
  int nsub = bonus->ivalue[0];
  int ned = bonus->ivalue[1];
  int nfac = bonus->ivalue[2];
  buf[0] = nsub;
  buf[1] = ned;
  buf[2] = nfac;
  int ndouble;
  if (nsub == 1 || nsub == 2) ndouble = 3*nsub+2+MAX_FACE_SIZE*nfac+1+1;
  else ndouble = 3*nsub+2*ned+MAX_FACE_SIZE*nfac+1+1;
  memcpy(&buf[3],bonus->dvalue,ndouble*sizeof(double));
  return 3+ndouble;
}

/* ---------------------------------------------------------------------- */

int BodyRoundedPolyhedron::unpack_border_body(AtomVecBody::Bonus *bonus,
                                           double *buf)
{
  int nsub = static_cast<int> (buf[0]);
  int ned = static_cast<int> (buf[1]);
  int nfac = static_cast<int> (buf[2]);
  bonus->ivalue[0] = nsub;
  bonus->ivalue[1] = ned;
  bonus->ivalue[2] = nfac;
  int ndouble;
  if (nsub == 1 || nsub == 2) ndouble = 3*nsub+2+MAX_FACE_SIZE*nfac+1+1;
  else ndouble = 3*nsub+2*ned+MAX_FACE_SIZE*nfac+1+1;
  memcpy(bonus->dvalue,&buf[3],ndouble*sizeof(double));
  return 3+ndouble;
}

/* ----------------------------------------------------------------------
   populate bonus data structure with data file values
------------------------------------------------------------------------- */

void BodyRoundedPolyhedron::data_body(int ibonus, int ninteger, int ndouble,
                             int *ifile, double *dfile)
{
  AtomVecBody::Bonus *bonus = &avec->bonus[ibonus];

  // set ninteger, ndouble in bonus and allocate 2 vectors of ints, doubles

  if (ninteger != 3)
    error->one(FLERR,"Incorrect # of integer values in "
               "Bodies section of data file");
  int nsub = ifile[0];
  int ned = ifile[1];
  int nfac = ifile[2];
  if (nsub < 1)
    error->one(FLERR,"Incorrect integer value in "
               "Bodies section of data file");

  // nentries = number of double entries to be read from Body section:
  // nsub == 1 || nsub == 2 || nsub == 3:
  //   6 for inertia + 3*nsub for vertex coords + 1 for rounded radius
  // nsub > 3:
  //   6 for inertia + 3*nsub for vertex coords + 2*nsub for edges +
  //   3*nfaces + 1 for rounded radius

  int nedges,nentries;
  if (nsub == 1 || nsub == 2) {
    nentries = 6 + 3*nsub + 1;
  } else {
    nedges = ned; //nsub + nfac - 2;
    nentries = 6 + 3*nsub + 2*nedges + MAX_FACE_SIZE*nfac + 1;
  }
  if (ndouble != nentries)
    error->one(FLERR,"Incorrect # of floating-point values in "
             "Bodies section of data file");

  bonus->ninteger = 3;
  bonus->ivalue = icp->get(bonus->iindex);
  bonus->ivalue[0] = nsub;
  bonus->ivalue[1] = ned;
  bonus->ivalue[2] = nfac;
  if (nsub == 1 || nsub == 2) bonus->ndouble = 3*nsub + 2*nsub + 1 + 1;
  else bonus->ndouble = 3*nsub + 2*nedges + MAX_FACE_SIZE*nfac + 1 + 1;
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

  // bonus->dvalue = the first 3*nsub elements are sub-particle displacements
  // find the enclosing radius of the body from the maximum displacement

  int i,m;
  double delta[3], rsq, erad, rrad;
  double erad2 = 0;
  int j = 6;
  int k = 0;
  for (i = 0; i < nsub; i++) {
    delta[0] = dfile[j];
    delta[1] = dfile[j+1];
    delta[2] = dfile[j+2];
    MathExtra::transpose_matvec(ex_space,ey_space,ez_space,
                                delta,&bonus->dvalue[k]);
    rsq = delta[0] * delta[0] + delta[1] * delta[1] +
      delta[2] * delta[2];
    if (rsq > erad2) erad2 = rsq;
    j += 3;
    k += 3;
  }

  // .. the next 2*nsub elements are edge ends

  if (nsub == 1) { // spheres
    nedges = 0;
    bonus->dvalue[k] = 0;
    *(&bonus->dvalue[k]+1) = 0;
    k += 2;

    rrad = 0.5 * dfile[j];
    bonus->dvalue[k] = rrad;
    erad = rrad; // enclosing radius = rounded_radius

    // the last element of bonus->dvalue is the rounded radius

    k++;
    bonus->dvalue[k] = rrad;

    atom->radius[bonus->ilocal] = erad;

  } else if (nsub == 2) { // rods
    nedges = 1;
    for (i = 0; i < nedges; i++) {
      bonus->dvalue[k] = 0;
      *(&bonus->dvalue[k]+1) = 1;
      k += 2;
    }

    erad = sqrt(erad2);
    bonus->dvalue[k] = erad;

    // the last element of bonus->dvalue is the rounded radius

    rrad = 0.5 * dfile[j];
    k++;
    bonus->dvalue[k] = rrad;

    atom->radius[bonus->ilocal] = erad + rrad;

  } else { // polyhedra

    // edges

    for (i = 0; i < nedges; i++) {
      bonus->dvalue[k] = dfile[j];
      *(&bonus->dvalue[k]+1) = dfile[j+1];
      k += 2;
      j += 2;
    }

    // faces

    for (i = 0; i < nfac; i++) {
      for (m = 0; m < MAX_FACE_SIZE; m++)
        *(&bonus->dvalue[k]+m) = dfile[j+m];
      k += MAX_FACE_SIZE;
      j += MAX_FACE_SIZE;
    }

    // the next to last element is the enclosing radius

    erad = sqrt(erad2);
    bonus->dvalue[k] = erad;

    // the last element bonus-> dvalue is the rounded radius

    rrad = 0.5 * dfile[j];
    k++;
    bonus->dvalue[k] = rrad;

    atom->radius[bonus->ilocal] = erad + rrad;
  }
}

/* ----------------------------------------------------------------------
   return radius of body particle defined by ifile/dfile params
   params are ordered as in data file
   called by Molecule class which needs single body size
------------------------------------------------------------------------- */

double BodyRoundedPolyhedron::radius_body(int /*ninteger*/, int ndouble,
				       int *ifile, double *dfile)
{
  int nsub = ifile[0];
  int ned = ifile[1];
  int nfac = ifile[2];
  int nedges = ned; //nsub + nfac - 2;

  int nentries;
  if (nsub == 1 || nsub == 2) nentries = 6 + 3*nsub + 1;
  else nentries = 6 + 3*nsub + 2*nedges + MAX_FACE_SIZE*nfac + 1;

  if (nsub < 1)
    error->one(FLERR,"Incorrect integer value in "
               "Bodies section of data file");
  if (ndouble != nentries)
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

  if (nsub > 2) offset += (2*nedges+MAX_FACE_SIZE*nfac);

  // add in radius of rounded corners

  return maxrad + 0.5*dfile[offset];
}

/* ---------------------------------------------------------------------- */

int BodyRoundedPolyhedron::noutcol()
{
  // the number of columns for the vertex coordinates

  return 3;
}

/* ---------------------------------------------------------------------- */

int BodyRoundedPolyhedron::noutrow(int ibonus)
{
  // only return the first nsub rows for the vertex coordinates

  return avec->bonus[ibonus].ivalue[0];
}

/* ---------------------------------------------------------------------- */

void BodyRoundedPolyhedron::output(int ibonus, int m, double *values)
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

int BodyRoundedPolyhedron::image(int ibonus, double flag1, double /*flag2*/,
                              int *&ivec, double **&darray)
{
  int nelements;
  double p[3][3];
  double *x, rrad;

  AtomVecBody::Bonus *bonus = &avec->bonus[ibonus];
  int nvertices = bonus->ivalue[0];

  if (nvertices == 1) { // spheres

    for (int i = 0; i < nvertices; i++) {
      imflag[i] = SPHERE;
      MathExtra::quat_to_mat(bonus->quat,p);
      MathExtra::matvec(p,&bonus->dvalue[3*i],imdata[i]);

      rrad = enclosing_radius(bonus);
      x = atom->x[bonus->ilocal];
      imdata[i][0] += x[0];
      imdata[i][1] += x[1];
      imdata[i][2] += x[2];
      if (flag1 <= 0) imdata[i][3] = 2*rrad;
      else imdata[i][3] = flag1;
    }

    nelements = nvertices;
  } else {
    //int nfaces = bonus->ivalue[2];
    int nedges = bonus->ivalue[1]; //nvertices + nfaces - 2;
    if (nvertices == 2) nedges = 1; // special case: rods
    double* edge_ends = &bonus->dvalue[3*nvertices];
    int pt1, pt2;

    for (int i = 0; i < nedges; i++) {
      imflag[i] = LINE;

      pt1 = static_cast<int>(edge_ends[2*i]);
      pt2 = static_cast<int>(edge_ends[2*i+1]);

      MathExtra::quat_to_mat(bonus->quat,p);
      MathExtra::matvec(p,&bonus->dvalue[3*pt1],imdata[i]);
      MathExtra::matvec(p,&bonus->dvalue[3*pt2],&imdata[i][3]);

      rrad = rounded_radius(bonus);
      x = atom->x[bonus->ilocal];
      imdata[i][0] += x[0];
      imdata[i][1] += x[1];
      imdata[i][2] += x[2];
      imdata[i][3] += x[0];
      imdata[i][4] += x[1];
      imdata[i][5] += x[2];

      if (flag1 <= 0) imdata[i][6] = 2*rrad;
      else imdata[i][6] = flag1;
    }

    nelements = nedges;
  }

  ivec = imflag;
  darray = imdata;
  return nelements;
}
