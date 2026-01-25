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

#include "body_rounded_polyhedron.h"

#include "atom.h"
#include "atom_vec_body.h"
#include "error.h"
#include "graphics.h"
#include "math_extra.h"
#include "math_eigen.h"
#include "memory.h"
#include "my_pool_chunk.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

static constexpr double EPSILON = 1.0e-7;
static constexpr int MAX_FACE_SIZE = 4;  // maximum number of vertices per face (for now)

/* ---------------------------------------------------------------------- */

BodyRoundedPolyhedron::BodyRoundedPolyhedron(LAMMPS *lmp, int narg, char **arg) :
  Body(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Invalid body rounded/polygon command");

  // nmin and nmax are minimum and maximum number of vertices

  int nmin = utils::inumeric(FLERR,arg[1],false,lmp);
  int nmax = utils::inumeric(FLERR,arg[2],false,lmp);
  if (nmin <= 0 || nmin > nmax) error->all(FLERR,"Invalid body rounded/polyhedron command");

  size_forward = 0;

  // 3 integers: nvertices, nedges, nfaces
  // doubles = 3*nvertices + 2*nedge + MAX_FACE_SIZE*nfaces +
  //           1 double for enclosing radius + 1 double for rounded radius

  size_border = 3 + 3*nmax + 2*nmax + MAX_FACE_SIZE*nmax + 1 + 1;

  // NOTE: need to set appropriate nnbin param for dcp

  icp = new MyPoolChunk<int>(1,3);
  dcp = new MyPoolChunk<double>(3*nmin+2+1+1,
                                3*nmax+2*nmax+MAX_FACE_SIZE*nmax+1+1);
  maxexchange = 3 + 3*nmax+2*nmax+MAX_FACE_SIZE*nmax+1+1;  // icp max + dcp max

  memory->create(imflag,3*nmax,"body/rounded/polyhedron:imflag");
  memory->create(imdata,3*nmax,9,"body/rounded/polyhedron:imdata");
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
  if (nvertices == 1) return 0;
  else if (nvertices == 2) return 1;
  return nedges;
}

/* ---------------------------------------------------------------------- */

double *BodyRoundedPolyhedron::edges(AtomVecBody::Bonus *bonus)
{
  return bonus->dvalue+3*nsub(bonus);
}

/* ---------------------------------------------------------------------- */

int BodyRoundedPolyhedron::nfaces(AtomVecBody::Bonus *bonus)
{
  int nvertices = bonus->ivalue[0];
  if (nvertices < 3) return 0;

  return bonus->ivalue[2];
}

/* ---------------------------------------------------------------------- */

double *BodyRoundedPolyhedron::faces(AtomVecBody::Bonus *bonus)
{
  int nvertices = bonus->ivalue[0];
  if (nvertices == 1 || nvertices == 2) return nullptr;
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

int BodyRoundedPolyhedron::unpack_border_body(AtomVecBody::Bonus *bonus, double *buf)
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
   populate bonus data structure with data file values for one body
------------------------------------------------------------------------- */

void BodyRoundedPolyhedron::data_body(int ibonus, int ninteger, int ndouble,
                                      int *ifile, double *dfile)
{
  AtomVecBody::Bonus *bonus = &avec->bonus[ibonus];

  // set ninteger, ndouble in bonus and allocate 2 vectors of ints, doubles

  if (ninteger != 3)
    error->one(FLERR,"Incorrect # of integer values in Bodies section of data file");
  int nsub = ifile[0];
  int ned = ifile[1];
  int nfac = ifile[2];
  if (nsub < 1) error->one(FLERR,"Incorrect integer value in Bodies section of data file");

  // nentries = number of double entries to be read from Body section:
  // nsub == 1,2:
  //   6 for inertia + 3*nsub + 1 for rounded radius
  // nsub > 2:
  //   6 for inertia + 3*nsub + 2*nedges + MAX_FACE_SIZE*nfaces + 1 for rounded radius

  int nedges,nentries;
  if (nsub < 3) {
    nentries = 6 + 3*nsub + 1;
  } else {
    nedges = ned; //nsub + nfac - 2;
    nentries = 6 + 3*nsub + 2*nedges + MAX_FACE_SIZE*nfac + 1;
  }

  if (ndouble != nentries)
    error->one(FLERR,"Incorrect # of floating-point values in Bodies section of data file");

  bonus->ninteger = 3;
  bonus->ivalue = icp->get(bonus->iindex);
  bonus->ivalue[0] = nsub;
  bonus->ivalue[1] = ned;
  bonus->ivalue[2] = nfac;
  if (nsub < 3) bonus->ndouble = 3*nsub + 2 + 1 + 1;
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
  int ierror = MathEigen::jacobi3(tensor,inertia,evectors);
  if (ierror) error->one(FLERR,"Insufficient Jacobi rotations for body nparticle");

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

  // the next 2*nsub elements are edge ends
  // the final two values are the enclosing radius and rounded radius
  // set atom->radius = enclosing + rounded radii (except for spheres)

  // spheres have just 1 edge

  if (nsub == 1) {
    bonus->dvalue[k] = 0;
    bonus->dvalue[k+1] = 0;
    k += 2;

    rrad = 0.5 * dfile[j];
    bonus->dvalue[k] = rrad;
    erad = rrad;

    k++;
    bonus->dvalue[k] = rrad;

    atom->radius[bonus->ilocal] = erad;

  // rods have just 1 edge

  } else if (nsub == 2) {
    bonus->dvalue[k] = 0;
    bonus->dvalue[k+1] = 1;
    k += 2;

    erad = sqrt(erad2);
    bonus->dvalue[k] = erad;

    rrad = 0.5 * dfile[j];
    k++;
    bonus->dvalue[k] = rrad;

    atom->radius[bonus->ilocal] = erad + rrad;

  // polyhedra have Nedges and Nfaces

  } else {

    // edges

    for (i = 0; i < nedges; i++) {
      bonus->dvalue[k] = dfile[j];
      bonus->dvalue[k+1] = dfile[j+1];
      k += 2;
      j += 2;
    }

    // faces

    for (i = 0; i < nfac; i++) {
      for (m = 0; m < MAX_FACE_SIZE; m++)
        bonus->dvalue[k+m] = dfile[j+m];
      k += MAX_FACE_SIZE;
      j += MAX_FACE_SIZE;
    }

    erad = sqrt(erad2);
    bonus->dvalue[k] = erad;

    rrad = 0.5 * dfile[j];
    k++;
    bonus->dvalue[k] = rrad;

    atom->radius[bonus->ilocal] = erad + rrad;
  }
}

/* ----------------------------------------------------------------------
   pack data struct for one body into buf for writing to data file
   if buf is a null pointer, just return buffer size
------------------------------------------------------------------------- */

int BodyRoundedPolyhedron::pack_data_body(tagint atomID, int ibonus, double *buf)
{
  int i,j,m;
  double values[3],p[3][3],pdiag[3][3],ispace[3][3];

  AtomVecBody::Bonus *bonus = &avec->bonus[ibonus];

  double *quat = bonus->quat;
  double *inertia = bonus->inertia;
  int *ivalue = bonus->ivalue;
  double *dvalue = bonus->dvalue;

  int nsub = ivalue[0];
  int nedge = ivalue[1];
  int nface = ivalue[2];

  if (buf) {

    // ID ninteger ndouble

    m = 0;
    buf[m++] = ubuf(atomID).d;
    buf[m++] = ubuf(3).d;
    if (nsub < 3) buf[m++] = ubuf(6 + 3*nsub + 1).d;
    else buf[m++] = ubuf(6 + 3*nsub + 2*nedge + MAX_FACE_SIZE*nface + 1).d;

    // 3 integers nsub,nedge,nface

    buf[m++] = ubuf(nsub).d;
    buf[m++] = ubuf(nedge).d;
    buf[m++] = ubuf(nface).d;

    // 6 moments of inertia

    MathExtra::quat_to_mat(quat,p);
    MathExtra::times3_diag(p,inertia,pdiag);
    MathExtra::times3_transpose(pdiag,p,ispace);

    buf[m++] = ispace[0][0];
    buf[m++] = ispace[1][1];
    buf[m++] = ispace[2][2];
    buf[m++] = ispace[0][1];
    buf[m++] = ispace[0][2];
    buf[m++] = ispace[1][2];

    // 3*nsub particle coords = displacement from COM in box frame

    for (i = 0; i < nsub; i++) {
      MathExtra::matvec(p,&dvalue[3*i],values);
      buf[m++] = values[0];
      buf[m++] = values[1];
      buf[m++] = values[2];
    }

    // 2*nedge edge indices
    // 4*nface face indices

    j = 3*nsub;

    if (nsub < 3) j += 2;
    else {
      for (i = 0; i < nedge; i++) {
        buf[m++] = static_cast<int> (dvalue[j++]);
        buf[m++] = static_cast<int> (dvalue[j++]);
      }
      for (i = 0; i < nface; i++) {
        buf[m++] = static_cast<int> (dvalue[j++]);
        buf[m++] = static_cast<int> (dvalue[j++]);
        buf[m++] = static_cast<int> (dvalue[j++]);
        buf[m++] = static_cast<int> (dvalue[j++]);
      }
    }

    // rounded diameter = 2 * last dvalue = rounded radius
    // j+1 to skip enclosing radius

    buf[m++] = 2.0 * dvalue[j+1];

  } else {
    m = 3 + 3 + 6 + 3*nsub + 1;
    if (nsub > 2) m += 2*nedge + MAX_FACE_SIZE*nface;
  }

  return m;
}

/* ----------------------------------------------------------------------
   write info for one body to data file
------------------------------------------------------------------------- */

int BodyRoundedPolyhedron::write_data_body(FILE *fp, double *buf)
{
  int m = 0;

  // atomID ninteger ndouble

  utils::print(fp,"{} {} {}\n",ubuf(buf[m]).i,ubuf(buf[m+1]).i,ubuf(buf[m+2]).i);
  m += 3;

  // nvert, nedge, nface

  const int nsub = (int) ubuf(buf[m++]).i;
  const int nedge = (int) ubuf(buf[m++]).i;
  const int nface = (int) ubuf(buf[m++]).i;
  utils::print(fp,"{} {} {}\n",nsub,nedge,nface);

  // inertia

  utils::print(fp,"{} {} {} {} {} {}\n",
             buf[m+0],buf[m+1],buf[m+2],buf[m+3],buf[m+4],buf[m+5]);
  m += 6;

  // nsub vertices

  for (int i = 0; i < nsub; i++, m+=3)
    utils::print(fp,"{} {} {}\n",buf[m],buf[m+1],buf[m+2]);

  // nedge 2-tuples and nface 4-tuples
  // unless nsub = 1 or 2

  if (nsub > 2) {
    for (int i = 0; i < nedge; i++, m+=2)
      utils::print(fp,"{} {}\n",static_cast<int> (buf[m]),static_cast<int> (buf[m+1]));
    for (int i = 0; i < nface; i++, m+=4)
      utils::print(fp,"{} {} {} {}\n",
                 static_cast<int> (buf[m]),static_cast<int> (buf[m+1]),
                 static_cast<int> (buf[m+2]),static_cast<int> (buf[m+3]));
  }

  // rounded diameter

  double diameter = buf[m++];
  utils::print(fp,"{}\n",diameter);

  return m;
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
    error->one(FLERR,"Incorrect integer value in Bodies section of data file");
  if (ndouble != nentries)
    error->one(FLERR,"Incorrect # of floating-point values in Bodies section of data file");

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

int BodyRoundedPolyhedron::image(int ibonus, double flag1, double flag2,
                                 int *&ivec, double **&darray)
{
  double p[3][3];

  int nelements = 0;
  AtomVecBody::Bonus *const bonus = &avec->bonus[ibonus];
  const double *const x = atom->x[bonus->ilocal];
  const double diam = (flag1 <= 0.0) ? 2.0 * rounded_radius(bonus) : flag1;

  MathExtra::quat_to_mat(bonus->quat,p); // get rotation matrix for body frame to box frame

  int nvertices = bonus->ivalue[0];
  if (nvertices == 1) { // special case: just one vertex -> one sphere
    imflag[0] = Graphics::SPHERE;
    // transform body frame position to box frame
    MathExtra::matvec(p,&bonus->dvalue[0],imdata[0]);
    // translate and set diameter
    imdata[0][0] += x[0];
    imdata[0][1] += x[1];
    imdata[0][2] += x[2];
    imdata[0][3] = diam;

    nelements = 1;
  } else {

    // select whether edges or faces or both should be drawn
    bool edgeflag = true;
    bool triflag = true;
    int flag = static_cast<int>(flag2);
    if (flag == 2) triflag = false;
    if (flag == 1) edgeflag = false;

    int nedges = bonus->ivalue[1];
    if (nvertices == 2) nedges = 1;                  // special case: just two vertices -> one rod
    double *edge_ends = &bonus->dvalue[3*nvertices]; // skip over vertex positions in body data
    if (edgeflag || (nedges == 1)) {                 // always draw edge for rod
      for (int i = 0; i < nedges; i++) {
        imflag[nelements] = Graphics::LINE;

        int pt1 = static_cast<int>(edge_ends[2*i]);
        int pt2 = static_cast<int>(edge_ends[2*i+1]);

        // transform body frame positions to box frame
        MathExtra::matvec(p,&bonus->dvalue[3*pt1],imdata[nelements]);
        MathExtra::matvec(p,&bonus->dvalue[3*pt2],&imdata[nelements][3]);

        // translate and set diameter
        imdata[nelements][0] += x[0];
        imdata[nelements][1] += x[1];
        imdata[nelements][2] += x[2];
        imdata[nelements][3] += x[0];
        imdata[nelements][4] += x[1];
        imdata[nelements][5] += x[2];
        imdata[nelements][6] = diam;

        ++nelements;
      }
    }

    if (triflag) {
      int nfaces = bonus->ivalue[2];
      // skip over vertex positions and edge indices in body data
      edge_ends = &bonus->dvalue[3*nvertices+2*nedges];
      for (int i = 0; i < nfaces; i++) {
        int pt1 = static_cast<int>(edge_ends[4*i]);
        int pt2 = static_cast<int>(edge_ends[4*i+1]);
        int pt3 = static_cast<int>(edge_ends[4*i+2]);
        int pt4 = static_cast<int>(edge_ends[4*i+3]);

        // quadrilateral face requires two triangles. triangle has fourth vertex index set to -1
        if (pt4 >= 0) {
          // first triangle
          imflag[nelements] = Graphics::TRI;
          MathExtra::matvec(p,&bonus->dvalue[3*pt1],imdata[nelements]);
          MathExtra::matvec(p,&bonus->dvalue[3*pt2],&imdata[nelements][3]);
          MathExtra::matvec(p,&bonus->dvalue[3*pt3],&imdata[nelements][6]);
          imdata[nelements][0] += x[0];
          imdata[nelements][1] += x[1];
          imdata[nelements][2] += x[2];
          imdata[nelements][3] += x[0];
          imdata[nelements][4] += x[1];
          imdata[nelements][5] += x[2];
          imdata[nelements][6] += x[0];
          imdata[nelements][7] += x[1];
          imdata[nelements][8] += x[2];
          ++nelements;

          // second triangle
          imflag[nelements] = Graphics::TRI;
          MathExtra::matvec(p,&bonus->dvalue[3*pt3],imdata[nelements]);
          MathExtra::matvec(p,&bonus->dvalue[3*pt4],&imdata[nelements][3]);
          MathExtra::matvec(p,&bonus->dvalue[3*pt1],&imdata[nelements][6]);
          imdata[nelements][0] += x[0];
          imdata[nelements][1] += x[1];
          imdata[nelements][2] += x[2];
          imdata[nelements][3] += x[0];
          imdata[nelements][4] += x[1];
          imdata[nelements][5] += x[2];
          imdata[nelements][6] += x[0];
          imdata[nelements][7] += x[1];
          imdata[nelements][8] += x[2];

          // shift triangles toward the outside of the body by half diameter when also drawing edges
          if (edgeflag) {
            double vec1[3];
            // get center of face
            double vec2[3] = {imdata[nelements][0],imdata[nelements][1],imdata[nelements][2]};
            vec2[0] += imdata[nelements][3] + imdata[nelements][6] + imdata[nelements-1][3];
            vec2[1] += imdata[nelements][4] + imdata[nelements][7] + imdata[nelements-1][4];
            vec2[2] += imdata[nelements][5] + imdata[nelements][8] + imdata[nelements-1][5];
            vec2[0] *= 0.25;
            vec2[1] *= 0.25;
            vec2[2] *= 0.25;
            // get direction from center of body to face and scale to half diameter length
            MathExtra::sub3(vec2,x,vec1);
            MathExtra::snormalize3(0.5*diam,vec1,vec2);
            // add shift to triangle corners
            imdata[nelements][0] += vec2[0];
            imdata[nelements][1] += vec2[1];
            imdata[nelements][2] += vec2[2];
            imdata[nelements][3] += vec2[0];
            imdata[nelements][4] += vec2[1];
            imdata[nelements][5] += vec2[2];
            imdata[nelements][6] += vec2[0];
            imdata[nelements][7] += vec2[1];
            imdata[nelements][8] += vec2[2];
            imdata[nelements-1][0] += vec2[0];
            imdata[nelements-1][1] += vec2[1];
            imdata[nelements-1][2] += vec2[2];
            imdata[nelements-1][3] += vec2[0];
            imdata[nelements-1][4] += vec2[1];
            imdata[nelements-1][5] += vec2[2];
            imdata[nelements-1][6] += vec2[0];
            imdata[nelements-1][7] += vec2[1];
            imdata[nelements-1][8] += vec2[2];
          }

          ++nelements;
        } else {
          imflag[nelements] = Graphics::TRI;
          MathExtra::matvec(p,&bonus->dvalue[3*pt1],imdata[nelements]);
          MathExtra::matvec(p,&bonus->dvalue[3*pt2],&imdata[nelements][3]);
          MathExtra::matvec(p,&bonus->dvalue[3*pt3],&imdata[nelements][6]);
          imdata[nelements][0] += x[0];
          imdata[nelements][1] += x[1];
          imdata[nelements][2] += x[2];
          imdata[nelements][3] += x[0];
          imdata[nelements][4] += x[1];
          imdata[nelements][5] += x[2];
          imdata[nelements][6] += x[0];
          imdata[nelements][7] += x[1];
          imdata[nelements][8] += x[2];

          // shift triangle toward the outside of the body by half diameter when also drawing edges
          if (edgeflag) {
            double vec1[3];
            // get center of face
            double vec2[3] = {imdata[nelements][0],imdata[nelements][1],imdata[nelements][2]};
            vec2[0] += imdata[nelements][3] + imdata[nelements][6];
            vec2[1] += imdata[nelements][4] + imdata[nelements][7];
            vec2[2] += imdata[nelements][5] + imdata[nelements][8];
            vec2[0] *= 1.0/3.0;
            vec2[1] *= 1.0/3.0;
            vec2[2] *= 1.0/3.0;
            // get direction from center of body to face and scale to half diameter length
            MathExtra::sub3(vec2,x,vec1);
            MathExtra::snormalize3(0.5*diam,vec1,vec2);
            // add shift to triangle corners
            imdata[nelements][0] += vec2[0];
            imdata[nelements][1] += vec2[1];
            imdata[nelements][2] += vec2[2];
            imdata[nelements][3] += vec2[0];
            imdata[nelements][4] += vec2[1];
            imdata[nelements][5] += vec2[2];
            imdata[nelements][6] += vec2[0];
            imdata[nelements][7] += vec2[1];
            imdata[nelements][8] += vec2[2];
          }
          ++nelements;
        }
      }
    }
  }

  ivec = imflag;
  darray = imdata;
  return nelements;
}
