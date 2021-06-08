// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_vec_tri.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "math_const.h"
#include "math_extra.h"
#include "math_eigen.h"
#include "memory.h"
#include "modify.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

#define EPSILON 0.001

/* ---------------------------------------------------------------------- */

AtomVecTri::AtomVecTri(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  bonus_flag = 1;

  size_forward_bonus = 4;
  size_border_bonus = 17;
  size_restart_bonus_one = 17;
  size_data_bonus = 10;

  atom->tri_flag = 1;
  atom->molecule_flag = atom->rmass_flag = 1;
  atom->radius_flag = atom->omega_flag = atom->angmom_flag = 1;
  atom->torque_flag = 1;
  atom->sphere_flag = 1;

  nlocal_bonus = nghost_bonus = nmax_bonus = 0;
  bonus = nullptr;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *) "molecule radius rmass omega angmom torque tri";
  fields_copy = (char *) "molecule radius rmass omega angmom";
  fields_comm = (char *) "";
  fields_comm_vel = (char *) "omega angmom";
  fields_reverse = (char *) "torque";
  fields_border = (char *) "molecule radius rmass";
  fields_border_vel = (char *) "molecule radius rmass omega";
  fields_exchange = (char *) "molecule radius rmass omega angmom";
  fields_restart = (char *) "molecule radius rmass omega angmom";
  fields_create = (char *) "molecule radius rmass omega angmom tri";
  fields_data_atom = (char *) "id molecule type tri rmass x";
  fields_data_vel = (char *) "id v omega angmom";

  setup_fields();
}

/* ---------------------------------------------------------------------- */

AtomVecTri::~AtomVecTri()
{
  memory->sfree(bonus);
}

/* ---------------------------------------------------------------------- */

void AtomVecTri::init()
{
  AtomVec::init();

  if (domain->dimension != 3)
    error->all(FLERR,"Atom_style tri can only be used in 3d simulations");
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecTri::grow_pointers()
{
  tri = atom->tri;
  radius = atom->radius;
  rmass = atom->rmass;
  omega = atom->omega;
  angmom = atom->angmom;
}

/* ----------------------------------------------------------------------
   grow bonus data structure
------------------------------------------------------------------------- */

void AtomVecTri::grow_bonus()
{
  nmax_bonus = grow_nmax_bonus(nmax_bonus);
  if (nmax_bonus < 0)
    error->one(FLERR,"Per-processor system is too big");

  bonus = (Bonus *) memory->srealloc(bonus,nmax_bonus*sizeof(Bonus),
                                     "atom:bonus");
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
   if delflag and atom J has bonus data, then delete it
------------------------------------------------------------------------- */

void AtomVecTri::copy_bonus(int i, int j, int delflag)
{
  // if deleting atom J via delflag and J has bonus data, then delete it

  if (delflag && tri[j] >= 0) {
    copy_bonus_all(nlocal_bonus-1,tri[j]);
    nlocal_bonus--;
  }

  // if atom I has bonus data, reset I's bonus.ilocal to loc J
  // do NOT do this if self-copy (I=J) since I's bonus data is already deleted

  if (tri[i] >= 0 && i != j) bonus[tri[i]].ilocal = j;
  tri[j] = tri[i];
}

/* ----------------------------------------------------------------------
   copy bonus data from I to J, effectively deleting the J entry
   also reset tri that points to I to now point to J
------------------------------------------------------------------------- */

void AtomVecTri::copy_bonus_all(int i, int j)
{
  tri[bonus[i].ilocal] = j;
  memcpy(&bonus[j],&bonus[i],sizeof(Bonus));
}

/* ----------------------------------------------------------------------
   clear ghost info in bonus data
   called before ghosts are recommunicated in comm and irregular
------------------------------------------------------------------------- */

void AtomVecTri::clear_bonus()
{
  nghost_bonus = 0;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->clear_bonus();
}

/* ---------------------------------------------------------------------- */

int AtomVecTri::pack_comm_bonus(int n, int *list, double *buf)
{
  int i,j,m;
  double *quat;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (tri[j] >= 0) {
      quat = bonus[tri[j]].quat;
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecTri::unpack_comm_bonus(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (tri[i] >= 0) {
      quat = bonus[tri[i]].quat;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecTri::pack_border_bonus(int n, int *list, double *buf)
{
  int i,j,m;
  double *quat,*c1,*c2,*c3,*inertia;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (tri[j] < 0) buf[m++] = ubuf(0).d;
    else {
      buf[m++] = ubuf(1).d;
      quat = bonus[tri[j]].quat;
      c1 = bonus[tri[j]].c1;
      c2 = bonus[tri[j]].c2;
      c3 = bonus[tri[j]].c3;
      inertia = bonus[tri[j]].inertia;
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
      buf[m++] = c1[0];
      buf[m++] = c1[1];
      buf[m++] = c1[2];
      buf[m++] = c2[0];
      buf[m++] = c2[1];
      buf[m++] = c2[2];
      buf[m++] = c3[0];
      buf[m++] = c3[1];
      buf[m++] = c3[2];
      buf[m++] = inertia[0];
      buf[m++] = inertia[1];
      buf[m++] = inertia[2];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecTri::unpack_border_bonus(int n, int first, double *buf)
{
  int i,j,m,last;
  double *quat,*c1,*c2,*c3,*inertia;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    tri[i] = (int) ubuf(buf[m++]).i;
    if (tri[i] == 0) tri[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      quat = bonus[j].quat;
      c1 = bonus[j].c1;
      c2 = bonus[j].c2;
      c3 = bonus[j].c3;
      inertia = bonus[j].inertia;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      c1[0] = buf[m++];
      c1[1] = buf[m++];
      c1[2] = buf[m++];
      c2[0] = buf[m++];
      c2[1] = buf[m++];
      c2[2] = buf[m++];
      c3[0] = buf[m++];
      c3[1] = buf[m++];
      c3[2] = buf[m++];
      inertia[0] = buf[m++];
      inertia[1] = buf[m++];
      inertia[2] = buf[m++];
      bonus[j].ilocal = i;
      tri[i] = j;
      nghost_bonus++;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecTri::pack_exchange_bonus(int i, double *buf)
{
  int m = 0;

  if (tri[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = tri[i];
    double *quat = bonus[j].quat;
    double *c1 = bonus[j].c1;
    double *c2 = bonus[j].c2;
    double *c3 = bonus[j].c3;
    double *inertia = bonus[j].inertia;
    buf[m++] = quat[0];
    buf[m++] = quat[1];
    buf[m++] = quat[2];
    buf[m++] = quat[3];
    buf[m++] = c1[0];
    buf[m++] = c1[1];
    buf[m++] = c1[2];
    buf[m++] = c2[0];
    buf[m++] = c2[1];
    buf[m++] = c2[2];
    buf[m++] = c3[0];
    buf[m++] = c3[1];
    buf[m++] = c3[2];
    buf[m++] = inertia[0];
    buf[m++] = inertia[1];
    buf[m++] = inertia[2];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecTri::unpack_exchange_bonus(int ilocal, double *buf)
{
  int m = 0;

  tri[ilocal] = (int) ubuf(buf[m++]).i;
  if (tri[ilocal] == 0) tri[ilocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *quat = bonus[nlocal_bonus].quat;
    double *c1 = bonus[nlocal_bonus].c1;
    double *c2 = bonus[nlocal_bonus].c2;
    double *c3 = bonus[nlocal_bonus].c3;
    double *inertia = bonus[nlocal_bonus].inertia;
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    c1[0] = buf[m++];
    c1[1] = buf[m++];
    c1[2] = buf[m++];
    c2[0] = buf[m++];
    c2[1] = buf[m++];
    c2[2] = buf[m++];
    c3[0] = buf[m++];
    c3[1] = buf[m++];
    c3[2] = buf[m++];
    inertia[0] = buf[m++];
    inertia[1] = buf[m++];
    inertia[2] = buf[m++];
    bonus[nlocal_bonus].ilocal = ilocal;
    tri[ilocal] = nlocal_bonus++;
  }

  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecTri::size_restart_bonus()
{
  int i;

  int n = 0;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) {
    if (tri[i] >= 0) n += size_restart_bonus_one;
    else n++;
  }

  return n;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including bonus data
------------------------------------------------------------------------- */

int AtomVecTri::pack_restart_bonus(int i, double *buf)
{
  int m = 0;

  if (tri[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = tri[i];
    double *quat = bonus[j].quat;
    double *c1 = bonus[j].c1;
    double *c2 = bonus[j].c2;
    double *c3 = bonus[j].c3;
    double *inertia = bonus[j].inertia;
    buf[m++] = quat[0];
    buf[m++] = quat[1];
    buf[m++] = quat[2];
    buf[m++] = quat[3];
    buf[m++] = c1[0];
    buf[m++] = c1[1];
    buf[m++] = c1[2];
    buf[m++] = c2[0];
    buf[m++] = c2[1];
    buf[m++] = c2[2];
    buf[m++] = c3[0];
    buf[m++] = c3[1];
    buf[m++] = c3[2];
    buf[m++] = inertia[0];
    buf[m++] = inertia[1];
    buf[m++] = inertia[2];
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including bonus data
------------------------------------------------------------------------- */

int AtomVecTri::unpack_restart_bonus(int ilocal, double *buf)
{
  int m = 0;

  tri[ilocal] = (int) ubuf(buf[m++]).i;
  if (tri[ilocal] == 0) tri[ilocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *quat = bonus[nlocal_bonus].quat;
    double *c1 = bonus[nlocal_bonus].c1;
    double *c2 = bonus[nlocal_bonus].c2;
    double *c3 = bonus[nlocal_bonus].c3;
    double *inertia = bonus[nlocal_bonus].inertia;
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    c1[0] = buf[m++];
    c1[1] = buf[m++];
    c1[2] = buf[m++];
    c2[0] = buf[m++];
    c2[1] = buf[m++];
    c2[2] = buf[m++];
    c3[0] = buf[m++];
    c3[1] = buf[m++];
    c3[2] = buf[m++];
    inertia[0] = buf[m++];
    inertia[1] = buf[m++];
    inertia[2] = buf[m++];
    bonus[nlocal_bonus].ilocal = ilocal;
    tri[ilocal] = nlocal_bonus++;
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack one line from Tris section of data file
------------------------------------------------------------------------- */

void AtomVecTri::data_atom_bonus(int m, char **values)
{
  if (tri[m]) error->one(FLERR,"Assigning tri parameters to non-tri atom");

  if (nlocal_bonus == nmax_bonus) grow_bonus();

  double c1[3],c2[3],c3[3];
  c1[0] = utils::numeric(FLERR,values[0],true,lmp);
  c1[1] = utils::numeric(FLERR,values[1],true,lmp);
  c1[2] = utils::numeric(FLERR,values[2],true,lmp);
  c2[0] = utils::numeric(FLERR,values[3],true,lmp);
  c2[1] = utils::numeric(FLERR,values[4],true,lmp);
  c2[2] = utils::numeric(FLERR,values[5],true,lmp);
  c3[0] = utils::numeric(FLERR,values[6],true,lmp);
  c3[1] = utils::numeric(FLERR,values[7],true,lmp);
  c3[2] = utils::numeric(FLERR,values[8],true,lmp);

  // check for duplicate points

  if (c1[0] == c2[0] && c1[1] == c2[1] && c1[2] == c2[2])
    error->one(FLERR,"Invalid shape in Triangles section of data file");
  if (c1[0] == c3[0] && c1[1] == c3[1] && c1[2] == c3[2])
    error->one(FLERR,"Invalid shape in Triangles section of data file");
  if (c2[0] == c3[0] && c2[1] == c3[1] && c2[2] == c3[2])
    error->one(FLERR,"Invalid shape in Triangles section of data file");

  // size = length of one edge

  double c2mc1[3],c3mc1[3];
  MathExtra::sub3(c2,c1,c2mc1);
  MathExtra::sub3(c3,c1,c3mc1);
  double size = MAX(MathExtra::len3(c2mc1),MathExtra::len3(c3mc1));

  // centroid = 1/3 of sum of vertices

  double centroid[3];
  centroid[0] = (c1[0]+c2[0]+c3[0]) / 3.0;
  centroid[1] = (c1[1]+c2[1]+c3[1]) / 3.0;
  centroid[2] = (c1[2]+c2[2]+c3[2]) / 3.0;

  double dx = centroid[0] - x[m][0];
  double dy = centroid[1] - x[m][1];
  double dz = centroid[2] - x[m][2];
  double delta = sqrt(dx*dx + dy*dy + dz*dz);

  if (delta/size > EPSILON)
    error->one(FLERR,"Inconsistent triangle in data file");

  x[m][0] = centroid[0];
  x[m][1] = centroid[1];
  x[m][2] = centroid[2];

  // reset tri radius and mass
  // rmass currently holds density
  // tri area = 0.5 len(U x V), where U,V are edge vectors from one vertex

  double c4[3];
  MathExtra::sub3(c1,centroid,c4);
  radius[m] = MathExtra::lensq3(c4);
  MathExtra::sub3(c2,centroid,c4);
  radius[m] = MAX(radius[m],MathExtra::lensq3(c4));
  MathExtra::sub3(c3,centroid,c4);
  radius[m] = MAX(radius[m],MathExtra::lensq3(c4));
  radius[m] = sqrt(radius[m]);

  double norm[3];
  MathExtra::cross3(c2mc1,c3mc1,norm);
  double area = 0.5 * MathExtra::len3(norm);
  rmass[m] *= area;

  // inertia = inertia tensor of triangle as 6-vector in Voigt ordering

  double inertia[6];
  MathExtra::inertia_triangle(c1,c2,c3,rmass[m],inertia);

  // diagonalize inertia tensor via Jacobi rotations
  // bonus[].inertia = 3 eigenvalues = principal moments of inertia
  // evectors and exzy_space = 3 evectors = principal axes of triangle

  double tensor[3][3],evectors[3][3];
  tensor[0][0] = inertia[0];
  tensor[1][1] = inertia[1];
  tensor[2][2] = inertia[2];
  tensor[1][2] = tensor[2][1] = inertia[3];
  tensor[0][2] = tensor[2][0] = inertia[4];
  tensor[0][1] = tensor[1][0] = inertia[5];

  int ierror = MathEigen::jacobi3(tensor,bonus[nlocal_bonus].inertia,evectors);
  if (ierror) error->one(FLERR,"Insufficient Jacobi rotations for triangle");

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

  // enforce 3 orthogonal vectors as a right-handed coordinate system
  // flip 3rd vector if needed

  MathExtra::cross3(ex_space,ey_space,norm);
  if (MathExtra::dot3(norm,ez_space) < 0.0) MathExtra::negate3(ez_space);

  // create initial quaternion

  MathExtra::exyz_to_q(ex_space,ey_space,ez_space,bonus[nlocal_bonus].quat);

  // bonus c1,c2,c3 = displacement of c1,c2,c3 from centroid
  // in basis of principal axes

  double disp[3];
  MathExtra::sub3(c1,centroid,disp);
  MathExtra::transpose_matvec(ex_space,ey_space,ez_space,
                              disp,bonus[nlocal_bonus].c1);
  MathExtra::sub3(c2,centroid,disp);
  MathExtra::transpose_matvec(ex_space,ey_space,ez_space,
                              disp,bonus[nlocal_bonus].c2);
  MathExtra::sub3(c3,centroid,disp);
  MathExtra::transpose_matvec(ex_space,ey_space,ez_space,
                              disp,bonus[nlocal_bonus].c3);

  bonus[nlocal_bonus].ilocal = m;
  tri[m] = nlocal_bonus++;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated bonus memory
------------------------------------------------------------------------- */

double AtomVecTri::memory_usage_bonus()
{
  double bytes = 0;
  bytes += (double)nmax_bonus*sizeof(Bonus);
  return bytes;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecTri::create_atom_post(int ilocal)
{
  double radius_one = 0.5;
  radius[ilocal] = radius_one;
  rmass[ilocal] = 4.0*MY_PI/3.0 * radius_one*radius_one*radius_one;
  tri[ilocal] = -1;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecTri::data_atom_post(int ilocal)
{
  tri_flag = tri[ilocal];
  if (tri_flag == 0) tri_flag = -1;
  else if (tri_flag == 1) tri_flag = 0;
  else error->one(FLERR,"Invalid tri flag in Atoms section of data file");
  tri[ilocal] = tri_flag;

  if (rmass[ilocal] <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  if (tri_flag < 0) {
    double radius_one = 0.5;
    radius[ilocal] = radius_one;
    rmass[ilocal] *= 4.0*MY_PI/3.0 * radius_one*radius_one*radius_one;
  } else radius[ilocal] = 0.0;

  omega[ilocal][0] = 0.0;
  omega[ilocal][1] = 0.0;
  omega[ilocal][2] = 0.0;
  angmom[ilocal][0] = 0.0;
  angmom[ilocal][1] = 0.0;
  angmom[ilocal][2] = 0.0;
}

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_data() to pack
------------------------------------------------------------------------- */

void AtomVecTri::pack_data_pre(int ilocal)
{
  tri_flag = tri[ilocal];
  rmass_one = rmass[ilocal];

  if (tri_flag < 0) tri[ilocal] = 0;
  else tri[ilocal] = 1;

  if (tri_flag < 0) {
    double radius_one = radius[ilocal];
    rmass[ilocal] /= 4.0*MY_PI/3.0 * radius_one*radius_one*radius_one;
  } else {
    double c2mc1[3],c3mc1[3],norm[3];
    MathExtra::sub3(bonus[tri_flag].c2,bonus[tri_flag].c1,c2mc1);
    MathExtra::sub3(bonus[tri_flag].c3,bonus[tri_flag].c1,c3mc1);
    MathExtra::cross3(c2mc1,c3mc1,norm);
    double area = 0.5 * MathExtra::len3(norm);
    rmass[ilocal] /= area;
  }
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_data()
------------------------------------------------------------------------- */

void AtomVecTri::pack_data_post(int ilocal)
{
  tri[ilocal] = tri_flag;
  rmass[ilocal] = rmass_one;
}

/* ----------------------------------------------------------------------
   pack bonus tri info for writing to data file
   if buf is nullptr, just return buffer size
------------------------------------------------------------------------- */

int AtomVecTri::pack_data_bonus(double *buf, int /*flag*/)
{
  int i,j;
  double xc,yc,zc;
  double dc1[3],dc2[3],dc3[3];
  double p[3][3];

  double **x = atom->x;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int m = 0;
  for (i = 0; i < nlocal; i++) {
    if (tri[i] < 0) continue;
    if (buf) {
      buf[m++] = ubuf(tag[i]).d;
      j = tri[i];
      MathExtra::quat_to_mat(bonus[j].quat,p);
      MathExtra::matvec(p,bonus[j].c1,dc1);
      MathExtra::matvec(p,bonus[j].c2,dc2);
      MathExtra::matvec(p,bonus[j].c3,dc3);
      xc = x[i][0];
      yc = x[i][1];
      zc = x[i][2];
      buf[m++] = xc + dc1[0];
      buf[m++] = yc + dc1[1];
      buf[m++] = zc + dc1[2];
      buf[m++] = xc + dc2[0];
      buf[m++] = yc + dc2[1];
      buf[m++] = zc + dc2[2];
      buf[m++] = xc + dc3[0];
      buf[m++] = yc + dc3[1];
      buf[m++] = zc + dc3[2];
    } else m += size_data_bonus;
  }
  return m;
}

/* ----------------------------------------------------------------------
   write bonus tri info to data file
------------------------------------------------------------------------- */

void AtomVecTri::write_data_bonus(FILE *fp, int n, double *buf, int /*flag*/)
{
  int i = 0;
  while (i < n) {
    fmt::print(fp,"{} {} {} {} {} {} {} {} {} {}\n", ubuf(buf[i]).i,
               buf[i+1],buf[i+2],buf[i+3],buf[i+4],buf[i+5],buf[i+6],
               buf[i+7],buf[i+8],buf[i+9]);
    i += size_data_bonus;
  }
}

/* ----------------------------------------------------------------------
   set equilateral tri of size in bonus data for particle I
   oriented symmetrically in xy plane
   this may create or delete entry in bonus data
------------------------------------------------------------------------- */

void AtomVecTri::set_equilateral(int i, double size)
{
  // also set radius = distance from center to corner-pt = len(c1)
  // unless size = 0.0, then set diameter = 1.0

  if (tri[i] < 0) {
    if (size == 0.0) return;
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *quat = bonus[nlocal_bonus].quat;
    double *c1 = bonus[nlocal_bonus].c1;
    double *c2 = bonus[nlocal_bonus].c2;
    double *c3 = bonus[nlocal_bonus].c3;
    double *inertia = bonus[nlocal_bonus].inertia;
    quat[0] = 1.0;
    quat[1] = 0.0;
    quat[2] = 0.0;
    quat[3] = 0.0;
    c1[0] = -size/2.0;
    c1[1] = -sqrt(3.0)/2.0 * size / 3.0;
    c1[2] = 0.0;
    c2[0] = size/2.0;
    c2[1] = -sqrt(3.0)/2.0 * size / 3.0;
    c2[2] = 0.0;
    c3[0] = 0.0;
    c3[1] = sqrt(3.0)/2.0 * size * 2.0/3.0;
    c3[2] = 0.0;
    inertia[0] = sqrt(3.0)/96.0 * size*size*size*size;
    inertia[1] = sqrt(3.0)/96.0 * size*size*size*size;
    inertia[2] = sqrt(3.0)/48.0 * size*size*size*size;
    radius[i] = MathExtra::len3(c1);
    bonus[nlocal_bonus].ilocal = i;
    tri[i] = nlocal_bonus++;
  } else if (size == 0.0) {
    radius[i] = 0.5;
    copy_bonus_all(nlocal_bonus-1,tri[i]);
    nlocal_bonus--;
    tri[i] = -1;
  } else {
    double *c1 = bonus[tri[i]].c1;
    double *c2 = bonus[tri[i]].c2;
    double *c3 = bonus[tri[i]].c3;
    double *inertia = bonus[tri[i]].inertia;
    c1[0] = -size/2.0;
    c1[1] = -sqrt(3.0)/2.0 * size / 3.0;
    c1[2] = 0.0;
    c2[0] = size/2.0;
    c2[1] = -sqrt(3.0)/2.0 * size / 3.0;
    c2[2] = 0.0;
    c3[0] = 0.0;
    c3[1] = sqrt(3.0)/2.0 * size * 2.0/3.0;
    c3[2] = 0.0;
    inertia[0] = sqrt(3.0)/96.0 * size*size*size*size;
    inertia[1] = sqrt(3.0)/96.0 * size*size*size*size;
    inertia[2] = sqrt(3.0)/48.0 * size*size*size*size;
    radius[i] = MathExtra::len3(c1);
  }
}
