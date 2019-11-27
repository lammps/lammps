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

#include "atom_vec_body.h"
#include <cstring>
#include <string>
#include "my_pool_chunk.h"
#include "style_body.h"
#include "body.h"
#include "atom.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecBody::AtomVecBody(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;
  bonus_flag = 1;

  // first 3 sizes do not include values from body itself
  // 1st,2nd body counts are added in process_args() via body style
  // 3rd body count is added in size_restart_bonus()
  // size_data_bonus is not used by Atom class for body style

  size_forward_bonus = 4;
  size_border_bonus = 9;
  size_restart_bonus_one = 9;
  size_data_bonus = 0;

  atom->body_flag = 1;
  atom->rmass_flag = 1;
  atom->angmom_flag = atom->torque_flag = 1;
  atom->radius_flag = 1;

  nlocal_bonus = nghost_bonus = nmax_bonus = 0;
  bonus = NULL;

  bptr = NULL;

  if (sizeof(double) == sizeof(int)) intdoubleratio = 1;
  else if (sizeof(double) == 2*sizeof(int)) intdoubleratio = 2;
  else error->all(FLERR,"Internal error in atom_style body");

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in the string does not matter
  //   except fields_data_atom and fields_data_vel which must match data file

  fields_grow = (char *) "radius rmass angmom torque body";
  fields_copy = (char *) "radius rmass angmom";
  fields_comm = NULL;
  fields_comm_vel = (char *) "angmom";
  fields_reverse = (char *) "torque";
  fields_border = (char *) "radius rmass";
  fields_border_vel = (char *) "radius rmass angmom";
  fields_exchange = (char *) "radius rmass angmom";
  fields_restart = (char *) "radius rmass angmom";
  fields_create = (char *) "radius rmass angmom tri";
  fields_data_atom = (char *) "id type body rmass x";
  fields_data_vel = (char *) "angmom";

  setup_fields();
}

/* ---------------------------------------------------------------------- */

AtomVecBody::~AtomVecBody()
{
  int nall = nlocal_bonus + nghost_bonus;
  for (int i = 0; i < nall; i++) {
    icp->put(bonus[i].iindex);
    dcp->put(bonus[i].dindex);
  }
  memory->sfree(bonus);

  delete bptr;
}

/* ----------------------------------------------------------------------
   process additional args
   instantiate Body class
   set size_forward and size_border to max sizes
------------------------------------------------------------------------- */

void AtomVecBody::process_args(int narg, char **arg)
{
  // suppress unused parameter warning dependent on style_body.h
  (void)(arg);

  if (narg < 1) error->all(FLERR,"Invalid atom_style body command");

  if (0) bptr = NULL;

#define BODY_CLASS
#define BodyStyle(key,Class) \
  else if (strcmp(arg[0],#key) == 0) bptr = new Class(lmp,narg,arg);
#include "style_body.h"
#undef BodyStyle
#undef BODY_CLASS

  else error->all(FLERR,utils::
                  check_packages_for_style("body",arg[0],lmp).c_str());

  bptr->avec = this;
  icp = bptr->icp;
  dcp = bptr->dcp;

  // max size of forward/border comm
  // bptr values = max number of additional ivalues/dvalues from Body class

  size_forward_bonus += bptr->size_forward;
  size_border_bonus += bptr->size_border;
}

/* ----------------------------------------------------------------------
   grow bonus data structure
------------------------------------------------------------------------- */

void AtomVecBody::grow_bonus()
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

void AtomVecBody::copy_bonus(int i, int j, int delflag)
{
  int *body = atom->body;

  // if deleting atom J via delflag and J has bonus data, then delete it

  if (delflag && body[j] >= 0) {
    int k = body[j];
    icp->put(bonus[k].iindex);
    dcp->put(bonus[k].dindex);
    copy_bonus_all(nlocal_bonus-1,k);
    nlocal_bonus--;
  }

  // if atom I has bonus data, reset I's bonus.ilocal to loc J
  // do NOT do this if self-copy (I=J) since I's bonus data is already deleted

  if (body[i] >= 0 && i != j) bonus[body[i]].ilocal = j;
  body[j] = body[i];
}

/* ----------------------------------------------------------------------
   copy bonus data from I to J, effectively deleting the J entry
   also reset body that points to I to now point to J
------------------------------------------------------------------------- */

void AtomVecBody::copy_bonus_all(int i, int j)
{
  atom->body[bonus[i].ilocal] = j;
  memcpy(&bonus[j],&bonus[i],sizeof(Bonus));
}

/* ----------------------------------------------------------------------
   clear ghost info in bonus data
   called before ghosts are recommunicated in comm and irregular
------------------------------------------------------------------------- */

void AtomVecBody::clear_bonus()
{
  int nall = nlocal_bonus + nghost_bonus;
  for (int i = nlocal_bonus; i < nall; i++) {
    icp->put(bonus[i].iindex);
    dcp->put(bonus[i].dindex);
  }
  nghost_bonus = 0;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->clear_bonus();
}

/* ---------------------------------------------------------------------- */

int AtomVecBody::pack_comm_bonus(int n, int *list, double *buf)
{
  int i,j,m;
  double *quat;

  int *body = atom->body;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (body[j] >= 0) {
      quat = bonus[body[j]].quat;
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
      m += bptr->pack_comm_body(&bonus[body[j]],&buf[m]);
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecBody::unpack_comm_bonus(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  int *body = atom->body;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (body[i] >= 0) {
      quat = bonus[body[i]].quat;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      m += bptr->unpack_comm_body(&bonus[body[i]],&buf[m]);
    }
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecBody::pack_border_bonus(int n, int *list, double *buf)
{
  int i,j,m;
  double *quat,*inertia;

  int *body = atom->body;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (body[j] < 0) buf[m++] = ubuf(0).d;
    else {
      buf[m++] = ubuf(1).d;
      quat = bonus[body[j]].quat;
      inertia = bonus[body[j]].inertia;
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
      buf[m++] = inertia[0];
      buf[m++] = inertia[1];
      buf[m++] = inertia[2];
      buf[m++] = ubuf(bonus[body[j]].ninteger).d;
      buf[m++] = ubuf(bonus[body[j]].ndouble).d;
      m += bptr->pack_border_body(&bonus[body[j]],&buf[m]);
    }
  } 

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBody::unpack_border_bonus(int n, int first, double *buf)
{
  int i,j,m,last;
  double *quat,*inertia;

  int *body = atom->body;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    body[i] = (int) ubuf(buf[m++]).i;
    if (body[i] == 0) body[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      quat = bonus[j].quat;
      inertia = bonus[j].inertia;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      inertia[0] = buf[m++];
      inertia[1] = buf[m++];
      inertia[2] = buf[m++];
      bonus[j].ninteger = (int) ubuf(buf[m++]).i;
      bonus[j].ndouble = (int) ubuf(buf[m++]).i;
      // corresponding put() calls are in clear_bonus()
      bonus[j].ivalue = icp->get(bonus[j].ninteger,bonus[j].iindex);
      bonus[j].dvalue = dcp->get(bonus[j].ndouble,bonus[j].dindex);
      m += bptr->unpack_border_body(&bonus[j],&buf[m]);
      bonus[j].ilocal = i;
      body[i] = j;
      nghost_bonus++;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecBody::pack_exchange_bonus(int i, double *buf)
{
  int m = 0;

  int *body = atom->body;

  if (body[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = body[i];
    double *quat = bonus[j].quat;
    double *inertia = bonus[j].inertia;
    buf[m++] = quat[0];
    buf[m++] = quat[1];
    buf[m++] = quat[2];
    buf[m++] = quat[3];
    buf[m++] = inertia[0];
    buf[m++] = inertia[1];
    buf[m++] = inertia[2];
    buf[m++] = ubuf(bonus[j].ninteger).d;
    buf[m++] = ubuf(bonus[j].ndouble).d;
    memcpy(&buf[m],bonus[j].ivalue,bonus[j].ninteger*sizeof(int));
    if (intdoubleratio == 1) m += bonus[j].ninteger;
    else m += (bonus[j].ninteger+1)/2;
    memcpy(&buf[m],bonus[j].dvalue,bonus[j].ndouble*sizeof(double));
    m += bonus[j].ndouble;
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBody::unpack_exchange_bonus(int ilocal, double *buf)
{
  int m = 0;

  int *body = atom->body;

  body[ilocal] = (int) ubuf(buf[m++]).i;
  if (body[ilocal] == 0) body[ilocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *quat = bonus[nlocal_bonus].quat;
    double *inertia = bonus[nlocal_bonus].inertia;
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    inertia[0] = buf[m++];
    inertia[1] = buf[m++];
    inertia[2] = buf[m++];
    bonus[nlocal_bonus].ninteger = (int) ubuf(buf[m++]).i;
    bonus[nlocal_bonus].ndouble = (int) ubuf(buf[m++]).i;
    // corresponding put() calls are in copy()
    bonus[nlocal_bonus].ivalue = icp->get(bonus[nlocal_bonus].ninteger,
                                          bonus[nlocal_bonus].iindex);
    bonus[nlocal_bonus].dvalue = dcp->get(bonus[nlocal_bonus].ndouble,
                                          bonus[nlocal_bonus].dindex);
    memcpy(bonus[nlocal_bonus].ivalue,&buf[m],
           bonus[nlocal_bonus].ninteger*sizeof(int));
    if (intdoubleratio == 1) m += bonus[nlocal_bonus].ninteger;
    else m += (bonus[nlocal_bonus].ninteger+1)/2;
    memcpy(bonus[nlocal_bonus].dvalue,&buf[m],
           bonus[nlocal_bonus].ndouble*sizeof(double));
    m += bonus[nlocal_bonus].ndouble;

    bonus[nlocal_bonus].ilocal = ilocal;
    body[ilocal] = nlocal_bonus++;
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including bonus data
------------------------------------------------------------------------- */

int AtomVecBody::size_restart_bonus()
{
  int i;

  int *body = atom->body;

  int n = 0;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) {
    if (body[i] >= 0) {
      n += size_restart_bonus_one;
      if (intdoubleratio == 1) n += bonus[body[i]].ninteger;
      else n += (bonus[body[i]].ninteger+1)/2;
      n += bonus[body[i]].ndouble;
    }
    n++;
  }

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecBody::pack_restart_bonus(int i, double *buf)
{
  int m = 0;

  int *body = atom->body;

  if (body[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = body[i];
    double *quat = bonus[j].quat;
    double *inertia = bonus[j].inertia;
    buf[m++] = quat[0];
    buf[m++] = quat[1];
    buf[m++] = quat[2];
    buf[m++] = quat[3];
    buf[m++] = inertia[0];
    buf[m++] = inertia[1];
    buf[m++] = inertia[2];
    buf[m++] = ubuf(bonus[j].ninteger).d;
    buf[m++] = ubuf(bonus[j].ndouble).d;
    memcpy(&buf[m],bonus[j].ivalue,bonus[j].ninteger*sizeof(int));
    if (intdoubleratio == 1) m += bonus[j].ninteger;
    else m += (bonus[j].ninteger+1)/2;
    memcpy(&buf[m],bonus[j].dvalue,bonus[j].ndouble*sizeof(double));
    m += bonus[j].ndouble;
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including bonus data
------------------------------------------------------------------------- */

int AtomVecBody::unpack_restart_bonus(int ilocal, double *buf)
{
  int m = 0;

  int *body = atom->body;

  body[ilocal] = (int) ubuf(buf[m++]).i;
  if (body[ilocal] == 0) body[ilocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *quat = bonus[nlocal_bonus].quat;
    double *inertia = bonus[nlocal_bonus].inertia;
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    inertia[0] = buf[m++];
    inertia[1] = buf[m++];
    inertia[2] = buf[m++];
    bonus[nlocal_bonus].ninteger = (int) ubuf(buf[m++]).i;
    bonus[nlocal_bonus].ndouble = (int) ubuf(buf[m++]).i;
    bonus[nlocal_bonus].ivalue = icp->get(bonus[nlocal_bonus].ninteger,
                                          bonus[nlocal_bonus].iindex);
    bonus[nlocal_bonus].dvalue = dcp->get(bonus[nlocal_bonus].ndouble,
                                          bonus[nlocal_bonus].dindex);
    memcpy(bonus[nlocal_bonus].ivalue,&buf[m],
           bonus[nlocal_bonus].ninteger*sizeof(int));
    if (intdoubleratio == 1) m += bonus[nlocal_bonus].ninteger;
    else m += (bonus[nlocal_bonus].ninteger+1)/2;
    memcpy(bonus[nlocal_bonus].dvalue,&buf[m],
           bonus[nlocal_bonus].ndouble*sizeof(double));
    m += bonus[nlocal_bonus].ndouble;
    bonus[nlocal_bonus].ilocal = ilocal;
    body[ilocal] = nlocal_bonus++;
  }

  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecBody::create_atom_post(int ilocal)
{
  atom->radius[ilocal] = 0.5;
  atom->rmass[ilocal] = 1.0;
  atom->body[ilocal] = -1;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecBody::data_atom_post(int ilocal)
{
  body_flag = atom->body[ilocal];
  if (body_flag == 0) body_flag = -1;
  else if (body_flag == 1) body_flag = 0;
  else error->one(FLERR,"Invalid body flag in Atoms section of data file");
  atom->body[ilocal] = body_flag;

  if (atom->rmass[ilocal] <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  atom->radius[ilocal] = 0.5;
  atom->angmom[ilocal][0] = 0.0;
  atom->angmom[ilocal][1] = 0.0;
  atom->angmom[ilocal][2] = 0.0;
}

/* ----------------------------------------------------------------------
   unpack one body from Bodies section of data file
------------------------------------------------------------------------- */

void AtomVecBody::data_body(int m, int ninteger, int ndouble,
                            int *ivalues, double *dvalues)
{
  if (atom->body[m]) 
    error->one(FLERR,"Assigning body parameters to non-body atom");
  if (nlocal_bonus == nmax_bonus) grow_bonus();
  bonus[nlocal_bonus].ilocal = m;
  bptr->data_body(nlocal_bonus,ninteger,ndouble,ivalues,dvalues);
  atom->body[m] = nlocal_bonus++;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecBody::memory_usage_bonus()
{
  bigint bytes = 0;
  bytes += nmax_bonus*sizeof(Bonus);
  bytes += icp->size + dcp->size;

  int nall = nlocal_bonus + nghost_bonus;
  for (int i = 0; i < nall; i++) {
    bytes += bonus[i].ninteger * sizeof(int);
    bytes += bonus[i].ndouble * sizeof(double);
  }

  return bytes;
}

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_data() to pack
------------------------------------------------------------------------- */

void AtomVecBody::pack_data_pre(int ilocal)
{ 
  body_flag = atom->body[ilocal];

  if (body_flag < 0) atom->body[ilocal] = 0;
  else atom->body[ilocal] = 1;
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_data()
------------------------------------------------------------------------- */

void AtomVecBody::pack_data_post(int ilocal)
{ 
  atom->body[ilocal] = body_flag;
}

/* ----------------------------------------------------------------------
   body computes its size based on ivalues/dvalues and returns it
------------------------------------------------------------------------- */

double AtomVecBody::radius_body(int ninteger, int ndouble,
                                int *ivalues, double *dvalues)
{
  return bptr->radius_body(ninteger,ndouble,ivalues,dvalues);
}

/* ----------------------------------------------------------------------
   reset quat orientation for atom M to quat_external
   called by Atom:add_molecule_atom()
------------------------------------------------------------------------- */

void AtomVecBody::set_quat(int m, double *quat_external)
{
  if (atom->body[m] < 0) error->one(FLERR,"Assigning quat to non-body atom");
  double *quat = bonus[atom->body[m]].quat;
  quat[0] = quat_external[0]; quat[1] = quat_external[1];
  quat[2] = quat_external[2]; quat[3] = quat_external[3];
}

/* ----------------------------------------------------------------------
   debug method for sanity checking of own/bonus data pointers
------------------------------------------------------------------------- */

/*
void AtomVecBody::check(int flag)
{
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->body[i] >= 0 && atom->body[i] >= nlocal_bonus) {
      printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
      errorx->one(FLERR,"BAD AAA");
    }
  }
  for (int i = atom->nlocal; i < atom->nlocal+atom->nghost; i++) {
    if (atom->body[i] >= 0 &&
        (atom->body[i] < nlocal_bonus ||
         atom->body[i] >= nlocal_bonus+nghost_bonus)) {
      printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
      errorx->one(FLERR,"BAD BBB");
    }
  }
  for (int i = 0; i < nlocal_bonus; i++) {
    if (bonus[i].ilocal < 0 || bonus[i].ilocal >= atom->nlocal) {
      printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
      errorx->one(FLERR,"BAD CCC");
    }
  }
  for (int i = 0; i < nlocal_bonus; i++) {
    if (atom->body[bonus[i].ilocal] != i) {
      printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
      errorx->one(FLERR,"BAD DDD");
    }
  }
  for (int i = nlocal_bonus; i < nlocal_bonus+nghost_bonus; i++) {
    if (bonus[i].ilocal < atom->nlocal ||
        bonus[i].ilocal >= atom->nlocal+atom->nghost) {
      printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
      errorx->one(FLERR,"BAD EEE");
    }
  }
  for (int i = nlocal_bonus; i < nlocal_bonus+nghost_bonus; i++) {
    if (atom->body[bonus[i].ilocal] != i) {
      printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
      errorx->one(FLERR,"BAD FFF");
    }
  }
}
*/
