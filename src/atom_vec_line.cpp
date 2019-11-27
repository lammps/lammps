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

#include "atom_vec_line.h"
#include <cmath>
#include <cstring>
#include "atom.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define EPSILON 0.001

/* ---------------------------------------------------------------------- */

AtomVecLine::AtomVecLine(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;
  bonus_flag = 1;

  size_forward_bonus = 1;
  size_border_bonus = 3;
  size_restart_bonus_one = 2;
  size_data_bonus = 5;

  atom->line_flag = 1;
  atom->molecule_flag = atom->rmass_flag = 1;
  atom->radius_flag = atom->omega_flag = atom->torque_flag = 1;
  atom->sphere_flag = 1;

  nlocal_bonus = nghost_bonus = nmax_bonus = 0;
  bonus = NULL;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in the string does not matter
  //   except fields_data_atom and fields_data_vel which must match data file

  fields_grow = (char *) "molecule radius rmass omega torque line";
  fields_copy = (char *) "molecule radius rmass omega";
  fields_comm = NULL;
  fields_comm_vel = (char *) "omega";
  fields_reverse = (char *) "torque";
  fields_border = (char *) "molecule radius rmass";
  fields_border_vel = (char *) "molecule radius rmass omega";
  fields_exchange = (char *) "molecule radius rmass omega";
  fields_restart = (char *) "molecule radius rmass omega";
  fields_create = (char *) "molecule radius rmass omega line";
  fields_data_atom = (char *) "id molecule type line rmass x";
  fields_data_vel = (char *) "omega";

  setup_fields();
}

/* ---------------------------------------------------------------------- */

AtomVecLine::~AtomVecLine()
{
  memory->sfree(bonus);
}

/* ---------------------------------------------------------------------- */

void AtomVecLine::init()
{
  AtomVec::init();

  if (domain->dimension != 2)
    error->all(FLERR,"Atom_style line can only be used in 2d simulations");
}

/* ----------------------------------------------------------------------
   grow bonus data structure
------------------------------------------------------------------------- */

void AtomVecLine::grow_bonus()
{
  nmax_bonus = grow_nmax_bonus(nmax_bonus);
  if (nmax_bonus < 0)
    error->one(FLERR,"Per-processor system is too big");

  bonus = (Bonus *) memory->srealloc(bonus,nmax_bonus*sizeof(Bonus),
                                     "atom:bonus");
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecLine::copy_bonus(int i, int j, int delflag)
{
  int *line = atom->line;

  // if deleting atom J via delflag and J has bonus data, then delete it

  if (delflag && line[j] >= 0) {
    copy_bonus_all(nlocal_bonus-1,line[j]);
    nlocal_bonus--;
  }

  // if atom I has bonus data, reset I's bonus.ilocal to loc J
  // do NOT do this if self-copy (I=J) since I's bonus data is already deleted

  if (line[i] >= 0 && i != j) bonus[line[i]].ilocal = j;
  line[j] = line[i];
}

/* ----------------------------------------------------------------------
   copy bonus data from I to J, effectively deleting the J entry
   also reset line that points to I to now point to J
------------------------------------------------------------------------- */

void AtomVecLine::copy_bonus_all(int i, int j)
{
  atom->line[bonus[i].ilocal] = j;
  memcpy(&bonus[j],&bonus[i],sizeof(Bonus));
}

/* ----------------------------------------------------------------------
   clear ghost info in bonus data
   called before ghosts are recommunicated in comm and irregular
------------------------------------------------------------------------- */

void AtomVecLine::clear_bonus()
{
  nghost_bonus = 0;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->clear_bonus();
}

/* ---------------------------------------------------------------------- */

int AtomVecLine::pack_comm_bonus(int n, int *list, double *buf)
{
  int i,j,m;

  int *line = atom->line;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (line[j] >= 0) buf[m++] = bonus[line[j]].theta;
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecLine::unpack_comm_bonus(int n, int first, double *buf)
{
  int i,m,last;

  int *line = atom->line;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (line[i] >= 0) bonus[line[i]].theta = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecLine::pack_border_bonus(int n, int *list, double *buf)
{
  int i,j,m;

  int *line = atom->line;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (line[j] < 0) buf[m++] = ubuf(0).d;
    else {
      buf[m++] = ubuf(1).d;
      buf[m++] = bonus[line[j]].length;
      buf[m++] = bonus[line[j]].theta;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecLine::unpack_border_bonus(int n, int first, double *buf)
{
  int i,j,m,last;

  int *line = atom->line;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    line[i] = (int) ubuf(buf[m++]).i;
    if (line[i] == 0) line[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      bonus[j].length = buf[m++];
      bonus[j].theta = buf[m++];
      bonus[j].ilocal = i;
      line[i] = j;
      nghost_bonus++;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecLine::pack_exchange_bonus(int i, double *buf)
{
  int m = 0;

  int *line = atom->line;

  if (line[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = line[i];
    buf[m++] = bonus[j].length;
    buf[m++] = bonus[j].theta;
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecLine::unpack_exchange_bonus(int ilocal, double *buf)
{
  int m = 0;

  int *line = atom->line;

  line[ilocal] = (int) ubuf(buf[m++]).i;
  if (line[ilocal] == 0) line[ilocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    bonus[nlocal_bonus].length = buf[m++];
    bonus[nlocal_bonus].theta = buf[m++];
    bonus[nlocal_bonus].ilocal = ilocal;
    line[ilocal] = nlocal_bonus++;
  }

  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecLine::size_restart_bonus()
{
  int i;

  int *line = atom->line;

  int n = 0;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) {
    if (line[i] >= 0) n += size_restart_bonus_one;
    n++;
  }

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecLine::pack_restart_bonus(int i, double *buf)
{
  int m = 0;

  int *line = atom->line;

  if (line[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = line[i];
    buf[m++] = bonus[j].length;
    buf[m++] = bonus[j].theta;
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including bonus data
------------------------------------------------------------------------- */

int AtomVecLine::unpack_restart_bonus(int ilocal, double *buf)
{
  int m = 0;

  int *line = atom->line;

  line[ilocal] = (int) ubuf(buf[m++]).i;
  if (line[ilocal] == 0) line[ilocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    bonus[nlocal_bonus].length = buf[m++];
    bonus[nlocal_bonus].theta = buf[m++];
    bonus[nlocal_bonus].ilocal = ilocal;
    line[ilocal] = nlocal_bonus++;
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack one line from Lines section of data file
------------------------------------------------------------------------- */

void AtomVecLine::data_atom_bonus(int m, char **values)
{
  int *line = atom->line;

  if (line[m]) error->one(FLERR,"Assigning line parameters to non-line atom");

  if (nlocal_bonus == nmax_bonus) grow_bonus();

  double x1 = utils::numeric(FLERR,values[0],true,lmp);
  double y1 = utils::numeric(FLERR,values[1],true,lmp);
  double x2 = utils::numeric(FLERR,values[2],true,lmp);
  double y2 = utils::numeric(FLERR,values[3],true,lmp);
  double dx = x2 - x1;
  double dy = y2 - y1;
  double length = sqrt(dx*dx + dy*dy);

  bonus[nlocal_bonus].length = length;
  if (dy >= 0.0) bonus[nlocal_bonus].theta = acos(dx/length);
  else bonus[nlocal_bonus].theta = -acos(dx/length);

  double xc = 0.5*(x1+x2);
  double yc = 0.5*(y1+y2);
  dx = xc - x[m][0];
  dy = yc - x[m][1];
  double delta = sqrt(dx*dx + dy*dy);

  if (delta/length > EPSILON)
    error->one(FLERR,"Inconsistent line segment in data file");

  x[m][0] = xc;
  x[m][1] = yc;

  // reset line radius and mass
  // rmass currently holds density

  atom->radius[m] = 0.5 * length;
  atom->rmass[m] *= length;

  bonus[nlocal_bonus].ilocal = m;
  line[m] = nlocal_bonus++;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated bonus memory
------------------------------------------------------------------------- */

bigint AtomVecLine::memory_usage_bonus()
{
  bigint bytes = 0;
  bytes += nmax_bonus*sizeof(Bonus);
  return bytes;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecLine::create_atom_post(int ilocal)
{
  double radius = 0.5;
  atom->radius[ilocal] = radius;
  atom->rmass[ilocal] = 4.0*MY_PI/3.0 * radius*radius*radius;
  atom->line[ilocal] = -1;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecLine::data_atom_post(int ilocal)
{
  line_flag = atom->line[ilocal];
  if (line_flag == 0) line_flag = -1;
  else if (line_flag == 1) line_flag = 0;
  else error->one(FLERR,"Invalid line flag in Atoms section of data file");
  atom->line[ilocal] = line_flag;

  if (atom->rmass[ilocal] <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  if (line_flag < 0) {
    double radius = 0.5;
    atom->radius[ilocal] = radius;
    atom->rmass[ilocal] *= 4.0*MY_PI/3.0 * radius*radius*radius;
  } else atom->radius[ilocal] = 0.0;

  atom->omega[ilocal][0] = 0.0;
  atom->omega[ilocal][1] = 0.0;
  atom->omega[ilocal][2] = 0.0;
}

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_data() to pack
------------------------------------------------------------------------- */

void AtomVecLine::pack_data_pre(int ilocal)
{ 
  line_flag = atom->line[ilocal];
  rmass = atom->rmass[ilocal];

  if (line_flag < 0) atom->line[ilocal] = 0;
  else atom->line[ilocal] = 1;

  if (line_flag < 0) {
    double radius = atom->radius[ilocal];
    atom->rmass[ilocal] /= 4.0*MY_PI/3.0 * radius*radius*radius;
  } else atom->rmass[ilocal] /= bonus[line_flag].length;
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_data()
------------------------------------------------------------------------- */

void AtomVecLine::pack_data_post(int ilocal)
{ 
  atom->line[ilocal] = line_flag;
  atom->rmass[ilocal] = rmass;
}

/* ----------------------------------------------------------------------
   set length value in bonus data for particle I
   oriented along x axis
   this may create or delete entry in bonus data
------------------------------------------------------------------------- */

void AtomVecLine::set_length(int i, double value)
{
  int *line = atom->line;

  if (line[i] < 0) {
    if (value == 0.0) return;
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    bonus[nlocal_bonus].length = value;
    bonus[nlocal_bonus].theta = 0.0;
    bonus[nlocal_bonus].ilocal = i;
    line[i] = nlocal_bonus++;
  } else if (value == 0.0) {
    copy_bonus_all(nlocal_bonus-1,line[i]);
    nlocal_bonus--;
    line[i] = -1;
  } else bonus[line[i]].length = value;

  // also set radius = half of length
  // unless value = 0.0, then set diameter = 1.0

  atom->radius[i] = 0.5 * value;
  if (value == 0.0) atom->radius[i] = 0.5;
}

/* ----------------------------------------------------------------------
   check consistency of internal Bonus data structure
   n = # of atoms in regular structure to check against
------------------------------------------------------------------------- */

/*
void AtomVecLine::consistency_check(int n, char *str)
{
  int iflag = 0;
  int count = 0;
  for (int i = 0; i < n; i++) {

    if (line[i] >= 0) {
      count++;
      if (line[i] >= nlocal_bonus) iflag++;
      if (bonus[line[i]].ilocal != i) iflag++;
      //if (comm->me == 1 && update->ntimestep == 873)
      //        printf("CCHK %s: %d %d: %d %d: %d %d\n",
      //       str,i,n,line[i],nlocal_bonus,bonus[line[i]].ilocal,iflag);
    }
  }

  if (iflag) {
    printf("BAD vecline ptrs: %s: %d %d: %d\n",str,comm->me,
           update->ntimestep,iflag);
    MPI_Abort(world,1);
  }

  if (count != nlocal_bonus) {
    char msg[128];
    printf("BAD vecline count: %s: %d %d: %d %d\n",
           str,comm->me,update->ntimestep,count,nlocal_bonus);
    MPI_Abort(world,1);
  }
}
*/
