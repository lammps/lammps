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

#include <cmath>
#include <cstring>
#include "compute_bond_local.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "update.h"
#include "domain.h"
#include "force.h"
#include "bond.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000
#define EPSILON 1.0e-12

enum{DIST,VELVIB,OMEGA,ENGTRANS,ENGVIB,ENGROT,ENGPOT,FORCE,VARIABLE};

/* ---------------------------------------------------------------------- */

ComputeBondLocal::ComputeBondLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  bstyle(NULL), vvar(NULL), dstr(NULL), vstr(NULL), vlocal(NULL), alocal(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal compute bond/local command");

  if (atom->avec->bonds_allow == 0)
    error->all(FLERR,"Compute bond/local used when bonds are not allowed");

  local_flag = 1;
  comm_forward = 3;

  // style args

  nvalues = narg - 3;
  bstyle = new int[nvalues];
  vstr = new char*[nvalues];
  vvar = new int[nvalues];

  nvalues = 0;
  nvar = 0;

  int iarg;
  for (iarg = 3; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"dist") == 0) bstyle[nvalues++] = DIST;
    else if (strcmp(arg[iarg],"engpot") == 0) bstyle[nvalues++] = ENGPOT;
    else if (strcmp(arg[iarg],"force") == 0) bstyle[nvalues++] = FORCE;
    else if (strcmp(arg[iarg],"engvib") == 0) bstyle[nvalues++] = ENGVIB;
    else if (strcmp(arg[iarg],"engrot") == 0) bstyle[nvalues++] = ENGROT;
    else if (strcmp(arg[iarg],"engtrans") == 0) bstyle[nvalues++] = ENGTRANS;
    else if (strcmp(arg[iarg],"omega") == 0) bstyle[nvalues++] = OMEGA;
    else if (strcmp(arg[iarg],"velvib") == 0) bstyle[nvalues++] = VELVIB;
    else if (strncmp(arg[iarg],"v_",2) == 0) {
      bstyle[nvalues++] = VARIABLE;
      int n = strlen(arg[iarg]);
      vstr[nvar] = new char[n];
      strcpy(vstr[nvar],&arg[iarg][2]);
      nvar++;
    } else break;
  }

  // optional args

  setflag = 0;
  dstr = NULL;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"set") == 0) {
      setflag = 1;
      if (iarg+3 > narg) error->all(FLERR,"Illegal compute bond/local command");
      if (strcmp(arg[iarg+1],"dist") == 0) {
        delete [] dstr;
        int n = strlen(arg[iarg+2]) + 1;
        dstr = new char[n];
        strcpy(dstr,arg[iarg+2]);
      } else error->all(FLERR,"Illegal compute bond/local command");
      iarg += 3;
    } else error->all(FLERR,"Illegal compute bond/local command");
  }

  // error check

  if (nvar) {
    if (!setflag)
      error->all(FLERR,"Compute bond/local variable requires a set variable");
    for (int i = 0; i < nvar; i++) {
      vvar[i] = input->variable->find(vstr[i]);
      if (vvar[i] < 0)
        error->all(FLERR,"Variable name for copute bond/local does not exist");
      if (!input->variable->equalstyle(vvar[i]))
        error->all(FLERR,"Variable for compute bond/local is invalid style");
    }

    if (dstr) {
      dvar = input->variable->find(dstr);
      if (dvar < 0)
        error->all(FLERR,"Variable name for compute bond/local does not exist");
      if (!input->variable->internalstyle(dvar))
        error->all(FLERR,"Variable for compute bond/local is invalid style");
    }
  } else if (setflag)
    error->all(FLERR,"Compute bond/local set with no variable");


  // set singleflag if need to call bond->single()
  // set velflag if compute any quantities based on velocities

  singleflag = 0;
  velflag = 0;
  for (int i = 0; i < nvalues; i++) {
    if (bstyle[i] == ENGPOT || bstyle[i] == FORCE) singleflag = 1;
    if (bstyle[i] == VELVIB || bstyle[i] == OMEGA || bstyle[i] == ENGTRANS ||
        bstyle[i] == ENGVIB || bstyle[i] == ENGROT) velflag = 1;
  }

  // initialize output

  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  nmax = 0;
  vlocal = NULL;
  alocal = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeBondLocal::~ComputeBondLocal()
{
  delete [] bstyle;
  for (int i = 0; i < nvar; i++) delete [] vstr[i];
  delete [] vstr;
  delete [] vvar;

  delete [] dstr;

  memory->destroy(vlocal);
  memory->destroy(alocal);
}

/* ---------------------------------------------------------------------- */

void ComputeBondLocal::init()
{
  if (force->bond == NULL)
    error->all(FLERR,"No bond style is defined for compute bond/local");

  if (nvar) {
    for (int i = 0; i < nvar; i++) {
      vvar[i] = input->variable->find(vstr[i]);
      if (vvar[i] < 0)
        error->all(FLERR,"Variable name for compute bond/local does not exist");
    }

    if (dstr) {
      dvar = input->variable->find(dstr);
      if (dvar < 0)
        error->all(FLERR,"Variable name for compute bond/local does not exist");
    }
  }

  // set ghostvelflag if need to acquire ghost atom velocities

  if (velflag && !comm->ghost_velocity) ghostvelflag = 1;
  else ghostvelflag = 0;

  // do initial memory allocation so that memory_usage() is correct

  initflag = 1;
  ncount = compute_bonds(0);
  initflag = 0;

  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
}

/* ---------------------------------------------------------------------- */

void ComputeBondLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and compute bond info

  ncount = compute_bonds(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
  ncount = compute_bonds(1);
}

/* ----------------------------------------------------------------------
   count bonds and compute bond info on this proc
   only count bond once if newton_bond is off
   all atoms in interaction must be in group
   all atoms in interaction must be known to proc
   if bond is deleted (type = 0), do not count
   if bond is turned off (type < 0), still count
   if flag is set, compute requested info about bond
   if bond is turned off (type < 0), energy = 0.0
------------------------------------------------------------------------- */

int ComputeBondLocal::compute_bonds(int flag)
{
  int i,m,n,nb,atom1,atom2,imol,iatom,btype,ivar;
  tagint tagprev;
  double dx,dy,dz,rsq;
  double mass1,mass2,masstotal,invmasstotal;
  double xcm[3],vcm[3];
  double delr1[3],delr2[3],delv1[3],delv2[3];
  double r12[3],vpar1,vpar2;
  double vvib,vrotsq;
  double inertia,omegasq;
  double mvv2e;
  double engpot,engtrans,engvib,engrot,fbond;
  double *ptr;

  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  tagint *tag = atom->tag;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *mask = atom->mask;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;

  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  int molecular = atom->molecular;

  Bond *bond = force->bond;

  // communicate ghost velocities if needed

  if (ghostvelflag && !initflag) comm->forward_comm_compute(this);

  // loop over all atoms and their bonds

  m = n = 0;
  for (atom1 = 0; atom1 < nlocal; atom1++) {
    if (!(mask[atom1] & groupbit)) continue;

    if (molecular == 1) nb = num_bond[atom1];
    else {
      if (molindex[atom1] < 0) continue;
      imol = molindex[atom1];
      iatom = molatom[atom1];
      nb = onemols[imol]->num_bond[iatom];
    }

    for (i = 0; i < nb; i++) {
      if (molecular == 1) {
        btype = bond_type[atom1][i];
        atom2 = atom->map(bond_atom[atom1][i]);
      } else {
        tagprev = tag[atom1] - iatom - 1;
        btype = onemols[imol]->bond_type[iatom][i];
        atom2 = atom->map(onemols[imol]->bond_atom[iatom][i]+tagprev);
      }

      if (atom2 < 0 || !(mask[atom2] & groupbit)) continue;
      if (newton_bond == 0 && tag[atom1] > tag[atom2]) continue;
      if (btype == 0) continue;

      if (!flag) {
        m++;
        continue;
      }

      dx = x[atom1][0] - x[atom2][0];
      dy = x[atom1][1] - x[atom2][1];
      dz = x[atom1][2] - x[atom2][2];
      domain->minimum_image(dx,dy,dz);
      rsq = dx*dx + dy*dy + dz*dz;

      if (btype == 0) {
        engpot = fbond = 0.0;
        engvib = engrot = engtrans = omegasq = vvib = 0.0;
      } else {

        if (singleflag) engpot = bond->single(btype,rsq,atom1,atom2,fbond);

        if (velflag) {
          if (rmass) {
            mass1 = rmass[atom1];
            mass2 = rmass[atom2];
          }
          else {
            mass1 = mass[type[atom1]];
            mass2 = mass[type[atom2]];
          }
          masstotal = mass1+mass2;
          invmasstotal = 1.0 / (masstotal);
          xcm[0] = (mass1*x[atom1][0] + mass2*x[atom2][0]) * invmasstotal;
          xcm[1] = (mass1*x[atom1][1] + mass2*x[atom2][1]) * invmasstotal;
          xcm[2] = (mass1*x[atom1][2] + mass2*x[atom2][2]) * invmasstotal;
          vcm[0] = (mass1*v[atom1][0] + mass2*v[atom2][0]) * invmasstotal;
          vcm[1] = (mass1*v[atom1][1] + mass2*v[atom2][1]) * invmasstotal;
          vcm[2] = (mass1*v[atom1][2] + mass2*v[atom2][2]) * invmasstotal;

          engtrans = 0.5 * masstotal * MathExtra::lensq3(vcm);

          // r12 = unit bond vector from atom1 to atom2

          MathExtra::sub3(x[atom2],x[atom1],r12);
          MathExtra::norm3(r12);

          // delr = vector from COM to each atom
          // delv = velocity of each atom relative to COM

          MathExtra::sub3(x[atom1],xcm,delr1);
          MathExtra::sub3(x[atom2],xcm,delr2);
          MathExtra::sub3(v[atom1],vcm,delv1);
          MathExtra::sub3(v[atom2],vcm,delv2);

          // vpar = component of delv parallel to bond vector

          vpar1 = MathExtra::dot3(delv1,r12);
          vpar2 = MathExtra::dot3(delv2,r12);
          engvib = 0.5 * (mass1*vpar1*vpar1 + mass2*vpar2*vpar2);

          // vvib = relative velocity of 2 atoms along bond direction
          // vvib < 0 for 2 atoms moving towards each other
          // vvib > 0 for 2 atoms moving apart

          vvib = vpar2 - vpar1;

          // vrotsq = tangential speed squared of atom1 only
          // omegasq = omega squared, and is the same for atom1 and atom2

          inertia = mass1*MathExtra::lensq3(delr1) +
            mass2*MathExtra::lensq3(delr2);
          vrotsq = MathExtra::lensq3(delv1) - vpar1*vpar1;
          omegasq = vrotsq / MathExtra::lensq3(delr1);

          engrot = 0.5 * inertia * omegasq;

          // sanity check: engtotal = engtrans + engvib + engrot

          //engtot = 0.5 * (mass1*MathExtra::lensq3(v[atom1]) +
          //                mass2*MathExtra::lensq3(v[atom2]));
          //if (fabs(engtot-engtrans-engvib-engrot) > EPSILON)
          //  error->one(FLERR,"Sanity check on 3 energy components failed");

          // scale energies by units

          mvv2e = force->mvv2e;
          engtrans *= mvv2e;
          engvib *= mvv2e;
          engrot *= mvv2e;
        }

        if (nvalues == 1) ptr = &vlocal[m];
        else ptr = alocal[m];

        if (nvar) {
          ivar = 0;
          if (dstr) input->variable->internal_set(dvar,sqrt(rsq));
        }

        for (n = 0; n < nvalues; n++) {
          switch (bstyle[n]) {
          case DIST:
            ptr[n] = sqrt(rsq);
            break;
          case ENGPOT:
            ptr[n] = engpot;
            break;
          case FORCE:
            ptr[n] = sqrt(rsq)*fbond;
            break;
          case ENGVIB:
            ptr[n] = engvib;
            break;
          case ENGROT:
            ptr[n] = engrot;
            break;
          case ENGTRANS:
            ptr[n] = engtrans;
            break;
          case OMEGA:
            ptr[n] = sqrt(omegasq);
            break;
          case VELVIB:
            ptr[n] = vvib;
            break;
          case VARIABLE:
            ptr[n] = input->variable->compute_equal(vvar[ivar]);
            ivar++;
            break;
          }
        }
      }

      m++;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int ComputeBondLocal::pack_forward_comm(int n, int *list, double *buf,
                                        int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  double **v = atom->v;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = v[j][0];
    buf[m++] = v[j][1];
    buf[m++] = v[j][2];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeBondLocal::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  double **v = atom->v;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void ComputeBondLocal::reallocate(int n)
{
  // grow vector_local or array_local

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->destroy(vlocal);
    memory->create(vlocal,nmax,"bond/local:vector_local");
    vector_local = vlocal;
  } else {
    memory->destroy(alocal);
    memory->create(alocal,nmax,nvalues,"bond/local:array_local");
    array_local = alocal;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeBondLocal::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  return bytes;
}
