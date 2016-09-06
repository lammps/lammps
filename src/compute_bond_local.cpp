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

#include <math.h>
#include <string.h>
#include "compute_bond_local.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "update.h"
#include "domain.h"
#include "force.h"
#include "bond.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000
#define SMALL 1.0e-15

enum{DIST,ENGPOT,FORCE,VELVIB,VELROT,ENGTRANS,ENGVIB,ENGROT};

/* ---------------------------------------------------------------------- */

ComputeBondLocal::ComputeBondLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  bstyle(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal compute bond/local command");

  if (atom->avec->bonds_allow == 0)
    error->all(FLERR,"Compute bond/local used when bonds are not allowed");

  local_flag = 1;
  nvalues = narg - 3;
  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  bstyle = new int[nvalues];

  nvalues = 0;
  for (int iarg = 3; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"dist") == 0) bstyle[nvalues++] = DIST;
    else if (strcmp(arg[iarg],"engpot") == 0) bstyle[nvalues++] = ENGPOT;
    else if (strcmp(arg[iarg],"force") == 0) bstyle[nvalues++] = FORCE;
    else if (strcmp(arg[iarg],"velvib") == 0) bstyle[nvalues++] = VELVIB;
    else if (strcmp(arg[iarg],"velrot") == 0) bstyle[nvalues++] = VELROT;
    else if (strcmp(arg[iarg],"engtrans") == 0) bstyle[nvalues++] = ENGTRANS;
    else if (strcmp(arg[iarg],"engvib") == 0) bstyle[nvalues++] = ENGVIB;
    else if (strcmp(arg[iarg],"engrot") == 0) bstyle[nvalues++] = ENGROT;
    else error->all(FLERR,"Invalid keyword in compute bond/local command");
  }

  // set singleflag if need to call bond->single()

  singleflag = 0;
  for (int i = 0; i < nvalues; i++)
    if (bstyle[i] != DIST) singleflag = 1;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeBondLocal::~ComputeBondLocal()
{
  memory->destroy(vector);
  memory->destroy(array);
  delete [] bstyle;
}

/* ---------------------------------------------------------------------- */

void ComputeBondLocal::init()
{
  if (force->bond == NULL)
    error->all(FLERR,"No bond style is defined for compute bond/local");

  // do initial memory allocation so that memory_usage() is correct

  ncount = compute_bonds(0);
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
  int i,m,n,nb,atom1,atom2,imol,iatom,btype;
  tagint tagprev;
  double dx,dy,dz,rsq;
  double dvx,dvy,dvz,vvib,vrotsq;
  double vcmx,vcmy,vcmz;
  double masstotal,massreduced;
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
  double engpot,engtrans,engvib,engrot,fbond;

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
        btype = atom->map(onemols[imol]->bond_type[iatom][i]);
        atom2 = atom->map(onemols[imol]->bond_atom[iatom][i]+tagprev);
      }

      if (atom2 < 0 || !(mask[atom2] & groupbit)) continue;
      if (newton_bond == 0 && tag[atom1] > tag[atom2]) continue;
      if (btype == 0) continue;

      if (flag) {
        dx = x[atom1][0] - x[atom2][0];
        dy = x[atom1][1] - x[atom2][1];
        dz = x[atom1][2] - x[atom2][2];
        domain->minimum_image(dx,dy,dz);
        rsq = dx*dx + dy*dy + dz*dz;

        if (singleflag && (btype > 0))
          engpot = bond->single(btype,rsq,atom1,atom2,fbond);
        else engpot = fbond = 0.0;

        dvx = v[atom1][0] - v[atom2][0];
        dvy = v[atom1][1] - v[atom2][1];
        dvz = v[atom1][2] - v[atom2][2];

        if (rmass) {
          masstotal = rmass[atom1]+rmass[atom2];
          vcmx = (rmass[atom1]*v[atom1][0] + rmass[atom2]*v[atom2][0]) / masstotal;
          vcmy = (rmass[atom1]*v[atom1][1] + rmass[atom2]*v[atom2][1]) / masstotal;
          vcmz = (rmass[atom1]*v[atom1][2] + rmass[atom2]*v[atom2][2]) / masstotal;
        }
        else {
          masstotal = mass[type[atom1]]+mass[type[atom2]];
          vcmx = (mass[type[atom1]]*v[atom1][0] + mass[type[atom2]]*v[atom2][0]) / masstotal;
          vcmy = (mass[type[atom1]]*v[atom1][1] + mass[type[atom2]]*v[atom2][1]) / masstotal;
          vcmz = (mass[type[atom1]]*v[atom1][2] + mass[type[atom2]]*v[atom2][2]) / masstotal;
        }
        engtrans=0.5*masstotal*(vcmx*vcmx+vcmy*vcmy+vcmz*vcmz)*force->mvv2e;

        for (int i = 0; i < nvalues; i++) {
          if (bstyle[i] == VELVIB || bstyle[i] == VELROT || bstyle[i] == ENGVIB || bstyle[i] == ENGROT) {
            // compute velocity for each bond by changing basis from x,y,z to that
            // along the bond vector from v'=inv(M)v, where the columns of M are
            // the bond vector and two other vectors that make up an orthonormal
            // basis
            double ione[3][3],inverse[3][3],norm;
            ione[0][0] = dx;
            ione[1][0] = dy;
            ione[2][0] = dz;
            // normalize
            norm = sqrt(ione[0][0]*ione[0][0] + ione[1][0]*ione[1][0] + ione[2][0]*ione[2][0]);
            ione[0][0] /= norm;
            ione[1][0] /= norm;
            ione[2][0] /= norm;
            // get vector that is perpendicular to the bond vector
            if(fabs(dz)>=SMALL)
              {
                ione[0][1] = 0.0;
                ione[1][1] = 1.0;
                ione[2][1] = (-ione[0][0] * ione[0][1] - ione[1][0] * ione[1][1]) / ione[2][0];
              }
              else if(fabs(dx)>=SMALL)
              {
                ione[1][1] = 0.0;
                ione[2][1] = 1.0;
                ione[0][1] = (-ione[1][0] * ione[1][1] - ione[2][0] * ione[2][1]) / ione[0][0];
              }
              else if(fabs(dy)>=SMALL)
              {
                ione[2][1] = 0.0;
                ione[0][1] = 1.0;
                ione[1][1] = (-ione[2][0] * ione[2][1] - ione[0][0] * ione[0][1]) / ione[1][0];
              }
            // normalize
            norm = sqrt(ione[0][1]*ione[0][1] + ione[1][1]*ione[1][1] + ione[2][1]*ione[2][1]);
            ione[0][1] /= norm;
            ione[1][1] /= norm;
            ione[2][1] /= norm;
            // find the last vector from the cross product
            ione[0][2] = ione[1][0] * ione[2][1] - ione[1][1] * ione[2][0];
            ione[1][2] = -(ione[0][0] * ione[2][1] - ione[0][1] * ione[2][0]);
            ione[2][2] = ione[0][0] * ione[1][1] - ione[0][1] * ione[1][0];
            // normalize
            norm = sqrt(ione[0][2]*ione[0][2] + ione[1][2]*ione[1][2] + ione[2][2]*ione[2][2]);
            ione[0][2] /= norm;
            ione[1][2] /= norm;
            ione[2][2] /= norm;
            // compute inverse
            double invdet = ione[0][0]*ione[1][1]*ione[2][2] +
              ione[0][1]*ione[1][2]*ione[2][0] + ione[0][2]*ione[1][0]*ione[2][1] -
              ione[0][0]*ione[1][2]*ione[2][1] - ione[0][1]*ione[1][0]*ione[2][2] -
              ione[2][0]*ione[1][1]*ione[0][2];
            invdet = 1.0/invdet; // determinant should always be 1, so this doesn't really matter
            inverse[0][0] = invdet*(ione[1][1]*ione[2][2] - ione[1][2]*ione[2][1]);
            inverse[0][1] = -invdet*(ione[0][1]*ione[2][2] - ione[0][2]*ione[2][1]);
            inverse[0][2] = invdet*(ione[0][1]*ione[1][2] - ione[0][2]*ione[1][1]);
            inverse[1][0] = -invdet*(ione[1][0]*ione[2][2] - ione[1][2]*ione[2][0]);
            inverse[1][1] = invdet*(ione[0][0]*ione[2][2] - ione[0][2]*ione[2][0]);
            inverse[1][2] = -invdet*(ione[0][0]*ione[1][2] - ione[0][2]*ione[1][0]);
            inverse[2][0] = invdet*(ione[1][0]*ione[2][1] - ione[1][1]*ione[2][0]);
            inverse[2][1] = -invdet*(ione[0][0]*ione[2][1] - ione[0][1]*ione[2][0]);
            inverse[2][2] = invdet*(ione[0][0]*ione[1][1] - ione[0][1]*ione[1][0]);
            vvib = inverse[0][0]*dvx + inverse[0][1]*dvy + inverse[0][2]*dvz;
            vrotsq = (inverse[1][0]*dvx + inverse[1][1]*dvy + inverse[1][2]*dvz) *
                     (inverse[1][0]*dvx + inverse[1][1]*dvy + inverse[1][2]*dvz) +
                     (inverse[2][0]*dvx + inverse[2][1]*dvy + inverse[2][2]*dvz) *
                     (inverse[2][0]*dvx + inverse[2][1]*dvy + inverse[2][2]*dvz);
            if (rmass) massreduced = rmass[atom1]*rmass[atom2]/(rmass[atom1]+rmass[atom2]);
            else massreduced = mass[type[atom1]]*mass[type[atom2]]/(mass[type[atom1]]+mass[type[atom2]]);
            engvib=0.5*massreduced*vvib*vvib*force->mvv2e;
            engrot=0.5*massreduced*vrotsq*force->mvv2e;
            break;
          }
        }

        if (nvalues == 1) ptr = &vector[m];
        else ptr = array[m];

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
          case VELVIB:
            ptr[n] = vvib;
            break;
          case VELROT:
            ptr[n] = sqrt(vrotsq);
            break;
          case ENGTRANS:
            ptr[n] = engtrans;
            break;
          case ENGVIB:
            ptr[n] = engvib;
            break;
          case ENGROT:
            ptr[n] = engrot;
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

void ComputeBondLocal::reallocate(int n)
{
  // grow vector or array and indices array

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->destroy(vector);
    memory->create(vector,nmax,"bond/local:vector");
    vector_local = vector;
  } else {
    memory->destroy(array);
    memory->create(array,nmax,nvalues,"bond/local:array");
    array_local = array;
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
