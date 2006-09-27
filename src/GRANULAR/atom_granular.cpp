/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "atom_granular.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"

/* ---------------------------------------------------------------------- */

AtomGranular::AtomGranular(int narg, char **arg) : Atom(narg, arg) {}

/* ---------------------------------------------------------------------- */

void AtomGranular::copy(int i, int j)
{
  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];
  phix[j][0] = phix[i][0];
  phix[j][1] = phix[i][1];
  phix[j][2] = phix[i][2];
  phiv[j][0] = phiv[i][0];
  phiv[j][1] = phiv[i][1];
  phiv[j][2] = phiv[i][2];
  radius[j] = radius[i];
  density[j] = density[i];
  rmass[j] = rmass[i];

  if (nextra_grow)
    for (int iextra = 0; iextra < nextra_grow; iextra++) 
      modify->fix[extra_grow[iextra]]->copy_arrays(i,j);
}

/* ---------------------------------------------------------------------- */

void AtomGranular::pack_comm(int n, int *list, double *buf, int *pbc_flags)
{
  int i,j,m;

  m = 0;
  if (pbc_flags[0] == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = phiv[j][0];
      buf[m++] = phiv[j][1];
      buf[m++] = phiv[j][2];
    }
  } else {
    double xprd = domain->xprd;
    double yprd = domain->yprd;
    double zprd = domain->zprd;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + pbc_flags[1]*xprd;
      buf[m++] = x[j][1] + pbc_flags[2]*yprd;
      buf[m++] = x[j][2] + pbc_flags[3]*zprd;
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = phiv[j][0];
      buf[m++] = phiv[j][1];
      buf[m++] = phiv[j][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void AtomGranular::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    phiv[i][0] = buf[m++];
    phiv[i][1] = buf[m++];
    phiv[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void AtomGranular::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
    buf[m++] = phia[i][0];
    buf[m++] = phia[i][1];
    buf[m++] = phia[i][2];
  }
}

/* ---------------------------------------------------------------------- */

void AtomGranular::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
    phia[j][0] += buf[m++];
    phia[j][1] += buf[m++];
    phia[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void AtomGranular::pack_border(int n, int *list, double *buf, int *pbc_flags)
{
  int i,j,m;

  m = 0;
  if (pbc_flags[0] == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = phiv[j][0];
      buf[m++] = phiv[j][1];
      buf[m++] = phiv[j][2];
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
    }
  } else {
    double xprd = domain->xprd;
    double yprd = domain->yprd;
    double zprd = domain->zprd;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + pbc_flags[1]*xprd;
      buf[m++] = x[j][1] + pbc_flags[2]*yprd;
      buf[m++] = x[j][2] + pbc_flags[3]*zprd;
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = phiv[j][0];
      buf[m++] = phiv[j][1];
      buf[m++] = phiv[j][2];
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
    }
  }
}

/* ---------------------------------------------------------------------- */

void AtomGranular::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = static_cast<int> (buf[m++]);
    type[i] = static_cast<int> (buf[m++]);
    mask[i] = static_cast<int> (buf[m++]);
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    phiv[i][0] = buf[m++];
    phiv[i][1] = buf[m++];
    phiv[i][2] = buf[m++];
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
  }
}

/* ----------------------------------------------------------------------
   pack all atom quantities for shipping to another proc
   xyz must be 1st 3 values, so that comm::exchange can test on them
------------------------------------------------------------------------- */

int AtomGranular::pack_exchange(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = phix[i][0];
  buf[m++] = phix[i][1];
  buf[m++] = phix[i][2];
  buf[m++] = phiv[i][0];
  buf[m++] = phiv[i][1];
  buf[m++] = phiv[i][2];
  buf[m++] = radius[i];
  buf[m++] = density[i];
  buf[m++] = rmass[i];
  buf[m++] = tag[i];
  buf[m++] = type[i];
  buf[m++] = mask[i];
  buf[m++] = image[i];

  if (nextra_grow)
    for (int iextra = 0; iextra < nextra_grow; iextra++) 
      m += modify->fix[extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomGranular::unpack_exchange(double *buf)
{
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  phix[nlocal][0] = buf[m++];
  phix[nlocal][1] = buf[m++];
  phix[nlocal][2] = buf[m++];
  phiv[nlocal][0] = buf[m++];
  phiv[nlocal][1] = buf[m++];
  phiv[nlocal][2] = buf[m++];
  radius[nlocal] = buf[m++];
  density[nlocal] = buf[m++];
  rmass[nlocal] = buf[m++];
  tag[nlocal] = static_cast<int> (buf[m++]);
  type[nlocal] = static_cast<int> (buf[m++]);
  mask[nlocal] = static_cast<int> (buf[m++]);
  image[nlocal] = static_cast<int> (buf[m++]);

  if (nextra_grow)
    for (int iextra = 0; iextra < nextra_grow; iextra++) 
      m += modify->fix[extra_grow[iextra]]->unpack_exchange(nlocal,&buf[m]);

  nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   unpack vels & phis specific to granular data file
------------------------------------------------------------------------- */

void AtomGranular::unpack_vels(int n, char *buf)
{
  int m,tagtmp;
  double vxtmp,vytmp,vztmp,phivxtmp,phivytmp,phivztmp;
  char *next;

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    sscanf(buf,"%d %lg %lg %lg %lg %lg %lg",
	   &tagtmp,&vxtmp,&vytmp,&vztmp,
	   &phivxtmp,&phivytmp,&phivztmp);
    if ((m = map(tagtmp)) >= 0) {
      v[m][0] = vxtmp;
      v[m][1] = vytmp;
      v[m][2] = vztmp;
      phiv[m][0] = phivxtmp;
      phiv[m][1] = phivytmp;
      phiv[m][2] = phivztmp;
    }
    buf = next + 1;
  }
}
