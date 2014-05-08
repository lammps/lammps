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
   Contributing author: Ilya Valuev (JIHT, Moscow, Russia)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "atom_vec_wavepacket.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "force.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecWavepacket::AtomVecWavepacket(LAMMPS *lmp) : AtomVec(lmp)
{
  comm_x_only = comm_f_only = 0;

  mass_type = 1;
  molecular = 0;

  size_forward = 4; // coords[3]+radius[1]
  size_reverse = 10; // force[3]+erforce[1]+ervelforce[1]+vforce[3]+csforce[2]
  size_border = 10; // coords[3]+tag[1]+type[1]+mask[1]+q[1]+spin[1]+eradius[1]+etag[1]
  size_velocity = 6; // +velocities[3]+ ervel[1]+cs[2]
  size_data_atom = 11; // for input file: 1-tag 2-type 3-q 4-spin 5-eradius 6-etag 7-cs_re 8-cs_im 9-x 10-y 11-z
  size_data_vel = 5; // for input file: vx vy vz ervel <??>
  xcol_data = 9; // starting column for x data

  atom->wavepacket_flag = 1;
  atom->electron_flag = 1; // compatible with eff
  atom->q_flag = atom->spin_flag = atom->eradius_flag =
    atom->ervel_flag = atom->erforce_flag = 1;

  atom->cs_flag = atom->csforce_flag = atom->vforce_flag = atom->ervelforce_flag = atom->etag_flag = 1;
}

/* ----------------------------------------------------------------------
   grow atom-electron arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecWavepacket::grow(int n)
{
  if (n == 0) grow_nmax();
  else nmax = n;
  atom->nmax = nmax;

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");

  q = memory->grow(atom->q,nmax,"atom:q");
  spin = memory->grow(atom->spin,nmax,"atom:spin");
  eradius = memory->grow(atom->eradius,nmax,"atom:eradius");
  ervel = memory->grow(atom->ervel,nmax,"atom:ervel");
  erforce = memory->grow(atom->erforce,nmax*comm->nthreads,"atom:erforce");

  cs = memory->grow(atom->cs,2*nmax,"atom:cs");
  csforce = memory->grow(atom->csforce,2*nmax,"atom:csforce");
  vforce = memory->grow(atom->vforce,3*nmax,"atom:vforce");
  ervelforce = memory->grow(atom->ervelforce,nmax,"atom:ervelforce");
  etag = memory->grow(atom->etag,nmax,"atom:etag");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecWavepacket::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  q = atom->q;
  eradius = atom->eradius; ervel = atom->ervel; erforce = atom->erforce;

  cs = atom->cs;
  csforce = atom->csforce;
  vforce = atom->vforce;
  ervelforce = atom->ervelforce;
  etag = atom->etag;

}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecWavepacket::copy(int i, int j, int delflag)
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

  q[j] = q[i];
  spin[j] = spin[i];
  eradius[j] = eradius[i];
  ervel[j] = ervel[i];

  cs[2*j] = cs[2*i];
  cs[2*j+1] = cs[2*i+1];
  etag[j] = etag[i];


  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */
// this will be used as partial pack for unsplit Hartree packets (v, ervel not regarded as separate variables)

int AtomVecWavepacket::pack_comm(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = eradius[j];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = eradius[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */
// this is a complete pack of all 'position' variables of AWPMD

int AtomVecWavepacket::pack_comm_vel(int n, int *list, double *buf,
                                   int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = eradius[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];

      buf[m++] = ervel[j];
      buf[m++] = cs[2*j];
      buf[m++] = cs[2*j+1];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = eradius[j];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];

        buf[m++] = ervel[j];
        buf[m++] = cs[2*j];
        buf[m++] = cs[2*j+1];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = eradius[j];
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
        buf[m++] = ervel[j];
        buf[m++] = cs[2*j];
        buf[m++] = cs[2*j+1];
      }
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecWavepacket::pack_comm_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = eradius[j];
    buf[m++] = ervel[j];
    buf[m++] = cs[2*j];
    buf[m++] = cs[2*j+1];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecWavepacket::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    eradius[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecWavepacket::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    eradius[i] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];

    ervel[i] = buf[m++];
    cs[2*i] =  buf[m++];
    cs[2*i+1] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecWavepacket::unpack_comm_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++){
    eradius[i] = buf[m++];
    ervel[i] = buf[m++];
    cs[2*i] =  buf[m++];
    cs[2*i+1] = buf[m++];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecWavepacket::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) { //10
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
    buf[m++] = erforce[i];

    buf[m++] = ervelforce[i];
    buf[m++] = vforce[3*i];
    buf[m++] = vforce[3*i+1];
    buf[m++] = vforce[3*i+2];
    buf[m++] = csforce[2*i];
    buf[m++] = csforce[2*i+1];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecWavepacket::pack_reverse_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++){
    buf[m++] = erforce[i];

    buf[m++] = ervelforce[i];
    buf[m++] = vforce[3*i];
    buf[m++] = vforce[3*i+1];
    buf[m++] = vforce[3*i+2];
    buf[m++] = csforce[2*i];
    buf[m++] = csforce[2*i+1];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecWavepacket::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
    erforce[j] += buf[m++];

    ervelforce[j] += buf[m++];
    vforce[3*j] += buf[m++];
    vforce[3*j+1] += buf[m++];
    vforce[3*j+2] += buf[m++];
    csforce[2*j] += buf[m++];
    csforce[2*j+1] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecWavepacket::unpack_reverse_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    erforce[j] += buf[m++];

    ervelforce[j] += buf[m++];
    vforce[3*j] += buf[m++];
    vforce[3*j+1] += buf[m++];
    vforce[3*j+2] += buf[m++];
    csforce[2*j] += buf[m++];
    csforce[2*j+1] += buf[m++];
  }
  return m;
}

/* ---------------------------------------------------------------------- */
// will be used for Hartree unsplit version (the etag is added however)
int AtomVecWavepacket::pack_border(int n, int *list, double *buf,
                                 int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = q[j];
      buf[m++] = ubuf(spin[j]).d;
      buf[m++] = eradius[j];
      buf[m++] = ubuf(etag[j]).d;
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = q[j];
      buf[m++] = ubuf(spin[j]).d;
      buf[m++] = eradius[j];
      buf[m++] = ubuf(etag[j]).d;
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecWavepacket::pack_border_vel(int n, int *list, double *buf,
                                     int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = q[j];
      buf[m++] = ubuf(spin[j]).d;
      buf[m++] = eradius[j];
      buf[m++] = ubuf(etag[j]).d;

      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];

      buf[m++] = ervel[j];
      buf[m++] = cs[2*j];
      buf[m++] = cs[2*j+1];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (domain->triclinic == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = q[j];
        buf[m++] = ubuf(spin[j]).d;
        buf[m++] = eradius[j];
        buf[m++] = ubuf(etag[j]).d;

        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];


        buf[m++] = ervel[j];
        buf[m++] = cs[2*j];
        buf[m++] = cs[2*j+1];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = q[j];
        buf[m++] = ubuf(spin[j]).d;
        buf[m++] = eradius[j];
        buf[m++] = ubuf(etag[j]).d;

        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }

        buf[m++] = ervel[j];
        buf[m++] = cs[2*j];
        buf[m++] = cs[2*j+1];
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecWavepacket::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = q[j];
    buf[m++] = ubuf(spin[j]).d;
    buf[m++] = eradius[j];

    buf[m++] = ubuf(etag[j]).d;
    buf[m++] = ervel[j];
    buf[m++] = cs[2*j];
    buf[m++] = cs[2*j+1];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecWavepacket::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    q[i] = buf[m++];
    spin[i] = (int) ubuf(buf[m++]).i;
    eradius[i] = buf[m++];
    etag[i] = (int) ubuf(buf[m++]).i;
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecWavepacket::unpack_border_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    q[i] = buf[m++];
    spin[i] = (int) ubuf(buf[m++]).i;
    eradius[i] = buf[m++];
    etag[i] = (int) ubuf(buf[m++]).i;
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    ervel[i] = buf[m++];
    cs[2*i] = buf[m++];
    cs[2*i+1] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecWavepacket::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    q[i] = buf[m++];
    spin[i] = (int) ubuf(buf[m++]).i;
    eradius[i] = buf[m++];
    etag[i] = (int) ubuf(buf[m++]).i;
    ervel[i] = buf[m++];
    cs[2*i] = buf[m++];
    cs[2*i+1] = buf[m++];
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecWavepacket::pack_exchange(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;

  buf[m++] = q[i];
  buf[m++] = ubuf(spin[i]).d;
  buf[m++] = eradius[i];
  buf[m++] = ervel[i];

  buf[m++] = ubuf(etag[i]).d;
  buf[m++] = cs[2*i];
  buf[m++] = cs[2*i+1];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecWavepacket::unpack_exchange(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;

  q[nlocal] = buf[m++];
  spin[nlocal] = (int) ubuf(buf[m++]).i;
  eradius[nlocal] = buf[m++];
  ervel[nlocal] = buf[m++];

  etag[nlocal] = (int) ubuf(buf[m++]).i;
  cs[2*nlocal] = buf[m++];
  cs[2*nlocal+1] = buf[m++];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecWavepacket::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = 18 * nlocal;        // Associated with pack_restart

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecWavepacket::pack_restart(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  buf[m++] = q[i];
  buf[m++] = ubuf(spin[i]).d;
  buf[m++] = eradius[i];
  buf[m++] = ervel[i];

  buf[m++] = ubuf(etag[i]).d;
  buf[m++] = cs[2*i];
  buf[m++] = cs[2*i+1];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecWavepacket::unpack_restart(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  q[nlocal] = buf[m++];
  spin[nlocal] = (int) ubuf(buf[m++]).i;
  eradius[nlocal] = buf[m++];
  ervel[nlocal] = buf[m++];

  etag[nlocal] = (int) ubuf(buf[m++]).i;
  cs[2*nlocal] = buf[m++];
  cs[2*nlocal+1] = buf[m++];

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
   AWPMD: creates a proton
------------------------------------------------------------------------- */

void AtomVecWavepacket::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) |
    ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  q[nlocal] = 1.;
  spin[nlocal] = 0;
  eradius[nlocal] = 0.0;
  ervel[nlocal] = 0.0;

  etag[nlocal] = 0;
  cs[2*nlocal] = 0.;
  cs[2*nlocal+1] = 0.;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
   AWPMD: 0-tag 1-type 2-q 3-spin 4-eradius 5-etag 6-cs_re 7-cs_im
------------------------------------------------------------------------- */

void AtomVecWavepacket::data_atom(double *coord, imageint imagetmp, 
                                  char **values)
{
  int nlocal = atom->nlocal;

  if (nlocal == nmax) grow(0);

  tag[nlocal] = ATOTAGINT(values[0]);
  type[nlocal] = atoi(values[1]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  q[nlocal] = atof(values[2]);
  spin[nlocal] = atoi(values[3]);
  eradius[nlocal] = atof(values[4]);
  if (eradius[nlocal] < 0.0)
    error->one(FLERR,"Invalid eradius in Atoms section of data file");

  etag[nlocal] = atoi(values[5]);
  cs[2*nlocal] = atoi(values[6]);
  cs[2*nlocal+1] = atof(values[7]);

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  ervel[nlocal] = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecWavepacket::data_atom_hybrid(int nlocal, char **values)
{
  q[nlocal] = atof(values[0]);
  spin[nlocal] = atoi(values[1]);
  eradius[nlocal] = atof(values[2]);
  if (eradius[nlocal] < 0.0)
    error->one(FLERR,"Invalid eradius in Atoms section of data file");

  etag[nlocal] = atoi(values[3]);
  cs[2*nlocal] = atoi(values[4]);
  cs[2*nlocal+1] = atof(values[5]);

  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  ervel[nlocal] = 0.0;

  return 3;
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecWavepacket::data_vel(int m, char **values)
{
  v[m][0] = atof(values[0]);
  v[m][1] = atof(values[1]);
  v[m][2] = atof(values[2]);
  ervel[m] = atof(values[3]);
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Velocities section of data file
------------------------------------------------------------------------- */

int AtomVecWavepacket::data_vel_hybrid(int m, char **values)
{
  ervel[m] = atof(values[0]);
  return 1;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecWavepacket::pack_data(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(type[i]).d;
    buf[i][2] = q[i];
    buf[i][3] = ubuf(spin[i]).d;
    buf[i][4] = eradius[i];
    buf[i][5] = ubuf(etag[i]).d;
    buf[i][6] = cs[2*i];
    buf[i][7] = cs[2*i+1];
    buf[i][8] = x[i][0];
    buf[i][9] = x[i][1];
    buf[i][10] = x[i][2];
    buf[i][11] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][12] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][13] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecWavepacket::pack_data_hybrid(int i, double *buf)
{
  buf[0] = q[i];
  buf[1] = ubuf(spin[i]).d;
  buf[2] = eradius[i];
  buf[3] = ubuf(etag[i]).d;
  buf[4] = cs[2*i];
  buf[5] = cs[2*i+1];
  return 6;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecWavepacket::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT 
            " %d %-1.16e %d %-1.16e %d %-1.16e %-1.16e %-1.16e "
            "%-1.16e %-1.16e %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            buf[i][2],(int) ubuf(buf[i][3]).i,buf[i][4],
            (int) ubuf(buf[i][5]).i,buf[i][6],buf[i][8],
            buf[i][8],buf[i][9],buf[i][10],
            (int) ubuf(buf[i][11]).i,(int) ubuf(buf[i][12]).i,
            (int) ubuf(buf[i][13]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecWavepacket::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %d %-1.16e %d %-1.16e %-1.16e",
          buf[0],(int) ubuf(buf[1]).i,buf[2],(int) ubuf(buf[3]).i,
          buf[4],buf[5]);
  return 6;
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

void AtomVecWavepacket::pack_vel(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = v[i][0];
    buf[i][2] = v[i][1];
    buf[i][3] = v[i][2];
    buf[i][4] = ervel[i];
  }
}

/* ----------------------------------------------------------------------
   pack hybrid velocity info for data file
------------------------------------------------------------------------- */

int AtomVecWavepacket::pack_vel_hybrid(int i, double *buf)
{
  buf[0] = ervel[i];
  return 1;
}

/* ----------------------------------------------------------------------
   write velocity info to data file
------------------------------------------------------------------------- */

void AtomVecWavepacket::write_vel(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " %-1.16e %-1.16e %-1.16e %-1.16e\n",
            (tagint) ubuf(buf[i][0]).i,buf[i][1],buf[i][2],buf[i][3],buf[i][4]);
}

/* ----------------------------------------------------------------------
   write hybrid velocity info to data file
------------------------------------------------------------------------- */

int AtomVecWavepacket::write_vel_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e",buf[0]);
  return 1;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecWavepacket::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("q")) bytes += memory->usage(q,nmax);
  if (atom->memcheck("spin")) bytes += memory->usage(spin,nmax);
  if (atom->memcheck("eradius")) bytes += memory->usage(eradius,nmax);
  if (atom->memcheck("ervel")) bytes += memory->usage(ervel,nmax);
  if (atom->memcheck("erforce"))
    bytes += memory->usage(erforce,nmax*comm->nthreads);

  if (atom->memcheck("ervelforce")) bytes += memory->usage(ervelforce,nmax);
  if (atom->memcheck("cs")) bytes += memory->usage(cs,2*nmax);
  if (atom->memcheck("csforce")) bytes += memory->usage(csforce,2*nmax);
  if (atom->memcheck("vforce")) bytes += memory->usage(vforce,3*nmax);
  if (atom->memcheck("etag")) bytes += memory->usage(etag,nmax);

  return bytes;
}
