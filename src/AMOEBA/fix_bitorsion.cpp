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

#include "fix_bitorsion.h"

#include <cmath>

#include <cstring>
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define BITORSIONMAX 6   // max # of BiTorsion terms stored by one atom
#define LISTDELTA 10000
#define LB_FACTOR 1.5
#define MAXLINE 1024

// NOTE: extra until figure things out

int tnx(int k) {
  return 0;
}
int tny(int k) {
  return 0;
}
int ttx(int i, int j) {
  return 0;
}
int tty(int i, int j) {
  return 0;
}
int ttxy(int i,int j) {
  return 0;
}
double tbf(int i,int j) {
  return 0.0;
}
double tbx(int i,int j) {
  return 0.0;
}
double tby(int i,int j) {
  return 0.0;
}
double tbxy(int i,int j) {
  return 0.0;
}

void chkttor(int ib, int ic, int id, int sign, double value1, double value2) {
}

void bcuint1(double *ftt, double *ft1, double *ft2, double *ft12,
             double x1l, double x1u, double y1l, double y1u, 
             double value1, double value2, 
             double e, double dedang1, double dedang2) {
}

/* ---------------------------------------------------------------------- */

FixBiTorsion::FixBiTorsion(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  bitorsion_list(nullptr), num_bitorsion(nullptr), bitorsion_type(nullptr),
  bitorsion_atom1(nullptr), bitorsion_atom2(nullptr), bitorsion_atom3(nullptr),
  bitorsion_atom4(nullptr), bitorsion_atom5(nullptr)
{
  if (narg != 4) error->all(FLERR,"Illegal fix bitorsion command");

  restart_global = 1;
  restart_peratom = 1;
  energy_global_flag = energy_peratom_flag = 1;
  virial_global_flag = virial_peratom_flag = 1;
  thermo_energy = thermo_virial = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  peratom_freq = 1;
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  wd_header = 1;
  wd_section = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // read and setup BiTorsion grid data

  read_grid_data(arg[3]);

  // perform initial allocation of atom-based arrays

  num_bitorsion = nullptr;
  bitorsion_type = nullptr;
  bitorsion_atom1 = nullptr;
  bitorsion_atom2 = nullptr;
  bitorsion_atom3 = nullptr;
  bitorsion_atom4 = nullptr;
  bitorsion_atom5 = nullptr;

  // register with Atom class

  nmax_previous = 0;
  grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART);

  // list of all bitorsions to compute on this proc

  nbitorsion_list = 0;
  max_bitorsion_list = 0;
  bitorsion_list = nullptr;

  // zero thermo energy

  ebitorsion = 0.0;
}

/* --------------------------------------------------------------------- */

FixBiTorsion::~FixBiTorsion()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,Atom::GROW);
  atom->delete_callback(id,Atom::RESTART);

  // per-atom data

  memory->destroy(num_bitorsion);
  memory->destroy(bitorsion_type);
  memory->destroy(bitorsion_atom1);
  memory->destroy(bitorsion_atom2);
  memory->destroy(bitorsion_atom3);
  memory->destroy(bitorsion_atom4);
  memory->destroy(bitorsion_atom5);

  // local list of bitorsions to compute

  memory->destroy(bitorsion_list);

  // BiTorsion grid data

  delete [] nxgrid;
  delete [] nygrid;
  for (int itype = 1; itype <= ntypes; itype++)
    memory->destroy(btgrid[itype]);
  delete [] btgrid;
}

/* ---------------------------------------------------------------------- */

int FixBiTorsion::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;
  mask |= PRE_REVERSE;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBiTorsion::init()
{
  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }

  // check if PairAmoeba disabled bitorsion terms

  Pair *pair = force->pair_match("amoeba",1,0);

  if (!pair) disable = 0;
  else {
    int tmp;
    int flag = *((int *) pair->extract("bitorsion_flag",tmp));
    disable = flag ? 0 : 1;
  }

  // constant

  onefifth = 0.2;
}

/* --------------------------------------------------------------------- */

void FixBiTorsion::setup(int vflag)
{
  pre_neighbor();

  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* --------------------------------------------------------------------- */

void FixBiTorsion::setup_pre_neighbor()
{
  pre_neighbor();
}

/* --------------------------------------------------------------------- */

void FixBiTorsion::setup_pre_reverse(int eflag, int vflag)
{
  pre_reverse(eflag,vflag);
}

/* --------------------------------------------------------------------- */

void FixBiTorsion::min_setup(int vflag)
{
  pre_neighbor();
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   create local list of bitorsions
   if one or more atoms in bitorsion are on this proc,
     this proc lists the bitorsion exactly once
------------------------------------------------------------------------- */

void FixBiTorsion::pre_neighbor()
{
  int i,m,atom1,atom2,atom3,atom4,atom5;

  // guesstimate initial length of local bitorsion list
  // if nbitorsions was not set (due to read_restart, no read_data),
  //   then list will grow by LISTDELTA chunks

  if (max_bitorsion_list == 0) {
    if (nprocs == 1) max_bitorsion_list = nbitorsions;
    else max_bitorsion_list = 
           static_cast<int> (LB_FACTOR*nbitorsions/nprocs);
    memory->create(bitorsion_list,max_bitorsion_list,6,
                   "bitorsion:bitorsion_list");
  }

  int nlocal = atom->nlocal;
  nbitorsion_list = 0;

  for (i = 0; i < nlocal; i++) {
    for (m = 0; m < num_bitorsion[i]; m++) {
      atom1 = atom->map(bitorsion_atom1[i][m]);
      atom2 = atom->map(bitorsion_atom2[i][m]);
      atom3 = atom->map(bitorsion_atom3[i][m]);
      atom4 = atom->map(bitorsion_atom4[i][m]);
      atom5 = atom->map(bitorsion_atom5[i][m]);

      if (atom1 == -1 || atom2 == -1 || atom3 == -1 ||
          atom4 == -1 || atom5 == -1)
        error->one(FLERR,"BiTorsion atoms {} {} {} {} {} {} missing on "
                   "proc {} at step {}",
                   bitorsion_atom1[i][m],bitorsion_atom2[i][m],
                   bitorsion_atom3[i][m],bitorsion_atom4[i][m],
                   bitorsion_atom5[i][m],me,update->ntimestep);
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      atom4 = domain->closest_image(i,atom4);
      atom5 = domain->closest_image(i,atom5);

      if (i <= atom1 && i <= atom2 && i <= atom3 &&
          i <= atom4 && i <= atom5) {
        if (nbitorsion_list == max_bitorsion_list) {
          max_bitorsion_list += LISTDELTA;
          memory->grow(bitorsion_list,max_bitorsion_list,6,
                       "bitorsion:bitorsion_list");
        }
        bitorsion_list[nbitorsion_list][0] = atom1;
        bitorsion_list[nbitorsion_list][1] = atom2;
        bitorsion_list[nbitorsion_list][2] = atom3;
        bitorsion_list[nbitorsion_list][3] = atom4;
        bitorsion_list[nbitorsion_list][4] = atom5;
        bitorsion_list[nbitorsion_list][5] = bitorsion_type[i][m];
        nbitorsion_list++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   store eflag, so can use it in post_force to tally per-atom energies
------------------------------------------------------------------------- */

void FixBiTorsion::pre_reverse(int eflag, int /*vflag*/)
{
  eflag_caller = eflag;
}

/* ---------------------------------------------------------------------- */

void FixBiTorsion::post_force(int vflag)
{
  if (disable) return;

  int ia,ib,ic,id,ie;
  int nlo,nhi,nt;
  int xlo,ylo;
  int pos1,pos2;
  double e,fgrp,sign;
  double angle1,angle2;
  double value1,value2;
  double cosine1,cosine2;
  double xt,yt,zt,rt2;
  double xu,yu,zu,ru2;
  double xv,yv,zv,rv2;
  double rtru,rurv;
  double x1l,x1u;
  double y1l,y1u;
  double xia,yia,zia;
  double xib,yib,zib;
  double xic,yic,zic;
  double xid,yid,zid;
  double xie,yie,zie;
  double xba,yba,zba;
  double xdc,ydc,zdc;
  double xcb,ycb,zcb;
  double xed,yed,zed;
  double rcb,rdc;
  double xca,yca,zca;
  double xdb,ydb,zdb;
  double xec,yec,zec;
  double dedang1,dedang2;
  double dedxt,dedyt,dedzt;
  double dedxu,dedyu,dedzu;
  double dedxu2,dedyu2,dedzu2;
  double dedxv2,dedyv2,dedzv2;
  double dedxia,dedyia,dedzia;
  double dedxib,dedyib,dedzib;
  double dedxic,dedyic,dedzic;
  double dedxid,dedyid,dedzid;
  double dedxib2,dedyib2,dedzib2;
  double dedxic2,dedyic2,dedzic2;
  double dedxid2,dedyid2,dedzid2;
  double dedxie2,dedyie2,dedzie2;
  double vxx,vyy,vzz;
  double vyx,vzx,vzy;
  double vxx2,vyy2,vzz2;
  double vyx2,vzx2,vzy2;
  double ftt[4],ft12[4];
  double ft1[4],ft2[4];

  // NOTE: extra until figure everything out

  int k,btype;
  double radian;
  double engfraction;
  int nlist,list[6];
  double v[6];

  // END of NOTE

  ebitorsion = 0.0;
  int eflag = eflag_caller;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;

  for (int n = 0; n < nbitorsion_list; n++) {

    ia = bitorsion_list[n][0];
    ib = bitorsion_list[n][1];
    ic = bitorsion_list[n][2];
    id = bitorsion_list[n][3];
    ie = bitorsion_list[n][4];

    // NOTE: is a btype ever used, i.e. as index into spline tables?

    btype = bitorsion_list[n][5];

    xia = x[ia][0];
    yia = x[ia][1];
    zia = x[ia][2];
    xib = x[ib][0];
    yib = x[ib][1];
    zib = x[ib][2];
    xic = x[ic][0];
    yic = x[ic][1];
    zic = x[ic][2];
    xid = x[id][0];
    yid = x[id][1];
    zid = x[id][2];
    xie = x[ie][0];
    yie = x[ie][1];
    zie = x[ie][2];

    xba = xib - xia;
    yba = yib - yia;
    zba = zib - zia;
    xcb = xic - xib;
    ycb = yic - yib;
    zcb = zic - zib;
    xdc = xid - xic;
    ydc = yid - yic;
    zdc = zid - zic;
    xed = xie - xid;
    yed = yie - yid;
    zed = zie - zid;

    xt = yba*zcb - ycb*zba;
    yt = zba*xcb - zcb*xba;
    zt = xba*ycb - xcb*yba;
    xu = ycb*zdc - ydc*zcb;
    yu = zcb*xdc - zdc*xcb;
    zu = xcb*ydc - xdc*ycb;

    rt2 = xt*xt + yt*yt + zt*zt;
    ru2 = xu*xu + yu*yu + zu*zu;
    rtru = sqrt(rt2 * ru2);
    xv = ydc*zed - yed*zdc;
    yv = zdc*xed - zed*xdc;
    zv = xdc*yed - xed*ydc;
    rv2 = xv*xv + yv*yv + zv*zv;
    rurv = sqrt(ru2 * rv2);

    if (rtru <= 0.0 || rurv <= 0.0) continue;

    rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
    cosine1 = (xt*xu + yt*yu + zt*zu) / rtru;
    cosine1 = MIN(1.0,MAX(-1.0,cosine1));
    angle1 = radian * acos(cosine1);
    sign = xba*xu + yba*yu + zba*zu;
    if (sign < 0.0)  angle1 = -angle1;
    value1 = angle1;

    rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc);
    cosine2 = (xu*xv + yu*yv + zu*zv) / rurv;
    cosine2 = MIN(1.0,MAX(-1.0,cosine2));
    angle2 = radian * acos(cosine2);
    sign = xcb*xv + ycb*yv + zcb*zv;
    if (sign < 0.0)  angle2 = -angle2;
    value2 = angle2;

    // check for inverted chirality at the central atom

    chkttor(ib,ic,id,sign,value1,value2);

    // use bicubic interpolation to compute spline values
    // NOTE: need to worry about C vs Fortran indexing
    //       both here and in the methods called
    // NOTE: make sure pos1 and pos2 are ints

    nlo = 1;
    nhi = tnx(k);
    while (nhi-nlo > 1) {
      nt = (nhi+nlo) / 2;
      if (ttx(nt,k) > value1) nhi = nt;
      else nlo = nt;
    }

    xlo = nlo;
    nlo = 1;
    nhi = tny(k);
    while (nhi-nlo > 1) {
      nt = (nhi + nlo)/2;
      if (tty(nt,k) > value2) nhi = nt;
      else nlo = nt;
    }
       
    ylo = nlo;
    x1l = ttx(xlo,k);
    x1u = ttx(xlo+1,k);
    y1l = tty(ylo,k);
    y1u = tty(ylo+1,k);
    pos2 = ylo*tnx(k) + xlo;
    pos1 = pos2 - tnx(k);

    ftt[0] = tbf(pos1,k);
    ftt[1] = tbf(pos1+1,k);
    ftt[2] = tbf(pos2+1,k);
    ftt[3] = tbf(pos2,k);

    ft1[0] = tbx(pos1,k);
    ft1[1] = tbx(pos1+1,k);
    ft1[2] = tbx(pos2+1,k);
    ft1[3] = tbx(pos2,k);

    ft2[0] = tby(pos1,k);
    ft2[1] = tby(pos1+1,k);
    ft2[2] = tby(pos2+1,k);
    ft2[3] = tby(pos2,k);

    ft12[0] = tbxy(pos1,k);
    ft12[1] = tbxy(pos1+1,k);
    ft12[2] = tbxy(pos2+1,k);
    ft12[3] = tbxy(pos2,k);

    bcuint1(ftt,ft1,ft2,ft12,x1l,x1u,y1l,y1u,value1,value2,
            e,dedang1,dedang2);

    // NOTE: remove ttorunit if 1.0 ?
    // NOTE: what are radian units ?
    // NOTE: value of sign is set twice above ??

    double ttorunit = 1.0;
    e *= ttorunit;
    dedang1 = sign * ttorunit * radian * dedang1;
    dedang2 = sign * ttorunit * radian * dedang2;

    // fraction of energy for each atom

    engfraction = e * onefifth;

    // chain rule terms for first angle derivative components

    xca = xic - xia;
    yca = yic - yia;
    zca = zic - zia;
    xdb = xid - xib;
    ydb = yid - yib;
    zdb = zid - zib;

    dedxt = dedang1 * (yt*zcb - ycb*zt) / (rt2*rcb);
    dedyt = dedang1 * (zt*xcb - zcb*xt) / (rt2*rcb);
    dedzt = dedang1 * (xt*ycb - xcb*yt) / (rt2*rcb);
    dedxu = -dedang1 * (yu*zcb - ycb*zu) / (ru2*rcb);
    dedyu = -dedang1 * (zu*xcb - zcb*xu) / (ru2*rcb);
    dedzu = -dedang1 * (xu*ycb - xcb*yu) / (ru2*rcb);
    
    // compute first derivative components for first angle

    dedxia = zcb*dedyt - ycb*dedzt;
    dedyia = xcb*dedzt - zcb*dedxt;
    dedzia = ycb*dedxt - xcb*dedyt;
    dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu;
    dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu;
    dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu;
    dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu;
    dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu;
    dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu;
    dedxid = zcb*dedyu - ycb*dedzu;
    dedyid = xcb*dedzu - zcb*dedxu;
    dedzid = ycb*dedxu - xcb*dedyu;

    // chain rule terms for second angle derivative components
    
    xec = xie - xic;
    yec = yie - yic;
    zec = zie - zic;
    
    dedxu2 = dedang2 * (yu*zdc - ydc*zu) / (ru2*rdc);
    dedyu2 = dedang2 * (zu*xdc - zdc*xu) / (ru2*rdc);
    dedzu2 = dedang2 * (xu*ydc - xdc*yu) / (ru2*rdc);
    dedxv2 = -dedang2 * (yv*zdc - ydc*zv) / (rv2*rdc);
    dedyv2 = -dedang2 * (zv*xdc - zdc*xv) / (rv2*rdc);
    dedzv2 = -dedang2 * (xv*ydc - xdc*yv) / (rv2*rdc);
      
    // compute first derivative components for second angle

    dedxib2 = zdc*dedyu2 - ydc*dedzu2;
    dedyib2 = xdc*dedzu2 - zdc*dedxu2;
    dedzib2 = ydc*dedxu2 - xdc*dedyu2;
    dedxic2 = ydb*dedzu2 - zdb*dedyu2 + zed*dedyv2 - yed*dedzv2;
    dedyic2 = zdb*dedxu2 - xdb*dedzu2 + xed*dedzv2 - zed*dedxv2;
    dedzic2 = xdb*dedyu2 - ydb*dedxu2 + yed*dedxv2 - xed*dedyv2;
    dedxid2 = zcb*dedyu2 - ycb*dedzu2 + yec*dedzv2 - zec*dedyv2;
    dedyid2 = xcb*dedzu2 - zcb*dedxu2 + zec*dedxv2 - xec*dedzv2;
    dedzid2 = ycb*dedxu2 - xcb*dedyu2 + xec*dedyv2 - yec*dedxv2;
    dedxie2 = zdc*dedyv2 - ydc*dedzv2;
    dedyie2 = xdc*dedzv2 - zdc*dedxv2;
    dedzie2 = ydc*dedxv2 - xdc*dedyv2;

    // increment the torsion-torsion energy and gradient

    if (ia < nlocal) {
      ebitorsion += engfraction;
      f[ia][0] += dedxia;
      f[ia][1] += dedyia;
      f[ia][2] += dedzia;
    }

    if (ib < nlocal) {
      ebitorsion += engfraction;
      f[ib][0] += dedxib + dedxib2;
      f[ib][1] += dedyib + dedyib2;
      f[ib][2] += dedzib + dedzib2;
    }

    if (ic < nlocal) {
      ebitorsion += engfraction;
      f[ic][0] += dedxic + dedxic2;
      f[ic][1] += dedyic + dedyic2;
      f[ic][2] += dedzic + dedzic2;
    }

    if (id < nlocal) {
      ebitorsion += engfraction;
      f[id][0] += dedxid + dedxid2;
      f[id][1] += dedyid + dedyid2;
      f[id][2] += dedzid + dedzid2;
    }

    if (ie < nlocal) {
      ebitorsion += engfraction;
      f[ie][0] += dedxie2;
      f[ie][1] += dedyie2;
      f[ie][2] += dedzie2;
    }

    // increment the internal virial tensor components

    if (evflag) {
      nlist = 0;
      if (ia < nlocal) list[nlist++] = ia;
      if (ib < nlocal) list[nlist++] = ib;
      if (ic < nlocal) list[nlist++] = ic;
      if (id < nlocal) list[nlist++] = id;
      if (ie < nlocal) list[nlist++] = ie;

      vxx = xcb*(dedxic+dedxid) - xba*dedxia + xdc*dedxid;
      vyx = ycb*(dedxic+dedxid) - yba*dedxia + ydc*dedxid;
      vzx = zcb*(dedxic+dedxid) - zba*dedxia + zdc*dedxid;
      vyy = ycb*(dedyic+dedyid) - yba*dedyia + ydc*dedyid;
      vzy = zcb*(dedyic+dedyid) - zba*dedyia + zdc*dedyid;
      vzz = zcb*(dedzic+dedzid) - zba*dedzia + zdc*dedzid;
      vxx2 = xdc*(dedxid2+dedxie2) - xcb*dedxib2 + xed*dedxie2;
      vyx2 = ydc*(dedxid2+dedxie2) - ycb*dedxib2 + yed*dedxie2;
      vzx2 = zdc*(dedxid2+dedxie2) - zcb*dedxib2 + zed*dedxie2;
      vyy2 = ydc*(dedyid2+dedyie2) - ycb*dedyib2 + yed*dedyie2;
      vzy2 = zdc*(dedyid2+dedyie2) - zcb*dedyib2 + zed*dedyie2;
      vzz2 = zdc*(dedzid2+dedzie2) - zcb*dedzib2 + zed*dedzie2;
    
      v[0] += vxx + vxx2;
      v[1] += vyy + vyy2;
      v[2] += vzz + vzz2;
      v[3] += vyx + vyx2;
      v[4] += vzx + vzx2;
      v[5] += vzy + vzy2;
      
      ev_tally(nlist,list,5.0,e,v);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBiTorsion::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixBiTorsion::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of BiTorision term
------------------------------------------------------------------------- */

double FixBiTorsion::compute_scalar()
{
  double all;
  MPI_Allreduce(&ebitorsion,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// methods to read BiTorsion grid file, perform interpolation
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

void FixBiTorsion::read_grid_data(char *bitorsion_file)
{
  char line[MAXLINE];
  char *eof;

  FILE *fp = nullptr;
  if (me == 0) {
    fp = utils::open_potential(bitorsion_file,lmp,nullptr);
    if (fp == nullptr)
      error->one(FLERR,"Cannot open fix bitorsion file {}: {}",
                 bitorsion_file, utils::getsyserror());

    eof = fgets(line,MAXLINE,fp);
    eof = fgets(line,MAXLINE,fp);
    eof = fgets(line,MAXLINE,fp);
    if (eof == nullptr) error->one(FLERR,"Unexpected end of fix bitorsion file");

    sscanf(line,"%d",&ntypes);
  }

  MPI_Bcast(&ntypes,1,MPI_INT,0,world);
  if (ntypes == 0) error->all(FLERR,"Fix bitorsion file has no types");

  btgrid = new double***[ntypes+1];
  nxgrid = new int[ntypes+1];
  nygrid = new int[ntypes+1];

  // read one array for each BiTorsion type from file

  int tmp,nx,ny;
  double xgrid,ygrid,value;

  for (int itype = 1; itype <= ntypes; itype++) {
    if (me == 0) {
      eof = fgets(line,MAXLINE,fp);
      eof = fgets(line,MAXLINE,fp);
      if (eof == nullptr) 
        error->one(FLERR,"Unexpected end of fix bitorsion file");
      sscanf(line,"%d %d %d",&tmp,&nx,&ny);
    }

    MPI_Bcast(&nx,1,MPI_INT,0,world);
    MPI_Bcast(&ny,1,MPI_INT,0,world);
    nxgrid[itype] = nx;
    nygrid[itype] = ny;

    memory->create(btgrid[itype],nx,ny,3,"bitorsion:btgrid");

    // NOTE: should read this chunk of lines with utils in single read

    if (me == 0) {
      for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
          eof = fgets(line,MAXLINE,fp);
          if (eof == nullptr) 
            error->one(FLERR,"Unexpected end of fix bitorsion file");
          sscanf(line,"%lg %lg %lg",&xgrid,&ygrid,&value);
          btgrid[itype][ix][iy][0] = xgrid;
          btgrid[itype][ix][iy][1] = ygrid;
          btgrid[itype][ix][iy][2] = value;
        }
      }
    }

    MPI_Bcast(&btgrid[itype][0][0][0],nx*ny*3,MPI_DOUBLE,0,world);
  }

  if (me == 0) fclose(fp);
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// methods to read and write data file
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

void FixBiTorsion::read_data_header(char *line)
{
  if (strstr(line,"bitorsions")) {
    sscanf(line,BIGINT_FORMAT,&nbitorsions);
  } else error->all(FLERR,"Invalid read data header line for fix bitorsion");
}

/* ----------------------------------------------------------------------
   unpack N lines in buf from section of data file labeled by keyword
   id_offset is applied to atomID fields if multiple data files are read
------------------------------------------------------------------------- */

void FixBiTorsion::read_data_section(char *keyword, int n, char *buf,
                                     tagint id_offset)
{
  int m,tmp,itype;
  tagint atom1,atom2,atom3,atom4,atom5;
  char *next;

  next = strchr(buf,'\n');
  *next = '\0';
  int nwords = utils::count_words(utils::trim_comment(buf));
  *next = '\n';

  if (nwords != 7)
    error->all(FLERR,"Incorrect {} format in data file",keyword);

  // loop over lines of BiTorsions
  // tokenize the line into values
  // each atom in the BiTorsion stores it

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    sscanf(buf,"%d %d " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT
           " " TAGINT_FORMAT " " TAGINT_FORMAT,
           &tmp,&itype,&atom1,&atom2,&atom3,&atom4,&atom5);

    atom1 += id_offset;
    atom2 += id_offset;
    atom3 += id_offset;
    atom4 += id_offset;
    atom5 += id_offset;

    if ((m = atom->map(atom1)) >= 0) {
      if (num_bitorsion[m] == BITORSIONMAX)
        error->one(FLERR,"Too many BIORSIONS for one atom");
      bitorsion_type[m][num_bitorsion[m]] = itype;
      bitorsion_atom1[m][num_bitorsion[m]] = atom1;
      bitorsion_atom2[m][num_bitorsion[m]] = atom2;
      bitorsion_atom3[m][num_bitorsion[m]] = atom3;
      bitorsion_atom4[m][num_bitorsion[m]] = atom4;
      bitorsion_atom5[m][num_bitorsion[m]] = atom5;
      num_bitorsion[m]++;
    }

    if ((m = atom->map(atom2)) >= 0) {
      if (num_bitorsion[m] == BITORSIONMAX)
        error->one(FLERR,"Too many BIORSIONS for one atom");
      bitorsion_type[m][num_bitorsion[m]] = itype;
      bitorsion_atom1[m][num_bitorsion[m]] = atom1;
      bitorsion_atom2[m][num_bitorsion[m]] = atom2;
      bitorsion_atom3[m][num_bitorsion[m]] = atom3;
      bitorsion_atom4[m][num_bitorsion[m]] = atom4;
      bitorsion_atom5[m][num_bitorsion[m]] = atom5;
      num_bitorsion[m]++;
    }

    if ((m = atom->map(atom3)) >= 0) {
      if (num_bitorsion[m] == BITORSIONMAX)
        error->one(FLERR,"Too many BIORSIONS for one atom");
      bitorsion_type[m][num_bitorsion[m]] = itype;
      bitorsion_atom1[m][num_bitorsion[m]] = atom1;
      bitorsion_atom2[m][num_bitorsion[m]] = atom2;
      bitorsion_atom3[m][num_bitorsion[m]] = atom3;
      bitorsion_atom4[m][num_bitorsion[m]] = atom4;
      bitorsion_atom5[m][num_bitorsion[m]] = atom5;
      num_bitorsion[m]++;
    }

    if ((m = atom->map(atom4)) >= 0) {
      if (num_bitorsion[m] == BITORSIONMAX)
        error->one(FLERR,"Too many BIORSIONS for one atom");
      bitorsion_type[m][num_bitorsion[m]] = itype;
      bitorsion_atom1[m][num_bitorsion[m]] = atom1;
      bitorsion_atom2[m][num_bitorsion[m]] = atom2;
      bitorsion_atom3[m][num_bitorsion[m]] = atom3;
      bitorsion_atom4[m][num_bitorsion[m]] = atom4;
      bitorsion_atom5[m][num_bitorsion[m]] = atom5;
      num_bitorsion[m]++;
    }

    if ((m = atom->map(atom5)) >= 0) {
      if (num_bitorsion[m] == BITORSIONMAX)
        error->one(FLERR,"Too many BIORSIONS for one atom");
      bitorsion_type[m][num_bitorsion[m]] = itype;
      bitorsion_atom1[m][num_bitorsion[m]] = atom1;
      bitorsion_atom2[m][num_bitorsion[m]] = atom2;
      bitorsion_atom3[m][num_bitorsion[m]] = atom3;
      bitorsion_atom4[m][num_bitorsion[m]] = atom4;
      bitorsion_atom5[m][num_bitorsion[m]] = atom5;
      num_bitorsion[m]++;
    }

    buf = next + 1;
  }
}

/* ---------------------------------------------------------------------- */

bigint FixBiTorsion::read_data_skip_lines(char * /*keyword*/)
{
  return nbitorsions;
}

/* ----------------------------------------------------------------------
   write Mth header line to file
   only called by proc 0
------------------------------------------------------------------------- */

void FixBiTorsion::write_data_header(FILE *fp, int /*mth*/)
{
  fprintf(fp,BIGINT_FORMAT " bitorsions\n",nbitorsions);
}

/* ----------------------------------------------------------------------
   return size I own for Mth data section
   # of data sections = 1 for this fix
  // nx = # of BiTorsions owned by my local atoms
  //   only atom3 owns the BiTorsion in this context
  // ny = 6 columns = type + 5 atom IDs
------------------------------------------------------------------------- */

void FixBiTorsion::write_data_section_size(int /*mth*/, int &nx, int &ny)
{
  int i,m;

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  nx = 0;
  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_bitorsion[i]; m++)
      if (bitorsion_atom3[i][m] == tag[i]) nx++;
  
  ny = 6;
}

/* ----------------------------------------------------------------------
   pack values for Mth data section into 2d buf
   buf allocated by caller as owned BiTorsions by 6
------------------------------------------------------------------------- */

void FixBiTorsion::write_data_section_pack(int /*mth*/, double **buf)
{
  int i,m;

  // 1st column = BiTorsion type
  // 2nd-6th columns = 5 atom IDs

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int n = 0;
  for (i = 0; i < nlocal; i++) {
    for (m = 0; m < num_bitorsion[i]; m++) {
      if (bitorsion_atom3[i][m] != tag[i]) continue;
      buf[n][0] = ubuf(bitorsion_type[i][m]).d;
      buf[n][1] = ubuf(bitorsion_atom1[i][m]).d;
      buf[n][2] = ubuf(bitorsion_atom2[i][m]).d;
      buf[n][3] = ubuf(bitorsion_atom3[i][m]).d;
      buf[n][4] = ubuf(bitorsion_atom4[i][m]).d;
      buf[n][5] = ubuf(bitorsion_atom5[i][m]).d;
      n++;
    }
  }
}

/* ----------------------------------------------------------------------
   write section keyword for Mth data section to file
   use Molecules or Charges if that is only field, else use fix ID
   only called by proc 0
------------------------------------------------------------------------- */

void FixBiTorsion::write_data_section_keyword(int /*mth*/, FILE *fp)
{
  fprintf(fp,"\nBiTorsions\n\n");
}

/* ----------------------------------------------------------------------
   write N lines from buf to file
   convert buf fields to int or double depending on styles
   index can be used to prepend global numbering
   only called by proc 0
------------------------------------------------------------------------- */

void FixBiTorsion::write_data_section(int /*mth*/, FILE *fp,
                                  int n, double **buf, int index)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,"%d %d " TAGINT_FORMAT " " TAGINT_FORMAT
            " " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT "\n",
            index+i,(int) ubuf(buf[i][0]).i,(tagint) ubuf(buf[i][1]).i,
            (tagint) ubuf(buf[i][2]).i,(tagint) ubuf(buf[i][3]).i,
            (tagint) ubuf(buf[i][4]).i,(tagint) ubuf(buf[i][5]).i);
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// methods for restart and communication
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixBiTorsion::write_restart(FILE *fp)
{
  if (comm->me == 0) {
    int size = sizeof(bigint);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(&nbitorsions,sizeof(bigint),1,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixBiTorsion::restart(char *buf)
{
  nbitorsions = *((bigint *) buf);
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixBiTorsion::pack_restart(int i, double *buf)
{
  int n = 1;
  for (int m = 0; m < num_bitorsion[i]; m++) {
    buf[n++] = ubuf(MAX(bitorsion_type[i][m],-bitorsion_type[i][m])).d;
    buf[n++] = ubuf(bitorsion_atom1[i][m]).d;
    buf[n++] = ubuf(bitorsion_atom2[i][m]).d;
    buf[n++] = ubuf(bitorsion_atom3[i][m]).d;
    buf[n++] = ubuf(bitorsion_atom4[i][m]).d;
    buf[n++] = ubuf(bitorsion_atom5[i][m]).d;
  }
  buf[0] = n;

  return n;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixBiTorsion::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

   int n = 0;
   for (int i = 0; i < nth; i++) n += static_cast<int> (extra[nlocal][n]);

   int count = static_cast<int> (extra[nlocal][n++]);
   num_bitorsion[nlocal] = (count-1)/6;

   for (int m = 0; m < num_bitorsion[nlocal]; m++) {
     bitorsion_type[nlocal][m] = (int) ubuf(extra[nlocal][n++]).i;
     bitorsion_atom1[nlocal][m] = (tagint) ubuf(extra[nlocal][n++]).i;
     bitorsion_atom2[nlocal][m] = (tagint) ubuf(extra[nlocal][n++]).i;
     bitorsion_atom3[nlocal][m] = (tagint) ubuf(extra[nlocal][n++]).i;
     bitorsion_atom4[nlocal][m] = (tagint) ubuf(extra[nlocal][n++]).i;
     bitorsion_atom5[nlocal][m] = (tagint) ubuf(extra[nlocal][n++]).i;
   }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixBiTorsion::maxsize_restart()
{
  return 1 + BITORSIONMAX*6;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixBiTorsion::size_restart(int nlocal)
{
  return 1 + num_bitorsion[nlocal]*6;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixBiTorsion::grow_arrays(int nmax)
{
  num_bitorsion = memory->grow(num_bitorsion,nmax,"cmap:num_bitorsion");
  bitorsion_type = memory->grow(bitorsion_type,nmax,BITORSIONMAX,
                                "cmap:bitorsion_type");
  bitorsion_atom1 = memory->grow(bitorsion_atom1,nmax,BITORSIONMAX,
                                 "cmap:bitorsion_atom1");
  bitorsion_atom2 = memory->grow(bitorsion_atom2,nmax,BITORSIONMAX,
                                 "cmap:bitorsion_atom2");
  bitorsion_atom3 = memory->grow(bitorsion_atom3,nmax,BITORSIONMAX,
                                 "cmap:bitorsion_atom3");
  bitorsion_atom4 = memory->grow(bitorsion_atom4,nmax,BITORSIONMAX,
                                 "cmap:bitorsion_atom4");
  bitorsion_atom5 = memory->grow(bitorsion_atom5,nmax,BITORSIONMAX,
                                 "cmap:bitorsion_atom5");

  // must initialize num_bitorsion to 0 for added atoms
  // may never be set for some atoms when data file is read

  for (int i = nmax_previous; i < nmax; i++) num_bitorsion[i] = 0;
  nmax_previous = nmax;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixBiTorsion::copy_arrays(int i, int j, int /*delflag*/)
{
  num_bitorsion[j] = num_bitorsion[i];

  for (int k = 0; k < num_bitorsion[j]; k++) {
    bitorsion_type[j][k] = bitorsion_type[i][k];
    bitorsion_atom1[j][k] = bitorsion_atom1[i][k];
    bitorsion_atom2[j][k] = bitorsion_atom2[i][k];
    bitorsion_atom3[j][k] = bitorsion_atom3[i][k];
    bitorsion_atom4[j][k] = bitorsion_atom4[i][k];
    bitorsion_atom5[j][k] = bitorsion_atom5[i][k];
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixBiTorsion::set_arrays(int i)
{
  num_bitorsion[i] = 0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixBiTorsion::pack_exchange(int i, double *buf)
{
  int n = 0;
  buf[n++] = ubuf(num_bitorsion[i]).d;
  for (int m = 0; m < num_bitorsion[i]; m++) {
    buf[n++] = ubuf(bitorsion_type[i][m]).d;
    buf[n++] = ubuf(bitorsion_atom1[i][m]).d;
    buf[n++] = ubuf(bitorsion_atom2[i][m]).d;
    buf[n++] = ubuf(bitorsion_atom3[i][m]).d;
    buf[n++] = ubuf(bitorsion_atom4[i][m]).d;
    buf[n++] = ubuf(bitorsion_atom5[i][m]).d;
  }
  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixBiTorsion::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  num_bitorsion[nlocal] = (int) ubuf(buf[n++]).i;
  for (int m = 0; m < num_bitorsion[nlocal]; m++) {
    bitorsion_type[nlocal][m] = (int) ubuf(buf[n++]).i;
    bitorsion_atom1[nlocal][m] = (tagint) ubuf(buf[n++]).i;
    bitorsion_atom2[nlocal][m] = (tagint) ubuf(buf[n++]).i;
    bitorsion_atom3[nlocal][m] = (tagint) ubuf(buf[n++]).i;
    bitorsion_atom4[nlocal][m] = (tagint) ubuf(buf[n++]).i;
    bitorsion_atom5[nlocal][m] = (tagint) ubuf(buf[n++]).i;
  }
  return n;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixBiTorsion::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = (double)nmax * sizeof(int);              // num_bitorsion
  bytes += (double)nmax*BITORSIONMAX * sizeof(int);      // bitorsion_type
  bytes += (double)5*nmax*BITORSIONMAX * sizeof(int);     // bitorsion_atom 12345
  bytes += (double)6*max_bitorsion_list * sizeof(int);    // bitorsion_list
  return bytes;
}
