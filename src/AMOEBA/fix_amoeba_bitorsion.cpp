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

#include "fix_amoeba_bitorsion.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "pair.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

static constexpr int BITORSIONMAX = 6;   // max # of BiTorsion terms stored by one atom
static constexpr int LISTDELTA = 10000;
static constexpr double LB_FACTOR = 1.5;
static constexpr int MAXLINE = 1024;

// spline weighting factors

static constexpr double WT[16][16] = {
  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {-3.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0,-2.0, 0.0, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0},
  { 2.0, 0.0, 0.0,-2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
  { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
  { 0.0, 0.0, 0.0, 0.0,-3.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0,-2.0, 0.0, 0.0,-1.0},
  { 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0,-2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0},
  {-3.0, 3.0, 0.0, 0.0,-2.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,-3.0, 3.0, 0.0, 0.0,-2.0,-1.0, 0.0, 0.0},
  { 9.0,-9.0, 9.0,-9.0, 6.0, 3.0,-3.0,-6.0, 6.0,-6.0,-3.0, 3.0, 4.0, 2.0, 1.0, 2.0},
  {-6.0, 6.0,-6.0, 6.0,-4.0,-2.0, 2.0, 4.0,-3.0, 3.0, 3.0,-3.0,-2.0,-1.0,-1.0,-2.0},
  { 2.0,-2.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0,-2.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
  {-6.0, 6.0,-6.0, 6.0,-3.0,-3.0, 3.0, 3.0,-4.0, 4.0, 2.0,-2.0,-2.0,-2.0,-1.0,-1.0},
  { 4.0,-4.0, 4.0,-4.0, 2.0, 2.0,-2.0,-2.0, 2.0,-2.0,-2.0, 2.0, 1.0, 1.0, 1.0, 1.0}
};

/* ---------------------------------------------------------------------- */

FixAmoebaBiTorsion::FixAmoebaBiTorsion(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), num_bitorsion(nullptr), bitorsion_type(nullptr), bitorsion_atom1(nullptr),
  bitorsion_atom2(nullptr), bitorsion_atom3(nullptr), bitorsion_atom4(nullptr),
  bitorsion_atom5(nullptr), bitorsion_list(nullptr)
{
  if (narg != 4) error->all(FLERR,"Illegal fix amoeba/bitorsion command");

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
  create_splines();

  // border comm of 1st neighbors in special list
  // so chkttor() can use them when central atom = ghost to check chirality
  // comm_border = max # of bonds per atom + 1 for count

  comm_border = 7;
  atom->add_callback(Atom::BORDER);

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

FixAmoebaBiTorsion::~FixAmoebaBiTorsion()
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

  delete[] nxgrid;
  delete[] nygrid;
  for (int itype = 1; itype <= nbitypes; itype++) {
    memory->destroy(ttx[itype]);
    memory->destroy(tty[itype]);
    memory->destroy(tbf[itype]);
    memory->destroy(tbx[itype]);
    memory->destroy(tby[itype]);
    memory->destroy(tbxy[itype]);
  }
  delete[] ttx;
  delete[] tty;
  delete[] tbf;
  delete[] tbx;
  delete[] tby;
  delete[] tbxy;
}

/* ---------------------------------------------------------------------- */

int FixAmoebaBiTorsion::setmask()
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

void FixAmoebaBiTorsion::init()
{
  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }

  // error check that PairAmoeba or PairHiippo exist

  pair = nullptr;
  pair = force->pair_match("^amoeba",0,0);
  if (!pair) pair = force->pair_match("^hippo",0,0);

  if (!pair)
    error->all(FLERR,"Cannot use fix amoeba/bitorsion w/out pair amoeba/hippo");

  // check if PairAmoeba or PairHippo disabled bitorsion terms

  int tmp;
  int flag = *((int *) pair->extract("bitorsion_flag",tmp));
  disable = flag ? 0 : 1;

  // constant

  onefifth = 0.2;
}

/* --------------------------------------------------------------------- */

void FixAmoebaBiTorsion::setup(int vflag)
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

void FixAmoebaBiTorsion::setup_pre_neighbor()
{
  pre_neighbor();
}

/* --------------------------------------------------------------------- */

void FixAmoebaBiTorsion::setup_pre_reverse(int eflag, int vflag)
{
  pre_reverse(eflag,vflag);
}

/* --------------------------------------------------------------------- */

void FixAmoebaBiTorsion::min_setup(int vflag)
{
  pre_neighbor();
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   create local list of bitorsions
   if one or more atoms in bitorsion are on this proc,
     this proc lists the bitorsion exactly once
------------------------------------------------------------------------- */

void FixAmoebaBiTorsion::pre_neighbor()
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

void FixAmoebaBiTorsion::pre_reverse(int eflag, int /*vflag*/)
{
  eflag_caller = eflag;
}

/* ---------------------------------------------------------------------- */

void FixAmoebaBiTorsion::post_force(int vflag)
{
  if (disable) return;

  int ia,ib,ic,id,ie,btype;
  int nx,ny,nlo,nhi,nt;
  int xlo,ylo;
  int pos1,pos2;
  double e,sign,dot;
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

  double engfraction;
  int nlist,list[6];
  double v[6];

  double radian2degree = 180.0 / MY_PI;

  ebitorsion = 0.0;
  int eflag = eflag_caller;
  ev_init(eflag,vflag);

  // set current ptrs to PairAmoeba amtype and atomic_num

  int tmp;
  amtype = (int *) pair->extract("amtype",tmp);
  atomic_num = (int *) pair->extract("atomic_num",tmp);

  // loop over local bitorsions

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;

  for (int n = 0; n < nbitorsion_list; n++) {

    ia = bitorsion_list[n][0];
    ib = bitorsion_list[n][1];
    ic = bitorsion_list[n][2];
    id = bitorsion_list[n][3];
    ie = bitorsion_list[n][4];
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
    angle1 = radian2degree * acos(cosine1);
    dot = xba*xu + yba*yu + zba*zu;
    if (dot < 0.0) angle1 = -angle1;
    value1 = angle1;

    rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc);
    cosine2 = (xu*xv + yu*yv + zu*zv) / rurv;
    cosine2 = MIN(1.0,MAX(-1.0,cosine2));
    angle2 = radian2degree * acos(cosine2);
    dot = xcb*xv + ycb*yv + zcb*zv;
    if (dot < 0.0) angle2 = -angle2;
    value2 = angle2;

    // check for inverted chirality at the central atom
    // inputs = ib,ic,id
    // outputs = sign,value1,value2

    chkttor(ib,ic,id,sign,value1,value2);

    // 2 binary searches to find location of angles 1,2 in grid
    // ttx,tty are 0-indexed here, 1-indexed in Tinker
    // xlo,ylo = final location, each one less than in Tinker

    nx = nxgrid[btype];
    ny = nygrid[btype];

    nlo = 0;
    nhi = nx-1;

    while (nhi-nlo > 1) {
      nt = (nhi+nlo) / 2;
      if (ttx[btype][nt] > value1) nhi = nt;
      else nlo = nt;
    }
    xlo = nlo;

    nlo = 0;
    nhi = ny-1;
    while (nhi-nlo > 1) {
      nt = (nhi + nlo)/2;
      if (tty[btype][nt] > value2) nhi = nt;
      else nlo = nt;
    }
    ylo = nlo;

    // fill ftt,ft1,ft2,ft12 vecs with spline coeffs near xlo,ylo grid pt
    // ttx,tty,tbf,tbx,tby,tbxy are 0-indexed here, 1-indexed in Tinker
    // xlo,ylo,pos1,pos2 are all one less than in Tinker

    x1l = ttx[btype][xlo];
    x1u = ttx[btype][xlo+1];
    y1l = tty[btype][ylo];
    y1u = tty[btype][ylo+1];

    pos2 = (ylo+1)*nx + xlo;
    pos1 = pos2 - nx;

    ftt[0] = tbf[btype][pos1];
    ftt[1] = tbf[btype][pos1+1];
    ftt[2] = tbf[btype][pos2+1];
    ftt[3] = tbf[btype][pos2];

    ft1[0] = tbx[btype][pos1];
    ft1[1] = tbx[btype][pos1+1];
    ft1[2] = tbx[btype][pos2+1];
    ft1[3] = tbx[btype][pos2];

    ft2[0] = tby[btype][pos1];
    ft2[1] = tby[btype][pos1+1];
    ft2[2] = tby[btype][pos2+1];
    ft2[3] = tby[btype][pos2];

    ft12[0] = tbxy[btype][pos1];
    ft12[1] = tbxy[btype][pos1+1];
    ft12[2] = tbxy[btype][pos2+1];
    ft12[3] = tbxy[btype][pos2];

    // bicuint1() uses bicubic interpolation to compute interpolated values
    // outputs = e,dedang1,dedang2

    bcuint1(ftt,ft1,ft2,ft12,x1l,x1u,y1l,y1u,value1,value2,
            e,dedang1,dedang2);

    dedang1 = sign * radian2degree * dedang1;
    dedang2 = sign * radian2degree * dedang2;

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
      f[ia][0] -= dedxia;
      f[ia][1] -= dedyia;
      f[ia][2] -= dedzia;
    }

    if (ib < nlocal) {
      ebitorsion += engfraction;
      f[ib][0] -= dedxib + dedxib2;
      f[ib][1] -= dedyib + dedyib2;
      f[ib][2] -= dedzib + dedzib2;
    }

    if (ic < nlocal) {
      ebitorsion += engfraction;
      f[ic][0] -= dedxic + dedxic2;
      f[ic][1] -= dedyic + dedyic2;
      f[ic][2] -= dedzic + dedzic2;
    }

    if (id < nlocal) {
      ebitorsion += engfraction;
      f[id][0] -= dedxid + dedxid2;
      f[id][1] -= dedyid + dedyid2;
      f[id][2] -= dedzid + dedzid2;
    }

    if (ie < nlocal) {
      ebitorsion += engfraction;
      f[ie][0] -= dedxie2;
      f[ie][1] -= dedyie2;
      f[ie][2] -= dedzie2;
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

      v[0] = -vxx - vxx2;
      v[1] = -vyy - vyy2;
      v[2] = -vzz - vzz2;
      v[3] = -vyx - vyx2;
      v[4] = -vzx - vzx2;
      v[5] = -vzy - vzy2;

      ev_tally(nlist,list,5.0,e,v);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixAmoebaBiTorsion::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAmoebaBiTorsion::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of BiTorsion term
------------------------------------------------------------------------- */

double FixAmoebaBiTorsion::compute_scalar()
{
  double all;
  MPI_Allreduce(&ebitorsion,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// methods to read BiTorsion grid file, spline grids, perform interpolation
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   read grid data from bitorsion_file produced by tinker2lmp.py
   one entry for each biotorsion type
   when complete:
     nbitypes = # of bitorsion types
     nxgrid,nygrid = x,y dimensions of grid for each type
     ttx,tty = vectors of x,y angles for each type
               length = nx or ny, 0-indexed
     tbf = vector of 2d grid values for each type
           length = nx*ny, 0-indexed, x varies fastest
     ttx,tty,tbf are similar to Tinker data structs
     here they are 0-indexed, in Tinker they are 1-indexed
------------------------------------------------------------------------- */

void FixAmoebaBiTorsion::read_grid_data(char *bitorsion_file)
{
  char line[MAXLINE] = {'\0'};
  char *eof;

  FILE *fp = nullptr;
  if (me == 0) {
    fp = utils::open_potential(bitorsion_file,lmp,nullptr);
    if (fp == nullptr)
      error->one(FLERR,"Cannot open fix amoeba/bitorsion file {}: {}",
                 bitorsion_file, utils::getsyserror());

    eof = fgets(line,MAXLINE,fp);
    eof = fgets(line,MAXLINE,fp);
    eof = fgets(line,MAXLINE,fp);
    if (eof == nullptr)
      error->one(FLERR,"Unexpected end of fix amoeba/bitorsion file");

    sscanf(line,"%d",&nbitypes);
  }

  MPI_Bcast(&nbitypes,1,MPI_INT,0,world);
  if (nbitypes == 0) error->all(FLERR,"Fix amoeba/bitorsion file has no types");

  // allocate data structs
  // type index ranges from 1 to Nbitypes, so allocate one larger

  nxgrid = new int[nbitypes+1];
  nygrid = new int[nbitypes+1];
  ttx = new double*[nbitypes+1];
  tty = new double*[nbitypes+1];
  tbf = new double*[nbitypes+1];

  // read one array for each BiTorsion type from file

  int tmp,nx,ny;
  double xgrid,ygrid,value;

  for (int itype = 1; itype <= nbitypes; itype++) {
    if (me == 0) {
      eof = fgets(line,MAXLINE,fp);
      eof = fgets(line,MAXLINE,fp);
      if (eof == nullptr)
        error->one(FLERR,"Unexpected end of fix amoeba/bitorsion file");
      sscanf(line,"%d %d %d",&tmp,&nx,&ny);
    }

    MPI_Bcast(&nx,1,MPI_INT,0,world);
    MPI_Bcast(&ny,1,MPI_INT,0,world);
    nxgrid[itype] = nx;
    nygrid[itype] = ny;

    memory->create(ttx[itype],nx,"bitorsion:ttx");
    memory->create(tty[itype],ny,"bitorsion:tty");
    memory->create(tbf[itype],nx*ny,"bitorsion:tbf");

    // NOTE: could read this chunk of lines with utils in single read?

    if (me == 0) {
      for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
          eof = fgets(line,MAXLINE,fp);
          if (eof == nullptr)
            error->one(FLERR,"Unexpected end of fix amoeba/bitorsion file");
          sscanf(line,"%lg %lg %lg",&xgrid,&ygrid,&value);
          if (iy == 0) ttx[itype][ix] = xgrid;
          if (ix == 0) tty[itype][iy] = ygrid;
          tbf[itype][iy*nx+ix] = value;
        }
      }
    }

    MPI_Bcast(ttx[itype],nx,MPI_DOUBLE,0,world);
    MPI_Bcast(tty[itype],ny,MPI_DOUBLE,0,world);
    MPI_Bcast(tbf[itype],nx*ny,MPI_DOUBLE,0,world);
  }

  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   create spline data structs for each bitorsion type
   based on Tinker ktortor.f
   when complete:
     tbx,tby = spline coeffs
     tbxy = vector of 2d grid values for each type, x varies fastest
     tbx,tby,tbxy are similar to Tinker data structs
     here they are 0-indexed, in Tinker they are 1-indexed
------------------------------------------------------------------------- */

void FixAmoebaBiTorsion::create_splines()
{
  int nx,ny;

  // allocate work vectors for cspline() and nspline() methods
  // all are 0-indexed here and in Tinker
  // tmp1,tmp2 = (x,y) inputs to spline methods
  // bs = retained output from spline methods
  // cs,ds,tmp3-7 = additional outputs from spline methods, not retained
  // allocate to max length of any grid dimension for all types

  double *bs,*cs,*ds;
  double *tmp1,*tmp2;
  double *tmp3,*tmp4,*tmp5,*tmp6,*tmp7;

  int maxdim = 0;
  for (int itype = 1; itype <= nbitypes; itype++) {
    maxdim = MAX(maxdim,nxgrid[itype]);
    maxdim = MAX(maxdim,nygrid[itype]);
  }

  memory->create(bs,maxdim,"bitorsion:bs");
  memory->create(cs,maxdim,"bitorsion:cs");
  memory->create(ds,maxdim,"bitorsion:ds");
  memory->create(tmp1,maxdim,"bitorsion:tmp1");
  memory->create(tmp2,maxdim,"bitorsion:tmp2");
  memory->create(tmp3,maxdim,"bitorsion:tmp3");
  memory->create(tmp4,maxdim,"bitorsion:tmp4");
  memory->create(tmp5,maxdim,"bitorsion:tmp5");
  memory->create(tmp6,maxdim,"bitorsion:tmp6");
  memory->create(tmp7,maxdim,"bitorsion:tmp7");

  // allocate data structs
  // type index ranges from 1 to Nbitypes, so allocate one larger

  tbx = new double*[nbitypes+1];
  tby = new double*[nbitypes+1];
  tbxy = new double*[nbitypes+1];

  for (int itype = 1; itype <= nbitypes; itype++) {

    nx = nxgrid[itype];
    ny = nygrid[itype];

    // cyclic = 1 if angle range is -180.0 to 180.0 in both dims
    // error if cyclic and x,y pairs of edges of 2d array values do not match
    // equality comparisons are within eps

    int cyclic = 1;
    double eps = 1.0e-6;

    if (fabs(fabs(ttx[itype][0]-ttx[itype][nx-1]) - 360.0) > eps) cyclic = 0;
    if (fabs(fabs(tty[itype][0]-tty[itype][ny-1]) - 360.0) > eps) cyclic = 0;

    if (cyclic) {
      // do error check on matching edge values
    }

    // allocate nx*ny vectors for itype

    memory->create(tbx[itype],nx*ny,"bitorsion:tbx");
    memory->create(tby[itype],nx*ny,"bitorsion:tby");
    memory->create(tbxy[itype],nx*ny,"bitorsion:tbxy");

    // spline fit of derivatives about 1st torsion

    for (int i = 0; i < nx; i++)
      tmp1[i] = ttx[itype][i];

    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++)
        tmp2[i] = tbf[itype][j*nx+i];

      if (cyclic)
        cspline(nx-1,tmp1,tmp2,bs,cs,ds,tmp3,tmp4,tmp5,tmp6,tmp7);
      else
        nspline(nx-1,tmp1,tmp2,bs,cs,tmp3,tmp4,tmp5,tmp6,tmp7);

      for (int i = 0; i < nx; i++)
        tbx[itype][j*nx+i] = bs[i];
    }

    // spline fit of derivatives about 2nd torsion

    for (int j = 0; j < ny; j++)
      tmp1[j] = ttx[itype][j];

    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++)
        tmp2[j] = tbf[itype][j*nx+i];

      if (cyclic)
        cspline(ny-1,tmp1,tmp2,bs,cs,ds,tmp3,tmp4,tmp5,tmp6,tmp7);
      else
        nspline(ny-1,tmp1,tmp2,bs,cs,tmp3,tmp4,tmp5,tmp6,tmp7);

      for (int j = 0; j < ny; j++)
        tby[itype][j*nx+i] = bs[j];
    }

    // spline fit of cross derivatives about both torsions

    for (int j = 0; j < ny; j++)
      tmp1[j] = ttx[itype][j];

    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++)
        tmp2[j] = tbx[itype][j*nx+i];

      if (cyclic)
        cspline(ny-1,tmp1,tmp2,bs,cs,ds,tmp3,tmp4,tmp5,tmp6,tmp7);
      else
        nspline(ny-1,tmp1,tmp2,bs,cs,tmp3,tmp4,tmp5,tmp6,tmp7);

      for (int j = 0; j < ny; j++)
        tbxy[itype][j*nx+i] = bs[j];
    }
  }

  // free work vectors local to this method

  memory->destroy(bs);
  memory->destroy(cs);
  memory->destroy(ds);
  memory->destroy(tmp1);
  memory->destroy(tmp2);
  memory->destroy(tmp3);
  memory->destroy(tmp4);
  memory->destroy(tmp5);
  memory->destroy(tmp6);
  memory->destroy(tmp7);
}

/* ----------------------------------------------------------------------
   Tinker method nspline()
   computes coefficients for an nonperiodic cubic spline
     with natural boundary conditions where the first and last second
     derivatives are already known
   all vectors are of length n+1 and are indexed from 0 to n inclusive
   n,x0,y0 are inputs
   rest of args are outputs
------------------------------------------------------------------------- */

void FixAmoebaBiTorsion::nspline(int n, double *x0, double *y0,
                                 double *s1, double *s2,
                                 double *h, double *g, double *dy,
                                 double *dla, double *dmu)
{
  int i;
  double t;

  // set first and last second deriviatives to zero

  double y21 = 0.0;
  double y2n = 0.0;

  // find the intervals to be used

  for (i = 0; i <= n-1; i++) {
    h[i] = x0[i+1] - x0[i];
    dy[i] = (y0[i+1]-y0[i]) / h[i];
  }

  // calculate the spline coeffcients

  for (i = 1; i <= n-1; i++) {
    dla[i] = h[i] / (h[i]+h[i-1]);
    dmu[i] = 1.0 - dla[i];
    g[i] = 3.0 * (dla[i]*dy[i-1]+dmu[i]*dy[i]);
  }

  // set the initial value via natural boundary condition

  dla[n] = 1.0;
  dla[0] = 0.0;
  dmu[n] = 0.0;
  dmu[0] = 1.0;
  g[0] = 3.0*dy[0] - 0.5*h[0]*y21;
  g[n] = 3.0*dy[n-1] + 0.5*h[n-1]*y2n;

  // solve the triagonal system of linear equations

  dmu[0] = 0.5 * dmu[0];
  g[0] = 0.5 * g[0];

  for (i = 1; i <= n; i++) {
    t = 2.0 - dmu[i-1]*dla[i];
    dmu[i] = dmu[i] / t;
    g[i] = (g[i]-g[i-1]*dla[i]) / t;
  }
  for (i = n-1; i >= 0; i--)
    g[i] = g[i] - dmu[i]*g[i+1];

  // get the first derivative at each grid point

  for (i = 0; i <= n; i++)
    s1[i] = g[i];

  // get the second derivative at each grid point

  s2[0] = y21;
  s2[n] = y2n;
  for (i = 1; i <= n-1; i++)
    s2[i] = 6.0*(y0[i+1]-y0[i])/(h[i]*h[i]) - 4.0*s1[i]/h[i] - 2.0*s1[i+1]/h[i];
}

/* ----------------------------------------------------------------------
   Tinker method cspline()
   computes coefficients for an periodic interpolating cubic spline
   all vectors are of length n+1 and are indexed from 0 to n inclusive
   n,xn,fn are inputs
   rest of args are outputs
------------------------------------------------------------------------- */

void FixAmoebaBiTorsion::cspline(int n, double *xn, double *fn,
                                 double *b, double *c, double *d,
                                 double *h, double *du, double *dm,
                                 double *rc, double *rs)
{
  int i;
  double temp1,temp2;

  double average = 0.5 * (fn[0] + fn[n]);
  fn[0] = average;
  fn[n] = average;

  // get auxiliary variables and matrix elements on first call

  for (i = 0; i < n; i++)
    h[i] = xn[i+1] - xn[i];
  h[n] = h[0];

  for (i = 1; i < n; i++)
    du[i] = h[i];
  du[n] = h[0];

  for (i = 1; i <= n; i++)
    dm[i] = 2.0 * (h[i-1]+h[i]);

  // compute the right hand side

  temp1 = (fn[1]-fn[0]) / h[0];
  for (i = 1; i < n; i++) {
    temp2 = (fn[i+1]-fn[i]) / h[i];
    rs[i]  = 3.0 * (temp2-temp1);
    temp1 = temp2;
  }
  rs[n] = 3.0 * ((fn[1]-fn[0])/h[0]-temp1);

  // solve the linear system with factorization

  int iflag;
  cytsy(n,dm,du,rc,rs,c,iflag);
  if (iflag != 1) return;

  // compute remaining spline coefficients

  c[0] = c[n];
  for (i = 0; i < n; i++) {
    b[i] = (fn[i+1]-fn[i])/h[i] - h[i]/3.0*(c[i+1]+2.0*c[i]);
    d[i] = (c[i+1]-c[i]) / (3.0*h[i]);
  }
  b[n] = (fn[1]-fn[n])/h[n] - h[n]/3.0*(c[1]+2.0*c[n]);
}

/* ----------------------------------------------------------------------
   Tinker method cytsy()
   solve cyclic tridiagonal system
   all vectors are of length n+1 and are indexed from 0 to n inclusive
   n,dm,du are inputs
   du,cr,rs,x,iflag are outputs
------------------------------------------------------------------------- */

void FixAmoebaBiTorsion::cytsy(int n, double *dm, double *du,
                               double *cr, double *rs, double *x, int &iflag)
{
  iflag = -2;
  if (n < 3) return;
  cytsyp(n,dm,du,cr,iflag);

  // update and back substitute as necessary

  if (iflag == 1) cytsys(n,dm,du,cr,rs,x);
}

/* ----------------------------------------------------------------------
   Tinker method ctsys()
   tridiagonal Cholesky factorization
   all vectors are of length n+1 and are indexed from 0 to n inclusive
   n,dm,du are inputs
   du,cr,iflag are outputs
------------------------------------------------------------------------- */

void FixAmoebaBiTorsion::cytsyp(int n, double *dm, double *du,
                                double *cr, int &iflag)
{
  int i;
  double temp1,temp2;

  // set error bound and test for condition n greater than 2

  double eps = 1.0e-8;

  iflag = -2;
  if (n < 3) return;

  // checking to see if matrix is positive definite

  double row = fabs(dm[1]) + fabs(du[1]) + fabs(du[n]);

  if (row == 0.0) {
    iflag = 0;
    return;
  }

  double d = 1.0 / row;

  if (dm[1] < 0.0) {
    iflag = -1;
    return;
  } else if (fabs(dm[1])*d <= eps) {
    iflag = 0;
    return;
  }

  // factoring a while checking for a positive definite and
  //   strong nonsingular matrix a

  temp1 = du[1];
  du[1] = du[1] / dm[1];
  cr[1] = du[n] / dm[1];

  for (i = 2; i < n; i++) {
    row = fabs(dm[i]) + fabs(du[i]) + fabs(temp1);
    if (row == 0.0) {
      iflag = 0;
      return;
    }
    d = 1.0 / row;
    dm[i] = dm[i] - temp1*du[i-1];
    if (dm[i] < 0.0) {
      iflag = -1;
      return;
    } else if (fabs(dm[i])*d <= eps) {
      iflag = 0;
      return;
    }
    if (i < n-1) {
      cr[i] = -temp1 * cr[i-1] / dm[i];
      temp1 = du[i];
      du[i] = du[i] / dm[i];
    } else {
      temp2 = du[i];
      du[i] = (du[i] - temp1*cr[i-1]) / dm[i];
    }
  }

  row = fabs(du[n]) + fabs(dm[n]) + fabs(temp2);
  if (row == 0.0) {
    iflag = 0;
    return;
  }

  d = 1.0 / row;
  dm[n] = dm[n] - dm[n-1]*du[n-1]*du[n-1];
  temp1 = 0.0;

  for (i = 1; i < n-1; i++)
    temp1 += dm[i]*cr[i]*cr[i];

  dm[n] = dm[n] - temp1;

  if (dm[n] < 0.0) {
    iflag = -1;
    return;
  } else if (fabs(dm[n])*d <= eps) {
    iflag = 0;
    return;
  }

  iflag = 1;
}

/* ----------------------------------------------------------------------
   Tinker method cytsys()
   tridiagonal solution from factors
   all vectors are of length n+1 and are indexed from 0 to n inclusive
   n,dm,du,cr are inputs
   rs,x are outputs
------------------------------------------------------------------------- */

void FixAmoebaBiTorsion::cytsys(int n, double *dm, double *du,
                                double *cr, double *rs, double *x)
{
  int i;

  // updating phase

  double temp = rs[1];
  rs[1] = temp / dm[1];
  double sum = cr[1] * temp;

  for (i = 2; i < n; i++) {
    temp = rs[i] - du[i-1]*temp;
    rs[i] = temp / dm[i];
    if (i != n-1)  sum += cr[i]*temp;
  }

  temp = rs[n] - du[n-1]*temp;
  temp = temp - sum;
  rs[n] = temp / dm[n];

  // back substitution phase

  x[n] = rs[n];
  x[n-1] = rs[n-1] - du[n-1]*x[n];
  for (i = n-2; i > 0; i--)
    x[i] = rs[i] - du[i]*x[i+1] - cr[i]*x[n];
}

/* ----------------------------------------------------------------------
   chkttor tests the attached atoms at a torsion-torsion central
   site and inverts the angle values if the site is chiral
------------------------------------------------------------------------- */

void FixAmoebaBiTorsion::chkttor(int ib, int ic, int id,
                                 double &sign, double &value1, double &value2)
{
  int i,ia,ilocal,jlocal,klocal,jtype,ktype;
  tagint j,k,m;
  double xac,yac,zac;
  double xbc,ybc,zbc;
  double xdc,ydc,zdc;
  double c1,c2,c3,vol;

  // test for chirality at the central torsion-torsion site

  sign = 1.0;
  if (atom->nspecial[ic][0] != 4) return;

  tagint **special = atom->special;
  tagint *tag = atom->tag;

  // j,k,m are atom IDs

  j = 0;
  for (i = 0; i < 4; i++) {
    m = special[ic][i];
    if (m != tag[ib] && m != tag[id]) {
      if (j == 0) j = m;
      else k = m;
    }
  }

  // convert atom IDs j,k to local indices jlocal,klocal closest to atom IC

  jlocal = atom->map(j);
  jlocal = domain->closest_image(ic,jlocal);
  jtype = amtype[jlocal];

  klocal = atom->map(k);
  klocal = domain->closest_image(ic,klocal);
  ktype = amtype[klocal];

  // atom ilocal = jlocal or klocal (or -1)
  // set atom IA = ilocal

  ilocal = -1;
  if (jtype > ktype) ilocal = jlocal;
  if (ktype > jtype) ilocal = klocal;
  if (atomic_num[jtype] > atomic_num[ktype]) ilocal = jlocal;
  if (atomic_num[ktype] > atomic_num[jtype]) ilocal = klocal;
  if (ilocal < 0) return;
  ia = ilocal;

  // compute the signed parallelpiped volume at central site

  double **x = atom->x;

  xac = x[ia][0] - x[ic][0];
  yac = x[ia][1] - x[ic][1];
  zac = x[ia][2] - x[ic][2];
  xbc = x[ib][0] - x[ic][0];
  ybc = x[ib][1] - x[ic][1];
  zbc = x[ib][2] - x[ic][2];
  xdc = x[id][0] - x[ic][0];
  ydc = x[id][1] - x[ic][1];
  zdc = x[id][2] - x[ic][2];
  c1 = ybc*zdc - zbc*ydc;
  c2 = ydc*zac - zdc*yac;
  c3 = yac*zbc - zac*ybc;
  vol = xac*c1 + xbc*c2 + xdc*c3;

  // invert the angle values if chirality has an inverted sign

  if (vol < 0.0) {
    sign = -1.0;
    value1 = -value1;
    value2 = -value2;
  }
}

/* ----------------------------------------------------------------------
   perform bicubic spline interpolation
   based on Tinker bcuint1.f
   input = all args except last 3
   output = ansy,ansy1,ansy2
------------------------------------------------------------------------- */

void FixAmoebaBiTorsion::bcuint1(double *y, double *y1,
                                 double *y2, double *y12,
                                 double x1l, double x1u, double x2l, double x2u,
                                 double x1, double x2,
                                 double &ansy, double &ansy1, double &ansy2)
{
  double c[4][4];

  // get coefficients, then perform bicubic interpolation

  bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c);

  double t = (x1-x1l) / (x1u-x1l);
  double u = (x2-x2l) / (x2u-x2l);

  ansy = ansy1 = ansy2 = 0.0;

  for (int i = 3; i >= 0; i--) {
    ansy = t*ansy + ((c[i][3]*u+c[i][2])*u+c[i][1])*u + c[i][0];
    ansy1 = u*ansy1 + (3.0*c[3][i]*t+2.0*c[2][i])*t + c[1][i];
    ansy2 = t*ansy2 + (3.0*c[i][3]*u+2.0*c[i][2])*u + c[i][1];
  }

  ansy1 /= x1u-x1l;
  ansy2 /= x2u-x2l;
}

/* ----------------------------------------------------------------------
   compute bicubic spline coeffs
   based on Tinker bcucof.f
   input = all args except c
   output = c
------------------------------------------------------------------------- */

void FixAmoebaBiTorsion::bcucof(double *y, double *y1, double *y2, double *y12,
                                double d1, double d2, double c[4][4])
{
  int i,j,k;
  double xx;
  double x[16],cl[16];

  // pack a temporary vector of corner values

  double d1d2 = d1 * d2;
  for (i = 0; i < 4; i++) {
    x[i] = y[i];
    x[i+4] = y1[i] * d1;
    x[i+8] = y2[i] * d2;
    x[i+12] = y12[i] * d1d2;
  }

  // matrix multiply by the stored weight table

  for (i = 0; i < 16; i++) {
    xx = 0.0;
    for (k = 0; k < 16; k++)
      xx += WT[i][k]*x[k];
    cl[i] = xx;
  }

  // unpack the result into the coefficient table

  j = 0;
  for (i = 0; i < 4; i++) {
    for (k = 0; k < 4; k++) {
      c[i][k] = cl[j++];
    }
  }
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// methods to read and write data file
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

void FixAmoebaBiTorsion::read_data_header(char *line)
{
  if (strstr(line,"bitorsions")) {
    sscanf(line,BIGINT_FORMAT,&nbitorsions);
  } else error->all(FLERR,
                    "Invalid read data header line for fix amoeba/bitorsion");
}

/* ----------------------------------------------------------------------
   unpack N lines in buf from section of data file labeled by keyword
   id_offset is applied to atomID fields if multiple data files are read
------------------------------------------------------------------------- */

void FixAmoebaBiTorsion::read_data_section(char *keyword, int n, char *buf,
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

bigint FixAmoebaBiTorsion::read_data_skip_lines(char * /*keyword*/)
{
  return nbitorsions;
}

/* ----------------------------------------------------------------------
   write Mth header line to file
   only called by proc 0
------------------------------------------------------------------------- */

void FixAmoebaBiTorsion::write_data_header(FILE *fp, int /*mth*/)
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

void FixAmoebaBiTorsion::write_data_section_size(int /*mth*/, int &nx, int &ny)
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

void FixAmoebaBiTorsion::write_data_section_pack(int /*mth*/, double **buf)
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

void FixAmoebaBiTorsion::write_data_section_keyword(int /*mth*/, FILE *fp)
{
  fprintf(fp,"\nBiTorsions\n\n");
}

/* ----------------------------------------------------------------------
   write N lines from buf to file
   convert buf fields to int or double depending on styles
   index can be used to prepend global numbering
   only called by proc 0
------------------------------------------------------------------------- */

void FixAmoebaBiTorsion::write_data_section(int /*mth*/, FILE *fp,
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

void FixAmoebaBiTorsion::write_restart(FILE *fp)
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

void FixAmoebaBiTorsion::restart(char *buf)
{
  nbitorsions = *((bigint *) buf);
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixAmoebaBiTorsion::pack_restart(int i, double *buf)
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

void FixAmoebaBiTorsion::unpack_restart(int nlocal, int nth)
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

int FixAmoebaBiTorsion::maxsize_restart()
{
  return 1 + BITORSIONMAX*6;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixAmoebaBiTorsion::size_restart(int nlocal)
{
  return 1 + num_bitorsion[nlocal]*6;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixAmoebaBiTorsion::grow_arrays(int nmax)
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

void FixAmoebaBiTorsion::copy_arrays(int i, int j, int /*delflag*/)
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

void FixAmoebaBiTorsion::set_arrays(int i)
{
  num_bitorsion[i] = 0;
}

/* ----------------------------------------------------------------------
   pack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixAmoebaBiTorsion::pack_border(int n, int *list, double *buf)
{
  int i, j, k;

  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  int m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = ubuf(nspecial[j][0]).d;
    for (k = 0; k < nspecial[j][0]; k++)
      buf[m++] = ubuf(special[j][k]).d;
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixAmoebaBiTorsion::unpack_border(int n, int first, double *buf)
{
  int i, k, last;

  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  int m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    nspecial[i][0] = (int) ubuf(buf[m++]).i;
    for (k = 0; k < nspecial[i][0]; k++)
      special[i][k] = (tagint) ubuf(buf[m++]).i;
  }

  return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixAmoebaBiTorsion::pack_exchange(int i, double *buf)
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

int FixAmoebaBiTorsion::unpack_exchange(int nlocal, double *buf)
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

double FixAmoebaBiTorsion::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = (double)nmax * sizeof(int);              // num_bitorsion
  bytes += (double)nmax*BITORSIONMAX * sizeof(int);       // bitorsion_type
  bytes += (double)5*nmax*BITORSIONMAX * sizeof(int);     // bitorsion_atom 12345
  bytes += (double)6*max_bitorsion_list * sizeof(int);    // bitorsion_list
  return bytes;
}
