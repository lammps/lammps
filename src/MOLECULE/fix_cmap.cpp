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

/* ----------------------------------------------------------------------
   Contributing authors:
   Xiaohu Hu, CMB/ORNL (hux2@ornl.gov)
   David Hyde-Volpe, Tigran Abramyan, and Robert A. Latour (Clemson University)
   Chris Lorenz (Kings College-London)

   Implementation of the CHARMM CMAP; adds an extra energy term for the
   peptide backbone dihedrals.  The tools/ch2lmp/charmm2lammps.pl
   conversion script, which generates an extra section in the LAMMPS data
   file, is needed in order to generate the info used by this fix style.

   References:
   - MacKerell et al., J. Am. Chem. Soc. 126(2004):698-699.
   - MacKerell et al., J. Comput. Chem. 25(2004):1400-1415.
------------------------------------------------------------------------- */

#include "fix_cmap.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "potential_file_reader.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <exception>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define MAXLINE 256
#define LISTDELTA 10000
#define LB_FACTOR 1.5

#define CMAPMAX 6   // max # of CMAP terms stored by one atom
#define CMAPDIM 24  // grid map dimension is 24 x 24
#define CMAPXMIN -360.0
#define CMAPXMIN2 -180.0
#define CMAPDX 15.0 // 360/CMAPDIM

/* ---------------------------------------------------------------------- */

FixCMAP::FixCMAP(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  crosstermlist(nullptr), num_crossterm(nullptr), crossterm_type(nullptr),
  crossterm_atom1(nullptr), crossterm_atom2(nullptr), crossterm_atom3(nullptr),
  crossterm_atom4(nullptr), crossterm_atom5(nullptr),
  g_axis(nullptr), cmapgrid(nullptr), d1cmapgrid(nullptr), d2cmapgrid(nullptr),
  d12cmapgrid(nullptr)
{
  if (narg != 4) error->all(FLERR,"Illegal fix cmap command");

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

  // allocate memory for CMAP data

  memory->create(g_axis,CMAPDIM,"cmap:g_axis");
  memory->create(cmapgrid,6,CMAPDIM,CMAPDIM,"cmap:grid");
  memory->create(d1cmapgrid,6,CMAPDIM,CMAPDIM,"cmap:d1grid");
  memory->create(d2cmapgrid,6,CMAPDIM,CMAPDIM,"cmap:d2grid");
  memory->create(d12cmapgrid,6,CMAPDIM,CMAPDIM,"cmap:d12grid");

  // read and setup CMAP data

  read_grid_map(arg[3]);

  // perform initial allocation of atom-based arrays
  // register with Atom class

  num_crossterm = nullptr;
  crossterm_type = nullptr;
  crossterm_atom1 = nullptr;
  crossterm_atom2 = nullptr;
  crossterm_atom3 = nullptr;
  crossterm_atom4 = nullptr;
  crossterm_atom5 = nullptr;

  nmax_previous = 0;
  FixCMAP::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART);

  // local list of crossterms

  ncmap = 0;
  maxcrossterm = 0;
  crosstermlist = nullptr;
}

/* --------------------------------------------------------------------- */

FixCMAP::~FixCMAP()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,Atom::GROW);
  atom->delete_callback(id,Atom::RESTART);

  memory->destroy(g_axis);
  memory->destroy(cmapgrid);
  memory->destroy(d1cmapgrid);
  memory->destroy(d2cmapgrid);
  memory->destroy(d12cmapgrid);

  memory->destroy(crosstermlist);

  memory->destroy(num_crossterm);
  memory->destroy(crossterm_type);
  memory->destroy(crossterm_atom1);
  memory->destroy(crossterm_atom2);
  memory->destroy(crossterm_atom3);
  memory->destroy(crossterm_atom4);
  memory->destroy(crossterm_atom5);
}

/* ---------------------------------------------------------------------- */

int FixCMAP::setmask()
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

void FixCMAP::init()
{
  int i;
  double angle;

  i = 0;
  angle = -180.0;
  while (angle < 180.0) {
    g_axis[i] = angle;
    angle += CMAPDX;
    i++;
  }

  // pre-compute the derivatives of the maps

  for (i = 0; i < 6; i++)
    set_map_derivatives(cmapgrid[i],d1cmapgrid[i],d2cmapgrid[i],d12cmapgrid[i]);

  // define newton_bond here in case restart file was read (not data file)

  newton_bond = force->newton_bond;

  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* --------------------------------------------------------------------- */

void FixCMAP::setup(int vflag)
{
  pre_neighbor();

  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel_respa);
  }
}

/* --------------------------------------------------------------------- */

void FixCMAP::setup_pre_neighbor()
{
  pre_neighbor();
}

/* --------------------------------------------------------------------- */

void FixCMAP::setup_pre_reverse(int eflag, int vflag)
{
  pre_reverse(eflag,vflag);
}

/* --------------------------------------------------------------------- */

void FixCMAP::min_setup(int vflag)
{
  pre_neighbor();
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   store local neighbor list as if newton_bond = OFF, even if actually ON
------------------------------------------------------------------------- */

void FixCMAP::pre_neighbor()
{
  int i,m,atom1,atom2,atom3,atom4,atom5;

  // guesstimate initial length of local crossterm list
  // if ncmap was not set (due to read_restart, no read_data),
  //   then list will grow by LISTDELTA chunks

  if (maxcrossterm == 0) {
    if (nprocs == 1) maxcrossterm = ncmap;
    else maxcrossterm = static_cast<int> (LB_FACTOR*ncmap/nprocs);
    memory->create(crosstermlist,maxcrossterm,6,"cmap:crosstermlist");
  }

  int nlocal = atom->nlocal;

  ncrosstermlist = 0;

  for (i = 0; i < nlocal; i++) {
    for (m = 0; m < num_crossterm[i]; m++) {
      atom1 = atom->map(crossterm_atom1[i][m]);
      atom2 = atom->map(crossterm_atom2[i][m]);
      atom3 = atom->map(crossterm_atom3[i][m]);
      atom4 = atom->map(crossterm_atom4[i][m]);
      atom5 = atom->map(crossterm_atom5[i][m]);

      if (atom1 == -1 || atom2 == -1 || atom3 == -1 ||
          atom4 == -1 || atom5 == -1)
        error->one(FLERR,"CMAP atoms {} {} {} {} {} missing on "
                                     "proc {} at step {}",
                                     crossterm_atom1[i][m],crossterm_atom2[i][m],
                                     crossterm_atom3[i][m],crossterm_atom4[i][m],
                                     crossterm_atom5[i][m],me,update->ntimestep);
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      atom4 = domain->closest_image(i,atom4);
      atom5 = domain->closest_image(i,atom5);

      if (i <= atom1 && i <= atom2 && i <= atom3 &&
          i <= atom4 && i <= atom5) {
        if (ncrosstermlist == maxcrossterm) {
          maxcrossterm += LISTDELTA;
          memory->grow(crosstermlist,maxcrossterm,6,"cmap:crosstermlist");
        }
        crosstermlist[ncrosstermlist][0] = atom1;
        crosstermlist[ncrosstermlist][1] = atom2;
        crosstermlist[ncrosstermlist][2] = atom3;
        crosstermlist[ncrosstermlist][3] = atom4;
        crosstermlist[ncrosstermlist][4] = atom5;
        crosstermlist[ncrosstermlist][5] = crossterm_type[i][m];
        ncrosstermlist++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   store eflag, so can use it in post_force to tally per-atom energies
------------------------------------------------------------------------- */

void FixCMAP::pre_reverse(int eflag, int /*vflag*/)
{
  eflag_caller = eflag;
}

/* ----------------------------------------------------------------------
   compute CMAP terms as if newton_bond = OFF, even if actually ON
------------------------------------------------------------------------- */

void FixCMAP::post_force(int vflag)
{
  int n,i1,i2,i3,i4,i5,type,nlist;
  int li1, li2, mli1,mli2,mli11,mli21,t1,li3,li4,mli3,mli4,mli31,mli41;
  int list[5];
  // vectors needed to calculate the cross-term dihedral angles
  double vb21x,vb21y,vb21z,vb32x,vb32y,vb32z,vb34x,vb34y,vb34z;
  double vb23x,vb23y,vb23z;
  double vb43x,vb43y,vb43z,vb45x,vb45y,vb45z,a1x,a1y,a1z,b1x,b1y,b1z;
  double a2x,a2y,a2z,b2x,b2y,b2z,r32,a1sq,b1sq,a2sq,b2sq,dpr21r32,dpr34r32;
  double dpr32r43,dpr45r43,r43,vb12x,vb12y,vb12z,vb54x,vb54y,vb54z;
  // cross-term dihedral angles
  double phi,psi,phi1,psi1;
  double f1[3],f2[3],f3[3],f4[3],f5[3],vcmap[6];
  double gs[4],d1gs[4],d2gs[4],d12gs[4];
  double engfraction;
  // vectors needed for the gradient/force calculation
  double dphidr1x,dphidr1y,dphidr1z,dphidr2x,dphidr2y,dphidr2z;
  double dphidr3x,dphidr3y,dphidr3z,dphidr4x,dphidr4y,dphidr4z;
  double dpsidr1x,dpsidr1y,dpsidr1z,dpsidr2x,dpsidr2y,dpsidr2z;
  double dpsidr3x,dpsidr3y,dpsidr3z,dpsidr4x,dpsidr4y,dpsidr4z;

  // Definition of cross-term dihedrals

  //         phi dihedral
  //   |--------------------|
  //   a1-----a2-----a3-----a4-----a5    cross-term atoms
  //   C      N      CA     C      N     cross-term atom types
  //          |--------------------|
  //               psi dihedral

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;

  ecmap = 0.0;
  int eflag = eflag_caller;
  ev_init(eflag,vflag);

  for (n = 0; n < ncrosstermlist; n++) {
    i1 = crosstermlist[n][0];
    i2 = crosstermlist[n][1];
    i3 = crosstermlist[n][2];
    i4 = crosstermlist[n][3];
    i5 = crosstermlist[n][4];

    type = crosstermlist[n][5];
    if (type == 0) continue;

    // calculate bond vectors for both dihedrals

    // phi
    // vb21 = r2 - r1

      vb21x = x[i2][0] - x[i1][0];
      vb21y = x[i2][1] - x[i1][1];
      vb21z = x[i2][2] - x[i1][2];
      vb12x = -1.0*vb21x;
      vb12y = -1.0*vb21y;
      vb12z = -1.0*vb21z;
      vb32x = x[i3][0] - x[i2][0];
      vb32y = x[i3][1] - x[i2][1];
      vb32z = x[i3][2] - x[i2][2];
      vb23x = -1.0*vb32x;
      vb23y = -1.0*vb32y;
      vb23z = -1.0*vb32z;

      vb34x = x[i3][0] - x[i4][0];
      vb34y = x[i3][1] - x[i4][1];
      vb34z = x[i3][2] - x[i4][2];

      // psi
      // bond vectors same as for phi: vb32

      vb43x = -1.0*vb34x;
      vb43y = -1.0*vb34y;
      vb43z = -1.0*vb34z;

      vb45x = x[i4][0] - x[i5][0];
      vb45y = x[i4][1] - x[i5][1];
      vb45z = x[i4][2] - x[i5][2];
      vb54x = -1.0*vb45x;
      vb54y = -1.0*vb45y;
      vb54z = -1.0*vb45z;

      // calculate normal vectors for planes that define the dihedral angles

      a1x = vb12y*vb23z - vb12z*vb23y;
      a1y = vb12z*vb23x - vb12x*vb23z;
      a1z = vb12x*vb23y - vb12y*vb23x;

      b1x = vb43y*vb23z - vb43z*vb23y;
      b1y = vb43z*vb23x - vb43x*vb23z;
      b1z = vb43x*vb23y - vb43y*vb23x;

      a2x = vb23y*vb34z - vb23z*vb34y;
      a2y = vb23z*vb34x - vb23x*vb34z;
      a2z = vb23x*vb34y - vb23y*vb34x;

      b2x = vb45y*vb43z - vb45z*vb43y;
      b2y = vb45z*vb43x - vb45x*vb43z;
      b2z = vb45x*vb43y - vb45y*vb43x;

      // calculate terms used later in calculations

      r32 = sqrt(vb32x*vb32x + vb32y*vb32y + vb32z*vb32z);
      a1sq = a1x*a1x + a1y*a1y + a1z*a1z;
      b1sq = b1x*b1x + b1y*b1y + b1z*b1z;

      r43 = sqrt(vb43x*vb43x + vb43y*vb43y + vb43z*vb43z);
      a2sq = a2x*a2x + a2y*a2y + a2z*a2z;
      b2sq = b2x*b2x + b2y*b2y + b2z*b2z;
      //if (a1sq<0.0001 || b1sq<0.0001 || a2sq<0.0001 || b2sq<0.0001)
      //  printf("a1sq b1sq a2sq b2sq: %f %f %f %f \n",a1sq,b1sq,a2sq,b2sq);
      if (a1sq<0.0001 || b1sq<0.0001 || a2sq<0.0001 || b2sq<0.0001) continue;
      dpr21r32 = vb21x*vb32x + vb21y*vb32y + vb21z*vb32z;
      dpr34r32 = vb34x*vb32x + vb34y*vb32y + vb34z*vb32z;
      dpr32r43 = vb32x*vb43x + vb32y*vb43y + vb32z*vb43z;
      dpr45r43 = vb45x*vb43x + vb45y*vb43y + vb45z*vb43z;

      // calculate the backbone dihedral angles as VMD and GROMACS

      phi = dihedral_angle_atan2(vb21x,vb21y,vb21z,a1x,a1y,a1z,b1x,b1y,b1z,r32);
      psi = dihedral_angle_atan2(vb32x,vb32y,vb32z,a2x,a2y,a2z,b2x,b2y,b2z,r43);

      if (phi == 180.0) phi= -180.0;
      if (psi == 180.0) psi= -180.0;

      phi1 = phi;
      if (phi1 < 0.0) phi1 += 360.0;
      psi1 = psi;
      if (psi1 < 0.0) psi1 += 360.0;

      // find the neighbor grid point index

      li1 = int(((phi1+CMAPXMIN2)/CMAPDX)+((CMAPDIM*1.0)/2.0));
      li2 = int(((psi1+CMAPXMIN2)/CMAPDX)+((CMAPDIM*1.0)/2.0));

      li3 = int((phi-CMAPXMIN2)/CMAPDX);
      li4 = int((psi-CMAPXMIN2)/CMAPDX);
      mli3 = li3 % CMAPDIM;
      mli4 = li4 % CMAPDIM;
      mli31 = (li3+1) % CMAPDIM;
      mli41 = (li4+1)  %CMAPDIM;
      mli1 = li1 % CMAPDIM;
      mli2 = li2 % CMAPDIM;
      mli11 = (li1+1) % CMAPDIM;
      mli21 = (li2+1)  %CMAPDIM;
      t1 = type-1;
      if (t1 < 0 || t1 > 5) error->all(FLERR,"Invalid CMAP crossterm_type");

      // determine the values and derivatives for the grid square points

      gs[0] = cmapgrid[t1][mli3][mli4];
      gs[1] = cmapgrid[t1][mli31][mli4];
      gs[2] = cmapgrid[t1][mli31][mli41];
      gs[3] = cmapgrid[t1][mli3][mli41];
      d1gs[0] = d1cmapgrid[t1][mli1][mli2];
      d1gs[1] = d1cmapgrid[t1][mli11][mli2];
      d1gs[2] = d1cmapgrid[t1][mli11][mli21];
      d1gs[3] = d1cmapgrid[t1][mli1][mli21];

      d2gs[0] = d2cmapgrid[t1][mli1][mli2];
      d2gs[1] = d2cmapgrid[t1][mli11][mli2];
      d2gs[2] = d2cmapgrid[t1][mli11][mli21];
      d2gs[3] = d2cmapgrid[t1][mli1][mli21];

      d12gs[0] = d12cmapgrid[t1][mli1][mli2];
      d12gs[1] = d12cmapgrid[t1][mli11][mli2];
      d12gs[2] = d12cmapgrid[t1][mli11][mli21];
      d12gs[3] = d12cmapgrid[t1][mli1][mli21];

      // calculate the cmap energy and the gradient (dE/dphi,dE/dpsi)

      bc_interpol(phi,psi,li3,li4,gs,d1gs,d2gs,d12gs);

      // sum up cmap energy contributions

      engfraction = 0.2 * E;
      if (i1 < nlocal) ecmap += engfraction;
      if (i2 < nlocal) ecmap += engfraction;
      if (i3 < nlocal) ecmap += engfraction;
      if (i4 < nlocal) ecmap += engfraction;
      if (i5 < nlocal) ecmap += engfraction;

      // calculate the derivatives dphi/dr_i

      dphidr1x = 1.0*r32/a1sq*a1x;
      dphidr1y = 1.0*r32/a1sq*a1y;
      dphidr1z = 1.0*r32/a1sq*a1z;

      dphidr2x = -1.0*r32/a1sq*a1x - dpr21r32/a1sq/r32*a1x +
        dpr34r32/b1sq/r32*b1x;
      dphidr2y = -1.0*r32/a1sq*a1y - dpr21r32/a1sq/r32*a1y +
        dpr34r32/b1sq/r32*b1y;
      dphidr2z = -1.0*r32/a1sq*a1z - dpr21r32/a1sq/r32*a1z +
        dpr34r32/b1sq/r32*b1z;

      dphidr3x = dpr34r32/b1sq/r32*b1x - dpr21r32/a1sq/r32*a1x - r32/b1sq*b1x;
      dphidr3y = dpr34r32/b1sq/r32*b1y - dpr21r32/a1sq/r32*a1y - r32/b1sq*b1y;
      dphidr3z = dpr34r32/b1sq/r32*b1z - dpr21r32/a1sq/r32*a1z - r32/b1sq*b1z;

      dphidr4x = r32/b1sq*b1x;
      dphidr4y = r32/b1sq*b1y;
      dphidr4z = r32/b1sq*b1z;

      // calculate the derivatives dpsi/dr_i

      dpsidr1x = 1.0*r43/a2sq*a2x;
      dpsidr1y = 1.0*r43/a2sq*a2y;
      dpsidr1z = 1.0*r43/a2sq*a2z;

      dpsidr2x = r43/a2sq*a2x + dpr32r43/a2sq/r43*a2x - dpr45r43/b2sq/r43*b2x;
      dpsidr2y = r43/a2sq*a2y + dpr32r43/a2sq/r43*a2y - dpr45r43/b2sq/r43*b2y;
      dpsidr2z = r43/a2sq*a2z + dpr32r43/a2sq/r43*a2z - dpr45r43/b2sq/r43*b2z;

      dpsidr3x = dpr45r43/b2sq/r43*b2x - dpr32r43/a2sq/r43*a2x - r43/b2sq*b2x;
      dpsidr3y = dpr45r43/b2sq/r43*b2y - dpr32r43/a2sq/r43*a2y - r43/b2sq*b2y;
      dpsidr3z = dpr45r43/b2sq/r43*b2z - dpr32r43/a2sq/r43*a2z - r43/b2sq*b2z;

      dpsidr4x = r43/b2sq*b2x;
      dpsidr4y = r43/b2sq*b2y;
      dpsidr4z = r43/b2sq*b2z;

      // calculate forces on cross-term atoms: F = -(dE/dPhi)*(dPhi/dr)

      f1[0] = dEdPhi*dphidr1x;
      f1[1] = dEdPhi*dphidr1y;
      f1[2] = dEdPhi*dphidr1z;
      f2[0] = dEdPhi*dphidr2x + dEdPsi*dpsidr1x;
      f2[1] = dEdPhi*dphidr2y + dEdPsi*dpsidr1y;
      f2[2] = dEdPhi*dphidr2z + dEdPsi*dpsidr1z;
      f3[0] = -dEdPhi*dphidr3x - dEdPsi*dpsidr2x;
      f3[1] = -dEdPhi*dphidr3y - dEdPsi*dpsidr2y;
      f3[2] = -dEdPhi*dphidr3z - dEdPsi*dpsidr2z;
      f4[0] = -dEdPhi*dphidr4x - dEdPsi*dpsidr3x;
      f4[1] = -dEdPhi*dphidr4y - dEdPsi*dpsidr3y;
      f4[2] = -dEdPhi*dphidr4z - dEdPsi*dpsidr3z;
      f5[0] = -dEdPsi*dpsidr4x;
      f5[1] = -dEdPsi*dpsidr4y;
      f5[2] = -dEdPsi*dpsidr4z;

      // apply force to each of the 5 atoms

      if (i1 < nlocal) {
        f[i1][0] += f1[0];
        f[i1][1] += f1[1];
        f[i1][2] += f1[2];
      }
      if (i2 < nlocal) {
        f[i2][0] += f2[0];
        f[i2][1] += f2[1];
        f[i2][2] += f2[2];
      }
      if (i3 < nlocal) {
        f[i3][0] += f3[0];
        f[i3][1] += f3[1];
        f[i3][2] += f3[2];
      }
      if (i4 < nlocal) {
        f[i4][0] += f4[0];
        f[i4][1] += f4[1];
        f[i4][2] += f4[2];
      }
      if (i5 < nlocal) {
        f[i5][0] += f5[0];
        f[i5][1] += f5[1];
        f[i5][2] += f5[2];
      }

      // tally energy and/or virial

      if (evflag) {
        nlist = 0;
        if (i1 < nlocal) list[nlist++] = i1;
        if (i2 < nlocal) list[nlist++] = i2;
        if (i3 < nlocal) list[nlist++] = i3;
        if (i4 < nlocal) list[nlist++] = i4;
        if (i5 < nlocal) list[nlist++] = i5;
        vcmap[0] = (vb12x*f1[0])+(vb32x*f3[0])+((vb43x+vb32x)*f4[0])+
          ((vb54x+vb43x+vb32x)*f5[0]);
        vcmap[1] = (vb12y*f1[1])+(vb32y*f3[1])+((vb43y+vb32y)*f4[1])+
          ((vb54y+vb43y+vb32y)*f5[1]);
        vcmap[2] = (vb12z*f1[2])+(vb32z*f3[2])+((vb43z+vb32z)*f4[2])+
          ((vb54z+vb43z+vb32z)*f5[2]);
        vcmap[3] = (vb12x*f1[1])+(vb32x*f3[1])+((vb43x+vb32x)*f4[1])+
          ((vb54x+vb43x+vb32x)*f5[1]);
        vcmap[4] = (vb12x*f1[2])+(vb32x*f3[2])+((vb43x+vb32x)*f4[2])+
          ((vb54x+vb43x+vb32x)*f5[2]);
        vcmap[5] = (vb12y*f1[2])+(vb32y*f3[2])+((vb43y+vb32y)*f4[2])+
          ((vb54y+vb43y+vb32y)*f5[2]);
        ev_tally(nlist,list,5.0,E,vcmap);
        //ev_tally(5,list,nlocal,newton_bond,E,vcmap);
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixCMAP::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCMAP::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of CMAP term
------------------------------------------------------------------------- */

double FixCMAP::compute_scalar()
{
  double all;
  MPI_Allreduce(&ecmap,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// methods to read CMAP potential file, perform interpolation
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

void FixCMAP::read_grid_map(char *cmapfile)
{
  if (comm->me == 0) {
    try {
      memset(&cmapgrid[0][0][0], 0, 6*CMAPDIM*CMAPDIM*sizeof(double));
      PotentialFileReader reader(lmp, cmapfile, "cmap grid");

      // there are six maps in this order.
      // alanine, alanine-proline, proline, proline-proline, glycine, glycine-proline.
      // read as one big blob of numbers while ignoring comments

      reader.next_dvector(&cmapgrid[0][0][0],6*CMAPDIM*CMAPDIM);

    } catch (std::exception &e) {
      error->one(FLERR,"Error reading CMAP potential file: {}", e.what());
    }
  }

  MPI_Bcast(&cmapgrid[0][0][0],6*CMAPDIM*CMAPDIM,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

void FixCMAP::spline(double *y, double *ddy, int n)
{
  // create the 2nd dervatives of a taublated function y_i(x_i)
  // at the tabulated points

  int i, j;
  double p, *u;

  memory->create(u,n-1,"cmap:u");

  ddy[0] = u[0] = 0.0;

  for (i = 1; i <= n-2; i++) {
    p = 1.0/(ddy[i-1]+4.0);
    ddy[i] = -p;
    u[i] = ((((6.0*y[i+1])-(12.0*y[i])+(6.0*y[i-1]))/(CMAPDX*CMAPDX))-u[i-1])*p;
  }

  ddy[n-1] = 0.0;

  for (j = n-2; j >= 0; j--)
    ddy[j] = ddy[j]*ddy[j+1] + u[j];

  memory->destroy(u);
}

/* ---------------------------------------------------------------------- */

void FixCMAP::spl_interpolate(double x, double *y, double *ddy, double &yo,
                              double &dyo)
{
  // perform a 1D cubic spline interpolation

  int ix;
  double a,b,a1,b1,a2,b2;

  ix = int((x-CMAPXMIN)/CMAPDX-(1./2.));

  a = (CMAPXMIN+(ix*1.0)*CMAPDX-x)/CMAPDX;
  b = (x-CMAPXMIN-(((ix-1)*1.0)*CMAPDX))/CMAPDX;

  a1 = a*a*a-a;
  b1 = b*b*b-b;

  a2 = 3.0*a*a-1.0;
  b2 = 3.0*b*b-1.0;
  yo = a*y[ix]+b*y[ix+1]+(a1*ddy[ix]+b1*ddy[ix+1])*(CMAPDX*CMAPDX)/6.0;
  dyo = (y[ix+1]-y[ix])/CMAPDX-a2/6.0*CMAPDX*ddy[ix]+b2/6.0*CMAPDX*ddy[ix+1];
}

/* ---------------------------------------------------------------------- */

void FixCMAP::set_map_derivatives(double **map, double **d1yo, double **d2yo,
                                  double **d12yo)
{
  // precompute the gradient and cross-derivatives of the map grid points.
  // use the bicubic spline to calculate the derivatives

  int i, j, k, ii, jj, xm, p;
  double phi, psi, d1y, d2y, d12y, tyyk,tdyk;
  double *tmp_y, *tmp_dy, *tmp_ddy, **tmap, **tddmap;
  int ix;
  double a,b,a1,b1,a2,b2;

  xm = CMAPDIM/2;
  p = CMAPDIM;

  d1y = 0.;
  d2y = 0.;
  d12y = 0.;

  memory->create(tmp_y,CMAPDIM*2,"cmap:tmp_y");
  memory->create(tmp_dy,CMAPDIM*2,"cmap:tmp_dy");
  memory->create(tmp_ddy,CMAPDIM*2,"cmap:tmp_ddy");
  memory->create(tmap,CMAPDIM*2,CMAPDIM*2,"cmap:tmap");
  memory->create(tddmap,CMAPDIM*2,CMAPDIM*2,"cmap:tddmap");

  // periodically expand the original map
  // use the expanded map for bicubic spline interpolation,
  //   which is used to obtain the derivatives
  // actual interpolation is done with bicubic interpolation

  for (i = 0; i < CMAPDIM*2; i++) {
    ii = ((i+CMAPDIM-xm)%CMAPDIM);
    for (j = 0; j < CMAPDIM*2; j++) {
      jj = ((j+CMAPDIM-xm)%CMAPDIM);
      tmap[i][j] = map[ii][jj];
    }
  }

  for (i = 0; i < CMAPDIM*2; i++)
    spline(tmap[i], tddmap[i], CMAPDIM*2);

  for (i = xm; i < CMAPDIM+xm; i++) {
    phi = (i-xm)*CMAPDX-180.0;
    for (j = xm; j < CMAPDIM+xm; j++) {
      psi = (j-xm)*CMAPDX-180.0;
      ix = int((psi-CMAPXMIN)/CMAPDX);
      a = (CMAPXMIN+((ix+1)*1.0)*CMAPDX-psi)/CMAPDX;
      b = (psi-CMAPXMIN-((ix)*1.0)*CMAPDX)/CMAPDX;
      a1 = a*a*a-a;
      b1 = b*b*b-b;
      a2 = 3.0*a*a-1.0;
      b2 = 3.0*b*b-1.0;
      for (k = 0; k < CMAPDIM*2; k++) {
        tyyk = tmp_y[k];
        tdyk = tmp_dy[k];
        tyyk = a*tmap[k][ix]+b*tmap[k][ix+1]+
          (a1*tddmap[k][ix]+b1*tddmap[k][ix+1])*(CMAPDX*CMAPDX)/6.0;
        tdyk = (tmap[k][ix+1]-tmap[k][ix])/CMAPDX-
          (a2/6.0*CMAPDX*tddmap[k][ix])+(b2/6.0*CMAPDX*tddmap[k][ix+1]);
        tmp_y[k] = tyyk;
        tmp_dy[k] = tdyk;
      }

      spline(tmp_y,tmp_ddy,CMAPDIM+xm+xm);
      ix = int((phi-CMAPXMIN)/CMAPDX);
      a = (CMAPXMIN+((ix+1)*1.0)*CMAPDX-phi)/CMAPDX;
      b = (phi-CMAPXMIN-(ix*1.0)*CMAPDX)/CMAPDX;
      a1 = a*a*a-a;
      b1 = b*b*b-b;
      a2 = 3.0*a*a-1.0;
      b2 = 3.0*b*b-1.0;
      d1y = (tmp_y[ix+1]-tmp_y[ix])/CMAPDX-
        a2/6.0*CMAPDX*tmp_ddy[ix]+b2/6.0*CMAPDX*tmp_ddy[ix+1];
      spline(tmp_dy,tmp_ddy,CMAPDIM+xm+xm);
      ix = int((phi-CMAPXMIN)/CMAPDX);
      a = (CMAPXMIN+((ix+1)*1.0)*CMAPDX-phi)/CMAPDX;
      b = (phi-CMAPXMIN-(ix*1.0)*CMAPDX)/CMAPDX;
      a1 = a*a*a-a;
      b1 = b*b*b-b;
      a2 = 3.0*a*a-1.0;
      b2 = 3.0*b*b-1.0;
      d2y = a*tmp_dy[ix]+b*tmp_dy[ix+1]+
        (a1*tmp_ddy[ix]+b1*tmp_ddy[ix+1])*(CMAPDX*CMAPDX)/6.0;
      d12y = (tmp_dy[ix+1]-tmp_dy[ix])/CMAPDX-
        a2/6.0*CMAPDX*tmp_ddy[ix]+b2/6.0*CMAPDX*tmp_ddy[ix+1];
      d1yo[i%p][j%p] = d1y;
      d2yo[i%p][j%p] = d2y;
      d12yo[i%p][j%p] = d12y;
    }
  }

  memory->destroy(tmp_y);
  memory->destroy(tmp_dy);
  memory->destroy(tmp_ddy);
  memory->destroy(tmap);
  memory->destroy(tddmap);
}

/* ---------------------------------------------------------------------- */

double FixCMAP::dihedral_angle_atan2(double fx, double fy, double fz,
                                      double ax, double ay, double az,
                                      double bx, double by, double bz,
                                      double absg)
{
  // calculate the dihedral angle

  double angle = 0.0, arg1, arg2;

  arg1 = absg*(fx*bx+fy*by+fz*bz);
  arg2 = ax*bx+ay*by+az*bz;

  if (arg1 == 0 && arg2 == 0)
    error->all(FLERR,"CMAP: atan2 function cannot take 2 zero arguments");
  else {
    angle = atan2(arg1,arg2);
    angle = angle*180.0/MY_PI;
  }

  return angle;
}

/* ---------------------------------------------------------------------- */

void FixCMAP::bc_coeff(double *gs, double *d1gs, double *d2gs, double *d12gs)
{
  // calculate the bicubic interpolation coefficients c_ij

  static int wt[16][16] =
    { {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
      {-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0},
      {2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
      {0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1},
      {0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1},
      {-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0},
      {9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2},
      {-6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2},
      {2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0},
      {-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1},
      {4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1}
    };

  int i, j, k, in;
  double xx, x[16];

  for (i = 0; i < 4; i++) {
    x[i] = gs[i];
    x[i+4] = d1gs[i]*CMAPDX;
    x[i+8] = d2gs[i]*CMAPDX;
    x[i+12] = d12gs[i]*CMAPDX*CMAPDX;
  }

  in = 0;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      xx = 0.0;
      for (k = 0; k < 16; k++) xx += wt[in][k]*x[k];
      in++;
      cij[i][j] = xx;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixCMAP::bc_interpol(double x1, double x2, int low1, int low2, double *gs,
                           double *d1gs, double *d2gs, double *d12gs)
{
  // for a given point of interest and its corresponding grid square values,
  //   gradients and cross-derivatives
  // calculate the interpolated value of the point of interest (POI)

  int i;
  double t, u, gs1l, gs2l;

  // set the interpolation coefficients

  bc_coeff(gs,d1gs,d2gs,d12gs);

  gs1l = g_axis[low1];
  gs2l = g_axis[low2];

  t = (x1-gs1l)/CMAPDX;
  u = (x2-gs2l)/CMAPDX;

  E = dEdPhi = dEdPsi = 0.0;

  for (i = 3; i >= 0; i--) {
    E = t*E + ((cij[i][3]*u+cij[i][2])*u+cij[i][1])*u+cij[i][0];
    dEdPhi = u*dEdPhi + (3.0*cij[3][i]*t+2.0*cij[2][i])*t+cij[1][i];
    dEdPsi = t*dEdPsi + (3.0*cij[i][3]*u+2.0*cij[i][2])*u+cij[i][1];
  }

  dEdPhi *= (180.0/MY_PI/CMAPDX);
  dEdPsi *= (180.0/MY_PI/CMAPDX);
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// methods to read and write data file
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

void FixCMAP::read_data_header(char *line)
{
  ValueTokenizer values(line);

  try {
    ncmap = values.next_bigint();
    if (values.count() == 2) {
      if (values.next_string() != "crossterms")
        throw TokenizerException("invalid format",utils::trim(line));
    } else if (values.count() == 3) {
      if ((values.next_string() != "cmap") || (values.next_string() != "crossterms"))
        throw TokenizerException("invalid format",utils::trim(line));
    } else {
      throw TokenizerException("valid format",utils::trim(line));
    }
  } catch (std::exception &e) {
    error->all(FLERR,"Invalid read data header line for fix cmap: {}", e.what());
  }

  // not set in constructor because this fix could be defined before newton command

  newton_bond = force->newton_bond;
}

/* ----------------------------------------------------------------------
   unpack N lines in buf from section of data file labeled by keyword
   id_offset is applied to atomID fields if multiple data files are read
   store CMAP interactions as if newton_bond = OFF, even if actually ON
------------------------------------------------------------------------- */

void FixCMAP::read_data_section(char * /*keyword*/, int /*n*/, char *buf,
                                 tagint id_offset)
{
  int m,itype;
  tagint atom1,atom2,atom3,atom4,atom5;

  auto lines = utils::split_lines(buf);
  if (lines.size() == 0) return;

  // loop over lines of CMAP crossterms
  // tokenize the line into values
  // add crossterm to one of my atoms, depending on newton_bond

  for (const auto &line : lines) {
    ValueTokenizer values(line);
    try {
      values.skip();
      itype = values.next_int();
      atom1 = values.next_tagint();
      atom2 = values.next_tagint();
      atom3 = values.next_tagint();
      atom4 = values.next_tagint();
      atom5 = values.next_tagint();
      if (values.has_next()) throw TokenizerException("too many items",line);
    } catch (std::exception &e) {
      error->all(FLERR,"Incorrect format of CMAP section: {}", e.what());
    }

    atom1 += id_offset;
    atom2 += id_offset;
    atom3 += id_offset;
    atom4 += id_offset;
    atom5 += id_offset;

    if ((m = atom->map(atom1)) >= 0) {
      if (num_crossterm[m] == CMAPMAX) error->one(FLERR,"Too many CMAP crossterms for one atom");
      crossterm_type[m][num_crossterm[m]] = itype;
      crossterm_atom1[m][num_crossterm[m]] = atom1;
      crossterm_atom2[m][num_crossterm[m]] = atom2;
      crossterm_atom3[m][num_crossterm[m]] = atom3;
      crossterm_atom4[m][num_crossterm[m]] = atom4;
      crossterm_atom5[m][num_crossterm[m]] = atom5;
      num_crossterm[m]++;
    }

    if ((m = atom->map(atom2)) >= 0) {
      if (num_crossterm[m] == CMAPMAX) error->one(FLERR,"Too many CMAP crossterms for one atom");
      crossterm_type[m][num_crossterm[m]] = itype;
      crossterm_atom1[m][num_crossterm[m]] = atom1;
      crossterm_atom2[m][num_crossterm[m]] = atom2;
      crossterm_atom3[m][num_crossterm[m]] = atom3;
      crossterm_atom4[m][num_crossterm[m]] = atom4;
      crossterm_atom5[m][num_crossterm[m]] = atom5;
      num_crossterm[m]++;
    }

    if ((m = atom->map(atom3)) >= 0) {
      if (num_crossterm[m] == CMAPMAX) error->one(FLERR,"Too many CMAP crossterms for one atom");
      crossterm_type[m][num_crossterm[m]] = itype;
      crossterm_atom1[m][num_crossterm[m]] = atom1;
      crossterm_atom2[m][num_crossterm[m]] = atom2;
      crossterm_atom3[m][num_crossterm[m]] = atom3;
      crossterm_atom4[m][num_crossterm[m]] = atom4;
      crossterm_atom5[m][num_crossterm[m]] = atom5;
      num_crossterm[m]++;
    }

    if ((m = atom->map(atom4)) >= 0) {
      if (num_crossterm[m] == CMAPMAX) error->one(FLERR,"Too many CMAP crossterms for one atom");
      crossterm_type[m][num_crossterm[m]] = itype;
      crossterm_atom1[m][num_crossterm[m]] = atom1;
      crossterm_atom2[m][num_crossterm[m]] = atom2;
      crossterm_atom3[m][num_crossterm[m]] = atom3;
      crossterm_atom4[m][num_crossterm[m]] = atom4;
      crossterm_atom5[m][num_crossterm[m]] = atom5;
      num_crossterm[m]++;
    }

    if ((m = atom->map(atom5)) >= 0) {
      if (num_crossterm[m] == CMAPMAX) error->one(FLERR,"Too many CMAP crossterms for one atom");
      crossterm_type[m][num_crossterm[m]] = itype;
      crossterm_atom1[m][num_crossterm[m]] = atom1;
      crossterm_atom2[m][num_crossterm[m]] = atom2;
      crossterm_atom3[m][num_crossterm[m]] = atom3;
      crossterm_atom4[m][num_crossterm[m]] = atom4;
      crossterm_atom5[m][num_crossterm[m]] = atom5;
      num_crossterm[m]++;
    }
  }
}

/* ---------------------------------------------------------------------- */

bigint FixCMAP::read_data_skip_lines(char * /*keyword*/)
{
  return ncmap;
}

/* ----------------------------------------------------------------------
   write Mth header line to file
   only called by proc 0
------------------------------------------------------------------------- */

void FixCMAP::write_data_header(FILE *fp, int /*mth*/)
{
  fmt::print(fp,"{} crossterms\n",ncmap);
}

/* ----------------------------------------------------------------------
   return size I own for Mth data section
   # of data sections = 1 for this fix
   nx = # of crossterms owned by my local atoms
     if newton_bond off, atom only owns crossterm if it is atom3
   ny = columns = type + 5 atom IDs
------------------------------------------------------------------------- */

void FixCMAP::write_data_section_size(int /*mth*/, int &nx, int &ny)
{
  int i,m;

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  nx = 0;
  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_crossterm[i]; m++)
      if (crossterm_atom3[i][m] == tag[i]) nx++;

  ny = 6;
}

/* ----------------------------------------------------------------------
   pack values for Mth data section into 2d buf
   buf allocated by caller as owned crossterms by 6
------------------------------------------------------------------------- */

void FixCMAP::write_data_section_pack(int /*mth*/, double **buf)
{
  int i,m;

  // 1st column = CMAP type
  // 2nd-6th columns = 5 atom IDs

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int n = 0;
  for (i = 0; i < nlocal; i++) {
    for (m = 0; m < num_crossterm[i]; m++) {
      if (crossterm_atom3[i][m] != tag[i]) continue;
      buf[n][0] = ubuf(crossterm_type[i][m]).d;
      buf[n][1] = ubuf(crossterm_atom1[i][m]).d;
      buf[n][2] = ubuf(crossterm_atom2[i][m]).d;
      buf[n][3] = ubuf(crossterm_atom3[i][m]).d;
      buf[n][4] = ubuf(crossterm_atom4[i][m]).d;
      buf[n][5] = ubuf(crossterm_atom5[i][m]).d;
      n++;
    }
  }
}

/* ----------------------------------------------------------------------
   write section keyword for Mth data section to file
   use Molecules or Charges if that is only field, else use fix ID
   only called by proc 0
------------------------------------------------------------------------- */

void FixCMAP::write_data_section_keyword(int /*mth*/, FILE *fp)
{
  fprintf(fp,"\nCMAP\n\n");
}

/* ----------------------------------------------------------------------
   write N lines from buf to file
   convert buf fields to int or double depending on styles
   index can be used to prepend global numbering
   only called by proc 0
------------------------------------------------------------------------- */

void FixCMAP::write_data_section(int /*mth*/, FILE *fp,
                                  int n, double **buf, int index)
{
  for (int i = 0; i < n; i++)
    fmt::print(fp,"{} {} {} {} {} {} {}\n",
               index+i,ubuf(buf[i][0]).i, ubuf(buf[i][1]).i, ubuf(buf[i][2]).i,
               ubuf(buf[i][3]).i,ubuf(buf[i][4]).i,ubuf(buf[i][5]).i);
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// methods for restart and communication
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixCMAP::write_restart(FILE *fp)
{
  if (comm->me == 0) {
    int size = sizeof(bigint);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(&ncmap,sizeof(bigint),1,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixCMAP::restart(char *buf)
{
  ncmap = *((bigint *) buf);
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixCMAP::pack_restart(int i, double *buf)
{
  int n = 1;
  for (int m = 0; m < num_crossterm[i]; m++) {
    buf[n++] = ubuf(MAX(crossterm_type[i][m],-crossterm_type[i][m])).d;
    buf[n++] = ubuf(crossterm_atom1[i][m]).d;
    buf[n++] = ubuf(crossterm_atom2[i][m]).d;
    buf[n++] = ubuf(crossterm_atom3[i][m]).d;
    buf[n++] = ubuf(crossterm_atom4[i][m]).d;
    buf[n++] = ubuf(crossterm_atom5[i][m]).d;
  }
  // pack buf[0] this way because other fixes unpack it
  buf[0] = n;

  return n;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixCMAP::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

   int n = 0;
   for (int i = 0; i < nth; i++) n += static_cast<int> (extra[nlocal][n]);

   int count = static_cast<int> (extra[nlocal][n++]);
   num_crossterm[nlocal] = (count-1)/6;

   for (int m = 0; m < num_crossterm[nlocal]; m++) {
     crossterm_type[nlocal][m] = (int) ubuf(extra[nlocal][n++]).i;
     crossterm_atom1[nlocal][m] = (tagint) ubuf(extra[nlocal][n++]).i;
     crossterm_atom2[nlocal][m] = (tagint) ubuf(extra[nlocal][n++]).i;
     crossterm_atom3[nlocal][m] = (tagint) ubuf(extra[nlocal][n++]).i;
     crossterm_atom4[nlocal][m] = (tagint) ubuf(extra[nlocal][n++]).i;
     crossterm_atom5[nlocal][m] = (tagint) ubuf(extra[nlocal][n++]).i;
   }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixCMAP::maxsize_restart()
{
  return 1 + CMAPMAX*6;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixCMAP::size_restart(int nlocal)
{
  return 1 + num_crossterm[nlocal]*6;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixCMAP::grow_arrays(int nmax)
{
  num_crossterm = memory->grow(num_crossterm,nmax,"cmap:num_crossterm");
  crossterm_type = memory->grow(crossterm_type,nmax,CMAPMAX,
                                "cmap:crossterm_type");
  crossterm_atom1 = memory->grow(crossterm_atom1,nmax,CMAPMAX,
                                 "cmap:crossterm_atom1");
  crossterm_atom2 = memory->grow(crossterm_atom2,nmax,CMAPMAX,
                                 "cmap:crossterm_atom2");
  crossterm_atom3 = memory->grow(crossterm_atom3,nmax,CMAPMAX,
                                 "cmap:crossterm_atom3");
  crossterm_atom4 = memory->grow(crossterm_atom4,nmax,CMAPMAX,
                                 "cmap:crossterm_atom4");
  crossterm_atom5 = memory->grow(crossterm_atom5,nmax,CMAPMAX,
                                 "cmap:crossterm_atom5");

  // must initialize num_crossterm to 0 for added atoms
  // may never be set for some atoms when data file is read

  for (int i = nmax_previous; i < nmax; i++) num_crossterm[i] = 0;
  nmax_previous = nmax;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixCMAP::copy_arrays(int i, int j, int /*delflag*/)
{
  num_crossterm[j] = num_crossterm[i];

  for (int k = 0; k < num_crossterm[j]; k++) {
    crossterm_type[j][k] = crossterm_type[i][k];
    crossterm_atom1[j][k] = crossterm_atom1[i][k];
    crossterm_atom2[j][k] = crossterm_atom2[i][k];
    crossterm_atom3[j][k] = crossterm_atom3[i][k];
    crossterm_atom4[j][k] = crossterm_atom4[i][k];
    crossterm_atom5[j][k] = crossterm_atom5[i][k];
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixCMAP::set_arrays(int i)
{
  num_crossterm[i] = 0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixCMAP::pack_exchange(int i, double *buf)
{
  int n = 0;
  buf[n++] = ubuf(num_crossterm[i]).d;
  for (int m = 0; m < num_crossterm[i]; m++) {
    buf[n++] = ubuf(crossterm_type[i][m]).d;
    buf[n++] = ubuf(crossterm_atom1[i][m]).d;
    buf[n++] = ubuf(crossterm_atom2[i][m]).d;
    buf[n++] = ubuf(crossterm_atom3[i][m]).d;
    buf[n++] = ubuf(crossterm_atom4[i][m]).d;
    buf[n++] = ubuf(crossterm_atom5[i][m]).d;
  }
  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixCMAP::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  num_crossterm[nlocal] = (int) ubuf(buf[n++]).i;
  for (int m = 0; m < num_crossterm[nlocal]; m++) {
    crossterm_type[nlocal][m] = (int) ubuf(buf[n++]).i;
    crossterm_atom1[nlocal][m] = (tagint) ubuf(buf[n++]).i;
    crossterm_atom2[nlocal][m] = (tagint) ubuf(buf[n++]).i;
    crossterm_atom3[nlocal][m] = (tagint) ubuf(buf[n++]).i;
    crossterm_atom4[nlocal][m] = (tagint) ubuf(buf[n++]).i;
    crossterm_atom5[nlocal][m] = (tagint) ubuf(buf[n++]).i;
  }
  return n;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixCMAP::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = (double)nmax * sizeof(int);        // num_crossterm
  bytes += (double)nmax*CMAPMAX * sizeof(int);      // crossterm_type
  bytes += (double)5*nmax*CMAPMAX * sizeof(int);    // crossterm_atom 12345
  bytes += (double)maxcrossterm*6 * sizeof(int);    // crosstermlist
  return bytes;
}
