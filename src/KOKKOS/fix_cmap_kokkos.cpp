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

/* ----------------------------------------------------------------------
   Contributing author: Mitch Murphy (alphataubio@gmail.com)
------------------------------------------------------------------------- */

#include "fix_cmap_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <iostream>

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr int LISTDELTA = 10000;
static constexpr double LB_FACTOR = 1.5;

static constexpr int CMAPMAX = 6;   // max # of CMAP terms stored by one atom
static constexpr int CMAPDIM = 24;  // grid map dimension is 24 x 24
static constexpr double CMAPXMIN2 = -180.0;
static constexpr double CMAPDX = 15.0; // 360/CMAPDIM


/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixCMAPKokkos<DeviceType>::FixCMAPKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixCMAP(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | F_MASK;
  datamask_modify = F_MASK;

  // allocate memory for CMAP data
  memoryKK->create_kokkos(k_g_axis,g_axis,CMAPDIM,"cmap:g_axis");
  memoryKK->create_kokkos(k_cmapgrid,cmapgrid,CMAPMAX,CMAPDIM,CMAPDIM,"cmap:grid");
  memoryKK->create_kokkos(k_d1cmapgrid,d1cmapgrid,CMAPMAX,CMAPDIM,CMAPDIM,"cmap:d1grid");
  memoryKK->create_kokkos(k_d2cmapgrid,d2cmapgrid,CMAPMAX,CMAPDIM,CMAPDIM,"cmap:d2grid");
  memoryKK->create_kokkos(k_d12cmapgrid,d12cmapgrid,CMAPMAX,CMAPDIM,CMAPDIM,"cmap:d12grid");

  d_g_axis = k_g_axis.template view<DeviceType>();
  d_cmapgrid = k_cmapgrid.template view<DeviceType>();
  d_d1cmapgrid = k_d1cmapgrid.template view<DeviceType>();
  d_d2cmapgrid = k_d2cmapgrid.template view<DeviceType>();
  d_d12cmapgrid = k_d12cmapgrid.template view<DeviceType>();

  // read and setup CMAP data
  read_grid_map(arg[3]);

  int i = 0;
  double angle = -180.0;

  while (angle < 180.0) {
    g_axis[i] = angle;
    angle += CMAPDX;
    i++;
  }

  FixCMAPKokkos::grow_arrays(atom->nmax);

  for( int i=0 ; i<CMAPDIM ; i++ )
    k_g_axis.h_view(i) = g_axis[i];

  for( int i=0 ; i<CMAPMAX ; i++ ) {

    // pre-compute the derivatives of the maps
    set_map_derivatives(cmapgrid[i],d1cmapgrid[i],d2cmapgrid[i],d12cmapgrid[i]);

    for( int j=0 ; j<CMAPDIM ; j++ ) {
      for( int k=0 ; k<CMAPDIM ; k++ ) {
        k_cmapgrid.h_view(i,j,k) = cmapgrid[i][j][k];
        k_d1cmapgrid.h_view(i,j,k) = d1cmapgrid[i][j][k];
        k_d2cmapgrid.h_view(i,j,k) = d2cmapgrid[i][j][k];
        k_d12cmapgrid.h_view(i,j,k) = d12cmapgrid[i][j][k];

      }
    }
  }

  k_g_axis.modify_host();
  k_cmapgrid.modify_host();
  k_d1cmapgrid.modify_host();
  k_d2cmapgrid.modify_host();
  k_d12cmapgrid.modify_host();

  k_g_axis.template sync<DeviceType>();
  k_cmapgrid.template sync<DeviceType>();
  k_d1cmapgrid.template sync<DeviceType>();
  k_d2cmapgrid.template sync<DeviceType>();
  k_d12cmapgrid.template sync<DeviceType>();

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixCMAPKokkos<DeviceType>::~FixCMAPKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_g_axis,g_axis);
  memoryKK->destroy_kokkos(k_cmapgrid,cmapgrid);
  memoryKK->destroy_kokkos(k_d1cmapgrid,d1cmapgrid);
  memoryKK->destroy_kokkos(k_d2cmapgrid,d2cmapgrid);
  memoryKK->destroy_kokkos(k_d12cmapgrid,d12cmapgrid);

  memoryKK->destroy_kokkos(k_crosstermlist,crosstermlist);

  memoryKK->destroy_kokkos(k_num_crossterm,num_crossterm);
  memoryKK->destroy_kokkos(k_crossterm_type,crossterm_type);
  memoryKK->destroy_kokkos(k_crossterm_atom1,crossterm_atom1);
  memoryKK->destroy_kokkos(k_crossterm_atom2,crossterm_atom2);
  memoryKK->destroy_kokkos(k_crossterm_atom3,crossterm_atom3);
  memoryKK->destroy_kokkos(k_crossterm_atom4,crossterm_atom4);
  memoryKK->destroy_kokkos(k_crossterm_atom5,crossterm_atom5);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixCMAPKokkos<DeviceType>::init()
{
  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot yet use respa with Kokkos");
}

/* ----------------------------------------------------------------------
   store local neighbor list as if newton_bond = OFF, even if actually ON
------------------------------------------------------------------------- */

template<class DeviceType>
void FixCMAPKokkos<DeviceType>::pre_neighbor()
{
  int i,m,atom1,atom2,atom3,atom4,atom5;
  const int me = comm->me;
  const int nprocs = comm->nprocs;

  // guesstimate initial length of local crossterm list
  // if ncmap was not set (due to read_restart, no read_data),
  //   then list will grow by LISTDELTA chunks

  if (maxcrossterm == 0) {
    if (nprocs == 1) maxcrossterm = ncmap;
    else maxcrossterm = static_cast<int> (LB_FACTOR*ncmap/nprocs);
    memoryKK->create_kokkos(k_crosstermlist,crosstermlist,maxcrossterm,CMAPMAX,"cmap:crosstermlist");
    d_crosstermlist = k_crosstermlist.template view<DeviceType>();
  }

  atomKK->sync(execution_space,X_MASK);
  d_x = atomKK->k_x.view<DeviceType>();
  int nlocal = atomKK->nlocal;

  map_style = atom->map_style;
  if (map_style == Atom::MAP_ARRAY) {
    k_map_array = atomKK->k_map_array;
    k_map_array.template sync<DeviceType>();
  } else if (map_style == Atom::MAP_HASH) {
    k_map_hash = atomKK->k_map_hash;
    k_map_hash.template sync<DeviceType>();
  }

  atomKK->k_sametag.sync<DeviceType>();
  d_sametag = atomKK->k_sametag.view<DeviceType>();

  ncrosstermlist = 0;

  for (i = 0; i < nlocal; i++) {
    for (m = 0; m < k_num_crossterm.h_view(i); m++) {

      atom1 = AtomKokkos::map_kokkos<DeviceType>(k_crossterm_atom1.h_view(i,m),map_style,k_map_array,k_map_hash);
      atom2 = AtomKokkos::map_kokkos<DeviceType>(k_crossterm_atom2.h_view(i,m),map_style,k_map_array,k_map_hash);
      atom3 = AtomKokkos::map_kokkos<DeviceType>(k_crossterm_atom3.h_view(i,m),map_style,k_map_array,k_map_hash);
      atom4 = AtomKokkos::map_kokkos<DeviceType>(k_crossterm_atom4.h_view(i,m),map_style,k_map_array,k_map_hash);
      atom5 = AtomKokkos::map_kokkos<DeviceType>(k_crossterm_atom5.h_view(i,m),map_style,k_map_array,k_map_hash);

      if (atom1 == -1 || atom2 == -1 || atom3 == -1 ||
          atom4 == -1 || atom5 == -1)
        error->one(FLERR,"CMAP atoms {} {} {} {} {} missing on "
                                     "proc {} at step {}",
                                     d_crossterm_atom1(i,m),d_crossterm_atom2(i,m),
                                     d_crossterm_atom3(i,m),d_crossterm_atom4(i,m),
                                     d_crossterm_atom5(i,m),me,update->ntimestep);
      atom1 = closest_image(i,atom1);
      atom2 = closest_image(i,atom2);
      atom3 = closest_image(i,atom3);
      atom4 = closest_image(i,atom4);
      atom5 = closest_image(i,atom5);

      if (i <= atom1 && i <= atom2 && i <= atom3 &&
          i <= atom4 && i <= atom5) {
        if (ncrosstermlist == maxcrossterm) {
          maxcrossterm += LISTDELTA;
          memoryKK->grow_kokkos(k_crosstermlist,crosstermlist,maxcrossterm,CMAPMAX,"cmap:crosstermlist");

          d_crosstermlist = k_crosstermlist.template view<DeviceType>();
        }
        d_crosstermlist(ncrosstermlist,0) = atom1;
        d_crosstermlist(ncrosstermlist,1) = atom2;
        d_crosstermlist(ncrosstermlist,2) = atom3;
        d_crosstermlist(ncrosstermlist,3) = atom4;
        d_crosstermlist(ncrosstermlist,4) = atom5;
        d_crosstermlist(ncrosstermlist,5) = d_crossterm_type(i,m);
        ncrosstermlist++;
      }
    }
  }
}


/* ----------------------------------------------------------------------
   compute CMAP terms as if newton_bond = OFF, even if actually ON
------------------------------------------------------------------------- */

template<class DeviceType>
void FixCMAPKokkos<DeviceType>::post_force(int vflag)
{

  d_x = atomKK->k_x.template view<DeviceType>();
  d_f = atomKK->k_f.template view<DeviceType>();
  atomKK->sync(execution_space,X_MASK|F_MASK);
  k_crosstermlist.template sync<DeviceType>();

  ecmap = 0.0;
  int eflag = eflag_caller;
  ev_init(eflag,vflag);

  copymode = 1;
  Kokkos::parallel_for(ncrosstermlist, *this);
  copymode = 0;
  atomKK->modified(execution_space,F_MASK);

       std::cerr << fmt::format("*** post_force ncrosstermlist {} vflag {} ecmap {}\n",ncrosstermlist,vflag,ecmap);

}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixCMAPKokkos<DeviceType>::operator()(const int n) const
{

  int i1,i2,i3,i4,i5,type;
  int li1, li2, mli1,mli2,mli11,mli21,t1,li3,li4,mli3,mli4,mli31,mli41;

  // vectors needed to calculate the cross-term dihedral angles
  double vb21x,vb21y,vb21z,vb32x,vb32y,vb32z,vb34x,vb34y,vb34z;
  double vb23x,vb23y,vb23z;
  double vb43x,vb43y,vb43z,vb45x,vb45y,vb45z,a1x,a1y,a1z,b1x,b1y,b1z;
  double a2x,a2y,a2z,b2x,b2y,b2z,r32,a1sq,b1sq,a2sq,b2sq,dpr21r32,dpr34r32;
  double dpr32r43,dpr45r43,r43,vb12x,vb12y,vb12z;
  // cross-term dihedral angles
  double phi,psi,phi1,psi1;
  double f1[3],f2[3],f3[3],f4[3],f5[3];
  double gs[4],d1gs[4],d2gs[4],d12gs[4];

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

  int nlocal = atomKK->nlocal;

  i1 = d_crosstermlist(n,0);
  i2 = d_crosstermlist(n,1);
  i3 = d_crosstermlist(n,2);
  i4 = d_crosstermlist(n,3);
  i5 = d_crosstermlist(n,4);
  type = d_crosstermlist(n,5);
  if (type == 0) return;

  // calculate bond vectors for both dihedrals

  // phi
  // vb21 = r2 - r1

  vb21x = d_x(i2,0) - d_x(i1,0);
  vb21y = d_x(i2,1) - d_x(i1,1);
  vb21z = d_x(i2,2) - d_x(i1,2);
  vb12x = -1.0*vb21x;
  vb12y = -1.0*vb21y;
  vb12z = -1.0*vb21z;
  vb32x = d_x(i3,0) - d_x(i2,0);
  vb32y = d_x(i3,1) - d_x(i2,1);
  vb32z = d_x(i3,2) - d_x(i2,2);
  vb23x = -1.0*vb32x;
  vb23y = -1.0*vb32y;
  vb23z = -1.0*vb32z;

  vb34x = d_x(i3,0) - d_x(i4,0);
  vb34y = d_x(i3,1) - d_x(i4,1);
  vb34z = d_x(i3,2) - d_x(i4,2);

  // psi
  // bond vectors same as for phi: vb32

  vb43x = -1.0*vb34x;
  vb43y = -1.0*vb34y;
  vb43z = -1.0*vb34z;

  vb45x = d_x(i4,0) - d_x(i5,0);
  vb45y = d_x(i4,1) - d_x(i5,1);
  vb45z = d_x(i4,2) - d_x(i5,2);

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
  if (a1sq<0.0001 || b1sq<0.0001 || a2sq<0.0001 || b2sq<0.0001) return;
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
  if (t1 < 0 || t1 > 5) Kokkos::abort("Invalid CMAP crossterm_type");

      // determine the values and derivatives for the grid square points

      gs[0] = d_cmapgrid(t1,mli3,mli4);
      gs[1] = d_cmapgrid(t1,mli31,mli4);
      gs[2] = d_cmapgrid(t1,mli31,mli41);
      gs[3] = d_cmapgrid(t1,mli3,mli41);
      d1gs[0] = d_d1cmapgrid(t1,mli1,mli2);
      d1gs[1] = d_d1cmapgrid(t1,mli11,mli2);
      d1gs[2] = d_d1cmapgrid(t1,mli11,mli21);
      d1gs[3] = d_d1cmapgrid(t1,mli1,mli21);
      d2gs[0] = d_d2cmapgrid(t1,mli1,mli2);
      d2gs[1] = d_d2cmapgrid(t1,mli11,mli2);
      d2gs[2] = d_d2cmapgrid(t1,mli11,mli21);
      d2gs[3] = d_d2cmapgrid(t1,mli1,mli21);
      d12gs[0] = d_d12cmapgrid(t1,mli1,mli2);
      d12gs[1] = d_d12cmapgrid(t1,mli11,mli2);
      d12gs[2] = d_d12cmapgrid(t1,mli11,mli21);
      d12gs[3] = d_d12cmapgrid(t1,mli1,mli21);

      // calculate the cmap energy and the gradient (dE/dphi,dE/dpsi)

      double E, dEdPhi, dEdPsi;
      bc_interpol(phi,psi,li3,li4,gs,d1gs,d2gs,d12gs,E,dEdPhi,dEdPsi);

      // sum up cmap energy contributions

      double ecmapKK = 0.0;

// FIXME: needed for compute_scalar()
      double engfraction = 0.2 * E;
      if (i1 < nlocal) ecmapKK += engfraction;
      if (i2 < nlocal) ecmapKK += engfraction;
      if (i3 < nlocal) ecmapKK += engfraction;
      if (i4 < nlocal) ecmapKK += engfraction;
      if (i5 < nlocal) ecmapKK += engfraction;

      //std::cerr << fmt::format("*** i {} {} {} {} {} nlocal {} E {} ecmapKK {}\n",
        //i1,i2,i3,i4,i5,nlocal,E,ecmapKK);

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
        d_f(i1,0) += f1[0];
        d_f(i1,1) += f1[1];
        d_f(i1,2) += f1[2];
      }
      if (i2 < nlocal) {
        d_f(i2,0) += f2[0];
        d_f(i2,1) += f2[1];
        d_f(i2,2) += f2[2];
      }
      if (i3 < nlocal) {
        d_f(i3,0) += f3[0];
        d_f(i3,1) += f3[1];
        d_f(i3,2) += f3[2];
      }
      if (i4 < nlocal) {
        d_f(i4,0) += f4[0];
        d_f(i4,1) += f4[1];
        d_f(i4,2) += f4[2];
      }
      if (i5 < nlocal) {
        d_f(i5,0) += f5[0];
        d_f(i5,1) += f5[1];
        d_f(i5,2) += f5[2];
      }

      // tally energy and/or virial

/*
      if (evflag) {
        //std::cerr << "******** tally energy and/or virial\n";
        int nlist = 0;
        int list[5];
        double vb54x = -1.0*vb45x;
        double vb54y = -1.0*vb45y;
        double vb54z = -1.0*vb45z;
        double vcmap[CMAPMAX];

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
*/

  }

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

template<class DeviceType>
void FixCMAPKokkos<DeviceType>::grow_arrays(int nmax)
{
  k_num_crossterm.template sync<LMPHostType>();
  k_crossterm_type.template sync<LMPHostType>();
  k_crossterm_atom1.template sync<LMPHostType>();
  k_crossterm_atom2.template sync<LMPHostType>();
  k_crossterm_atom3.template sync<LMPHostType>();
  k_crossterm_atom4.template sync<LMPHostType>();
  k_crossterm_atom5.template sync<LMPHostType>();

  // force reallocation on host
  k_num_crossterm.template modify<LMPHostType>();
  k_crossterm_type.template modify<LMPHostType>();
  k_crossterm_atom1.template modify<LMPHostType>();
  k_crossterm_atom2.template modify<LMPHostType>();
  k_crossterm_atom3.template modify<LMPHostType>();
  k_crossterm_atom4.template modify<LMPHostType>();
  k_crossterm_atom5.template modify<LMPHostType>();

  memoryKK->grow_kokkos(k_num_crossterm,num_crossterm,nmax,"cmap:num_crossterm");
  memoryKK->grow_kokkos(k_crossterm_type,crossterm_type,nmax,CMAPMAX,"cmap:crossterm_type");
  memoryKK->grow_kokkos(k_crossterm_atom1,crossterm_atom1,nmax,CMAPMAX,"cmap:crossterm_atom1");
  memoryKK->grow_kokkos(k_crossterm_atom2,crossterm_atom2,nmax,CMAPMAX,"cmap:crossterm_atom2");
  memoryKK->grow_kokkos(k_crossterm_atom3,crossterm_atom3,nmax,CMAPMAX,"cmap:crossterm_atom3");
  memoryKK->grow_kokkos(k_crossterm_atom4,crossterm_atom4,nmax,CMAPMAX,"cmap:crossterm_atom4");
  memoryKK->grow_kokkos(k_crossterm_atom5,crossterm_atom5,nmax,CMAPMAX,"cmap:crossterm_atom5");

  d_num_crossterm = k_num_crossterm.template view<DeviceType>();
  d_crossterm_type = k_crossterm_type.template view<DeviceType>();
  d_crossterm_atom1 = k_crossterm_atom1.template view<DeviceType>();
  d_crossterm_atom2 = k_crossterm_atom2.template view<DeviceType>();
  d_crossterm_atom3 = k_crossterm_atom3.template view<DeviceType>();
  d_crossterm_atom4 = k_crossterm_atom4.template view<DeviceType>();
  d_crossterm_atom5 = k_crossterm_atom5.template view<DeviceType>();

  // must initialize num_crossterm to 0 for added atoms
  // may never be set for some atoms when data file is read

  for (int i = nmax_previous; i < nmax; i++) k_num_crossterm.h_view(i) = 0;
  nmax_previous = nmax;

  k_num_crossterm.template modify<LMPHostType>();
  k_crossterm_type.template modify<LMPHostType>();
  k_crossterm_atom1.template modify<LMPHostType>();
  k_crossterm_atom2.template modify<LMPHostType>();
  k_crossterm_atom3.template modify<LMPHostType>();
  k_crossterm_atom4.template modify<LMPHostType>();
  k_crossterm_atom5.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

template<class DeviceType>
void FixCMAPKokkos<DeviceType>::copy_arrays(int i, int j, int delflag)
{
  k_num_crossterm.template sync<LMPHostType>();
  k_crossterm_type.template sync<LMPHostType>();
  k_crossterm_atom1.template sync<LMPHostType>();
  k_crossterm_atom2.template sync<LMPHostType>();
  k_crossterm_atom3.template sync<LMPHostType>();
  k_crossterm_atom4.template sync<LMPHostType>();
  k_crossterm_atom5.template sync<LMPHostType>();

  FixCMAP::copy_arrays(i,j,delflag);

  k_num_crossterm.template modify<LMPHostType>();
  k_crossterm_type.template modify<LMPHostType>();
  k_crossterm_atom1.template modify<LMPHostType>();
  k_crossterm_atom2.template modify<LMPHostType>();
  k_crossterm_atom3.template modify<LMPHostType>();
  k_crossterm_atom4.template modify<LMPHostType>();
  k_crossterm_atom5.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

template<class DeviceType>
void FixCMAPKokkos<DeviceType>::set_arrays(int i)
{
  k_num_crossterm.sync_host();
  num_crossterm[i] = 0;
  k_num_crossterm.modify_host();
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixCMAPKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  k_num_crossterm.sync_host();
  k_crossterm_type.sync_host();
  k_crossterm_atom1.sync_host();
  k_crossterm_atom2.sync_host();
  k_crossterm_atom3.sync_host();
  k_crossterm_atom4.sync_host();
  k_crossterm_atom5.sync_host();

  int m = FixCMAP::pack_exchange(i,buf);

  k_num_crossterm.modify_host();
  k_crossterm_type.modify_host();
  k_crossterm_atom1.modify_host();
  k_crossterm_atom2.modify_host();
  k_crossterm_atom3.modify_host();
  k_crossterm_atom4.modify_host();
  k_crossterm_atom5.modify_host();

  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixCMAPKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf)
{
  k_num_crossterm.sync_host();
  k_crossterm_type.sync_host();
  k_crossterm_atom1.sync_host();
  k_crossterm_atom2.sync_host();
  k_crossterm_atom3.sync_host();
  k_crossterm_atom4.sync_host();
  k_crossterm_atom5.sync_host();

  int m = FixCMAP::unpack_exchange(nlocal,buf);

  k_num_crossterm.modify_host();
  k_crossterm_type.modify_host();
  k_crossterm_atom1.modify_host();
  k_crossterm_atom2.modify_host();
  k_crossterm_atom3.modify_host();
  k_crossterm_atom4.modify_host();
  k_crossterm_atom5.modify_host();

  return m;
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixCMAPKokkos<DeviceType>::dihedral_angle_atan2(double fx, double fy, double fz,
                                      double ax, double ay, double az,
                                      double bx, double by, double bz,
                                      double absg) const
{
  // calculate the dihedral angle

  double angle = 0.0, arg1, arg2;

  arg1 = absg*(fx*bx+fy*by+fz*bz);
  arg2 = ax*bx+ay*by+az*bz;

  if (arg1 == 0 && arg2 == 0)
    Kokkos::abort("CMAP: atan2 function cannot take 2 zero arguments");
  else {
    angle = Kokkos::atan2(arg1,arg2);
    angle = angle*180.0/MY_PI;
  }

  return angle;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixCMAPKokkos<DeviceType>::bc_interpol(double x1, double x2, int low1, int low2, double *gs,
                           double *d1gs, double *d2gs, double *d12gs,
                           double &E, double &dEdPhi, double &dEdPsi ) const
{

  // FUSE bc_coeff() and bc_interpol() inline functions for
  // KOKKOS version to avoid passing cij[][] array back and forth

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
  double xx, x[16], cij[4][4];

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

  // for a given point of interest and its corresponding grid square values,
  //   gradients and cross-derivatives
  // calculate the interpolated value of the point of interest (POI)

  double t, u, gs1l, gs2l;

  // set the interpolation coefficients
  // bc_coeff(gs,d1gs,d2gs,d12gs,&cij[0]);

  gs1l = d_g_axis(low1);
  gs2l = d_g_axis(low2);

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



/* ----------------------------------------------------------------------
   return local index of atom J or any of its images that is closest to atom I
   if J is not a valid index like -1, just return it
   copied from domain.cpp
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int FixCMAPKokkos<DeviceType>::closest_image(const int i, int j) const
{
  if (j < 0) return j;

  const X_FLOAT xi0 = d_x(i,0);
  const X_FLOAT xi1 = d_x(i,1);
  const X_FLOAT xi2 = d_x(i,2);

  int closest = j;
  X_FLOAT delx = xi0 - d_x(j,0);
  X_FLOAT dely = xi1 - d_x(j,1);
  X_FLOAT delz = xi2 - d_x(j,2);
  X_FLOAT rsqmin = delx*delx + dely*dely + delz*delz;
  X_FLOAT rsq;

  while (d_sametag[j] >= 0) {
    j = d_sametag[j];
    delx = xi0 - d_x(j,0);
    dely = xi1 - d_x(j,1);
    delz = xi2 - d_x(j,2);
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq < rsqmin) {
      rsqmin = rsq;
      closest = j;
    }
  }

  return closest;
}


namespace LAMMPS_NS {
template class FixCMAPKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixCMAPKokkos<LMPHostType>;
#endif
}
