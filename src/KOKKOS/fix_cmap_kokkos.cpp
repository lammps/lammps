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
   Contributing author: Mitch Murphy (alphataubio at gmail)
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

using namespace LAMMPS_NS;
using namespace MathConst;

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
  exchange_comm_device = sort_device = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

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

  d_count = typename AT::t_int_scalar("fix_cmap:count");
  h_count = Kokkos::create_mirror_view(d_count);

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

  memoryKK->destroy_kokkos(k_num_crossterm,num_crossterm);
  memoryKK->destroy_kokkos(k_crossterm_type,crossterm_type);
  memoryKK->destroy_kokkos(k_crossterm_atom1,crossterm_atom1);
  memoryKK->destroy_kokkos(k_crossterm_atom2,crossterm_atom2);
  memoryKK->destroy_kokkos(k_crossterm_atom3,crossterm_atom3);
  memoryKK->destroy_kokkos(k_crossterm_atom4,crossterm_atom4);
  memoryKK->destroy_kokkos(k_crossterm_atom5,crossterm_atom5);

  memoryKK->destroy_kokkos(d_crosstermlist);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixCMAPKokkos<DeviceType>::init()
{
  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot yet use respa with Kokkos");

  // on KOKKOS, allocate enough for all crossterms on each GPU to avoid grow operation in device code

  maxcrossterm = ncmap;
  memoryKK->create_kokkos(d_crosstermlist,maxcrossterm,CMAPMAX,"cmap:crosstermlist");
}

/* ----------------------------------------------------------------------
   store local neighbor list as if newton_bond = OFF, even if actually ON
------------------------------------------------------------------------- */

template<class DeviceType>
void FixCMAPKokkos<DeviceType>::pre_neighbor()
{
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

  copymode = 1;
  Kokkos::parallel_scan(Kokkos::RangePolicy<DeviceType,TagFixCmapPreNeighbor>(0,nlocal),*this,ncrosstermlist);
  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixCMAPKokkos<DeviceType>::operator()(TagFixCmapPreNeighbor, const int i, int &l_ncrosstermlist, const bool is_final ) const
{
  for( int m = 0; m < d_num_crossterm(i); m++) {

    int atom1 = AtomKokkos::map_kokkos<DeviceType>(d_crossterm_atom1(i,m),map_style,k_map_array,k_map_hash);
    int atom2 = AtomKokkos::map_kokkos<DeviceType>(d_crossterm_atom2(i,m),map_style,k_map_array,k_map_hash);
    int atom3 = AtomKokkos::map_kokkos<DeviceType>(d_crossterm_atom3(i,m),map_style,k_map_array,k_map_hash);
    int atom4 = AtomKokkos::map_kokkos<DeviceType>(d_crossterm_atom4(i,m),map_style,k_map_array,k_map_hash);
    int atom5 = AtomKokkos::map_kokkos<DeviceType>(d_crossterm_atom5(i,m),map_style,k_map_array,k_map_hash);

    if( atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1 || atom5 == -1)
      Kokkos::abort("CMAP atoms missing on proc");

    atom1 = closest_image(i,atom1);
    atom2 = closest_image(i,atom2);
    atom3 = closest_image(i,atom3);
    atom4 = closest_image(i,atom4);
    atom5 = closest_image(i,atom5);

    if( i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4 && i <= atom5) {
      if (l_ncrosstermlist > maxcrossterm) Kokkos::abort("l_ncrosstermlist > maxcrossterm");
      if(is_final) {
        d_crosstermlist(l_ncrosstermlist,0) = atom1;
        d_crosstermlist(l_ncrosstermlist,1) = atom2;
        d_crosstermlist(l_ncrosstermlist,2) = atom3;
        d_crosstermlist(l_ncrosstermlist,3) = atom4;
        d_crosstermlist(l_ncrosstermlist,4) = atom5;
        d_crosstermlist(l_ncrosstermlist,5) = d_crossterm_type(i,m);
      }
      l_ncrosstermlist++;
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

  int eflag = eflag_caller;
  ev_init(eflag,vflag);

  copymode = 1;
  nlocal = atomKK->nlocal;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixCmapPostForce>(0,ncrosstermlist),*this,ecmap);
  copymode = 0;
  atomKK->modified(execution_space,F_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixCMAPKokkos<DeviceType>::operator()(TagFixCmapPostForce, const int n, double &ecmapKK) const
{
  // Definition of cross-term dihedrals

  //         phi dihedral
  //   |--------------------|
  //   a1-----a2-----a3-----a4-----a5    cross-term atoms
  //   C      N      CA     C      N     cross-term atom types
  //          |--------------------|
  //               psi dihedral

  int i1 = d_crosstermlist(n,0);
  int i2 = d_crosstermlist(n,1);
  int i3 = d_crosstermlist(n,2);
  int i4 = d_crosstermlist(n,3);
  int i5 = d_crosstermlist(n,4);
  int type = d_crosstermlist(n,5);
  if (type == 0) return;

  // calculate bond vectors for both dihedrals

  // phi
  // vb21 = r2 - r1

  double vb21x = d_x(i2,0) - d_x(i1,0);
  double vb21y = d_x(i2,1) - d_x(i1,1);
  double vb21z = d_x(i2,2) - d_x(i1,2);
  double vb12x = -1.0*vb21x;
  double vb12y = -1.0*vb21y;
  double vb12z = -1.0*vb21z;
  double vb32x = d_x(i3,0) - d_x(i2,0);
  double vb32y = d_x(i3,1) - d_x(i2,1);
  double vb32z = d_x(i3,2) - d_x(i2,2);
  double vb23x = -1.0*vb32x;
  double vb23y = -1.0*vb32y;
  double vb23z = -1.0*vb32z;

  double vb34x = d_x(i3,0) - d_x(i4,0);
  double vb34y = d_x(i3,1) - d_x(i4,1);
  double vb34z = d_x(i3,2) - d_x(i4,2);

  // psi
  // bond vectors same as for phi: vb32

  double vb43x = -1.0*vb34x;
  double vb43y = -1.0*vb34y;
  double vb43z = -1.0*vb34z;

  double vb45x = d_x(i4,0) - d_x(i5,0);
  double vb45y = d_x(i4,1) - d_x(i5,1);
  double vb45z = d_x(i4,2) - d_x(i5,2);

  // calculate normal vectors for planes that define the dihedral angles

  double a1x = vb12y*vb23z - vb12z*vb23y;
  double a1y = vb12z*vb23x - vb12x*vb23z;
  double a1z = vb12x*vb23y - vb12y*vb23x;

  double b1x = vb43y*vb23z - vb43z*vb23y;
  double b1y = vb43z*vb23x - vb43x*vb23z;
  double b1z = vb43x*vb23y - vb43y*vb23x;

  double a2x = vb23y*vb34z - vb23z*vb34y;
  double a2y = vb23z*vb34x - vb23x*vb34z;
  double a2z = vb23x*vb34y - vb23y*vb34x;

  double b2x = vb45y*vb43z - vb45z*vb43y;
  double b2y = vb45z*vb43x - vb45x*vb43z;
  double b2z = vb45x*vb43y - vb45y*vb43x;

  // calculate terms used later in calculations

  double r32 = sqrt(vb32x*vb32x + vb32y*vb32y + vb32z*vb32z);
  double a1sq = a1x*a1x + a1y*a1y + a1z*a1z;
  double b1sq = b1x*b1x + b1y*b1y + b1z*b1z;

  double r43 = sqrt(vb43x*vb43x + vb43y*vb43y + vb43z*vb43z);
  double a2sq = a2x*a2x + a2y*a2y + a2z*a2z;
  double b2sq = b2x*b2x + b2y*b2y + b2z*b2z;
  if (a1sq<0.0001 || b1sq<0.0001 || a2sq<0.0001 || b2sq<0.0001) return;

  // vectors needed to calculate the cross-term dihedral angles

  double dpr21r32 = vb21x*vb32x + vb21y*vb32y + vb21z*vb32z;
  double dpr34r32 = vb34x*vb32x + vb34y*vb32y + vb34z*vb32z;
  double dpr32r43 = vb32x*vb43x + vb32y*vb43y + vb32z*vb43z;
  double dpr45r43 = vb45x*vb43x + vb45y*vb43y + vb45z*vb43z;

  // cross-term dihedral angles
  // calculate the backbone dihedral angles as VMD and GROMACS

  double phi = dihedral_angle_atan2(vb21x,vb21y,vb21z,a1x,a1y,a1z,b1x,b1y,b1z,r32);
  double psi = dihedral_angle_atan2(vb32x,vb32y,vb32z,a2x,a2y,a2z,b2x,b2y,b2z,r43);

  if (phi == 180.0) phi= -180.0;
  if (psi == 180.0) psi= -180.0;

  double phi1 = phi;
  if (phi1 < 0.0) phi1 += 360.0;
  double psi1 = psi;
  if (psi1 < 0.0) psi1 += 360.0;

  // find the neighbor grid point index

  int li1 = int(((phi1+CMAPXMIN2)/CMAPDX)+((CMAPDIM*1.0)/2.0));
  int li2 = int(((psi1+CMAPXMIN2)/CMAPDX)+((CMAPDIM*1.0)/2.0));
  int li3 = int((phi-CMAPXMIN2)/CMAPDX);
  int li4 = int((psi-CMAPXMIN2)/CMAPDX);
  int mli3 = li3 % CMAPDIM;
  int mli4 = li4 % CMAPDIM;
  int mli31 = (li3+1) % CMAPDIM;
  int mli41 = (li4+1)  %CMAPDIM;
  int mli1 = li1 % CMAPDIM;
  int mli2 = li2 % CMAPDIM;
  int mli11 = (li1+1) % CMAPDIM;
  int mli21 = (li2+1)  %CMAPDIM;
  int t1 = type-1;
  if (t1 < 0 || t1 > 5) Kokkos::abort("Invalid CMAP crossterm_type");

  // determine the values and derivatives for the grid square points

  double gs[4],d1gs[4],d2gs[4],d12gs[4];

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
  // needed for compute_scalar()

  double engfraction = 0.2 * E;
  if (i1 < nlocal) ecmapKK += engfraction;
  if (i2 < nlocal) ecmapKK += engfraction;
  if (i3 < nlocal) ecmapKK += engfraction;
  if (i4 < nlocal) ecmapKK += engfraction;
  if (i5 < nlocal) ecmapKK += engfraction;

  // calculate the derivatives dphi/dr_i

  double dphidr1x = 1.0*r32/a1sq*a1x;
  double dphidr1y = 1.0*r32/a1sq*a1y;
  double dphidr1z = 1.0*r32/a1sq*a1z;

  double dphidr2x = -1.0*r32/a1sq*a1x - dpr21r32/a1sq/r32*a1x + dpr34r32/b1sq/r32*b1x;
  double dphidr2y = -1.0*r32/a1sq*a1y - dpr21r32/a1sq/r32*a1y + dpr34r32/b1sq/r32*b1y;
  double dphidr2z = -1.0*r32/a1sq*a1z - dpr21r32/a1sq/r32*a1z + dpr34r32/b1sq/r32*b1z;

  double dphidr3x = dpr34r32/b1sq/r32*b1x - dpr21r32/a1sq/r32*a1x - r32/b1sq*b1x;
  double dphidr3y = dpr34r32/b1sq/r32*b1y - dpr21r32/a1sq/r32*a1y - r32/b1sq*b1y;
  double dphidr3z = dpr34r32/b1sq/r32*b1z - dpr21r32/a1sq/r32*a1z - r32/b1sq*b1z;

  double dphidr4x = r32/b1sq*b1x;
  double dphidr4y = r32/b1sq*b1y;
  double dphidr4z = r32/b1sq*b1z;

  // calculate the derivatives dpsi/dr_i

  double dpsidr1x = 1.0*r43/a2sq*a2x;
  double dpsidr1y = 1.0*r43/a2sq*a2y;
  double dpsidr1z = 1.0*r43/a2sq*a2z;

  double dpsidr2x = r43/a2sq*a2x + dpr32r43/a2sq/r43*a2x - dpr45r43/b2sq/r43*b2x;
  double dpsidr2y = r43/a2sq*a2y + dpr32r43/a2sq/r43*a2y - dpr45r43/b2sq/r43*b2y;
  double dpsidr2z = r43/a2sq*a2z + dpr32r43/a2sq/r43*a2z - dpr45r43/b2sq/r43*b2z;

  double dpsidr3x = dpr45r43/b2sq/r43*b2x - dpr32r43/a2sq/r43*a2x - r43/b2sq*b2x;
  double dpsidr3y = dpr45r43/b2sq/r43*b2y - dpr32r43/a2sq/r43*a2y - r43/b2sq*b2y;
  double dpsidr3z = dpr45r43/b2sq/r43*b2z - dpr32r43/a2sq/r43*a2z - r43/b2sq*b2z;

  double dpsidr4x = r43/b2sq*b2x;
  double dpsidr4y = r43/b2sq*b2y;
  double dpsidr4z = r43/b2sq*b2z;

  // calculate forces on cross-term atoms: F = -(dE/dPhi)*(dPhi/dr)
  // apply force to each of the 5 atoms

  if (i1 < nlocal) {
    Kokkos::atomic_add(&d_f(i1,0), dEdPhi*dphidr1x);
    Kokkos::atomic_add(&d_f(i1,1), dEdPhi*dphidr1y);
    Kokkos::atomic_add(&d_f(i1,2), dEdPhi*dphidr1z);
  }
  if (i2 < nlocal) {
    Kokkos::atomic_add(&d_f(i2,0), dEdPhi*dphidr2x + dEdPsi*dpsidr1x);
    Kokkos::atomic_add(&d_f(i2,1), dEdPhi*dphidr2y + dEdPsi*dpsidr1y);
    Kokkos::atomic_add(&d_f(i2,2), dEdPhi*dphidr2z + dEdPsi*dpsidr1z);
  }
  if (i3 < nlocal) {
    Kokkos::atomic_add(&d_f(i3,0), -dEdPhi*dphidr3x - dEdPsi*dpsidr2x);
    Kokkos::atomic_add(&d_f(i3,1), -dEdPhi*dphidr3y - dEdPsi*dpsidr2y);
    Kokkos::atomic_add(&d_f(i3,2), -dEdPhi*dphidr3z - dEdPsi*dpsidr2z);
  }
  if (i4 < nlocal) {
    Kokkos::atomic_add(&d_f(i4,0), -dEdPhi*dphidr4x - dEdPsi*dpsidr3x);
    Kokkos::atomic_add(&d_f(i4,1), -dEdPhi*dphidr4y - dEdPsi*dpsidr3y);
    Kokkos::atomic_add(&d_f(i4,2), -dEdPhi*dphidr4z - dEdPsi*dpsidr3z);
  }
  if (i5 < nlocal) {
    Kokkos::atomic_add(&d_f(i5,0), -dEdPsi*dpsidr4x);
    Kokkos::atomic_add(&d_f(i5,1), -dEdPsi*dpsidr4y);
    Kokkos::atomic_add(&d_f(i5,2), -dEdPsi*dpsidr4z);
  }
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
   sort local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixCMAPKokkos<DeviceType>::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  // always sort on the device

  k_num_crossterm.sync_device();
  k_crossterm_type.sync_device();
  k_crossterm_atom1.sync_device();
  k_crossterm_atom2.sync_device();
  k_crossterm_atom3.sync_device();
  k_crossterm_atom4.sync_device();
  k_crossterm_atom5.sync_device();

  Sorter.sort(LMPDeviceType(), k_num_crossterm.d_view);
  Sorter.sort(LMPDeviceType(), k_crossterm_type.d_view);
  Sorter.sort(LMPDeviceType(), k_crossterm_atom1.d_view);
  Sorter.sort(LMPDeviceType(), k_crossterm_atom2.d_view);
  Sorter.sort(LMPDeviceType(), k_crossterm_atom3.d_view);
  Sorter.sort(LMPDeviceType(), k_crossterm_atom4.d_view);
  Sorter.sort(LMPDeviceType(), k_crossterm_atom5.d_view);

  k_num_crossterm.modify_device();
  k_crossterm_type.modify_device();
  k_crossterm_atom1.modify_device();
  k_crossterm_atom2.modify_device();
  k_crossterm_atom3.modify_device();
  k_crossterm_atom4.modify_device();
  k_crossterm_atom5.modify_device();
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

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange
------------------------------------------------------------------------- */

template<class DeviceType>
int FixCMAPKokkos<DeviceType>::pack_exchange_kokkos(
   const int &nsend, DAT::tdual_xfloat_2d &k_buf,
   DAT::tdual_int_1d k_exchange_sendlist, DAT::tdual_int_1d k_copylist,
   ExecutionSpace space)
{
  k_buf.template sync<DeviceType>();
  k_copylist.template sync<DeviceType>();
  k_exchange_sendlist.template sync<DeviceType>();

  k_num_crossterm.template sync<DeviceType>();
  k_crossterm_type.template sync<DeviceType>();
  k_crossterm_atom1.template sync<DeviceType>();
  k_crossterm_atom2.template sync<DeviceType>();
  k_crossterm_atom3.template sync<DeviceType>();
  k_crossterm_atom4.template sync<DeviceType>();
  k_crossterm_atom5.template sync<DeviceType>();

  auto d_buf = typename ArrayTypes<DeviceType>::t_xfloat_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));
  auto d_copylist = k_copylist.template view<DeviceType>();
  auto d_exchange_sendlist = k_exchange_sendlist.template view<DeviceType>();

  Kokkos::deep_copy(d_count,0);

  auto l_num_crossterm = d_num_crossterm;
  auto l_crossterm_type = d_crossterm_type;
  auto l_crossterm_atom1 = d_crossterm_atom1;
  auto l_crossterm_atom2 = d_crossterm_atom2;
  auto l_crossterm_atom3 = d_crossterm_atom3;
  auto l_crossterm_atom4 = d_crossterm_atom4;
  auto l_crossterm_atom5 = d_crossterm_atom5;
  auto l_count = d_count;

  copymode = 1;

  Kokkos::parallel_scan(Kokkos::RangePolicy<DeviceType>(0,nsend), KOKKOS_LAMBDA(const int &mysend, int &offset, const bool &final) {

    const int i = d_exchange_sendlist(mysend);

    if (!final) offset += (1+l_num_crossterm(i)*6);
    else {

      int m = nsend + offset;
      d_buf(mysend) = d_ubuf(m).d;
      d_buf(m++) = d_ubuf(l_num_crossterm(i)).d;

      for (int k = 0; k < l_num_crossterm(i); k++) {
        d_buf(m++) = d_ubuf(l_crossterm_type(i,k)).d;
        d_buf(m++) = d_ubuf(l_crossterm_atom1(i,k)).d;
        d_buf(m++) = d_ubuf(l_crossterm_atom2(i,k)).d;
        d_buf(m++) = d_ubuf(l_crossterm_atom3(i,k)).d;
        d_buf(m++) = d_ubuf(l_crossterm_atom4(i,k)).d;
        d_buf(m++) = d_ubuf(l_crossterm_atom5(i,k)).d;
      }

      if (mysend == nsend-1) l_count() = m;
      offset = m - nsend;

      const int j = d_copylist(mysend);
      if (j > -1) {
        l_num_crossterm(i) = l_num_crossterm(j);
        for (int k = 0; k < l_num_crossterm(i); k++) {
          l_crossterm_type(i,k) = l_crossterm_type(j,k);
          l_crossterm_atom1(i,k) = l_crossterm_atom1(j,k);
          l_crossterm_atom2(i,k) = l_crossterm_atom2(j,k);
          l_crossterm_atom3(i,k) = l_crossterm_atom3(j,k);
          l_crossterm_atom4(i,k) = l_crossterm_atom4(j,k);
          l_crossterm_atom5(i,k) = l_crossterm_atom5(j,k);
        }
      }
    }
  });

  copymode = 0;

  k_buf.template modify<DeviceType>();
  if (space == Host) k_buf.template sync<LMPHostType>();
  else k_buf.template sync<LMPDeviceType>();

  k_num_crossterm.template modify<DeviceType>();
  k_crossterm_type.template modify<DeviceType>();
  k_crossterm_atom1.template modify<DeviceType>();
  k_crossterm_atom2.template modify<DeviceType>();
  k_crossterm_atom3.template modify<DeviceType>();
  k_crossterm_atom4.template modify<DeviceType>();
  k_crossterm_atom5.template modify<DeviceType>();

  Kokkos::deep_copy(h_count,d_count);
  return h_count();
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange
------------------------------------------------------------------------- */

template <class DeviceType>
void FixCMAPKokkos<DeviceType>::unpack_exchange_kokkos(
  DAT::tdual_xfloat_2d &k_buf, DAT::tdual_int_1d &k_indices, int nrecv,
  int nrecv1, int nextrarecv1, ExecutionSpace /*space*/)
{
  k_buf.template sync<DeviceType>();
  k_indices.template sync<DeviceType>();

  k_num_crossterm.template sync<DeviceType>();
  k_crossterm_type.template sync<DeviceType>();
  k_crossterm_atom1.template sync<DeviceType>();
  k_crossterm_atom2.template sync<DeviceType>();
  k_crossterm_atom3.template sync<DeviceType>();
  k_crossterm_atom4.template sync<DeviceType>();
  k_crossterm_atom5.template sync<DeviceType>();

  auto d_buf = typename ArrayTypes<DeviceType>::t_xfloat_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));

  auto d_indices = k_indices.template view<DeviceType>();

  auto l_num_crossterm = d_num_crossterm;
  auto l_crossterm_type = d_crossterm_type;
  auto l_crossterm_atom1 = d_crossterm_atom1;
  auto l_crossterm_atom2 = d_crossterm_atom2;
  auto l_crossterm_atom3 = d_crossterm_atom3;
  auto l_crossterm_atom4 = d_crossterm_atom4;
  auto l_crossterm_atom5 = d_crossterm_atom5;

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nrecv), KOKKOS_LAMBDA(const int &i) {
    int index = d_indices(i);
    if (index > -1) {
      int m = d_ubuf(d_buf(i)).i;
      if (i >= nrecv1) m = nextrarecv1 + d_ubuf(d_buf(nextrarecv1 + i - nrecv1)).i;
      l_num_crossterm(index) = static_cast<int> (d_ubuf(d_buf(m++)).i);
      for (int k = 0; k < l_num_crossterm(index); k++) {
        l_crossterm_type(index,k) = static_cast<int> (d_ubuf(d_buf(m++)).i);
        l_crossterm_atom1(index,k) = static_cast<tagint> (d_ubuf(d_buf(m++)).i);
        l_crossterm_atom2(index,k) = static_cast<tagint> (d_ubuf(d_buf(m++)).i);
        l_crossterm_atom3(index,k) = static_cast<tagint> (d_ubuf(d_buf(m++)).i);
        l_crossterm_atom4(index,k) = static_cast<tagint> (d_ubuf(d_buf(m++)).i);
        l_crossterm_atom5(index,k) = static_cast<tagint> (d_ubuf(d_buf(m++)).i);
      }
    }
  });

  copymode = 0;

  k_num_crossterm.template modify<DeviceType>();
  k_crossterm_type.template modify<DeviceType>();
  k_crossterm_atom1.template modify<DeviceType>();
  k_crossterm_atom2.template modify<DeviceType>();
  k_crossterm_atom3.template modify<DeviceType>();
  k_crossterm_atom4.template modify<DeviceType>();
  k_crossterm_atom5.template modify<DeviceType>();
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

  const int wt[16][16] =
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
