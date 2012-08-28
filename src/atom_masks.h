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

#ifndef LMP_ATOM_MASK_H
#define LMP_ATOM_MASK_H

// per-atom data masks

#define EMPTY_MASK     0x00000000
#define ALL_MASK       0xffffffff
#define SAMETAG_MASK   0x01000000
#define EXTENDED_MASK  0x80000000

// standard

#define X_MASK         0x00000001
#define V_MASK         0x00000002
#define F_MASK         0x00000004
#define TAG_MASK       0x00000008
#define TYPE_MASK      0x00000010
#define MASK_MASK      0x00000020
#define IMAGE_MASK     0x00000040
#define Q_MASK         0x00000080
#define MOLECULE_MASK  0x00000100
#define RMASS_MASK     0x00000200
#define BOND_MASK      0x00000400
#define ANGLE_MASK     0x00000800
#define DIHEDRAL_MASK  0x00001000
#define IMPROPER_MASK  0x00002000
#define SPECIAL_MASK   0x00004000
#define MAP_MASK       0x00008000

// granular

#define RADIUS_MASK    0x00010000
#define DENSITY_MASK   0x00020000
#define OMEGA_MASK     0x00040000
#define TORQUE_MASK    0x00080000
#define ANGMOM_MASK    0x00100000
#define GRANULAR_MASK  0x001f0000

// peridynamics

#define VFRAC_MASK     0x00200000
#define S0_MASK        0x00400000
#define X0_MASK        0x00800000
#define PERI_MASK      0x00e00000

#define ELLIPSOID_MASK 0x00000001
#define LINE_MASK      0x00000002
#define TRI_MASK       0x00000004

// electron

#define SPIN_MASK      0x00000010
#define ERADIUS_MASK   0x00000020
#define ERVEL_MASK     0x00000040
#define ERFORCE_MASK   0x00000080
#define ERVELFORCE_MASK 0x00000100

#define CS_MASK        0x00000200
#define CSFORCE_MASK   0x00000400
#define VFORCE_MASK    0x00000800

#define ELECTRON_MASK  0x00000f00

// SPH

#define ETAG_MASK      0x00001000
#define RHO_MASK       0x00002000
#define DRHO_MASK      0x00004000
#define E_MASK         0x00008000
#define DE_MASK        0x00010000
#define VEST_MASK      0x00020000
#define CV_MASK        0x00040000

#endif
