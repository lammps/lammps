// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_ATOM_MASK_H
#define LMP_ATOM_MASK_H

// per-atom data masks

#define EMPTY_MASK     0x0000000000000000
#define ALL_MASK       0xffffffffffffffff

// standard

#define X_MASK         0x0000000000000001
#define V_MASK         0x0000000000000002
#define F_MASK         0x0000000000000004
#define TAG_MASK       0x0000000000000008
#define TYPE_MASK      0x0000000000000010
#define MASK_MASK      0x0000000000000020
#define IMAGE_MASK     0x0000000000000040
#define Q_MASK         0x0000000000000080
#define MOLECULE_MASK  0x0000000000000100
#define RMASS_MASK     0x0000000000000200
#define BOND_MASK      0x0000000000000400
#define ANGLE_MASK     0x0000000000000800
#define DIHEDRAL_MASK  0x0000000000001000
#define IMPROPER_MASK  0x0000000000002000
#define SPECIAL_MASK   0x0000000000004000
#define ENERGY_MASK    0x0000000000008000
#define VIRIAL_MASK    0x0000000000010000
#define MU_MASK        0x0000000000020000

// SPIN

#define SP_MASK        0x0000000000040000
#define FM_MASK        0x0000000000080000
#define FML_MASK       0x0000000000100000

// DPD

#define DPDRHO_MASK    0x0000000000200000
#define DPDTHETA_MASK  0x0000000000400000
#define UCOND_MASK     0x0000000000800000
#define UMECH_MASK     0x0000000001000000
#define UCHEM_MASK     0x0000000002000000
#define UCG_MASK       0x0000000004000000
#define UCGNEW_MASK    0x0000000008000000
#define DUCHEM_MASK    0x0000000010000000
#define DVECTOR_MASK   0x0000000020000000

// granular

#define RADIUS_MASK    0x0000000040000000
#define OMEGA_MASK     0x0000000080000000
#define TORQUE_MASK    0x0000000100000000
#define ANGMOM_MASK    0x0000000200000000

#endif
