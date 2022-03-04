/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Hasan Metin Aktulga, Purdue University
   (now at Lawrence Berkeley National Laboratory, hmaktulga@lbl.gov)

   Heavily modified and adapted for LAMMPS by the LAMMPS developers.
------------------------------------------------------------------------- */

#ifndef LMP_REAXFF_DEFS_H
#define LMP_REAXFF_DEFS_H

#if defined(__IBMC__)
#define inline __inline__
#endif /*IBMC*/

#define SQR(x) ((x) * (x))
#define CUBE(x) ((x) * (x) * (x))
#define DEG2RAD(a) ((a) *constPI / 180.0)
#define RAD2DEG(a) ((a) *180.0 / constPI)
#define MAX3(x, y, z) MAX(MAX(x, y), z)

#define constPI 3.14159265
#define C_ele 332.06371
#define EV_to_KCALpMOL 14.400000    // ElectronVolt --> KCAL per MOLe
#define KCALpMOL_to_EV 23.02        // 23.060549 //KCAL per MOLe --> ElectronVolt

#define HB_THRESHOLD 1e-2    // 0.01

#define REAX_MIN_CAP 50
#define REAX_MIN_NBRS 100
#define MIN_HENTRIES 100
#define MAX_BONDS 30
#define MIN_BONDS 25
#define REAX_MIN_HBONDS 25
#define MIN_3BODIES 1000
#define REAX_SAFE_ZONE 1.2
#define REAX_SAFER_ZONE 1.4
#define DANGER_ZONE 0.90
#define LOOSE_ZONE 0.75

#define MAXREAXBOND 24 /* used in fix_reaxff_bonds.cpp and pair_reaxff.cpp */
#define MAXSPECBOND 24 /* used in fix_reaxff_species.cpp and pair_reaxff.cpp */

#define REAX_MAX_3BODY_PARAM 5
#define REAX_MAX_4BODY_PARAM 5

namespace ReaxFF {
/******************* ENUMERATORS *************************/
enum { BONDS, THREE_BODIES, HBONDS, FAR_NBRS, LIST_N };
enum { TYP_BOND, TYP_THREE_BODY, TYP_HBOND, TYP_FAR_NEIGHBOR, TYP_N };
}    // namespace ReaxFF

#endif
