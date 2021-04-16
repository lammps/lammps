

/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#ifndef LMP_REAXC_DEFS_H
#define LMP_REAXC_DEFS_H

#if defined(__IBMC__)
#define inline __inline__
#endif /*IBMC*/

#ifndef SUCCESS
#define SUCCESS  1
#endif
#ifndef FAILURE
#define FAILURE  0
#endif

#define SQR(x)        ((x)*(x))
#define CUBE(x)       ((x)*(x)*(x))
#define DEG2RAD(a)    ((a)*constPI/180.0)
#define RAD2DEG(a)    ((a)*180.0/constPI)
#define MAX3(x,y,z)   MAX( MAX(x,y), z)

#define constPI        3.14159265
#define C_ele          332.06371
//#define K_B         503.398008   // kcal/mol/K
#define K_B             0.831687   // amu A^2 / ps^2 / K
#define F_CONV          1e6 / 48.88821291 / 48.88821291   // --> amu A / ps^2
#define E_CONV          0.002391   // amu A^2 / ps^2 --> kcal/mol
#define EV_to_KCALpMOL 14.400000   // ElectronVolt --> KCAL per MOLe
#define KCALpMOL_to_EV 23.02       // 23.060549 //KCAL per MOLe --> ElectronVolt
#define ECxA_to_DEBYE   4.803204   // elem. charge * Ang -> debye
#define CAL_to_JOULES   4.184000   // CALories --> JOULES
#define JOULES_to_CAL   1/4.184000 // JOULES --> CALories
#define AMU_to_GRAM     1.6605e-24
#define ANG_to_CM       1e-8
#define AVOGNR          6.0221367e23
#define P_CONV          1e-24 * AVOGNR * JOULES_to_CAL

#define MAX_STR             1024
#define MAX_LINE            1024
#define MAX_TOKENS          1024
#define MAX_TOKEN_LEN       1024

#define NUM_INTRS      10
#define ALMOST_ZERO    1e-10
#define NEG_INF       -1e10
#define NO_BOND        1e-3  // 0.001
#define HB_THRESHOLD   1e-2  // 0.01

#define REAX_MIN_CAP    50
#define REAX_MIN_NBRS   100
#define MIN_HENTRIES    100
#define MAX_BONDS       30
#define MIN_BONDS       25
#define REAX_MIN_HBONDS 25
#define MIN_3BODIES     1000
#define REAX_SAFE_ZONE  1.2
#define REAX_SAFER_ZONE 1.4
#define DANGER_ZONE     0.90
#define LOOSE_ZONE      0.75
#define MAX_3BODY_PARAM 5
#define MAX_4BODY_PARAM 5

#define MASTER_NODE 0

#define MAXREAXBOND 24 /* used in fix_reaxc_bonds.cpp and pair_reaxc.cpp */
#define MAXSPECBOND 24 /* used in fix_reaxc_species.cpp and pair_reaxc.cpp */

/******************* ENUMERATIONS *************************/

enum lists { BONDS, OLD_BONDS, THREE_BODIES,
             HBONDS, FAR_NBRS, DBOS, DDELTAS, LIST_N };

enum message_tags {NONE, INIT_DESCS, ATOM_LINES, BOND_LINES, ANGLE_LINES};

enum interactions {TYP_VOID, TYP_BOND, TYP_THREE_BODY,
  TYP_HBOND, TYP_FAR_NEIGHBOR, TYP_DBO, TYP_DDELTA, TYP_N};

#endif
