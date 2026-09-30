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
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include "pair_oxdna3_hbond.h"

#include "atom.h"
#include "comm.h"
#include "constants_oxdna.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "mf_oxdna.h"
#include "neigh_list.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>
#include <cassert>

using namespace LAMMPS_NS;
using namespace MFOxdna;

/* ---------------------------------------------------------------------- */

PairOxdna3Hbond::PairOxdna3Hbond(LAMMPS *lmp) : PairOxdnaHbond(lmp)
{
  single_enable = 0;
  writedata = 0;
  trim_flag = 0;

  // sequence-specific base-pairing strength
  // A:0 C:1 G:2 T:3, 5'- [i][j] -3'

  alpha_hb[0][0] = 1.00000;
  alpha_hb[0][1] = 1.00000;
  alpha_hb[0][2] = 1.00000;
  alpha_hb[0][3] = 0.6493620379646540;

  alpha_hb[1][0] = 1.00000;
  alpha_hb[1][1] = 1.00000;
  alpha_hb[1][2] = 1.1999420813642658;
  alpha_hb[1][3] = 1.00000;

  alpha_hb[2][0] = 1.00000;
  alpha_hb[2][1] = 1.1999420813642658;
  alpha_hb[2][2] = 1.00000;
  alpha_hb[2][3] = 1.00000;

  alpha_hb[3][0] = 0.6493620379646540;
  alpha_hb[3][1] = 1.00000;
  alpha_hb[3][2] = 1.00000;
  alpha_hb[3][3] = 1.00000;

}

/* ----------------------------------------------------------------------
   set coeffs
------------------------------------------------------------------------- */
void PairOxdna3Hbond::coeff(int narg, char **arg)
{
  int count;

  if (narg != 3) error->all(FLERR,"Incorrect args for pair coefficients in oxdna3/hbond, use potential file" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi,imod4,jmod4;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  if (ihi>4 || jhi>4) error->all(FLERR, "pair oxdna3/hbond does not support more than 4 atom types for A, C, G and T");

  // h-bonding interaction
  count = 0;

  double epsilon_hb_one, a_hb_one, cut_hb_0_one, cut_hb_c_one, cut_hb_lo_one, cut_hb_hi_one;
  double b_hb_lo_one, b_hb_hi_one, cut_hb_lc_one, cut_hb_hc_one, tmp, shift_hb_one;

  double a_hb1_one, theta_hb1_0_one, dtheta_hb1_ast_one;
  double b_hb1_one, dtheta_hb1_c_one;

  double a_hb2_one, theta_hb2_0_one, dtheta_hb2_ast_one;
  double b_hb2_one, dtheta_hb2_c_one;

  double a_hb3_one, theta_hb3_0_one, dtheta_hb3_ast_one;
  double b_hb3_one, dtheta_hb3_c_one;

  double a_hb4_at, theta_hb4_0_at, dtheta_hb4_ast_at;
  double a_hb4_cg, theta_hb4_0_cg, dtheta_hb4_ast_cg;
  double b_hb4_at, dtheta_hb4_c_at;
  double b_hb4_cg, dtheta_hb4_c_cg;

  double a_hb7_one, theta_hb7_0_one, dtheta_hb7_ast_one;
  double b_hb7_one, dtheta_hb7_c_one;

  double a_hb8_one, theta_hb8_0_one, dtheta_hb8_ast_one;
  double b_hb8_one, dtheta_hb8_c_one;

  // read values from potential file
  if (comm->me == 0) {
    PotentialFileReader reader(lmp, arg[2], "oxdna3 potential", " (hbond)");
    reader.set_bufsize(65336);
    char * line;
    std::string iloc, jloc, potential_name;

    while ((line = reader.next_line())) {
     try {
        ValueTokenizer values(line);
        iloc = values.next_string();
        jloc = values.next_string();
        potential_name = values.next_string();
        if (iloc == arg[0] && jloc == arg[1] && potential_name == "hbond") {

          epsilon_hb_one = values.next_double();
          a_hb_one = values.next_double();
          cut_hb_0_one = values.next_double();
          cut_hb_c_one = values.next_double();
          cut_hb_lo_one = values.next_double();
          cut_hb_hi_one = values.next_double();

          a_hb1_one = values.next_double();
          theta_hb1_0_one = values.next_double();
          dtheta_hb1_ast_one = values.next_double();

          a_hb2_one = values.next_double();
          theta_hb2_0_one = values.next_double();
          dtheta_hb2_ast_one = values.next_double();

          a_hb3_one = values.next_double();
          theta_hb3_0_one = values.next_double();
          dtheta_hb3_ast_one = values.next_double();

          a_hb4_at = values.next_double();
          theta_hb4_0_at = values.next_double();
          dtheta_hb4_ast_at = values.next_double();

          a_hb4_cg = values.next_double();
          theta_hb4_0_cg = values.next_double();
          dtheta_hb4_ast_cg = values.next_double();

          a_hb7_one = values.next_double();
          theta_hb7_0_one = values.next_double();
          dtheta_hb7_ast_one = values.next_double();

          a_hb8_one = values.next_double();
          theta_hb8_0_one = values.next_double();
          dtheta_hb8_ast_one = values.next_double();

          break;
        } else continue;
      } catch (std::exception &e) {
        error->one(FLERR, "Problem parsing oxDNA3 potential file: {}", e.what());
      }
    }
    if ((iloc != arg[0]) || (jloc != arg[1]) || (potential_name != "hbond"))
      error->one(FLERR, "No corresponding hbond potential found in file {} for pair type {} {}",
                 arg[2], arg[0], arg[1]);
  }

  MPI_Bcast(&epsilon_hb_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&a_hb_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_hb_0_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_hb_c_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_hb_lo_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_hb_hi_one, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&a_hb1_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_hb1_0_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_hb1_ast_one, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&a_hb2_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_hb2_0_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_hb2_ast_one, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&a_hb3_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_hb3_0_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_hb3_ast_one, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&a_hb4_at, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_hb4_0_at, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_hb4_ast_at, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&a_hb4_cg, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_hb4_0_cg, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_hb4_ast_cg, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&a_hb7_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_hb7_0_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_hb7_ast_one, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&a_hb8_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_hb8_0_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_hb8_ast_one, 1, MPI_DOUBLE, 0, world);

  b_hb_lo_one = 2*a_hb_one*exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one))*
        2*a_hb_one*exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one))*
        (1-exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one)))*
        (1-exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one)))/
        (4*((1-exp(-a_hb_one*(cut_hb_lo_one -cut_hb_0_one)))*
        (1-exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one)))-
        (1-exp(-a_hb_one*(cut_hb_c_one -cut_hb_0_one)))*
        (1-exp(-a_hb_one*(cut_hb_c_one-cut_hb_0_one)))));

  cut_hb_lc_one = cut_hb_lo_one - a_hb_one*exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one))*
        (1-exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one)))/b_hb_lo_one;

  b_hb_hi_one = 2*a_hb_one*exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one))*
        2*a_hb_one*exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one))*
        (1-exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one)))*
        (1-exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one)))/
        (4*((1-exp(-a_hb_one*(cut_hb_hi_one -cut_hb_0_one)))*
        (1-exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one)))-
        (1-exp(-a_hb_one*(cut_hb_c_one -cut_hb_0_one)))*
        (1-exp(-a_hb_one*(cut_hb_c_one-cut_hb_0_one)))));

  cut_hb_hc_one = cut_hb_hi_one - a_hb_one*exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one))*
        (1-exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one)))/b_hb_hi_one;

  tmp = 1 - exp(-(cut_hb_c_one-cut_hb_0_one) * a_hb_one);
  shift_hb_one = epsilon_hb_one * tmp * tmp;

  b_hb1_one = a_hb1_one*a_hb1_one*dtheta_hb1_ast_one*dtheta_hb1_ast_one/(1-a_hb1_one*dtheta_hb1_ast_one*dtheta_hb1_ast_one);
  dtheta_hb1_c_one = 1/(a_hb1_one*dtheta_hb1_ast_one);

  b_hb2_one = a_hb2_one*a_hb2_one*dtheta_hb2_ast_one*dtheta_hb2_ast_one/(1-a_hb2_one*dtheta_hb2_ast_one*dtheta_hb2_ast_one);
  dtheta_hb2_c_one = 1/(a_hb2_one*dtheta_hb2_ast_one);

  b_hb3_one = a_hb3_one*a_hb3_one*dtheta_hb3_ast_one*dtheta_hb3_ast_one/(1-a_hb3_one*dtheta_hb3_ast_one*dtheta_hb3_ast_one);
  dtheta_hb3_c_one = 1/(a_hb3_one*dtheta_hb3_ast_one);

  b_hb4_at = a_hb4_at*a_hb4_at*dtheta_hb4_ast_at*dtheta_hb4_ast_at/(1-a_hb4_at*dtheta_hb4_ast_at*dtheta_hb4_ast_at);
  dtheta_hb4_c_at = 1/(a_hb4_at*dtheta_hb4_ast_at);

  b_hb4_cg = a_hb4_cg*a_hb4_cg*dtheta_hb4_ast_cg*dtheta_hb4_ast_cg/(1-a_hb4_cg*dtheta_hb4_ast_cg*dtheta_hb4_ast_cg);
  dtheta_hb4_c_cg = 1/(a_hb4_cg*dtheta_hb4_ast_cg);

  b_hb7_one = a_hb7_one*a_hb7_one*dtheta_hb7_ast_one*dtheta_hb7_ast_one/(1-a_hb7_one*dtheta_hb7_ast_one*dtheta_hb7_ast_one);
  dtheta_hb7_c_one = 1/(a_hb7_one*dtheta_hb7_ast_one);

  b_hb8_one = a_hb8_one*a_hb8_one*dtheta_hb8_ast_one*dtheta_hb8_ast_one/(1-a_hb8_one*dtheta_hb8_ast_one*dtheta_hb8_ast_one);
  dtheta_hb8_c_one = 1/(a_hb8_one*dtheta_hb8_ast_one);

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {

      imod4 = i%4;
      if (imod4 == 0) imod4 = 4;
      jmod4 = j%4;
      if (jmod4 == 0) jmod4 = 4;

      epsilon_hb[i][j] = epsilon_hb_one;
      epsilon_hb[i][j] *= alpha_hb[imod4-1][jmod4-1];
      a_hb[i][j] = a_hb_one;
      cut_hb_0[i][j] = cut_hb_0_one;
      cut_hb_c[i][j] = cut_hb_c_one;
      cut_hb_lo[i][j] = cut_hb_lo_one;
      cut_hb_hi[i][j] = cut_hb_hi_one;
      cut_hb_lc[i][j] = cut_hb_lc_one;
      cut_hb_hc[i][j] = cut_hb_hc_one;
      b_hb_lo[i][j] = b_hb_lo_one;
      b_hb_hi[i][j] = b_hb_hi_one;
      shift_hb[i][j] = shift_hb_one;
      shift_hb[i][j] *= alpha_hb[imod4-1][jmod4-1];

      a_hb1[i][j] = a_hb1_one;
      theta_hb1_0[i][j] = theta_hb1_0_one;
      dtheta_hb1_ast[i][j] = dtheta_hb1_ast_one;
      b_hb1[i][j] = b_hb1_one;
      dtheta_hb1_c[i][j] = dtheta_hb1_c_one;

      a_hb2[i][j] = a_hb2_one;
      theta_hb2_0[i][j] = theta_hb2_0_one;
      dtheta_hb2_ast[i][j] = dtheta_hb2_ast_one;
      b_hb2[i][j] = b_hb2_one;
      dtheta_hb2_c[i][j] = dtheta_hb2_c_one;

      a_hb3[i][j] = a_hb3_one;
      theta_hb3_0[i][j] = theta_hb3_0_one;
      dtheta_hb3_ast[i][j] = dtheta_hb3_ast_one;
      b_hb3[i][j] = b_hb3_one;
      dtheta_hb3_c[i][j] = dtheta_hb3_c_one;

      if((imod4==1 && jmod4==4) || (imod4==4 && jmod4==1)){
        a_hb4[i][j] = a_hb4_at;
        theta_hb4_0[i][j] = theta_hb4_0_at;
        dtheta_hb4_ast[i][j] = dtheta_hb4_ast_at;
        b_hb4[i][j] = b_hb4_at;
        dtheta_hb4_c[i][j] = dtheta_hb4_c_at;
      }

      if((imod4==2 && jmod4==3) || (imod4==3 && jmod4==2)){
        a_hb4[i][j] = a_hb4_cg;
        theta_hb4_0[i][j] = theta_hb4_0_cg;
        dtheta_hb4_ast[i][j] = dtheta_hb4_ast_cg;
        b_hb4[i][j] = b_hb4_cg;
        dtheta_hb4_c[i][j] = dtheta_hb4_c_cg;
      }

      a_hb7[i][j] = a_hb7_one;
      theta_hb7_0[i][j] = theta_hb7_0_one;
      dtheta_hb7_ast[i][j] = dtheta_hb7_ast_one;
      b_hb7[i][j] = b_hb7_one;
      dtheta_hb7_c[i][j] = dtheta_hb7_c_one;

      a_hb8[i][j] = a_hb8_one;
      theta_hb8_0[i][j] = theta_hb8_0_one;
      dtheta_hb8_ast[i][j] = dtheta_hb8_ast_one;
      b_hb8[i][j] = b_hb8_one;
      dtheta_hb8_c[i][j] = dtheta_hb8_c_one;

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna3/hbond" + utils::errorurl(21));

}
