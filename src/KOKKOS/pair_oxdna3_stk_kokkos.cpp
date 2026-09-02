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

#include "pair_oxdna3_stk_kokkos.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "mf_oxdna.h"
#include "potential_file_reader.h"
#include "math_special.h"

#include <cassert>
#include <cmath>

using namespace LAMMPS_NS;
using namespace MathSpecial;
using namespace MFOxdna;

/* ----------------------------------------------------------------------
   IMPORTANT NOTE ! We entirely code duplicate PairOxdna3StkKokkos::coeff
   into PairOxdna3Stk::coeff. So any edits made in one need to manually be
   made to the other !
   The vanilla version is in: src/CG-DNA/pair_oxdna3_stk.cpp
------------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdna3StkKokkos<DeviceType>::coeff(int narg, char **arg)
{
  // Due to class templating of DeviceType, we need this-> on everything. We use local variables
  // so that we can as much as possible just copy-paste the vanilla code (it's cleaner this way also).
  auto *error = this->error;
  auto *atom = this->atom;
  auto *comm = this->comm;
  auto *lmp = this->lmp;
  MPI_Comm world = this->world;

  int count;

  if (narg != 4) error->all(FLERR,"Incorrect args for pair coefficients in oxdna3/stk, use potential file" + utils::errorurl(21));
  if (!this->allocated) this->allocate();

  // NOTE: allocate() needs this-> still, but otherwise this is a direct copy and paste from the
  // vanilla code. These pointer aliases must be taken AFTER allocate(), since allocate() is what
  // assigns the underlying member pointers.

  auto *setflag = this->setflag;

  auto *epsilon_st = this->epsilon_st;
  auto *a_st = this->a_st;
  auto *cut_st_0 = this->cut_st_0;
  auto *cut_st_c = this->cut_st_c;
  auto *cut_st_lo = this->cut_st_lo;
  auto *cut_st_hi = this->cut_st_hi;
  auto *cut_st_lc = this->cut_st_lc;
  auto *cut_st_hc = this->cut_st_hc;
  auto *b_st_lo = this->b_st_lo;
  auto *b_st_hi = this->b_st_hi;
  auto *shift_st = this->shift_st;
  auto *cutsq_st_hc = this->cutsq_st_hc;

  auto *a_st4 = this->a_st4;
  auto *theta_st4_0 = this->theta_st4_0;
  auto *dtheta_st4_ast = this->dtheta_st4_ast;
  auto *b_st4 = this->b_st4;
  auto *dtheta_st4_c = this->dtheta_st4_c;

  auto *a_st5 = this->a_st5;
  auto *theta_st5_0 = this->theta_st5_0;
  auto *dtheta_st5_ast = this->dtheta_st5_ast;
  auto *b_st5 = this->b_st5;
  auto *dtheta_st5_c = this->dtheta_st5_c;

  auto *a_st6 = this->a_st6;
  auto *theta_st6_0 = this->theta_st6_0;
  auto *dtheta_st6_ast = this->dtheta_st6_ast;
  auto *b_st6 = this->b_st6;
  auto *dtheta_st6_c = this->dtheta_st6_c;

  auto *a_st1 = this->a_st1;
  auto *cosphi_st1_ast = this->cosphi_st1_ast;
  auto *b_st1 = this->b_st1;
  auto *cosphi_st1_c = this->cosphi_st1_c;
  auto *a_st2 = this->a_st2;
  auto *cosphi_st2_ast = this->cosphi_st2_ast;
  auto *b_st2 = this->b_st2;
  auto *cosphi_st2_c = this->cosphi_st2_c;

  // START OF VANILLA CODE DUPLICATION

  int ilo,ihi,jlo,jhi,nlo,nhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  assert((ilo == jlo) & (ihi == jhi));
  nlo = ilo;
  nhi = ihi;

  if (nhi > 4) error->all(FLERR, "pair oxdna3/stk does not support more than 4 atom types for A, C, G and T");

  // stacking interaction
  count = 0;

  double tmp, theta_st4_0_one;
  double T, epsilon_st_one, xi_st_one, kappa_st_one, a_st_one, b_st_lo_one, b_st_hi_one;

  double a_st5_one, theta_st5_0_one, dtheta_st5_ast_one;
  double b_st5_one, dtheta_st5_c_one;

  double a_st6_one, theta_st6_0_one, dtheta_st6_ast_one;
  double b_st6_one, dtheta_st6_c_one;

  double a_st1_one, cosphi_st1_ast_one, b_st1_one, cosphi_st1_c_one;
  double a_st2_one, cosphi_st2_ast_one, b_st2_one, cosphi_st2_c_one;

  T = utils::numeric(FLERR,arg[2],false,lmp);

  for (int i = 0; i <= nhi; i++) { // type 0 for terminal j
    for (int j = 0; j <= nhi; j++) {
      for (int k = 0; k <= nhi; k++) {
        for (int l = 0; l <= nhi; l++) { // type 0 for terminal k
          cut_st_0[i][j][k][l] = 0.0;
          cut_st_c[i][j][k][l] = 0.0;
          cut_st_lo[i][j][k][l] = 0.0;
          cut_st_hi[i][j][k][l] = 0.0;
          a_st4[i][j][k][l] = 0.0;
          dtheta_st4_ast[i][j][k][l] = 0.0;
        }
      }
    }
  }

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, arg[3], "oxdna3 potential", " (stk)");
    reader.set_bufsize(65336);
    char * line;
    std::string iloc, jloc, potential_name;

    while ((line = reader.next_line())) {
      try {
        ValueTokenizer values(line);
        iloc = values.next_string();
        jloc = values.next_string();
        potential_name = values.next_string();
        if (iloc == arg[0] && jloc == arg[1] && potential_name == "stk") {

          xi_st_one = values.next_double();
          kappa_st_one = values.next_double();
          epsilon_st_one = this->stacking_strength(xi_st_one, kappa_st_one, T);

          a_st_one = values.next_double();

          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  cut_st_0[i][j][k][l] = values.next_double();
                  cut_st_0[i][j][k][0] += cut_st_0[i][j][k][l];
                  cut_st_0[0][j][k][l] += cut_st_0[i][j][k][l];
                  cut_st_0[0][j][k][0] += cut_st_0[i][j][k][l];
                }
              }
            }
          }
          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  cut_st_c[i][j][k][l] = values.next_double();
                  cut_st_c[i][j][k][0] += cut_st_c[i][j][k][l];
                  cut_st_c[0][j][k][l] += cut_st_c[i][j][k][l];
                  cut_st_c[0][j][k][0] += cut_st_c[i][j][k][l];
                }
              }
            }
          }
          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  cut_st_lo[i][j][k][l] = values.next_double();
                  cut_st_lo[i][j][k][0] += cut_st_lo[i][j][k][l];
                  cut_st_lo[0][j][k][l] += cut_st_lo[i][j][k][l];
                  cut_st_lo[0][j][k][0] += cut_st_lo[i][j][k][l];
                }
              }
            }
          }
          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  cut_st_hi[i][j][k][l] = values.next_double();
                  cut_st_hi[i][j][k][0] += cut_st_hi[i][j][k][l];
                  cut_st_hi[0][j][k][l] += cut_st_hi[i][j][k][l];
                  cut_st_hi[0][j][k][0] += cut_st_hi[i][j][k][l];
                }
              }
            }
          }
          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  a_st4[i][j][k][l] = values.next_double();
                  a_st4[i][j][k][0] += a_st4[i][j][k][l];
                  a_st4[0][j][k][l] += a_st4[i][j][k][l];
                  a_st4[0][j][k][0] += a_st4[i][j][k][l];
                }
              }
            }
          }

          theta_st4_0_one = values.next_double();

          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                  dtheta_st4_ast[i][j][k][l] = values.next_double();
                  dtheta_st4_ast[i][j][k][0] += dtheta_st4_ast[i][j][k][l];
                  dtheta_st4_ast[0][j][k][l] += dtheta_st4_ast[i][j][k][l];
                  dtheta_st4_ast[0][j][k][0] += dtheta_st4_ast[i][j][k][l];
                }
              }
            }
          }

          a_st5_one = values.next_double();
          theta_st5_0_one = values.next_double();
          dtheta_st5_ast_one = values.next_double();
          a_st6_one = values.next_double();
          theta_st6_0_one = values.next_double();
          dtheta_st6_ast_one = values.next_double();
          a_st1_one = values.next_double();
          cosphi_st1_ast_one = values.next_double();
          a_st2_one = values.next_double();
          cosphi_st2_ast_one = values.next_double();

          break;
        } else continue;
      } catch (std::exception &e) {
        error->one(FLERR, "Problem parsing oxDNA3 potential file: {}", e.what());
      }
    }
    if ((iloc != arg[0]) || (jloc != arg[1]) || (potential_name != "stk"))
      error->one(FLERR, "No corresponding stk potential found in file {} for pair type {} {}",
                 arg[3], arg[0], arg[1]);

    // calculate sequence-averaged parameters for terminal base step j-k
    for (int i = nlo; i <= nhi; i++) {
      for (int j = nlo; j <= nhi; j++) {
        for (int k = nlo; k <= nhi; k++) {
          cut_st_0[i][j][k][0] /= nhi;
          cut_st_c[i][j][k][0] /= nhi;
          cut_st_lo[i][j][k][0] /= nhi;
          cut_st_hi[i][j][k][0] /= nhi;
          a_st4[i][j][k][0] /= nhi;
          dtheta_st4_ast[i][j][k][0] /= nhi;
        }
      }
    }
    for (int j = nlo; j <= nhi; j++) {
      for (int k = nlo; k <= nhi; k++) {
        for (int l = nlo; l <= nhi; l++) {
          cut_st_0[0][j][k][l] /= nhi;
          cut_st_c[0][j][k][l] /= nhi;
          cut_st_lo[0][j][k][l] /= nhi;
          cut_st_hi[0][j][k][l] /= nhi;
          a_st4[0][j][k][l] /= nhi;
          dtheta_st4_ast[0][j][k][l] /= nhi;
        }
      }
    }
    for (int j = nlo; j <= nhi; j++) {
      for (int k = nlo; k <= nhi; k++) {
        cut_st_0[0][j][k][0] /= powint(nhi,2);
        cut_st_c[0][j][k][0] /= powint(nhi,2);
        cut_st_lo[0][j][k][0] /= powint(nhi,2);
        cut_st_hi[0][j][k][0] /= powint(nhi,2);
        a_st4[0][j][k][0] /= powint(nhi,2);
        dtheta_st4_ast[0][j][k][0] /= powint(nhi,2);
      }
    }

  }

  MPI_Bcast(&epsilon_st_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&a_st_one, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&cut_st_0[0][0][0][0], 625, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_st_c[0][0][0][0], 625, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_st_lo[0][0][0][0], 625, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_st_hi[0][0][0][0], 625, MPI_DOUBLE, 0, world);
  MPI_Bcast(&a_st4[0][0][0][0], 625, MPI_DOUBLE, 0, world);

  MPI_Bcast(&theta_st4_0_one, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&dtheta_st4_ast[0][0][0][0], 625, MPI_DOUBLE, 0, world);

  MPI_Bcast(&a_st5_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_st5_0_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_st5_ast_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&a_st6_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_st6_0_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_st6_ast_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&a_st1_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cosphi_st1_ast_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&a_st2_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cosphi_st2_ast_one, 1, MPI_DOUBLE, 0, world);

  // smoothing - determined through continuity and differentiability

  // smoothing strength coincidentally identical for all pairs ij, hence use AAAA tetramer value below
  b_st_lo_one = 2*a_st_one*exp(-a_st_one*(cut_st_lo[1][1][1][1]-cut_st_0[1][1][1][1]))*
      2*a_st_one*exp(-a_st_one*(cut_st_lo[1][1][1][1]-cut_st_0[1][1][1][1]))*
      (1-exp(-a_st_one*(cut_st_lo[1][1][1][1]-cut_st_0[1][1][1][1])))*
      (1-exp(-a_st_one*(cut_st_lo[1][1][1][1]-cut_st_0[1][1][1][1])))/
      (4*((1-exp(-a_st_one*(cut_st_lo[1][1][1][1] -cut_st_0[1][1][1][1])))*
      (1-exp(-a_st_one*(cut_st_lo[1][1][1][1]-cut_st_0[1][1][1][1])))-
      (1-exp(-a_st_one*(cut_st_c[1][1][1][1] -cut_st_0[1][1][1][1])))*
      (1-exp(-a_st_one*(cut_st_c[1][1][1][1]-cut_st_0[1][1][1][1])))));

  // smoothing strength coincidentally identical for all pairs ij, hence use AAAA tetramer value below
  b_st_hi_one = 2*a_st_one*exp(-a_st_one*(cut_st_hi[1][1][1][1]-cut_st_0[1][1][1][1]))*
      2*a_st_one*exp(-a_st_one*(cut_st_hi[1][1][1][1]-cut_st_0[1][1][1][1]))*
      (1-exp(-a_st_one*(cut_st_hi[1][1][1][1]-cut_st_0[1][1][1][1])))*
      (1-exp(-a_st_one*(cut_st_hi[1][1][1][1]-cut_st_0[1][1][1][1])))/
      (4*((1-exp(-a_st_one*(cut_st_hi[1][1][1][1] -cut_st_0[1][1][1][1])))*
      (1-exp(-a_st_one*(cut_st_hi[1][1][1][1]-cut_st_0[1][1][1][1])))-
      (1-exp(-a_st_one*(cut_st_c[1][1][1][1] -cut_st_0[1][1][1][1])))*
      (1-exp(-a_st_one*(cut_st_c[1][1][1][1]-cut_st_0[1][1][1][1])))));

  b_st5_one = a_st5_one*a_st5_one*dtheta_st5_ast_one*dtheta_st5_ast_one/
      (1-a_st5_one*dtheta_st5_ast_one*dtheta_st5_ast_one);
  dtheta_st5_c_one = 1/(a_st5_one*dtheta_st5_ast_one);

  b_st6_one = a_st6_one*a_st6_one*dtheta_st6_ast_one*dtheta_st6_ast_one/
      (1-a_st6_one*dtheta_st6_ast_one*dtheta_st6_ast_one);
  dtheta_st6_c_one = 1/(a_st6_one*dtheta_st6_ast_one);

  b_st1_one = a_st1_one*a_st1_one*cosphi_st1_ast_one*cosphi_st1_ast_one/
      (1-a_st1_one*cosphi_st1_ast_one*cosphi_st1_ast_one);
  cosphi_st1_c_one = 1/(a_st1_one*cosphi_st1_ast_one);

  b_st2_one = a_st2_one*a_st2_one*cosphi_st2_ast_one*cosphi_st2_ast_one/
      (1-a_st2_one*cosphi_st2_ast_one*cosphi_st2_ast_one);
  cosphi_st2_c_one = 1/(a_st2_one*cosphi_st2_ast_one);


  // parameters, uniform or depending on base step
  for (int i = nlo; i <= nhi; i++) {
    for (int j = nlo; j <= nhi; j++) {

      epsilon_st[i][j] = epsilon_st_one * this->eta_st[i-1][j-1];

      a_st[i][j] = a_st_one;
      b_st_lo[i][j] = b_st_lo_one;
      b_st_hi[i][j] = b_st_hi_one;
      theta_st4_0[i][j] = theta_st4_0_one;

      a_st5[i][j] = a_st5_one;
      theta_st5_0[i][j] = theta_st5_0_one;
      dtheta_st5_ast[i][j] = dtheta_st5_ast_one;
      b_st5[i][j] = b_st5_one;
      dtheta_st5_c[i][j] = dtheta_st5_c_one;

      a_st6[i][j] = a_st6_one;
      theta_st6_0[i][j] = theta_st6_0_one;
      dtheta_st6_ast[i][j] = dtheta_st6_ast_one;
      b_st6[i][j] = b_st6_one;
      dtheta_st6_c[i][j] = dtheta_st6_c_one;

      a_st1[i][j] = a_st1_one;
      cosphi_st1_ast[i][j] = cosphi_st1_ast_one;
      b_st1[i][j] = b_st1_one;
      cosphi_st1_c[i][j] = cosphi_st1_c_one;

      a_st2[i][j] = a_st2_one;
      cosphi_st2_ast[i][j] = cosphi_st2_ast_one;
      b_st2[i][j] = b_st2_one;
      cosphi_st2_c[i][j] = cosphi_st2_c_one;

    }
  }

  // parameters depending on tetramer
  for (int i = 0; i <= nhi; i++) { // type 0 for terminal j
    for (int j = nlo; j <= nhi; j++) {
      for (int k = nlo; k <= nhi; k++) {
        for (int l = 0; l <= nhi; l++) { // type 0 for terminal k

          cut_st_lc[i][j][k][l] = cut_st_lo[i][j][k][l]
              - a_st_one*exp(-a_st_one*(cut_st_lo[i][j][k][l]-cut_st_0[i][j][k][l]))*
              (1-exp(-a_st_one*(cut_st_lo[i][j][k][l]-cut_st_0[i][j][k][l])))/b_st_lo_one;

          cut_st_hc[i][j][k][l] = cut_st_hi[i][j][k][l]
              - a_st_one*exp(-a_st_one*(cut_st_hi[i][j][k][l]-cut_st_0[i][j][k][l]))*
              (1-exp(-a_st_one*(cut_st_hi[i][j][k][l]-cut_st_0[i][j][k][l])))/b_st_hi_one;

          cutsq_st_hc[i][j][k][l] = cut_st_hc[i][j][k][l]*cut_st_hc[i][j][k][l];

          tmp = 1 - exp(-(cut_st_c[i][j][k][l]-cut_st_0[i][j][k][l]) * a_st_one);
          shift_st[i][j][k][l] = epsilon_st_one * this->eta_st[j-1][k-1] * tmp * tmp;

          b_st4[i][j][k][l] = a_st4[i][j][k][l]*a_st4[i][j][k][l]*dtheta_st4_ast[i][j][k][l]*
              dtheta_st4_ast[i][j][k][l]/(1-a_st4[i][j][k][l]*dtheta_st4_ast[i][j][k][l]*dtheta_st4_ast[i][j][k][l]);
              dtheta_st4_c[i][j][k][l] = 1/(a_st4[i][j][k][l]*dtheta_st4_ast[i][j][k][l]);

        }
      }
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna3/stk" + utils::errorurl(21));

  // END OF VANILLA CODE DUPLICATION HERE - now we just need to sync the tetramer arrays to device
  // (the non-tetramer Kokkos views are synced within ::init_one)

  this->coeff_set_tetramers_kokkos(narg, arg);
}

namespace LAMMPS_NS {
template class PairOxdna3StkKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairOxdna3StkKokkos<LMPHostType>;
#endif
}
