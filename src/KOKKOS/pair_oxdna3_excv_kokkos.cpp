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

#include "pair_oxdna3_excv_kokkos.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "potential_file_reader.h"
#include "math_special.h"

#include <cmath>
#include <cassert>

using namespace LAMMPS_NS;
using namespace MathSpecial;

/* ----------------------------------------------------------------------
   IMPORTANT NOTE ! We entirely code duplicate PairOxdna3ExcvKokkos::coeff
   into PairOxdna3Excv::coeff. So any edits made in one need to manually be
   made to the other !
   The vanilla version is in: src/CG-DNA/pair_oxdna3_excv.cpp
------------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdna3ExcvKokkos<DeviceType>::coeff(int narg, char **arg)
{
  // Due to class templating of DeviceType, we need this-> on everything. We use local variables
  // so that we can as much as possible just copy-paste the vanilla code (it's cleaner this way also).
  auto *error = this->error;
  auto *atom = this->atom;
  auto *comm = this->comm;
  auto *lmp = this->lmp;
  MPI_Comm world = this->world;

  int count;

  if (narg != 3)
    error->all(FLERR,"Incorrect args for pair coefficients in oxdna3/excv, use potential file" + utils::errorurl(21));

  if (!this->allocated) this->allocate();

  // NOTE: allocate() needs this-> still, but otherwise this is a direct copy and paste from the
  // vanilla code. These pointer aliases must be taken AFTER allocate(), since allocate() is what
  // assigns the underlying member pointers.

  auto *setflag = this->setflag;

  auto *epsilon_bkbk = this->epsilon_bkbk;
  auto *sigma_bkbk = this->sigma_bkbk;
  auto *cut_bkbk_ast = this->cut_bkbk_ast;
  auto *b_bkbk = this->b_bkbk;
  auto *cut_bkbk_c = this->cut_bkbk_c;
  auto *lj1_bkbk = this->lj1_bkbk;
  auto *lj2_bkbk = this->lj2_bkbk;
  auto *cutsq_bkbk_ast = this->cutsq_bkbk_ast;
  auto *cutsq_bkbk_c = this->cutsq_bkbk_c;

  auto *epsilon_bkbs = this->epsilon_bkbs;
  auto *sigma_bkbs = this->sigma_bkbs;
  auto *cut_bkbs_ast = this->cut_bkbs_ast;
  auto *b_bkbs = this->b_bkbs;
  auto *cut_bkbs_c = this->cut_bkbs_c;
  auto *lj1_bkbs = this->lj1_bkbs;
  auto *lj2_bkbs = this->lj2_bkbs;
  auto *cutsq_bkbs_ast = this->cutsq_bkbs_ast;
  auto *cutsq_bkbs_c = this->cutsq_bkbs_c;

  auto *epsilon_bsbs = this->epsilon_bsbs;
  auto *sigma_bsbs = this->sigma_bsbs;
  auto *cut_bsbs_ast = this->cut_bsbs_ast;
  auto *b_bsbs = this->b_bsbs;
  auto *cut_bsbs_c = this->cut_bsbs_c;
  auto *lj1_bsbs = this->lj1_bsbs;
  auto *lj2_bsbs = this->lj2_bsbs;
  auto *cutsq_bsbs_ast = this->cutsq_bsbs_ast;
  auto *cutsq_bsbs_c = this->cutsq_bsbs_c;

  auto *sigma4_bsbs = this->sigma4_bsbs;
  auto *cut4_bsbs_ast = this->cut4_bsbs_ast;
  auto *cut4sq_bsbs_ast = this->cut4sq_bsbs_ast;
  auto *lj14_bsbs = this->lj14_bsbs;
  auto *lj24_bsbs = this->lj24_bsbs;
  auto *b4_bsbs = this->b4_bsbs;
  auto *cut4_bsbs_c = this->cut4_bsbs_c;
  auto *cut4sq_bsbs_c = this->cut4sq_bsbs_c;

  // START OF VANILLA CODE DUPLICATION

  int ilo,ihi,jlo,jhi,nlo,nhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  assert((ilo == jlo) & (ihi == jhi));
  nlo = ilo;
  nhi = ihi;

  if (nhi > 4) error->all(FLERR, "pair oxdna3/excv does not support more than 4 atom types for A, C, G and T");

  double epsilon_bkbk_one, sigma_bkbk_one;
  double cut_bkbk_ast_one, cut_bkbk_c_one, b_bkbk_one;

  double epsilon_bkbs_one, sigma_bkbs_one;
  double cut_bkbs_ast_one, cut_bkbs_c_one, b_bkbs_one;

  double epsilon_bsbs_one, sigma_bsbs_one;
  double cut_bsbs_ast_one, cut_bsbs_c_one, b_bsbs_one;

  for (int i = 0; i <= nhi; i++) { // type 0 for terminal j
    for (int j = 0; j <= nhi; j++) {
      for (int k = 0; k <= nhi; k++) {
        for (int l = 0; l <= nhi; l++) { // type 0 for terminal k
          sigma4_bsbs[i][j][k][l] = 0.0;
          cut4_bsbs_ast[i][j][k][l] = 0.0;
        }
      }
    }
  }

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, arg[2], "oxdna3 potential", " (excv)");
    reader.set_bufsize(65336);
    char * line;
    std::string iloc, jloc, potential_name;

    while ((line = reader.next_line())) {
      try {
        ValueTokenizer values(line);
        iloc = values.next_string();
        jloc = values.next_string();
        potential_name = values.next_string();
        if (iloc == arg[0] && jloc == arg[1] && potential_name == "excv") {
          // Excluded volume interaction
          // LJ backbone-backbone parameters
          epsilon_bkbk_one = values.next_double();
          sigma_bkbk_one = values.next_double();
          cut_bkbk_ast_one = values.next_double();

          // LJ backbone-base parameters
          epsilon_bkbs_one = values.next_double();
          sigma_bkbs_one = values.next_double();
          cut_bkbs_ast_one = values.next_double();

          // LJ base-base parameters
          epsilon_bsbs_one = values.next_double();
          sigma_bsbs_one = values.next_double();
          cut_bsbs_ast_one = values.next_double();

          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                sigma4_bsbs[i][j][k][l] = values.next_double();
                sigma4_bsbs[i][j][k][0] += sigma4_bsbs[i][j][k][l];
                sigma4_bsbs[0][j][k][l] += sigma4_bsbs[i][j][k][l];
                sigma4_bsbs[0][j][k][0] += sigma4_bsbs[i][j][k][l];
                }
              }
            }
          }

          for (int i = nlo; i <= nhi; i++) {
            for (int j = nlo; j <= nhi; j++) {
              for (int k = nlo; k <= nhi; k++) {
                for (int l = nlo; l <= nhi; l++) {
                cut4_bsbs_ast[i][j][k][l] = values.next_double();
                cut4_bsbs_ast[i][j][k][0] += cut4_bsbs_ast[i][j][k][l];
                cut4_bsbs_ast[0][j][k][l] += cut4_bsbs_ast[i][j][k][l];
                cut4_bsbs_ast[0][j][k][0] += cut4_bsbs_ast[i][j][k][l];
                }
              }
            }
          }

          break;
        } else continue;
      } catch (std::exception &e) {
        error->one(FLERR, "Problem parsing oxdna3 potential file: {}", e.what());
      }
    }
    if ((iloc != arg[0]) || (jloc != arg[1]) || (potential_name != "excv"))
      error->one(FLERR, "No corresponding excv potential found in file {} for pair type {} {}",
                 arg[2], arg[0], arg[1]);


    // calculate sequence-averaged parameters for terminal base step j-k
    for (int i = nlo; i <= nhi; i++) {
      for (int j = nlo; j <= nhi; j++) {
        for (int k = nlo; k <= nhi; k++) {
          sigma4_bsbs[i][j][k][0] /= nhi;
          cut4_bsbs_ast[i][j][k][0] /= nhi;
        }
      }
    }
    for (int j = nlo; j <= nhi; j++) {
      for (int k = nlo; k <= nhi; k++) {
        for (int l = nlo; l <= nhi; l++) {
          sigma4_bsbs[0][j][k][l] /= nhi;
          cut4_bsbs_ast[0][j][k][l] /= nhi;

        }
      }
    }
    for (int j = nlo; j <= nhi; j++) {
      for (int k = nlo; k <= nhi; k++) {
        sigma4_bsbs[0][j][k][0] /= powint(nhi,2);
        cut4_bsbs_ast[0][j][k][0] /= powint(nhi,2);
      }
    }

  }

  // The 3x3 MPI broadcasts below are indifferent to the version of oxDNA that is simulated at
  // compile/runtime in the KOKKOS build/case.
  MPI_Bcast(&epsilon_bkbk_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&sigma_bkbk_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_bkbk_ast_one, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&epsilon_bkbs_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&sigma_bkbs_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_bkbs_ast_one, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&epsilon_bsbs_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&sigma_bsbs_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_bsbs_ast_one, 1, MPI_DOUBLE, 0, world);

  // But for the tetramers, we put in the prefix
  MPI_Bcast(&sigma4_bsbs[0][0][0][0], 625, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut4_bsbs_ast[0][0][0][0], 625, MPI_DOUBLE, 0, world);

  // backbone-backbone
  count = 0;

  // smoothing - determined through continuity and differentiability
  b_bkbk_one = 4.0/sigma_bkbk_one
      *(6.0*pow(sigma_bkbk_one/cut_bkbk_ast_one,7)-12.0*pow(sigma_bkbk_one/cut_bkbk_ast_one,13))
      *4.0/sigma_bkbk_one*(6.0*pow(sigma_bkbk_one/cut_bkbk_ast_one,7)-12.0*pow(sigma_bkbk_one/cut_bkbk_ast_one,13))
      /4.0/(4.0*(pow(sigma_bkbk_one/cut_bkbk_ast_one,12)-pow(sigma_bkbk_one/cut_bkbk_ast_one,6)));

  cut_bkbk_c_one = cut_bkbk_ast_one
      - 2.0*4.0*(pow(sigma_bkbk_one/cut_bkbk_ast_one,12)-pow(sigma_bkbk_one/cut_bkbk_ast_one,6))
      /(4.0/sigma_bkbk_one*(6.0*pow(sigma_bkbk_one/cut_bkbk_ast_one,7)-12.0*pow(sigma_bkbk_one/cut_bkbk_ast_one,13)));

  // backbone-backbone parameters depending on base step
  for (int i = nlo; i <= nhi; i++) {
    for (int j = nlo; j <= nhi; j++) {
      epsilon_bkbk[i][j] = epsilon_bkbk_one;
      sigma_bkbk[i][j] = sigma_bkbk_one;
      cut_bkbk_ast[i][j] = cut_bkbk_ast_one;
      b_bkbk[i][j] = b_bkbk_one;
      cut_bkbk_c[i][j] = cut_bkbk_c_one;
      lj1_bkbk[i][j] = 4.0 * epsilon_bkbk[i][j] * pow(sigma_bkbk[i][j],12.0);
      lj2_bkbk[i][j] = 4.0 * epsilon_bkbk[i][j] * pow(sigma_bkbk[i][j],6.0);
      cutsq_bkbk_ast[i][j] = cut_bkbk_ast[i][j]*cut_bkbk_ast[i][j];
      cutsq_bkbk_c[i][j]  = cut_bkbk_c[i][j]*cut_bkbk_c[i][j];
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/excv" + utils::errorurl(21));

  // backbone-base
  count = 0;

  // smoothing - determined through continuity and differentiability
  b_bkbs_one = 4.0/sigma_bkbs_one
      *(6.0*pow(sigma_bkbs_one/cut_bkbs_ast_one,7)-12.0*pow(sigma_bkbs_one/cut_bkbs_ast_one,13))
      *4.0/sigma_bkbs_one*(6.0*pow(sigma_bkbs_one/cut_bkbs_ast_one,7)-12.0*pow(sigma_bkbs_one/cut_bkbs_ast_one,13))
      /4.0/(4.0*(pow(sigma_bkbs_one/cut_bkbs_ast_one,12)-pow(sigma_bkbs_one/cut_bkbs_ast_one,6)));

  cut_bkbs_c_one = cut_bkbs_ast_one
      - 2.0*4.0*(pow(sigma_bkbs_one/cut_bkbs_ast_one,12)-pow(sigma_bkbs_one/cut_bkbs_ast_one,6))
      /(4.0/sigma_bkbs_one*(6.0*pow(sigma_bkbs_one/cut_bkbs_ast_one,7)-12.0*pow(sigma_bkbs_one/cut_bkbs_ast_one,13)));

  // backbone-base parameters depending on base step
  for (int i = nlo; i <= nhi; i++) {
    for (int j = nlo; j <= nhi; j++) {
      epsilon_bkbs[i][j] = epsilon_bkbs_one;
      sigma_bkbs[i][j] = sigma_bkbs_one;
      cut_bkbs_ast[i][j] = cut_bkbs_ast_one;
      b_bkbs[i][j] = b_bkbs_one;
      cut_bkbs_c[i][j] = cut_bkbs_c_one;
      lj1_bkbs[i][j] = 4.0 * epsilon_bkbs[i][j] * pow(sigma_bkbs[i][j],12.0);
      lj2_bkbs[i][j] = 4.0 * epsilon_bkbs[i][j] * pow(sigma_bkbs[i][j],6.0);
      cutsq_bkbs_ast[i][j] = cut_bkbs_ast[i][j]*cut_bkbs_ast[i][j];
      cutsq_bkbs_c[i][j]  = cut_bkbs_c[i][j]*cut_bkbs_c[i][j];
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/excv" + utils::errorurl(21));

  // base-base
  count = 0;

  // smoothing - determined through continuity and differentiability
  b_bsbs_one = 4.0/sigma_bsbs_one
      *(6.0*pow(sigma_bsbs_one/cut_bsbs_ast_one,7)-12.0*pow(sigma_bsbs_one/cut_bsbs_ast_one,13))
      *4.0/sigma_bsbs_one*(6.0*pow(sigma_bsbs_one/cut_bsbs_ast_one,7)-12.0*pow(sigma_bsbs_one/cut_bsbs_ast_one,13))
      /4.0/(4.0*(pow(sigma_bsbs_one/cut_bsbs_ast_one,12)-pow(sigma_bsbs_one/cut_bsbs_ast_one,6)));

  cut_bsbs_c_one = cut_bsbs_ast_one
      - 2.0*4.0*(pow(sigma_bsbs_one/cut_bsbs_ast_one,12)-pow(sigma_bsbs_one/cut_bsbs_ast_one,6))
      /(4.0/sigma_bsbs_one*(6.0*pow(sigma_bsbs_one/cut_bsbs_ast_one,7)-12.0*pow(sigma_bsbs_one/cut_bsbs_ast_one,13)));

  // base-base parameters depending on base step
  for (int i = nlo; i <= nhi; i++) {
    for (int j = nlo; j <= nhi; j++) {
      epsilon_bsbs[i][j] = epsilon_bsbs_one;
      sigma_bsbs[i][j] = sigma_bsbs_one;
      cut_bsbs_ast[i][j] = cut_bsbs_ast_one;
      b_bsbs[i][j] = b_bsbs_one;
      cut_bsbs_c[i][j] = cut_bsbs_c_one;
      lj1_bsbs[i][j] = 4.0 * epsilon_bsbs[i][j] * pow(sigma_bsbs[i][j],12.0);
      lj2_bsbs[i][j] = 4.0 * epsilon_bsbs[i][j] * pow(sigma_bsbs[i][j],6.0);
      cutsq_bsbs_ast[i][j] = cut_bsbs_ast[i][j]*cut_bsbs_ast[i][j];
      cutsq_bsbs_c[i][j]  = cut_bsbs_c[i][j]*cut_bsbs_c[i][j];
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/excv");

  // base-base parameters depending on tetramer
  count = 0;

  for (int i = 0; i <= nhi; i++) { // type 0 for terminal j
    for (int j = nlo; j <= nhi; j++) {
      for (int k = nlo; k <= nhi; k++) {
        for (int l = 0; l <= nhi; l++) { // type 0 for terminal k

          b4_bsbs[i][j][k][l] = 4.0/sigma4_bsbs[i][j][k][l]
              *(6.0*pow(sigma4_bsbs[i][j][k][l]/cut4_bsbs_ast[i][j][k][l],7)
              -12.0*pow(sigma4_bsbs[i][j][k][l]/cut4_bsbs_ast[i][j][k][l],13))
              *4.0/sigma4_bsbs[i][j][k][l]*(6.0*pow(sigma4_bsbs[i][j][k][l]/cut4_bsbs_ast[i][j][k][l],7)
              -12.0*pow(sigma4_bsbs[i][j][k][l]/cut4_bsbs_ast[i][j][k][l],13))
              /4.0/(4.0*(pow(sigma4_bsbs[i][j][k][l]/cut4_bsbs_ast[i][j][k][l],12)
              -pow(sigma4_bsbs[i][j][k][l]/cut4_bsbs_ast[i][j][k][l],6)));

          cut4_bsbs_c[i][j][k][l] = cut4_bsbs_ast[i][j][k][l]
              - 2.0*4.0*(pow(sigma4_bsbs[i][j][k][l]/cut4_bsbs_ast[i][j][k][l],12)
              -pow(sigma4_bsbs[i][j][k][l]/cut4_bsbs_ast[i][j][k][l],6))
              /(4.0/sigma4_bsbs[i][j][k][l]*(6.0*pow(sigma4_bsbs[i][j][k][l]/cut4_bsbs_ast[i][j][k][l],7)
              -12.0*pow(sigma4_bsbs[i][j][k][l]/cut4_bsbs_ast[i][j][k][l],13)));

          cut4sq_bsbs_ast[i][j][k][l] = cut4_bsbs_ast[i][j][k][l]*cut4_bsbs_ast[i][j][k][l];
          cut4sq_bsbs_c[i][j][k][l]  = cut4_bsbs_c[i][j][k][l]*cut4_bsbs_c[i][j][k][l];
          lj14_bsbs[i][j][k][l] = 4.0 * epsilon_bsbs[j][k] * pow(sigma4_bsbs[i][j][k][l],12.0);
          lj24_bsbs[i][j][k][l] = 4.0 * epsilon_bsbs[j][k] * pow(sigma4_bsbs[i][j][k][l],6.0);
          count++;
       }
      }
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/excv");

  // END OF VANILLA CODE DUPLICATION HERE - now we just need to sync the tetramer arrays to device
  // (the non-tetramer Kokkos views are synced within ::init_one)

  this->coeff_set_tetramers_kokkos(narg, arg);
}

namespace LAMMPS_NS {
template class PairOxdna3ExcvKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairOxdna3ExcvKokkos<LMPHostType>;
#endif
}
