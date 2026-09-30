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

#include "bond_oxdna3_fene.h"
#include "constants_oxdna.h"
#include "nucleotide_oxdna.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "math_special.h"
#include "potential_file_reader.h"

using namespace LAMMPS_NS;
using namespace MathSpecial;

/* ----------------------------------------------------------------------
   set coeffs
------------------------------------------------------------------------- */
void BondOxdna3Fene::coeff(int narg, char **arg)
{
  if (narg != 2)
    error->all(FLERR, "Incorrect args for bond coefficients in oxdna3/fene, use potential file" + utils::errorurl(21));

  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  int n = atom->ntypes;
  if (n > 4)
    error->all(FLERR, "bond oxdna3/fene does not support more than 4 atom types for A, C, G and T");

  for (int i = 0; i <= n; i++) {
    for (int j = 0; j <= n; j++) {
      for (int k = 0; k <= n; k++) {
        for (int l = 0; l <= n; l++) {
          Delta[ilo][i][j][k][l] = 0.0;
          r0[ilo][i][j][k][l] = 0.0;
        }
      }
    }
  }

  if (comm->me == 0) {    // read values from potential file
    PotentialFileReader reader(lmp, arg[1], "oxdna3 potential", " (fene)");
    reader.set_bufsize(65336);
    char *line;
    std::string iloc, potential_name;

    while ((line = reader.next_line())) {
      try {
        ValueTokenizer values(line);
        iloc = values.next_string();
        potential_name = values.next_string();
        if (iloc == arg[0] && potential_name == "fene") {
          k[ilo] = values.next_double();
          for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
              for (int k = 1; k <= n; k++) {
                for (int l = 1; l <= n; l++) {
                  Delta[ilo][i][j][k][l] = values.next_double();
                  Delta[ilo][i][j][k][0] += Delta[ilo][i][j][k][l];
                  Delta[ilo][0][j][k][l] += Delta[ilo][i][j][k][l];
                  Delta[ilo][0][j][k][0] += Delta[ilo][i][j][k][l];
                }
              }
            }
          }
          for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
              for (int k = 1; k <= n; k++) {
                for (int l = 1; l <= n; l++) {
                  r0[ilo][i][j][k][l] = values.next_double();
                  r0[ilo][i][j][k][0] += r0[ilo][i][j][k][l];
                  r0[ilo][0][j][k][l] += r0[ilo][i][j][k][l];
                  r0[ilo][0][j][k][0] += r0[ilo][i][j][k][l];
                }
              }
            }
          }
          break;
        } else
          continue;
      } catch (std::exception &e) {
        error->one(FLERR, "Problem parsing oxdna3 potential file: {}", e.what());
      }
    }
    if ((iloc != arg[0]) || (potential_name != "fene"))
      error->one(FLERR, "No corresponding fene potential found in file {} for bond type {}", arg[1], arg[0]);

    // calculate sequence-averaged parameters for terminal base step j-k
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
        for (int k = 1; k <= n; k++) {
          Delta[ilo][i][j][k][0] /= n;
          r0[ilo][i][j][k][0] /= n;
        }
      }
    }
    for (int j = 1; j <= n; j++) {
      for (int k = 1; k <= n; k++) {
        for (int l = 1; l <= n; l++) {
          Delta[ilo][0][j][k][l] /= n;
          r0[ilo][0][j][k][l] /= n;
        }
      }
    }
    for (int j = 1; j <= n; j++) {
      for (int k = 1; k <= n; k++) {
        Delta[ilo][0][j][k][0] /= powint(n, 2);
        r0[ilo][0][j][k][0] /= powint(n, 2);
      }
    }
  }

  // communicate parameters for bond type ilo
  MPI_Bcast(&k[ilo], 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Delta[ilo][0][0][0][0], 625, MPI_DOUBLE, 0, world);
  MPI_Bcast(&r0[ilo][0][0][0][0], 625, MPI_DOUBLE, 0, world);

  // set parameters for all other bond types
  int count = 0;
  for (int ib = ilo; ib <= ihi; ib++) {
    k[ib] = k[ilo];
    for (int i = 0; i <= n; i++) {    // type 0 for terminal j
      for (int j = 0; j <= n; j++) {
        for (int k = 0; k <= n; k++) {
          for (int l = 0; l <= n; l++) {    // type 0 for terminal k
            Delta[ib][i][j][k][l] = Delta[ilo][i][j][k][l];
            r0[ib][i][j][k][l] = r0[ilo][i][j][k][l];
          }
        }
      }
    }
    setflag[ib] = 1;
    count++;
  }

  if (count == 0)
    error->all(FLERR, "Incorrect args for bond coefficients in oxdna3/fene" + utils::errorurl(21));
}
