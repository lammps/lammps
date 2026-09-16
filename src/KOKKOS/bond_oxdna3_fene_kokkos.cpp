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

#include "bond_oxdna3_fene_kokkos.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "math_special.h"
#include "potential_file_reader.h"

using namespace LAMMPS_NS;
using namespace MathSpecial;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
BondOxdna3FENEKokkos<DeviceType>::BondOxdna3FENEKokkos(LAMMPS *lmp) : BondOxdnaFENEKokkos<DeviceType>(lmp)
{
  this->oxdnaflag = BondOxdnaFENEKokkos<DeviceType>::EnabledOXDNAFlag::OXDNA2; // oxDNA3 uses same as OXDNA2 here
}

/* ----------------------------------------------------------------------
   set coeffs
   IMPORTANT NOTE ! We entirely code duplicate BondOxdna3Fene::coeff into
   BondOxdna3FENEKokkos::coeff. So any edits made in one needs to manually
   be made to the other ! We did it this way to avoid messy workarounds in
   KOKKOS due to its inheritance structure.
   The vanilla version is in: src/CG-DNA/bond_oxdna3_fene.cpp
------------------------------------------------------------------------- */

template<class DeviceType>
void BondOxdna3FENEKokkos<DeviceType>::coeff(int narg, char **arg)
{
  // Due to class templating of DeviceType, we need this-> on everything. We use local variables
  // so that we can as much as possible just copy-paste the vanilla code (it's cleaner this way also).
  auto *error = this->error;
  auto *atom = this->atom;
  auto *comm = this->comm;
  auto *lmp = this->lmp;
  auto &Delta = this->Delta;
  auto &r0 = this->r0;
  MPI_Comm world = this->world;

  // START OF VANILLA CODE DUPLICATION
  // NOTE: allocate() needs this-> still, and some of the local variables need to come after allocate().
  // Apart from that, this is a vanilla copy-paste.

  if (narg != 2)
    error->all(FLERR, "Incorrect args for bond coefficients in oxdna3/fene, use potential file" + utils::errorurl(21));

  if (!this->allocated) this->allocate();

  // These auto's are not vanilla copy-paste, but they need to come after the allocate() call.
  auto *setflag = this->setflag;
  auto *k = this->k;
  auto &k_k = this->k_k;
  auto &k_r0 = this->k_r0;
  auto &k_Delta = this->k_Delta;

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

  // END OF VANILLA CODE DUPICATION HERE - now we just need to copy the data into the Kokkos views and sync to device

  int m = atom->nbondtypes;
  for (int i = 1; i <= m; i++) {
    k_k.view_host()[i] = k[i];
    for (int n1 = 0; n1 <= n; n1++) {
      for (int n2 = 0; n2 <= n; n2++) {
        for (int n3 = 0; n3 <= n; n3++) {
          for (int n4 = 0; n4 <= n; n4++) {
            k_r0.view_host()(i,n1,n2,n3,n4) = r0[i][n1][n2][n3][n4];
            k_Delta.view_host()(i,n1,n2,n3,n4) = Delta[i][n1][n2][n3][n4];
          }
        }
      }
    }
  }

  k_k.template modify<LMPHostType>();
  k_r0.template modify<LMPHostType>();
  k_Delta.template modify<LMPHostType>();

  // sync to device
  k_k.template sync<DeviceType>();
  k_r0.template sync<DeviceType>();
  k_Delta.template sync<DeviceType>();
}

namespace LAMMPS_NS {
template class BondOxdna3FENEKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class BondOxdna3FENEKokkos<LMPHostType>;
#endif
}
