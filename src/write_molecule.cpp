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
   Contributing author:  Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "write_molecule.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "json.h"
#include "molecule.h"
#include "safe_pointers.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void WriteMolecule::command(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "write_molecule", error);

  auto imol = atom->find_molecule(arg[0]);
  if (imol == -1)
    error->all(FLERR, Error::ARGZERO, "Molecule template ID {} for write_molecule does not exist",
               arg[0]);
  auto *mol = atom->molecules[imol];
  int nset = mol->nset;
  if (nset < 1)
    error->all(FLERR, Error::ARGZERO, "Molecule template {} has invalid data (nset = {})", nset);

  if ((nset > 1) && !strchr(arg[1], '%'))
    error->all(FLERR, 1,
               "Output filename must contain a '%' for molecule template {} with "
               "multiple data sets",
               arg[0]);

  // from here on, we only need MPI rank 0
  if (comm->me) return;

  if (utils::strmatch(arg[1], "\\.json$")) {
    for (int i = 0; i < nset; ++i) {
      mol = atom->molecules[imol];
      std::string filename = arg[1];
      auto idx = filename.rfind('%');
      if (nset > 1) filename.replace(idx, 1, std::to_string(i + 1));
      utils::logmesg(lmp, "Writing molecule {} to file {} in JSON format\n", arg[0], filename);
      SafeFilePtr fp = fopen(filename.c_str(), "w");
      if (fp == nullptr)
        error->one(FLERR, 1, "Could not open file {} for writing: {}", filename,
                   utils::getsyserror());
      auto moldata = mol->to_json();
      fputs(moldata.dump(2).c_str(), fp);
      fputs("\n", fp);
      ++imol;
    }
  } else {
    for (int i = 0; i < nset; ++i) {
      mol = atom->molecules[imol];
      std::string filename = arg[1];
      auto idx = filename.rfind('%');
      if (nset > 1) filename.replace(idx, 1, std::to_string(i + 1));
      utils::logmesg(lmp, "Writing molecule {} to file {} in native format\n", arg[0], filename);
      SafeFilePtr fp = fopen(filename.c_str(), "w");
      if (fp == nullptr)
        error->one(FLERR, 1, "Could not open file {} for writing: {}", filename,
                   utils::getsyserror());
      utils::print(fp, "# MOLECULE {}, units = {}, set {} of {}, {}\n\n", arg[0],
                   update->unit_style, i + 1, nset, mol->title);
      mol->print(fp);
      fputs("\n", fp);
      ++imol;
    }
  }
}
