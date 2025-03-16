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
   Contributing author: Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "write_pdb.h"

#include "atom.h"
#include "atom_vec.h"
//#include "bond.h"
#include "comm.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "label_map.h"
#include "memory.h"
#include "output.h"

#include <cstring>
#include <string>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   Constructor
------------------------------------------------------------------------- */
WritePDB::WritePDB(LAMMPS *lmp) : Command(lmp) {}

/* ----------------------------------------------------------------------
   Command handler
------------------------------------------------------------------------- */
void WritePDB::command(int narg, char **arg) {
  if (narg < 1) error->all(FLERR, "Write_pdb command requires a filename");

  std::string filename = arg[0];
  write_pdb(filename);
}

/* ----------------------------------------------------------------------
   Function to write PDB file
------------------------------------------------------------------------- */
void WritePDB::write_pdb(const std::string &filename) {
  FILE *fp = fopen(filename.c_str(), "w");
  if (!fp) error->one(FLERR, "Cannot open PDB file for writing");

  fprintf(fp, "HEADER    LAMMPS PDB OUTPUT\n");

  int *molecule = atom->molecule;
  double **x = atom->x;
  int *type = atom->type;
  char **name = atom->name;
  char **segment = atom->segment;
  char **residue = atom->residue;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    fprintf(fp, "ATOM  %5d  %-4s %-4s %-4s   %8.3f%8.3f%8.3f\n",
            i+1, name[i], residue[i], segment[i], x[i][0], x[i][1], x[i][2]);
  }

  fprintf(fp, "END\n");
  fclose(fp);
}
