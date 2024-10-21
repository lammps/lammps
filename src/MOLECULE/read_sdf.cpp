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

#include "read_sdf.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "label_map.h"
#include "memory.h"
#include "molecule.h"
#include "text_file_reader.h"
#include "tokenizer.h"

#include <cstring>

using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

ReadSdf::ReadSdf(LAMMPS *lmp) : Command(lmp) {}

/* ---------------------------------------------------------------------- */

void ReadSdf::command(int narg, char **arg)
{

  atom->molecular = Atom::ATOMIC;

  if (atom->find_molecule(arg[0]) >= 0)
    error->all(FLERR,"Reuse of molecule template ID {}", arg[0]);

  atom->molecules = (Molecule **)
      memory->srealloc(atom->molecules,(atom->nmolecule+1)*sizeof(Molecule *), "read_sdf::molecules");
  Molecule *molecule = new Molecule(lmp);
  atom->molecules[atom->nmolecule] = molecule;
  atom->molecules[atom->nmolecule]->nset = 1;
  atom->nmolecule++;

  molecule->xflag = molecule->typeflag = 1;

  molecule->id = utils::strdup(arg[0]);
  if (!utils::is_id(molecule->id))
    error->all(FLERR, "Molecule template ID must have only alphanumeric or underscore characters");

  if (comm->me == 0) {
    try {
      FILE *fp = fopen(arg[1], "r");
      utils::logmesg(lmp, "Reading SDF file: {}\n", arg[1]);
      TextFileReader reader(fp, "structure-data format (SDF)");
      reader.skip_line();
      reader.skip_line();
      reader.skip_line();
      molecule->natoms = reader.next_int();
      memory->create(molecule->x, molecule->natoms, 3, "molecule:x");
      memory->create(molecule->type, molecule->natoms, "molecule:type");

      for( int i=0; i<molecule->natoms ; i++) {
        char *line = reader.next_line(4);
        ValueTokenizer values(line);

        // x, y, z, element
        molecule->x[i][0] = values.next_double();
        molecule->x[i][1] = values.next_double();
        molecule->x[i][2] = values.next_double();
        std::string element = values.next_string();
        molecule->type[i] = atom->lmap->find(element, Atom::ATOM);
        molecule->ntypes = MAX(molecule->ntypes, molecule->type[i]);

        // FIXME: read partial charges if available

      }

      fclose(fp);
      fp = nullptr;

    } catch (EOFException &) {
      // reached end of file
      printf("reached EOF\n");
    } catch (std::exception &e) {
      error->one(FLERR, "Error reading sdf file: {}", e.what());
    }
  }
}
