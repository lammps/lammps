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
   Contributing author: Mitch Murphy (alphataubio@gmail.com)
------------------------------------------------------------------------- */

#include "read_psf.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "label_map.h"
#include "memory.h"
#include "reader.h"
#include "text_file_reader.h"
#include "potential_file_reader.h"

#include <cstring>

#define MAX_PSF_LABEL_SIZE 8 // psf EXT format

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ReadPsf::ReadPsf(LAMMPS *lmp) :
    Command(lmp)
{

  int flag,cols;
  int index_atom_iarray = atom->find_custom("psf_segment_residue_name",flag,cols);

  // if atom custom psf_segment_residue_name doesn't exist, add it
  if( index_atom_iarray == -1 )
    index_atom_iarray = atom->add_custom("psf_segment_residue_name",0,3,0);

  atom_iarray_psf = atom->iarray[index_atom_iarray];

}

/* ---------------------------------------------------------------------- */

void ReadPsf::command(int narg, char **arg)
{

  if (domain->box_exist == 0)
    error->all(FLERR,"Read_psf command before simulation box is defined");

  if (!atom->labelmapflag)
    atom->add_label_map();

  char **lmap_arg;
  memory->create(lmap_arg,3,MAX_PSF_LABEL_SIZE+1,"read_psf:lmap_arg");

  if (comm->me == 0) {
    try {
      PotentialFileReader reader(lmp, arg[0], "Protein Structure Format (PSF)");
      reader.skip_line();
      reader.skip_line();
      int ntitle = reader.next_int();

      for( int i=0; i<ntitle ; i++)
        reader.skip_line();

      // FIXME: check same number of atoms in psf file as in lammps
      int natom = reader.next_int();

      for( int i=0; i<natom ; i++) {
        char *line = reader.next_line(9);
        ValueTokenizer values(line);

        // atom tag
        tagint atom_tag = values.next_tagint();

        // atom segment
        std::string segment = values.next_string();
        int segment_type = atom->lmap->find_or_add_psf(segment, Atom::SEGMENT);
        atom_iarray_psf[i][0] = segment_type;

        // skip molecule id
        values.skip(1);

        // atom residue
        std::string residue = values.next_string();
        int residue_type = atom->lmap->find_or_add_psf(residue, Atom::RESIDUE);
        atom_iarray_psf[i][1] = residue_type;

        // atom name
        std::string name = values.next_string();
        int name_type = atom->lmap->find_or_add_psf(name, Atom::NAME);
        atom_iarray_psf[i][2] = name_type;

        // atom type
        int atom_type = atom->type[atom_tag-1];
        strcpy(lmap_arg[0], "atom");
        strcpy(lmap_arg[1],std::to_string(atom_type).c_str());
        strcpy(lmap_arg[2],values.next_string().c_str());
        atom->lmap->modify_lmap(3,lmap_arg);
      }
    } catch (EOFException &) {
      // reached end of file
      printf("reached EOF\n");
    } catch (std::exception &e) {
      error->one(FLERR, "Error reading psf file: {}", e.what());
    }
  }

  memory->destroy(lmap_arg);
}

