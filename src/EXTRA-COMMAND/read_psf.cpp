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

#include "read_psf.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "label_map.h"
#include "memory.h"
#include "platform.h"
#include "reader.h"
#include "text_file_reader.h"

#include <cstring>
#include <string>
#include <vector>

#define MAX_PSF_LABEL_SIZE 8 // psf EXT format

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ReadPsf::ReadPsf(LAMMPS *lmp) : Command(lmp) {}

/* ---------------------------------------------------------------------- */

void ReadPsf::command(int narg, char **arg)
{

  if (domain->box_exist == 0)
    error->all(FLERR,"Read_psf command before simulation box is defined");
  
  // Set flags for CHARMM string properties
  atom->segment_flag = 1;
  atom->residue_flag = 1;
  atom->name_flag = 1;
  
  atom->segment = new std::string[atom->nlocal];
  atom->residue = new std::string[atom->nlocal];
  atom->name = new std::string[atom->nlocal];
  
  if (!atom->labelmapflag) atom->add_label_map();
  LabelMap *lmap = atom->lmap;
  char **lmap_arg;
  memory->create(lmap_arg,3,MAX_PSF_LABEL_SIZE+1,"read_psf:lmap_arg");

  try {
    open(arg[0]);
    utils::logmesg(lmp, "Reading PSF file: {}\n", arg[0]);
    TextFileReader reader(fp, "Protein Structure Format (PSF)");
    reader.skip_line();
    reader.skip_line();
    int ntitle = reader.next_int();

    for( int i=0; i<ntitle ; i++)
      reader.skip_line();

    int natom = reader.next_int();

    for( int i=0; i<natom ; i++) {
      char *line = reader.next_line(9);
      ValueTokenizer values(line);

      // atom tag
      tagint atom_tag = values.next_tagint();
      int atom_index = atom->map(atom_tag);
      
      if( atom_index != -1 ) {
        // atom segment
        atom->segment[atom_index] = values.next_string();
        // skip molecule id
        values.skip(1);
        // residue
        atom->residue[atom_index] = values.next_string();
        // name
        atom->name[atom_index] = values.next_string();
        
        int type_id = atom->type[atom_index];
        strncpy(lmap_arg[0], "atom", 5);
        strncpy(lmap_arg[1],std::to_string(type_id).c_str(),9);
        strncpy(lmap_arg[2],values.next_string().c_str(),9);
        //utils::logmesg(lmp, " *** lmap_arg {} atom_tag {} atom_index {} type_id {} {}\n", lmap_arg[0], atom_tag, atom_index, lmap_arg[1], lmap_arg[2]);
        lmap->modify_lmap(3,lmap_arg);
      }
    }
    
    memory->destroy(lmap_arg);

    // close file
    if (compressed) platform::pclose(fp);
    else fclose(fp);
    fp = nullptr;

  } catch (EOFException &) {
    // reached end of file
    printf("reached EOF\n");
  } catch (std::exception &e) {
    error->one(FLERR, "Error reading psf file: {}", e.what());
  }
}

/* ----------------------------------------------------------------------
   proc 0 opens data file
   test if compressed
------------------------------------------------------------------------- */

void ReadPsf::open(const std::string &file)
{
  if (platform::has_compress_extension(file)) {
    compressed = 1;
    fp = platform::compressed_read(file);
    if (!fp) error->one(FLERR, "Cannot open compressed file {}", file);
  } else {
    compressed = 0;
    fp = fopen(file.c_str(), "r");
    if (!fp) error->one(FLERR, "Cannot open file {}: {}", file, utils::getsyserror());
  }
}
