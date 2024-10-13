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

#include <iostream>

#define MAX_PSF_LABEL_SIZE 8 // psf EXT format

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ReadPsf::ReadPsf(LAMMPS *lmp) :
    Command(lmp)
{

  int flag,cols;
  int index_atom_iarray = atom->find_custom("psf",flag,cols);

  // if atom custom psf doesn't exist, add it
  if( index_atom_iarray == -1 )
    index_atom_iarray = atom->add_custom("psf",0,3,0);

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

  int sendsize;
  int **sendbuf;

  if (comm->me == 0) {
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
      memory->create(sendbuf,natom,4,"read_psf:sendbuf");

      for( int i=0; i<natom ; i++) {
        char *line = reader.next_line(9);
        ValueTokenizer values(line);

        // atom tag
        tagint atom_tag = values.next_tagint();
        int atom_index = atom->map(atom_tag);

        // atom segment
        std::string segment = values.next_string();
        int segment_id = atom->lmap->find_or_add_psf(segment, Atom::SEGMENT);

        // skip molecule id
        values.skip(1);

        // residue
        std::string residue = values.next_string();
        int residue_id = atom->lmap->find_or_add_psf(residue, Atom::RESIDUE);

        // name
        std::string name = values.next_string();
        int name_id = atom->lmap->find_or_add_psf(name, Atom::NAME);

        // determine if this proc owns the atom
        if( atom_index != -1 ) {
          atom_iarray_psf[atom_index][0] = segment_id;
          atom_iarray_psf[atom_index][1] = residue_id;
          atom_iarray_psf[atom_index][2] = name_id;

          // type
          int type_id = atom->type[atom_index];

          strcpy(lmap_arg[0], "atom");
          strcpy(lmap_arg[1],std::to_string(type_id).c_str());
          strcpy(lmap_arg[2],values.next_string().c_str());
          //utils::logmesg(lmp, "read_psf... lmap_arg {} {} {}\n", lmap_arg[0], lmap_arg[1], lmap_arg[2]);
          atom->lmap->modify_lmap(3,lmap_arg);

        } else {
          sendbuf[sendsize][0] = atom_tag;
          sendbuf[sendsize][1] = segment_id;
          sendbuf[sendsize][2] = residue_id;
          sendbuf[sendsize][3] = name_id;
          sendsize++;
        }
      }
      memory->destroy(lmap_arg);

      // close file
      if (compressed)
        platform::pclose(fp);
      else
        fclose(fp);
      fp = nullptr;
      
    } catch (EOFException &) {
      // reached end of file
      printf("reached EOF\n");
    } catch (std::exception &e) {
      error->one(FLERR, "Error reading psf file: {}", e.what());
    }

  }

  MPI_Bcast(&sendsize,1,MPI_INT,0,world);
  if (comm->me > 0) memory->create(sendbuf,sendsize,4,"read_psf:sendbuf");
  MPI_Bcast(&sendbuf[0][0],sendsize*4,MPI_INT,0,world);

  if (comm->me > 0) {

    for( int i=0; i<sendsize ; i++) {

      int atom_index = atom->map(sendbuf[i][0]);

      if( atom_index != -1 ) {
        atom_iarray_psf[atom_index][0] = sendbuf[i][1]; // segment_id
        atom_iarray_psf[atom_index][1] = sendbuf[i][2]; // residue_id
        atom_iarray_psf[atom_index][2] = sendbuf[i][3]; // name_id
      }

    }
  }

  memory->destroy(sendbuf);

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


