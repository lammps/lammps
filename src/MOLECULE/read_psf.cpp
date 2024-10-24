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
#include <vector>

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

  int nprocs;
  MPI_Comm_size(world,&nprocs);

  if (domain->box_exist == 0)
    error->all(FLERR,"Read_psf command before simulation box is defined");

  if (!atom->labelmapflag)
    atom->add_label_map();

  LabelMap *lmap = atom->lmap;
  char **lmap_arg;
  memory->create(lmap_arg,3,MAX_PSF_LABEL_SIZE+1,"read_psf:lmap_arg");

  int sendsize = 0;
  int **sendbuf;
  std::vector<std::string> atom_types;

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
        int segment_id = lmap->find_or_add_psf(segment, Atom::SEGMENT);

        // skip molecule id
        values.skip(1);

        // residue
        std::string residue = values.next_string();
        int residue_id = lmap->find_or_add_psf(residue, Atom::RESIDUE);

        // name
        std::string name = values.next_string();
        int name_id = lmap->find_or_add_psf(name, Atom::NAME);

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
          //utils::logmesg(lmp, " *** lmap_arg (me==0) {} {} {}\n", lmap_arg[0], lmap_arg[1], lmap_arg[2]);
          lmap->modify_lmap(3,lmap_arg);
        } else {
          sendbuf[sendsize][0] = atom_tag;
          sendbuf[sendsize][1] = segment_id;
          sendbuf[sendsize][2] = residue_id;
          sendbuf[sendsize][3] = name_id;
          atom_types.push_back(values.next_string());
          sendsize++;
        }
      }

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

  if( nprocs>1 ) {

    // SEND TAG/SEGMENT/RESIDUE/NAME NOT OWNED BY PROC 0 TO OTHER PROCS

    MPI_Bcast(&sendsize,1,MPI_INT,0,world);
    if (comm->me > 0) memory->create(sendbuf,sendsize,4,"read_psf:sendbuf");
    MPI_Bcast(&sendbuf[0][0],sendsize*4,MPI_INT,0,world);

    int recvsize = 0;
    int **recvbuf;
    memory->create(recvbuf,sendsize,2,"read_psf:recvbuf");

    if (comm->me > 0) {
      for( int i=0; i<sendsize ; i++) {
        int atom_index = atom->map(sendbuf[i][0]);
        if( atom_index != -1 ) {
          atom_iarray_psf[atom_index][0] = sendbuf[i][1]; // segment_id
          atom_iarray_psf[atom_index][1] = sendbuf[i][2]; // residue_id
          atom_iarray_psf[atom_index][2] = sendbuf[i][3]; // name_id
          recvbuf[recvsize][0] = i;
          recvbuf[recvsize][1] = atom->type[atom_index]; //type_id
          recvsize++;
        }
      }
    }

    // GET BACK ATOM TYPES NOT OWNED BY PROC 0 AND UPDATE MISSING LMAP TYPES

    int *recvcounts = new int[nprocs];
    int *displs = new int[nprocs];
    recvsize *= 2;
    MPI_Gather( &recvsize, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, world);
    recvcounts[0] = 0;
    displs[0] = 0;

    if (comm->me == 0) {
      for ( int i=1; i<nprocs; i++) displs[i] = displs[i-1]+recvcounts[i-1];
      MPI_Gatherv(MPI_IN_PLACE,0,MPI_INT,&recvbuf[0][0],recvcounts,displs,MPI_INT,0,world);
    } else
      MPI_Gatherv(&recvbuf[0][0],recvsize,MPI_INT,0,0,0,MPI_INT,0,world);

    if (comm->me == 0) {
      for( int j=0; j<sendsize ; j++) {
        int i = recvbuf[j][0];
        if( lmap->find(atom_types[i], Atom::ATOM) == -1 ) {
          int type_id = recvbuf[j][1];
          strcpy(lmap_arg[0], "atom");
          strcpy(lmap_arg[1],std::to_string(type_id).c_str());
          strcpy(lmap_arg[2],atom_types[i].c_str());
          lmap->modify_lmap(3,lmap_arg);
        }
      }
    }

    // SEND COMPLETE LMAP TO OTHER PROCS

    int lmapsizes[4];
    lmapsizes[0] = lmap->natomtypes;
    lmapsizes[1] = lmap->nsegmenttypes;
    lmapsizes[2] = lmap->nresiduetypes;
    lmapsizes[3] = lmap->nnametypes;
    MPI_Bcast(&lmapsizes[0],4,MPI_INT,0,world);
    int lmapsize = lmapsizes[0] + lmapsizes[1] + lmapsizes[2] + lmapsizes[3];
    char **lmapbuf;
    memory->create(lmapbuf,lmapsize,9,"read_psf:lmapbuf");

    if (comm->me == 0) {
      for ( int i=0; i<lmapsizes[0]; i++)
        strcpy(lmapbuf[i],lmap->typelabel[i].c_str());
      for ( int i=0; i<lmapsizes[1]; i++)
        strcpy(lmapbuf[lmapsizes[0]+i],lmap->stypelabel[i].c_str());
      for ( int i=0; i<lmapsizes[2]; i++)
        strcpy(lmapbuf[lmapsizes[0]+lmapsizes[1]+i],lmap->rtypelabel[i].c_str());
      for ( int i=0; i<lmapsizes[3]; i++)
        strcpy(lmapbuf[lmapsizes[0]+lmapsizes[1]+lmapsizes[2]+i],lmap->ntypelabel[i].c_str());
    }

    MPI_Bcast(&lmapbuf[0][0],lmapsize*9,MPI_CHAR,0,world);

    if (comm->me > 0) {
      for ( int i=0; i<lmapsizes[0]; i++)
        lmap->find_or_create(std::string(lmapbuf[i]), lmap->typelabel, lmap->typelabel_map);
      for ( int i=0; i<lmapsizes[1]; i++)
        lmap->find_or_add_psf(std::string(lmapbuf[lmapsizes[0]+i]), Atom::SEGMENT);
      for ( int i=0; i<lmapsizes[2]; i++)
        lmap->find_or_add_psf(std::string(lmapbuf[lmapsizes[0]+lmapsizes[1]+i]), Atom::RESIDUE);
      for ( int i=0; i<lmapsizes[3]; i++)
        lmap->find_or_add_psf(std::string(lmapbuf[lmapsizes[0]+lmapsizes[1]+lmapsizes[2]+i]), Atom::NAME);
    }
    // CLEANUP
    memory->destroy(recvbuf);
    memory->destroy(lmapbuf);
    delete [] recvcounts;
    delete [] displs;
  }
  memory->destroy(sendbuf);
  memory->destroy(lmap_arg);
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
