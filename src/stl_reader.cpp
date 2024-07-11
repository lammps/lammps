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

#include "stl_reader.h"

#include "comm.h"
#include "error.h"
#include "memory.h"
#include "text_file_reader.h"

using namespace LAMMPS_NS;

#define DELTA 16384

/* ---------------------------------------------------------------------- */

STLReader::STLReader(LAMMPS *lmp) : Pointers(lmp)
{
  ntris = maxtris = 0;
  tris = nullptr;
}

/* ---------------------------------------------------------------------- */

STLReader::~STLReader()
{
  memory->destroy(tris);
}

/* ---------------------------------------------------------------------- */

int STLReader::read_file(const char *filename, double **&caller_tris)
{
  int me = comm->me;

  if (me == 0) {

    // open file (text or binary)

    FILE *fp = fopen(filename, "rb");
    if (fp == nullptr) error->one(FLERR, "Cannot open file {}: {}", filename, utils::getsyserror());

    // try reading the file in ASCII format

    TextFileReader reader(fp, "STL mesh");

    try {

      char *line = reader.next_line();
      if (!line || !utils::strmatch(line, "^solid"))
        throw TokenizerException("Invalid STL mesh file format", "");

      line += 6;
      if (utils::strmatch(line, "^binary"))
        throw TokenizerException("Invalid STL mesh file format", "");

      utils::logmesg(lmp, "Reading STL object {} from text file {}\n", utils::trim(line), filename);

      while ((line = reader.next_line())) {

        // next line is facet line with 5 words

        auto values = utils::split_words(line);

        // otherwise next line should be endsolid and are done reading file

        if ((values.size() != 5) || !utils::strmatch(values[0], "^facet")) {
          if (!utils::strmatch(values[0], "^endsolid"))
            throw TokenizerException("Error reading endsolid", "");
          break;
        }

        // ignore normal

        line = reader.next_line(2);
        if (!line || !utils::strmatch(line, "^ *outer *loop"))
          throw TokenizerException("Error reading outer loop", "");

        // three corner points of a single triangle

        if (ntris == maxtris) {
          maxtris += DELTA;
          memory->grow(tris,maxtris,9,"stl_reader:tris");
        }

        for (int k = 0; k < 3; ++k) {
          line = reader.next_line(4);
          values = utils::split_words(line);
          if ((values.size() != 4) || !utils::strmatch(values[0], "^vertex"))
            throw TokenizerException("Error reading vertex", "");

          tris[ntris][3*k+0] = utils::numeric(FLERR, values[1], false, lmp);
          tris[ntris][3*k+1] = utils::numeric(FLERR, values[2], false, lmp);
          tris[ntris][3*k+2] = utils::numeric(FLERR, values[3], false, lmp);
        }
        
        line = reader.next_line(1);
        if (!line || !utils::strmatch(line, "^ *endloop"))
          throw TokenizerException("Error reading endloop", "");
        line = reader.next_line(1);
        if (!line || !utils::strmatch(line, "^ *endfacet"))
          throw TokenizerException("Error reading endfacet", "");

        ntris++;
      }

    // if read as text file failed, try reading as binary file

    } catch (std::exception &e) {

      if (utils::strmatch(e.what(), "^Invalid STL mesh file format")) {
        char title[80];
        float triangle[12];
        uint32_t untri;
        uint16_t attr;
        size_t count;

        rewind(fp);
        count = fread(title, sizeof(char), 80, fp);
        title[79] = '\0';
        count = fread(&untri, sizeof(untri), 1, fp);
        if (count <= 0) {
          error->one(FLERR, "Error reading STL file {}: {}", filename, utils::getsyserror());
        } else {
          utils::logmesg(lmp, "Reading STL object {} from binary file {}\n", utils::trim(title),
                         filename);
        }

        if (untri > MAXSMALLINT)
          error->one(FLERR, "Triangle counts in STL file exceed integer limit");

        ntris = untri;
        memory->create(tris,ntris,9,"stl_reader:tris");

        for (int i = 0; i < ntris; ++i) {
          count = fread(triangle, sizeof(float), 12, fp);
          if (count != 12)
            error->one(FLERR, "Error reading STL file {}: {}", filename, utils::getsyserror());
          count = fread(&attr, sizeof(attr), 1, fp);

          // ignore normal
          
          int m = 0;
          for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
              tris[i][m++] = triangle[3*j + k + 3];
        }

      } else {
        error->all(FLERR, "Error reading STL file {}: {}", filename, e.what());
      }
    }

    if (fp) fclose(fp);
  }

  // MPI broadcast of triangle vertex coords to all procs
  // allow for 9*ntris to exceed max-allowed size of MPI_Bcast()

  MPI_Bcast(&ntris,1,MPI_INT,0,world);
  if (ntris == 0) error->all(FLERR,"STL file has no triangles");
  if (me) memory->create(tris,ntris,9,"stl_reader:tris");

  bigint ntotal = (bigint) ntris * 9;
  if (ntotal < MAXSMALLINT)
    MPI_Bcast(&tris[0][0],9*ntris,MPI_DOUBLE,0,world);
  else {
    double *source = &tris[0][0];
    bigint n = 0;
    while (n < ntotal) {
      int nsize = MIN(MAXSMALLINT,ntotal-n);
      MPI_Bcast(&source[n],nsize,MPI_DOUBLE,0,world);
      n += nsize;
    }
  }

  // return values

  caller_tris = tris;
  return ntris;
}

