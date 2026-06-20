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
#include "tokenizer.h"

#include <cstdint>

using namespace LAMMPS_NS;

// byte layout of a binary STL file:
//   80-byte header + 4-byte (uint32) triangle count, then per triangle
//   12 floats (normal + 3 vertices) + a 2-byte attribute count

static constexpr long STL_BIN_HEADER = 80 + sizeof(uint32_t);              // 84
static constexpr long STL_BIN_PER_TRI = 12 * sizeof(float) + sizeof(uint16_t);  // 50

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

/* ----------------------------------------------------------------------
   read triangles on rank 0, broadcast vertex coords to all ranks
   caller_tris is filled with 9 doubles per triangle (3 vertices, no normal)
   allow for 9*ntris to exceed the max allowed size of a single MPI_Bcast()
------------------------------------------------------------------------- */

int STLReader::read_file(const char *filename, double **&caller_tris)
{
  int me = comm->me;

  if (me == 0) {
    try {
      std::string title;
      std::vector<Triangle> parsed = parse(filename, &title);
      ntris = (int) parsed.size();
      maxtris = ntris;
      utils::logmesg(lmp, "Reading STL object {} with {} triangles from file {}\n",
                     title.empty() ? "(unnamed)" : title, ntris, filename);

      memory->create(tris, (ntris > 0) ? ntris : 1, 9, "stl_reader:tris");
      for (int i = 0; i < ntris; i++) {
        int m = 0;
        for (int j = 0; j < 3; j++)
          for (int k = 0; k < 3; k++)
            tris[i][m++] = parsed[i].vert[j][k];
      }
    } catch (std::exception &e) {
      error->one(FLERR, "{}", e.what());
    }
  }

  MPI_Bcast(&ntris, 1, MPI_INT, 0, world);
  if (ntris == 0) error->all(FLERR, "STL file {} has no triangles", filename);
  if (me) memory->create(tris, ntris, 9, "stl_reader:tris");

  bigint ntotal = (bigint) ntris * 9;
  if (ntotal < MAXSMALLINT)
    MPI_Bcast(&tris[0][0], 9*ntris, MPI_DOUBLE, 0, world);
  else {
    double *source = &tris[0][0];
    bigint n = 0;
    while (n < ntotal) {
      int nsize = MIN(MAXSMALLINT, ntotal-n);
      MPI_Bcast(&source[n], nsize, MPI_DOUBLE, 0, world);
      n += nsize;
    }
  }

  caller_tris = tris;
  return ntris;
}

/* ----------------------------------------------------------------------
   stand-alone STL parser: auto-detect ASCII vs binary and return triangles
   a binary STL file is exactly STL_BIN_HEADER + STL_BIN_PER_TRI*ntri bytes,
   which is the robust way to distinguish it from an ASCII file (whose header
   may legally also start with the word "solid")
------------------------------------------------------------------------- */

std::vector<STLReader::Triangle> STLReader::parse(const std::string &filename, std::string *title_out)
{
  FILE *fp = fopen(filename.c_str(), "rb");
  if (!fp)
    throw STLReaderException(fmt::format("Cannot open STL file {}: {}", filename,
                                         utils::getsyserror()));

  // determine the file size and the triangle count claimed by a binary header

  bool is_binary = false;
  if ((fseek(fp, 0, SEEK_END) == 0)) {
    long filesize = ftell(fp);
    if (filesize >= STL_BIN_HEADER) {
      uint32_t ntri_claim = 0;
      if ((fseek(fp, 80, SEEK_SET) == 0) && (fread(&ntri_claim, sizeof(ntri_claim), 1, fp) == 1)) {
        bigint expected = (bigint) STL_BIN_HEADER + (bigint) ntri_claim * STL_BIN_PER_TRI;
        if (expected == (bigint) filesize) is_binary = true;
      }
    }
  }
  rewind(fp);

  std::string title;
  std::vector<Triangle> triangles;

  try {
    if (is_binary) {
      triangles = parse_binary(fp, title);
    } else {
      fclose(fp);
      fp = nullptr;
      triangles = parse_text(filename, title);
    }
  } catch (...) {
    if (fp) fclose(fp);
    throw;
  }
  if (fp) fclose(fp);

  if (title_out) *title_out = title;
  return triangles;
}

/* ----------------------------------------------------------------------
   parse an ASCII STL file using the LAMMPS text file reader and tokenizer
------------------------------------------------------------------------- */

std::vector<STLReader::Triangle> STLReader::parse_text(const std::string &filename, std::string &title)
{
  TextFileReader reader(filename, "STL mesh");

  char *line = reader.next_line();
  if (!line || !utils::strmatch(line, "^ *solid"))
    throw STLReaderException(fmt::format("File {} is not a valid ASCII STL file", filename));

  // solid name may be empty; use a std::string so there is no risk of running
  // past the end of the buffer (as bare pointer arithmetic on "solid" would)

  std::string header(line);
  auto pos = header.find("solid");
  title = utils::trim(header.substr(pos + 5));

  std::vector<Triangle> triangles;

  try {
    while ((line = reader.next_line())) {
      auto words = utils::split_words(line);
      if (words.empty()) continue;
      if (utils::strmatch(words[0], "^endsolid")) break;
      if (!utils::strmatch(words[0], "^facet"))
        throw STLReaderException(
            fmt::format("Expected 'facet' or 'endsolid' in STL file, got: {}", utils::trim(line)));

      Triangle tri;
      tri.normal[0] = tri.normal[1] = tri.normal[2] = 0.0;

      // facet line is "facet normal nx ny nz"; the normal is optional and
      // tolerated to be absent or unparsable (it is recomputed when needed)

      if ((words.size() >= 5) && utils::strmatch(words[1], "^normal")) {
        try {
          ValueTokenizer nvals(fmt::format("{} {} {}", words[2], words[3], words[4]));
          tri.normal[0] = nvals.next_double();
          tri.normal[1] = nvals.next_double();
          tri.normal[2] = nvals.next_double();
        } catch (TokenizerException &) {
          tri.normal[0] = tri.normal[1] = tri.normal[2] = 0.0;
        }
      }

      line = reader.next_line();
      if (!line || !utils::strmatch(line, "^ *outer *loop"))
        throw STLReaderException("Error reading 'outer loop' in STL file");

      for (int k = 0; k < 3; k++) {
        ValueTokenizer values = reader.next_values(4);
        if (values.next_string() != "vertex")
          throw STLReaderException(fmt::format("Error reading vertex {} of facet in STL file", k+1));
        tri.vert[k][0] = values.next_double();
        tri.vert[k][1] = values.next_double();
        tri.vert[k][2] = values.next_double();
      }

      line = reader.next_line();
      if (!line || !utils::strmatch(line, "^ *endloop"))
        throw STLReaderException("Error reading 'endloop' in STL file");
      line = reader.next_line();
      if (!line || !utils::strmatch(line, "^ *endfacet"))
        throw STLReaderException("Error reading 'endfacet' in STL file");

      triangles.push_back(tri);
    }
  } catch (TokenizerException &e) {
    throw STLReaderException(fmt::format("Error parsing STL file {}: {}", filename, e.what()));
  }

  return triangles;
}

/* ----------------------------------------------------------------------
   parse a binary STL file: 80-byte header, uint32 triangle count, then
   per triangle 12 little-endian floats (normal + 3 vertices) + 2-byte attr
------------------------------------------------------------------------- */

std::vector<STLReader::Triangle> STLReader::parse_binary(FILE *fp, std::string &title)
{
  rewind(fp);

  char head[80];
  if (fread(head, 1, 80, fp) != 80)
    throw STLReaderException("Unexpected end of binary STL file while reading header");

  // the header is a fixed 80-byte field that need not be null terminated

  std::size_t len = 0;
  while ((len < sizeof(head)) && (head[len] != '\0')) ++len;
  title = utils::trim(std::string(head, len));

  uint32_t ntri = 0;
  if (fread(&ntri, sizeof(ntri), 1, fp) != 1)
    throw STLReaderException("Unexpected end of binary STL file while reading triangle count");
  if (ntri > (uint32_t) MAXSMALLINT)
    throw STLReaderException("Number of triangles in STL file exceeds integer limit");

  std::vector<Triangle> triangles;
  triangles.reserve(ntri);

  float buf[12];
  uint16_t attr;
  for (uint32_t i = 0; i < ntri; i++) {
    if (fread(buf, sizeof(float), 12, fp) != 12)
      throw STLReaderException(
          fmt::format("Unexpected end of binary STL file at triangle {} of {}", i+1, ntri));
    if (fread(&attr, sizeof(attr), 1, fp) != 1)
      throw STLReaderException(
          fmt::format("Unexpected end of binary STL file reading attributes of triangle {}", i+1));

    Triangle tri;
    for (int k = 0; k < 3; k++) tri.normal[k] = buf[k];
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        tri.vert[j][k] = buf[3 + 3*j + k];
    triangles.push_back(tri);
  }

  return triangles;
}
