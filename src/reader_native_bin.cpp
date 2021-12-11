// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "reader_native_bin.h"

#include "error.h"
#include "memory.h"
#include "tokenizer.h"

#include <cstring>
#include <utility>
#include <iostream>

using namespace LAMMPS_NS;

enum{ID,TYPE,X,Y,Z,VX,VY,VZ,Q,IX,IY,IZ,FX,FY,FZ};
enum{UNSET,NOSCALE_NOWRAP,NOSCALE_WRAP,SCALE_NOWRAP,SCALE_WRAP};

/* ---------------------------------------------------------------------- */

ReaderNativeBin::ReaderNativeBin(LAMMPS *lmp) : ReaderNative(lmp)
{
  fieldindex = nullptr;
  buf = new double[maxbuf];
}

/* ---------------------------------------------------------------------- */

ReaderNativeBin::~ReaderNativeBin()
{
  delete [] buf;
  memory->destroy(fieldindex);
}

/* ----------------------------------------------------------------------
   overload the open_file function to open the binary file
------------------------------------------------------------------------- */
void ReaderNativeBin::open_file(const std::string &file)
{
  if (fp != nullptr) close_file();

  if (platform::has_compress_extension(file)) {
    fp = platform::compressed_read(file, true);
    if (!fp) error->one(FLERR,"Cannot open compressed file for reading");
  } else {
    compressed = 0;
    fp = fopen(file.c_str(), "rb");
  }

  if (!fp) error->one(FLERR,"Cannot open file {}: {}", file, utils::getsyserror());
}

/* ----------------------------------------------------------------------
   read and return time stamp from dump file
   if first read reaches end-of-file, return 1 so caller can open next file
   only called by proc 0
------------------------------------------------------------------------- */

int ReaderNativeBin::read_time(bigint &ntimestep)
{

  int endian = 0x0001;
  revision = 0x0001;
  magic_string = nullptr;
  unit_style = nullptr;

  fread(&ntimestep, sizeof(bigint), 1, fp);

  // detect end-of-file
  if (feof(fp)) return 1;

  // detect newer format
  if (ntimestep < 0) {
    // first bigint encodes negative format name length
    bigint magic_string_len = -ntimestep;

    delete[] magic_string;
    magic_string = new char[magic_string_len + 1];
    read_buf(magic_string, sizeof(char), magic_string_len);
    magic_string[magic_string_len] = '\0';

    // read endian flag
    read_buf(&endian, sizeof(int), 1);

    // read revision number
    read_buf(&revision, sizeof(int), 1);

    // read the real ntimestep
    read_buf(&ntimestep, sizeof(bigint), 1);
  }

  return 0;
}

/* ----------------------------------------------------------------------
   skip snapshot from timestamp onward
   only called by proc 0
------------------------------------------------------------------------- */

void ReaderNativeBin::skip()
{
  bigint natoms;
  int triclinic;
  skip_buf(sizeof(bigint));
  read_buf(&triclinic, sizeof(int), 1);
  skip_buf((sizeof(int)+sizeof(double))*6);
  if (triclinic) {
    skip_buf(sizeof(double)*3);
  }
  skip_buf(sizeof(int));

  skip_reading_magic_str();

  // read chunk and skip them

  int nchunk;
  read_buf(&nchunk, sizeof(int), 1);

  int n;
  for (int i = 0; i < nchunk; i++) {
    read_buf(&n, sizeof(int), 1);
    read_double_chunk(n);
  }

  delete[] magic_string;
  delete[] unit_style;
}

void ReaderNativeBin::skip_reading_magic_str()
{
  if (magic_string && revision > 0x0001) {
    int len;
    read_buf(&len, sizeof(int), 1);

    if (len > 0) {
      // has units
      skip_buf(sizeof(char)*len);
    }

    char flag = 0;
    read_buf(&flag, sizeof(char), 1);

    if (flag) {
      skip_buf(sizeof(double));
    }

    read_buf(&len, sizeof(int), 1);
    skip_buf(sizeof(char)*len);
  }
}

/* ----------------------------------------------------------------------
   read remaining header info:
     return natoms
     box bounds, triclinic (inferred), fieldflag (1 if any fields not found),
     xyz flags = UNSET (not a requested field), SCALE/WRAP as in enum
   if fieldflag set:
     match Nfield fields to per-atom column labels
     allocate and set fieldindex = which column each field maps to
     fieldtype = X,VX,IZ etc
     fieldlabel = user-specified label or nullptr if use fieldtype default
   xyz flags = scaleflag+wrapflag if has fieldlabel name,
     else set by x,xs,xu,xsu
   only called by proc 0
------------------------------------------------------------------------- */

bigint ReaderNativeBin::read_header(double box[3][3], int &boxinfo, int &triclinic,
                                 int fieldinfo, int nfield,
                                 int *fieldtype, char **fieldlabel,
                                 int scaleflag, int wrapflag, int &fieldflag,
                                 int &xflag, int &yflag, int &zflag)
{
  bigint natoms;

  read_buf(&natoms, sizeof(bigint), 1);

  boxinfo = 1;
  triclinic = 0;
  box[0][2] = box[1][2] = box[2][2] = 0.0;

  int boundary[3][2];
  read_buf(&triclinic, sizeof(int), 1);
  read_buf(&boundary[0][0], sizeof(int), 6);
  read_buf(box[0], sizeof(double), 2);
  read_buf(box[1], sizeof(double), 2);
  read_buf(box[2], sizeof(double), 2);
  if (triclinic) {
    read_buf(&box[0][2], sizeof(double), 1);
    read_buf(&box[1][2], sizeof(double), 1);
    read_buf(&box[2][2], sizeof(double), 1);
  }

  // exatract column labels and match to requested fields
  read_buf(&size_one, sizeof(int), 1);

  if (!fieldinfo) {
    skip_reading_magic_str();
    return natoms;
  }

  int len = 0;
  char *labelline;

  if (magic_string && revision > 0x0001) {
    // newer format includes units string, columns string
    // and time
    read_buf(&len, sizeof(int), 1);
    labelline = new char[len + 1];

    if (len > 0) {
      // has units
      delete[] unit_style;
      unit_style = new char[len + 1];
      read_buf(unit_style, sizeof(char), len);
      unit_style[len] = '\0';
    }

    char flag = 0;
    read_buf(&flag, sizeof(char), 1);

    if (flag) {
      double time;
      read_buf(&time, sizeof(double), 1);
    }

    read_buf(&len, sizeof(int), 1);
    read_buf(labelline, sizeof(char), len);
    labelline[len] = '\0';
  }

  std::map<std::string, int> labels;
  Tokenizer tokens(labelline);
  int nwords = 0;

  while (tokens.has_next()) {
    labels[tokens.next()] = nwords++;
  }

  if (nwords == 0) {
    return 1;
  }

  match_fields(nfield, xflag, yflag, zflag, fieldtype, fieldlabel, scaleflag, wrapflag, fieldflag, labels);


  return natoms;
}

/* ----------------------------------------------------------------------
   read N atom lines from dump file
   stores appropriate values in fields array
   return 0 if success, 1 if error
   only called by proc 0
------------------------------------------------------------------------- */

void ReaderNativeBin::read_atoms(int n, int nfield, double **fields)
{
  if (feof(fp)) {
    error->one(FLERR,"Unexpected end of dump file");
  }

  int i_atom = 0;
  int nchunk;
  read_buf(&nchunk, sizeof(int), 1);
  for (int i = 0; i < nchunk; i++) {

    read_buf(&n, sizeof(int), 1);

    // read chunk and write as size_one values per line
    read_double_chunk(n);
    n /= size_one;
    int m = 0;
    for (int j = 0; j < n; j++)
    {
      double *words = &buf[m];
      for (int k = 0; k < nfield; k++)
        fields[i_atom][k] = words[fieldindex[k]];
      i_atom += 1;
      m+=size_one;
    }
  }

  delete[] magic_string;
  delete[] unit_style;
}

void ReaderNativeBin::read_buf(void * ptr, size_t size, size_t count)
{
  fread(ptr, size, count, fp);

  // detect end-of-file
  if (feof(fp)) error->one(FLERR,"Unexpected end of dump file");
}

void ReaderNativeBin::read_double_chunk(size_t count)
{
  // extend buffer to fit chunk size
  if (count > maxbuf) {
    if (buf) delete[] buf;
    buf = new double[count];
    maxbuf = count;
  }
  read_buf(buf, sizeof(double), count);
}

void ReaderNativeBin::skip_buf(size_t size)
{
      char tmp[size];
      read_buf(tmp, 1, size);
}