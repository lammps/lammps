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

using namespace LAMMPS_NS;

#define MAXLINE 1024        // max line length in dump file

// also in read_dump.cpp

enum{ID,TYPE,X,Y,Z,VX,VY,VZ,Q,IX,IY,IZ,FX,FY,FZ};
enum{UNSET,NOSCALE_NOWRAP,NOSCALE_WRAP,SCALE_NOWRAP,SCALE_WRAP};

/* ---------------------------------------------------------------------- */

ReaderNativeBin::ReaderNativeBin(LAMMPS *lmp) : ReaderNative(lmp)
{
  line = new char[MAXLINE];
  fieldindex = nullptr;
}

/* ---------------------------------------------------------------------- */

ReaderNativeBin::~ReaderNativeBin()
{
  delete [] line;
  memory->destroy(fieldindex);
}

/* ----------------------------------------------------------------------
   read and return time stamp from dump file
   if first read reaches end-of-file, return 1 so caller can open next file
   only called by proc 0
------------------------------------------------------------------------- */

int ReaderNativeBin::read_time(bigint &ntimestep)
{

  endian = 0x0001;
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
    fread(magic_string, sizeof(char), magic_string_len, fp);
    magic_string[magic_string_len] = '\0';

    // read endian flag
    fread(&endian, sizeof(int), 1, fp);

    // read revision number
    fread(&revision, sizeof(int), 1, fp);

    // read the real ntimestep
    fread(&ntimestep, sizeof(bigint), 1, fp);
  }

  return 0;
}

/* ----------------------------------------------------------------------
   skip snapshot from timestamp onward
   only called by proc 0
------------------------------------------------------------------------- */

void ReaderNativeBin::skip()
{
  int nchunk;
  fread(&nchunk, sizeof(int), 1, fp);
  if (feof(fp)) {
    error->one(FLERR,"Unexpected end of dump file");
  }
  for (int i = 0; i < nchunk; i++) {

    int n;
    fread(&n, sizeof(int), 1, fp);

    // extend buffer to fit chunk size

    if (n > maxbuf) {
      if (buf) delete[] buf;
      buf = new double[n];
      maxbuf = n;
    }
    // read chunk and write as size_one values per line
    fread(buf, sizeof(double), n, fp);
    if (feof(fp)) {
      error->one(FLERR,"Unexpected end of dump file");
    }
  }

  delete[] magic_string;
  delete[] unit_style;
  magic_string = nullptr;
  unit_style = nullptr;
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

  fread(&natoms, sizeof(bigint), 1, fp);

  boxinfo = 1;
  triclinic = 0;
  box[0][2] = box[1][2] = box[2][2] = 0.0;

  int boundary[3][2];
  fread(&triclinic, sizeof(int), 1, fp);
  fread(&boundary[0][0], 6 * sizeof(int), 1, fp);
  fread(&box[0][0], sizeof(double), 1, fp);
  fread(&box[0][1], sizeof(double), 1, fp);
  fread(&box[1][0], sizeof(double), 1, fp);
  fread(&box[1][1], sizeof(double), 1, fp);
  fread(&box[2][0], sizeof(double), 1, fp);
  fread(&box[2][1], sizeof(double), 1, fp);
  if (triclinic) {
    fread(&box[0][2], sizeof(double), 1, fp);
    fread(&box[1][2], sizeof(double), 1, fp);
    fread(&box[2][2], sizeof(double), 1, fp);
  }

  if (!fieldinfo) return natoms;

  // exatract column labels and match to requested fields
  fread(&size_one, sizeof(int), 1, fp);

  int len = 0;
  char *labelline;

  if (magic_string && revision > 0x0001) {
    // newer format includes units string, columns string
    // and time
    fread(&len, sizeof(int), 1, fp);
    labelline = new char[len + 1];

    if (len > 0) {
      // has units
      delete[] unit_style;
      unit_style = new char[len + 1];
      fread(unit_style, sizeof(char), len, fp);
      unit_style[len] = '\0';
    }

    char flag = 0;
    fread(&flag, sizeof(char), 1, fp);

    if (flag) {
      double time;
      fread(&time, sizeof(double), 1, fp);
    }

    fread(&len, sizeof(int), 1, fp);
    fread(labelline, sizeof(char), len, fp);
    labelline[len] = '\0';
  }

  std::map<std::string, int> labels;
  Tokenizer tokens(labelline);
  nwords = 0;

  while (tokens.has_next()) {
    labels[tokens.next()] = nwords++;
  }

  if (nwords == 0) {
    return 1;
  }

  match_field(nfield, xflag, yflag, zflag, fieldtype, fieldlabel, scaleflag, wrapflag, fieldflag, labels);

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
  fread(&nchunk, sizeof(int), 1, fp);
  for (int i = 0; i < nchunk; i++) {

    fread(&n, sizeof(int), 1, fp);

    // extend buffer to fit chunk size

    if (n > maxbuf) {
      if (buf) delete[] buf;
      buf = new double[n];
      maxbuf = n;
    }

    // read chunk and write as size_one values per line

    fread(buf, sizeof(double), n, fp);
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
  magic_string = nullptr;
  unit_style = nullptr;
}