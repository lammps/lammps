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

#include "reader_native.h"

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

ReaderNative::ReaderNative(LAMMPS *lmp) : Reader(lmp)
{
  line = new char[MAXLINE];
  fieldindex = nullptr;
  maxbuf = 0;
  databuf = nullptr;
  magic_string = nullptr;
  unit_style = nullptr;
}

/* ---------------------------------------------------------------------- */

ReaderNative::~ReaderNative()
{
  delete[] line;
  delete[] magic_string;
  delete[] unit_style;
  memory->destroy(fieldindex);
  memory->destroy(databuf);
}

/* ----------------------------------------------------------------------
   read and return time stamp from dump file
   if first read reaches end-of-file, return 1 so caller can open next file
   only called by proc 0
------------------------------------------------------------------------- */

int ReaderNative::read_time(bigint &ntimestep)
{
  if (binary) {
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
  } else {
    char *eof = fgets(line,MAXLINE,fp);
    if (eof == nullptr) return 1;

    // skip over unit and time information, if present.

    if (utils::strmatch(line,"^\\s*ITEM: UNITS\\s*$"))
      read_lines(2);

    if (utils::strmatch(line,"^\\s*ITEM: TIME\\s*$"))
      read_lines(2);

    if (!utils::strmatch(line,"^\\s*ITEM: TIMESTEP\\s*$"))
      error->one(FLERR,"Dump file is incorrectly formatted");

    read_lines(1);
    int rv = sscanf(line,BIGINT_FORMAT,&ntimestep);
    if (rv != 1)
      error->one(FLERR,"Dump file is incorrectly formatted");
  }
  return 0;
}

/* ----------------------------------------------------------------------
   skip snapshot from timestamp onward
   only called by proc 0
------------------------------------------------------------------------- */

void ReaderNative::skip()
{
  if (binary) {
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
      skip_buf(n*sizeof(double));
    }

    delete[] magic_string;
    delete[] unit_style;
    magic_string = nullptr;
    unit_style = nullptr;
  } else {
    read_lines(2);
    bigint natoms;
    int rv = sscanf(line,BIGINT_FORMAT,&natoms);
    if (rv != 1)
      error->one(FLERR,"Dump file is incorrectly formatted");

    read_lines(5);

    // invoke read_lines() in chunks no larger than MAXSMALLINT

    int nchunk;
    bigint nremain = natoms;
    while (nremain) {
      nchunk = MIN(nremain,MAXSMALLINT);
      read_lines(nchunk);
      nremain -= nchunk;
    }
  }
}

void ReaderNative::skip_reading_magic_str()
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

bigint ReaderNative::read_header(double box[3][3], int &boxinfo, int &triclinic,
                                 int fieldinfo, int nfield,
                                 int *fieldtype, char **fieldlabel,
                                 int scaleflag, int wrapflag, int &fieldflag,
                                 int &xflag, int &yflag, int &zflag)
{
  bigint natoms = 0;
  int len = 0;
  char *labelline;

  if (binary) {
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

    // extract column labels and match to requested fields
    read_buf(&size_one, sizeof(int), 1);

    if (!fieldinfo) {
      skip_reading_magic_str();
      return natoms;
    }

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
  } else {
    int rv;

    read_lines(2);
    rv = sscanf(line,BIGINT_FORMAT,&natoms);
    if (rv != 1)
      error->one(FLERR,"Dump file is incorrectly formatted");

    boxinfo = 1;
    triclinic = 0;
    box[0][2] = box[1][2] = box[2][2] = 0.0;
    read_lines(1);
    if (line[strlen("ITEM: BOX BOUNDS ")] == 'x') triclinic = 1;

    read_lines(1);
    if (!triclinic) rv = 2 - sscanf(line,"%lg %lg",&box[0][0],&box[0][1]);
    else rv = 3 - sscanf(line,"%lg %lg %lg",&box[0][0],&box[0][1],&box[0][2]);
    if (rv != 0) error->one(FLERR,"Dump file is incorrectly formatted");

    read_lines(1);
    if (!triclinic) rv = 2 - sscanf(line,"%lg %lg",&box[1][0],&box[1][1]);
    else rv = 3 - sscanf(line,"%lg %lg %lg",&box[1][0],&box[1][1],&box[1][2]);
    if (rv != 0) error->one(FLERR,"Dump file is incorrectly formatted");

    read_lines(1);
    if (!triclinic) rv = 2 - sscanf(line,"%lg %lg",&box[2][0],&box[2][1]);
    else rv = 3 - sscanf(line,"%lg %lg %lg",&box[2][0],&box[2][1],&box[2][2]);
    if (rv != 0) error->one(FLERR,"Dump file is incorrectly formatted");

    read_lines(1);

    // if no field info requested, just return

    if (!fieldinfo) return natoms;

    // extract column labels and match to requested fields

    labelline = &line[strlen("ITEM: ATOMS ")];
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

  // match each field with a column of per-atom data
  // if fieldlabel set, match with explicit column
  // else infer one or more column matches from fieldtype
  // xyz flag set by scaleflag + wrapflag (if fieldlabel set) or column label

  memory->create(fieldindex,nfield,"read_dump:fieldindex");

  int s_index,u_index,su_index;
  xflag = UNSET;
  yflag = UNSET;
  zflag = UNSET;

  for (int i = 0; i < nfield; i++) {
    if (fieldlabel[i]) {
      fieldindex[i] = find_label(fieldlabel[i], labels);
      if (fieldtype[i] == X) xflag = 2*scaleflag + wrapflag + 1;
      else if (fieldtype[i] == Y) yflag = 2*scaleflag + wrapflag + 1;
      else if (fieldtype[i] == Z) zflag = 2*scaleflag + wrapflag + 1;
    }

    else if (fieldtype[i] == ID)
      fieldindex[i] = find_label("id", labels);
    else if (fieldtype[i] == TYPE)
      fieldindex[i] = find_label("type", labels);

    else if (fieldtype[i] == X) {
      fieldindex[i] = find_label("x", labels);
      xflag = NOSCALE_WRAP;
      if (fieldindex[i] < 0) {
        fieldindex[i] = nwords;
        s_index = find_label("xs", labels);
        u_index = find_label("xu", labels);
        su_index = find_label("xsu", labels);
        if (s_index >= 0 && s_index < fieldindex[i]) {
          fieldindex[i] = s_index;
          xflag = SCALE_WRAP;
        }
        if (u_index >= 0 && u_index < fieldindex[i]) {
          fieldindex[i] = u_index;
          xflag = NOSCALE_NOWRAP;
        }
        if (su_index >= 0 && su_index < fieldindex[i]) {
          fieldindex[i] = su_index;
          xflag = SCALE_NOWRAP;
        }
      }
      if (fieldindex[i] == nwords) fieldindex[i] = -1;

    } else if (fieldtype[i] == Y) {
      fieldindex[i] = find_label("y", labels);
      yflag = NOSCALE_WRAP;
      if (fieldindex[i] < 0) {
        fieldindex[i] = nwords;
        s_index = find_label("ys", labels);
        u_index = find_label("yu", labels);
        su_index = find_label("ysu", labels);
        if (s_index >= 0 && s_index < fieldindex[i]) {
          fieldindex[i] = s_index;
          yflag = SCALE_WRAP;
        }
        if (u_index >= 0 && u_index < fieldindex[i]) {
          fieldindex[i] = u_index;
          yflag = NOSCALE_NOWRAP;
        }
        if (su_index >= 0 && su_index < fieldindex[i]) {
          fieldindex[i] = su_index;
          yflag = SCALE_NOWRAP;
        }
      }
      if (fieldindex[i] == nwords) fieldindex[i] = -1;

    } else if (fieldtype[i] == Z) {
      fieldindex[i] = find_label("z", labels);
      zflag = NOSCALE_WRAP;
      if (fieldindex[i] < 0) {
        fieldindex[i] = nwords;
        s_index = find_label("zs", labels);
        u_index = find_label("zu", labels);
        su_index = find_label("zsu", labels);
        if (s_index >= 0 && s_index < fieldindex[i]) {
          fieldindex[i] = s_index;
          zflag = SCALE_WRAP;
        }
        if (u_index >= 0 && u_index < fieldindex[i]) {
          fieldindex[i] = u_index;
          zflag = NOSCALE_NOWRAP;
        }
        if (su_index >= 0 && su_index < fieldindex[i]) {
          fieldindex[i] = su_index;
          zflag = SCALE_NOWRAP;
        }
      }
      if (fieldindex[i] == nwords) fieldindex[i] = -1;

    } else if (fieldtype[i] == VX)
      fieldindex[i] = find_label("vx", labels);
    else if (fieldtype[i] == VY)
      fieldindex[i] = find_label("vy", labels);
    else if (fieldtype[i] == VZ)
      fieldindex[i] = find_label("vz", labels);

    else if (fieldtype[i] == FX)
      fieldindex[i] = find_label("fx", labels);
    else if (fieldtype[i] == FY)
      fieldindex[i] = find_label("fy", labels);
    else if (fieldtype[i] == FZ)
      fieldindex[i] = find_label("fz", labels);

    else if (fieldtype[i] == Q)
      fieldindex[i] = find_label("q", labels);

    else if (fieldtype[i] == IX)
      fieldindex[i] = find_label("ix", labels);
    else if (fieldtype[i] == IY)
      fieldindex[i] = find_label("iy", labels);
    else if (fieldtype[i] == IZ)
      fieldindex[i] = find_label("iz", labels);
  }

  // set fieldflag = -1 if any unfound fields

  fieldflag = 0;
  for (int i = 0; i < nfield; i++)
    if (fieldindex[i] < 0) fieldflag = -1;

  return natoms;
}

/* ----------------------------------------------------------------------
   read N atom lines from dump file
   stores appropriate values in fields array
   return 0 if success, 1 if error
   only called by proc 0
------------------------------------------------------------------------- */

void ReaderNative::read_atoms(int n, int nfield, double **fields)
{
  if (binary) {
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
          double *words = &databuf[m];
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

  } else {
    int i,m;
    char *eof;

    for (i = 0; i < n; i++) {
      eof = fgets(line,MAXLINE,fp);
      if (eof == nullptr) error->one(FLERR,"Unexpected end of dump file");

      // tokenize the line
      std::vector<std::string> words = Tokenizer(line).as_vector();

      if ((int)words.size() < nwords) error->one(FLERR,"Insufficient columns in dump file");

      // convert selected fields to floats

      for (m = 0; m < nfield; m++)
        fields[i][m] = atof(words[fieldindex[m]].c_str());
    }
  }
}

/* ----------------------------------------------------------------------
   match label to any of N labels
   return index of match or -1 if no match
------------------------------------------------------------------------- */

int ReaderNative::find_label(const std::string &label, const std::map<std::string, int> & labels)
{
  auto it = labels.find(label);
  if (it != labels.end()) {
      return it->second;
  }
  return -1;
}

/* ----------------------------------------------------------------------
   read N lines from dump file
   only last one is saved in line
   only called by proc 0
------------------------------------------------------------------------- */

void ReaderNative::read_lines(int n)
{
  char *eof = nullptr;
  if (n <= 0) return;
  for (int i = 0; i < n; i++) eof = fgets(line,MAXLINE,fp);
  if (eof == nullptr) error->one(FLERR,"Unexpected end of dump file");
}

void ReaderNative::read_buf(void * ptr, size_t size, size_t count)
{
  fread(ptr, size, count, fp);

  // detect end-of-file
  if (feof(fp)) error->one(FLERR,"Unexpected end of dump file");
}

void ReaderNative::read_double_chunk(size_t count)
{
  // extend buffer to fit chunk size
  if (count > maxbuf) {
    memory->grow(databuf,count,"reader:databuf");
    maxbuf = count;
  }
  read_buf(databuf, sizeof(double), count);
}

void ReaderNative::skip_buf(size_t size)
{
  bigint pos = platform::ftell(fp);
  pos += size;
  platform::fseek(fp,pos);
}
