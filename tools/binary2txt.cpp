/* -----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------ */

// Convert one or more LAMMPS binary dump files to ASCII text files
//
// this serial code must be compiled on a platform that can read the binary
//   dump files since binary formats are not compatible across all platforms
//
// Syntax: binary2txt file1 file2 ...
// Creates:           file1.txt file2.txt ...

#include <cstdio>
#include <cstring>

// these must match settings in src/lmptype.h which builds LAMMPS with
//   -DLAMMPS_SMALLBIG (the default), -DLAMMPS_BIGBIG, or -DLAMMPS_SMALLSMALL
// you can edit the tools/Makefile to enforce the same setting
//   for the build of binary2txt, e.g.
//   g++ -g -DLAMMPS_BIGBIG binarytxt.o -o binary2txt
//   again -DLAMMPS_SMALLBIG is the default

#define __STDC_FORMAT_MACROS
#include <cinttypes>

#ifndef PRId64
#define PRId64 "ld"
#endif

#if !defined(LAMMPS_SMALLSMALL) && !defined(LAMMPS_BIGBIG) && !defined(LAMMPS_SMALLBIG)
#define LAMMPS_SMALLBIG
#endif

#if defined(LAMMPS_SMALLBIG)
typedef int tagint;
typedef int64_t bigint;
#define BIGINT_FORMAT "%" PRId64
#elif defined(LAMMPS_SMALLSMALL)
typedef int tagint;
typedef int bigint;
#define BIGINT_FORMAT "%d"
#else /* LAMMPS_BIGBIG */
typedef int64_t tagint;
typedef int64_t bigint;
#define BIGINT_FORMAT "%" PRId64
#endif

int main(int narg, char **arg)
{
  int i, j, k, m, n;
  bigint ntimestep, natoms;
  int size_one, nchunk, triclinic;
  double xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz;
  int boundary[3][2];
  char boundstr[9];

  int maxbuf = 0;
  double *buf = nullptr;

  if (narg == 1) {
    printf("Syntax: binary2txt file1 file2 ...\n");
    return 1;
  }

  // loop over files

  for (int iarg = 1; iarg < narg; iarg++) {
    printf("%s:", arg[iarg]);
    fflush(stdout);
    FILE *fp = fopen(arg[iarg], "rb");
    if (!fp) {
      printf("ERROR: Could not open %s\n", arg[iarg]);
      return 1;
    }

    n = strlen(arg[iarg]) + 1 + 4;
    auto filetxt = new char[n];
    strcpy(filetxt, arg[iarg]);
    strcat(filetxt, ".txt");
    FILE *fptxt = fopen(filetxt, "w");
    delete[] filetxt;

    // detect newer format
    char *magic_string = nullptr;
    char *columns = nullptr;
    char *unit_style = nullptr;

    // loop over snapshots in file

    while (true) {
      int endian = 0x0001;
      int revision = 0x0001;

      fread(&ntimestep, sizeof(bigint), 1, fp);

      // detect end-of-file

      if (feof(fp)) {
        fclose(fp);
        fclose(fptxt);
        break;
      }

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

      fread(&natoms, sizeof(bigint), 1, fp);
      fread(&triclinic, sizeof(int), 1, fp);
      fread(&boundary[0][0], 6 * sizeof(int), 1, fp);
      fread(&xlo, sizeof(double), 1, fp);
      fread(&xhi, sizeof(double), 1, fp);
      fread(&ylo, sizeof(double), 1, fp);
      fread(&yhi, sizeof(double), 1, fp);
      fread(&zlo, sizeof(double), 1, fp);
      fread(&zhi, sizeof(double), 1, fp);
      if (triclinic) {
        fread(&xy, sizeof(double), 1, fp);
        fread(&xz, sizeof(double), 1, fp);
        fread(&yz, sizeof(double), 1, fp);
      }
      fread(&size_one, sizeof(int), 1, fp);

      if (magic_string && revision > 0x0001) {
        // newer format includes units string, columns string
        // and time
        int len = 0;
        fread(&len, sizeof(int), 1, fp);

        if (len > 0) {
          // has units
          delete[] unit_style;
          unit_style = new char[len + 1];
          fread(unit_style, sizeof(char), len, fp);
          unit_style[len] = '\0';
          fprintf(fptxt, "ITEM: UNITS\n");
          fprintf(fptxt, "%s\n", unit_style);
        }

        char flag = 0;
        fread(&flag, sizeof(char), 1, fp);

        if (flag) {
          double time;
          fread(&time, sizeof(double), 1, fp);
          fprintf(fptxt, "ITEM: TIME\n%.16g\n", time);
        }

        fread(&len, sizeof(int), 1, fp);
        delete[] columns;
        columns = new char[len + 1];
        fread(columns, sizeof(char), len, fp);
        columns[len] = '\0';
      }

      fread(&nchunk, sizeof(int), 1, fp);

      fprintf(fptxt, "ITEM: TIMESTEP\n");
      fprintf(fptxt, BIGINT_FORMAT "\n", ntimestep);
      fprintf(fptxt, "ITEM: NUMBER OF ATOMS\n");
      fprintf(fptxt, BIGINT_FORMAT "\n", natoms);

      m = 0;
      for (int idim = 0; idim < 3; idim++) {
        for (int iside = 0; iside < 2; iside++) {
          if (boundary[idim][iside] == 0)
            boundstr[m++] = 'p';
          else if (boundary[idim][iside] == 1)
            boundstr[m++] = 'f';
          else if (boundary[idim][iside] == 2)
            boundstr[m++] = 's';
          else if (boundary[idim][iside] == 3)
            boundstr[m++] = 'm';
        }
        boundstr[m++] = ' ';
      }
      boundstr[8] = '\0';

      if (!triclinic) {
        fprintf(fptxt, "ITEM: BOX BOUNDS %s\n", boundstr);
        fprintf(fptxt, "%-1.16e %-1.16e\n", xlo, xhi);
        fprintf(fptxt, "%-1.16e %-1.16e\n", ylo, yhi);
        fprintf(fptxt, "%-1.16e %-1.16e\n", zlo, zhi);
      } else {
        fprintf(fptxt, "ITEM: BOX BOUNDS xy xz yz %s\n", boundstr);
        fprintf(fptxt, "%-1.16e %-1.16e %-1.16e\n", xlo, xhi, xy);
        fprintf(fptxt, "%-1.16e %-1.16e %-1.16e\n", ylo, yhi, xz);
        fprintf(fptxt, "%-1.16e %-1.16e %-1.16e\n", zlo, zhi, yz);
      }

      if (columns)
        fprintf(fptxt, "ITEM: ATOMS %s\n", columns);
      else
        fprintf(fptxt, "ITEM: ATOMS\n");

      // loop over processor chunks in file

      for (i = 0; i < nchunk; i++) {
        fread(&n, sizeof(int), 1, fp);

        // extend buffer to fit chunk size

        if (n > maxbuf) {
          delete[] buf;
          buf = new double[n];
          maxbuf = n;
        }

        // read chunk and write as size_one values per line

        fread(buf, sizeof(double), n, fp);
        n /= size_one;
        m = 0;
        for (j = 0; j < n; j++) {
          for (k = 0; k < size_one; k++) {
            if (k + 1 < size_one) {
              fprintf(fptxt, "%g ", buf[m++]);
            } else {
              fprintf(fptxt, "%g", buf[m++]);
            }
          }
          fprintf(fptxt, "\n");
        }
      }

      printf(" " BIGINT_FORMAT, ntimestep);
      fflush(stdout);
    }
    printf("\n");
    delete[] columns;
    delete[] magic_string;
    delete[] unit_style;
    columns = nullptr;
    magic_string = nullptr;
    unit_style = nullptr;
  }

  delete[] buf;
  return 0;
}
