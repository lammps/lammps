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

/* ----------------------------------------------------------------------
   Contributing author: Paul Coffman (IBM)
------------------------------------------------------------------------- */

#include "dump_custom_mpiio.h"

#include "domain.h"
#include "error.h"
#include "fix.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

#include "omp_compat.h"
#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;

#define DUMP_BUF_CHUNK_SIZE 16384
#define DUMP_BUF_INCREMENT_SIZE 4096

/* ---------------------------------------------------------------------- */

DumpCustomMPIIO::DumpCustomMPIIO(LAMMPS *lmp, int narg, char **arg) : DumpCustom(lmp, narg, arg)
{
  if (me == 0)
    error->warning(FLERR, "MPI-IO output is unmaintained and unreliable. Use with caution.");
}

/* ---------------------------------------------------------------------- */

DumpCustomMPIIO::~DumpCustomMPIIO()
{
  if (multifile == 0) MPI_File_close(&mpifh);
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::openfile()
{
  if (singlefile_opened) {    // single file already opened, so just return after resetting filesize
    mpifo = currentFileSize;
    MPI_File_set_size(mpifh, mpifo + headerSize + sumFileSize);
    currentFileSize = mpifo + headerSize + sumFileSize;
    return;
  }
  if (multifile == 0) singlefile_opened = 1;

  // if one file per timestep, replace '*' with current timestep

  filecurrent = filename;

  if (multifile) {
    filecurrent = utils::strdup(utils::star_subst(filecurrent, update->ntimestep, padflag));
    if (maxfiles > 0) {
      if (numfiles < maxfiles) {
        nameslist[numfiles] = utils::strdup(filecurrent);
        ++numfiles;
      } else {
        remove(nameslist[fileidx]);
        delete[] nameslist[fileidx];
        nameslist[fileidx] = utils::strdup(filecurrent);
        fileidx = (fileidx + 1) % maxfiles;
      }
    }
  }

  if (append_flag) {    // append open
    int err = MPI_File_open(world, filecurrent, MPI_MODE_CREATE | MPI_MODE_APPEND | MPI_MODE_WRONLY,
                            MPI_INFO_NULL, &mpifh);
    if (err != MPI_SUCCESS)
      error->one(FLERR, "Cannot open dump file {}: {}", filecurrent, utils::getsyserror());

    int myrank;
    MPI_Comm_rank(world, &myrank);
    if (myrank == 0) MPI_File_get_size(mpifh, &mpifo);
    MPI_Bcast(&mpifo, 1, MPI_LMP_BIGINT, 0, world);
    MPI_File_set_size(mpifh, mpifo + headerSize + sumFileSize);
    currentFileSize = mpifo + headerSize + sumFileSize;

  } else {    // replace open

    int err =
        MPI_File_open(world, filecurrent, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifh);
    if (err != MPI_SUCCESS)
      error->one(FLERR, "Cannot open dump file {}: {}", filecurrent, utils::getsyserror());

    mpifo = 0;

    MPI_File_set_size(mpifh, (MPI_Offset) (headerSize + sumFileSize));
    currentFileSize = (headerSize + sumFileSize);
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::write()
{
  if (domain->triclinic == 0) {
    boxxlo = domain->boxlo[0];
    boxxhi = domain->boxhi[0];
    boxylo = domain->boxlo[1];
    boxyhi = domain->boxhi[1];
    boxzlo = domain->boxlo[2];
    boxzhi = domain->boxhi[2];
  } else {
    boxxlo = domain->boxlo_bound[0];
    boxxhi = domain->boxhi_bound[0];
    boxylo = domain->boxlo_bound[1];
    boxyhi = domain->boxhi_bound[1];
    boxzlo = domain->boxlo_bound[2];
    boxzhi = domain->boxhi_bound[2];
    boxxy = domain->xy;
    boxxz = domain->xz;
    boxyz = domain->yz;
  }

  // nme = # of dump lines this proc contributes to dump

  nme = count();

  // ntotal = total # of dump lines in snapshot
  // nmax = max # of dump lines on any proc

  bigint bnme = nme;
  MPI_Allreduce(&bnme, &ntotal, 1, MPI_LMP_BIGINT, MPI_SUM, world);

  int nmax;
  MPI_Allreduce(&nme, &nmax, 1, MPI_INT, MPI_MAX, world);

  // write timestep header
  // for multiproc,
  //   nheader = # of lines in this file via Allreduce on clustercomm

  bigint nheader = ntotal;

  // insure filewriter proc can receive everyone's info
  // limit nmax*size_one to int since used as arg in MPI_Rsend() below
  // pack my data into buf
  // if sorting on IDs also request ID list from pack()
  // sort buf as needed

  if (nmax > maxbuf) {
    if ((bigint) nmax * size_one > MAXSMALLINT)
      error->all(FLERR, "Too much per-proc info for dump");
    maxbuf = nmax;
    memory->destroy(buf);
    memory->create(buf, (maxbuf * size_one), "dump:buf");
  }
  if (sort_flag && sortcol == 0 && nmax > maxids) {
    maxids = nmax;
    memory->destroy(ids);
    memory->create(ids, maxids, "dump:ids");
  }

  if (sort_flag && sortcol == 0)
    pack(ids);
  else
    pack(nullptr);
  if (sort_flag) sort();

  // determine how much data needs to be written for setting the file size and prepocess it prior to writing
  performEstimate = 1;
  write_header(nheader);
  write_data(nme, buf);
  MPI_Bcast(&sumFileSize, 1, MPI_LMP_BIGINT, (nprocs - 1), world);

  openfile();

  performEstimate = 0;
  write_header(nheader);    // mpifo now points to end of header info

  // now actually write the data
  performEstimate = 0;
  write_data(nme, buf);

  if (multifile) MPI_File_close(&mpifh);
  if (multifile) delete[] filecurrent;
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::init_style()
{
  // assemble ITEMS: column string from defaults and user values

  delete[] columns;
  std::string combined;
  int icol = 0;
  for (const auto &item : utils::split_words(columns_default)) {
    if (combined.size()) combined += " ";
    if (keyword_user[icol].size())
      combined += keyword_user[icol];
    else
      combined += item;
    ++icol;
  }
  columns = utils::strdup(combined);

  // format = copy of default or user-specified line format

  delete[] format;
  if (format_line_user)
    format = utils::strdup(format_line_user);
  else
    format = utils::strdup(format_default);

  // tokenize the format string and add space at end of each format element
  // if user-specified int/float format exists, use it instead
  // if user-specified column format exists, use it instead
  // lo priority = line, medium priority = int/float, hi priority = column

  auto words = utils::split_words(format);
  if ((int) words.size() < nfield)
    error->all(FLERR, "Dump_modify format line is too short: {}", format);

  int i = 0;
  for (const auto &word : words) {
    if (i >= nfield) break;
    delete[] vformat[i];

    if (format_column_user[i])
      vformat[i] = utils::strdup(std::string(format_column_user[i]) + " ");
    else if (vtype[i] == Dump::INT && format_int_user)
      vformat[i] = utils::strdup(std::string(format_int_user) + " ");
    else if (vtype[i] == Dump::DOUBLE && format_float_user)
      vformat[i] = utils::strdup(std::string(format_float_user) + " ");
    else if (vtype[i] == Dump::BIGINT && format_bigint_user)
      vformat[i] = utils::strdup(std::string(format_bigint_user) + " ");
    else
      vformat[i] = utils::strdup(word + " ");

    // remove trailing blank on last column's format
    if (i == nfield - 1) vformat[i][strlen(vformat[i]) - 1] = '\0';

    ++i;
  }

  // setup boundary string

  domain->boundary_string(boundstr);

  // setup function ptrs

  if (binary && domain->triclinic == 0)
    header_choice = &DumpCustomMPIIO::header_binary;
  else if (binary && domain->triclinic == 1)
    header_choice = &DumpCustomMPIIO::header_binary_triclinic;
  else if (!binary && domain->triclinic == 0)
    header_choice = &DumpCustomMPIIO::header_item;
  else if (!binary && domain->triclinic == 1)
    header_choice = &DumpCustomMPIIO::header_item_triclinic;

  if (binary)
    write_choice = &DumpCustomMPIIO::write_binary;
  else
    write_choice = &DumpCustomMPIIO::write_string;

  // find current ptr for each compute,fix,variable
  // check that fix frequency is acceptable

  for (i = 0; i < ncompute; i++) {
    compute[i] = modify->get_compute_by_id(id_compute[i]);
    if (!compute[i])
      error->all(FLERR, "Could not find dump custom/mpiio compute ID {}", id_compute[i]);
  }

  for (i = 0; i < nfix; i++) {
    fix[i] = modify->get_fix_by_id(id_fix[i]);
    if (!fix[i]) error->all(FLERR, "Could not find dump custom/mpiio fix ID {}", id_fix[i]);
    if (nevery % fix[i]->peratom_freq)
      error->all(FLERR, "dump custom/mpiio and fix not computed at compatible times");
  }

  for (i = 0; i < nvariable; i++) {
    int ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0)
      error->all(FLERR, "Could not find dump custom/mpiio variable name {}", id_variable[i]);
    variable[i] = ivariable;
  }

  // set index and check validity of region

  if (idregion && !domain->get_region_by_id(idregion))
    error->all(FLERR, "Region {} for dump custom/mpiio does not exist", idregion);
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::write_header(bigint ndump)
{
  (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::header_binary(bigint ndump)
{
  if (performEstimate) {

    headerBuffer = (char *) malloc((2 * sizeof(bigint)) + (9 * sizeof(int)) + (6 * sizeof(double)));

    headerSize = 0;
    memcpy(headerBuffer + headerSize, &update->ntimestep, sizeof(bigint));
    headerSize += sizeof(bigint);

    memcpy(headerBuffer + headerSize, &ndump, sizeof(bigint));
    headerSize += sizeof(bigint);

    memcpy(headerBuffer + headerSize, &domain->triclinic, sizeof(int));
    headerSize += sizeof(int);

    memcpy(headerBuffer + headerSize, &domain->boundary[0][0], 6 * sizeof(int));
    headerSize += 6 * sizeof(int);

    memcpy(headerBuffer + headerSize, &boxxlo, sizeof(double));
    headerSize += sizeof(double);

    memcpy(headerBuffer + headerSize, &boxxhi, sizeof(double));
    headerSize += sizeof(double);

    memcpy(headerBuffer + headerSize, &boxylo, sizeof(double));
    headerSize += sizeof(double);

    memcpy(headerBuffer + headerSize, &boxyhi, sizeof(double));
    headerSize += sizeof(double);

    memcpy(headerBuffer + headerSize, &boxzlo, sizeof(double));
    headerSize += sizeof(double);

    memcpy(headerBuffer + headerSize, &boxzhi, sizeof(double));
    headerSize += sizeof(double);

    memcpy(headerBuffer + headerSize, &size_one, sizeof(int));
    headerSize += sizeof(int);

    memcpy(headerBuffer + headerSize, &nprocs, sizeof(int));
    headerSize += sizeof(int);
  } else {    // write data
    if (me == 0)
      MPI_File_write_at(mpifh, mpifo, headerBuffer, headerSize, MPI_BYTE, MPI_STATUS_IGNORE);
    mpifo += headerSize;
    free(headerBuffer);
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::header_binary_triclinic(bigint ndump)
{
  if (performEstimate) {

    headerBuffer = (char *) malloc((2 * sizeof(bigint)) + (9 * sizeof(int)) + (9 * sizeof(double)));

    headerSize = 0;
    memcpy(headerBuffer + headerSize, &update->ntimestep, sizeof(bigint));
    headerSize += sizeof(bigint);

    memcpy(headerBuffer + headerSize, &ndump, sizeof(bigint));
    headerSize += sizeof(bigint);

    memcpy(headerBuffer + headerSize, &domain->triclinic, sizeof(int));
    headerSize += sizeof(int);

    memcpy(headerBuffer + headerSize, &domain->boundary[0][0], 6 * sizeof(int));
    headerSize += 6 * sizeof(int);

    memcpy(headerBuffer + headerSize, &boxxlo, sizeof(double));
    headerSize += sizeof(double);

    memcpy(headerBuffer + headerSize, &boxxhi, sizeof(double));
    headerSize += sizeof(double);

    memcpy(headerBuffer + headerSize, &boxylo, sizeof(double));
    headerSize += sizeof(double);

    memcpy(headerBuffer + headerSize, &boxyhi, sizeof(double));
    headerSize += sizeof(double);

    memcpy(headerBuffer + headerSize, &boxzlo, sizeof(double));
    headerSize += sizeof(double);

    memcpy(headerBuffer + headerSize, &boxzhi, sizeof(double));
    headerSize += sizeof(double);

    memcpy(headerBuffer + headerSize, &boxxy, sizeof(double));
    headerSize += sizeof(double);

    memcpy(headerBuffer + headerSize, &boxxz, sizeof(double));
    headerSize += sizeof(double);

    memcpy(headerBuffer + headerSize, &boxyz, sizeof(double));
    headerSize += sizeof(double);

    memcpy(headerBuffer + headerSize, &size_one, sizeof(int));
    headerSize += sizeof(int);

    memcpy(headerBuffer + headerSize, &nprocs, sizeof(int));
    headerSize += sizeof(int);

  } else {    // write data

    if (me == 0)
      MPI_File_write_at(mpifh, mpifo, headerBuffer, headerSize, MPI_BYTE, MPI_STATUS_IGNORE);
    mpifo += headerSize;
    free(headerBuffer);
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::header_item(bigint ndump)
{
  if (performEstimate) {

    auto itemtxt = fmt::format("ITEM: TIMESTEP\n{}\n", update->ntimestep);
    itemtxt += fmt::format("ITEM: NUMBER OF ATOMS\n{}\n", ndump);
    itemtxt += fmt::format("ITEM: BOX BOUNDS {}\n", boundstr);
    itemtxt += fmt::format("{} {}\n{} {}\n{} {}\n", boxxlo, boxxhi, boxylo, boxyhi, boxzlo, boxzhi);
    itemtxt += fmt::format("ITEM: ATOMS {}\n", columns);

    headerSize = itemtxt.size();
    headerBuffer = utils::strdup(itemtxt);

  } else {    // write data

    if (me == 0)
      MPI_File_write_at(mpifh, mpifo, headerBuffer, headerSize, MPI_CHAR, MPI_STATUS_IGNORE);
    mpifo += headerSize;
    delete[] headerBuffer;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::header_item_triclinic(bigint ndump)
{
  if (performEstimate) {

    auto itemtxt = fmt::format("ITEM: TIMESTEP\n{}\n", update->ntimestep);
    itemtxt += fmt::format("ITEM: NUMBER OF ATOMS\n{}\n", ndump);
    itemtxt += fmt::format("ITEM: BOX BOUNDS xy xz yz {}\n", boundstr);
    itemtxt += fmt::format("{} {} {}\n{} {} {}\n{} {} {}\n", boxxlo, boxxhi, boxxy, boxylo, boxyhi,
                           boxxz, boxzlo, boxzhi, boxyz);
    itemtxt += fmt::format("ITEM: ATOMS {}\n", columns);

    headerSize = itemtxt.size();
    headerBuffer = utils::strdup(itemtxt);

  } else {    // write data

    if (me == 0)
      MPI_File_write_at(mpifh, mpifo, headerBuffer, headerSize, MPI_CHAR, MPI_STATUS_IGNORE);
    mpifo += headerSize;
    delete[] headerBuffer;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n, mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::write_binary(int n, double *mybuf)
{
  n *= size_one;

  if (performEstimate) {

    bigint incPrefix = 0;
    bigint bigintNme = (bigint) nme;
    MPI_Scan(&bigintNme, &incPrefix, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    sumFileSize = (incPrefix * size_one * sizeof(double)) + (nprocs * sizeof(int));
    offsetFromHeader = ((incPrefix - bigintNme) * size_one * sizeof(double)) + (me * sizeof(int));
  } else {
    int byteBufSize = (n * sizeof(double)) + sizeof(int);

    char *bufWithSize;
    memory->create(bufWithSize, byteBufSize, "dump:bufWithSize");
    memcpy(bufWithSize, (char *) (&n), sizeof(int));
    memcpy(&((char *) bufWithSize)[sizeof(int)], mybuf, (n * sizeof(double)));
    MPI_File_write_at_all(mpifh, mpifo + offsetFromHeader, bufWithSize, byteBufSize, MPI_BYTE,
                          MPI_STATUS_IGNORE);
    memory->destroy(bufWithSize);

    if (flush_flag) MPI_File_sync(mpifh);
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::write_string(int n, double *mybuf)
{
  if (performEstimate) {

#if defined(_OPENMP)
    int nthreads = omp_get_max_threads();
    if ((nthreads > 1) && !(lmp->kokkos))
      nsme = convert_string_omp(n, mybuf);    // not (yet) compatible with Kokkos
    else
      nsme = convert_string(n, mybuf);
#else
    nsme = convert_string(n, mybuf);
#endif
    bigint incPrefix = 0;
    bigint bigintNsme = (bigint) nsme;
    MPI_Scan(&bigintNsme, &incPrefix, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    sumFileSize = (incPrefix * sizeof(char));
    offsetFromHeader = ((incPrefix - bigintNsme) * sizeof(char));
  } else {
    MPI_File_write_at_all(mpifh, mpifo + offsetFromHeader, sbuf, nsme, MPI_CHAR, MPI_STATUS_IGNORE);
    if (flush_flag) MPI_File_sync(mpifh);
  }
}

#if defined(_OPENMP)

/* ----------------------------------------------------------------------
   multithreaded version - convert mybuf of doubles to one big formatted string in sbuf
   return -1 if strlen exceeds an int, since used as arg in MPI calls in Dump
------------------------------------------------------------------------- */

int DumpCustomMPIIO::convert_string_omp(int n, double *mybuf)
{
  char **mpifh_buffer_line_per_thread;
  int mpifhStringCount;
  int *mpifhStringCountPerThread, *bufOffset, *bufRange, *bufLength;

  mpifhStringCount = 0;

  int nthreads = omp_get_max_threads();
  if (nthreads > n) {    // call serial version
    convert_string(n, mybuf);

  } else {
    memory->create(mpifhStringCountPerThread, nthreads, "dump:mpifhStringCountPerThread");
    mpifh_buffer_line_per_thread = (char **) malloc(nthreads * sizeof(char *));
    memory->create(bufOffset, nthreads, "dump:bufOffset");
    memory->create(bufRange, nthreads, "dump:bufRange");
    memory->create(bufLength, nthreads, "dump:bufLength");

    int i = 0;
    for (i = 0; i < (nthreads - 1); i++) {
      mpifhStringCountPerThread[i] = 0;
      bufOffset[i] = (int) (i * (int) (floor((double) n / (double) nthreads)) * size_one);
      bufRange[i] = (int) (floor((double) n / (double) nthreads));
      bufLength[i] = DUMP_BUF_CHUNK_SIZE;
      mpifh_buffer_line_per_thread[i] = (char *) malloc(DUMP_BUF_CHUNK_SIZE * sizeof(char));
      mpifh_buffer_line_per_thread[i][0] = '\0';
    }
    mpifhStringCountPerThread[i] = 0;
    bufOffset[i] = (int) (i * (int) (floor((double) n / (double) nthreads)) * size_one);
    bufRange[i] = n - (i * (int) (floor((double) n / (double) nthreads)));
    bufLength[i] = DUMP_BUF_CHUNK_SIZE;
    mpifh_buffer_line_per_thread[i] = (char *) malloc(DUMP_BUF_CHUNK_SIZE * sizeof(char));
    mpifh_buffer_line_per_thread[i][0] = '\0';

#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(bufOffset, bufRange, bufLength, \
                                                 mpifhStringCountPerThread,      \
                                                 mpifh_buffer_line_per_thread, mybuf)
    {
      int tid = omp_get_thread_num();
      int m = 0;

      for (int i = 0; i < bufRange[tid]; i++) {

        if ((bufLength[tid] - mpifhStringCountPerThread[tid]) < DUMP_BUF_INCREMENT_SIZE) {
          mpifh_buffer_line_per_thread[tid] = (char *) realloc(
              mpifh_buffer_line_per_thread[tid],
              (mpifhStringCountPerThread[tid] + DUMP_BUF_CHUNK_SIZE) * sizeof(char));
          bufLength[tid] = (mpifhStringCountPerThread[tid] + DUMP_BUF_CHUNK_SIZE) * sizeof(char);
        }
        for (int j = 0; j < size_one; j++) {

          if (vtype[j] == Dump::INT)
            mpifhStringCountPerThread[tid] +=
                sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),
                        vformat[j], static_cast<int>(mybuf[bufOffset[tid] + m]));
          else if (vtype[j] == Dump::DOUBLE)
            mpifhStringCountPerThread[tid] +=
                sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),
                        vformat[j], mybuf[bufOffset[tid] + m]);
          else if (vtype[j] == Dump::STRING)
            mpifhStringCountPerThread[tid] +=
                sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),
                        vformat[j], typenames[(int) mybuf[bufOffset[tid] + m]]);
          m++;
        }
        mpifhStringCountPerThread[tid] +=
            sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]), "\n");
      }
    }

#pragma omp barrier
    mpifhStringCount = 0;
    for (i = 0; i < nthreads; i++) { mpifhStringCount += mpifhStringCountPerThread[i]; }

    memory->destroy(bufOffset);
    memory->destroy(bufRange);
    memory->destroy(bufLength);

    if (mpifhStringCount > 0) {
      if (mpifhStringCount > maxsbuf) {
        if (mpifhStringCount > MAXSMALLINT) return -1;
        maxsbuf = mpifhStringCount + 1;
        memory->grow(sbuf, maxsbuf, "dump:sbuf");
      }
      sbuf[0] = '\0';
    }

    for (int i = 0; i < nthreads; i++) {
      strcat(sbuf, mpifh_buffer_line_per_thread[i]);
      free(mpifh_buffer_line_per_thread[i]);
    }

    memory->destroy(mpifhStringCountPerThread);
    free(mpifh_buffer_line_per_thread);
  }

  return mpifhStringCount;
}
#endif
