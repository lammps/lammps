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
   Contributing author: Paul Coffman (IBM)
------------------------------------------------------------------------- */

#include "dump_xyz_mpiio.h"

#include "compute.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "update.h"

#include <cmath>
#include <cstring>

#include "omp_compat.h"
#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;

#define MAX_TEXT_HEADER_SIZE 4096
#define DUMP_BUF_CHUNK_SIZE 16384
#define DUMP_BUF_INCREMENT_SIZE 4096

/* ---------------------------------------------------------------------- */

DumpXYZMPIIO::DumpXYZMPIIO(LAMMPS *lmp, int narg, char **arg) :
  DumpXYZ(lmp, narg, arg) {
  if (me == 0)
    error->warning(FLERR,"MPI-IO output is unmaintained and unreliable. Use with caution.");
}

/* ---------------------------------------------------------------------- */

DumpXYZMPIIO::~DumpXYZMPIIO()
{
  if (multifile == 0) MPI_File_close(&mpifh);
}

/* ---------------------------------------------------------------------- */

void DumpXYZMPIIO::openfile()
{
  if (singlefile_opened) { // single file already opened, so just return after resetting filesize
    mpifo = currentFileSize;
    MPI_File_set_size(mpifh,mpifo+headerSize+sumFileSize);
    currentFileSize = mpifo+headerSize+sumFileSize;
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

  if (append_flag) { // append open
    int err = MPI_File_open( world, filecurrent, MPI_MODE_CREATE | MPI_MODE_APPEND |
                             MPI_MODE_WRONLY  , MPI_INFO_NULL, &mpifh);
    if (err != MPI_SUCCESS) error->one(FLERR, "Cannot open dump file {}", filecurrent);

    int myrank;
    MPI_Comm_rank(world,&myrank);
    if (myrank == 0) MPI_File_get_size(mpifh,&mpifo);
    MPI_Bcast(&mpifo, 1, MPI_LMP_BIGINT, 0, world);
    MPI_File_set_size(mpifh,mpifo+headerSize+sumFileSize);
    currentFileSize = mpifo+headerSize+sumFileSize;
  } else { // replace open

    int err = MPI_File_open( world, filecurrent, MPI_MODE_CREATE |  MPI_MODE_WRONLY,
                             MPI_INFO_NULL, &mpifh);
    if (err != MPI_SUCCESS) error->one(FLERR, "Cannot open dump file {}",filecurrent);
    mpifo = 0;

    MPI_File_set_size(mpifh,(MPI_Offset) (headerSize+sumFileSize));
    currentFileSize = (headerSize+sumFileSize);
  }
}

/* ---------------------------------------------------------------------- */

void DumpXYZMPIIO::write()
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
  MPI_Allreduce(&bnme,&ntotal,1,MPI_LMP_BIGINT,MPI_SUM,world);

  int nmax;
  MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);

  // write timestep header
  // for multiproc,
  //   nheader = # of lines in this file via Allreduce on clustercomm

  bigint nheader = ntotal;

  // ensure filewriter proc can receive everyone's info
  // limit nmax*size_one to int since used as arg in MPI_Rsend() below
  // pack my data into buf
  // if sorting on IDs also request ID list from pack()
  // sort buf as needed

  if (nmax > maxbuf) {
    if ((bigint) nmax * size_one > MAXSMALLINT)
      error->all(FLERR,"Too much per-proc info for dump");
    maxbuf = nmax;
    memory->destroy(buf);
    memory->create(buf,(maxbuf*size_one),"dump:buf");
  }
  if (sort_flag && sortcol == 0 && nmax > maxids) {
    maxids = nmax;
    memory->destroy(ids);
    memory->create(ids,maxids,"dump:ids");
  }

  if (sort_flag && sortcol == 0) pack(ids);
  else pack(nullptr);
  if (sort_flag) sort();

  // determine how much data needs to be written for setting the file size and prepocess it prior to writing
  performEstimate = 1;
  write_header(nheader);
  write_data(nme,buf);
  MPI_Bcast(&sumFileSize, 1, MPI_LMP_BIGINT, (nprocs-1), world);

  openfile();

  performEstimate = 0;
  write_header(nheader); // actually write the header - mpifo now points to end of header info

  // now actually write the data
  performEstimate = 0;
  write_data(nme,buf);

  if (multifile) MPI_File_close(&mpifh);
  if (multifile) delete[] filecurrent;
}

/* ---------------------------------------------------------------------- */

void DumpXYZMPIIO::init_style()
{
  // format = copy of default or user-specified line format

  delete[] format;
  char *str;
  if (format_line_user) str = format_line_user;
  else str = format_default;

  int n = strlen(str) + 2;
  format = new char[n];
  strcpy(format,str);
  strcat(format,"\n");

  // initialize typenames array to be backward compatible by default
  // a 32-bit int can be maximally 10 digits plus sign

  if (typenames == nullptr) {
    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++)
      typenames[itype] = utils::strdup(std::to_string(itype));
  }

  // setup function ptr
  write_choice = &DumpXYZMPIIO::write_string;
}

/* ---------------------------------------------------------------------- */

void DumpXYZMPIIO::write_header(bigint n)
{
  auto header = fmt::format("{}\n Atoms. Timestep: {}", n, update->ntimestep);
  if (time_flag) header += fmt::format(" Time: {:.6f}", compute_time());
  header += "\n";

  if (performEstimate) {
    headerSize = header.size();
  } else { // write data
    if (me == 0)
      MPI_File_write_at(mpifh,mpifo,(void *)header.c_str(),header.size(),MPI_CHAR,MPI_STATUS_IGNORE);
    mpifo += header.size();
  }
}

/* ---------------------------------------------------------------------- */

void DumpXYZMPIIO::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpXYZMPIIO::write_string(int n, double *mybuf)
{
  if (performEstimate) {

#if defined(_OPENMP)
    int nthreads = omp_get_max_threads();
    if (nthreads > 1)
      nsme = convert_string_omp(n,mybuf);
    else
      nsme = convert_string(n,mybuf);
#else

    nsme = convert_string(n,mybuf);
#endif
    bigint incPrefix = 0;
    bigint bigintNsme = (bigint) nsme;
    MPI_Scan(&bigintNsme,&incPrefix,1,MPI_LMP_BIGINT,MPI_SUM,world);
    sumFileSize = (incPrefix*sizeof(char));
    offsetFromHeader = ((incPrefix-bigintNsme)*sizeof(char));
  }
  else { // write data
    MPI_File_write_at_all(mpifh,mpifo+offsetFromHeader,sbuf,nsme,MPI_CHAR,MPI_STATUS_IGNORE);
    if (flush_flag)
      MPI_File_sync(mpifh);
  }
}

#if defined(_OPENMP)

/* ----------------------------------------------------------------------
   multithreaded version - convert mybuf of doubles to one big formatted string in sbuf
   return -1 if strlen exceeds an int, since used as arg in MPI calls in Dump
------------------------------------------------------------------------- */

int DumpXYZMPIIO::convert_string_omp(int n, double *mybuf)
{
  char **mpifh_buffer_line_per_thread;
  int mpifhStringCount;
  int *mpifhStringCountPerThread, *bufOffset, *bufRange, *bufLength;

  mpifhStringCount = 0;

  int nthreads = omp_get_max_threads();
  if (nthreads > n) { // call serial version
    convert_string(n,mybuf);

  }
  else {
    memory->create(mpifhStringCountPerThread,nthreads,"dump:mpifhStringCountPerThread");
    mpifh_buffer_line_per_thread = (char **) malloc(nthreads*sizeof(char*));
    memory->create(bufOffset,nthreads,"dump:bufOffset");
    memory->create(bufRange,nthreads,"dump:bufRange");
    memory->create(bufLength,nthreads,"dump:bufLength");

    int i=0;
    for (i=0;i<(nthreads-1);i++) {
      mpifhStringCountPerThread[i] = 0;
      bufOffset[i] = (int) (i*(int)(floor((double)n/(double)nthreads))*size_one);
      bufRange[i] = (int)(floor((double)n/(double)nthreads));
      bufLength[i] = DUMP_BUF_CHUNK_SIZE;
      mpifh_buffer_line_per_thread[i] = (char *) malloc(DUMP_BUF_CHUNK_SIZE * sizeof(char));
      mpifh_buffer_line_per_thread[i][0] = '\0';
    }
    mpifhStringCountPerThread[i] = 0;
    bufOffset[i] = (int) (i*(int)(floor((double)n/(double)nthreads))*size_one);
    bufRange[i] =  n-(i*(int)(floor((double)n/(double)nthreads)));
    bufLength[i] = DUMP_BUF_CHUNK_SIZE;
    mpifh_buffer_line_per_thread[i] = (char *) malloc(DUMP_BUF_CHUNK_SIZE * sizeof(char));
    mpifh_buffer_line_per_thread[i][0] = '\0';

#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(bufOffset, bufRange, bufLength, mpifhStringCountPerThread, mpifh_buffer_line_per_thread, mybuf)
    {
      int tid = omp_get_thread_num();
      int m=0;

      for (int i = 0; i < bufRange[tid]; i++) {

        if ((bufLength[tid] - mpifhStringCountPerThread[tid]) < DUMP_BUF_INCREMENT_SIZE) {
          mpifh_buffer_line_per_thread[tid] = (char *) realloc(mpifh_buffer_line_per_thread[tid],(mpifhStringCountPerThread[tid]+DUMP_BUF_CHUNK_SIZE) * sizeof(char));
          bufLength[tid] = (mpifhStringCountPerThread[tid]+DUMP_BUF_CHUNK_SIZE) * sizeof(char);
        }

        mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),format,typenames[static_cast<int> (mybuf[bufOffset[tid]+m+1])],mybuf[bufOffset[tid]+m+2],mybuf[bufOffset[tid]+m+3],mybuf[bufOffset[tid]+m+4]);
          m += size_one;

      }
    }

#pragma omp barrier
    mpifhStringCount = 0;
    for (i=0;i<nthreads;i++) {
      mpifhStringCount += mpifhStringCountPerThread[i];
    }

    memory->destroy(bufOffset);
    memory->destroy(bufRange);
    memory->destroy(bufLength);

    if (mpifhStringCount > 0) {
      if (mpifhStringCount > maxsbuf) {
        if (mpifhStringCount > MAXSMALLINT) return -1;
        maxsbuf = mpifhStringCount+1;
        memory->grow(sbuf,maxsbuf,"dump:sbuf");
      }
      sbuf[0] = '\0';
    }


    for (int i=0;i<nthreads;i++) {
      strcat(sbuf,mpifh_buffer_line_per_thread[i]);
      free(mpifh_buffer_line_per_thread[i]);
    }

    memory->destroy(mpifhStringCountPerThread);
    free(mpifh_buffer_line_per_thread);

  }

  return mpifhStringCount;
}
#endif
