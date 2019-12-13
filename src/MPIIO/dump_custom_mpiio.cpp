/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
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
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "domain.h"
#include "input.h"
#include "variable.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;

#define MAX_TEXT_HEADER_SIZE 4096
#define DUMP_BUF_CHUNK_SIZE 16384
#define DUMP_BUF_INCREMENT_SIZE 4096

enum{ID,MOL,TYPE,ELEMENT,MASS,
     X,Y,Z,XS,YS,ZS,XSTRI,YSTRI,ZSTRI,XU,YU,ZU,XUTRI,YUTRI,ZUTRI,
     XSU,YSU,ZSU,XSUTRI,YSUTRI,ZSUTRI,
     IX,IY,IZ,
     VX,VY,VZ,FX,FY,FZ,
     Q,MUX,MUY,MUZ,MU,RADIUS,DIAMETER,
     OMEGAX,OMEGAY,OMEGAZ,ANGMOMX,ANGMOMY,ANGMOMZ,
     TQX,TQY,TQZ,SPIN,ERADIUS,ERVEL,ERFORCE,
     COMPUTE,FIX,VARIABLE};
enum{LT,LE,GT,GE,EQ,NEQ};

/* ---------------------------------------------------------------------- */

DumpCustomMPIIO::DumpCustomMPIIO(LAMMPS *lmp, int narg, char **arg) :
  DumpCustom(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

DumpCustomMPIIO::~DumpCustomMPIIO()
{
  if (multifile == 0) MPI_File_close(&mpifh);
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::openfile()
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
    char *filestar = filecurrent;
    filecurrent = new char[strlen(filestar) + 16];
    char *ptr = strchr(filestar,'*');
    *ptr = '\0';
    if (padflag == 0)
      sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
          filestar,update->ntimestep,ptr+1);
    else {
      char bif[8],pad[16];
      strcpy(bif,BIGINT_FORMAT);
      sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
      sprintf(filecurrent,pad,filestar,update->ntimestep,ptr+1);
    }
    *ptr = '*';
    if (maxfiles > 0) {
      if (numfiles < maxfiles) {
        nameslist[numfiles] = new char[strlen(filecurrent)+1];
        strcpy(nameslist[numfiles],filecurrent);
        ++numfiles;
      } else {
        remove(nameslist[fileidx]);
        delete[] nameslist[fileidx];
        nameslist[fileidx] = new char[strlen(filecurrent)+1];
        strcpy(nameslist[fileidx],filecurrent);
        fileidx = (fileidx + 1) % maxfiles;
      }
    }
  }

  if (append_flag) { // append open
    int err = MPI_File_open( world, filecurrent, MPI_MODE_CREATE | MPI_MODE_APPEND | MPI_MODE_WRONLY  , MPI_INFO_NULL, &mpifh);
    if (err != MPI_SUCCESS) {
      char str[128];
      sprintf(str,"Cannot open dump file %s",filecurrent);
      error->one(FLERR,str);
    }
    int myrank;
    MPI_Comm_rank(world,&myrank);
    if (myrank == 0)
      MPI_File_get_size(mpifh,&mpifo);
    MPI_Bcast(&mpifo, 1, MPI_LMP_BIGINT, 0, world);
    MPI_File_set_size(mpifh,mpifo+headerSize+sumFileSize);
    currentFileSize = mpifo+headerSize+sumFileSize;

  }
  else { // replace open

    int err = MPI_File_open( world, filecurrent, MPI_MODE_CREATE | MPI_MODE_WRONLY  , MPI_INFO_NULL, &mpifh);
    if (err != MPI_SUCCESS) {
      char str[128];
      sprintf(str,"Cannot open dump file %s",filecurrent);
      error->one(FLERR,str);
    }
    mpifo = 0;

    MPI_File_set_size(mpifh,(MPI_Offset) (headerSize+sumFileSize));
    currentFileSize = (headerSize+sumFileSize);

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
  MPI_Allreduce(&bnme,&ntotal,1,MPI_LMP_BIGINT,MPI_SUM,world);

  int nmax;
  MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);

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
  else pack(NULL);
  if (sort_flag) sort();

  // determine how much data needs to be written for setting the file size and prepocess it prior to writing
  performEstimate = 1;
  write_header(nheader);
  write_data(nme,buf);
  MPI_Bcast(&sumFileSize, 1, MPI_LMP_BIGINT, (nprocs-1), world);

  openfile();

  performEstimate = 0;
  write_header(nheader); // mpifo now points to end of header info

  // now actually write the data
  performEstimate = 0;
  write_data(nme,buf);

  if (multifile)    MPI_File_close(&mpifh);
  if (multifile) delete [] filecurrent;
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::init_style()
{
  // format = copy of default or user-specified line format

  delete [] format;
  char *str;
  if (format_line_user) str = format_line_user;
  else str = format_default;
  char *r_token;

  int n = strlen(str) + 1;
  format = new char[n];
  strcpy(format,str);

  // tokenize the format string and add space at end of each format element
  // if user-specified int/float format exists, use it instead
  // if user-specified column format exists, use it instead
  // lo priority = line, medium priority = int/float, hi priority = column

  char *ptr;
  for (int i = 0; i < size_one; i++) {
    r_token = format;
    if (i == 0) ptr = utils::strtok_r(r_token," \0",&r_token);
    else ptr = utils::strtok_r(NULL," \0",&r_token);
    if (ptr == NULL) error->all(FLERR,"Dump_modify format line is too short");
    delete [] vformat[i];

    if (format_column_user[i]) {
      vformat[i] = new char[strlen(format_column_user[i]) + 2];
      strcpy(vformat[i],format_column_user[i]);
    } else if (vtype[i] == Dump::INT && format_int_user) {
      vformat[i] = new char[strlen(format_int_user) + 2];
      strcpy(vformat[i],format_int_user);
    } else if (vtype[i] == Dump::DOUBLE && format_float_user) {
      vformat[i] = new char[strlen(format_float_user) + 2];
      strcpy(vformat[i],format_float_user);
    } else if (vtype[i] == Dump::BIGINT && format_bigint_user) {
      vformat[i] = new char[strlen(format_bigint_user) + 2];
      strcpy(vformat[i],format_bigint_user);
    } else {
      vformat[i] = new char[strlen(ptr) + 2];
      strcpy(vformat[i],ptr);
    }

    vformat[i] = strcat(vformat[i]," ");
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

  if (binary) write_choice = &DumpCustomMPIIO::write_binary;
  else write_choice = &DumpCustomMPIIO::write_string;

  // find current ptr for each compute,fix,variable
  // check that fix frequency is acceptable

  int icompute;
  for (int i = 0; i < ncompute; i++) {
    icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0) error->all(FLERR,"Could not find dump custom compute ID");
    compute[i] = modify->compute[icompute];
  }

  int ifix;
  for (int i = 0; i < nfix; i++) {
    ifix = modify->find_fix(id_fix[i]);
    if (ifix < 0) error->all(FLERR,"Could not find dump custom fix ID");
    fix[i] = modify->fix[ifix];
    if (nevery % modify->fix[ifix]->peratom_freq)
      error->all(FLERR,"Dump custom and fix not computed at compatible times");
  }

  int ivariable;
  for (int i = 0; i < nvariable; i++) {
    ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0)
      error->all(FLERR,"Could not find dump custom variable name");
    variable[i] = ivariable;
  }

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for dump custom does not exist");
  }
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

    headerBuffer = (char *) malloc((2*sizeof(bigint)) + (9*sizeof(int)) + (6*sizeof(double)));

    headerSize = 0;
    memcpy(&((char*)headerBuffer)[headerSize],&update->ntimestep,sizeof(bigint));
    headerSize += sizeof(bigint);

    memcpy(&((char*)headerBuffer)[headerSize],&ndump,sizeof(bigint));
    headerSize += sizeof(bigint);

    memcpy(&((char*)headerBuffer)[headerSize],&domain->triclinic,sizeof(int));
    headerSize += sizeof(int);

    memcpy(&((char*)headerBuffer)[headerSize],&domain->boundary[0][0],6*sizeof(int));
    headerSize += 6*sizeof(int);

    memcpy(&((char*)headerBuffer)[headerSize],&boxxlo,sizeof(double));
    headerSize += sizeof(double);

    memcpy(&((char*)headerBuffer)[headerSize],&boxxhi,sizeof(double));
    headerSize += sizeof(double);

    memcpy(&((char*)headerBuffer)[headerSize],&boxylo,sizeof(double));
    headerSize += sizeof(double);

    memcpy(&((char*)headerBuffer)[headerSize],&boxyhi,sizeof(double));
    headerSize += sizeof(double);

    memcpy(&((char*)headerBuffer)[headerSize],&boxzlo,sizeof(double));
    headerSize += sizeof(double);

    memcpy(&((char*)headerBuffer)[headerSize],&boxzhi,sizeof(double));
    headerSize += sizeof(double);

    memcpy(&((char*)headerBuffer)[headerSize],&size_one,sizeof(int));
    headerSize += sizeof(int);

    memcpy(&((char*)headerBuffer)[headerSize],&nprocs,sizeof(int));
    headerSize += sizeof(int);
  }
  else { // write data
    if (me == 0)
      MPI_File_write_at(mpifh,mpifo,headerBuffer,headerSize,MPI_BYTE,MPI_STATUS_IGNORE);
    mpifo += headerSize;
    free(headerBuffer);
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::header_binary_triclinic(bigint ndump)
{
  if (performEstimate) {

    headerBuffer = (char *) malloc((2*sizeof(bigint)) + (9*sizeof(int)) + (6*sizeof(double)));

    headerSize = 0;
    memcpy(&((char*)headerBuffer)[headerSize],&update->ntimestep,sizeof(bigint));
    headerSize += sizeof(bigint);

    memcpy(&((char*)headerBuffer)[headerSize],&ndump,sizeof(bigint));
    headerSize += sizeof(bigint);

    memcpy(&((char*)headerBuffer)[headerSize],&domain->triclinic,sizeof(int));
    headerSize += sizeof(int);

    memcpy(&((char*)headerBuffer)[headerSize],&domain->boundary[0][0],6*sizeof(int));
    headerSize += 6*sizeof(int);

    memcpy(&((char*)headerBuffer)[headerSize],&boxxlo,sizeof(double));
    headerSize += sizeof(double);

    memcpy(&((char*)headerBuffer)[headerSize],&boxxhi,sizeof(double));
    headerSize += sizeof(double);

    memcpy(&((char*)headerBuffer)[headerSize],&boxylo,sizeof(double));
    headerSize += sizeof(double);

    memcpy(&((char*)headerBuffer)[headerSize],&boxyhi,sizeof(double));
    headerSize += sizeof(double);

    memcpy(&((char*)headerBuffer)[headerSize],&boxzlo,sizeof(double));
    headerSize += sizeof(double);

    memcpy(&((char*)headerBuffer)[headerSize],&boxzhi,sizeof(double));
    headerSize += sizeof(double);

    memcpy(&((char*)headerBuffer)[headerSize],&boxxy,sizeof(double));
    headerSize += sizeof(double);

    memcpy(&((char*)headerBuffer)[headerSize],&boxxz,sizeof(double));
    headerSize += sizeof(double);

    memcpy(&((char*)headerBuffer)[headerSize],&boxyz,sizeof(double));
    headerSize += sizeof(double);

    memcpy(&((char*)headerBuffer)[headerSize],&size_one,sizeof(int));
    headerSize += sizeof(int);

    memcpy(&((char*)headerBuffer)[headerSize],&nprocs,sizeof(int));
    headerSize += sizeof(int);

  }
  else { // write data

    if (me == 0)
      MPI_File_write_at(mpifh,mpifo,headerBuffer,headerSize,MPI_BYTE,MPI_STATUS_IGNORE);
    mpifo += headerSize;
    free(headerBuffer);
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::header_item(bigint ndump)
{
  if (performEstimate) {

    headerBuffer = (char *) malloc(MAX_TEXT_HEADER_SIZE);

    headerSize = 0;
    headerSize += sprintf(((char*)&((char*)headerBuffer)[headerSize]),"ITEM: TIMESTEP\n");
    headerSize += sprintf(&((char*)headerBuffer)[headerSize],BIGINT_FORMAT "\n",update->ntimestep);
    headerSize += sprintf(&((char*)headerBuffer)[headerSize],"ITEM: NUMBER OF ATOMS\n");
    headerSize += sprintf(&((char*)headerBuffer)[headerSize],BIGINT_FORMAT "\n",ndump);
    headerSize += sprintf(&((char*)headerBuffer)[headerSize],"ITEM: BOX BOUNDS %s\n",boundstr);
    headerSize += sprintf(&((char*)headerBuffer)[headerSize],"%g %g\n",boxxlo,boxxhi);
    headerSize += sprintf(&((char*)headerBuffer)[headerSize],"%g %g\n",boxylo,boxyhi);
    headerSize += sprintf(&((char*)headerBuffer)[headerSize],"%g %g\n",boxzlo,boxzhi);
    headerSize += sprintf(&((char*)headerBuffer)[headerSize],"ITEM: ATOMS %s\n",columns);
  }
  else { // write data

    if (me == 0)
      MPI_File_write_at(mpifh,mpifo,headerBuffer,headerSize,MPI_CHAR,MPI_STATUS_IGNORE);
    mpifo += headerSize;
    free(headerBuffer);
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::header_item_triclinic(bigint ndump)
{
  if (performEstimate) {

    headerBuffer = (char *) malloc(MAX_TEXT_HEADER_SIZE);

    headerSize = 0;
    headerSize += sprintf(&((char*)headerBuffer)[headerSize],"ITEM: TIMESTEP\n");
    headerSize += sprintf(&((char*)headerBuffer)[headerSize],BIGINT_FORMAT "\n",update->ntimestep);
    headerSize += sprintf(&((char*)headerBuffer)[headerSize],"ITEM: NUMBER OF ATOMS\n");
    headerSize += sprintf(&((char*)headerBuffer)[headerSize],BIGINT_FORMAT "\n",ndump);
    headerSize += sprintf(&((char*)headerBuffer)[headerSize],"ITEM: BOX BOUNDS xy xz yz %s\n",boundstr);
    headerSize += sprintf(&((char*)headerBuffer)[headerSize],"%g %g %g\n",boxxlo,boxxhi,boxxy);
    headerSize += sprintf(&((char*)headerBuffer)[headerSize],"%g %g %g\n",boxylo,boxyhi,boxxz);
    headerSize += sprintf(&((char*)headerBuffer)[headerSize],"%g %g %g\n",boxzlo,boxzhi,boxyz);
    headerSize += sprintf(&((char*)headerBuffer)[headerSize],"ITEM: ATOMS %s\n",columns);
  }
  else { // write data

    if (me == 0)
      MPI_File_write_at(mpifh,mpifo,headerBuffer,headerSize,MPI_CHAR,MPI_STATUS_IGNORE);
    mpifo += headerSize;
    free(headerBuffer);
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::write_binary(int n, double *mybuf)
{
  n *= size_one;

  if (performEstimate) {

    bigint incPrefix = 0;
    bigint bigintNme = (bigint) nme;
    MPI_Scan(&bigintNme,&incPrefix,1,MPI_LMP_BIGINT,MPI_SUM,world);
    sumFileSize = (incPrefix*size_one*sizeof(double)) + (nprocs * sizeof(int));
    offsetFromHeader = ((incPrefix-bigintNme)*size_one*sizeof(double)) + (me * sizeof(int));
  }
  else {
    int byteBufSize = (n*sizeof(double)) + sizeof(int);

    char *bufWithSize;
    memory->create(bufWithSize,byteBufSize,"dump:bufWithSize");
    memcpy(bufWithSize,(char*)(&n),sizeof(int));
    memcpy(&((char*)bufWithSize)[sizeof(int)],mybuf,(n*sizeof(double)));
    MPI_File_write_at_all(mpifh,mpifo+offsetFromHeader,bufWithSize,byteBufSize,MPI_BYTE,MPI_STATUS_IGNORE);
    memory->destroy(bufWithSize);

    if (flush_flag)
      MPI_File_sync(mpifh);
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomMPIIO::write_string(int n, double *mybuf)
{
  if (performEstimate) {

#if defined(_OPENMP)
    int nthreads = omp_get_max_threads();
    if ((nthreads > 1) && !(lmp->kokkos))
      nsme = convert_string_omp(n,mybuf); // not (yet) compatible with Kokkos
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
  else {
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

int DumpCustomMPIIO::convert_string_omp(int n, double *mybuf)
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

#pragma omp parallel default(none) shared(bufOffset, bufRange, bufLength, mpifhStringCountPerThread, mpifh_buffer_line_per_thread, mybuf)
    {
      int tid = omp_get_thread_num();
      int m=0;

      for (int i = 0; i < bufRange[tid]; i++) {

        if ((bufLength[tid] - mpifhStringCountPerThread[tid]) < DUMP_BUF_INCREMENT_SIZE) {
          mpifh_buffer_line_per_thread[tid] = (char *) realloc(mpifh_buffer_line_per_thread[tid],(mpifhStringCountPerThread[tid]+DUMP_BUF_CHUNK_SIZE) * sizeof(char));
          bufLength[tid] = (mpifhStringCountPerThread[tid]+DUMP_BUF_CHUNK_SIZE) * sizeof(char);
        }
        for (int j = 0; j < size_one; j++) {

          if (vtype[j] == Dump::INT)
            mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),vformat[j],static_cast<int> (mybuf[bufOffset[tid]+m]));
          else if (vtype[j] == Dump::DOUBLE)
            mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),vformat[j],mybuf[bufOffset[tid]+m]);
          else if (vtype[j] == Dump::STRING)
            mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),vformat[j],typenames[(int) mybuf[bufOffset[tid]+m]]);
          m ++;
        }
        mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),"\n");

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
        maxsbuf = mpifhStringCount;
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
