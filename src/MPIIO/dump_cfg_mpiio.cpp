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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "dump_cfg_mpiio.h"
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "fix.h"
#include "variable.h"
#include "update.h"
#include "memory.h"
#include "error.h"

#ifdef LMP_USER_IO_TIMER
#include <sys/times.h>
#include <hwi/include/bqc/A2_inlines.h>
#include <stdlib.h>
long dumpCFGTimestamps[10];
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;

#define MAX_TEXT_HEADER_SIZE 4096
#define DUMP_BUF_CHUNK_SIZE 16384
#define DUMP_BUF_INCREMENT_SIZE 4096

enum{INT,DOUBLE,STRING,BIGINT};   // same as in DumpCustom

#define UNWRAPEXPAND 10.0
#define ONEFIELD 32
#define DELTA 1048576

/* ---------------------------------------------------------------------- */

DumpCFGMPIIO::DumpCFGMPIIO(LAMMPS *lmp, int narg, char **arg) :
  DumpCFG(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

DumpCFGMPIIO::~DumpCFGMPIIO()
{
  if (multifile == 0) MPI_File_close(&mpifh);
}

/* ---------------------------------------------------------------------- */

void DumpCFGMPIIO::openfile()
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

void DumpCFGMPIIO::write()
{

#ifdef LMP_USER_IO_TIMER
	long startTimeBase, endTimeBase;
	MPI_Barrier(world); // timestamp barrier
	if (me == 0)
		startTimeBase = GetTimeBase();
#endif

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

#ifdef LMP_USER_IO_TIMER
	MPI_Barrier(world); // timestamp barrier
	dumpCFGTimestamps[0] = GetTimeBase();
#endif

  openfile();

#ifdef LMP_USER_IO_TIMER
	MPI_Barrier(world); // timestamp barrier
	dumpCFGTimestamps[1] = GetTimeBase();
#endif

  performEstimate = 0;
  write_header(nheader); // mpifo now points to end of header info

#ifdef LMP_USER_IO_TIMER
	MPI_Barrier(world); // timestamp barrier
	dumpCFGTimestamps[2] = GetTimeBase();
#endif

  // now actually write the data
  performEstimate = 0;
  write_data(nme,buf);

#ifdef LMP_USER_IO_TIMER
	MPI_Barrier(world); // timestamp barrier
	dumpCFGTimestamps[3] = GetTimeBase();
#endif

  if (multifile)    MPI_File_close(&mpifh);
  if (multifile) delete [] filecurrent;

#ifdef LMP_USER_IO_TIMER
	MPI_Barrier(world); // timestamp barrier
	dumpCFGTimestamps[4] = GetTimeBase();
	if (me == 0) {
		endTimeBase = GetTimeBase();
		printf("total dump cycles: %ld - estimates and setup: %ld openfile: %ld write header: %ld write data: %ld close file: %ld\n",(long) (endTimeBase-startTimeBase),(long) (dumpCFGTimestamps[0]-startTimeBase),(long) (dumpCFGTimestamps[1]-dumpCFGTimestamps[0]),(long) (dumpCFGTimestamps[2]-dumpCFGTimestamps[1]),(long) (dumpCFGTimestamps[3]-dumpCFGTimestamps[2]),(long) (dumpCFGTimestamps[4]-dumpCFGTimestamps[3]));
	}
#endif

}

/* ---------------------------------------------------------------------- */

void DumpCFGMPIIO::init_style()
{
  if (multifile == 0 && !multifile_override)
    error->all(FLERR,"Dump cfg requires one snapshot per file");

  DumpCustom::init_style();

  // setup function ptrs

  write_choice = &DumpCFGMPIIO::write_string;

}

/* ---------------------------------------------------------------------- */

void DumpCFGMPIIO::write_header(bigint n)
{
  // set scale factor used by AtomEye for CFG viz
  // default = 1.0
  // for peridynamics, set to pre-computed PD scale factor
  //   so PD particles mimic C atoms
  // for unwrapped coords, set to UNWRAPEXPAND (10.0)
  //   so molecules are not split across periodic box boundaries

  if (performEstimate) {

    headerBuffer = (char *) malloc(MAX_TEXT_HEADER_SIZE);

    headerSize = 0;

    double scale = 1.0;
    if (atom->peri_flag) scale = atom->pdscale;
    else if (unwrapflag == 1) scale = UNWRAPEXPAND;

    char str[64];
  
    sprintf(str,"Number of particles = %s\n",BIGINT_FORMAT);
    headerSize += sprintf(((char*)&((char*)headerBuffer)[headerSize]),str,n);
    headerSize += sprintf(((char*)&((char*)headerBuffer)[headerSize]),"A = %g Angstrom (basic length-scale)\n",scale);
    headerSize += sprintf(((char*)&((char*)headerBuffer)[headerSize]),"H0(1,1) = %g A\n",domain->xprd);
    headerSize += sprintf(((char*)&((char*)headerBuffer)[headerSize]),"H0(1,2) = 0 A \n");
    headerSize += sprintf(((char*)&((char*)headerBuffer)[headerSize]),"H0(1,3) = 0 A \n");
    headerSize += sprintf(((char*)&((char*)headerBuffer)[headerSize]),"H0(2,1) = %g A \n",domain->xy);
    headerSize += sprintf(((char*)&((char*)headerBuffer)[headerSize]),"H0(2,2) = %g A\n",domain->yprd);
    headerSize += sprintf(((char*)&((char*)headerBuffer)[headerSize]),"H0(2,3) = 0 A \n");
    headerSize += sprintf(((char*)&((char*)headerBuffer)[headerSize]),"H0(3,1) = %g A \n",domain->xz);
    headerSize += sprintf(((char*)&((char*)headerBuffer)[headerSize]),"H0(3,2) = %g A \n",domain->yz);
    headerSize += sprintf(((char*)&((char*)headerBuffer)[headerSize]),"H0(3,3) = %g A\n",domain->zprd);
    headerSize += sprintf(((char*)&((char*)headerBuffer)[headerSize]),".NO_VELOCITY.\n");
    headerSize += sprintf(((char*)&((char*)headerBuffer)[headerSize]),"entry_count = %d\n",nfield-2);
    for (int i = 0; i < nfield-5; i++)
      headerSize += sprintf(((char*)&((char*)headerBuffer)[headerSize]),"auxiliary[%d] = %s\n",i,auxname[i]);
  }
  else { // write data

    if (me == 0)
      MPI_File_write_at(mpifh,mpifo,headerBuffer,headerSize,MPI_CHAR,MPI_STATUS_IGNORE);
    mpifo += headerSize;
    free(headerBuffer);
  }
}

#if defined(_OPENMP)

/* ----------------------------------------------------------------------
   convert mybuf of doubles to one big formatted string in sbuf
   return -1 if strlen exceeds an int, since used as arg in MPI calls in Dump
------------------------------------------------------------------------- */

int DumpCFGMPIIO::convert_string_omp(int n, double *mybuf)
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

      if (unwrapflag == 0) {
        for (int i = 0; i < bufRange[tid]; i++) {
          if ((bufLength[tid] - mpifhStringCountPerThread[tid]) < DUMP_BUF_INCREMENT_SIZE) {
            mpifh_buffer_line_per_thread[tid] = (char *) realloc(mpifh_buffer_line_per_thread[tid],(mpifhStringCountPerThread[tid]+DUMP_BUF_CHUNK_SIZE) * sizeof(char));
            bufLength[tid] = (mpifhStringCountPerThread[tid]+DUMP_BUF_CHUNK_SIZE) * sizeof(char);
          }
          for (int j = 0; j < size_one; j++) {  
            if (j == 0) {
              mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),"%f \n",(mybuf[bufOffset[tid]+m]));
            } else if (j == 1) {
              mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),"%s \n",typenames[(int) mybuf[bufOffset[tid]+m]]);
            } else if (j >= 2) {
            if (vtype[j] == INT) 
              mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),vformat[j],static_cast<int> (mybuf[bufOffset[tid]+m]));
            else if (vtype[j] == DOUBLE) 
              mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),vformat[j],mybuf[bufOffset[tid]+m]);
            else if (vtype[j] == STRING) 
              mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),vformat[j],typenames[(int) mybuf[bufOffset[tid]+m]]);
            else if (vtype[j] == BIGINT) 
              mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),vformat[j],static_cast<bigint> (mybuf[bufOffset[tid]+m]));
          }
          m++;
        } // for j
        mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),"\n");        
      } // for i
    } // wrap flag
    else if (unwrapflag == 1) {
      for (int i = 0; i < bufRange[tid]; i++) {
        if ((bufLength[tid] - mpifhStringCountPerThread[tid]) < DUMP_BUF_INCREMENT_SIZE) {
          mpifh_buffer_line_per_thread[tid] = (char *) realloc(mpifh_buffer_line_per_thread[tid],(mpifhStringCountPerThread[tid]+DUMP_BUF_CHUNK_SIZE) * sizeof(char));
          bufLength[tid] = (mpifhStringCountPerThread[tid]+DUMP_BUF_CHUNK_SIZE) * sizeof(char);
        }
        for (int j = 0; j < size_one; j++) {  
          double unwrap_coord;
          if (j == 0) {
	  //offset += sprintf(&sbuf[offset],"%f \n",mybuf[m]);
            mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),"%f \n",mybuf[bufOffset[tid]+m]);
          } else if (j == 1) {
	 // offset += sprintf(&sbuf[offset],"%s \n",typenames[(int) mybuf[m]]);
            mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),"%s \n",typenames[(int) mybuf[bufOffset[tid]+m]]);
          } else if (j >= 2 && j <= 4) {
            unwrap_coord = (mybuf[bufOffset[tid]+m] - 0.5)/UNWRAPEXPAND + 0.5;
          //offset += sprintf(&sbuf[offset],vformat[j],unwrap_coord);
            mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),vformat[j],unwrap_coord);
          } else if (j >= 5 ) {
            if (vtype[j] == INT) 
            //offset += 
            //  sprintf(&sbuf[offset],vformat[j],static_cast<int> (mybuf[m]));
              mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),vformat[j],static_cast<int> (mybuf[bufOffset[tid]+m]));
            else if (vtype[j] == DOUBLE) 
            // offset += sprintf(&sbuf[offset],vformat[j],mybuf[m]);
              mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),vformat[j],mybuf[bufOffset[tid]+m]);
            else if (vtype[j] == STRING) 
            // offset += 
            //  sprintf(&sbuf[offset],vformat[j],typenames[(int) mybuf[m]]);
              mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),vformat[j],typenames[(int) mybuf[bufOffset[tid]+m]]);
            else if (vtype[j] == BIGINT) 
            // offset += 
            //  sprintf(&sbuf[offset],vformat[j],static_cast<bigint> (mybuf[m]));
              mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),vformat[j],static_cast<bigint> (mybuf[bufOffset[tid]+m]));
          }
          m++;
        } // for j
        mpifhStringCountPerThread[tid] += sprintf(&(mpifh_buffer_line_per_thread[tid][mpifhStringCountPerThread[tid]]),"\n");        
      } // for i
    } // unwrap flag
  } // pragma omp parallel
  
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
  
  } // else omp

  return mpifhStringCount;
 
}

#endif

/* ---------------------------------------------------------------------- */

void DumpCFGMPIIO::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpCFGMPIIO::write_string(int n, double *mybuf)
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
  else {
    MPI_File_write_at_all(mpifh,mpifo+offsetFromHeader,sbuf,nsme,MPI_CHAR,MPI_STATUS_IGNORE);
    if (flush_flag)
      MPI_File_sync(mpifh);
  }
}

