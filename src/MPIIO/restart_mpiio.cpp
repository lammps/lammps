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

#include "mpi.h"
#include "restart_mpiio.h"
#include "error.h"
#include "limits.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RestartMPIIO::RestartMPIIO(LAMMPS *lmp) : Pointers(lmp)
{
  mpiio_exists = 1;
  MPI_Comm_size(world,&nprocs);
  MPI_Comm_rank(world,&myrank);
}

/* ----------------------------------------------------------------------
   calls MPI_File_open in read-only mode, read_restart should call this
   for some file servers it is most efficient to only read or only write
------------------------------------------------------------------------- */

void RestartMPIIO::openForRead(char *filename)
{
  int err = MPI_File_open(world, filename, MPI_MODE_RDONLY ,
                          MPI_INFO_NULL, &mpifh);
  if (err != MPI_SUCCESS) {
    char str[MPI_MAX_ERROR_STRING+128];
    char mpiErrorString[MPI_MAX_ERROR_STRING];
    int mpiErrorStringLength;
    MPI_Error_string(err, mpiErrorString, &mpiErrorStringLength);
    sprintf(str,"Cannot open restart file for reading - MPI error: %s",
            mpiErrorString);
    error->one(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   calls MPI_File_open in write-only mode, write_restart should call this
   for some file servers it is most efficient to only read or only write
------------------------------------------------------------------------- */

void RestartMPIIO::openForWrite(char *filename)
{
  int err = MPI_File_open(world, filename, MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &mpifh);
  if (err != MPI_SUCCESS) {
    char str[MPI_MAX_ERROR_STRING+128];
    char mpiErrorString[MPI_MAX_ERROR_STRING];
    int mpiErrorStringLength;
    MPI_Error_string(err, mpiErrorString, &mpiErrorStringLength);
    sprintf(str,"Cannot open restart file for writing - MPI error: %s",
            mpiErrorString);
    error->one(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   determine the absolute offset for the data to be written with
     MPI_Scan of the send sizes
   compute the file size based off the MPI_Scan send size value on the last rank
   set the filesize with ftruncate via MPI_File_set_size
   write the data via collective MPI-IO by calling MPI_File_write_at_all
------------------------------------------------------------------------- */

void RestartMPIIO::write(MPI_Offset headerOffset, int send_size, double *buf)
{
  MPI_Status mpiStatus;
  bigint incPrefix = 0;
  bigint bigintSendSize = (bigint) send_size;
  MPI_Scan(&bigintSendSize,&incPrefix,1,MPI_LMP_BIGINT,MPI_SUM,world);

  bigint largestIncPrefix = incPrefix;
  MPI_Bcast(&largestIncPrefix, 1, MPI_LMP_BIGINT, (nprocs-1), world);

  int err = MPI_File_set_size(mpifh,
                              (headerOffset+(largestIncPrefix*sizeof(double))));
  if (err != MPI_SUCCESS) {
    char str[MPI_MAX_ERROR_STRING+128];
    char mpiErrorString[MPI_MAX_ERROR_STRING];
    int mpiErrorStringLength;
    MPI_Error_string(err, mpiErrorString, &mpiErrorStringLength);
    sprintf(str,"Cannot set restart file size - MPI error: %s",
            mpiErrorString);
    error->one(FLERR,str);
  }

  err = MPI_File_write_at_all(mpifh,headerOffset +
                              ((incPrefix-bigintSendSize)*sizeof(double)),
                              buf,send_size,MPI_DOUBLE,&mpiStatus);
  if (err != MPI_SUCCESS) {
    char str[MPI_MAX_ERROR_STRING+128];
    char mpiErrorString[MPI_MAX_ERROR_STRING];
    int mpiErrorStringLength;
    MPI_Error_string(err, mpiErrorString, &mpiErrorStringLength);
    sprintf(str,"Cannot write to restart file - MPI error: %s",
            mpiErrorString);
    error->one(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   read the data into buf via collective MPI-IO by calling MPI_File_read_at_all
   with the chunkOffset and chunkSize provided
   if the consolidated chunksize is greater than INT_MAX
     can only happen in extreme situation of reading restart file on
     much fewer ranks than written and with relatively large data sizes
   follow the collective IO call with rank independant IO to read remaining data
------------------------------------------------------------------------- */

void RestartMPIIO::read(MPI_Offset chunkOffset, bigint chunkSize, double *buf)
{
  MPI_Status mpiStatus;

  int intChunkSize;
  bigint remainingSize = 0;
  if (chunkSize > INT_MAX) {
    intChunkSize = INT_MAX;
    remainingSize = chunkSize - INT_MAX;
  }
  else intChunkSize = (int) chunkSize;

  int err = MPI_File_read_at_all(mpifh,chunkOffset,buf,intChunkSize,
                                 MPI_DOUBLE,&mpiStatus);
  if (err != MPI_SUCCESS) {
    char str[MPI_MAX_ERROR_STRING+128];
    char mpiErrorString[MPI_MAX_ERROR_STRING];
    int mpiErrorStringLength;
    MPI_Error_string(err, mpiErrorString, &mpiErrorStringLength);
    sprintf(str,"Cannot read from restart file - MPI error: %s",
            mpiErrorString);
    error->one(FLERR,str);
  }

  MPI_Offset currentOffset = chunkOffset+intChunkSize;
  MPI_Offset bufOffset = intChunkSize;
  while (remainingSize > 0) {
    int currentChunkSize;
    if (remainingSize > INT_MAX) {
      currentChunkSize = INT_MAX;
      remainingSize -= INT_MAX;
    }
    else {
      currentChunkSize = remainingSize;
      remainingSize = 0;
    }
    int err = MPI_File_read_at(mpifh,currentOffset,&buf[bufOffset],
                               currentChunkSize,MPI_DOUBLE,&mpiStatus);
    if (err != MPI_SUCCESS) {
      char str[MPI_MAX_ERROR_STRING+128];
      char mpiErrorString[MPI_MAX_ERROR_STRING];
      int mpiErrorStringLength;
      MPI_Error_string(err, mpiErrorString, &mpiErrorStringLength);
      sprintf(str,"Cannot read from restart file - MPI error: %s",
              mpiErrorString);
      error->one(FLERR,str);
    }
    currentOffset += currentChunkSize;
    bufOffset += currentChunkSize;
  }
}

/* ----------------------------------------------------------------------
   calls MPI_File_close
------------------------------------------------------------------------- */

void RestartMPIIO::close()
{
  int err = MPI_File_close(&mpifh);
  if (err != MPI_SUCCESS) {
    char str[MPI_MAX_ERROR_STRING+128];
    char mpiErrorString[MPI_MAX_ERROR_STRING];
    int mpiErrorStringLength;
    MPI_Error_string(err, mpiErrorString, &mpiErrorStringLength);
    sprintf(str,"Cannot close restart file - MPI error: %s",mpiErrorString);
    error->one(FLERR,str);
  }
}
