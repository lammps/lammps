/* ----------------------------------------------------------------------
   CSlib - Client/server library for code coupling
   https://cslib.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright 2018 National Technology & Engineering Solutions of
   Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.
   This software is distributed under the modified Berkeley Software
   Distribution (BSD) License.

   See the README file in the top-level CSlib directory.
------------------------------------------------------------------------- */

// C style library interface to CSlib class

#ifdef MPI_YES
#include <mpi.h>
#else
#include <mpi_dummy.h>
#endif
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cslib_wrap.h"
#include "cslib.h"

using namespace CSLIB_NS;

// ----------------------------------------------------------------------

void cslib_open(int csflag, const char *mode, const void *ptr,
                const void *pcomm, void **csptr)
{
  CSlib *cs = new CSlib(csflag,mode,ptr,pcomm);
  *csptr = (void *) cs;
}

// ----------------------------------------------------------------------

void cslib_open_fortran(int csflag, const char *mode, const char *str,
                        const void *pcomm, void **csptr)
{
  MPI_Comm ccomm;
  void *pccomm = nullptr;

  if (pcomm) {
    MPI_Fint *fcomm = (MPI_Fint *) pcomm;
    ccomm = MPI_Comm_f2c(*fcomm);
    pccomm = &ccomm;
  }

  CSlib *cs = new CSlib(csflag,mode,str,pccomm);
  *csptr = (void *) cs;
}

// ----------------------------------------------------------------------

void cslib_open_fortran_mpi_one(int csflag, const char *mode,
                                const void *pboth, const void *pcomm,
                                void **csptr)
{
  MPI_Comm ccomm,cboth;
  void *pccomm,*pcboth;

  MPI_Fint *fcomm = (MPI_Fint *) pcomm;
  ccomm = MPI_Comm_f2c(*fcomm);
  pccomm = &ccomm;

  MPI_Fint *fboth = (MPI_Fint *) pboth;
  cboth = MPI_Comm_f2c(*fboth);
  pcboth = &cboth;

  CSlib *cs = new CSlib(csflag,mode,pcboth,pccomm);
  *csptr = (void *) cs;
}

// ----------------------------------------------------------------------

void cslib_close(void *ptr)
{
  CSlib *cs = (CSlib *) ptr;
  delete cs;
}

// ----------------------------------------------------------------------

void cslib_send(void *ptr, int msgID, int nfield)
{
  CSlib *cs = (CSlib *) ptr;
  cs->send(msgID,nfield);
}

// ----------------------------------------------------------------------

void cslib_pack_int(void *ptr, int id, int value)
{
  CSlib *cs = (CSlib *) ptr;
  cs->pack_int(id,value);
}

// ----------------------------------------------------------------------

void cslib_pack_int64(void *ptr, int id, int64_t value)
{
  CSlib *cs = (CSlib *) ptr;
  cs->pack_int64(id,value);
}

// ----------------------------------------------------------------------

void cslib_pack_float(void *ptr, int id, float value)
{
  CSlib *cs = (CSlib *) ptr;
  cs->pack_float(id,value);
}

// ----------------------------------------------------------------------

void cslib_pack_double(void *ptr, int id, double value)
{
  CSlib *cs = (CSlib *) ptr;
  cs->pack_double(id,value);
}

// ----------------------------------------------------------------------

void cslib_pack_string(void *ptr, int id, char *value)
{
  CSlib *cs = (CSlib *) ptr;
  cs->pack_string(id,value);
}

// ----------------------------------------------------------------------

void cslib_pack(void *ptr, int id, int ftype, int flen, void *data)
{
  CSlib *cs = (CSlib *) ptr;
  cs->pack(id,ftype,flen,data);
}

// ----------------------------------------------------------------------

void cslib_pack_parallel(void *ptr, int id, int ftype,
			 int nlocal, int *ids, int nper, void *data)
{
  CSlib *cs = (CSlib *) ptr;
  cs->pack_parallel(id,ftype,nlocal,ids,nper,data);
}

// ----------------------------------------------------------------------

int cslib_recv(void *ptr, int *nfield_caller,
	       int **fieldID_caller, int **fieldtype_caller,
	       int **fieldlen_caller)
{
  CSlib *cs = (CSlib *) ptr;

  int nfield;
  int *fieldID,*fieldtype,*fieldlen;
  int msgID = cs->recv(nfield,fieldID,fieldtype,fieldlen);

  *nfield_caller = nfield;
  *fieldID_caller = fieldID;
  *fieldtype_caller = fieldtype;
  *fieldlen_caller = fieldlen;

  return msgID;
}

// ----------------------------------------------------------------------

int cslib_unpack_int(void *ptr, int id)
{
  CSlib *cs = (CSlib *) ptr;
  return cs->unpack_int(id);
}
// ----------------------------------------------------------------------

int64_t cslib_unpack_int64(void *ptr, int id)
{
  CSlib *cs = (CSlib *) ptr;
  return cs->unpack_int64(id);
}

// ----------------------------------------------------------------------

float cslib_unpack_float(void *ptr, int id)
{
  CSlib *cs = (CSlib *) ptr;
  return cs->unpack_float(id);
}

// ----------------------------------------------------------------------

double cslib_unpack_double(void *ptr, int id)
{
  CSlib *cs = (CSlib *) ptr;
  return cs->unpack_double(id);
}

// ----------------------------------------------------------------------

char *cslib_unpack_string(void *ptr, int id)
{
  CSlib *cs = (CSlib *) ptr;
  return cs->unpack_string(id);
}

// ----------------------------------------------------------------------

void *cslib_unpack(void *ptr, int id)
{
  CSlib *cs = (CSlib *) ptr;
  return cs->unpack(id);
}

// ----------------------------------------------------------------------

void cslib_unpack_data(void *ptr, int id, void *data)
{
  CSlib *cs = (CSlib *) ptr;
  cs->unpack(id,data);
}

// ----------------------------------------------------------------------

void cslib_unpack_parallel(void *ptr, int id, int nlocal, int *ids,
			   int nper, void *data)
{
  CSlib *cs = (CSlib *) ptr;
  cs->unpack_parallel(id,nlocal,ids,nper,data);
}

// ----------------------------------------------------------------------

int cslib_extract(void *ptr, int flag)
{
  CSlib *cs = (CSlib *) ptr;
  return cs->extract(flag);
}
