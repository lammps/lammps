/* ----------------------------------------------------------------------
   CSlib - Client/server library for code coupling
   http://cslib.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright 2018 National Technology & Engineering Solutions of
   Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.
   This software is distributed under the modified Berkeley Software
   Distribution (BSD) License.

   See the README file in the top-level CSlib directory.
------------------------------------------------------------------------- */

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>

#include "cslib.h"
#include "msg_file.h"
#include "msg_zmq.h"
#include "msg_mpi_one.h"
#include "msg_mpi_two.h"

using namespace CSLIB_NS;

#define MAXTYPE 5       // # of defined field data types

/* ---------------------------------------------------------------------- */

CSlib::CSlib(int csflag, const char *mode, const void *ptr, const void *pcomm)
{
  if (pcomm) myworld = (uint64_t) *((MPI_Comm *) pcomm);
  else myworld = 0;

#ifdef MPI_NO
  if (pcomm) 
    error_all("constructor(): CSlib invoked with MPI_Comm "
              "but built w/out MPI support");
#endif
#ifdef MPI_YES              // NOTE: this could be OK to allow ??
                            // would allow a parallel app to invoke CSlib
                            //   in parallel and/or in serial
  if (!pcomm) 
    error_all("constructor(): CSlib invoked w/out MPI_Comm "
              "but built with MPI support");
#endif

  client = server = 0;
  if (csflag == 0) client = 1;
  else if (csflag == 1) server = 1;
  else error_all("constructor(): Invalid client/server arg");

  if (pcomm == NULL) {
    me = 0;
    nprocs = 1;

    if (strcmp(mode,"file") == 0) msg = new MsgFile(csflag,ptr);
    else if (strcmp(mode,"zmq") == 0) msg = new MsgZMQ(csflag,ptr);
    else if (strcmp(mode,"mpi/one") == 0) 
      error_all("constructor(): No mpi/one mode for serial lib usage");
    else if (strcmp(mode,"mpi/two") == 0)
      error_all("constructor(): No mpi/two mode for serial lib usage");
    else error_all("constructor(): Unknown mode");

  } else if (pcomm) {
    MPI_Comm world = (MPI_Comm) myworld;
    MPI_Comm_rank(world,&me);
    MPI_Comm_size(world,&nprocs);

    if (strcmp(mode,"file") == 0) msg = new MsgFile(csflag,ptr,world);
    else if (strcmp(mode,"zmq") == 0) msg = new MsgZMQ(csflag,ptr,world);
    else if (strcmp(mode,"mpi/one") == 0) msg = new MsgMPIOne(csflag,ptr,world);
    else if (strcmp(mode,"mpi/two") == 0) msg = new MsgMPITwo(csflag,ptr,world);
    else error_all("constructor(): Unknown mode");
  }

  maxfield = 0;
  fieldID = fieldtype = fieldlen = fieldoffset = NULL;
  maxheader = 0;
  header = NULL;
  maxbuf = 0;
  buf = NULL;

  recvcounts = displs = NULL;
  maxglobal = 0;
  allids = NULL;
  maxfieldbytes = 0;
  fielddata = NULL;
  
  pad = "\0\0\0\0\0\0\0";    // just length 7 since will have trailing NULL
  
  nsend = nrecv = 0;
}

/* ---------------------------------------------------------------------- */

CSlib::~CSlib()
{
  deallocate_fields();
  sfree(header);
  sfree(buf);
  
  sfree(recvcounts);
  sfree(displs);
  sfree(allids);
  sfree(fielddata);

  delete msg;
}

/* ---------------------------------------------------------------------- */

void CSlib::send(int msgID_caller, int nfield_caller)
{
  if (nfield_caller < 0) error_all("send(): Invalid nfield");

  msgID = msgID_caller;
  nfield = nfield_caller;
  allocate_fields();

  fieldcount = 0;
  nbuf = 0;
  
  if (fieldcount == nfield) send_message();
}

/* ---------------------------------------------------------------------- */

void CSlib::pack_int(int id, int value)
{
  pack(id,1,1,&value);
}

/* ---------------------------------------------------------------------- */

void CSlib::pack_int64(int id, int64_t value)
{
  pack(id,2,1,&value);
}

/* ---------------------------------------------------------------------- */

void CSlib::pack_float(int id, float value)
{
  pack(id,3,1,&value);
}

/* ---------------------------------------------------------------------- */

void CSlib::pack_double(int id, double value)
{
  pack(id,4,1,&value);
}

/* ---------------------------------------------------------------------- */

void CSlib::pack_string(int id, char *value)
{
  pack(id,5,strlen(value)+1,value);
}

/* ---------------------------------------------------------------------- */

void CSlib::pack(int id, int ftype, int flen, void *data)
{
  if (find_field(id,fieldcount) >= 0)
    error_all("pack(): Reuse of field ID");
  if (ftype < 1 || ftype > MAXTYPE) error_all("pack(): Invalid ftype");
  if (flen < 0) error_all("pack(): Invalid flen");
    
  fieldID[fieldcount] = id;
  fieldtype[fieldcount] = ftype;
  fieldlen[fieldcount] = flen;

  int nbytes,nbytesround;
  onefield(ftype,flen,nbytes,nbytesround);

  memcpy(&buf[nbuf],data,nbytes);
  memcpy(&buf[nbuf+nbytes],pad,nbytesround-nbytes);
  nbuf += nbytesround;
  
  fieldcount++;
  if (fieldcount == nfield) send_message();
}

/* ---------------------------------------------------------------------- */

void CSlib::pack_parallel(int id, int ftype,
			  int nlocal, int *ids, int nper, void *data)
{
  int i,j,k,m;

  if (find_field(id,fieldcount) >= 0)
    error_all("pack_parallel(): Reuse of field ID");
  if (ftype < 1 || ftype > MAXTYPE) error_all("pack_parallel(): Invalid ftype");
  if (nlocal < 0) error_all("pack_parallel(): Invalid nlocal");
  if (nper < 1) error_all("pack_parallel(): Invalid nper");

  MPI_Comm world = (MPI_Comm) myworld;

  // NOTE: check for overflow of maxglobal and flen

  int nglobal;
  MPI_Allreduce(&nlocal,&nglobal,1,MPI_INT,MPI_SUM,world);
  int flen = nper*nglobal;

  fieldID[fieldcount] = id;
  fieldtype[fieldcount] = ftype;
  fieldlen[fieldcount] = flen;
  
  // nlocal datums, each of nper length, from all procs
  // final data in buf = datums for all natoms, ordered by ids

  if (recvcounts == NULL) {
    recvcounts = (int *) smalloc(nprocs*sizeof(int));
    displs = (int *) smalloc(nprocs*sizeof(int));
  }

  MPI_Allgather(&nlocal,1,MPI_INT,recvcounts,1,MPI_INT,world);

  displs[0] = 0;
  for (int iproc = 1; iproc < nprocs; iproc++)
    displs[iproc] = displs[iproc-1] + recvcounts[iproc-1];

  if (ids && nglobal > maxglobal) {
    sfree(allids);
    maxglobal = nglobal;
    // NOTE: maxglobal*sizeof(int) could overflow int
    allids = (int *) smalloc(maxglobal*sizeof(int));
  }

  MPI_Allgatherv(ids,nlocal,MPI_INT,allids,
                 recvcounts,displs,MPI_INT,world);
  
  int nlocalsize = nper*nlocal;
  MPI_Allgather(&nlocalsize,1,MPI_INT,recvcounts,1,MPI_INT,world);

  displs[0] = 0;
  for (int iproc = 1; iproc < nprocs; iproc++)
    displs[iproc] = displs[iproc-1] + recvcounts[iproc-1];

  int nbytes,nbytesround;
  onefield(ftype,flen,nbytes,nbytesround);

  if (ftype == 1) {
    int *alldata;
    if (ids) {
      if (nbytes > maxfieldbytes) {
        sfree(fielddata);
        maxfieldbytes = nbytes;   
        fielddata = (char *) smalloc(maxfieldbytes);
      }
      alldata = (int *) fielddata;
    } else alldata = (int *) &buf[nbuf];
    MPI_Allgatherv(data,nlocalsize,MPI_INT,alldata,
		   recvcounts,displs,MPI_INT,world);
    if (ids) {
      int *bufptr = (int *) &buf[nbuf];
      m = 0;
      for (i = 0; i < nglobal; i++) {
	j = (allids[i]-1) * nper;
	if (nper == 1) bufptr[j] = alldata[m++];
	else
	  for (k = 0; k < nper; k++)
	    bufptr[j++] = alldata[m++];
      }
    }

  } else if (ftype == 2) {
    int64_t *alldata;
    if (ids) {
      if (nbytes > maxfieldbytes) {
        sfree(fielddata);
        maxfieldbytes = nbytes;   
        fielddata = (char *) smalloc(maxfieldbytes);
      }
      alldata = (int64_t *) fielddata;
    } else alldata = (int64_t *) &buf[nbuf];
    // NOTE: may be just MPI_LONG on some machines
    MPI_Allgatherv(data,nlocalsize,MPI_LONG_LONG,alldata,
		   recvcounts,displs,MPI_LONG_LONG,world);
    if (ids) {
      int64_t *bufptr = (int64_t *) &buf[nbuf];
      m = 0;
      for (i = 0; i < nglobal; i++) {
	j = (allids[i]-1) * nper;
	if (nper == 1) bufptr[j] = alldata[m++];
	else
	  for (k = 0; k < nper; k++)
	    bufptr[j++] = alldata[m++];
      }
    }
    
  } else if (ftype == 3) {
    float *alldata;
    if (ids) {
      if (nbytes > maxfieldbytes) {
        sfree(fielddata);
        maxfieldbytes = nbytes;   
        fielddata = (char *) smalloc(maxfieldbytes);
      }
      alldata = (float *) fielddata;
    } else alldata = (float *) &buf[nbuf];
    MPI_Allgatherv(data,nlocalsize,MPI_FLOAT,alldata,
                   recvcounts,displs,MPI_FLOAT,world);
    if (ids) {
      float *bufptr = (float *) &buf[nbuf];
      m = 0;
      for (i = 0; i < nglobal; i++) {
	j = (allids[i]-1) * nper;
	if (nper == 1) bufptr[j] = alldata[m++];
	else
	  for (k = 0; k < nper; k++)
	    bufptr[j++] = alldata[m++];
      }
    }

  } else if (ftype == 4) {
    double *alldata;
    if (ids) {
      if (nbytes > maxfieldbytes) {
        sfree(fielddata);
        maxfieldbytes = nbytes;   
        fielddata = (char *) smalloc(maxfieldbytes);
      }
      alldata = (double *) fielddata;
    } else alldata = (double *) &buf[nbuf];
    MPI_Allgatherv(data,nlocalsize,MPI_DOUBLE,alldata,
                   recvcounts,displs,MPI_DOUBLE,world);
    if (ids) {
      double *bufptr = (double *) &buf[nbuf];
      m = 0;
      for (i = 0; i < nglobal; i++) {
	j = (allids[i]-1) * nper;
	if (nper == 1) bufptr[j] = alldata[m++];
	else
	  for (k = 0; k < nper; k++)
	    bufptr[j++] = alldata[m++];
      }
    }

    /* eventually ftype = BYTE, but not yet
  } else if (ftype == 5) {
    char *alldata;
    if (ids) {
      if (nbytes > maxfieldbytes) {
        sfree(fielddata);
        maxfieldbytes = nbytes;   
        fielddata = (char *) smalloc(maxfieldbytes);
      }
      alldata = (char *) fielddata;
    } else alldata = (char *) &buf[nbuf];
    MPI_Allgatherv(data,nlocalsize,MPI_CHAR,alldata,
                   recvcounts,displs,MPI_CHAR,world);
    if (ids) {
      char *bufptr = (char *) &buf[nbuf];
      m = 0;
      for (i = 0; i < nglobal; i++) {
	j = (allids[i]-1) * nper;
	memcpy(&bufptr[j],&alldata[m],nper);
	m += nper;
      }
    }
    */
  }

  memcpy(&buf[nbuf+nbytes],pad,nbytesround-nbytes);
  nbuf += nbytesround;

  fieldcount++;
  if (fieldcount == nfield) send_message();
}

/* ---------------------------------------------------------------------- */

void CSlib::send_message()
{
  // setup header message

  int m = 0;
  header[m++] = msgID;
  header[m++] = nfield;
  for (int ifield = 0; ifield < nfield; ifield++) {
    header[m++] = fieldID[ifield];
    header[m++] = fieldtype[ifield];
    header[m++] = fieldlen[ifield];
  }

  msg->send(nheader,header,nbuf,buf);
  nsend++;
}

/* ---------------------------------------------------------------------- */

int CSlib::recv(int &nfield_caller, int *&fieldID_caller, 
		int *&fieldtype_caller, int *&fieldlen_caller)
{
  msg->recv(maxheader,header,maxbuf,buf);
  nrecv++;

  // unpack header message
  
  int m = 0;
  msgID = header[m++];
  nfield = header[m++];
  allocate_fields();

  int nbytes,nbytesround;

  nbuf = 0;
  for (int ifield = 0; ifield < nfield; ifield++) {
    fieldID[ifield] = header[m++];
    fieldtype[ifield] = header[m++];
    fieldlen[ifield] = header[m++];
    fieldoffset[ifield] = nbuf;
    onefield(fieldtype[ifield],fieldlen[ifield],nbytes,nbytesround);
    nbuf += nbytesround;
  }
  
  // return message parameters

  nfield_caller = nfield;
  fieldID_caller = fieldID;
  fieldtype_caller = fieldtype;
  fieldlen_caller = fieldlen;

  return msgID;
}

/* ---------------------------------------------------------------------- */

int CSlib::unpack_int(int id)
{
  int ifield = find_field(id,nfield);
  if (ifield < 0) error_all("unpack_int(): Unknown field ID");
  if (fieldtype[ifield] != 1) error_all("unpack_int(): Mis-match of ftype");
  if (fieldlen[ifield] != 1) error_all("unpack_int(): Flen is not 1");

  int *ptr = (int *) unpack(id);
  return *ptr;
}

/* ---------------------------------------------------------------------- */

int64_t CSlib::unpack_int64(int id)
{
  int ifield = find_field(id,nfield);
  if (ifield < 0) error_all("unpack_int64(): Unknown field ID");
  if (fieldtype[ifield] != 2) error_all("unpack_int64(): Mis-match of ftype");
  if (fieldlen[ifield] != 1) error_all("unpack_int64(): Flen is not 1");

  int64_t *ptr = (int64_t *) unpack(id);
  return *ptr;
}

/* ---------------------------------------------------------------------- */

float CSlib::unpack_float(int id)
{
  int ifield = find_field(id,nfield);
  if (ifield < 0) error_all("unpack_float(): Unknown field ID");
  if (fieldtype[ifield] != 3) error_all("unpack_float(): Mis-match of ftype");
  if (fieldlen[ifield] != 1) error_all("unpack_float(): Flen is not 1");

  float *ptr = (float *) unpack(id);
  return *ptr;
}

/* ---------------------------------------------------------------------- */

double CSlib::unpack_double(int id)
{
  int ifield = find_field(id,nfield);
  if (ifield < 0) error_all("unpack_double(): Unknown field ID");
  if (fieldtype[ifield] != 4) error_all("unpack_double(): Mis-match of ftype");
  if (fieldlen[ifield] != 1) error_all("unpack_double(): Flen is not 1");

  double *ptr = (double *) unpack(id);
  return *ptr;
}

/* ---------------------------------------------------------------------- */

char *CSlib::unpack_string(int id)
{
  int ifield = find_field(id,nfield);
  if (ifield < 0) error_all("unpack_string(): Unknown field ID");
  if (fieldtype[ifield] != 5) error_all("unpack_string(): Mis-match of ftype");

  char *ptr = (char *) unpack(id);
  return ptr;
}

/* ---------------------------------------------------------------------- */

void *CSlib::unpack(int id)
{
  int ifield = find_field(id,nfield);
  if (ifield < 0) error_all("unpack(): Unknown field ID");
  return &buf[fieldoffset[ifield]];
}

/* ---------------------------------------------------------------------- */

void CSlib::unpack(int id, void *data)
{
  int ifield = find_field(id,nfield);
  if (ifield < 0) error_all("unpack(): Unknown field ID");
  
  int ftype = fieldtype[ifield];
  int nbytes = fieldlen[ifield];
  if (ftype == 1) nbytes *= sizeof(int);
  else if (ftype == 2) nbytes *= sizeof(int64_t);
  else if (ftype == 3) nbytes *= sizeof(float);
  else if (ftype == 4) nbytes *= sizeof(double);
  memcpy(data,&buf[fieldoffset[ifield]],nbytes);
}

/* ---------------------------------------------------------------------- */

void CSlib::unpack_parallel(int id, int nlocal, int *ids, int nper, void *data)
{
  int i,j,k,m;

  int ifield = find_field(id,nfield);
  if (ifield < 0) error_all("unpack_parallel(): Unknown field ID");
  if (nlocal < 0) error_all("unpack_parallel(): Invalid nlocal");
  if (nper < 1) error_all("pack_parallel(): Invalid nper");

  MPI_Comm world = (MPI_Comm) myworld;

  int upto;
  if (!ids) {
    MPI_Scan(&nlocal,&upto,1,MPI_INT,MPI_SUM,world);
    upto -= nlocal;
  }
  
  if (fieldtype[ifield] == 1) {
    int *local = (int *) data;
    int *global = (int *) &buf[fieldoffset[ifield]];
    if (!ids) memcpy(local,&global[nper*upto],nper*nlocal*sizeof(int));
    else {
      m = 0;
      for (i = 0; i < nlocal; i++) {
	j = (ids[i]-1) * nper;
	if (nper == 1) local[m++] = global[j];
	else
	  for (k = 0; k < nper; k++)
	    local[m++] = global[j++];
      }
    } 

  } else if (fieldtype[ifield] == 2) {
    int64_t *local = (int64_t *) data;
    int64_t *global = (int64_t *) &buf[fieldoffset[ifield]];
    if (!ids) memcpy(local,&global[nper*upto],nper*nlocal*sizeof(int64_t));
    else {
      m = 0;
      for (i = 0; i < nlocal; i++) {
	j = (ids[i]-1) * nper;
	if (nper == 1) local[m++] = global[j];
	else
	  for (k = 0; k < nper; k++)
	    local[m++] = global[j++];
      }
    }

  } else if (fieldtype[ifield] == 3) {
    float *local = (float *) data;
    float *global = (float *) &buf[fieldoffset[ifield]];
    if (!ids) memcpy(local,&global[nper*upto],nper*nlocal*sizeof(float));
    else {
      m = 0;
      for (i = 0; i < nlocal; i++) {
	j = (ids[i]-1) * nper;
	if (nper == 1) local[m++] = global[j];
	else
	  for (k = 0; k < nper; k++)
	    local[m++] = global[j++];
      }
    }
    
  } else if (fieldtype[ifield] == 4) {
    double *local = (double *) data;
    double *global = (double *) &buf[fieldoffset[ifield]];
    if (!ids) memcpy(local,&global[nper*upto],nper*nlocal*sizeof(double));
    else {
      m = 0;
      for (i = 0; i < nlocal; i++) {
	j = (ids[i]-1) * nper;
	if (nper == 1) local[m++] = global[j];
	else
	  for (k = 0; k < nper; k++)
	    local[m++] = global[j++];
      }
    }
    
    /* eventually ftype = BYTE, but not yet
  } else if (fieldtype[ifield] == 5) {
    char *local = (char *) data;
    char *global = (char *) &buf[fieldoffset[ifield]];
    if (!ids) memcpy(local,&global[nper*upto],nper*nlocal*sizeof(char));
    else {
      m = 0;
      for (i = 0; i < nlocal; i++) {
	j = (ids[i]-1) * nper;
	memcpy(&local[m],&global[j],nper);
	m += nper;
      }
    }
    */
  }
}

/* ---------------------------------------------------------------------- */

int CSlib::extract(int flag)
{
  if (flag == 1) return nsend;
  if (flag == 2) return nrecv;
  error_all("extract(): Invalid flag");
  return 0;
}

/* ---------------------------------------------------------------------- */

void CSlib::onefield(int ftype, int flen, int &nbytes, int &nbytesround)
{
  int64_t bigbytes,bigbytesround;
  int64_t biglen = flen;
  
  if (ftype == 1) bigbytes = biglen * sizeof(int);
  else if (ftype == 2) bigbytes = biglen * sizeof(int64_t);
  else if (ftype == 3) bigbytes = biglen * sizeof(float);
  else if (ftype == 4) bigbytes = biglen * sizeof(double);
  else if (ftype == 5) bigbytes = biglen * sizeof(char);
  bigbytesround = roundup(bigbytes,8);

  if (nbuf + bigbytesround > INT_MAX)
    error_all("pack(): Message size exceeds 32-bit integer limit");

  nbytes = (int) bigbytes;
  nbytesround = (int) bigbytesround;
  if (nbuf + nbytesround > maxbuf) {
    maxbuf = nbuf + nbytesround;
    buf = (char *) srealloc(buf,maxbuf);
  }
}

/* ---------------------------------------------------------------------- */

int CSlib::find_field(int id, int n)
{
  int ifield;
  for (ifield = 0; ifield < n; ifield++)
    if (id == fieldID[ifield]) return ifield;
  return -1;
}

/* ---------------------------------------------------------------------- */

void CSlib::allocate_fields()
{
  int64_t bigbytes = (2 + 3*((int64_t) nfield)) * sizeof(int);
  if (bigbytes > INT_MAX)
    error_all("send(): Message header size exceeds 32-bit integer limit");

  nheader = 2;
  nheader += 3 * nfield;
  
  if (nfield > maxfield) {
    deallocate_fields();
    maxfield = nfield;
    fieldID = new int[maxfield];
    fieldtype = new int[maxfield];
    fieldlen = new int[maxfield];
    fieldoffset = new int[maxfield];
  }
  
  if (nheader > maxheader) {
    sfree(header);
    maxheader = nheader;
    header = (int *) smalloc(maxheader*sizeof(int));
  }
}

/* ---------------------------------------------------------------------- */

void CSlib::deallocate_fields()
{
  delete [] fieldID;
  delete [] fieldtype;
  delete [] fieldlen;
  delete [] fieldoffset;
}

/* ---------------------------------------------------------------------- */

void *CSlib::smalloc(int nbytes)
{
  if (nbytes == 0) return NULL;
  void *ptr = malloc(nbytes);
  if (ptr == NULL) {
    char str[128];
    sprintf(str,"malloc(): Failed to allocate %d bytes",nbytes);
    error_one(str);
  }
  return ptr;
}

/* ---------------------------------------------------------------------- */

void *CSlib::srealloc(void *ptr, int nbytes)
{
  if (nbytes == 0) {
    sfree(ptr);
    return NULL;
  }
  
  ptr = realloc(ptr,nbytes);
  if (ptr == NULL) {
    char str[128];
    sprintf(str,"realloc(): Failed to reallocate %d bytes",nbytes);
    error_one(str);
  }
  return ptr;
}

/* ---------------------------------------------------------------------- */

void CSlib::sfree(void *ptr)
{
  if (ptr == NULL) return;
  free(ptr);
}

/* ---------------------------------------------------------------------- */

void CSlib::error_all(const char *str)
{
  if (me == 0) printf("CSlib ERROR: %s\n",str);
  MPI_Comm world = (MPI_Comm) myworld;
  MPI_Abort(world,1);
}

/* ---------------------------------------------------------------------- */

void CSlib::error_one(const char *str)
{
  printf("CSlib ERROR: %s\n",str);
  MPI_Comm world = (MPI_Comm) myworld;
  MPI_Abort(world,1);
}

/* ----------------------------------------------------------------------
   round N up to multiple of nalign and return it
   NOTE: see mapreduce/src/keyvalue.cpp for doing this as uint64_t
------------------------------------------------------------------------- */

int64_t CSlib::roundup(int64_t n, int nalign)
{
  if (n % nalign == 0) return n;
  n = (n/nalign + 1) * nalign;
  return n;
}
