#include "xdr_compat.h"

#include <cstdlib>
#include <cstring>

/*
 * This file contains an implementation of the Sun External Data Representation (XDR)
 * routines.  They have been adapted specifically for the use with the LAMMPS xtc dump
 * style to produce compressed trajectory files in the Gromacs XTC format.
 *
 * The XDR sources are avaiable under the BSD 3-clause license for example in
 * the MIT Kerberos 5 distribution with the following copyright notice and license.
 *
 * @(#)xdr.h    2.2 88/07/29 4.0 RPCSRC
 *
 * Copyright (c) 2010, Oracle America, Inc.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 *     * Neither the name of the "Oracle America, Inc." nor the names of
 *       its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef FALSE
#define FALSE (0)
#endif
#ifndef TRUE
#define TRUE (1)
#endif

#define BYTES_PER_XDR_UNIT (4)

/*
 * for unit alignment
 */
static char xdr_zero[BYTES_PER_XDR_UNIT] = {0, 0, 0, 0};

static xdr_uint32_t xdr_swapbytes(xdr_uint32_t x)
{
  xdr_uint32_t y;
  int i;
  char *px = (char *) &x;
  char *py = (char *) &y;

  for (i = 0; i < 4; i++) py[i] = px[3 - i];

  return y;
}

static xdr_uint32_t xdr_htonl(xdr_uint32_t x)
{
  short s = 0x0F00;
  if (*((char *) &s) == (char) 0x0F) {
    /* bigendian, do nothing */
    return x;
  } else {
    /* smallendian,swap bytes */
    return xdr_swapbytes(x);
  }
}

static xdr_uint32_t xdr_ntohl(xdr_uint32_t x)
{
  short s = 0x0F00;
  if (*((char *) &s) == (char) 0x0F) {
    /* bigendian, do nothing */
    return x;
  } else {
    /* smallendian, swap bytes */
    return xdr_swapbytes(x);
  }
}

/*
 * XDR integers
 */
bool_t xdr_int(XDR *xdrs, int *ip)
{
  xdr_int32_t l;

  switch (xdrs->x_op) {

    case XDR_ENCODE:
      l = (xdr_int32_t) (*ip);
      return xdr_putint32(xdrs, &l);
      break;

    case XDR_DECODE:
      if (!xdr_getint32(xdrs, &l)) return FALSE;
      *ip = (int) l;
      return TRUE;
      break;

    case XDR_FREE:
      return TRUE;
      break;
  }
  return FALSE;
}

/*
 * XDR opaque data
 * Allows the specification of a fixed size sequence of opaque bytes.
 * cp points to the opaque object and cnt gives the byte length.
 */
bool_t xdr_opaque(XDR *xdrs, char *cp, unsigned int cnt)
{
  unsigned int rndup;
  static char crud[BYTES_PER_XDR_UNIT];

  /*
   * if no data we are done
   */
  if (cnt == 0) return TRUE;

  /*
   * round byte count to full xdr units
   */
  rndup = cnt % BYTES_PER_XDR_UNIT;
  if (rndup > 0) rndup = BYTES_PER_XDR_UNIT - rndup;

  switch (xdrs->x_op) {

    case XDR_DECODE:
      if (!xdr_getbytes(xdrs, cp, cnt)) { return FALSE; }
      if (rndup == 0) return TRUE;
      return xdr_getbytes(xdrs, (char *) crud, rndup);
      break;

    case XDR_ENCODE:
      if (!xdr_putbytes(xdrs, cp, cnt)) { return FALSE; }
      if (rndup == 0) return TRUE;
      return xdr_putbytes(xdrs, xdr_zero, rndup);
      break;

    case XDR_FREE:
      return TRUE;
      break;
  }
  return FALSE;
}

/* Floating-point stuff */

bool_t xdr_float(XDR *xdrs, float *fp)
{
  xdr_int32_t tmp;

  switch (xdrs->x_op) {

    case XDR_ENCODE:
      tmp = *(xdr_int32_t *) fp;
      return (xdr_putint32(xdrs, &tmp));
      break;

    case XDR_DECODE:
      if (xdr_getint32(xdrs, &tmp)) {
        *(xdr_int32_t *) fp = tmp;
        return TRUE;
      }
      break;

    case XDR_FREE:
      return TRUE;
      break;
  }
  return FALSE;
}

/* Array routines */

/*
 * xdr_vector():
 *
 * XDR a fixed length array. Unlike variable-length arrays,
 * the storage of fixed length arrays is static and unfreeable.
 * > basep: base of the array
 * > size: size of the array
 * > elemsize: size of each element
 * > xdr_elem: routine to XDR each element
 */
bool_t xdr_vector(XDR *xdrs, char *basep, unsigned int nelem, unsigned int elemsize,
                  xdrproc_t xdr_elem)
{
#define LASTUNSIGNED ((unsigned int) 0 - 1)
  unsigned int i;
  char *elptr;

  elptr = basep;
  for (i = 0; i < nelem; i++) {
    if (!(*xdr_elem)(xdrs, elptr, LASTUNSIGNED)) { return FALSE; }
    elptr += elemsize;
  }
  return TRUE;
#undef LASTUNSIGNED
}

static bool_t xdrstdio_getbytes(XDR *, char *, unsigned int);
static bool_t xdrstdio_putbytes(XDR *, char *, unsigned int);
static void xdrstdio_destroy(XDR *);
static bool_t xdrstdio_getint32(XDR *, xdr_int32_t *);
static bool_t xdrstdio_putint32(XDR *, xdr_int32_t *);

/*
 * Ops vector for stdio type XDR
 */
static const struct xdr_ops xdrstdio_ops = {
    xdrstdio_getbytes, /* deserialize counted bytes */
    xdrstdio_putbytes, /* serialize counted bytes */
    xdrstdio_destroy,  /* destroy stream */
    xdrstdio_getint32, /* deserialize a int */
    xdrstdio_putint32, /* serialize a int */
};

/*
 * Initialize a stdio xdr stream.
 * Sets the xdr stream handle xdrs for use on the stream file.
 * Operation flag is set to op.
 */
void xdrstdio_create(XDR *xdrs, FILE *file, enum xdr_op op)
{
  xdrs->x_op = op;
  /* We have to add the const since the `struct xdr_ops' in `struct XDR'
     is not `const'.  */
  xdrs->x_ops = (struct xdr_ops *) &xdrstdio_ops;
  xdrs->x_private = (char *) file;
  xdrs->x_handy = 0;
  xdrs->x_base = nullptr;
}

/*
 * Destroy a stdio xdr stream.
 * Cleans up the xdr stream handle xdrs previously set up by xdrstdio_create.
 */
static void xdrstdio_destroy(XDR *xdrs)
{
  (void) fflush((FILE *) xdrs->x_private);
  /* xx should we close the file ?? */
}

static bool_t xdrstdio_getbytes(XDR *xdrs, char *addr, unsigned int len)
{
  if ((len != 0) && (fread(addr, (int) len, 1, (FILE *) xdrs->x_private) != 1)) return FALSE;
  return TRUE;
}

static bool_t xdrstdio_putbytes(XDR *xdrs, char *addr, unsigned int len)
{
  if ((len != 0) && (fwrite(addr, (int) len, 1, (FILE *) xdrs->x_private) != 1)) return FALSE;
  return TRUE;
}

static bool_t xdrstdio_getint32(XDR *xdrs, xdr_int32_t *ip)
{
  xdr_int32_t mycopy;

  if (fread((char *) &mycopy, 4, 1, (FILE *) xdrs->x_private) != 1) return FALSE;
  *ip = xdr_ntohl(mycopy);
  return TRUE;
}

static bool_t xdrstdio_putint32(XDR *xdrs, xdr_int32_t *ip)
{
  xdr_int32_t mycopy = xdr_htonl(*ip);

  ip = &mycopy;
  if (fwrite((char *) ip, 4, 1, (FILE *) xdrs->x_private) != 1) return FALSE;
  return TRUE;
}

#ifdef __cplusplus
}
#endif
