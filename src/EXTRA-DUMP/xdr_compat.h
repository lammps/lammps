#ifndef LMP_XDR_COMPAT_H
#define LMP_XDR_COMPAT_H

#include <cstdint>
#include <cstdio>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * This file contains the definitions for Sun External Data Representation (XDR).
 * They have been adapted specifically for the use with the LAMMPS xtc dump style
 * to produce compressed trajectory files in the Gromacs XTC format.
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

/* compatibility typedefs */
typedef int bool_t;

typedef int32_t xdr_int32_t;
typedef uint32_t xdr_uint32_t;

/*
 * Xdr operations.  XDR_ENCODE causes the type to be encoded into the
 * stream.  XDR_DECODE causes the type to be extracted from the stream.
 * XDR_FREE can be used to release the space allocated by an
 * XDR_DECODE request.
 */

enum xdr_op { XDR_ENCODE = 0, XDR_DECODE = 1, XDR_FREE = 2 };

/*
 * The XDR handle.
 * Contains operation which is being applied to the stream,
 * an operations vector for the particular implementation (e.g. see xdr_mem.c),
 * and two private fields for the use of the particular implementation.
 */

typedef struct XDR XDR;
struct XDR {
  enum xdr_op x_op; /* operation; fast additional param */
  struct xdr_ops *x_ops;
  char *x_public;  /* users' data */
  char *x_private; /* pointer to private data */
  char *x_base;    /* private used for position info */
  int x_handy;     /* extra private word */
};

struct xdr_ops {
  /* get some bytes from XDR stream */
  bool_t (*x_getbytes)(XDR *__xdrs, char *__addr, unsigned int __len);
  /* put some bytes to XDR stream */
  bool_t (*x_putbytes)(XDR *__xdrs, char *__addr, unsigned int __len);
  /* free privates of this xdr_stream */
  void (*x_destroy)(XDR *__xdrs);
  /* get a int from XDR stream */
  bool_t (*x_getint32)(XDR *__xdrs, xdr_int32_t *__ip);
  /* put a int to XDR stream */
  bool_t (*x_putint32)(XDR *__xdrs, xdr_int32_t *__ip);
};

/*
 * A xdrproc_t exists for each data type which is to be encoded or decoded.
 *
 * The second argument to the xdrproc_t is a pointer to an opaque pointer.
 * The opaque pointer generally points to a structure of the data type
 * to be decoded.  If this pointer is 0, then the type routines should
 * allocate dynamic storage of the appropriate size and return it.
 */

typedef bool_t (*xdrproc_t)(XDR *, void *, ...);

/*
 * Operations defined on a XDR handle
 *
 * XDR          *xdrs;
 * xdr_int32_t  *int32p;
 * unsigned int  len;
 */

#define xdr_getint32(xdrs, int32p) (*(xdrs)->x_ops->x_getint32)(xdrs, int32p)
#define xdr_putint32(xdrs, int32p) (*(xdrs)->x_ops->x_putint32)(xdrs, int32p)
#define xdr_getbytes(xdrs, addr, len) (*(xdrs)->x_ops->x_getbytes)(xdrs, addr, len)
#define xdr_putbytes(xdrs, addr, len) (*(xdrs)->x_ops->x_putbytes)(xdrs, addr, len)
#define xdr_destroy(xdrs)                                            \
  do {                                                               \
    if ((xdrs)->x_ops->x_destroy) (*(xdrs)->x_ops->x_destroy)(xdrs); \
  } while (0)

extern bool_t xdr_int(XDR *__xdrs, int *__ip);
extern bool_t xdr_opaque(XDR *__xdrs, char *__cp, unsigned int __cnt);
extern bool_t xdr_vector(XDR *__xdrs, char *__basep, unsigned int __nelem, unsigned int __elemsize,
                         xdrproc_t __xdr_elem);
extern bool_t xdr_float(XDR *__xdrs, float *__fp);
extern void xdrstdio_create(XDR *__xdrs, FILE *__file, enum xdr_op __xop);

#ifdef __cplusplus
}
#endif

#endif /* XDR_COMPAT_H */
