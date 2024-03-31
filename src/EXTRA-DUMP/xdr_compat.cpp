// clang-format off
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
 * @(#)xdr.h	2.2 88/07/29 4.0 RPCSRC
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

/* compatibility typedefs */
#if defined(_WIN32) || defined(__APPLE__) || defined(__FreeBSD__) || defined(__DragonFly__) || \
  defined(__OpenBSD__) || defined(__NetBSD__) || (defined(__linux__) && !defined(__GLIBC_MINOR__))
typedef char *caddr_t;
typedef unsigned int u_int;
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
  char *px=(char *)&x;
  char *py=(char *)&y;

  for (i=0;i<4;i++)
    py[i]=px[3-i];

  return y;
}

static xdr_uint32_t xdr_htonl(xdr_uint32_t x)
{
  short s=0x0F00;
  if (*((char *)&s)==(char)0x0F) {
    /* bigendian, do nothing */
    return x;
  } else {
    /* smallendian,swap bytes */
    return xdr_swapbytes(x);
  }
}

static xdr_uint32_t xdr_ntohl(xdr_uint32_t x)
{
  short s=0x0F00;
  if (*((char *)&s)==(char)0x0F) {
    /* bigendian, do nothing */
    return x;
  } else {
    /* smallendian, swap bytes */
    return xdr_swapbytes(x);
  }
}


/*
 * Free a data structure using XDR
 * Not a filter, but a convenient utility nonetheless
 */
void xdr_free(xdrproc_t proc, char *objp)
{
  XDR x;

  x.x_op = XDR_FREE;
  (*proc)(&x, objp);
}

/*
 * XDR nothing
 */
bool_t xdr_void(void)
{
  return TRUE;
}

/*
 * XDR integers
 */
bool_t xdr_int(XDR *xdrs, int *ip)
{
  xdr_int32_t l;

  switch (xdrs->x_op)
  {
    case XDR_ENCODE:
      l = (xdr_int32_t)(*ip);
      return xdr_putint32(xdrs, &l);

    case XDR_DECODE:
      if (!xdr_getint32(xdrs, &l))
            return FALSE;
      *ip = (int) l;
      return TRUE;

    case XDR_FREE:
      return TRUE;
  }
  return FALSE;
}


/*
 * XDR unsigned integers
 */
bool_t xdr_u_int(XDR *xdrs, unsigned int *up)
{
  xdr_uint32_t l;

  switch (xdrs->x_op)
  {
    case XDR_ENCODE:
      l = (xdr_uint32_t)(*up);
      return xdr_putuint32(xdrs, &l);

    case XDR_DECODE:
      if (!xdr_getuint32(xdrs, &l))
            return FALSE;
      *up = (unsigned int) l;
      return TRUE;

    case XDR_FREE:
      return TRUE;
  }
  return FALSE;
}




/*
 * XDR short integers
 */
bool_t
xdr_short(XDR *xdrs, short *sp)
{
  xdr_int32_t l;

  switch (xdrs->x_op)
    {
    case XDR_ENCODE:
      l = (xdr_int32_t) *sp;
      return xdr_putint32(xdrs, &l);

    case XDR_DECODE:
      if (!xdr_getint32(xdrs, &l))
        {
          return FALSE;
        }
      *sp = (short) l;
      return TRUE;

    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
}


/*
 * XDR unsigned short integers
 */
bool_t
xdr_u_short(XDR *xdrs, unsigned short *usp)
{
  xdr_uint32_t l;

  switch (xdrs->x_op)
    {
    case XDR_ENCODE:
      l = (xdr_uint32_t) *usp;
      return xdr_putuint32(xdrs, &l);

    case XDR_DECODE:
      if (!xdr_getuint32(xdrs, &l))
        {
          return FALSE;
        }
          *usp = (unsigned short) l;
      return TRUE;

    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
}


/*
 * XDR a char
 */
bool_t
xdr_char(XDR *xdrs, char *cp)
{
  int i;

  i = (*cp);
  if (!xdr_int(xdrs, &i))
    {
      return FALSE;
    }
  *cp = i;
  return TRUE;
}

/*
 * XDR an unsigned char
 */
bool_t
xdr_u_char(XDR *xdrs, unsigned char *cp)
{
  unsigned int u;

  u = (*cp);
  if (!xdr_u_int(xdrs, &u))
    {
      return FALSE;
    }
  *cp = u;
  return TRUE;
}

/*
 * XDR booleans
 */
bool_t
xdr_bool(XDR *xdrs, int *bp)
{
#define XDR_FALSE        ((xdr_int32_t) 0)
#define XDR_TRUE        ((xdr_int32_t) 1)

  xdr_int32_t lb;

  switch (xdrs->x_op)
    {
    case XDR_ENCODE:
      lb = *bp ? XDR_TRUE : XDR_FALSE;
      return xdr_putint32(xdrs, &lb);

    case XDR_DECODE:
      if (!xdr_getint32(xdrs, &lb))
        {
          return FALSE;
        }
      *bp = (lb == XDR_FALSE) ? FALSE : TRUE;
      return TRUE;

    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
#undef XDR_FALSE
#undef XDR_TRUE
}



/*
 * XDR opaque data
 * Allows the specification of a fixed size sequence of opaque bytes.
 * cp points to the opaque object and cnt gives the byte length.
 */
bool_t
xdr_opaque(XDR *xdrs, char *cp, unsigned int cnt)
{
  unsigned int rndup;
  static char crud[BYTES_PER_XDR_UNIT];

  /*
   * if no data we are done
   */
  if (cnt == 0)
    return TRUE;

  /*
   * round byte count to full xdr units
   */
  rndup = cnt % BYTES_PER_XDR_UNIT;
  if (rndup > 0)
    rndup = BYTES_PER_XDR_UNIT - rndup;

  switch (xdrs->x_op)
    {
    case XDR_DECODE:
      if (!xdr_getbytes(xdrs, cp, cnt))
        {
          return FALSE;
        }
      if (rndup == 0)
        return TRUE;
      return xdr_getbytes(xdrs, (char *)crud, rndup);

    case XDR_ENCODE:
      if (!xdr_putbytes(xdrs, cp, cnt))
        {
          return FALSE;
        }
      if (rndup == 0)
        return TRUE;
      return xdr_putbytes(xdrs, xdr_zero, rndup);

    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
}


/*
 * XDR null terminated ASCII strings
 * xdr_string deals with "C strings" - arrays of bytes that are
 * terminated by a nullptr character.  The parameter cpp references a
 * pointer to storage; If the pointer is null, then the necessary
 * storage is allocated.  The last parameter is the max allowed length
 * of the string as specified by a protocol.
 */
bool_t
xdr_string(XDR *xdrs, char **cpp, unsigned int maxsize)
{
  char *sp = *cpp;        /* sp is the actual string pointer */
  unsigned int size = 0;
  unsigned int nodesize = 0;

  /*
   * first deal with the length since xdr strings are counted-strings
   */
  switch (xdrs->x_op)
    {
    case XDR_FREE:
      if (sp == nullptr)
        {
          return TRUE;                /* already free */
        }
      /* fall through... */
    case XDR_ENCODE:
      if (sp == nullptr)
            return FALSE;
      size = strlen(sp);
      break;
    case XDR_DECODE:
      break;
    }

  if (!xdr_u_int(xdrs, &size))
    {
      return FALSE;
    }
  if (size > maxsize)
    {
      return FALSE;
    }
  nodesize = size + 1;

  /*
   * now deal with the actual bytes
   */
  switch (xdrs->x_op)
    {
    case XDR_DECODE:
      if (nodesize == 0)
        {
          return TRUE;
        }
      if (sp == nullptr)
        *cpp = sp = (char *) malloc(nodesize);
      if (sp == nullptr)
        {
          (void) fputs("xdr_string: out of memory\n", stderr);
          return FALSE;
        }
      sp[size] = 0;
      return xdr_opaque(xdrs, sp, size);

    case XDR_ENCODE:
      return xdr_opaque(xdrs, sp, size);

    case XDR_FREE:
      free(sp);
      *cpp = nullptr;
      return TRUE;
    }
  return FALSE;
}



/* Floating-point stuff */

bool_t
xdr_float(XDR *xdrs, float *fp)
{
        xdr_int32_t tmp;

        switch (xdrs->x_op) {

        case XDR_ENCODE:
                tmp = *(xdr_int32_t *)fp;
               return (xdr_putint32(xdrs, &tmp));

                break;

        case XDR_DECODE:
                        if (xdr_getint32(xdrs, &tmp)) {
                                *(xdr_int32_t *)fp = tmp;
                                return TRUE;
                        }

                break;

        case XDR_FREE:
                return TRUE;
        }
        return FALSE;
}


bool_t
xdr_double(XDR *xdrs, double *dp)
{

  /* Windows and some other systems dont define double-precision
   * word order in the header files, so unfortunately we have
   * to calculate it!
   */
  static int LSW=-1; /* Least significant fp word */
  int *ip;
  xdr_int32_t tmp[2];

  if (LSW<0) {
    double x=0.987654321; /* Just a number */

    /* Possible representations in IEEE double precision:
     * (S=small endian, B=big endian)
     *
     * Byte order, Word order, Hex
     *     S           S       b8 56 0e 3c dd 9a ef 3f
     *     B           S       3c 0e 56 b8 3f ef 9a dd
     *     S           B       dd 9a ef 3f b8 56 0e 3c
     *     B           B       3f ef 9a dd 3c 0e 56 b8
     */

    unsigned char ix = *((char *)&x);

    if (ix==0xdd || ix==0x3f)
      LSW=1;  /* Big endian word order */
    else if (ix==0xb8 || ix==0x3c)
      LSW=0;  /* Small endian word order */
    else { /* Catch strange errors */
      printf("Error when detecting floating-point word order.\n"
             "Do you have a non-IEEE system?\n"
             "If possible, use the XDR libraries provided with your system,\n"
             "instead of the Gromacs fallback XDR source.\n");
      exit(0);
    }
  }

  switch (xdrs->x_op) {

  case XDR_ENCODE:
    ip = (int *)dp;
    tmp[0] = ip[!LSW];
    tmp[1] = ip[LSW];
    return (xdr_putint32(xdrs, tmp) &&
               xdr_putint32(xdrs, tmp+1));

    break;

  case XDR_DECODE:
    ip = (int *)dp;
    if (xdr_getint32(xdrs, tmp+!LSW) &&
           xdr_getint32(xdrs, tmp+LSW)) {
        ip[0] = tmp[0];
        ip[1] = tmp[1];
        return TRUE;
    }

    break;

  case XDR_FREE:
    return TRUE;
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
bool_t
xdr_vector(XDR *xdrs, char *basep, unsigned int nelem,
            unsigned int elemsize, xdrproc_t xdr_elem)
{
#define LASTUNSIGNED        ((unsigned int)0-1)
  unsigned int i;
  char *elptr;

  elptr = basep;
  for (i = 0; i < nelem; i++)
    {
      if (!(*xdr_elem)(xdrs, elptr, LASTUNSIGNED))
        {
          return FALSE;
        }
      elptr += elemsize;
    }
  return TRUE;
#undef LASTUNSIGNED
}

static bool_t xdrstdio_getbytes(XDR *, char *, unsigned int);
static bool_t xdrstdio_putbytes(XDR *, char *, unsigned int);
static unsigned int xdrstdio_getpos(XDR *);
static bool_t xdrstdio_setpos(XDR *, unsigned int);
static xdr_int32_t *xdrstdio_inline(XDR *, int);
static void xdrstdio_destroy(XDR *);
static bool_t xdrstdio_getint32(XDR *, xdr_int32_t *);
static bool_t xdrstdio_putint32(XDR *, xdr_int32_t *);
static bool_t xdrstdio_getuint32(XDR *, xdr_uint32_t *);
static bool_t xdrstdio_putuint32(XDR *, xdr_uint32_t *);

/*
 * Ops vector for stdio type XDR
 */
static const struct xdr_ops xdrstdio_ops =
{
  xdrstdio_getbytes,               /* deserialize counted bytes */
  xdrstdio_putbytes,             /* serialize counted bytes */
  xdrstdio_getpos,                /* get offset in the stream */
  xdrstdio_setpos,                /* set offset in the stream */
  xdrstdio_inline,                /* prime stream for inline macros */
  xdrstdio_destroy,                /* destroy stream */
  xdrstdio_getint32,        /* deserialize a int */
  xdrstdio_putint32,        /* serialize a int */
  xdrstdio_getuint32,        /* deserialize a int */
  xdrstdio_putuint32                /* serialize a int */
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
  if ((len != 0) && (fread(addr, (int) len, 1, (FILE *) xdrs->x_private) != 1))
    return FALSE;
  return TRUE;
}

static bool_t xdrstdio_putbytes(XDR *xdrs, char *addr, unsigned int len)
{
  if ((len != 0) && (fwrite(addr, (int) len, 1, (FILE *) xdrs->x_private) != 1))
    return FALSE;
  return TRUE;
}

static unsigned int xdrstdio_getpos(XDR *xdrs)
{
  return (unsigned int) ftell((FILE *) xdrs->x_private);
}

static bool_t xdrstdio_setpos(XDR *xdrs, unsigned int pos)
{
  return fseek((FILE *) xdrs->x_private, (xdr_int32_t) pos, 0) < 0 ? FALSE : TRUE;
}

static xdr_int32_t *xdrstdio_inline(XDR * /*xdrs*/, int /*len*/)
{
  /*
   * Must do some work to implement this: must ensure
   * enough data in the underlying stdio buffer,
   * that the buffer is aligned so that we can indirect through a
   * long *, and stuff this pointer in xdrs->x_buf.  Doing
   * a fread or fwrite to a scratch buffer would defeat
   * most of the gains to be had here and require storage
   * management on this buffer, so we don't do this.
   */
  return nullptr;
}

static bool_t xdrstdio_getint32(XDR *xdrs, xdr_int32_t *ip)
{
  xdr_int32_t mycopy;

  if (fread((char *) &mycopy, 4, 1, (FILE *) xdrs->x_private) != 1)
    return FALSE;
  *ip = xdr_ntohl(mycopy);
  return TRUE;
}

static bool_t xdrstdio_putint32(XDR *xdrs, xdr_int32_t *ip)
{
  xdr_int32_t mycopy = xdr_htonl(*ip);

  ip = &mycopy;
  if (fwrite((char *) ip, 4, 1, (FILE *) xdrs->x_private) != 1)
    return FALSE;
  return TRUE;
}

static bool_t xdrstdio_getuint32(XDR *xdrs, xdr_uint32_t *ip)
{
        xdr_uint32_t mycopy;

        if (fread((char *) &mycopy, 4, 1, (FILE *) xdrs->x_private) != 1)
                return FALSE;
        *ip = xdr_ntohl (mycopy);
        return TRUE;
}

static bool_t xdrstdio_putuint32(XDR *xdrs, xdr_uint32_t *ip)
{
        xdr_uint32_t mycopy = xdr_htonl (*ip);

        ip = &mycopy;
        if (fwrite((char *) ip, 4, 1, (FILE *) xdrs->x_private) != 1)
                return FALSE;
        return TRUE;
}

#ifdef __cplusplus
}
#endif
