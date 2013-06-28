#ifndef LMP_XDR_COMPAT_H
#define LMP_XDR_COMPAT_H

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * This file is needed for systems, that do not provide XDR support
 * in their system libraries. It was written for windows, but will
 * most probably work on other platforms too. better make sure you
 * test that the xtc files produced are ok before using it.
 *
 * It is also needed on BG/L, BG/P and Cray XT3/XT4/XT5 as we don't
 * have XDR support in the lightweight kernel runtimes either.
 *
 * This file contains the definitions for Sun External Data
 * Representation (XDR) headers and routines.
 *
 * Although the rest of LAMPPS is GPL, you can copy and use the XDR
 * routines in any way you want as long as you obey Sun's license:
 * Sun RPC is a product of Sun Microsystems, Inc. and is provided for
 * unrestricted use provided that this legend is included on all tape
 * media and as a part of the software program in whole or part.  Users
 * may copy or modify Sun RPC without charge, but are not authorized
 * to license or distribute it to anyone else except as part of a product or
 * program developed by the user.
 *
 * SUN RPC IS PROVIDED AS IS WITH NO WARRANTIES OF ANY KIND INCLUDING THE
 * WARRANTIES OF DESIGN, MERCHANTIBILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE, OR ARISING FROM A COURSE OF DEALING, USAGE OR TRADE PRACTICE.
 *
 * Sun RPC is provided with no support and without any obligation on the
 * part of Sun Microsystems, Inc. to assist in its use, correction,
 * modification or enhancement.
 *
 * SUN MICROSYSTEMS, INC. SHALL HAVE NO LIABILITY WITH RESPECT TO THE
 * INFRINGEMENT OF COPYRIGHTS, TRADE SECRETS OR ANY PATENTS BY SUN RPC
 * OR ANY PART THEREOF.
 *
 * In no event will Sun Microsystems, Inc. be liable for any lost revenue
 * or profits or other special, indirect and consequential damages, even if
 * Sun has been advised of the possibility of such damages.
 *
 * Sun Microsystems, Inc.
 * 2550 Garcia Avenue
 * Mountain View, California  94043
 */

/*
 * Xdr operations.  XDR_ENCODE causes the type to be encoded into the
 * stream.  XDR_DECODE causes the type to be extracted from the stream.
 * XDR_FREE can be used to release the space allocated by an
 * XDR_DECODE request.
 */

typedef int bool_t;

#if defined(__MINGW32_VERSION)
typedef char * caddr_t;
typedef unsigned int u_int;
#endif

/*
 * Aninteger type that is 32 bits wide. Check if int,
 * long or short is 32 bits and die if none of them is :-)
 */
#if (INT_MAX == 2147483647)
   typedef int xdr_int32_t;
   typedef unsigned int xdr_uint32_t;
#elif (LONG_MAX == 2147483647L)
   typedef long xdr_int32_t;
   typedef unsigned long xdr_uint32_t;
#elif (SHRT_MAX == 2147483647)
   typedef short xdr_int32_t;
   typedef unsigned short xdr_uint32_t;
#else
#  error ERROR: No 32 bit wide integer type found!
#endif

enum xdr_op {
  XDR_ENCODE = 0,
  XDR_DECODE = 1,
  XDR_FREE = 2
};

#ifndef FALSE
#      define  FALSE   (0)
#endif
#ifndef TRUE
#      define  TRUE    (1)
#endif

#define BYTES_PER_XDR_UNIT        (4)
/* Macro to round up to units of 4. */
#define XDR_RNDUP(x)  (((x) + BYTES_PER_XDR_UNIT - 1) & ~(BYTES_PER_XDR_UNIT - 1))


/*
 * The XDR handle.
 * Contains operation which is being applied to the stream,
 * an operations vector for the particular implementation (e.g. see xdr_mem.c),
 * and two private fields for the use of the particular implementation.
 */

typedef struct XDR XDR;
struct XDR
  {
    enum xdr_op x_op;                /* operation; fast additional param */
    struct xdr_ops  *x_ops;
    char *x_public;                /* users' data */
    char *x_private;                /* pointer to private data */
    char *x_base;                /* private used for position info */
    int x_handy;                /* extra private word */
  };

struct xdr_ops
      {
        bool_t (*x_getbytes) (XDR *__xdrs, char *__addr, unsigned int __len);
        /* get some bytes from " */
        bool_t (*x_putbytes) (XDR *__xdrs, char *__addr, unsigned int __len);
        /* put some bytes to " */
        unsigned int (*x_getpostn) (XDR *__xdrs);
        /* returns bytes off from beginning */
        bool_t (*x_setpostn) (XDR *__xdrs, unsigned int __pos);
        /* lets you reposition the stream */
        xdr_int32_t *(*x_inline) (XDR *__xdrs, int __len);
        /* buf quick ptr to buffered data */
        void (*x_destroy) (XDR *__xdrs);
        /* free privates of this xdr_stream */
        bool_t (*x_getint32) (XDR *__xdrs, xdr_int32_t *__ip);
        /* get a int from underlying stream */
        bool_t (*x_putint32) (XDR *__xdrs, xdr_int32_t *__ip);
        /* put a int to " */
        bool_t (*x_getuint32) (XDR *__xdrs, xdr_uint32_t *__ip);
        /* get a unsigned int from underlying stream */
        bool_t (*x_putuint32) (XDR *__xdrs, xdr_uint32_t *__ip);
        /* put a int to " */
};

/*
 * A xdrproc_t exists for each data type which is to be encoded or decoded.
 *
 * The second argument to the xdrproc_t is a pointer to an opaque pointer.
 * The opaque pointer generally points to a structure of the data type
 * to be decoded.  If this pointer is 0, then the type routines should
 * allocate dynamic storage of the appropriate size and return it.
 */

typedef bool_t (*xdrproc_t) (XDR *, void *,...);

/*
 * Operations defined on a XDR handle
 *
 * XDR          *xdrs;
 * xdr_int32_t  *int32p;
 * long         *longp;
 * char         *addr;
 * unsigned int  len;
 * unsigned int  pos;
 */


#define xdr_getint32(xdrs, int32p)                      \
        (*(xdrs)->x_ops->x_getint32)(xdrs, int32p)

#define xdr_putint32(xdrs, int32p)                      \
        (*(xdrs)->x_ops->x_putint32)(xdrs, int32p)

#define xdr_getuint32(xdrs, uint32p)                      \
        (*(xdrs)->x_ops->x_getuint32)(xdrs, uint32p)

#define xdr_putuint32(xdrs, uint32p)                      \
        (*(xdrs)->x_ops->x_putuint32)(xdrs, uint32p)

#define xdr_getbytes(xdrs, addr, len)                        \
        (*(xdrs)->x_ops->x_getbytes)(xdrs, addr, len)

#define xdr_putbytes(xdrs, addr, len)                        \
        (*(xdrs)->x_ops->x_putbytes)(xdrs, addr, len)

#define xdr_getpos(xdrs)                                \
        (*(xdrs)->x_ops->x_getpostn)(xdrs)

#define xdr_setpos(xdrs, pos)                                \
        (*(xdrs)->x_ops->x_setpostn)(xdrs, pos)

#define        xdr_inline(xdrs, len)                                \
        (*(xdrs)->x_ops->x_inline)(xdrs, len)

#define        xdr_destroy(xdrs)                                        \
        do {                                                        \
                if ((xdrs)->x_ops->x_destroy)                        \
                        (*(xdrs)->x_ops->x_destroy)(xdrs);        \
        } while (0)


extern bool_t xdr_int (XDR *__xdrs, int *__ip);
extern bool_t xdr_u_int (XDR *__xdrs, unsigned int *__ip);
extern bool_t xdr_short (XDR *__xdrs, short *__ip);
extern bool_t xdr_u_short (XDR *__xdrs, unsigned short *__ip);
extern bool_t xdr_bool (XDR *__xdrs, int *__bp);
extern bool_t xdr_opaque (XDR *__xdrs, char *__cp, unsigned int __cnt);
extern bool_t xdr_string (XDR *__xdrs, char **__cpp, unsigned int __maxsize);
extern bool_t xdr_char (XDR *__xdrs, char *__cp);
extern bool_t xdr_u_char (XDR *__xdrs, unsigned char *__cp);
extern bool_t xdr_vector (XDR *__xdrs, char *__basep, unsigned int __nelem,
                          unsigned int __elemsize, xdrproc_t __xdr_elem);
extern bool_t xdr_float (XDR *__xdrs, float *__fp);
extern bool_t xdr_double (XDR *__xdrs, double *__dp);
extern void xdrstdio_create (XDR *__xdrs, FILE *__file, enum xdr_op __xop);

/* free memory buffers for xdr */
extern void xdr_free (xdrproc_t __proc, char *__objp);

#ifdef __cplusplus
}
#endif


#endif /* XDR_COMPAT_H */
