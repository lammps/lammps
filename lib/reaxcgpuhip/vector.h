/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, haktulga@cs.purdue.edu
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#ifndef __VECTOR_H_
#define __VECTOR_H_

#include "reaxc_types.h"

#include <assert.h>

#include "random.h"

#ifdef __cplusplus
extern "C"  {
#endif

#if defined(LAMMPS_REAX) || defined(PURE_REAX)
CUDA_HOST_DEVICE static inline int Vector_isZero( const real * const v, int k )
{
    assert( k >= 0 );

    for ( --k; k >= 0; --k )
    {
        if ( FABS( v[k] ) > ALMOST_ZERO )
        {
            return FALSE;
        }
    }

    return TRUE;
}


CUDA_HOST_DEVICE static inline void Vector_MakeZero( real * const v, int k )
{
    assert( k >= 0 );

    for ( --k; k >= 0; --k )
    {
        v[k] = 0;
    }
}


CUDA_HOST_DEVICE static inline void Vector_Copy( real * const dest, const real * const v, int k )
{
    assert( k >= 0 );

    for ( --k; k >= 0; --k )
    {
        dest[k] = v[k];
    }
}


CUDA_HOST_DEVICE static inline void Vector_Scale( real * const dest, real c,
        const real * const v, int k )
{
    assert( k >= 0 );

    for ( --k; k >= 0; --k )
    {
        dest[k] = c * v[k];
    }
}


CUDA_HOST_DEVICE static inline void Vector_Sum( real * const dest, real c,
        const real * const v, real d, const real * const y, int k )
{
    assert( k >= 0 );

    for ( --k; k >= 0; --k )
    {
        dest[k] = c * v[k] + d * y[k];
    }
}


CUDA_HOST_DEVICE static inline void Vector_Add( real * const dest, real c,
        const real * const v, int k )
{
    assert( k >= 0 );

    for ( --k; k >= 0; --k )
    {
        dest[k] += c * v[k];
    }
}


CUDA_HOST_DEVICE static inline real Dot( const real * const v1,
        const real * const v2, int k )
{
    real ret;

    assert( k >= 0 );
    ret = 0.0;

    for ( --k; k >= 0; --k )
    {
        ret +=  v1[k] * v2[k];
    }

    return ret;
}


CUDA_HOST_DEVICE static inline real Norm( const real * const v1, int k )
{
    real ret;

    assert( k >= 0 );
    ret = 0.0;

    for ( --k; k >= 0; --k )
    {
        ret +=  SQR( v1[k] );
    }

    return SQRT( ret );
}


CUDA_HOST_DEVICE static inline void Vector_Print( FILE * const fout,
        const char * const vname, const real * const v, int k )
{
    int i;

    assert( k >= 0 );

    fprintf( fout, "%s:", vname );
    for ( i = 0; i < k; ++i )
    {
        fprintf( fout, "%8.3f\n", v[i] );
    }
    fprintf( fout, "\n" );
}


CUDA_HOST_DEVICE static inline void rvec_Copy( rvec dest, const rvec src )
{
    dest[0] = src[0];
    dest[1] = src[1];
    dest[2] = src[2];
}


CUDA_HOST_DEVICE static inline void rvec_Scale( rvec ret, real c, const rvec v )
{
    ret[0] = c * v[0];
    ret[1] = c * v[1];
    ret[2] = c * v[2];
}


CUDA_HOST_DEVICE static inline void rvec_Add( rvec ret, const rvec v )
{
    ret[0] += v[0];
    ret[1] += v[1];
    ret[2] += v[2];
}


CUDA_HOST_DEVICE static inline void rvec_ScaledAdd( rvec ret, real c, const rvec v )
{
    ret[0] += c * v[0];
    ret[1] += c * v[1];
    ret[2] += c * v[2];
}


CUDA_HOST_DEVICE static inline void rvec_Sum( rvec ret, const rvec v1, const rvec v2 )
{
    ret[0] = v1[0] + v2[0];
    ret[1] = v1[1] + v2[1];
    ret[2] = v1[2] + v2[2];
}


CUDA_HOST_DEVICE static inline void rvec_ScaledSum( rvec ret, real c1, const rvec v1,
        real c2, const rvec v2 )
{
    ret[0] = c1 * v1[0] + c2 * v2[0];
    ret[1] = c1 * v1[1] + c2 * v2[1];
    ret[2] = c1 * v1[2] + c2 * v2[2];
}


CUDA_HOST_DEVICE static inline real rvec_Dot( const rvec v1, const rvec v2 )
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}


CUDA_HOST_DEVICE static inline real rvec_ScaledDot( real c1, const rvec v1,
        real c2, const rvec v2 )
{
    return (c1 * c2) * (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
}


CUDA_HOST_DEVICE static inline void rvec_Multiply( rvec r, const rvec v1, const rvec v2 )
{
    r[0] = v1[0] * v2[0];
    r[1] = v1[1] * v2[1];
    r[2] = v1[2] * v2[2];
}


CUDA_HOST_DEVICE static inline void rvec_iMultiply( rvec r, const ivec v1, const rvec v2 )
{
    r[0] = v1[0] * v2[0];
    r[1] = v1[1] * v2[1];
    r[2] = v1[2] * v2[2];
}


CUDA_HOST_DEVICE static inline void rvec_Divide( rvec r, const rvec v1, const rvec v2 )
{
    r[0] = v1[0] / v2[0];
    r[1] = v1[1] / v2[1];
    r[2] = v1[2] / v2[2];
}


CUDA_HOST_DEVICE static inline void rvec_iDivide( rvec r, const rvec v1, const ivec v2 )
{
    r[0] = v1[0] / v2[0];
    r[1] = v1[1] / v2[1];
    r[2] = v1[2] / v2[2];
}


CUDA_HOST_DEVICE static inline void rvec_Invert( rvec r, const rvec v )
{
    r[0] = 1.0 / v[0];
    r[1] = 1.0 / v[1];
    r[2] = 1.0 / v[2];
}


CUDA_HOST_DEVICE static inline void rvec_Cross( rvec ret,
        const rvec v1, const rvec v2 )
{
    ret[0] = v1[1] * v2[2] - v1[2] * v2[1];
    ret[1] = v1[2] * v2[0] - v1[0] * v2[2];
    ret[2] = v1[0] * v2[1] - v1[1] * v2[0];
}


CUDA_HOST_DEVICE static inline void rvec_OuterProduct( rtensor r,
        const rvec v1, const rvec v2 )
{
    int i, j;

    for ( i = 0; i < 3; ++i )
    {
        for ( j = 0; j < 3; ++j )
        {
            r[i][j] = v1[i] * v2[j];
        }
    }
}


CUDA_HOST_DEVICE static inline real rvec_Norm_Sqr( rvec v )
{
    return SQR(v[0]) + SQR(v[1]) + SQR(v[2]);
}


CUDA_HOST_DEVICE static inline real rvec_Norm( rvec v )
{
    return SQRT( SQR(v[0]) + SQR(v[1]) + SQR(v[2]) );
}


CUDA_HOST_DEVICE static inline int rvec_isZero( rvec v )
{
    if ( FABS(v[0]) > ALMOST_ZERO ||
            FABS(v[1]) > ALMOST_ZERO ||
            FABS(v[2]) > ALMOST_ZERO )
    {
        return FALSE;
    }

    return TRUE;
}


CUDA_HOST_DEVICE static inline void rvec_MakeZero( rvec v )
{
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
}


static inline void rvec_Random( rvec v )
{
    v[0] = Random( 2.0 ) - 1.0;
    v[1] = Random( 2.0 ) - 1.0;
    v[2] = Random( 2.0 ) - 1.0;
}

CUDA_HOST_DEVICE static inline void rtensor_Multiply( rtensor ret,
        const rtensor m1, const rtensor m2 )
{
    int i, j, k;
    rtensor temp;

    // check if the result matrix is the same as one of m1, m2.
    // if so, we cannot modify the contents of m1 or m2, so
    // we have to use a temp matrix.
    if ( ret == m1 || ret == m2 )
    {
        for ( i = 0; i < 3; ++i )
        {
            for ( j = 0; j < 3; ++j )
            {
                temp[i][j] = 0;

                for ( k = 0; k < 3; ++k )
                {
                    temp[i][j] += m1[i][k] * m2[k][j];
                }
            }
        }

        for ( i = 0; i < 3; ++i )
        {
            for ( j = 0; j < 3; ++j )
            {
                ret[i][j] = temp[i][j];
            }
        }
    }
    else
    {
        for ( i = 0; i < 3; ++i )
        {
            for ( j = 0; j < 3; ++j )
            {
                ret[i][j] = 0;

                for ( k = 0; k < 3; ++k )
                {
                    ret[i][j] += m1[i][k] * m2[k][j];
                }
            }
        }
    }
}


CUDA_HOST_DEVICE static inline void rtensor_MatVec( rvec ret,
        const rtensor m, const rvec v )
{
    int i;
    rvec temp;

    // if ret is the same vector as v, we cannot modify the
    // contents of v until all computation is finished.
    if ( ret == v )
    {
        for ( i = 0; i < 3; ++i )
        {
            temp[i] = m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2];
        }

        for ( i = 0; i < 3; ++i )
        {
            ret[i] = temp[i];
        }
    }
    else
    {
        for ( i = 0; i < 3; ++i )
        {
            ret[i] = m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2];
        }
    }
}


CUDA_HOST_DEVICE static inline void rtensor_Scale( rtensor ret,
        real c, const rtensor m )
{
    int i, j;

    for ( i = 0; i < 3; ++i )
    {
        for ( j = 0; j < 3; ++j )
        {
            ret[i][j] = c * m[i][j];
        }
    }
}


CUDA_HOST_DEVICE static inline void rtensor_Add( rtensor ret,
        const rtensor t )
{
    int i, j;

    for ( i = 0; i < 3; ++i )
    {
        for ( j = 0; j < 3; ++j )
        {
            ret[i][j] += t[i][j];
        }
    }
}


CUDA_HOST_DEVICE static inline void rtensor_ScaledAdd( rtensor ret,
        real c, const rtensor t )
{
    int i, j;

    for ( i = 0; i < 3; ++i )
    {
        for ( j = 0; j < 3; ++j )
        {
            ret[i][j] += c * t[i][j];
        }
    }
}


CUDA_HOST_DEVICE static inline void rtensor_Sum( rtensor ret,
        const rtensor t1, const rtensor t2 )
{
    int i, j;

    for ( i = 0; i < 3; ++i )
    {
        for ( j = 0; j < 3; ++j )
        {
            ret[i][j] = t1[i][j] + t2[i][j];
        }
    }
}


CUDA_HOST_DEVICE static inline void rtensor_ScaledSum( rtensor ret,
        real c1, const rtensor t1, real c2, const rtensor t2 )
{
    int i, j;

    for ( i = 0; i < 3; ++i )
    {
        for ( j = 0; j < 3; ++j )
        {
            ret[i][j] = c1 * t1[i][j] + c2 * t2[i][j];
        }
    }
}


CUDA_HOST_DEVICE static inline void rtensor_Copy( rtensor ret,
        const rtensor t )
{
    int i, j;

    for ( i = 0; i < 3; ++i )
    {
        for ( j = 0; j < 3; ++j )
        {
            ret[i][j] = t[i][j];
        }
    }
}


CUDA_HOST_DEVICE static inline void rtensor_Identity( rtensor t )
{
    t[0][0] = 1.0;
    t[1][1] = 1.0;
    t[2][2] = 1.0;
    t[0][1] = 0.0;
    t[0][2] = 0.0;
    t[1][0] = 0.0;
    t[1][2] = 0.0;
    t[2][0] = 0.0;
    t[2][1] = 0.0;
}


CUDA_HOST_DEVICE static inline void rtensor_MakeZero( rtensor t )
{
    t[0][0] = 0.0;
    t[0][1] = 0.0;
    t[0][2] = 0.0;
    t[1][0] = 0.0;
    t[1][1] = 0.0;
    t[1][2] = 0.0;
    t[2][0] = 0.0;
    t[2][1] = 0.0;
    t[2][2] = 0.0;
}


CUDA_HOST_DEVICE static inline void rtensor_Transpose( rtensor ret,
        const rtensor t )
{
    ret[0][0] = t[0][0];
    ret[1][1] = t[1][1];
    ret[2][2] = t[2][2];
    ret[0][1] = t[1][0];
    ret[0][2] = t[2][0];
    ret[1][0] = t[0][1];
    ret[1][2] = t[2][1];
    ret[2][0] = t[0][2];
    ret[2][1] = t[1][2];
}


CUDA_HOST_DEVICE static inline real rtensor_Det( const rtensor t )
{
    return ( t[0][0] * (t[1][1] * t[2][2] - t[1][2] * t[2][1] ) +
            t[0][1] * (t[1][2] * t[2][0] - t[1][0] * t[2][2] ) +
            t[0][2] * (t[1][0] * t[2][1] - t[1][1] * t[2][0] ) );
}


CUDA_HOST_DEVICE static inline real rtensor_Trace( rtensor t )
{
    return (t[0][0] + t[1][1] + t[2][2]);
}


CUDA_HOST_DEVICE static inline void Print_rTensor( FILE * const fp,
        const rtensor t )
{
    int i, j;

    for ( i = 0; i < 3; i++ )
    {
        fprintf( fp, "[" );

        for ( j = 0; j < 3; j++ )
        {
            fprintf( fp, "%8.3f,\t", t[i][j] );
        }

        fprintf( fp, "]\n" );
    }
}


CUDA_HOST_DEVICE static inline void ivec_MakeZero( ivec v )
{
    v[0] = 0;
    v[1] = 0;
    v[2] = 0;
}


CUDA_HOST_DEVICE static inline void ivec_Copy( ivec dest, const ivec src )
{
    dest[0] = src[0];
    dest[1] = src[1];
    dest[2] = src[2];
}


CUDA_HOST_DEVICE static inline void ivec_Scale( ivec dest, real C, const ivec src )
{
    dest[0] = (int)(C * src[0]);
    dest[1] = (int)(C * src[1]);
    dest[2] = (int)(C * src[2]);
}


CUDA_HOST_DEVICE static inline void ivec_rScale( ivec dest, real C, const rvec src )
{
    dest[0] = (int)(C * src[0]);
    dest[1] = (int)(C * src[1]);
    dest[2] = (int)(C * src[2]);
}


CUDA_HOST_DEVICE static inline int ivec_isZero( const ivec v )
{
    if ( v[0] == 0 && v[1] == 0 && v[2] == 0 )
    {
        return TRUE;
    }

    return FALSE;
}


CUDA_HOST_DEVICE static inline int ivec_isEqual( const ivec v1, const ivec v2 )
{
    if ( v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2] )
    {
        return TRUE;
    }

    return FALSE;
}


CUDA_HOST_DEVICE static inline void ivec_Sum( ivec dest, const ivec v1, const ivec v2 )
{
    dest[0] = v1[0] + v2[0];
    dest[1] = v1[1] + v2[1];
    dest[2] = v1[2] + v2[2];
}


CUDA_HOST_DEVICE static inline void ivec_ScaledSum( ivec dest,
        int k1, const ivec v1, int k2, const ivec v2 )
{
    dest[0] = k1 * v1[0] + k2 * v2[0];
    dest[1] = k1 * v1[1] + k2 * v2[1];
    dest[2] = k1 * v1[2] + k2 * v2[2];
}


CUDA_HOST_DEVICE static inline void ivec_Add( ivec dest, const ivec v )
{
    dest[0] += v[0];
    dest[1] += v[1];
    dest[2] += v[2];
}


CUDA_HOST_DEVICE static inline void ivec_ScaledAdd( ivec dest,
        int k, const ivec v )
{
    dest[0] += k * v[0];
    dest[1] += k * v[1];
    dest[2] += k * v[2];
}



CUDA_HOST_DEVICE static inline void ivec_Max( ivec res,
        const ivec v1, const ivec v2 )
{
    res[0] = MAX( v1[0], v2[0] );
    res[1] = MAX( v1[1], v2[1] );
    res[2] = MAX( v1[2], v2[2] );
}


CUDA_HOST_DEVICE static inline void ivec_Max3( ivec res,
        const ivec v1, const ivec v2, const ivec v3 )
{
    res[0] = MAX3( v1[0], v2[0], v3[0] );
    res[1] = MAX3( v1[1], v2[1], v3[1] );
    res[2] = MAX3( v1[2], v2[2], v3[2] );
}
#endif


#ifdef __cplusplus
}
#endif


#endif
