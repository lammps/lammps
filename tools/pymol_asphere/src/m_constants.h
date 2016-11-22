/***************************************************************************
                                m_constants.h
                             -------------------
                               W. Michael Brown

  Misc constants

 __________________________________________________________________________
    This file is part of the Math Library
 __________________________________________________________________________

    begin                : Wed Aug 10 2005
    copyright            : (C) 2005 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

/*! \file */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <math.h>
using namespace std;

#define MATHLIB_VER "0.15"

#ifndef PI
///\def PI
/// pi
#define PI 3.1415926535897932384626433832795
#endif
///\def TWOPI
/// pi*2
#define TWOPI 6.2831853071795862
///\def HALFPI
/// pi/2
#define HALFPI 1.5707963267948966
///\def DEGTORAD
/// Convert Degrees to Radians (pi/180)
#define DEGTORAD 0.017453292519943295
///\def SQRT_TWO
/// sqrt(2.0)
#define SQRT_TWO 1.4142135623730951
///\def SQRT_PI
/// sqrt(PI)
#define SQRT_PI 1.7724538509055159
///\def INF
/// Infinity
#define INF 1e308
///\def MINUSINF
/// Negative infinity
#define MINUSINF -1e308

#ifndef EPS
///\def EPS
/// Small number
#define EPS 1e-100
#endif

/** \mainpage Math Library
  * \section intro Introduction
  * Math library with containers and operations for vectors, matrices, graphs,
  * cartesian coordinates, quaternions, Euler angles, support vector machine
  * models, etc. **/


#endif
