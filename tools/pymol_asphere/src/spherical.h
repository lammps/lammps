/***************************************************************************
                                spherical.h
														  W. Michael Brown
                             -------------------

    Stuff for working spherical coordinates

 __________________________________________________________________________

    Part of the Math Library
 __________________________________________________________________________

    begin                : Tue Aug 29 2006
    copyright            : (C) 2006 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

#ifndef SPHERICAL_H
#define SPHERICAL_H

#include "miscm.h"
#include "m_constants.h"
#include "cartesian.h"
#include <assert.h>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;

// Other coordinates
template<class numtyp> class ThreeD;

// Friends
template<class numtyp> class Ball;
template<class numtyp>
ostream & operator<< (ostream &out, const Ball<numtyp> &t);
template<class numtyp>
istream & operator>> (istream &in, Ball<numtyp> &t);

enum { THETA,		///<0
 			 PHI		///<1
};

/// Two dimensional spherical coordinates on a unit sphere
/**   The elements can be accessed directly .theta or .phi
  *     or by using the operator [] ( [THETA], [PHI] )
  *
  *   Input and output are overloaded for element I/O of the form "theta phi"
  *     <<, >>
  **/
template<class numtyp>
class Ball {
 public:
 	/// Empty construct. Not necessarily initialized to [0 0]
  Ball();
  /// Assignment Constructor
  Ball(numtyp theta, numtyp phi);
  /// Assign theta and phi to the value
  Ball(numtyp value);
  /// Convert from cartesian
  Ball(const ThreeD<numtyp> &pt);

  numtyp theta;
	numtyp phi;

  numtyp &operator[](unsigned i);

  friend ostream & operator<< <>(ostream &out, const Ball &t);
  friend istream & operator>> <>(istream &in, Ball &t);

  /// Add both angles
  void operator += (const Ball<numtyp> &two);
  /// Add to both angles
  Ball<numtyp> operator + (const numtyp two) const;
  /// Multiply both angles
  Ball<numtyp> operator * (const numtyp two) const;
  /// Divide both angles
  void operator /= (const numtyp two);
  /// Add both angles
  Ball<numtyp> operator + (const Ball<numtyp> &two);

  /// Distance between two points (along arc)
	/** \note The form of calculation used suffers from round off error
	  * when points are antipodal **/
  numtyp dist(const Ball &two) const;
  /// Distance squared between two points (along arc)
	/** \note The form of calculation used suffers from round off error
	  * when points are antipodal **/
  numtyp dist2(const Ball &two) const;

  /// Move coordinates into array
  void to_array(numtyp *array);
  /// Set coordinates from array
  void from_array(numtyp *array);

  // -------------- Weird functions that help with coord templating
  /// Returns 2
	unsigned dimensionality();
  // Returns true
  bool check_bounds(numtyp min,numtyp max);
 private:
};

///\var typedef Ball<double> BallD
/// Double unit sphere
typedef Ball<double> BallD;
///\var typedef Ball<double> BallF
/// Float unit sphere
typedef Ball<float> BallF;



#endif
