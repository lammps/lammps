/***************************************************************************
                                    miscm.h
                             -------------------
                               W. Michael Brown

  Miscellaneous functions that do not deserve their own class

 __________________________________________________________________________
    This file is part of the Math Library
 __________________________________________________________________________

    begin                : May 30 2003
    copyright            : (C) 2003 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

#ifndef MISCM_H
#define MISCM_H

#include <float.h>
#include <math.h>
#include <vector>
#include "error.h"

using namespace std;

/// Miscellaneous functions that do not deserve their own class
/** \e a contains functions for \n
  * - simple math functions
  * - finite precision stuff */

namespace am {
	/// Returns the square of a number
	double square(double);
	/// Rounds a number
	double round(double);

  /// Return the range of elements in a vector (0 if empty)
  double range(const vector<double> &v);
  /// Return the average of elements in a vector (0 if empty)
  double mean(const vector<double> &v);
  
  /// Return the max of two objects
  template <class numtyp>
  numtyp max(numtyp one,numtyp two);

  /// Return the min of two objects
  template <class numtyp>
  numtyp min(numtyp one,numtyp two);

  /// Return the -1 for negative, 0 for zero, and 1 for positive
  int sign(double v);

  /// Swap two objects
  void swap(double a, double b);

  // Finite Precision stuff

	/// Returns a number representing zero for finite checks
	double epsilon(double number);
  /// Returns number closer to zero by the smallest interval possible
	double minus_eps(double number);
	/// Returns number farther from zero by the smallest interval possible
	double plus_eps(double number);
 	/// Returns number farther from zero by smallest interval * \b m
	double plus_Meps(double m,double number);
  /// Returns number farther from zero by smallest interval * 10^8
	double plus_2eps(double number);
  /// Returns number closer to zero by smallest interval * 10^8
	double minus_2eps(double number);

	/// Returns false if the number is stored as NAN
	bool not_nan(double number);

  // Matrix stuff

  /// Invert a matrix (Gauss-Jordan)
  /** Generates error 303 L 9 for singular matrix
    * No checking for memory limitations **/
  void invert(double **matrix, unsigned size, Error &error);

  /// Move a value from a fraction of one range to a fraction of another
  /** am::rerange(0,1,0.3,0,100) will return 30. No checking to enforce
    * that the value actually lies within the range is made. **/
  double rerange(double one_start, double one_end, double value,
                 double two_start, double two_end);
}

#endif
