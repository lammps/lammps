/***************************************************************************
                               spherical.cpp
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

#include "spherical.h"

// Empty construct. Not necessarily initialized to [0 0]
template<class numtyp>
Ball<numtyp>::Ball() {
}

// Assign theta and phi to the value
template<class numtyp>
Ball<numtyp>::Ball(numtyp value) : theta(value), phi(value) {
}

// Assignment Constructor
template<class numtyp>
Ball<numtyp>::Ball(numtyp thet, numtyp ph) : theta(thet), phi(ph) {
}

// Convert from cartesian
template<class numtyp>
Ball<numtyp>::Ball(const ThreeD<numtyp> &pt) {
  theta=atan2(pt.y,pt.x);
  if (theta<0)
    theta+=TWOPI;
  phi=acos(pt.z/pt.norm());
}

template<class numtyp>
numtyp & Ball<numtyp>::operator[](unsigned i) {
  if (i==THETA)
	  return theta;
	else
	  return phi;
}

template<class numtyp>
ostream & operator<< (ostream &out, const Ball<numtyp> &t) {
  out << t.theta << " " << t.phi;
  return out;
}

template<class numtyp>
istream & operator>> (istream &in, Ball<numtyp> &t) {
  in >> t.theta >> t.phi;
  return in;
}

// Distance between two points (along arc)
template<class numtyp>
numtyp Ball<numtyp>::dist(const Ball &two) const {
  double dot=cPt(*this).dot(cPt(two));
  if (dot>1.0)
    dot=1.0;
  return acos(dot);
}

// Distance squared between two points
template<class numtyp>
numtyp Ball<numtyp>::dist2(const Ball &two) const {
  numtyp d=dist(two);
	return d*d;
}

// Add both angles
template<class numtyp>
void Ball<numtyp>::operator += (const Ball<numtyp> &two) {
  theta+=two.theta;
  phi+=two.phi;
}

// Add to both angles
template<class numtyp>
Ball<numtyp> Ball<numtyp>::operator + (const numtyp two) const {
  return Ball(theta+two,phi+two);
}

// Multiply both angles
template<class numtyp>
Ball<numtyp> Ball<numtyp>::operator * (const numtyp two) const {
  return Ball(theta*two,phi*two);
}

// Divide both angles
template<class numtyp>
void Ball<numtyp>::operator /= (const numtyp two) {
  theta/=two;
  phi/=two;
}

// Add both angles
template<class numtyp>
Ball<numtyp> Ball<numtyp>::operator + (const Ball<numtyp> &two) {
  return Ball(theta+two.theta,phi+two.phi);
}

// Move coordinates into array
template<class numtyp>
void Ball<numtyp>::to_array(numtyp *array) {
  *array=theta;
  array++;
  *array=phi;
}

// Set coordinates from array
template<class numtyp>
void Ball<numtyp>::from_array(numtyp *array) {
  theta=*array;
  array++;
  phi=*array;
}

// Returns 2
template<class numtyp>
unsigned Ball<numtyp>::dimensionality() {
  return 2;
}

// Returns true
template<class numtyp>
bool Ball<numtyp>::check_bounds(numtyp min,numtyp max) {
  return true;
}

template class Ball<double>;
template class Ball<float>;
template ostream & operator<< <double>(ostream &out, const Ball<double> &t);
template ostream & operator<< <float>(ostream &out, const Ball<float> &t);
template istream & operator>> <double>(istream &in, Ball<double> &t);
template istream & operator>> <float>(istream &in, Ball<float> &t);


