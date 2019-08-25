/****************************************************************************
 * cartesian.h
 *    Shared stuff for dealing with Cartesian coordinates 
 *		Also contains quaternion structures and operations
 *
 *    Cartesian point of doubles: cPt
 *    Cartesian point of integers: iPt
 *    Vector of doubles: vectorPt
 *    Color of doubles:   colorPt
 *    Quaternion: Quaternion
 *
 *    Can be accessed as cPt.x, cPt[X], colorPt[GREEN], iPt[I], etc.
 *
 *
 * W. Michael Brown
 ****************************************************************************/

/*! \file */

#ifndef CARTESIAN_H
#define CARTESIAN_H

#include "miscm.h"
#include "m_constants.h"
#include "spherical.h"
#include <assert.h>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;

enum { X,		///<0
 			 Y,		///<1
 	     Z    ///<2
};
enum { I,   ///<0
 			 J,   ///<1
   		 K,   ///<2
   		 W		///<3
};
enum { RED,   ///<0
			 GREEN, ///<1
     	 BLUE   ///<2
};
  
// Other coordinates
template<class numtyp> class Ball;
  
// Template friend declarations
template<class numtyp> class TwoD;
template<class numtyp>
ostream & operator<< (ostream &out, const TwoD<numtyp> &t);
template<class numtyp>
istream & operator>> (istream &in, TwoD<numtyp> &t);

template<class numtyp> class ThreeD;
template<class numtyp>
ostream & operator<< (ostream &out, const ThreeD<numtyp> &t);
template<class numtyp>
istream & operator>> (istream &in, ThreeD<numtyp> &t);
template<class numtyp>
ThreeD<numtyp> operator+ (const numtyp, const ThreeD<numtyp> &two);
template<class numtyp>
ThreeD<numtyp> operator- (const numtyp, const ThreeD<numtyp> &two);
template<class numtyp>
ThreeD<numtyp> operator* (const numtyp, const ThreeD<numtyp> &two);
template<class numtyp>
ThreeD<numtyp> operator/ (const numtyp, const ThreeD<numtyp> &two);

/// Two dimensional vector
/**   The elements can be accessed directly .x or .y
  *     or by using the operator [] ( [X], [Y] or [I], [J] )
  *
  *   The following operators are currently overloaded:
  *     +,-,*
  *
  *   operators *,/ returns a vector with each element multiplied
  *   times the corresponding element in the other vector (Matlab .* )
  *
  *   the member function dot can be used to compute dot products
  *
  *   Input and output are overloaded for element I/O of the form "x y"
  *     <<, >>
  *
  * \sa cartesian.h for typedefs and defines*/
template<class numtyp>
class TwoD {
 public:
 	/// Empty construct. Not necessarily initialized to [0 0]
  TwoD();
  /// Assign both x and y the value
  TwoD(numtyp x);
  /// Assignment Constructor
  TwoD(numtyp cx, numtyp cy);
  
  /// Type conversion
  TwoD(const TwoD<float> &two);
  /// Ball projection (onto xy-plane)
  TwoD(const Ball<float> &ball);

  numtyp x; ///< First element
  numtyp y; ///< Last element

  numtyp &operator[](unsigned i);

  friend ostream & operator<< <>(ostream &out, const TwoD<numtyp> &t);
  friend istream & operator>> <>(istream &in, TwoD<numtyp> &t);

  TwoD<numtyp> operator + (const TwoD<numtyp> &two) const;
  void operator += (const TwoD<numtyp> &two);
  TwoD<numtyp> operator + (const numtyp two) const;
  TwoD<numtyp> operator - (const TwoD<numtyp> &two) const;
  TwoD<numtyp> operator - (const numtyp two) const;
  TwoD<numtyp> operator * (const numtyp two) const;
  TwoD<numtyp> operator * (const TwoD<numtyp> &two) const;
  void operator /= (const numtyp two);

  bool operator != (const TwoD<numtyp> &two) const;
  
  /// Dot Product
  numtyp dot(const TwoD<numtyp> &two) const;
  
  /// Distance between two points
  numtyp dist(const TwoD<numtyp> &two) const;
  /// Distance squared between two points
  numtyp dist2(const TwoD<numtyp> &two) const;
  /// Returns one of two normals to a line represented by vector
  TwoD<numtyp> normal();

  /// Move coordinates into array
  void to_array(numtyp *array);
  /// Set coordinates from array
  void from_array(numtyp *array);
  
  // -------------- Weird functions that help with coord templating
  unsigned dimensionality();
  /// Returns the index of *this (for unsigned) in a square 3D array
  /** \param s_size s_size[0]=1Dsize **/
  numtyp array_index(vector<unsigned> &s_size);
  /// Increment a 2D index from min to max (same as nested for)
  /** Returns false when increment is complete **/
  bool increment_index(TwoD &minp,TwoD &maxp);
  /// Return false if x or y is not within the inclusive range
  bool check_bounds(numtyp min,numtyp max);
 private:
};

///\var typedef TwoD<double> c2DPt
/// Two dimensional vector of doubles
typedef TwoD<double> c2DPt;


/// Three dimensional vector
/**   The elements can be accessed directly .x or .y or .z
  *     or by using the operator [] ( [X], [Y], [Z] or [I], [J], [K] or
  *     [RED], [GREEN], [BLUE] )
  *
  *   The following operators are currently overloaded:
  *     +,-,*,/,+=,-=,*=,/=,==,!=
  *
  *   operators *,/,*=,/= returns a vector with each element multiplied
  *   times the corresponding element in the other vector (Matlab .* )
  *   or with each element multiplied by a scalar
  *
  *   the member function dot can be used to compute dot products

  *
  *   Input and output are overloaded for element I/O of the form "x y z"
  *     <<, >>
  *
  * \sa cartesian.h for typedefs and defines*/
template<class numtyp>
class ThreeD {
 public:
 	/// Assignment Constructor 
  ThreeD(numtyp cx, numtyp cy, numtyp cz);
  /// Assign all the value
  ThreeD(numtyp in);
  /// Type Conversion
  ThreeD(const ThreeD<unsigned>&);
  /// Type Conversion
  ThreeD(const ThreeD<int>&);
  /// Type Conversion
  ThreeD(const ThreeD<float>&);
  /// Type Conversion
  ThreeD(const ThreeD<double>&);
  /// TwoD with (z-coordinate set to zero)
  ThreeD(const TwoD<float>&);
  /// Spherical Conversion
  ThreeD(const Ball<double>&);
  /// Spherical Conversion
  ThreeD(const Ball<float>&);
  /// Empty construct (Not necessarily initialized to zero)
  ThreeD();

  numtyp x;   ///< First Element
  numtyp y;   ///< Second Element
  numtyp z;   ///< Last Element

  friend ostream & operator<< <>(ostream &out, const ThreeD<numtyp> &t);
  friend istream & operator>> <>(istream &in, ThreeD<numtyp> &t);
  friend ThreeD<numtyp> operator+ <>(const numtyp, const ThreeD<numtyp> &two);
  friend ThreeD<numtyp> operator- <>(const numtyp, const ThreeD<numtyp> &two);
  friend ThreeD<numtyp> operator* <>(const numtyp, const ThreeD<numtyp> &two);
  friend ThreeD<numtyp> operator/ <>(const numtyp, const ThreeD<numtyp> &two);

  inline numtyp &operator[](unsigned i) {
    switch(i) {
    case X: return x;
    case Y: return y;
    case Z: return z;
    }
    return x;
  }

  inline numtyp operator[](unsigned i) const {
    switch(i) {
    case X: return x;
    case Y: return y;
    case Z: return z;
    }
    return x;
  }

  bool operator == (const ThreeD<numtyp> &two) const;
  bool operator != (const ThreeD<numtyp> &two) const;
  ThreeD<numtyp> operator + (const ThreeD<numtyp> &two) const;
  ThreeD<numtyp> operator + (const numtyp &two) const;
  ThreeD<numtyp> operator - (const ThreeD<numtyp> &two) const;
  ThreeD<numtyp> operator - (const numtyp &two) const;
  ThreeD<numtyp> operator * (const numtyp &two) const;
  ThreeD<numtyp> operator * (const ThreeD<numtyp> &two) const;
  ThreeD<numtyp> operator / (const numtyp &two) const;
  ThreeD<numtyp> operator / (const ThreeD<numtyp> &two) const;
  void operator = (const ThreeD &two);
  void operator += (const numtyp &two);
  void operator += (const ThreeD &two);
  void operator -= (const numtyp &two);
  void operator -= (const ThreeD &two);
  inline void operator *= (const numtyp &two) {
    x*=two; y*=two; z*=two;
  }

  void operator /= (const numtyp &two);

  /// Move coordinates into array
  void to_array(numtyp *array);
  /// Set coordinates from array
  void from_array(numtyp *array);
  
  /// Returns the dot product of *this and two
  numtyp dot(const ThreeD<numtyp> &two) const;

  /// Returns the cross product of \b *this and \b two
  /** The input vectors do not need to be normalized, however, the output
    * vector will not be */
  ThreeD<numtyp> cross(const ThreeD<numtyp> &two) const;
  /// Returns the angle between \b *this and \b two
  numtyp angle(const ThreeD<numtyp> &two) const;
  /// Returns an arbitrary vector that is perpendicular to the input
  /** The output vector is not normalized */
  ThreeD<numtyp> perpendicular();
  /// Rotate a vector in xz-plane by t radians
  ThreeD<numtyp> rotatey(double t) const;  
  /// Rotate a vector in xy-plane by t radians
  ThreeD<numtyp> rotatez(double t) const;
  /// Magnitude of vector
  numtyp norm() const;
  /// Squared norm of a vector
  numtyp norm2() const;
  /// Magnitude of vector
  numtyp hypot() const;
  /// Distance between two points
  inline numtyp dist(const ThreeD<numtyp> &two) {
    return (*this-two).norm();
  }
  /// Distance squared between two points
  numtyp dist2(const ThreeD<numtyp> &two);
  /// Converts \b *this to the unit vector
  inline void normalize() {
    numtyp temp=norm();
    #ifdef NANCHECK
    assert(temp!=0);
    #endif
    *this/=temp;
  }

  /// Return the unit vector of \b *this
  ThreeD<numtyp> unit();

  /// Returns the projection of the point or vector onto Z-plane
  c2DPt projz();

  /// Set this to be the max of one and two for each dimension
  void max(ThreeD &one, ThreeD &two);
  /// Set this to be the min of one and two for each dimension
  void min(ThreeD &one, ThreeD &two);

  // -------------- Weird functions that help with coord templating
  /// Returns 3
  unsigned dimensionality();
  /// Returns the index of *this (for unsigned) in a square 3D array
  /** \param s_size s_size[0]=1Dsize, and s_size[1]=1D size*1Dsize **/
  numtyp array_index(vector<unsigned> &s_size);
  /// Increment a 3D index from min to max (same as nested for)
  /** \note This is currently only implemented for unsigned numbers
    * Returns false when increment is complete **/
  bool increment_index(ThreeD &minp,ThreeD &maxp);
  /// Return false if x,y, or z is not within the inclusive range
  bool check_bounds(numtyp min,numtyp max);
 private:
};

///\var typedef ThreeD<int> iPt;
/// Three dimensional vector of integers
typedef ThreeD<int> iPt;
///\var typedef ThreeD<unsigned> uPt;
/// Three dimensional vector of unsigned
typedef ThreeD<unsigned> uPt;
///\var typedef ThreeD<double> cPt;
/// Three dimensional vector of doubles
typedef ThreeD<double> cPt;
///\var typedef ThreeD<double> vectorPt;
/// Three dimensional vector of doubles
typedef ThreeD<double> vectorPt;
///\var typedef ThreeD<double> colorPt;
/// Three dimensional vector of doubles
typedef ThreeD<double> colorPt;

///\def ORIGIN
/// Point at origin
#define ORIGIN cPt(0.0,0.0,0.0)
///\def XAXIS
/// Unit vector for x-axis
#define XAXIS vectorPt(1.0,0.0,0.0)
///\def YAXIS
/// Unit vector for y-axis
#define YAXIS vectorPt(0.0,1.0,0.0)
///\def ZAXIS
/// Unit vector for z-axis
#define ZAXIS vectorPt(0.0,0.0,1.0)

class RotMat;

/// Euler Rotation
class EulerRot {
 public:
  EulerRot();
  EulerRot(double theta,double psi,double phi);

  /// Rotate the rotation axis by a given rotation matrix
  void rotate_axis(const RotMat &rotmat);
  
  void operator= (const EulerRot &two);
  friend ostream & operator<<(ostream &out, const EulerRot &rot);
  
	double theta;
	double psi;
	double phi;
};

///\def EU_NOROTATE
/// EulerRot with no rotation
#define EU_NOROTATE EulerRot(0.0,0.0,0.0)

//---------------------------Quaternion Stuff -------------------------

/// Class for handling Quaternion vectors
/**   The elements can be accessed directly .w .i .j or .k
  *     or by using the operator[] ( [W], [X], [Y], [Z] or [W], [I], [J], [K] )
  *
  *   The following operators are currently overloaded:
  *     +=,*,*=,=,==
  *
  *   Overloaded * concatenates two successive rotations \n
  *     \e It \e is \e not \e communative \n
  *     For \e join=q1*q2, \e join is equivalent to performing rotation \e q1
  *     \b followed by \e q2.
  *
  *   Input and output are overloaded for element I/O of the form "w i j k"
  *     <<, >>
  *
  * \sa cartesian.h for typedefs and defines*/
class Quaternion {
 public:
  Quaternion();
  ~Quaternion();

	double w; ///<Real Component
  double i; ///<Imaginary Component
  double j; ///<Imaginary Component
  double k; ///<Imaginary Component
  
	double &operator[](unsigned index);
	
  /// Copy Constructor
  Quaternion(const Quaternion &two);
  /// Assignment constructor
  Quaternion(double inw, double ini, double inj, double ink);
	/// Axis-Angle Rotation (\b v must be a \b UNIT vector)
	Quaternion(const vectorPt &v, double angle);
  /// Rotation from vector 1 to vector 2 (v1 and v2 need not be normalized)
  /** Angle is between \b v1 and \b v2 and axis is calculated as the cross
    * product of \b v1 and \b v2 */
	Quaternion(const vectorPt &v1, const vectorPt &v2);
	/// Spherical rotation
	Quaternion(double longitude, double latitude);
  /// From Euler angles
  Quaternion(const EulerRot &erot);
	
  void operator += (const Quaternion &two);
  Quaternion operator * (const Quaternion &two) const;
  void operator*= (const Quaternion &two);
  Quaternion operator * (const double two) const;
  void operator*= (const double two);
  void operator= (const Quaternion &two);
  bool operator== (const Quaternion &two) const;
  
  /// Returns the conjugate
  Quaternion conj() const;
  /// Returns the norm
  double norm();
  /// Form \b *this into the unit vector
  void normalize();
  /// Returns the unit vector of \b *this
  Quaternion unit();
  
  friend ostream & operator<<(ostream &out, const Quaternion &q);
  friend istream & operator>>(istream &in, Quaternion &q);
};

/// Rotation matrix (about origin)
/** Rotation of ThreeD can be applied via overloaded * operator **/
class RotMat {
 public:
  /// Unspecified matrix!
  RotMat();
  /// Initialize with quaternion
  RotMat(const Quaternion &q);
  /// Set based on quaternion
  void set(const Quaternion &q);
  /// Assert that the rotation is proper
  void proper();

  cPt operator*(const cPt &in) const;

  double x_x,x_y,x_z,y_x,y_y,y_z,z_x,z_y,z_z;
};  

///\def NOROTATE
/// Quaternion with no rotation
#define NOROTATE Quaternion(1.0,0.0,0.0,0.0)

/// Global geometry functions
namespace c {
	/// Returns point on a line closest to an outside point
	/** The line is described by a directional vector and a point on the line
	  *\param v Vector describing the line direction
	  *\param v_point Point on the line
	  *\param point The Outside point */
  cPt point_to_line(const vectorPt &v, const cPt &v_point,
									  const cPt &point);

  /// Return the closest distance between two line segments input as end points
  /**\param l1_1 End point of segment 1
    *\param l1_2 End point of segment 1
    *\param l2_1 End point of segment 2
    *\param l2_2 End point of segment 2 */
  double closest_approach(const cPt &l1_1, const cPt &l1_2,
													const cPt &l2_1, const cPt &l2_2);
  /// Calculates the points where two line segments are closest
  /**\param l1_1 End point of segment 1
    *\param l1_2 End point of segment 1
    *\param l2_1 End point of segment 2
    *\param l2_2 End point of segment 2
    *\param close_l1 Calculated closest point on line segment 1
    *\param close_l2 Calculated closest point on line segment 2 */
  void closest_approach_points(const cPt &l1_1, const cPt &l1_2,
															 const cPt &l2_1, const cPt &l2_2,
														   cPt &close_l1, cPt &close_l2);

  /// Returns true if two line segments intersect
  bool intersect(const c2DPt &line1_start, const c2DPt &line1_end,
								 const c2DPt &line2_start, const c2DPt &line2_end);
	
	/// Liang-Barsky Intersection between line segment and line
	bool sline_intersect(const c2DPt &line1_start, const c2DPt &line1_end,
											 const c2DPt &line2normal, const c2DPt &line2_point);
  /// Average position
  cPt mean(vector<cPt> &vec);
}
 
#endif
