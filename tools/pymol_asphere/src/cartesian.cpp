/****************************************************************************
 * cartesian.cpp
 *    Shared stuff for dealing with Cartesian coordinates
 *		Also contains quaternion structures and operations
 *
 *    Cartesian point of doubles: cPt
 *    Cartesian point of intergers: iPt
 *    Vector of doubles: vectorPt
 *    Color of doubles:   colorPt
 *    Quaternion: Quaternion
 *
 *    Can be accessed as cPt.x, cPt[X], colorPt[GREEN], iPt[I], etc.
 *
 *
 * W. Michael Brown
 * 5/26/03
 ****************************************************************************/

#include "cartesian.h"


// ------------------------ TwoD Stuff

template<class numtyp>
TwoD<numtyp>::TwoD() {
}

// Assign both x and y the value
template<class numtyp>
TwoD<numtyp>::TwoD(numtyp in) {
  x=in;
  y=in;
}

// Assignment Constructor
template<class numtyp>
TwoD<numtyp>::TwoD(numtyp cx, numtyp cy) {
	x=cx;
	y=cy;
}

// Type conversion
template<class numtyp>
TwoD<numtyp>::TwoD(const TwoD<float> &two) {
	x=numtyp(two.x);
	y=numtyp(two.y);
}

// Ball projection (onto xy-plane)
template<class numtyp>
TwoD<numtyp>::TwoD(const Ball<float> &ball) {
  x=numtyp(cos(ball.theta));
  y=numtyp(sin(ball.theta));
}  

template<class numtyp>
numtyp & TwoD<numtyp>::operator[](unsigned i) {
	switch(i) {
   case X: return x;
   case Y: return y;
  }
  return x;
}

template<class numtyp>
TwoD<numtyp> TwoD<numtyp>::operator + (const TwoD<numtyp> &two) const {
	return TwoD(x+two.x, y+two.y);
}

template<class numtyp>
void TwoD<numtyp>::operator += (const TwoD<numtyp> &two) {
  x+=two.x;
  y+=two.y;
}

template<class numtyp>
TwoD<numtyp> TwoD<numtyp>::operator + (const numtyp two) const {
	return TwoD(x+two, y+two);
}

template<class numtyp>
TwoD<numtyp> TwoD<numtyp>::operator - (const TwoD<numtyp> &two) const {
	return TwoD(x-two.x, y-two.y);
}

template<class numtyp>
TwoD<numtyp> TwoD<numtyp>::operator - (const numtyp two) const {
	return TwoD(x-two, y-two);
}

template<class numtyp>
TwoD<numtyp> TwoD<numtyp>::operator * (const numtyp two) const {
	return TwoD(x*two, y*two);
}

template<class numtyp>
TwoD<numtyp> TwoD<numtyp>::operator * (const TwoD<numtyp> &two) const {
	return TwoD<numtyp>(x*two.x,y*two.y);
}

template<class numtyp>
void TwoD<numtyp>::operator /= (const numtyp two) {
  x/=two;
  y/=two;
}

template<class numtyp>
bool TwoD<numtyp>::operator != (const TwoD<numtyp> &two) const {
  if (x!=two.x || y!=two.y)
  	return true;
  return false;
}

// Move coordinates into array
template<class numtyp>
void TwoD<numtyp>::to_array(numtyp *array) {
  *array=x;
  array++;
  *array=y;
}
  
// Set coordinates from array
template<class numtyp>
void TwoD<numtyp>::from_array(numtyp *array) {
  x=*array;
  array++;
  y=*array;
}

// Dot Product
template<class numtyp>
numtyp TwoD<numtyp>::dot(const TwoD<numtyp> &two) const {
	return x*two.x+y*two.y;
}

// Returns one of two normals to a line represented by vector
template<class numtyp>
TwoD<numtyp> TwoD<numtyp>::normal() {
	return TwoD(y,x*numtyp(-1));
}

template<class numtyp>
numtyp TwoD<numtyp>::dist2(const TwoD<numtyp> &two) const {
	TwoD diff=*this-two;
	return diff.dot(diff);
}

// Returns two
template<class numtyp>
unsigned TwoD<numtyp>::dimensionality() {
  return 2;
}

// Returns the index of *this (for unsigned) in a square 3D array
/* s_size[0]=1Dsize */
template<class numtyp>
numtyp TwoD<numtyp>::array_index(vector<unsigned> &s_size) {
  return y+x*s_size[0];
}

// Increment a 3D index from min to max (same as nested for)
template <>
bool TwoD<unsigned>::increment_index(TwoD &minp,TwoD &maxp) {
  x++;
  if (x!=maxp.x)
    return true;
  x=minp.x;
  y++;
  if (y!=maxp.y)
    return true;
  return false;
}

// Return false if x or y is not within the inclusive range
template<class numtyp>
bool TwoD<numtyp>::check_bounds(numtyp minp,numtyp maxp) {
  if (x<minp || y<minp || x>maxp || y>maxp)
    return false;
  return true;
}

template<class numtyp>
ostream & operator<<(ostream &out, const TwoD<numtyp> &t) {
	out << t.x << " " << t.y;
	return out;
}

template<class numtyp>
istream & operator>>(istream &in, TwoD<numtyp> &t) {
	in >> t.x >> t.y;
	return in;
}

// ------------------------ ThreeD Stuff

// Assignment Constructor
template<class numtyp>
ThreeD<numtyp>::ThreeD(numtyp cx, numtyp cy, numtyp cz) {
	x=cx; y=cy; z=cz;
}

// Assignment Constructor
template<class numtyp>
ThreeD<numtyp>::ThreeD(numtyp in) {
	x=in; y=in; z=in;
}

// Type Conversion
template<class numtyp>
ThreeD<numtyp>::ThreeD(const ThreeD<unsigned>& in) {
  x=numtyp(in.x);  y=numtyp(in.y);  z=numtyp(in.z);
}

// Type Conversion
template<class numtyp>
ThreeD<numtyp>::ThreeD(const ThreeD<int>& in) {
  x=numtyp(in.x);  y=numtyp(in.y);  z=numtyp(in.z);
}

// Type Conversion
template<class numtyp>
ThreeD<numtyp>::ThreeD(const ThreeD<float>& in) {
  x=numtyp(in.x);  y=numtyp(in.y);  z=numtyp(in.z);
}

// Type Conversion
template<class numtyp>
ThreeD<numtyp>::ThreeD(const ThreeD<double>& in) {
  x=numtyp(in.x);  y=numtyp(in.y);  z=numtyp(in.z);
}

// TwoD with (z-coordinate set to zero)
template<class numtyp>
ThreeD<numtyp>::ThreeD(const TwoD<float>& in) {
  x=numtyp(in.x);  y=numtyp(in.y);  z=0;
}

// Spherical Conversion
template<class numtyp>
ThreeD<numtyp>::ThreeD(const Ball<double>& ball) {
  double sp=sin(ball.phi);
  x=numtyp(sp*cos(ball.theta));
  y=numtyp(sp*sin(ball.theta));
  z=numtyp(cos(ball.phi));
}

// Spherical Conversion
template<class numtyp>
ThreeD<numtyp>::ThreeD(const Ball<float>& ball) {
  double sp=sin(ball.phi);
  x=numtyp(sp*cos(ball.theta));
  y=numtyp(sp*sin(ball.theta));
  z=numtyp(cos(ball.phi));
}

template<class numtyp>
ThreeD<numtyp>::ThreeD() {
}

template<class numtyp>
void ThreeD<numtyp>::operator = (const ThreeD &two) {
  #ifdef NANCHECK
  //assert(a::not_nan(two.x) && a::not_nan(two.y) && a::not_nan(two.z));
  #endif
  x=two.x; y=two.y; z=two.z;
}

template<class numtyp>
bool ThreeD<numtyp>::operator == (const ThreeD<numtyp> &two) const {
  if (x!=two.x || y!=two.y || z!=two.z)
    return false;
  return true;
}

template<class numtyp>
bool ThreeD<numtyp>::operator != (const ThreeD<numtyp> &two) const {
  if (x!=two.x || y!=two.y || z!=two.z)
    return true;
  return false;
}

template<class numtyp>
ThreeD<numtyp> ThreeD<numtyp>::operator + (const ThreeD<numtyp> &two) const {
  return (ThreeD<numtyp>(x+two.x,y+two.y,z+two.z));
}

template <class numtyp>
ThreeD<numtyp> ThreeD<numtyp>::operator + (const numtyp &two) const {
  return (ThreeD<numtyp>(x+two,y+two,z+two));
}

template<class numtyp>
void ThreeD<numtyp>::operator += (const numtyp &two) {
  x+=two; y+=two; z+=two;
}

template<class numtyp>
void ThreeD<numtyp>::operator += (const ThreeD &two) {
  x+=two.x; y+=two.y; z+=two.z;
}

template<class numtyp>
void ThreeD<numtyp>::operator -= (const ThreeD &two) {
  x-=two.x; y-=two.y; z-=two.z;
}

template <class numtyp>
ThreeD<numtyp> ThreeD<numtyp>::operator - (const ThreeD<numtyp> &two) const {
  return (ThreeD<numtyp>(x-two.x,y-two.y,z-two.z));
}

template <class numtyp>
ThreeD<numtyp> ThreeD<numtyp>::operator - (const numtyp &two) const {
  return (ThreeD<numtyp>(x-two,y-two,z-two));
}

template<class numtyp>
void ThreeD<numtyp>::operator -= (const numtyp &two) {
  x-=two; y-=two; z-=two;
}

template <class numtyp>
ThreeD<numtyp> ThreeD<numtyp>::operator * (const numtyp &two) const {
  return (ThreeD<numtyp>(x*two,y*two,z*two));
}

template <class numtyp>
ThreeD<numtyp> ThreeD<numtyp>::operator * (const ThreeD<numtyp> &two) const {
	return ThreeD<numtyp>(x*two.x,y*two.y,z*two.z);
}

template <class numtyp>
ThreeD<numtyp> ThreeD<numtyp>::operator / (const ThreeD<numtyp> &two) const {
	return ThreeD<numtyp>(x/two.x,y/two.y,z/two.z);
}

// Dot Product
template <class numtyp>
numtyp ThreeD<numtyp>::dot(const ThreeD<numtyp> &two) const {
	return (x*two.x+y*two.y+z*two.z);
}

// Move coordinates into array
template<class numtyp>
void ThreeD<numtyp>::to_array(numtyp *array) {
  *array=x;
  array++;
  *array=y;
  array++;
  *array=z;
}

// Set coordinates from array
template<class numtyp>
void ThreeD<numtyp>::from_array(numtyp *array) {
  x=*array;
  array++;
  y=*array;
  array++;
  z=*array;
}


// Return the cross product of *this and two
template<class numtyp>
ThreeD<numtyp> ThreeD<numtyp>::cross(const ThreeD<numtyp> &two) const {
	return ThreeD<numtyp> (y*two.z-z*two.y, z*two.x-x*two.z, x*two.y-y*two.x);
}

// Return the angle between two vectors
template<class numtyp>
numtyp ThreeD<numtyp>::angle(const ThreeD<numtyp> &two) const {
	#ifdef NANCHECK
	//assert(fabs(double((*this*two)/(norm()*two.norm())))<=1);
	#endif
	return(numtyp(acos( double(dot(two)/(norm()*two.norm())) )));
}

// Returns an arbitrary vector that is perpendicular to the input
template <class numtyp>
ThreeD<numtyp> ThreeD<numtyp>::perpendicular() {
  ThreeD<numtyp> two=*this;
	if (y==0 && z==0)
		two.y+=1;
	else
		two.x+=1;
	return cross(two);
}

template <class numtyp>
ThreeD<numtyp> ThreeD<numtyp>::rotatey(double t) const {
	ThreeD<numtyp> two;
  numtyp sint=numtyp(sin(t));
  numtyp cost=numtyp(cos(t));
	two.x=(x*cost)-(z*sint);
	two.z=(x*sint)+(z*cost);
	two.y=y;
	return two;
}

template <class numtyp>
ThreeD<numtyp> ThreeD<numtyp>::rotatez(double t) const {
	ThreeD<numtyp> two;
  numtyp sint=numtyp(sin(t));
  numtyp cost=numtyp(cos(t));
	two.x=(x*cost)-(y*sint);
	two.y=(x*sint)+(y*cost);
	two.z=z;
	return two;
}

template <class numtyp>
ThreeD<numtyp> ThreeD<numtyp>::operator / (const numtyp &two) const {
  return (ThreeD<numtyp>(x/two,y/two,z/two));
}

template<class numtyp>
void ThreeD<numtyp>::operator /= (const numtyp &two) {
  x/=two; y/=two; z/=two;
}

// Magnitude of vector
template <class numtyp>
numtyp ThreeD<numtyp>::hypot() const {
  return norm();
}

// Return the norm of a vector
template<class numtyp>
numtyp ThreeD<numtyp>::norm() const {
  numtyp xi=x*x; numtyp yi=y*y; numtyp zi=z*z;
  return numtyp(sqrt(double(xi+yi+zi)));
}

// Return the squared norm of a vector
template<class numtyp>
numtyp ThreeD<numtyp>::norm2() const {
  return x*x+y*y+z*z;
}

template <class numtyp>
numtyp ThreeD<numtyp>::dist2(const ThreeD<numtyp> &two) {
	ThreeD diff=*this-two;
	return diff.dot(diff);
}

// Return unit vector
template<class numtyp>
ThreeD<numtyp> ThreeD<numtyp>::unit() {
	ThreeD<numtyp> unit=*this;
	unit.normalize();
	return unit;
}

/// Returns 3
template<class numtyp>
unsigned ThreeD<numtyp>::dimensionality() {
  return 3;
}

template<class numtyp>
c2DPt ThreeD<numtyp>::projz() {
	return c2DPt(x,y);
}

template<class numtyp>
ostream & operator<< (ostream &out, const ThreeD<numtyp> &t) {
	out << t.x << " " << t.y << " " << t.z;
	return out;
}

template<class numtyp>
istream & operator>>(istream &in, ThreeD<numtyp> &t) {
	in >> t.x >> t.y >> t.z;
	return in;
}

template<class numtyp>
ThreeD<numtyp> operator+ (const numtyp one, const ThreeD<numtyp> &two) {
	return ThreeD<numtyp>(one+two.x,one+two.y,one+two.z);
}

template<class numtyp>
ThreeD<numtyp> operator- (const numtyp one, const ThreeD<numtyp> &two) {
	return ThreeD<numtyp>(one-two.x,one-two.y,one-two.z);
}

template<class numtyp>
ThreeD<numtyp> operator* (const numtyp one, const ThreeD<numtyp> &two) {
	return ThreeD<numtyp>(one*two.x,one*two.y,one*two.z);
}

template<class numtyp>
ThreeD<numtyp> operator/ (const numtyp one, const ThreeD<numtyp> &two) {
	return ThreeD<numtyp>(one/two.x,one/two.y,one/two.z);
}

// Set this to be the max of one and two for each dimension
template<class numtyp>
void ThreeD<numtyp>::max(ThreeD<numtyp> &one, ThreeD<numtyp> &two) {
	x=am::max(one.x,two.x);
	y=am::max(one.y,two.y);
	z=am::max(one.z,two.z);
}

// Set this to be the max of one and two for each dimension
template<class numtyp>
void ThreeD<numtyp>::min(ThreeD<numtyp> &one, ThreeD<numtyp> &two) {
	x=am::min(one.x,two.x);
	y=am::min(one.y,two.y);
	z=am::min(one.z,two.z);
}

// Returns the index of *this (for unsigned) in a square 3D array
/* s_size[0]=1Dsize, and s_size[1]=1D size*1Dsize */
template<class numtyp>
numtyp ThreeD<numtyp>::array_index(vector<unsigned> &s_size) {
  return z+y*s_size[0]+x*s_size[1];
}

// Increment a 3D index from min to max (same as nested for)
template <>
bool ThreeD<unsigned>::increment_index(ThreeD &minp,ThreeD &maxp) {
  x++;
  if (x!=maxp.x)
    return true;
  x=minp.x;
  y++;
  if (y!=maxp.y)
    return true;
  y=minp.y;
  z++;
  if (z!=maxp.z)
    return true;
  return false;
}

// Return false if x,y, or z is not within the inclusive range
template<class numtyp>
bool ThreeD<numtyp>::check_bounds(numtyp minp,numtyp maxp) {
  if (x<minp || y<minp || z<minp || x>maxp || y>maxp || z>maxp)
    return false;
  return true;
}

// --------------------- Quaternion Stuff -----------------------

Quaternion::Quaternion() {
	#ifdef NANCHECK
	w=0; i=0; j=0; k=0;
	#endif
}

Quaternion::~Quaternion() {
}

Quaternion::Quaternion(const Quaternion &two) {
  #ifdef NANCHECK
  //assert(a::not_nan(two.w) && a::not_nan(two.i) && a::not_nan(two.j) &&
  //			 a::not_nan(two.k));
  #endif
	w=two.w; i=two.i; j=two.j; k=two.k;
}

Quaternion::Quaternion(double inw, double ini, double inj, double ink) {
  #ifdef NANCHECK
  //assert(a::not_nan(inw) && a::not_nan(ini) && a::not_nan(inj) &&
  //			 a::not_nan(ink));
  #endif
	w=inw; i=ini; j=inj; k=ink;
}

// Set the quaterion according to axis-angle (must be a UNIT vector)
Quaternion::Quaternion(const vectorPt &v, double angle) {
  double halfa=angle/2;
	double sina=sin(halfa);
	w=cos(halfa);
	i=v.x*sina;
	j=v.y*sina;
	k=v.z*sina;
  #ifdef NANCHECK
  //assert(a::not_nan(w) && a::not_nan(i) && a::not_nan(j) && a::not_nan(k));
  #endif
}

// Set the quaterion according to rotation from vector 1 to vector 2
//    The input vectors do not need to be normalized
Quaternion::Quaternion(const vectorPt &v1, const vectorPt &v2) {
	double angle=acos( v1.dot(v2)/(v1.norm()*v2.norm()) );
	vectorPt axis=v1.cross(v2);
  double norm=axis.norm();
  if (norm<EPS) {
    *this=NOROTATE;
    return;
  }    
    
	axis/=norm;
  double halfa=angle/2;
	double s=sin(halfa);

	w=cos(halfa);
	i=axis.x*s;
	j=axis.y*s;
	k=axis.z*s;
  #ifdef NANCHECK
  //assert(a::not_nan(w) && a::not_nan(i) && a::not_nan(j) && a::not_nan(k));
  #endif
}

// From Euler angles
Quaternion::Quaternion(const EulerRot &erot) {
  double c1 = cos(erot.theta/2);
  double s1 = sin(erot.theta/2);
  double c2 = cos(erot.psi/2);
  double s2 = sin(erot.psi/2);
  double c3 = cos(erot.phi/2);
  double s3 = sin(erot.phi/2);
  double c1c2 = c1*c2;
  double s1s2 = s1*s2;
  w =c1c2*c3 - s1s2*s3;
  i =c1c2*s3 + s1s2*c3;
	j =s1*c2*c3 + c1*s2*s3;
	k =c1*s2*c3 - s1*c2*s3;
}

double & Quaternion::operator[](unsigned index) {
	switch(index) {
	 case I: return i;
   case J: return j;
   case K: return k;
   case W: return w;
  }
  return w;
}

Quaternion operator + (const Quaternion one, const Quaternion two) {
	Quaternion q(one.w+two.w,one.i+two.i,one.j+two.j,one.k+two.k);
	return q;
}

void Quaternion::operator += (const Quaternion &two) {
  w+=two.w; i+=two.i; j+=two.j; k+=two.k;
}

Quaternion Quaternion::operator* (const Quaternion &two) const {
	Quaternion q;
	q.w=two.w*w-two.i*i-two.j*j-two.k*k;
	q.i=two.w*i+two.i*w+two.j*k-two.k*j;
	q.j=two.w*j-two.i*k+two.j*w+two.k*i;
	q.k=two.w*k+two.i*j-two.j*i+two.k*w;

	return (q);
}

void Quaternion::operator*=(const Quaternion &two) {
	(*this)=(*this)*two;
}

Quaternion Quaternion::operator* (const double two) const {
	Quaternion q;
	q.w=two*w;
	q.i=two*i;
	q.j=two*j;
	q.k=two*k;

	return (q);
}

void Quaternion::operator*=(const double two) {
  w*=two;
  i*=two;
  j*=two;
  k*=two;
}

void Quaternion::operator = (const Quaternion &two) {
  #ifdef NANCHECK
  //assert(a::not_nan(two.w) && a::not_nan(two.i) && a::not_nan(two.j) &&
  //			 a::not_nan(two.k));
  #endif
  w=two.w; i=two.i; j=two.j; k=two.k; return;
}

bool Quaternion::operator == (const Quaternion &two) const {
  if (w!=two.w || i!=two.i || j!=two.j || k!=two.k)
		return false;
  return true;
}

Quaternion Quaternion::conj() const {
	Quaternion q;
	q.w=w; q.i=i*-1.0; q.j=j*-1.0; q.k=k*-1.0;
	return q;
}

double Quaternion::norm() {
  return (sqrt(w*w+i*i+j*j+k*k));
}

void Quaternion::normalize() {
  double temp=norm();

  if (temp!=0) {
    w/=temp; i/=temp; j/=temp; k/=temp;
  } else {
    w=1; i=0; j=0; k=0;
  }
}

Quaternion Quaternion::unit() {
	Quaternion unit=*this;
	unit.normalize();
	return unit;
}

ostream & operator<<(ostream &out, const Quaternion &q) {
	out << q.w << " " << q.i << " " << q.j << " " << q.k << " ";
	return out;
}

istream & operator>>(istream &in, Quaternion &q) {
	in >> q.w >> q.i >> q.j >> q.k;
	return in;
}

// ------------------- RotMat Stuff
RotMat::RotMat() {
}

RotMat::RotMat(const Quaternion &q) {
	set(q);
}

void RotMat::set(const Quaternion &q) {
 	// x'=x(w^2+i^2-k^2-j^2)+y(2ij-2kw)+z(2jw+2ik)
	// y'=x(2ij+2kw)+y(j^2-k^2+w^2-i^2)+z(2jk-2iw)
	// z'=x(2ik-2jw)+y(2jk+2iw)+z(k^2-j^2-i^2+w^2)
	double w2=q.w*q.w;
	double i2=q.i*q.i;
	double j2=q.j*q.j;
	double k2=q.k*q.k;
  double twoij=2.0*q.i*q.j;
  double twoik=2.0*q.i*q.k;
  double twojk=2.0*q.j*q.k;
  double twoiw=2.0*q.i*q.w;
  double twojw=2.0*q.j*q.w;
  double twokw=2.0*q.k*q.w;

  x_x=w2+i2-j2-k2;
  x_y=twoij-twokw;
  x_z=twojw+twoik;

  y_x=twoij+twokw;
  y_y=w2-i2+j2-k2;
	y_z=twojk-twoiw;

	z_x=twoik-twojw;
	z_y=twojk+twoiw;
	z_z=w2-i2-j2+k2;
}

cPt RotMat::operator*(const cPt &in) const {
  cPt out;
  out.x=x_x*in.x+x_y*in.y+x_z*in.z;
  out.y=y_x*in.x+y_y*in.y+y_z*in.z;
  out.z=z_x*in.x+z_y*in.y+z_z*in.z;
  return out;
}

void RotMat::proper() {
  double det=x_x*y_y*z_z-x_x*y_z*z_y-x_y*y_x*z_z+x_y*y_z*z_x+x_z*y_x*z_y-x_z*
             x_y*z_x;
  if (det<0) {
    x_x*=-1.0;
    x_y*=-1.0;
    x_z*=-1.0;
    y_x*=-1.0;
    y_y*=-1.0;
    y_z*=-1.0;
    z_x*=-1.0;
    z_y*=-1.0;
    z_z*=-1.0;
  }
}

EulerRot::EulerRot() {
}

EulerRot::EulerRot(double theta_i,double psi_i,double phi_i) {
	theta=theta_i;
 	psi=psi_i;
 	phi=phi_i;
}

// Rotate the rotation axis by a given rotation matrix
void EulerRot::rotate_axis(const RotMat &rotmat) {
	double c1=cos(theta/2.0);
  double s1=sin(theta/2.0);
	double c2=cos(psi/2.0);
	double s2=sin(psi/2.0);
	double c3=cos(phi/2.0);
	double s3=sin(phi/2.0);
	double c1c2=c1*c2;
	double s1s2=s1*s2;
  double w=c1c2*c3 - s1s2*s3;

  vectorPt axis;
	axis.x=c1c2*s3 + s1s2*c3;
	axis.y=s1*c2*c3 + c1*s2*s3;
	axis.z=c1*s2*c3 - s1*c2*s3;

  double angle=2.0*acos(w);
	double norm = axis.norm();
	if (norm < 0.00000001)
    return; // No rotation
  axis/=norm;

  axis=rotmat*axis;

  double s=sin(angle);
	double c=cos(angle);
	double t=1.0-c;
	if ((axis.x*axis.y*t + axis.z*s) > 0.998) { // north pole singularity detected
		theta = 2.0*atan2(axis.x*sin(angle/2.0),cos(angle/2.0));
		psi = PI/2.0;
		phi = 0;
		return;
	}
	if ((axis.x*axis.y*t + axis.z*s) < -0.998) { // south pole singularity detected
		theta = -2.0*atan2(axis.x*sin(angle/2.0),cos(angle/2.0));
		psi = -PI/2.0;
		phi = 0;
		return;
	}

  theta = atan2(axis.y*s-axis.x*axis.z*t, 1-(axis.y*axis.y+axis.z*axis.z)*t);
	psi = asin(axis.x*axis.y*t+axis.z*s) ;
	phi = atan2(axis.x*s-axis.y*axis.z*t, 1-(axis.x*axis.x+axis.z*axis.z)*t);
}


void EulerRot::operator= (const EulerRot &two) {
	theta=two.theta;
	psi=two.psi;
	phi=two.phi;
}

ostream & operator<<(ostream &out, const EulerRot &rot) {
	out << rot.theta << " " << rot.psi << " " << rot.phi;
	return out;
}

// ------------------- Global functions

// Find the point on the line in the direction of v closest to an outside point
//   The point v_point is known to lie on the line
cPt c::point_to_line(const vectorPt &v, const cPt &v_point, const cPt &point){
	vectorPt u=point-v_point;
	double uTv=u.dot(v);
	double vTv=v.dot(v);
	double qT=uTv/vTv;

	return ( (v*qT)+v_point);
}

// Return the closest distance between two line segments input as end points
double c::closest_approach(const cPt &l1_1, const cPt &l1_2, const cPt &l2_1,
													 const cPt &l2_2) {
  vectorPt u = l1_2 - l1_1;
  vectorPt v = l2_2 - l2_1;
  vectorPt w = l1_1-l2_1;
  double a = u.dot(u);
  double b = u.dot(v);
  double c = v.dot(v);
  double d = u.dot(w);
  double e = v.dot(w);
  double D = a*c - b*b;       
  double sc, sN, sD = D;      
  double tc, tN, tD = D;      

  // compute the closest points between the two lines
  if (D < 0.00000000001) { // parrallel lines
    sN = 0.0; sD = 1.0; tN = e; tD = c;
  } else {  // get the closest points on the infinite lines
    sN = (b*e - c*d);
    tN = (a*e - b*d);
    if (sN < 0.0) {       
      sN = 0.0; tN = e; tD = c;
    } else if (sN > sD) {  
			sN = sD; tN = e + b; tD = c;
		}
  }

  if (tN < 0.0) {
    tN = 0.0;
    if (-d < 0.0)
      sN = 0.0;
    else if (-d > a)
      sN = sD;
    else {
      sN = -d; sD = a;
    }
  } else if (tN > tD) {   
    tN = tD;
    if ((-d + b) < 0.0)
      sN = 0;
    else if ((-d + b) > a)
      sN = sD;
    else {
      sN = (-d + b); sD = a;
    }
  }

	sc=sN/sD;
	tc=tN/tD;

  // vector connecting points
  vectorPt dP=w+(u*sc)-(v*tc); 

  return dP.norm();
}

void c::closest_approach_points(const cPt &l1_1, const cPt &l1_2,
																const cPt &l2_1, const cPt &l2_2,
																cPt &close_l1, cPt &close_l2) {
  vectorPt u = l1_2 - l1_1;
  vectorPt v = l2_2 - l2_1;
  vectorPt w = l1_1-l2_1;
  double a = u.dot(u);
  double b = u.dot(v);
  double c = v.dot(v);
  double d = u.dot(w);
  double e = v.dot(w);
  double D = a*c - b*b;
  double sc, sN, sD = D;
  double tc, tN, tD = D;

  // compute the closest points between the two lines
  if (D < 0.00000000001) { // parrallel lines
    sN = 0.0; sD = 1.0; tN = e; tD = c;
  } else {  // get the closest points on the infinite lines
    sN = (b*e - c*d);
    tN = (a*e - b*d);
    if (sN < 0.0) {
      sN = 0.0; tN = e; tD = c;
    } else if (sN > sD) {
			sN = sD; tN = e + b; tD = c;
		}
  }

  if (tN < 0.0) {
    tN = 0.0;
    if (-d < 0.0)
      sN = 0.0;
    else if (-d > a)
      sN = sD;
    else {
      sN = -d; sD = a;
    }
  } else if (tN > tD) {
    tN = tD;
    if ((-d + b) < 0.0)
      sN = 0;
    else if ((-d + b) > a)
      sN = sD;
    else {
      sN = (-d + b); sD = a;
    }
  }

	sc=sN/sD;
	tc=tN/tD;

  close_l1=u*sc+l1_1;
  close_l2=v*tc+l2_1;

  return;
}

// Returns true if two line segments intersect
bool c::intersect(const c2DPt &line1_start, const c2DPt &line1_end,
									const c2DPt &line2_start, const c2DPt &line2_end) {
  // Check segment1, line2 intersection
	if (!sline_intersect(line1_start,line1_end,(line2_end-line2_start).normal(),
											 line2_start))
	  return false;

  // Check segment2, line1 intersection
	if (!sline_intersect(line2_start,line2_end,(line1_end-line1_start).normal(),
											 line1_start))
	  return false;

	return true;
}  
							 
// Liang-Barsky Intersection between line segment and line
bool c::sline_intersect(const c2DPt &line1_start, const c2DPt &line1_end,
												const c2DPt &line2normal, const c2DPt &line2_point) {
	double denom=(line1_end-line1_start).dot(line2normal);
	if (fabs(denom)<am::epsilon(denom))
		return false; // Parallel lines
	double numerator=(line1_start-line2_point).dot(line2normal)*-1.0;
	double t=numerator/denom;
	if (t<0 || t>1)
		return false;
	return true;
}

/// Average position
cPt c::mean(vector<cPt> &vec) {
  cPt average(vec[0]);
  for (unsigned i=1; i<vec.size(); i++) 
    average+=vec[i];
  return average/vec.size();
}

template class TwoD<double>;
template class TwoD<float>;
template class TwoD<unsigned>;
template ostream & operator<< <double>(ostream &out, const TwoD<double> &t);
template ostream & operator<< <float>(ostream &out, const TwoD<float> &t);
template istream & operator>> <double>(istream &in, TwoD<double> &t);
template istream & operator>> <float>(istream &in, TwoD<float> &t);
template class ThreeD<int>;
template class ThreeD<float>;
template class ThreeD<double>;
template class ThreeD<unsigned>;
template ostream & operator<< <float>(ostream &out, const ThreeD<float> &t);
template istream & operator>> <float>(istream &in, ThreeD<float> &t);
template ostream & operator<< <double>(ostream &out, const ThreeD<double> &t);
template istream & operator>> <double>(istream &in, ThreeD<double> &t);
template ostream & operator<< <unsigned>(ostream &out, const ThreeD<unsigned> &t);
template istream & operator>> <unsigned>(istream &in, ThreeD<unsigned> &t);
template ThreeD<unsigned> operator+ (const unsigned one,
																		 const ThreeD<unsigned> &two);
template ThreeD<double> operator+(const double one, const ThreeD<double> &two);
template ThreeD<unsigned> operator- (const unsigned one,
																		 const ThreeD<unsigned> &two);
template ThreeD<double> operator-(const double one, const ThreeD<double> &two);
template ThreeD<unsigned> operator* (const unsigned one,
																		 const ThreeD<unsigned> &two);
template ThreeD<double> operator*(const double one, const ThreeD<double> &two);
template ThreeD<unsigned> operator/ (const unsigned one,
																		 const ThreeD<unsigned> &two);
template ThreeD<double> operator/(const double one, const ThreeD<double> &two);


