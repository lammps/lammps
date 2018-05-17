//HGvector.h
//class listing for svector class

#ifndef HGVECTOR_H
#define HGVECTOR_H 1

//#include <iostream>
//#include <cmath>
//#include "consts.h" 

//using namespace std;


class HGvector
{
 public:
  double x,y,z;

  HGvector() {x=y=z=0;}
  HGvector(double X, double Y, double Z) {x=X; y=Y; z=Z;}
  void clear() {x=y=z=0;}
  bool is_zero() {return (x==0.0 && y==0.0 && z==0.0);}
  double sqmag() {return x*x+y*y+z*z;}
  double mag() {return sqrt(x*x+y*y+z*z);}
  // // double minsqmag(double& Lx, double& Ly)
  // //   {minimg(Lx, Ly); return sqmag();}
  // // double minmag(double& Lx, double& Ly)
  // //   {return sqrt(minsqmag(Lx, Ly));}
  // // void minimg(double& Lx, double& Ly){
  // //   while (abs(x/Lx)>0.5) x = (x<0) ? x+Lx : x-Lx;
  // //   while (abs(y/Ly)>0.5) y = (y<0) ? y+Ly : y-Ly;
  // // }
  void set(double X, double Y, double Z) {x=X; y=Y; z=Z;}
  void operator = (const HGvector& v) {x=v.x; y=v.y; z=v.z;}
  void operator *= (double f) {x*=f; y*=f; z*=f;}
  void operator /= (double f) {x/=f; y/=f; z/=f;}
  void operator += (double f) {x+=f; y+=f; z+=f;}
  void operator -= (double f) {x-=f; y-=f; z-=f;}
  void operator += (const HGvector& v) {x+=v.x; y+=v.y; z+=v.z;}
  void operator -= (const HGvector& v) {x-=v.x; y-=v.y; z-=v.z;}
  friend double operator ^ (const HGvector& V1, const HGvector& V2) //dot
    {return V1.x*V2.x+V1.y*V2.y+V1.z*V2.z;}
  friend HGvector operator * (const HGvector& V1, const HGvector& V2)  //cross
    {return HGvector(V1.y*V2.z-V1.z*V2.y, 
	                   V1.z*V2.x-V1.x*V2.z,
	                   V1.x*V2.y-V1.y*V2.x);}
  friend HGvector operator * (const HGvector& V1, const double f)
    {return HGvector(V1.x*f, V1.y*f, V1.z*f);}
  friend HGvector operator * (const double f, const HGvector& V1)
    {return HGvector(V1.x*f, V1.y*f, V1.z*f);}
  friend HGvector operator / (const HGvector& V1, const double f)
    {return HGvector(V1.x/f, V1.y/f, V1.z/f);}
   friend HGvector operator + (const HGvector& V1, const HGvector& V2)
    {return HGvector(V1.x+V2.x, V1.y+V2.y, V1.z+V2.z);}
  friend HGvector operator - (const HGvector& V1, const HGvector& V2)
    {return HGvector(V1.x-V2.x, V1.y-V2.y, V1.z-V2.z);}
  // void EulerTrans(double, double, double);
  // void Rec2Cyl();
  // void Cyl2Rec();
  // void Rec2Sph();
  // void Sph2Rec();
  // friend void Print(svector, ostream&);
  //friend std::ostream& operator<< (std::ostream&, const svector&); //output
  //friend std::istream& operator>> (std::istream&, svector&); //input
};

#endif
