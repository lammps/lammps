# ifndef WAVEPACKET_H
# define WAVEPACKET_H

# ifndef _USE_MATH_DEFINES
# define _USE_MATH_DEFINES
# endif

# include <cmath>
# include <complex>
# include <functional>
# include "cvector_3.h"

using namespace std;

/** @file wpmd.h 
    @brief Classes to handle Gaussian Wave Packets. */

// Constants
const double MIN_EXP_ARG = -15.;  // Minimum noticeable argument for exp function

class WavePacket;


///\en Template for v=der operation in \ref Wavepacket::int2phys_der()
template<class Type>
struct eq_second : public binary_function <Type, Type, Type> {
  Type operator()(const Type& _Left, const Type& _Right) const{
    return _Right;
  }
};

///\en Template for v=-der operation in \ref Wavepacket::int2phys_der()
template<class Type>
struct eq_minus_second : public binary_function <Type, Type, Type> {
  Type operator()(const Type& _Left, const Type& _Right) const{
    return -_Right;
  }
};

///\en Compares complex numbers on a per component basis.
///    \return  \retval 0 if all component differences are 0 within tolerance \a tol (EQUAL),
///             \retval -1 for LESS
///             \retval 2  for GREATER
template< class CT >
int compare_compl(const CT &a, const CT& b, double tol=0.){
  double dd=real(a)-real(b);
  if(dd<-tol)
    return -1;
  if(dd>tol)
    return 1;
  dd=imag(a)-imag(b);
  if(dd<-tol)
    return -1;
  if(dd>tol)
    return 1;
  return 0;
}


///\en Compares vectors on a per component basis.
///    \return  \retval 0 if all component differences are 0 within tolerance \a tol (EQUAL),
///             \retval -1 for LESS
///             \retval 2  for GREATER
inline int compare_vec(const Vector_3 &a, const Vector_3& b, double tol=0.){
  for(int i=0;i<3;i++){
    double dd=a[i]-b[i];
    if(dd<-tol)
      return -1;
    if(dd>tol)
      return 1;
  }
  return 0;
}


/// wavepacket is w(x)=exp(-a*x^2+b*x+lz)
class WavePacket{
  /// constructs a conjugate packet
  friend WavePacket conj(const WavePacket &wp);

public:
  cdouble a;
  cVector_3 b;
  cdouble lz;

  WavePacket(const cdouble &sa=cdouble(1.,0.),const cVector_3 &sb=cVector_3(0.,0.), const cdouble &slz=0.): a(sa), b(sb), lz(slz){
  }

  WavePacket operator*(const WavePacket& other) const {
    return WavePacket(a+other.a,b+other.b,lz+other.lz);
  }
 
  /// returns the integral of w(x) over 3D space
  cdouble integral() const {
    cdouble z = lz + b.norm2()/(4.*a);
    if(real(z) < MIN_EXP_ARG) return 0.;
    return pow(M_PI/a,3./2.)*exp(z);
  }

  /// init normalized packet with physical parameters: r0, p0, width, pw
  /// w(x)=(3/2pi width^(3/4)exp[-(3/(4 width^2)-i pw/(2*width) )(x-r)^2+i p (x-r)]
  void init(const double width=1., const Vector_3 &r=0., const Vector_3 &p=0.,  const double pw=0.){
    a  = (3./(2.*width) - cdouble(0.,1.)*pw)/(2.*width);
    b  = (2.*a)*r + cdouble(0.,1.)*p;
    lz = log(3./(2.*M_PI*width*width))*(3./4.) + (-a*r.norm2() - cdouble(0.,1.)*(p*r));
  }

  /// init normalized packet with complex parameters a and b
  /// w(x)=(3/2pi width^(3/4)exp[-(3/(4 width^2)-i pw/(2*width) )(x-r)^2+i p (x-r)]
  void init(const cdouble &a_, const cVector_3 &b_){
    a = a_;
    b = b_;
    Vector_3 r = get_r();
    double r2 = r.norm2();
    lz = cdouble( log( real(a)*(2./M_PI) )*(3./4.) - r2*real(a), r2*imag(a) - r*imag(b) );
  }

  cdouble operator()(const Vector_3 &x) const {
    return exp(lz - a*x.norm2() + b*x);
  }

  /// ajusts lz so that Integral[w*(conj(w))] is 1 after this operation
  WavePacket& normalize(){
    //lz=0.;
    //lz=overlap(conj(*this));
    //lz=(-1./2)*log(lz);
    Vector_3 r = get_r();
    double r2 = r.norm2();
    lz = cdouble( log( real(a)*(2./M_PI) )*(3./4.) - r2*real(a), r2*imag(a) - r*imag(b) );
    return *this;
  }

  /// computes 3D overlap of wavepackets Inegral[w*other]
  cdouble overlap(const WavePacket &other) const{
    WavePacket wp=(*this)*other;
    return wp.integral();
  }

  /// returns translated packet to the position of  r'=r+dr
  WavePacket translate(const Vector_3 &dr) const {
    WavePacket wp(a,b,lz);
    wp.b+=2*dr*a;
    Vector_3 r=get_r();
    double dr2=(r+dr).norm2()-r.norm2();
    wp.lz+=-a*dr2-cdouble(0.,(get_p()*dr));
    return wp;
  }

  ///  width
  double get_width() const{
    return sqrt(3./4./real(a));
  }
  /// width momentum
  double get_pwidth() const{
    return -imag(a)*sqrt(3./real(a));
  }
  /// both width and width momentum
  pair<double,double> get_width_pars() const{
    double c=sqrt(3./2./real(a));
    return make_pair(c/2,-imag(a)*c);
  }
  /// position
  Vector_3 get_r() const {
    return real(b)/(2*real(a));
  }
  /// momentum
  Vector_3 get_p() const {
    return imag(b) - real(b)*(imag(a)/real(a));
  }

  ///\en Transforms derivatives of a function whith respect to WP parameters
  ///    from internal into physical representation, i. e.:\n
  ///    from df/d{are,aim,b0re,b0im,b1re,b1im,b2re,b2im} (8 values accessed by input iterator d_it in the given order)\n
  ///    to   df/d{x0,x1,x2}, df/d{p0,p1,p2}, df/dw, df/dpw 
  ///    The supplied inputs (val) are modified by op: val=op(val,phys_der).
  ///    Use operation=eq_second for the supplied inputs to be replaced by new physical derivative values.
  ///    The inpput and output locations may coinside, an internal buffer is used for transformation.
  template<template<class A> class operation, class d_it, class dfdx_it, class dfdp_it, class dfdw_it, class dfdpw_it>
  void int2phys_der(d_it dfdi_,dfdx_it dfdx, dfdp_it dfdp,  dfdw_it dfdw, dfdpw_it dfdpw, double h_p=h_plank) const {
    operation<double> op;
    double dfdi[8], dfn[8];// internal buffer
    for(int i=0;i<8;i++)
      dfdi[i]=*dfdi_++;
    double w=get_width();
    Vector_3 r=get_r();
    double t=3/(2*w*w*w);
    dfn[6]= -t*dfdi[0]-imag(a)*dfdi[1]/w;  //dw
    dfn[7]=-dfdi[1]/(2*w*h_p);  // dpw
    for(int i=0;i<3;i++){
      dfn[i]= 2*real(a)*dfdi[2+2*i]+2*imag(a)*dfdi[2+2*i+1];
      dfn[3+i]= dfdi[2+2*i+1]*(/*m_electron*/1./h_p) ; //*(h_plank/m_electron);
      dfn[7]+=-(r[i]*dfdi[2+2*i+1]/w)/h_p;  
      dfn[6]+=-2*r[i]*(t*dfdi[2+2*i]+imag(a)*dfdi[2+2*i+1]/w);
    }
    int i=0;
    for(int k=0;k<3;k++)
      *dfdx++=op(*dfdx,dfn[i++]);
    for(int k=0;k<3;k++)
      *dfdp++=op(*dfdp,dfn[i++]);
    *dfdw=op(*dfdw,dfn[i++]);
    *dfdpw=op(*dfdpw,dfn[i++]);
  }

  
  ///\en Compares the wave packet to another on a per component basis.
  ///    \return  \retval 0 if all component differences are 0 within tolerance \a tol (EQUAL),
  ///             \retval -1 for LESS
  ///             \retval 2  for GREATER
  int compare(const WavePacket &other, double tol=0.) const {
    int res=compare_compl(a,other.a, tol);
    if(res)
      return res;
    for(int i=0;i<3;i++){
      res=compare_compl(b[i],other.b[i]);
      if(res)
        return res;
    }
    return 0;
  }

};

    /*double w=wk.get_width();
          Vector_3 r=wk.get_r();
          double t=3/(2*w*w*w);
          fe_w[ic1+k1]+= t*E_der[s1][indw1+8*k1]+imag(wk.a)*E_der[s1][indw1+8*k1+1]/w;
          fe_pw[ic1+k1]+=E_der[s1][indw1+8*k1+1]/(2*w*h_plank);
          for(int i=0;i<3;i++){
            fe_x[ic1+k1][i]+= -2*real(wk.a)*E_der[s1][indw1+8*k1+2+2*i]-2*imag(wk.a)*E_der[s1][indw1+8*k1+2+2*i+1];
            fe_p[ic1+k1][i]+= (-E_der[s1][indw1+8*k1+2+2*i+1])*(m_electron/h_plank); //*(h_plank/m_electron);
            fe_pw[ic1+k1]+=(r[i]*E_der[s1][indw1+8*k1+2+2*i+1]/w)/h_plank;  
            fe_w[ic1+k1]+=2*r[i]*(t*E_der[s1][indw1+8*k1+2+2*i]+imag(wk.a)*E_der[s1][indw1+8*k1+2+2*i+1]/w);
          }*/


/// constructs a conjugate packet
inline WavePacket conj(const WavePacket &wp){
  return WavePacket(conj(wp.a),conj(wp.b),conj(wp.lz));
}


# endif  // WAVEPACKET_H