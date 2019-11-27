/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: GridMD, ivutils
 *
 *****************************************************************************/
  
/*s****************************************************************************
 * $Log: wpmd.h,v $
 * Revision 1.4  2011/06/11 16:53:55  valuev
 * sync with LAMMPS
 *
 * Revision 1.3  2011/06/11 16:45:23  morozov
 * Fixed erf.c for Windows and Unix
 *
 * Revision 1.2  2011/06/10 19:25:17  morozov
 * *** empty log message ***
 *
 * Revision 1.43  2011/06/10 19:20:53  valuev
 * fixes
 *
 * Revision 1.42  2011/06/09 22:55:08  valuev
 * norm matrices
 *
 * Revision 1.41  2011/06/07 19:58:42  valuev
 * corrected partitions
 *
 * Revision 1.40  2011/06/07 17:43:00  valuev
 * added Y derivatives
 *
 * Revision 1.39  2011/06/03 08:13:33  valuev
 * added partitions to account for ghost atoms
 *
 * Revision 1.38  2011/06/02 22:11:17  morozov
 * Compatibility with LAPACK library
 *
 * Revision 1.37  2011/06/01 23:45:35  valuev
 * modified for LAMMPS compatibility
 *
 * Revision 1.36  2011/05/28 17:16:22  valuev
 * fixed template<>, some fixes to UHF
 *
 * Revision 1.35  2011/05/25 21:30:48  morozov
 * Compatibility with ICC 11.1
 *
 * Revision 1.34  2011/05/25 05:23:43  valuev
 * fixed variable transformation for norm matrix
 *
 * Revision 1.33  2011/05/24 19:54:32  valuev
 * fixed sqmatrix::iterator
 *
 * Revision 1.32  2011/05/21 23:06:49  valuev
 * Norm matrix transform to pysical variables
 *
 * Revision 1.31  2011/05/20 21:39:49  valuev
 * separated norm calculation
 *
 * Revision 1.30  2011/05/14 18:56:19  valuev
 * derivative for ee split interactions
 *
 * Revision 1.29  2011/05/05 08:56:02  valuev
 * working split wp version
 *
 * Revision 1.28  2011/05/04 16:48:52  valuev
 * fixed syntax
 *
 * Revision 1.27  2011/05/04 09:04:48  valuev
 * completed wp_split (except for ee forces)
 *
 * Revision 1.26  2011/04/29 03:07:20  valuev
 * new split wp features
 *
 * Revision 1.25  2011/04/22 09:52:49  valuev
 * working on split WP
 *
 * Revision 1.24  2011/04/20 08:43:09  valuev
 * started adding split packet WPMD
 *
 * Revision 1.23  2010/09/03 12:17:48  morozov
 * The order of parameters in Norm matrix is changed to mimic the order of the single WP parameter storage
 *
 * Revision 1.22  2009/08/27 00:01:36  morozov
 * First working MPI equilibration
 *
 * Revision 1.21  2009/04/14 14:44:10  valuev
 * fixed momentum calculation in hartree model, added "fix" constraint and model="hartree"
 *  to parameters
 *
 * Revision 1.20  2009/04/13 17:00:45  morozov
 * Fixed norm-matrix ratio in AWPMC algorithm
 *
 * Revision 1.19  2009/04/06 17:00:28  morozov
 * Fixed Hartree version of WPMC
 *
 * Revision 1.18  2009/04/01 10:06:37  valuev
 * added Hartee factorization to AWPMD
 *
 * Revision 1.17  2009/03/24 16:10:05  morozov
 * Fixed errors in Norm-matrix calculation related to PBC
 *
 * Revision 1.16  2009/03/17 11:40:04  morozov
 * The prefactor of NormMatrix is corrected
 *
 * Revision 1.15  2008/07/23 16:42:12  valuev
 * Added AWPMD Monte-Carlo
 *
 * Revision 1.14  2008/07/23 15:58:32  valuev
 * *** empty log message ***
 *
 * Revision 1.13  2008/07/21 02:23:22  morozov
 * *** empty log message ***
 *
 * Revision 1.12  2008/07/18 18:15:31  morozov
 * *** empty log message ***
 *
 * Revision 1.11  2008/05/29 13:33:05  valuev
 * VASP band structure
 *
 * Revision 1.10  2008/05/14 17:17:26  morozov
 * Passed 2- and 3-electron test. Added Norm matrix.
 *
 * Revision 1.9  2008/05/05 17:29:32  morozov
 * Fixed errors with Hermitian matrix indeces. Redesigned cVector_3 operations.
 *
 * Revision 1.8  2008/04/28 22:16:45  valuev
 * restructured coulomb term
 *
 * Revision 1.7  2008/04/28 09:54:13  valuev
 * corrected summation for Eee part
 *
*******************************************************************************/
# ifndef WPMD_H
# define WPMD_H
  
/** @file wpmd.h 
    @brief Classes for Wave Packet Molecular Dynamics of two component plasma. */

# ifndef _USE_MATH_DEFINES
# define _USE_MATH_DEFINES
# endif
# include <complex>
# include <vector>
# include <cmath> 
# include "logexc.h"
# include "cvector_3.h"
# include "pairhash.h"
# include "TCP/tcpdefs.h"
# include "wavepacket.h"
# include "erf.h"
# include "cerf.h"


using namespace std;

# include "lapack_inter.h"

typedef hmatrix<cdouble> chmatrix;

const cdouble i_unit = cdouble(0.,1.), i_unit1 = -i_unit;
const double h_plank2 = h_plank * 2.;

//cdouble ccerf(const cdouble &a);

//e calculates cerf(c*z)/z
inline cdouble cerf_div(const cdouble &z, const cdouble &c=i_unit){
  if((fabs(real(z))+fabs(imag(z)))<1e-8)
    return c*two_over_sqr_pi;
  else 
    return cerf(z*c)/z;
}

//e calculates erf(c*z)/z
inline double erf_div(const double &z, double c=1){
  if(fabs(z)<1e-8)
    return c*two_over_sqr_pi;
  else 
    return erf(z*c)/z;
}

template<class T1, class T2>
struct _myrefpair{
  T1& first;
  T2& second;
  _myrefpair(T1& a, T2 &b):first(a),second(b){}
  _myrefpair operator=(const pair<T1,T2> &other){
    first=other.first;
    second=other.second;
    return *this;
  }
};

template< class T1, class T2>
_myrefpair<T1,T2> _mytie(T1& var1, T2 &var2){
  return _myrefpair<T1,T2>(var1, var2);
}

inline pair<double,double> operator*(const pair<double,double> &right, double left){
  return make_pair(right.first*left,right.second*left);
}

// Auxilary class to handle the normalizing term derivatives
class NormDeriv
{
public:
  cdouble l;    // lambda = (f over a_re)
  double m;     // mu = (f over a_im) / i
  cVector_3 u;  // u = (f over b_re)
  Vector_3 v;   // v = (f over b_im) / i

  NormDeriv() {}
  NormDeriv(const WavePacket& wp) { set(wp); }

  //e Create NormDeriv object and calculate the derivatived for the given WP
  void set(const WavePacket& wp){
    Vector_3 br = real(wp.b), bi = imag(wp.b);
    double ar = real(wp.a), ai = imag(wp.a);
    double i_2ar = 0.5 / ar, ai_ar = ai / ar;

    v = (-i_2ar) * br;
    m = v.norm2();
    u = v * (i_unit1 * ai_ar) - wp.b * i_2ar;
    l = (1.5*i_2ar + m) + cdouble(0.,2.) * ( (br*bi)*i_2ar*i_2ar - m*ai_ar );
  }
};

inline NormDeriv conj(const NormDeriv& src){
  NormDeriv dst;
  dst.l = conj(src.l);
  dst.m = -src.m;
  dst.u = conj(src.u);
  dst.v = - src.v;
  return dst;
}

///\en Auxilary class to handle derivatives of overlaps
class OverlapDeriv{
public:
  WavePacket w1, w2, w12;
  NormDeriv d1, d2;
  cdouble I0, I2;
  cVector_3 I1;
  cdouble bb_4a;
  sqmatrix<cdouble> IDD;


  OverlapDeriv():I0(0),I1(0),IDD(10){}
  
  void set1(const WavePacket& w1_) {
    w1=w1_;
    d1.set(w1);
    d1=conj(d1);
  }

  //e Create NormDeriv object and calculate the derivatived for the given WP
  void set2(const WavePacket& w2_, const cdouble *I0_=NULL){
    w2=w2_;
    d2.set(w2);
    w12=conj(w1)*w2;
    if(!I0_)
      I0 = w12.integral();
    else
      I0=*I0_;
    I1 = w12.b * (I0 / w12.a / 2);
    bb_4a = w12.b.norm2() / w12.a / 4;
    I2 = I0 * (bb_4a + 1.5) / w12.a;
  }

  cdouble da2_re() const {
    return (d2.l*I0 - I2);
  }

  cdouble da2_im() const {
    return i_unit*(d2.m*I0 - I2);
  }

  cdouble da1_re() const {
    return (d1.l*I0 - I2);
  }

  cdouble da1_im() const {
    return -i_unit*(d1.m*I0 - I2);
  }

  cdouble db2_re(int i) const {
    return d2.u[i]*I0 + I1[i];
  }

  cdouble db2_im(int i) const {
    return i_unit*(d2.v[i]*I0 + I1[i]);
  }

  cdouble db1_re(int i) const {
    return d1.u[i]*I0 + I1[i];
  }

  cdouble db1_im(int i) const {
    return -i_unit*(d1.v[i]*I0 + I1[i]);
  }

  ///\en Calculates  derivative overlap matrix IDD
  void calc_der_overlap(bool self=false, cdouble cc1=0., cdouble c2=0.);


};


class AWPMD {
//protected:
public:
  int ne[2], ni;
  int nwp[2]; ///<\en number of wavepackets (same as ne for unsplit version)
  int nvar[2]; ///<\en full number of dynamic variables for each spin
  chmatrix O[2], Y[2], Te[2], Tei[2];
  smatrix<unsigned char> Oflg[2];  ///<\en equals 0 for non-overlaping packets
  sqmatrix<double> Norm[2];  ///<\en Norm matrix
  vector<WavePacket> wp[2]; ///<\en wave packets for electrons (spin up =0 and down =1)
  vector<double> qe[2];  ///<\en electron charges
  vector<double> qi;  ///<\en ion charges
  vector<Vector_3> xi;  ///<\en ion coordinates
  int pbc; ///<\en pbc flag
  Vector_3 cell; ///<\en cell coordinates (L1,L2,L3)
  double Lextra; ///<\en width PBC, unset if negative
  double harm_w0_4;
  double w0;
  int calc_ii; ///<\en flag indicating whether to calculate ion-ion interaction
  int norm_needed; ///<\en flag indicating whether to prepare norm matrix data in interaction function

  enum {NORM_UNDEFINED, NORM_CALCULATED, NORM_FACTORIZED, NORM_INVERTED};
  int norm_matrix_state[2];
  
  // Arrays for temporal data
  chmatrix IDD;  // Second derivatives of the overlap integral (used in Norm matrix)
  vector<cdouble> ID, IDYs;  // First derivatives of the overlap integral (used in Norm matrix)
  vector<int> ipiv;  // The pivot indices array

  recmatrix<cdouble> L[2]; ///<\en overlap derivative matrix for each spin
  recmatrix<cdouble> M[2]; ///<\en 'normalized' overlap derivative matrix for each spin: M=YL

public:
  enum {NONE=0, HARM, FIX, RELAX} constraint;

  ///\em Sets approximation level for quantum problem: \n
  ///    HARTREE Hartree product (no antisymmetrization) \n
  ///    DPRODUCT product of det0*det1 of antisymmetrized functions for spins 0, 1 \n
  ///    UHF unrestricted Hartree-Fock
  enum APPROX {HARTREE, DPRODUCT,  UHF } approx; 
  ///\em Sets overlap matrix element to zero if the overlap norm is less than this value
  double ovl_tolerance;

  double Ee[2]; ///<\en electron kinetic energy for each spin
  double Eei[2]; ///<\en electron-ion energy for each spin
  double Eee, Ew; ///<\en electron-electron energy
  double Eii; ///<\en ion-ion energy
  double Edk; ///<\en sum of diagonal kinetic energy terms
  double Edc; ///<\en sum of diagonal Coulomb energy terms

  vector<double> Eep[2]; ///<\en per particle electron kinetic energy for each spin
  vector<double> Eeip[2]; ///<\en per particle electron-ion energy for each spin
  vector<double> Eeep[2]; ///<\en per particle electron-electron energy for each spin
  vector<double> Ewp[2]; ///<\en per particle restrain energy for each spin
  vector<double> Eiep; ///<\en per particle ion-electron energy
  vector<double> Eiip; ///<\en per particle ion-ion energy


  ///\en \{ Conversion constants that depend on the unit system used (for LAMMPS compatibility).
  ///       Default is GRIDMD units. Change them according to your unit system.
  double me; ///<\en electron mass (LAMMPS: me in the appropriate unit system)
  double one_h; ///<\en inverse of Plancks constant (LAMMPS: conversion [(m*v)/h] to [distance]  )
  double h2_me; ///<\en Plancks constant squared divided by electron mass (LAMMPS: conversion [h^2/(m*r^2)] to [Energy]  )
  double coul_pref; ///<\en Coulomb prefactor (e2 for GRIDMD) (LAMMPS: conversion [q^2/r] to [Energy]  )
  ///    \}

  ///\en 0 -- indicates that the inter-partition force should be full, and energy half,\n
  ///    1 -- inter-partition force and energy counts one half (LAMMPS compatibility)
  int newton_pair; 

  //int myid; ///<\en id for partitions
  
  ///\en Partition arrays storing the tags of particles. The initial tags should be >0.
  ///    If the tag stored is <0, then the particle is ghost with -tag.
  ///    partition1[2] is for ions, 0, 1 for each electron spin
  vector<int> partition1[3]; 
  //vector<int> partition2[3]; ///<\en 2 for ions


  int tag_index(int i, int j){
    return i==j ? -1 : (i>j ? (i-2)*(i-1)/2+j : (j-2)*(j-1)/2+i );
  }


  ///\en 1 -- all my, -1 all other, 2 -- my mixed term, -2 -- other mixed term 
  int check_ee(int s1,int icj1,int ick2){
    //printf(" (%d %d) ",partition1[s1][icj1],partition1[s1][ick2]);
    int c1=(int)(partition1[s1][icj1]>0);
    int c2=(int)(partition1[s1][ick2]>0);
    int res;
    if(c1!=c2){ // mixed
      int tag1=abs(partition1[s1][icj1]);
      int tag2=abs(partition1[s1][ick2]);
      int num=tag_index(tag1-1,tag2-1);
      if(num<0){ // compare wave packets
        int cmp= s1<2 ? 
          wp[s1][icj1].compare(wp[s1][ick2],1e-15) : 
          compare_vec(xi[icj1],xi[ick2],1e-15);
        if((cmp>0 && c1) || (cmp<0 && c2))
          res= 2; // my mixed term
        else
          res= -2; // not my term
      }
      else // parity check
        res=num%2 ? 2 : -2;
    }
    else if(c1)
      res=1; // all my
    else
      res=-1; // all other
    return res;
  }

  ///\en Returns electron-electron inter-partition multipliers for energy (first) and force (second)
  ///    for a 4- and 2- electron additive terms (all inter-partition interactions are 
  ///    calculated only once based on particle tags)
  ///    If force multiplier is zero, then the term may be omitted (energy will also be zero).
  ///    NOW ASSIGNS BASED ON THE FIRST PAIR ONLY
  pair<double, double> check_part1(int s1,int icj1,int ick2){
    int res=check_ee(s1,icj1,ick2);    
    if(res==1){ // my term
      //printf(" *\n");
      return make_pair(1.,1.); // all at my partition
    }
    else if(res==-1){
      //printf(" \n");
      return make_pair(0.,0.); // all at other partition
    }
    else if(res==2){
      //printf(" *\n");
      return make_pair(1.,1.); // my inter-partition
    }
    else if(res==-2){ 
      //printf(" \n");
      return make_pair(0., newton_pair ? 0.0 : 1. ); // other inter-partition: must add force if newton comm is off
    }
    return make_pair(0.,0.); // nonsense
  }

  ///\en Returns elctron-ion inter-partition multipliers for energy (first) and force (second)
  ///    for ion-electron additive terms (all inter-partition interactions are 
  ///    calculated only once based on particle tags)
  ///    If force multiplier is zero, then the term may be omitted (energy will also be zero).
  ///    BASED ON ION ATTACHMENT
  pair<double,double> check_part1ei(int s1,int icj1,int ick2, int ion){
    //printf("%d ",partition1[2][ion]);
    int ci=(int)(partition1[2][ion]>0);
    
    if(!newton_pair){ // care about mixed terms
      int cee=check_ee(s1,icj1,ick2);
      if((cee==2 || cee==-2) || (ci && cee==-1) || (!ci && cee==1)) // all mixed variants
        make_pair(0., 1. ); // other inter-partition: must add force if newton comm is off 
    }
    if(ci){
      //printf(" *\n");
      return make_pair(1.,1.); // my term
    }
    else{
      //printf(" \n");
      return make_pair(0.,0.); // all at other partition
    }
  }

  ///\en Returns ion-ion inter-partition multipliers for energy (first) and force (second)
  ///    for ion-electron additive terms (all inter-partition interactions are 
  ///    calculated only once based on particle tags)
  ///    If force multiplier is zero, then the term may be omitted (energy will also be zero).
  pair<double,double> check_part1ii(int ion1, int ion2){
    return check_part1(2,ion1,ion2);
  }



  AWPMD():pbc(0),Lextra(-1),constraint(NONE),newton_pair(1) {
    nvar[0]=nvar[1]=ne[0]=ne[1]=ni=0;
    norm_matrix_state[0] = norm_matrix_state[1] = NORM_UNDEFINED;
    ovl_tolerance=0.;
    approx = DPRODUCT;
    
    me=m_electron;
    one_h=1./h_plank;
    h2_me=h_sq/me;
    coul_pref=::coul_pref;

    calc_ii=0;
    norm_needed=0;
  }

protected:
  //e translates wp2 to the nearest image postion relative to wp1
  //e gets the translation vector
  Vector_3 move_to_image(const WavePacket &wp1, WavePacket &wp2) const {
    Vector_3 r1=wp1.get_r();
    Vector_3 r2=wp2.get_r();
    Vector_3 dr=r2-r1;
    Vector_3 ndr=dr.rcell(cell,pbc); // [0,L)
    ndr=ndr.rcell1(cell,pbc); // [-L/2,L/2)
    ndr-=dr;
    wp2=wp2.translate(ndr);  // wln.b=wln.b+ndr.a
    return ndr;
  }

  //e gets the overlap packet taking PBC into account
  WavePacket pbc_mul(const WavePacket &wp1, const WavePacket &wp2) const {
    if(!pbc)
      return wp1*conj(wp2);
    Vector_3 r1=wp1.get_r();
    Vector_3 r2=wp2.get_r();
    Vector_3 dr=r2-r1; // distance
    Vector_3 drn=dr.rcell1(cell,pbc); // distance within PBC
    Vector_3 rtrans=drn-dr; // new location of wp2 according to PBC (nearest image)
    WavePacket wpn=wp2.translate(rtrans);
    wpn=wp1*(conj(wpn));
    // reducing the result to elementary cell
    //r1=wpn.get_r();
    //r2=r1.rcell(cell,pbc);
    //dr=r2-r1;
    //wpn=wpn.translate(dr);
    return wpn;
  }
  ///\en resizes all internal arrays according to new electrons added
  virtual void resize(int flag);
public:
  
  
  
  ///\en Prepares to setup a new system of particles using \ref add_ion() and add_electron().
  ///    There is no need to call this function when using 
  ///    \ref set_electrons() and \ref set_ions() to setup particles.
  virtual void reset(){
    for(int s=0;s<2;s++){
      nwp[s]=ne[s]=nvar[s]=0;
      wp[s].clear();
      qe[s].clear();
      partition1[s].clear();
      //partition2[s].clear();
    }
    partition1[2].clear();
    ni=0;
    xi.clear();
    qi.clear();
  }

  //e sets Periodic Boundary Conditions
  //e using bit flags: 0x1 -- PBC along X
  //e                  0x2 -- PBC along Y
  //e                  0x4 -- PBC along Z
  //e cell specifies the lengths of the simulation box in all directions
  //e if PBCs are used, the corresponding coordinates of electrons and ions
  //e in periodic directions must be within the range  [0, cell[per_dir]) 
  //e @returns 1 if OK
  int set_pbc(const Vector_3P pcell=NULL, int pbc_=0x7);
 
  ///\en Setup electrons: forms internal wave packet representations.
  ///    If PBCs are used the coords must be within a range [0, cell).
  ///    Default electron mass is AWPMD::me.
  ///    Default (q=NULL )electron charges are -1.
  int set_electrons(int spin, int n, Vector_3P x, Vector_3P v, double* w, double* pw, double mass=-1, double *q=NULL);

  //e setup ion charges and coordinates
  //e if PBCs are used the coords must be within a range [0, cell)
  int set_ions(int n, double* q, Vector_3P x);

  ///\en Adds an ion with charge q and position x,
  ///    \return id of the ion starting from 0
  ///    The tags must be nonzero, >0 for the local particle, <0 for ghost particle. 
  ///    Unique particle id  is abs(tag).
  ///    Default tag (0) means inserting the current particle id as local particle. 
  int add_ion(double q, const Vector_3 &x, int tag=0){
    qi.push_back(q);
    xi.push_back(x);
    ni=(int)xi.size();
    if(tag==0)
      tag=ni;
    partition1[2].push_back(tag);
    return ni-1;
  }


  //e calculates interaction in the system of ni ions + electrons 
  //e the electonic subsystem must be previously setup by set_electrons, ionic by set_ions
  //e the iterators are describing ionic system only
  // 0x1 -- give back ion forces
  // 0x2 -- add ion forces to the existing set
  // 0x4 -- calculate derivatives for electronic time step (NOT IMPLEMENTED)
  //e if PBCs are used the coords must be within a range [0, cell)
  virtual int interaction(int flag=0, Vector_3P fi=NULL, Vector_3P fe_x=NULL, 
                                      Vector_3P fe_p=NULL, double *fe_w=NULL, double *fe_pw=NULL, Vector_2P fe_c=NULL);

  //e same as interaction, but using Hartee factorization (no antisymmetrization)
  virtual int interaction_hartree(int flag=0, Vector_3P fi=NULL, Vector_3P fe_x=NULL, 
                                      Vector_3P fe_p=NULL, double *fe_w=NULL, double *fe_pw=NULL, Vector_2P fe_c=NULL);

  ///\en Calculates ion-ion interactions and updates Eii and ion forces if requested. This function
  ///    is called form intaraction() and interaction_hartree if calc_ii is set.
  virtual int interaction_ii(int flag,Vector_3P fi=NULL);
  
  //e Calculates Norm matrix
  //e The result is saved in AWPMD::Norm[s]
  void norm_matrix(int s);

  //e Performs LU-factorization of the Norm matrix
  //e AWPMD::Norm[s] is replaced by the LU matrix
  void norm_factorize(int s);

  //e Invert Norm matrix
  //e AWPMD::Norm[s] is replaced by the inverted matrix
  void norm_invert(int s);

  //e Get the determinant of the norm-matrix for the particles with spin s
  double norm_matrix_det(int s);

  //e Get the determinant logarithm of the norm-matrix for the particles with spin s
  double norm_matrix_detl(int s);

  double get_energy();

  //e makes timestep of electronic component (NOT IMPLEMENTED)
  int step(double dt);

  ///\en Gets current electronic coordinates.
  ///    Transforms the momenta to velocity v according to mass setting (-1 means me)
  int get_electrons(int spin, Vector_3P x, Vector_3P v, double* w, double* pw, double mass=-1);

  void set_harm_constr(double w0) {
    constraint = HARM;
    harm_w0_4 = h_sq*9./(8.*m_electron)/(w0*w0*w0*w0);
  }

  void set_fix_constr(double w0_) {
    constraint = FIX;
    w0 = w0_;
  }

  ///\en Prepares force arrays according to \a flag setting for interaction()
  virtual void clear_forces(int flagi,Vector_3P fi, Vector_3P fe_x, 
                    Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c=NULL);


  ///\en Creates wave packet acording to the given physical parameters.
  ///    The function may change its arguments by applying existing constraints!
  ///    Default mass (-1) is the electron mass AWPMD::me.
  WavePacket create_wp(Vector_3 &x, Vector_3 &v, double &w, double &pw, double mass=-1);

  
};


# endif
