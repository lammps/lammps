/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// machine-specific C++ -> Fortran calling syntax
// It defines the FORTRAN macro for converting variable and function
// names from FORTRAN to C. Different compilers do this in different
// ways. The default is add an underscore to the lower case string.
// Other definitions of the macro can be invoked by defining the
// corresponding macro at compile time using -D e.g. -D_IBM

// CONS(a,b) should return ab, the concatenation of its arguments.
// If compiler is using strict ISO C standards, the ## works.
// Otherwise try the old /**/ trick and test.
// If that fails, you will need to figure out
// a definition for the FORTRAN macro that works on your machine.

#if  __STDC__
#define CONS(a,b) a##b
#elif defined(_IBM)
#define CONS(a,b) a##b
#else
#define CONS(a,b) a/**/b
#warning "The following declaration is a test of the CONS macro"
#warning "If it fails, pair_reax_fortran.h must be modified by hand"
static int my_apples_my_oranges = 1;
static int my_applesoroanges = CONS(my_apples,_my_oranges);
#endif


#ifdef _IBM
#define FORTRAN(lcname,ucname)  lcname
#endif

#ifdef _F2C_LINUX
#define FORTRAN(lcname,ucname)  CONS(lcname,__)
#endif

#ifndef FORTRAN
#define FORTRAN(lcname,ucname)  CONS(lcname,_)
#endif

// hard-wired array sizes set in Fortran library
// accesses include file from Fortran library

#include "reax_defs.h"

class ReaxParams {
 public:
  enum {nneighmax=NNEIGHMAXDEF,
        nat=NATDEF,
        nattot=NATTOTDEF,
        nsort=NSORTDEF,
        mbond=MBONDDEF,
        nbomax=NBOMAXDEF,
  };
};

// data structures corresponding to values in Fortran library

extern "C" struct {
  double abo[ReaxParams::nat];
} FORTRAN(cbkabo,CBKABO);

extern "C" struct {
  double bo[ReaxParams::nbomax];
} FORTRAN(cbkbo,CBKBO);

extern "C" struct {
  double c[3*ReaxParams::nat]; double cglobal[3*ReaxParams::nattot];
  int itag[ReaxParams::nat];
} FORTRAN(cbkc,CBKC);

extern "C" struct {
double ch[ReaxParams::nat];
} FORTRAN(cbkch,CBKCH);

extern "C" struct {
  double chi[ReaxParams::nsort];
  double eta[ReaxParams::nsort];
  double gam[ReaxParams::nsort];
} FORTRAN(cbkchb,CBKCHB);

extern "C" struct {
  double d[3*ReaxParams::nat]; double estrain[ReaxParams::nat];
} FORTRAN(cbkd,CBKD);

extern "C" struct {
  double atomvirial[6*ReaxParams::nat];
  double virial[6];
  int Lvirial;
  int Latomvirial;
} FORTRAN(cbkvirial,CBKVIRIAL);

extern "C" struct {
  int ia[ReaxParams::nat*(ReaxParams::mbond+3)];
  int iag[ReaxParams::nat*(ReaxParams::mbond+3)];
} FORTRAN(cbkia,CBKIA);

extern "C" struct {
  double vlp[ReaxParams::nat];
  double dvlpdsbo[ReaxParams::nat];
} FORTRAN(cbklonpar,CBKLONPAR);

extern "C" struct {
  int nubon1[ReaxParams::nat*(ReaxParams::mbond)];
  int nubon2[ReaxParams::nat*(ReaxParams::mbond)];
} FORTRAN(cbknubon2,CBKNUBON2);

extern "C" struct {
  int nvl1[ReaxParams::nneighmax * ReaxParams::nat];
  int nvl2[ReaxParams::nneighmax * ReaxParams::nat];
  int nvpair;
  int nvlself;
} FORTRAN(cbkpairs,CBKPAIRS);

extern "C" struct {
  int nvlbo[ReaxParams::nneighmax * ReaxParams::nat];
} FORTRAN(cbknvlbo,CBKNVLBO);

extern "C" struct {
  int nvlown[ReaxParams::nneighmax * ReaxParams::nat];
} FORTRAN(cbknvlown,CBKNVLOWN);

extern "C" struct {
  char qa[20*ReaxParams::nattot+10];
} FORTRAN(cbkqa,CBKQA);

extern "C" struct {
  double eb;
  double eoop;
  double epen;
  double estrc;
  double deda[3];
  double pressu;
  double efi;
  double elp;
  double emol;
  double ea;
  double eres;
  double et;
  double eradbo;
  double ev;
  double eco;
  double ecoa;
  double ehb;
  double sw;
  double ew;
  double ep;
  double ekin;
} FORTRAN(cbkenergies,CBKENERGIES);

extern "C" struct {
  double tset;
  double dseed;
  double tempmd;
  double ts2;
  double ts22;
  int nmolo;
  int nmolo5;
  int nbon;
  int na;
  int namov;
  int na_local;
} FORTRAN(rsmall,RSMALL);

// external routines provided by Fortran library

extern "C" void FORTRAN(readc,READC)();
extern "C" void FORTRAN(reaxinit,REAXINIT)();
extern "C" void FORTRAN(ffinpt,FFINPT)();
extern "C" void FORTRAN(tap7th,TAP7TH)();
extern "C" void FORTRAN(taper,TAPER)(double*,double*);
extern "C" void FORTRAN(readgeo,READGEO)();
extern "C" void FORTRAN(srtatom,SRTATOM)();
extern "C" void FORTRAN(vlist,VLIST) ();
extern "C" void FORTRAN(srtbon1,SRTBON1)(int*,int*,double*,int*,int*);
extern "C" void FORTRAN(molec,MOLEC)();
extern "C" void FORTRAN(encalc,ENCALC)();
extern "C" void FORTRAN(getswb,GETSWB)(double*);
extern "C" void FORTRAN(getswa,GETSWA)(double*);
extern "C" void FORTRAN(getvrange,GET_VRANGE)(double*);
extern "C" void FORTRAN(getnvlist,GET_NVLIST)(int*);
extern "C" void FORTRAN(getvlbora,GETVLBORA)(double*);
extern "C" void FORTRAN(cgsolve,CGSOLVE)
  (int*,double*,int*,double*,double*,int*);
extern "C" void FORTRAN(getnval,GETNVAL)(int*);
extern "C" void FORTRAN(getntor,GETNTOR)(int*);
extern "C" void FORTRAN(getnhb,GETNHB)(int*);
extern "C" void FORTRAN(getnbonall,GETNBONALL)(int*);
extern "C" void FORTRAN(getnneighmax,GETNNEIGHMAX)(int*);
extern "C" void FORTRAN(getnat,GETNAT)(int*);
extern "C" void FORTRAN(getnattot,GETNATTOT)(int*);
extern "C" void FORTRAN(getnsort,GETNSORT)(int*);
extern "C" void FORTRAN(getmbond,GETMBOND)(int*);
extern "C" void FORTRAN(getnso,GETNSO)(int*);
extern "C" void FORTRAN(setngeofor,SETNGEOFOR)(int*);
extern "C" void FORTRAN(mdsav,MDSAV)(int*);
extern "C" void FORTRAN(getnsbmax,GETNSBMAX)(int*);
extern "C" void FORTRAN(getnsbma2,GETNSBMA2)(int*);
extern "C" void FORTRAN(getcutof3,GETCUTOF3)(double*);
