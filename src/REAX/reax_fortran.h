/////:EOH~

/* CONS(a,b) should return ab, the concatenation
   of its arguments */
#if  __STDC__ 
#define CONS(a,b) a##b
#else
#define CONS(a,b) a/**/b
#endif

#ifdef _ARDENT
#define FORTRAN(lcname,ucname)  ucname
#endif

#ifdef _IBM
#define FORTRAN(lcname,ucname)  lcname
#endif

#ifdef _HP
#define FORTRAN(lcname,ucname)  lcname
#endif

#ifdef _F2C_LINUX
#define FORTRAN(lcname,ucname)  CONS(lcname,__)
#endif

#ifndef FORTRAN
#define FORTRAN(lcname,ucname)  CONS(lcname,_)
#endif


