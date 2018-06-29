/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_wrapper_Plumed_h
#define __PLUMED_wrapper_Plumed_h

/**
\page ReferencePlumedH Reference for interfacing MD codes with PLUMED

  Plumed.h and Plumed.c contain the external plumed interface, which is used to
  integrate it with MD engines. This interface is very general, and is expected
  not to change across plumed versions. Plumed.c also implements a dummy version
  of the interface, so as to allow a code to be fully linked even if the plumed
  library is not available yet. These files could be directly included in the official
  host MD distribution. In this manner, it will be sufficient to link the plumed
  library at link time (on all systems) or directly at runtime (on system where
  dynamic loading is enabled) to include plumed features.

  Why is Plumed.c written in C and not C++? The reason is that the resulting Plumed.o
  needs to be linked with the host MD code immediately (whereas the rest of plumed
  could be linked a posteriori). Imagine the MD code is written in FORTRAN: when we
  link the Plumed.o file we would like not to need any C++ library linked. In this
  manner, we do not need to know which C++ compiler will be used to compile plumed.
  The C++ library is only linked to the "rest" of plumed, which actually use it.
  Anyway, Plumed.c is written in such a manner to allow its compilation also in C++
  (C++ is a bit stricter than C; compatibility is checked when PlumedStatic.cpp,
  which basically includes Plumed.c, is compiled with the C++ compiler). This will
  allow e.g. MD codes written in C++ to just incorporate Plumed.c (maybe renamed into
  Plumed.cpp), without the need of configuring a plain C compiler.

  Plumed interface can be used from C, C++ and FORTRAN. Everything concerning plumed
  is hidden inside a single object type, which is described in C by a structure
  (struct \ref plumed), in C++ by a class (PLMD::Plumed) and in FORTRAN by a
  fixed-length string (CHARACTER(LEN=32)). Obviously C++ can use both struct
  and class interfaces, but the first should be preferred. The reference interface
  is the C one, whereas FORTRAN and C++ interfaces are implemented as wrappers
  around it.

  In the C++ interface, all the routines are implemented as methods of PLMD::Plumed.
  In the C and FORTRAN interfaces, all the routines are named plumed_*, to
  avoid potential name clashes. Notice that the entire plumed library
  is implemented in C++, and it is hidden inside the PLMD namespace.
  If the used C++ compiler supports C++11, PLMD::Plumed object defines move semantics
  so as to be usable in STL containers. That is, you can declare a std::vector<PLMD::Plumed>.

  Handlers to the plumed object can be converted among different representations,
  to allow inter-operability among languages. In C, there are tools to convert
  to/from FORTRAN, whereas in C++ there are tools to convert to/from FORTRAN and C.

  These handlers only contain a pointer to the real structure, so that
  when a plumed object is brought from one language to another,
  it brings a reference to the same environment.

  Moreover, to simplify life in all cases where a single Plumed object is
  required for the entire simulation (which covers most of the practical
  applications with conventional MD codes) it is possible to take advantage
  of a global interface, which is implicitly referring to a unique global instance.
  The global object should still be initialized and finalized properly.

  The basic method to send a message to plumed is
\verbatim
  (C) plumed_cmd
  (C++) PLMD::Plumed::cmd
  (FORTRAN)  PLUMED_F_CMD
\endverbatim

  To initialize a plumed object, use:
\verbatim
  (C)        plumed_create
  (C++)      (constructor of PLMD::Plumed)
  (FORTRAN)  PLUMED_F_CREATE
\endverbatim

  To finalize it, use
\verbatim
  (C)        plumed_finalize
  (C++)      (destructor of PLMD::Plumed)
  (FORTRAN)  PLUMED_F_FINALIZE
\endverbatim

  To access to the global-object, use
\verbatim
  (C)        plumed_gcreate, plumed_gfinalize, plumed_gcmd
  (C++)      PLMD::Plumed::gcreate, PLMD::Plumed::gfinalize, PLMD::Plumed::gcmd
  (FORTRAN)  PLUMED_F_GCREATE, PLUMED_F_GFINALIZE, PLUMED_F_GCMD
\endverbatim

  To check if the global object has been initialized, use
\verbatim
  (C)        plumed_ginitialized
  (C++)      PLMD::Plumed::ginitialized
  (FORTRAN)  PLUMED_F_GINITIALIZED
\endverbatim

  To check if plumed library is available (this is useful for runtime linking), use
\verbatim
  (C)        plumed_installed
  (C++)      PLMD::Plumed::installed
  (FORTRAN)  PLUMED_F_INSTALLED
\endverbatim

  To convert handlers use
\verbatim
  (C)        plumed_c2f                 (C to FORTRAN)
  (C)        plumed_f2c                 (FORTRAN to C)
  (C++)      Plumed(plumed) constructor (C to C++)
  (C++)      operator plumed() cast     (C++ to C)
  (C++)      Plumed(char*)  constructor (FORTRAN to C++)
  (C++)      toFortran(char*)           (C++ to FORTRAN)
\endverbatim

\verbatim
  FORTRAN interface
    SUBROUTINE PLUMED_F_INSTALLED(i)
      INTEGER,           INTENT(OUT)   :: i
    SUBROUTINE PLUMED_F_GINITIALIZED(i)
      INTEGER,           INTENT(OUT)   :: i
    SUBROUTINE PLUMED_F_GCREATE()
    SUBROUTINE PLUMED_F_GCMD(key,val)
      CHARACTER(LEN=*), INTENT(IN)     :: key
      UNSPECIFIED_TYPE, INTENT(INOUT)  :: val(*)
    SUBROUTINE PLUMED_F_GFINALIZE()
    SUBROUTINE PLUMED_F_GLOBAL(p)
      CHARACTER(LEN=32), INTENT(OUT)   :: p
    SUBROUTINE PLUMED_F_CREATE(p)
      CHARACTER(LEN=32), INTENT(OUT)   :: p
    SUBROUTINE PLUMED_F_CMD(p,key,val)
      CHARACTER(LEN=32), INTENT(IN)    :: p
      CHARACTER(LEN=*),  INTENT(IN)    :: key
      UNSPECIFIED_TYPE,  INTENT(INOUT) :: val(*)
    SUBROUTINE PLUMED_F_FINALIZE(p)
      CHARACTER(LEN=32), INTENT(IN)    :: p
\endverbatim

  The main routine is "cmd", which accepts two arguments:
  key is a string containing the name of the command
  val is the argument. it is declared const so as to use allow passing const objects, but in practice plumed
      is going to modify val in several cases (using a const_cast).
  In some cases val can be omitted: just pass a NULL pointer (in C++, val is optional and can be omitted).
  The set of possible keys is the real API of the plumed library, and will be expanded with time.
  New commands will be added, but backward compatibility will be retained as long as possible.

  To pass plumed a callback function use the following syntax (not available in FORTRAN yet)
\verbatim
    plumed_function_holder ff;
    ff.p=your_function;
    plumed_cmd(plumed,"xxxx",&ff);
\endverbatim
  (this is passing the your_function() function to the "xxxx" command)
*/

#ifdef __cplusplus
extern "C" {
#endif

/* Generic function pointer */
typedef void (*plumed_function_pointer)(void);

/**
  \brief Holder for function pointer.

  To pass plumed a callback function use the following syntax:
\verbatim
    plumed_function_holder ff;
    ff.p=your_function;
    plumed_cmd(plumed,"xxxx",&ff);
\endverbatim
  (this is going to pass the your_function() function to the "xxxx" command)
*/

typedef struct {
  plumed_function_pointer p;
} plumed_function_holder;

/**
  \brief Main plumed object

  This is an object containing a Plumed instance, which should be used in
  the MD engine. It should first be initialized with plumed_create(),
  then it communicates with the MD engine using plumed_cmd(). Finally,
  before the termination, it should be deallocated with plumed_finalize().
  Its interface is very simple and general, and is expected
  not to change across plumed versions. See \ref ReferencePlumedH.
*/
typedef struct {
  /**
    \private
    \brief Void pointer holding the real PlumedMain structure
  */
  void*p;
} plumed;

/** \relates plumed
    \brief Constructor

    \return The constructed plumed object
*/
plumed plumed_create(void);

/** \relates plumed
    \brief Tells p to execute a command

    \param p The plumed object on which command is acting
    \param key The name of the command to be executed
    \param val The argument. It is declared as const to allow calls like plumed_cmd(p,"A","B"),
               but for some choice of key it can change the content
*/
void plumed_cmd(plumed p,const char*key,const void*val);

/** \relates plumed
    \brief Destructor

    \param p The plumed object to be deallocated
*/
void plumed_finalize(plumed p);

/** \relates plumed
    \brief Check if plumed is installed (for runtime binding)

    \return 1 if plumed is installed, 0 otherwise
*/
int plumed_installed(void);

/** \relates plumed
    \brief Retrieves an handler to the global structure.
*/
plumed plumed_global(void);

/** \relates plumed
    \brief Check if the global interface has been initialized

    \return 1 if plumed has been initialized, 0 otherwise
*/
int plumed_ginitialized(void);

/* global C interface, working on a global object */

/** \relates plumed
    \brief Constructor for the global interface.

    \note Equivalent to plumed_create(), but initialize the static global plumed object
*/
void plumed_gcreate(void);

/** \relates plumed
    \brief Tells to the global interface to execute a command.

    \param key The name of the command to be executed
    \param val The argument. It is declared as const to allow calls like plumed_gcmd("A","B"),
               but for some choice of key it can change the content

    \note Equivalent to plumed_cmd(), but acting on the global plumed object.
          It thus does not require the plumed object to be specified.
*/
void plumed_gcmd(const char* key,const void* val);

/** \relates plumed
    \brief Destructor for the global interface.

    \note Equivalent to plumed_finalize(), but acting on the global plumed object.
          It thus does not require the plumed object to be specified.
*/
void plumed_gfinalize(void);

/* routines to convert char handler from/to plumed objects */

/** \related plumed
    \brief Converts a C handler to a FORTRAN handler

    \param p The C handler
    \param c The FORTRAN handler (a char[32])

    This function can be used to convert a plumed object created in C to
    a plumed handler that can be used in FORTRAN.
\verbatim
#include <plumed/wrapper/Plumed.h>
int main(int argc,char*argv[]){
  plumed p;
  p=plumed_create();
  char fortran_handler[32];
  plumed_c2f(p,fortran_handler);
  printf("DEBUG: this is a string representation for the plumed handler: %s\n",fortran_handler);
  fortran_routine(fortran_handler);
  plumed_finalize(p);
  return 0;
}
\endverbatim
  Here `fortran_routine` is a routine implemented in FORTRAN that manipulates the
  fortran_handler.
*/
void   plumed_c2f(plumed p,char* c);

/** \related plumed
    \brief Converts a FORTRAN handler to a C handler
    \param c The FORTRAN handler (a char[32])
    \return The C handler

    This function can be used to convert a plumed object created in FORTRAN
    to a plumed handler that can be used in C.
\verbatim
void c_routine(char handler[32]){
  plumed p;
  p=plumed_f2c(handler);
  plumed_cmd(p,"init",NULL);
}
\endverbatim
  Here `c_routine` is a C function that can be called from FORTRAN
  and interact with the provided plumed handler.
*/
plumed plumed_f2c(const char* c);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

/* this is to include the NULL pointer */
#include <cstdlib>

/* C++ interface is hidden in PLMD namespace (same as plumed library) */
namespace PLMD {

/**
  C++ wrapper for \ref plumed.

  This class provides a C++ interface to PLUMED.
*/

class Plumed {
  /**
    C structure.
  */
  plumed main;
  /**
     keeps track if the object was created from scratch using
     the defaults destructor (reference=false) or if it was imported
     from C or FORTRAN (reference=true). In the latter case, the
     plumed_finalize() method is not called when destructing the object,
     since it is expected to be finalized in the C/FORTRAN code
  */
  bool reference;
public:
  /**
     Check if plumed is installed (for runtime binding)
     \return true if plumed is installed, false otherwise
     \note Equivalent to plumed_installed() but returns a bool
  */
  static bool installed();
  /**
     Check if global-plumed has been initialized
     \return true if global plumed object (see global()) is initialized (i.e. if gcreate() has been
             called), false otherwise.
     \note Equivalent to plumed_ginitialized() but returns a bool
  */
  static bool ginitialized();
  /**
     Initialize global-plumed.
     \note Equivalent to plumed_gcreate()
  */
  static void gcreate();
  /**
     Send a command to global-plumed
      \param key The name of the command to be executed
      \param val The argument. It is declared as const to allow calls like gcmd("A","B"),
                 but for some choice of key it can change the content
     \note Equivalent to plumed_gcmd()
  */
  static void gcmd(const char* key,const void* val);
  /**
     Finalize global-plumed
  */
  static void gfinalize();
  /**
     Returns the Plumed global object
     \return The Plumed global object
  */
  static Plumed global();
  /**
     Constructor.
    \note Performs the same task a plumed_create()
  */
  Plumed();
  /**
     Clone a Plumed object from a FORTRAN char* handler
     \param c The FORTRAN handler (a char[32]).

   \attention The Plumed object created in this manner
              will not finalize the corresponding plumed structure.
              It is expected that the FORTRAN code calls plumed_c_finalize for it
  */
// to have maximum portability of this file I do not use the explicit keyword here
// I thus add a suppress command for cppcheck
// cppcheck-suppress noExplicitConstructor
  Plumed(const char*c);
  /**
     Clone a Plumed object from a C plumed structure
     \param p The C plumed structure.

   \attention The Plumed object created in this manner
              will not finalize the corresponding plumed structure.
              It is expected that the C code calls plumed_finalize for it
  */
// to have maximum portability of this file I do not use the explicit keyword here
// I thus add a suppress command for cppcheck
// cppcheck-suppress noExplicitConstructor
  Plumed(plumed p);
private:
  /** Copy constructor is disabled (private and unimplemented)
    The problem here is that after copying it will not be clear who is
    going to finalize the corresponding plumed structure.
  */
  Plumed(const Plumed&);
  /** Assignment operator is disabled (private and unimplemented)
    The problem here is that after copying it will not be clear who is
    going to finalize the corresponding plumed structure.
  */
  Plumed&operator=(const Plumed&);
public:
  /*
    PLUMED 2.4 requires a C++11 compiler.
    Anyway, since Plumed.h file might be redistributed with other codes
    and it should be possible to combine it with earlier PLUMED versions,
    we here explicitly check if C+11 is available before enabling move semantics.
    This could still create problems if a compiler 'cheats', setting  __cplusplus > 199711L
    but not supporting move semantics. Hopefully will not happen!
  */
#if __cplusplus > 199711L
  /** Move constructor.
    Only if move semantics is enabled.
    It allows storing PLMD::Plumed objects in STL containers.
  */
  Plumed(Plumed&&);
  /** Move assignment.
    Only if move semantics is enabled.
  */
  Plumed& operator=(Plumed&&);
#endif
  /**
     Retrieve the C plumed structure for this object
  */
  operator plumed()const;
  /**
     Retrieve a FORTRAN handler for this object
      \param c The FORTRAN handler (a char[32]).
  */
  void toFortran(char*c)const;
  /**
     Send a command to this plumed object
      \param key The name of the command to be executed
      \param val The argument. It is declared as const to allow calls like p.cmd("A","B"),
                 but for some choice of key it can change the content
      \note Equivalent to plumed_cmd()
  */
  void cmd(const char*key,const void*val=NULL);
  /**
     Destructor

     Destructor is virtual so as to allow correct inheritance from Plumed object.
     To avoid linking problems with g++, I specify "inline" also here (in principle
     it should be enough to specify it down in the definition of the function, but
     for some reason that I do not understand g++ does not inline it properly in that
     case and complains when Plumed.h is included but Plumed.o is not linked. Anyway, the
     way it is done here seems to work properly).
  */
  inline virtual ~Plumed();
};

/* All methods are inlined so as to avoid the compilation of an extra c++ file */

inline
bool Plumed::installed() {
  return plumed_installed();
}

inline
Plumed::Plumed():
  main(plumed_create()),
  reference(false)
{}

inline
Plumed::Plumed(const char*c):
  main(plumed_f2c(c)),
  reference(true)
{}

inline
Plumed::Plumed(plumed p):
  main(p),
  reference(true)
{}

#if __cplusplus > 199711L
inline
Plumed::Plumed(Plumed&& p):
  main(p.main),
  reference(p.reference)
{}

inline
Plumed& Plumed::operator=(Plumed&& p) {
  main=p.main;
  reference=p.reference;
  return *this;
}
#endif

inline
Plumed::operator plumed()const {
  return main;
}

inline
void Plumed::toFortran(char*c)const {
  plumed_c2f(main,c);
}

inline
void Plumed::cmd(const char*key,const void*val) {
  plumed_cmd(main,key,val);
}

inline
Plumed::~Plumed() {
  if(!reference)plumed_finalize(main);
}

inline
bool Plumed::ginitialized() {
  return plumed_ginitialized();
}

inline
void Plumed::gcreate() {
  plumed_gcreate();
}

inline
void Plumed::gcmd(const char* key,const void* val) {
  plumed_gcmd(key,val);
}

inline
void Plumed::gfinalize() {
  plumed_gfinalize();
}

inline
Plumed Plumed::global() {
  return plumed_global();
}

}

#endif


#endif
