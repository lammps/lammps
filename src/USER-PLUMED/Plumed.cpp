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
#include "Plumed.h"

#ifdef __PLUMED_HAS_DLOPEN
#include <dlfcn.h>
#endif
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <limits.h>

/* DECLARATION USED ONLY IN THIS FILE */

#ifdef __cplusplus
extern "C" {
#endif

/**
  Function pointer to plumed_create
*/

typedef void*(*plumed_create_pointer)(void);
/**
  Function pointer to plumed_cmd
*/
typedef void(*plumed_cmd_pointer)(void*,const char*,const void*);

/**
  Function pointer to plumed_finalize
*/
typedef void(*plumed_finalize_pointer)(void*);

/**
   Holder for plumedmain function pointers.
*/
typedef struct {
  plumed_create_pointer create;
  plumed_cmd_pointer cmd;
  plumed_finalize_pointer finalize;
} plumed_plumedmain_function_holder;

/**
   Holder for plumed symbol table.
*/
typedef struct {
  int version;
  plumed_plumedmain_function_holder functions;
} plumed_symbol_table_type;

/**
  Register for plumedmain function pointers
*/
plumed_plumedmain_function_holder* plumed_kernel_register(const plumed_plumedmain_function_holder*);

#ifdef __PLUMED_STATIC_KERNEL
/* Real interface */
void*plumed_plumedmain_create(void);
void plumed_plumedmain_cmd(void*,const char*,const void*);
void plumed_plumedmain_finalize(void*);
#else
/* dummy interface */
void*plumed_dummy_create(void);
void plumed_dummy_cmd(void*,const char*,const void*);
void plumed_dummy_finalize(void*);
#endif

#ifdef __cplusplus
}
#endif

/* END OF DECLARATION USED ONLY IN THIS FILE */

/* These are the dummy routines which are used when plumed is not available */

#ifdef __PLUMED_STATIC_KERNEL

static int installed=1;

#else

static int installed=0;

static int dummy;

void*plumed_dummy_create(void) {
  return (void*)&dummy;
}

void plumed_dummy_cmd(void*p,const char*key,const void*val) {
  (void) p;   /* avoid warning on unused parameter */
  (void) key; /* avoid warning on unused parameter */
  (void) val; /* avoid warning on unused parameter */
  fprintf(stderr,"+++ ERROR: you are trying to use plumed, but it is not available +++\n");
  fprintf(stderr,"+++ Check your PLUMED_KERNEL environment variable +++\n");
  exit(1);
}

void plumed_dummy_finalize(void*p) {
  (void) p; /* avoid warning on unused parameter */
}

#endif

plumed_plumedmain_function_holder* plumed_kernel_register(const plumed_plumedmain_function_holder* f) {
  /*
    Argument f is present for historical reasons but ignored in PLUMED>=2.5.
  */
  if(f) {
    if(getenv("PLUMED_LOAD_DEBUG")) {
      fprintf(stderr,"+++ Ignoring registration at %p (%p,%p,%p) +++\n",(void*)f,(void*)f->create,(void*)f->cmd,(void*)f->finalize);
    }
  }
#ifdef __PLUMED_STATIC_KERNEL
  /*
    When __PLUMED_STATIC_KERNEL is defined, the function holder is initialized
    to statically bound plumed_plumedmain_create, plumed_plumedmain_cmd, plumed_plumedmain_finalize and
    cannot be changed. This saves from mis-set values for PLUMED_KERNEL.
  */
  static plumed_plumedmain_function_holder g= {plumed_plumedmain_create,plumed_plumedmain_cmd,plumed_plumedmain_finalize};
#else
  /*
    On the other hand, for runtime binding, we use dlsym to find the relevant functions.
    Notice that as of PLUMED 2.5 self registration of the kernel is ignored, so argument f
    is not used anymore.
    Also notice that we should put some guard here for safe multithread calculations.
  */
  static plumed_plumedmain_function_holder g= {plumed_dummy_create,plumed_dummy_cmd,plumed_dummy_finalize};
  static int first=1;
#ifdef __PLUMED_HAS_DLOPEN
  const char* path;
  char* pathcopy;
  void* p;
  char* pc;
  plumed_symbol_table_type* plumed_symbol_table_ptr;
  plumed_plumedmain_function_holder functions;
  char* debug;
  size_t strlenpath;
  int dlopenmode;
  /*
    f==NULL is required here otherwise we would enter this block a second time
    when plumed_kernel_register is called by the just loaded shared library.
  */
  if(first && f==NULL) {
    path=getenv("PLUMED_KERNEL");
    debug=getenv("PLUMED_LOAD_DEBUG");
#ifdef __PLUMED_DEFAULT_KERNEL
    /*
      This variable allows a default path for the kernel to be hardcoded.
      Can be useful for hardcoding the predefined plumed location
      still allowing the user to override this choice setting PLUMED_KERNEL.
      The path should be chosen at compile time adding e.g.
      -D__PLUMED_DEFAULT_KERNEL=/opt/local/lib/libplumed.dylib
    */
    /* This is required to add quotes */
#define PLUMED_QUOTE_DIRECT(name) #name
#define PLUMED_QUOTE(macro) PLUMED_QUOTE_DIRECT(macro)
    if(! (path && (*path) )) path=PLUMED_QUOTE(__PLUMED_DEFAULT_KERNEL);
#endif
    if(path && (*path)) {
      fprintf(stderr,"+++ Loading the PLUMED kernel runtime +++\n");
      fprintf(stderr,"+++ PLUMED_KERNEL=\"%s\" +++\n",path);
      if(getenv("PLUMED_LOAD_NAMESPACE") && !strcmp(getenv("PLUMED_LOAD_NAMESPACE"),"LOCAL")) {
        dlopenmode=RTLD_NOW|RTLD_LOCAL;
        if(debug) fprintf(stderr,"+++ Loading with mode RTLD_NOW|RTLD_LOCAL +++\n");
      } else {
        dlopenmode=RTLD_NOW|RTLD_GLOBAL;
        if(debug) fprintf(stderr,"+++ Loading with mode RTLD_NOW|RTLD_GLOBAL +++\n");
      }
      p=dlopen(path,dlopenmode);
      if(!p) {
        /*
          Something went wrong. We try to remove "Kernel" string from the PLUMED_KERNEL variable
          and load directly the shared library. Notice that this particular path is only expected
          to be necessary when using PLUMED<=2.4 and the symbols in the main executable are
          not visible. All the other cases (either PLUMED>=2.5 or symbols in the main executable visible)
          should work correctly without entering here.
        */
        fprintf(stderr,"+++ An error occurred. Message from dlopen(): %s +++\n",dlerror());
        strlenpath=strlen(path);
        pathcopy=(char*) malloc(strlenpath+1);
        strncpy(pathcopy,path,strlenpath+1);
        pc=pathcopy+strlenpath-6;
        while(pc>=pathcopy && memcmp(pc,"Kernel",6)) pc--;
        if(pc>=pathcopy) {
          memmove(pc, pc+6, strlen(pc)-5);
          fprintf(stderr,"+++ Trying %s +++\n",pathcopy);
          p=dlopen(pathcopy,dlopenmode);
          if(!p) fprintf(stderr,"+++ An error occurred. Message from dlopen(): %s +++\n",dlerror());
        }
        free(pathcopy);
      }
      if(p) {
        functions.create=NULL;
        functions.cmd=NULL;
        functions.finalize=NULL;
      /*
        If the library was loaded, use dlsym to initialize pointers.
        Notice that as of PLUMED 2.5 we ignore self registrations.
        Pointers are searched in the form of a single pointer to a structure, which
        is the standard way in PLUMED 2.5, as well as using alternative names used in
        PLUMED 2.0 to 2.4 (e.g. plumedmain_create) and in some intermediate versions between
        PLUMED 2.4 and 2.5 (e.g. plumed_plumedmain_create). The last chance is probably
        unnecessary and might be removed at some point.
      */
        plumed_symbol_table_ptr=(plumed_symbol_table_type*) dlsym(p,"plumed_symbol_table");
        if(plumed_symbol_table_ptr) functions=plumed_symbol_table_ptr->functions;
        if(debug && plumed_symbol_table_ptr) {
          fprintf(stderr,"+++ plumed_symbol_table version %i found at %p +++\n",plumed_symbol_table_ptr->version,(void*)plumed_symbol_table_ptr);
          fprintf(stderr,"+++ plumed_function_pointers found at %p (%p,%p,%p) +++\n",(void*)&plumed_symbol_table_ptr->functions,(void*)functions.create,(void*)functions.cmd,(void*)functions.finalize);
        }

        if(!functions.create) {
          functions.create=(plumed_create_pointer) dlsym(p,"plumedmain_create");
          if(debug && functions.create) fprintf(stderr,"+++ %s found at %p +++\n","plumedmain_create",(void*)functions.create);
        }
        if(!functions.create) {
          functions.create=(plumed_create_pointer) dlsym(p,"plumed_plumedmain_create");
          if(debug && functions.create) fprintf(stderr,"+++ %s found at %p +++\n","plumed_plumedmain_create",(void*)functions.create);
        }

        if(!functions.cmd) {
          functions.cmd=(plumed_cmd_pointer) dlsym(p,"plumedmain_cmd");
          if(debug && functions.cmd) fprintf(stderr,"+++ %s found at %p +++\n","plumedmain_cmd",(void*)functions.cmd);
        }
        if(!functions.cmd) {
          functions.cmd=(plumed_cmd_pointer) dlsym(p,"plumed_plumedmain_cmd");
          if(debug && functions.cmd) fprintf(stderr,"+++ %s found at %p +++\n","plumed_plumedmain_cmd",(void*)functions.cmd);
        }

        if(!functions.finalize) {
          functions.finalize=(plumed_finalize_pointer) dlsym(p,"plumedmain_finalize");
          if(debug && functions.finalize) fprintf(stderr,"+++ %s found at %p +++\n","plumedmain_finalize",(void*)functions.finalize);
        }
        if(!functions.finalize) {
          functions.finalize=(plumed_finalize_pointer) dlsym(p,"plumed_plumedmain_finalize");
          if(debug && functions.finalize) fprintf(stderr,"+++ %s found at %p +++\n","plumed_plumedmain_finalize",(void*)functions.finalize);
        }

        if(functions.create && functions.cmd && functions.finalize) {
          g=functions;
          installed=1;
        } else {
          if(!functions.create) fprintf(stderr,"+++ pointer to (plumed_)plumedmain_create not found +++\n");
          if(!functions.cmd) fprintf(stderr,"+++ pointer to (plumed_)plumedmain_cmd not found +++\n");
          if(!functions.finalize) fprintf(stderr,"+++ pointer to (plumed_)plumedmain_finalize not found +++\n");
        }
      }
    }
    first=0;
  }
#endif
#endif
  return &g;
}

/* C wrappers: */

plumed plumed_create(void) {
  plumed p;
  plumed_plumedmain_function_holder*h=plumed_kernel_register(NULL);
  assert(h);
  assert(h->create);
  p.p=(*(h->create))();
  assert(p.p);
  return p;
}

void plumed_cmd(plumed p,const char*key,const void*val) {
  plumed_plumedmain_function_holder*h=plumed_kernel_register(NULL);
  assert(p.p);
  assert(h);
  assert(h->cmd);
  (*(h->cmd))(p.p,key,val);
}

void plumed_finalize(plumed p) {
  plumed_plumedmain_function_holder*h=plumed_kernel_register(NULL);
  assert(p.p);
  assert(h);
  assert(h->finalize);
  (*(h->finalize))(p.p);
  p.p=NULL;
}

int plumed_installed(void) {
  plumed_kernel_register(NULL);
  return installed;
}

/* we declare a Plumed_g_main object here, in such a way that it is always available */

static plumed gmain= {NULL};

plumed plumed_global(void) {
  return gmain;
}

void plumed_gcreate(void) {
  assert(gmain.p==NULL);
  gmain=plumed_create();
}

void plumed_gcmd(const char*key,const void*val) {
  assert(gmain.p);
  plumed_cmd(gmain,key,val);
}

void plumed_gfinalize(void) {
  assert(gmain.p);
  plumed_finalize(gmain);
  gmain.p=NULL;
}

int plumed_ginitialized(void) {
  if(gmain.p) return 1;
  else                return 0;
}

void plumed_c2f(plumed p,char*c) {
  unsigned i;
  unsigned char* cc;
  /*
    Convert the address stored in p.p into a proper FORTRAN string
    made of only ASCII characters. For this to work, the two following
    assertions should be satisfied:
  */
  assert(CHAR_BIT<=12);
  assert(sizeof(p.p)<=16);

  assert(c);
  cc=(unsigned char*)&p.p;
  for(i=0; i<sizeof(p.p); i++) {
    /*
      characters will range between '0' (ASCII 48) and 'o' (ASCII 111=48+63)
    */
    c[2*i]=cc[i]/64+48;
    c[2*i+1]=cc[i]%64+48;
  }
}

plumed plumed_f2c(const char*c) {
  plumed p;
  unsigned i;
  unsigned char* cc;

  assert(CHAR_BIT<=12);
  assert(sizeof(p.p)<=16);

  assert(c);
  cc=(unsigned char*)&p.p;
  for(i=0; i<sizeof(p.p); i++) {
    /*
      perform the reversed transform
    */
    cc[i]=(c[2*i]-48)*64 + (c[2*i+1]-48);
  }
  return p;
}


#ifdef __cplusplus
extern "C" {
#endif

/*
  Fortran wrappers
  These are just like the global C wrappers. They are
  just defined here and not declared in the .h file since they
  should not be used from c/c++ anyway.
*/

/*
  First we assume no name mangling
*/

void plumed_f_installed(int*i) {
  assert(i);
  *i=plumed_installed();
}

void plumed_f_ginitialized(int*i) {
  assert(i);
  *i=plumed_ginitialized();
}

void plumed_f_gcreate(void) {
  plumed_gcreate();
}

void plumed_f_gcmd(char*key,void*val) {
  plumed_gcmd(key,val);
}

void plumed_f_gfinalize(void) {
  plumed_gfinalize();
}

void plumed_f_create(char*c) {
  plumed p;
  p=plumed_create();
  plumed_c2f(p,c);
}

void plumed_f_cmd(char*c,char*key,void*val) {
  plumed p;
  p=plumed_f2c(c);
  plumed_cmd(p,key,val);
}

void plumed_f_finalize(char*c) {
  plumed p;
  p=plumed_f2c(c);
  plumed_finalize(p);
}

void plumed_f_global(char*c) {
  plumed_c2f(gmain,c);
}

/*
  Then we add wrappers for there functions to cover all
  the possible fortran mangling schemes, which should be:
  without underscore, with one underscore and with two underscores
  lower or upper case
*/

#define IMPLEMENT(lower,upper,implem) \
  void lower ##_  implem \
  void lower ##__ implem \
  void upper      implem \
  void upper ##_  implem \
  void upper ##__ implem

IMPLEMENT(plumed_f_gcreate,     PLUMED_F_GCREATE,     (void) {plumed_f_gcreate();})
IMPLEMENT(plumed_f_gcmd,        PLUMED_F_GCMD,        (char* key,void* val) {plumed_f_gcmd(key,val);})
IMPLEMENT(plumed_f_gfinalize,   PLUMED_F_GFINALIZE,   (void) {plumed_f_gfinalize();})
IMPLEMENT(plumed_f_ginitialized,PLUMED_F_GINITIALIZED,(int*i) {plumed_f_ginitialized(i);})
IMPLEMENT(plumed_f_create,      PLUMED_F_CREATE,      (char*c) {plumed_f_create(c);})
IMPLEMENT(plumed_f_cmd,         PLUMED_F_CMD,         (char*c,char* key,void* val) {plumed_f_cmd(c,key,val);})
IMPLEMENT(plumed_f_finalize,    PLUMED_F_FINALIZE,    (char*c) {plumed_f_finalize(c);})
IMPLEMENT(plumed_f_installed,   PLUMED_F_INSTALLED,   (int*i) {plumed_f_installed(i);})
IMPLEMENT(plumed_f_global,      PLUMED_F_GLOBAL,      (char*c) {plumed_f_global(c);})

#ifdef __cplusplus
}
#endif




