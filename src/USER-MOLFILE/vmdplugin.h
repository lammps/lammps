/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2006 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: vmdplugin.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.32 $       $Date: 2009/02/24 05:12:35 $
 *
 ***************************************************************************/

/** @file
 * This header must be included by every VMD plugin library.  It defines the
 * API for every plugin so that VMD can organize the plugins it finds.
 */

#ifndef VMD_PLUGIN_H
#define VMD_PLUGIN_H


/*
 * Preprocessor tricks to make it easier for us to redefine the names of
 * functions when building static plugins.
 */
#if !defined(VMDPLUGIN)
/**
  * macro defining VMDPLUGIN if it hasn't already been set to the name of
  * a static plugin that is being compiled.  This is the catch-all case.
  */
#define VMDPLUGIN vmdplugin
#endif
/** concatenation macro, joins args x and y together as a single string */
#define xcat(x, y) cat(x, y)
/** concatenation macro, joins args x and y together as a single string */
#define cat(x, y) x ## y

/*
 *  macros to correctly define plugin function names depending on whether
 *  the plugin is being compiled for static linkage or dynamic loading.
 *  When compiled for static linkage, each plugin needs to have unique
 *  function names for all of its entry points.  When compiled for dynamic
 *  loading, the plugins must name their entry points consistently so that
 *  the plugin loading mechanism can find the register, register_tcl, init,
 *  and fini routines via dlopen() or similar operating system interfaces.
 */
/*@{*/
/** Macro names entry points correctly for static linkage or dynamic loading */
#define VMDPLUGIN_register     xcat(VMDPLUGIN, _register)
#define VMDPLUGIN_register_tcl xcat(VMDPLUGIN, _register_tcl)
#define VMDPLUGIN_init         xcat(VMDPLUGIN, _init)
#define VMDPLUGIN_fini         xcat(VMDPLUGIN, _fini)
/*@}*/


/** "WIN32" is defined on both WIN32 and WIN64 platforms... */
#if (defined(WIN32))
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#if !defined(STATIC_PLUGIN)
#if defined(VMDPLUGIN_EXPORTS)
/**
 *  Only define DllMain for plugins, not in VMD or in statically linked plugins
 *  VMDPLUGIN_EXPORTS is only defined when compiling dynamically loaded plugins
 */
BOOL APIENTRY DllMain( HANDLE hModule,
                       DWORD ul_reason_for_call,
                       LPVOID lpReserved
                     )
{
  return TRUE;
}

#define VMDPLUGIN_API __declspec(dllexport)
#else
#define VMDPLUGIN_API __declspec(dllimport)
#endif /* VMDPLUGIN_EXPORTS */
#else  /* ! STATIC_PLUGIN */
#define VMDPLUGIN_API
#endif /* ! STATIC_PLUGIN */
#else
/** If we're not compiling on Windows, then this macro is defined empty */
#define VMDPLUGIN_API
#endif

/** define plugin linkage correctly for both C and C++ based plugins */
#ifdef __cplusplus
#define VMDPLUGIN_EXTERN extern "C" VMDPLUGIN_API
#else
#define VMDPLUGIN_EXTERN extern VMDPLUGIN_API
#endif  /* __cplusplus */

/*
 * Plugin API functions start here
 */


/**
 * Init routine: called the first time the library is loaded by the
 * application and before any other API functions are referenced.
 * Return 0 on success.
 */
VMDPLUGIN_EXTERN int VMDPLUGIN_init(void);

/**
 * Macro for creating a struct header used in all plugin structures.
 *
 * This header should be placed at the top of every plugin API definition
 * so that it can be treated as a subtype of the base plugin type.
 *
 * abiversion: Defines the ABI for the base plugin type (not for other plugins)
 * type: A string descriptor of the plugin type.
 * name: A name for the plugin.
 * author: A string identifier, possibly including newlines.
 * Major and minor version.
 * is_reentrant: Whether this library can be run concurrently with itself.
 */
#define vmdplugin_HEAD \
  int abiversion; \
  const char *type; \
  const char *name; \
  const char *prettyname; \
  const char *author; \
  int majorv; \
  int minorv; \
  int is_reentrant;

/**
  * Typedef for generic plugin header, individual plugins can
  * make their own structures as long as the header info remains
  * the same as the generic plugin header, most easily done by
  * using the vmdplugin_HEAD macro.
  */
typedef struct {
  vmdplugin_HEAD
} vmdplugin_t;

/**
 * Use this macro to initialize the abiversion member of each plugin
 */
#define vmdplugin_ABIVERSION  16

/*@{*/
/** Use this macro to indicate a plugin's thread-safety at registration time */
#define VMDPLUGIN_THREADUNSAFE 0
#define VMDPLUGIN_THREADSAFE   1
/*@}*/

/*@{*/
/** Error return code for use in the plugin registration and init functions */
#define VMDPLUGIN_SUCCESS      0
#define VMDPLUGIN_ERROR       -1
/*@}*/

/**
 * Function pointer typedef for register callback functions
 */
typedef int (*vmdplugin_register_cb)(void *, vmdplugin_t *);

/**
 * Allow the library to register plugins with the application.
 * The callback should be called using the passed-in void pointer, which
 * should not be interpreted in any way by the library.  Each vmdplugin_t
 * pointer passed to the application should point to statically-allocated
 * or heap-allocated memory and should never be later modified by the plugin.
 * Applications must be permitted to retain only a copy of the the plugin
 * pointer, without making any deep copy of the items in the struct.
 */
VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *, vmdplugin_register_cb);

/**
 * Allow the library to register Tcl extensions.
 * This API is optional; if found by dlopen, it will be called after first
 * calling init and register.
 */
VMDPLUGIN_EXTERN int VMDPLUGIN_register_tcl(void *, void *tcl_interp,
    vmdplugin_register_cb);

/**
 * The Fini method is called when the application will no longer use
 * any plugins in the library.
 */
VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void);

#endif   /* VMD_PLUGIN_H */
