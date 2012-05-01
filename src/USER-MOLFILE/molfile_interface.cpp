/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple)
------------------------------------------------------------------------- */

#include "molfile_interface.h"

#include <sys/types.h>
#include <stdio.h>
#include <dirent.h>
#include <string.h>

#if defined(_WIN32)
#include <windows.h>
#else
#include <dirent.h>
#include <dlfcn.h>
#endif

#include "molfile_plugin.h"

#define DEBUG 0

extern "C" {
  typedef int (*initfunc)(void);
  typedef int (*regfunc)(void *, vmdplugin_register_cb);
  typedef int (*finifunc)(void);

  // callback function for plugin registration.
  static int plugin_register_cb(void *v, vmdplugin_t *p) 
  {
    molfile_plugin_t *m = static_cast<molfile_plugin_t *>(v);
    if (strcmp(m->name, p->name) == 0)
      memcpy(v,p,sizeof(molfile_plugin_t));
    return 0;
  }

  // directory traversal helper functions

#if defined(_WIN32)

  // Win32 directory traversal handle
  typedef struct {
    HANDLE h;
    WIN32_FIND_DATA fd;
    char *name;
    char *searchname;
    int dlen;
  } dirhandle_t;

  // open a directory handle
  static dirhandle_t *my_opendir(const char *dirname)
  {
    dirhandle_t *d;
    int len;

    if (dirname == NULL) 
      return NULL;
    d = new dirhandle_t;

    len = 2 + strlen(dirname);
    d->name = new char[len];
    strcpy(d->name, dirname);
    strcat(d->name, "\\");
    d->dlen = len;

    len += 1;
    d->searchname = new char[len];
    strcpy(d->searchname, dirname);
    strcat(d->searchname, "\\*");

    d->h = FindFirstFile(d->searchname, &(d->fd));
    if (d->h == ((HANDLE)(-1))) {
      delete[] d->searchname;
      delete[] d->name;
      delete d;
      return NULL;
    }
    return d;
  }

  // get next file name from directory handle
  static char *my_readdir(dirhandle_t *d)
  {
    if (FindNextFile(d->h, &(d->fd))) {
      return d->fd.cFileName; 
    }
    return NULL;     
  }

  // close directory handle
  static void my_closedir(dirhandle_t *d)
  {
    if (d->h != NULL) {
      FindClose(d->h);
    }
    delete[] d->searchname;
    delete[] d->name;
    delete d;
  }

  // open a shared object file
  static void *my_dlopen(const char *fname) {
    return (void *)LoadLibrary(fname);
  }

  // report error message from dlopen
  static const char *my_dlerror(void) {
    static CHAR szBuf[80]; 
    DWORD dw = GetLastError();

    sprintf(szBuf, "my_dlopen failed: GetLastError returned %u\n", dw); 
    return szBuf;
  }

  // resolve a symbol in shared object
  static void *my_dlsym(void *h, const char *sym) {
    return (void *)GetProcAddress((HINSTANCE)h, sym);
  }

  // close a shared object
  static int my_dlclose(void *h) {
    /* FreeLibrary returns nonzero on success */
    return !FreeLibrary((HINSTANCE)h);
  }

#else

  // Unix directory traversal handle
  typedef struct {
    DIR *d;
    char *name;
    int dlen;
  } dirhandle_t;

  // open a directory handle
  static dirhandle_t *my_opendir(const char *dirname)
  {
    dirhandle_t *d;
    int len;

    if (dirname == NULL) return NULL;

    d = new dirhandle_t;
    len = 2 + strlen(dirname);
    d->name = new char[len];
    strcpy(d->name,dirname);
    strcat(d->name,"/");
    d->dlen = len;

    d->d = opendir(d->name);
    if (d->d == NULL) {
      delete[] d->name;
      delete d;
      return NULL;
    }
    return d;
  }

  // get next file name from directory handle
  static char *my_readdir(dirhandle_t *d)
  {
    struct dirent *p;

    if ((p = readdir(d->d)) != NULL) {
      return p->d_name;
    }

    return NULL;
  }

  // close directory handle
  static void my_closedir(dirhandle_t *d)
  {
    if (d->d != NULL) {
      closedir(d->d);
    }
    delete[] d->name;
    delete d;
    return;
  }

  // open a shared object file
  static void *my_dlopen(const char *fname) {
    return dlopen(fname, RTLD_NOW);
  }

  // report error message from dlopen
  static const char *my_dlerror(void) {
    return dlerror();
  }

  // resolve a symbol in shared object
  static void *my_dlsym(void *h, const char *sym) {
    return dlsym(h, sym);
  }

  // close a shared object
  static int my_dlclose(void *h) {
    return dlclose(h);
  }

#endif

} // end of extern "C" region

using namespace LAMMPS_NS;

// constructor.
MolfileInterface::MolfileInterface()
  : _plugin(0), _dso(0), _rptr(0), _wptr(0), _natoms(0), _mode(M_NONE) 
{
  _name = new char[5];
  strcpy(_name,"none");
}

// destructor.
MolfileInterface::~MolfileInterface()
{
  forget_plugin();
  _mode = M_NONE;
  delete[] _name;
}

// register the best matching plugin in a given directory
int MolfileInterface::find_plugin(const char *plugindir,
				  const char *filetype,
				  const int mode)
{
  dirhandle_t *dir;
  char *filename, *ext;
  molfile_plugin_t *plugin;
  int retval = E_NONE;

  // evict current plugin, if the mode changes
  if (mode != _mode)
    forget_plugin();

  dir = my_opendir(plugindir);
  if (dir == NULL) return E_NODIR;
#if DEBUG
    printf("searching for plugins in directory: %s\n",plugindir);
#endif

  // search for suitable file names.
  while(1) {
    void *dso;
    char *fullname;
    int len;

    filename = my_readdir(dir);
    if (filename == NULL) break;

    // only look at .so files
    ext = strrchr(filename, '.');
    if (ext == NULL) continue;
    if (strcasecmp(ext,".so") != 0) continue;
#if DEBUG
    printf("next DSO file: %s\n",filename);
#endif

    // construct full pathname of potential DSO
    len = dir->dlen;
    len += strlen(filename);
    fullname = new char[len];
    strcpy(fullname,dir->name);
    strcat(fullname,filename);

#if DEBUG
    printf("trying to dlopen file %s\n",fullname);
#endif
    dso = my_dlopen(fullname);
    if (dso == NULL) {
      delete[] fullname;
      continue;
    }
#if DEBUG
    printf("dlopen handle: 0x%0x\n",dso);
#endif

    // check for required plugin symbols
    void *ifunc = my_dlsym(dso,"vmdplugin_init");
    void *rfunc = my_dlsym(dso,"vmdplugin_register");
    void *ffunc = my_dlsym(dso,"vmdplugin_fini");
    if (ifunc == NULL || rfunc == NULL || ffunc == NULL) {
      delete[] fullname;
      my_dlclose(dso);
      continue;
    }

    // intialize plugin. skip plugin if it fails.
    if (((initfunc)(ifunc))()) {
      delete[] fullname;
      my_dlclose(dso);
      continue;
    }

    // temporarily register plugin, i.e. make copy of plugin struct.
    plugin = new molfile_plugin_t;
    plugin->name=filetype;
    ((regfunc)rfunc)(plugin, plugin_register_cb);

#if DEBUG
    printf("plugin name: %s v%d.%d\n", plugin->name,
	   plugin->majorv, plugin->minorv);
    printf("plugin abiversion: %d\n", plugin->abiversion);
    printf("%s molfile plugin by %s\n", plugin->prettyname, plugin->author);
#endif

    // make some checks, if the plugin is suitable
    int use_this_plugin = 1;

    // check if the callback found a matching plugin
    // and thus overwrote the struct.
    if (plugin->name == filetype)
      use_this_plugin = 0;

    // check if the ABI matches the one used to compile this code
    if (plugin->abiversion != vmdplugin_ABIVERSION) {
      use_this_plugin = 0;
      if (retval < E_ABI)
	retval = E_ABI;
    }

    // check if (basic) reading is supported
    if ((mode & M_READ) &&
	( (plugin->open_file_read == NULL) ||
	  (plugin->read_next_timestep  == NULL) ||
	  (plugin->close_file_read == NULL) )) {
      use_this_plugin = 0;
      if (retval < E_MODE)
	retval = E_MODE;
    }

    // check if (basic) writing is supported
    if ( (mode & M_WRITE) &&
	 ( (plugin->open_file_write == NULL) ||
	   (plugin->write_timestep  == NULL) ||
	   (plugin->close_file_write == NULL) )) {
      use_this_plugin = 0;
      if (retval < E_MODE)
	retval = E_MODE;
    }

    // check if we have an updated plugin.
    if (_dso && _plugin) {
      molfile_plugin_t *p;
      p = static_cast<molfile_plugin_t *>(_plugin);

      if (p->majorv > plugin->majorv) {
	use_this_plugin = 0;
	if (retval < E_VERSION)
	  retval = E_VERSION;
      }

      if ( (p->majorv == plugin->majorv) &&
	   (p->minorv >= plugin->minorv) ) {
	use_this_plugin = 0;
	if (retval < E_VERSION)
	  retval = E_VERSION;
      }
    }

    // store plugin info in class if it qualifies
    if (use_this_plugin) {
      forget_plugin();

      len = 16;
      len += strlen(plugin->prettyname);
      len += strlen(plugin->author);
      delete[] _name;
      _name = new char[len];
      sprintf(_name,"%s v%d.%d by %s",plugin->prettyname,
	      plugin->majorv, plugin->minorv, plugin->author);

      _plugin = plugin;
      _dso = dso;
      _mode = mode;
      retval = E_MATCH;
    } else {

      delete plugin;
      my_dlclose(dso);
    }

    delete[] fullname;
  }
  my_closedir(dir);
  return retval;
}

// deregister a plugin
void MolfileInterface::forget_plugin()
{
  if (_plugin) {
    molfile_plugin_t *p;
    p = static_cast<molfile_plugin_t *>(_plugin);
    delete p;
  }

  if (_dso) {
    void *ffunc = my_dlsym(_dso,"vmdplugin_fini");
    if (ffunc)
      ((finifunc)ffunc)();
    my_dlclose(_dso);
  }

  strcpy(_name,"none");
  _plugin = NULL;
  _dso = NULL;
  _mode = M_NONE;
}

// open file for writing
int MolfileInterface::open_write(const char *name, const int natoms)
{
  if (!_plugin || !_dso || !(_mode & M_WRITE))
    return 1;
  molfile_plugin_t *p = static_cast<molfile_plugin_t *>(_plugin);
  
  _wptr = p->open_file_write(name,_type,natoms);
  if (_wptr == NULL)
    return 1;
  
  return 0;
}

