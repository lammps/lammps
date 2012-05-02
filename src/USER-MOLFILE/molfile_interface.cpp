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
#include <stdlib.h>

#if defined(_WIN32)
#include <windows.h>
#else
#include <dirent.h>
#include <dlfcn.h>
#endif

#include "molfile_plugin.h"

#if vmdplugin_ABIVERSION < 16
#error "unsupported VMD molfile plugin ABI version"
#endif

#define DEBUG 0

extern "C" {
  typedef int (*initfunc)(void);
  typedef int (*regfunc)(void *, vmdplugin_register_cb);
  typedef int (*finifunc)(void);

  typedef struct {
    void *p;
    const char *name;
  } plugin_reginfo_t;

  // callback function for plugin registration.
  static int plugin_register_cb(void *v, vmdplugin_t *p) 
  {
    plugin_reginfo_t *r = static_cast<plugin_reginfo_t *>(v);
    // make sure we have the proper plugin type (native reader)
    // for the desired file type (called "name" at this level)
    if ((strcmp(MOLFILE_PLUGIN_TYPE,p->type) == 0)
	&& (strcmp(r->name, p->name) == 0) ) {
      r->p = static_cast<void *>(p);
    }
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
MolfileInterface::MolfileInterface(const char *type, const int mode)
  : _plugin(0), _dso(0), _ptr(0), _info(0), _natoms(0),
    _mode(mode), _caps(M_NONE)
{
  _name = new char[5];
  strcpy(_name,"none");
  _type = new char[1+strlen(type)];
  strcpy(_type,type);
}

// destructor.
MolfileInterface::~MolfileInterface()
{
  forget_plugin();

  if (_info) {
    molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);
    delete[] a;
    _info = NULL;
  }
  delete[] _name;
  delete[] _type;
}

// register the best matching plugin in a given directory
int MolfileInterface::find_plugin(const char *pluginpath)
{
  dirhandle_t *dir;
  char *filename, *ext, *next, *path, *plugindir;
  int retval = E_NONE;

#if defined(_WIN32)
#define MY_PATHSEP ';'
#else
#define MY_PATHSEP ':'
#endif
  if (pluginpath == NULL) return E_DIR;
  plugindir = path = strdup(pluginpath);

  while (plugindir) {
    // check if this a single directory or path.
    next = strchr(plugindir,MY_PATHSEP);
    if (next) {
      *next = NULL;
      ++next;
    }

    dir = my_opendir(plugindir);
    if (!dir)
      retval = (retval > E_DIR) ? retval : E_DIR;

    // search for suitable file names and try to inspect them
    while(dir) {
      char *fullname;
      int len;

      filename = my_readdir(dir);
      if (filename == NULL) break;

      // only look at .so files
      ext = strrchr(filename, '.');
      if (ext == NULL) continue;
      if (strcasecmp(ext,".so") != 0) continue;

      // construct full pathname of potential DSO
      len = dir->dlen;
      len += strlen(filename);
      fullname = new char[len];
      strcpy(fullname,dir->name);
      strcat(fullname,filename);

      // try to register plugin at file name.
      int rv = load_plugin(fullname);
      if (rv > retval) retval = rv;

      delete[] fullname;
    }
    if (dir)
      my_closedir(dir);

    plugindir = next;
  }
  free(path);
  return retval;
}

// register the best matching plugin in a given directory
int MolfileInterface::load_plugin(const char *filename)
{
  void *dso;
  int len, retval = E_NONE;

  // access shared object
  dso = my_dlopen(filename);
  if (dso == NULL)
    return E_FILE;

  // check for required plugin symbols
  void *ifunc = my_dlsym(dso,"vmdplugin_init");
  void *rfunc = my_dlsym(dso,"vmdplugin_register");
  void *ffunc = my_dlsym(dso,"vmdplugin_fini");
  if (ifunc == NULL || rfunc == NULL || ffunc == NULL) {
    my_dlclose(dso);
    return E_SYMBOL;
  }

  // intialize plugin. skip plugin if it fails.
  if (((initfunc)(ifunc))()) {
    my_dlclose(dso);
    return E_SYMBOL;
  }

  // pre-register plugin.
  // the callback will be called for each plugin in the DSO and
  // check the file type. plugin->name will change if successful.
  plugin_reginfo_t reginfo;
  reginfo.p = NULL;
  reginfo.name=_type;
  ((regfunc)rfunc)(&reginfo, plugin_register_cb);

  // make some checks to see if the plugin is suitable or not.
  molfile_plugin_t *plugin = static_cast<molfile_plugin_t *>(reginfo.p);

  // if the callback found a matching plugin and copied the struct,
  // its name element will point to a different location now.
  if (plugin == NULL) {
    retval = E_TYPE;

    // check if the ABI matches the one used to compile this code
  } else if (plugin->abiversion != vmdplugin_ABIVERSION) {
    retval = E_ABI;

    // check if (basic) reading is supported
  } else if ((_mode & M_READ) &&
	     ( (plugin->open_file_read == NULL) ||
	       (plugin->read_next_timestep  == NULL) ||
	       (plugin->close_file_read == NULL) )) {
    retval = E_MODE;

    // check if (basic) writing is supported
  } else if ( (_mode & M_WRITE) &&
	      ( (plugin->open_file_write == NULL) ||
		(plugin->write_timestep  == NULL) ||
		(plugin->close_file_write == NULL) )) {
    retval = E_MODE;

    // make some additional check, if we
    // already have a plugin registered.
    // NOTE: this has to come last.
  } else if (_dso && _plugin) {
    molfile_plugin_t *p;
    p = static_cast<molfile_plugin_t *>(_plugin);

    // check if the new plugin is of a newer major version
    if (p->majorv > plugin->majorv) {
      retval = E_VERSION;

    // check if the new plugin is of a newer minor version
    } else if ( (p->majorv == plugin->majorv) &&
		(p->minorv >= plugin->minorv) ) {
      retval = E_VERSION;
    }
  }

  // bingo! this one is a keeper.
  if (retval == E_NONE) {

    // make sure any existing plugin is wiped out
    forget_plugin();

    delete[] _name;
    len = 16;
    len += strlen(plugin->prettyname);
    len += strlen(plugin->author);
    _name = new char[len];
    sprintf(_name,"%s v%d.%d by %s",plugin->prettyname,
	    plugin->majorv, plugin->minorv, plugin->author);

    // determine plugin capabilities
    _caps = M_NONE;
    if (plugin->read_next_timestep)    _caps |= M_READ;
    if (plugin->write_timestep)        _caps |= M_WRITE;
    if (plugin->read_structure)        _caps |= M_RSTRUCT;
    if (plugin->write_structure)       _caps |= M_WSTRUCT;
    if (plugin->read_bonds)            _caps |= M_RBONDS;
    if (plugin->write_bonds)           _caps |= M_WBONDS;
    if (plugin->read_angles)           _caps |= M_RANGLES;
    if (plugin->write_angles)          _caps |= M_WANGLES;
    if (plugin->read_volumetric_data)  _caps |= M_RVOL;
    if (plugin->write_volumetric_data) _caps |= M_WVOL;

    if (_mode & M_WRITE)
      _mode |= (_caps & M_WSTRUCT);
    else if (_mode & M_READ)
      _mode |= (_caps & M_RSTRUCT);

    _plugin = plugin;
    _dso = dso;
    return E_MATCH;
  }

  // better luck next time. clean up and return.
  my_dlclose(dso);
  return retval;
}

// deregister a plugin and close or reset all associated objects.
void MolfileInterface::forget_plugin()
{
  if (_ptr)
    close();

  if (_plugin)
    _plugin = NULL;

  if (_dso) {
    void *ffunc = my_dlsym(_dso,"vmdplugin_fini");
    if (ffunc)
      ((finifunc)ffunc)();
    my_dlclose(_dso);
  }
  _dso = NULL;

  delete[] _name;
    _name = new char[5];
  strcpy(_name,"none");

  _caps = M_NONE;
}

// open file for writing
int MolfileInterface::open(const char *name, int *natoms)
{
  if (!_plugin || !_dso || !natoms)
    return E_FILE;
  molfile_plugin_t *p = static_cast<molfile_plugin_t *>(_plugin);
  
  if (_mode & M_WRITE)
    _ptr = p->open_file_write(name,_type,*natoms);
  else if (_mode & M_READ)
    _ptr = p->open_file_read(name,_type,natoms);

  if (_ptr == NULL)
    return E_FILE;

  _natoms = *natoms;
  // we need to deal with structure information
  if (_mode & (M_RSTRUCT|M_WSTRUCT)) {
    molfile_atom_t *a = new molfile_atom_t[_natoms];
    _info = a;
    memset(_info,0,_natoms*sizeof(molfile_atom_t));
  }

  return E_NONE;
}

// safely close file
int MolfileInterface::close()
{
  if (!_plugin || !_dso || !_ptr)
    return E_FILE;

  molfile_plugin_t *p = static_cast<molfile_plugin_t *>(_plugin);

  if (_mode & M_WRITE) {
    p->close_file_write(_ptr);
  } else if (_mode & M_READ) {
    p->close_file_read(_ptr);
  }

  if (_info) {
    molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);
    delete[] a;
    _info = NULL;
  }
  _ptr = NULL;
  _natoms = 0;

  return E_NONE;
}


// read or write timestep
int MolfileInterface::timestep(float *coords, float *vels,
			       float *cell, double *simtime)
{
  if (!_plugin || !_dso || !_ptr)
    return 1;

  molfile_plugin_t *p = static_cast<molfile_plugin_t *>(_plugin);
  molfile_timestep_t *t = new molfile_timestep_t;
  int rv;
  
  if (_mode & M_WRITE) {
    t->coords = coords;
    t->velocities = vels;
    if (cell != NULL) {
      t->A = cell[0];
      t->B = cell[1];
      t->C = cell[2];
      t->alpha = cell[3];
      t->beta = cell[4];
      t->gamma = cell[5];
    } else {
      t->A = 0.0;
      t->B = 0.0;
      t->C = 0.0;
      t->alpha = 90.0;
      t->beta = 90.0;
      t->gamma = 90.0;
    }

    if (simtime)
      t->physical_time = *simtime;
    else
      t->physical_time = 0.0;

    rv = p->write_timestep(_ptr,t);

  } else {
    t->coords = coords;
    t->velocities = vels;
    rv = p->read_next_timestep(_ptr, _natoms, t);
    if (cell != NULL) {
      cell[0] = t->A;
      cell[1] = t->B;
      cell[2] = t->C;
      cell[3] = t->alpha;
      cell[4] = t->beta;
      cell[5] = t->gamma;
    }

    if (simtime)
      *simtime = t->physical_time;
  }

  return 0;
}

// set/get per type floating point property
int MolfileInterface::property(int propid, int *types, float *prop)
{
  if ((_info == NULL) || (types == NULL) || (prop == NULL))
    return P_NONE;

  molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);
  char buf[64];

  if (_mode & M_WSTRUCT) {

    for (int i=0; i < _natoms; ++i) {
      if (propid & P_OCCP)
	a[i].occupancy  = prop[types[i]];
      if (propid & P_BFAC)
	a[i].bfactor    = prop[types[i]];
      if (propid & P_MASS)
	a[i].mass       = prop[types[i]];
      if (propid & P_CHRG)
	a[i].charge     = prop[types[i]];
      if (propid & P_RADS)
	a[i].radius     = prop[types[i]];

      if (propid & (P_NAME|P_TYPE|P_RESN|P_SEGN|P_CHAI)) {
	sprintf(buf,"%g",prop[types[i]]);
	buf[15] = 0;

	if (propid & P_NAME)
	  strcpy(a[i].name,buf);
	if (propid & P_TYPE)
	  strcpy(a[i].type,buf);

	buf[7] = 0;
	if (propid & P_RESN)
	  strcpy(a[i].resname,buf);
	if (propid & P_SEGN)
	  strcpy(a[i].segid,buf);

	buf[1] = 0;
	if (propid & P_CHAI)
	  strcpy(a[i].chain,buf);
      }
    }
  } 
}


// set/get per type integer property
int MolfileInterface::property(int propid, int *types, int *prop) 
{
  if ((_info == NULL) || (types == NULL) || (prop == NULL))
    return P_NONE;

  molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);
  char buf[64];
  
  if (_mode & M_WSTRUCT) {

    for (int i=0; i < _natoms; ++i) {
      if (propid & P_RESI)
	a[i].resid = prop[types[i]];
      if (propid & P_ATMN)
	a[i].atomicnumber = prop[types[i]];
      if (propid & P_MASS)
	a[i].mass = prop[types[i]];
      if (propid & P_CHRG)
	a[i].charge = prop[types[i]];
      if (propid & P_RADS)
	a[i].radius = prop[types[i]];

      if (propid & (P_NAME|P_TYPE|P_RESN|P_SEGN|P_CHAI)) {
	sprintf(buf,"%d",prop[types[i]]);
	buf[15] = 0;

	if (propid & P_NAME)
	  strcpy(a[i].name,buf);
	if (propid & P_TYPE)
	  strcpy(a[i].type,buf);

	buf[7] = 0;
	if (propid & P_RESN)
	  strcpy(a[i].resname,buf);
	if (propid & P_SEGN)
	  strcpy(a[i].segid,buf);

	buf[1] = 0;
	if (propid & P_CHAI)
	  strcpy(a[i].chain,buf);
      }
    }
  } 
}

// set/get per type string property
int MolfileInterface::property(int propid, int *types, char **prop)
{

}

// set/get per atom floating point property
int MolfileInterface::property(int propid, float *prop)
{
}

// set/get per atom integer property
int MolfileInterface::property(int propid, int *prop)
{
  
}

// set/get per atom string property
int MolfileInterface::property(int propid, char **prop)
{
}


