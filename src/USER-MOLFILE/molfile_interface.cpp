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
#include <ctype.h>

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

  /* periodic table of elements for translation of ordinal to atom type */
  static const char *pte_label[] = {
    "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P" , "S",  "Cl", "Ar", "K",  "Ca", "Sc",
    "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc",
    "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
    "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg"
  };
  static const int nr_pte_entries = sizeof(pte_label) / sizeof(char *);

  /* corresponding table of masses. */
  static const float pte_mass[] = {
    /* X  */ 0.00000, 1.00794, 4.00260, 6.941, 9.012182, 10.811,
    /* C  */ 12.0107, 14.0067, 15.9994, 18.9984032, 20.1797,
    /* Na */ 22.989770, 24.3050, 26.981538, 28.0855, 30.973761,
    /* S  */ 32.065, 35.453, 39.948, 39.0983, 40.078, 44.955910,
    /* Ti */ 47.867, 50.9415, 51.9961, 54.938049, 55.845, 58.9332,
    /* Ni */ 58.6934, 63.546, 65.409, 69.723, 72.64, 74.92160,
    /* Se */ 78.96, 79.904, 83.798, 85.4678, 87.62, 88.90585,
    /* Zr */ 91.224, 92.90638, 95.94, 98.0, 101.07, 102.90550,
    /* Pd */ 106.42, 107.8682, 112.411, 114.818, 118.710, 121.760,
    /* Te */ 127.60, 126.90447, 131.293, 132.90545, 137.327,
    /* La */ 138.9055, 140.116, 140.90765, 144.24, 145.0, 150.36,
    /* Eu */ 151.964, 157.25, 158.92534, 162.500, 164.93032,
    /* Er */ 167.259, 168.93421, 173.04, 174.967, 178.49, 180.9479,
    /* W  */ 183.84, 186.207, 190.23, 192.217, 195.078, 196.96655,
    /* Hg */ 200.59, 204.3833, 207.2, 208.98038, 209.0, 210.0, 222.0,
    /* Fr */ 223.0, 226.0, 227.0, 232.0381, 231.03588, 238.02891,
    /* Np */ 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0,
    /* Md */ 258.0, 259.0, 262.0, 261.0, 262.0, 266.0, 264.0, 269.0,
    /* Mt */ 268.0, 271.0, 272.0
  };

  /*
   * corresponding table of VDW radii.
   * van der Waals radii are taken from A. Bondi,
   * J. Phys. Chem., 68, 441 - 452, 1964,
   * except the value for H, which is taken from R.S. Rowland & R. Taylor,
   * J.Phys.Chem., 100, 7384 - 7391, 1996. Radii that are not available in
   * either of these publications have RvdW = 2.00 Å.
   * The radii for Ions (Na, K, Cl, Ca, Mg, and Cs are based on the CHARMM27
   * Rmin/2 parameters for (SOD, POT, CLA, CAL, MG, CES) by default.
   */
  static const float pte_vdw_radius[] = {
    /* X  */ 1.5, 1.2, 1.4, 1.82, 2.0, 2.0,
    /* C  */ 1.7, 1.55, 1.52, 1.47, 1.54,
    /* Na */ 1.36, 1.18, 2.0, 2.1, 1.8,
    /* S  */ 1.8, 2.27, 1.88, 1.76, 1.37, 2.0,
    /* Ti */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Ni */ 1.63, 1.4, 1.39, 1.07, 2.0, 1.85,
    /* Se */ 1.9, 1.85, 2.02, 2.0, 2.0, 2.0,
    /* Zr */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Pd */ 1.63, 1.72, 1.58, 1.93, 2.17, 2.0,
    /* Te */ 2.06, 1.98, 2.16, 2.1, 2.0,
    /* La */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Eu */ 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Er */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* W  */ 2.0, 2.0, 2.0, 2.0, 1.72, 1.66,
    /* Hg */ 1.55, 1.96, 2.02, 2.0, 2.0, 2.0, 2.0,
    /* Fr */ 2.0, 2.0, 2.0, 2.0, 2.0, 1.86,
    /* Np */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Md */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Mt */ 2.0, 2.0, 2.0
  };

  /* lookup functions */

  static const char *get_pte_label(const int idx)
  {
    if ((idx < 1) || (idx >= nr_pte_entries)) return pte_label[0];

    return pte_label[idx];
  }

  static float get_pte_mass(const int idx)
  {
    if ((idx < 1) || (idx >= nr_pte_entries)) return pte_mass[0];

    return pte_mass[idx];
  }

  static float get_pte_vdw_radius(const int idx)
  {
    if ((idx < 1) || (idx >= nr_pte_entries)) return pte_vdw_radius[0];

#if 1
    /* Replace with Hydrogen radius with an "all-atom" radius */
    if (idx == 1)
      return 1.0;    /* H  */
#else
    /* Replace with old VMD atom radii values */
    switch (idx) {
    case  1: return 1.0;    /* H  */
    case  6: return 1.5;    /* C  */
    case  7: return 1.4;    /* N  */
    case  8: return 1.3;    /* O  */
    case  9: return 1.2;    /* F  */
    case 15: return 1.5;    /* P  */
    case 16: return 1.9;    /* S  */
    }
#endif

    return pte_vdw_radius[idx];
  }

  static int get_pte_idx(const char *label)
  {
    int i;
    char atom[3];

    /* zap string */
    atom[0] = (char) 0;
    atom[1] = (char) 0;
    atom[2] = (char) 0;
    /* if we don't have a null-pointer, there must be at least two
     * chars, which is all we need. we convert to the capitalization
     * convention of the table above during assignment. */
    if (label != NULL) {
      atom[0] = (char) toupper((int) label[0]);
      atom[1] = (char) tolower((int) label[1]);
    }
    /* discard numbers in atom label */
    if (isdigit(atom[1])) atom[1] = (char) 0;

    for (i=0; i < nr_pte_entries; ++i) {
      if ( (pte_label[i][0] == atom[0])
           && (pte_label[i][1] == atom[1]) ) return i;
    }

    return 0;
  }

  static int get_pte_idx_from_string(const char *label) {
    int i, ind;
    char atom[3];

    if (label != NULL) {
      /* zap string */
      atom[0] = atom[1] = atom[2] = '\0';

      for (ind=0,i=0; (ind<2) && (label[i]!='\0'); i++) {
        if (label[i] != ' ') {
          atom[ind] = toupper(label[i]);
          ind++;
        }
      }

      if (ind < 1)
        return 0; /* no non-whitespace characters */

      for (i=0; i < nr_pte_entries; ++i) {
        if ((toupper(pte_label[i][0]) == atom[0]) && (toupper(pte_label[i][1]) == atom[1]))
          return i;
      }
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
      *next = '\0';
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
    if (plugin->read_next_timestep)      _caps |= M_READ;
    if (plugin->write_timestep)          _caps |= M_WRITE;
#if vmdplugin_ABIVERSION > 10
    // required to tell if velocities are present
    if (plugin->read_timestep_metadata)  _caps |= M_RVELS;
    // we can always offer velocities. we may not know if
    // they will be written by the plugin though.
    if (plugin->write_timestep)          _caps |= M_WVELS;
#endif
    if (plugin->read_structure)          _caps |= M_RSTRUCT;
    if (plugin->write_structure)         _caps |= M_WSTRUCT;
    if (plugin->read_bonds)              _caps |= M_RBONDS;
    if (plugin->write_bonds)             _caps |= M_WBONDS;
    if (plugin->read_angles)             _caps |= M_RANGLES;
    if (plugin->write_angles)            _caps |= M_WANGLES;
    if (plugin->read_volumetric_data)    _caps |= M_RVOL;
    if (plugin->write_volumetric_data)   _caps |= M_WVOL;

    if (_mode & M_WRITE) {
      _mode |= (_caps & M_WSTRUCT);
      _mode |= (_caps & M_WVELS);
    } else if (_mode & M_READ) {
      _mode |= (_caps & M_RSTRUCT);
      _mode |= (_caps & M_RVELS);
    }

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

// open file for reading or writing
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
  // we need to deal with structure information,
  // so we allocate and initialize storage for it.
  if (_mode & (M_RSTRUCT|M_WSTRUCT)) {
    molfile_atom_t *a = new molfile_atom_t[_natoms];
    _info = a;
    memset(_info,0,_natoms*sizeof(molfile_atom_t));
    for (int i=0; i < _natoms; ++i) {
      a[i].name[0] = 'X';
      a[i].type[0] = a[i].resname[0] = a[i].segid[0] = 'U';
      a[i].type[1] = a[i].resname[1] = a[i].segid[1] = 'N';
      a[i].type[2] = a[i].resname[2] = a[i].segid[2] = 'K';
      a[i].chain[0] = 'X';
    }
  }
  return E_NONE;
}

// get of set atom structure information
int MolfileInterface::structure()
{
  if (!_plugin || !_dso)
    return E_FILE;
  molfile_plugin_t *p = static_cast<molfile_plugin_t *>(_plugin);

  int optflags = MOLFILE_NOOPTIONS;

  if (_mode & M_WSTRUCT) {
    optflags |= (_props & P_BFAC) ? MOLFILE_BFACTOR : 0;
    optflags |= (_props & P_OCCP) ? MOLFILE_OCCUPANCY : 0;
    optflags |= (_props & P_MASS) ? MOLFILE_MASS : 0;
    optflags |= (_props & P_CHRG) ? MOLFILE_CHARGE : 0;
    optflags |= (_props & P_RADS) ? MOLFILE_RADIUS : 0;
    optflags |= (_props & P_ATMN) ? MOLFILE_ATOMICNUMBER : 0;

    molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);
    p->write_structure(_ptr,optflags,a);
  } else if (_mode & M_RSTRUCT) {
    molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);
    p->read_structure(_ptr,&optflags,a);
    // mandatory properties
    _props = P_NAME|P_TYPE|P_RESN|P_RESI|P_SEGN|P_CHAI;
    // optional properties
    _props |= (optflags & MOLFILE_BFACTOR) ? P_BFAC : 0;
    _props |= (optflags & MOLFILE_OCCUPANCY) ? P_OCCP : 0;
    _props |= (optflags & MOLFILE_MASS) ? P_MASS : 0;
    _props |= (optflags & MOLFILE_CHARGE) ? P_CHRG : 0;
    _props |= (optflags & MOLFILE_RADIUS) ? P_RADS : 0;
    _props |= (optflags & MOLFILE_ATOMICNUMBER) ? P_ATMN : 0;
  }
  return 0;
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
      t->A = 0.0f;
      t->B = 0.0f;
      t->C = 0.0f;
      t->alpha = 90.0f;
      t->beta = 90.0f;
      t->gamma = 90.0f;
    }

    if (simtime)
      t->physical_time = *simtime;
    else
      t->physical_time = 0.0;

    rv = p->write_timestep(_ptr,t);

  } else {
    // no coordinate storage => skip step
    if (coords == NULL) {
      rv = p->read_next_timestep(_ptr, _natoms, NULL);
    } else {
      t->coords = coords;
      t->velocities = vels;
      t->A = 0.0f;
      t->B = 0.0f;
      t->C = 0.0f;
      t->alpha = 90.0f;
      t->beta = 90.0f;
      t->gamma = 90.0f;
      t->physical_time = 0.0;
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

    if (rv == MOLFILE_EOF)
      return 1;
  }

  return 0;
}

// functions to read properties from molfile structure

#define PROPUPDATE(PROP,ENTRY,VAL)              \
  if (propid == PROP) { VAL = a.ENTRY; }

#define PROPSTRCPY(PROP,ENTRY,VAL)              \
  if (propid == PROP) { strcpy(VAL,a.ENTRY); }

// single precision floating point props
static float read_float_property(molfile_atom_t &a, const int propid)
{
  float prop = 0.0f;
  int iprop = 0;
  PROPUPDATE(MolfileInterface::P_OCCP,occupancy,prop);
  PROPUPDATE(MolfileInterface::P_BFAC,bfactor,prop);
  PROPUPDATE(MolfileInterface::P_MASS,mass,prop);
  PROPUPDATE(MolfileInterface::P_CHRG,charge,prop);
  PROPUPDATE(MolfileInterface::P_RADS,radius,prop);

  PROPUPDATE((MolfileInterface::P_ATMN|MolfileInterface::P_MASS),
             atomicnumber,iprop);
  PROPUPDATE((MolfileInterface::P_ATMN|MolfileInterface::P_RADS),
             atomicnumber,iprop);
  if (propid & MolfileInterface::P_ATMN) {
    if (propid & MolfileInterface::P_MASS)
      prop = get_pte_mass(iprop);
    if (propid & MolfileInterface::P_RADS)
      prop = get_pte_vdw_radius(iprop);
  }

  return prop;
}

// integer and derived props
static int read_int_property(molfile_atom_t &a, const int propid)
{
  int prop = 0;
  const char * sprop;

  PROPUPDATE(MolfileInterface::P_RESI,resid,prop);
  PROPUPDATE(MolfileInterface::P_ATMN,atomicnumber,prop);

  PROPUPDATE((MolfileInterface::P_ATMN|MolfileInterface::P_NAME),
             name,sprop);
  PROPUPDATE((MolfileInterface::P_ATMN|MolfileInterface::P_TYPE),
             type,sprop);

  if (propid & MolfileInterface::P_ATMN) {
    if (propid & (MolfileInterface::P_NAME|MolfileInterface::P_TYPE))
      prop = get_pte_idx_from_string(sprop);
  }

  return prop;
}

// string and derived props
static const char *read_string_property(molfile_atom_t &a,
                                        const int propid)
{
  const char *prop = NULL;
  int iprop = 0;
  PROPUPDATE(MolfileInterface::P_NAME,name,prop);
  PROPUPDATE(MolfileInterface::P_TYPE,type,prop);
  PROPUPDATE(MolfileInterface::P_RESN,resname,prop);
  PROPUPDATE(MolfileInterface::P_SEGN,segid,prop);

  PROPUPDATE((MolfileInterface::P_ATMN|MolfileInterface::P_NAME),
             atomicnumber,iprop);
  PROPUPDATE((MolfileInterface::P_ATMN|MolfileInterface::P_TYPE),
             atomicnumber,iprop);

  if (propid & MolfileInterface::P_ATMN) {
    if (propid & (MolfileInterface::P_NAME|MolfileInterface::P_TYPE))
      prop = get_pte_label(iprop);
  }

  return prop;
}
#undef PROPUPDATE
#undef PROPSTRCPY

// functions to store properties into molfile structure

#define PROPUPDATE(PROP,ENTRY,VAL)                                  \
  if ((propid & PROP) == PROP) { a.ENTRY = VAL; plist |= PROP; }

#define PROPSTRCPY(PROP,ENTRY,VAL)                                      \
  if ((propid & PROP) == PROP) { strcpy(a.ENTRY,VAL); plist |= PROP; }

// floating point props
static int write_atom_property(molfile_atom_t &a,
                               const int propid,
                               const float prop)
{
  int plist = MolfileInterface::P_NONE;
  PROPUPDATE(MolfileInterface::P_OCCP,occupancy,prop);
  PROPUPDATE(MolfileInterface::P_BFAC,bfactor,prop);
  PROPUPDATE(MolfileInterface::P_MASS,mass,prop);
  PROPUPDATE(MolfileInterface::P_CHRG,charge,prop);
  PROPUPDATE(MolfileInterface::P_RADS,radius,prop);
  return plist;
}

// double precision floating point props
static int write_atom_property(molfile_atom_t &a,
                               const int propid,
                               const double prop)
{
  return write_atom_property(a,propid,static_cast<float>(prop));
}

// integer and derived props
static int write_atom_property(molfile_atom_t &a,
                               const int propid,
                               const int prop)
{
  int plist = MolfileInterface::P_NONE;
  PROPUPDATE(MolfileInterface::P_RESI,resid,prop);
  PROPUPDATE(MolfileInterface::P_ATMN,atomicnumber,prop);
  PROPUPDATE((MolfileInterface::P_ATMN|MolfileInterface::P_MASS),
             mass,get_pte_mass(prop));
  PROPSTRCPY((MolfileInterface::P_ATMN|MolfileInterface::P_NAME),
             name,get_pte_label(prop));
  PROPSTRCPY((MolfileInterface::P_ATMN|MolfileInterface::P_TYPE),
             type,get_pte_label(prop));
  return plist;
}

// integer and derived props
static int write_atom_property(molfile_atom_t &a,
                               const int propid,
                               const char *prop)
{
  int plist = MolfileInterface::P_NONE;
  PROPSTRCPY(MolfileInterface::P_NAME,name,prop);
  PROPSTRCPY(MolfileInterface::P_TYPE,type,prop);
  PROPSTRCPY(MolfileInterface::P_RESN,resname,prop);
  PROPSTRCPY(MolfileInterface::P_SEGN,segid,prop);
  return plist;
}
#undef PROPUPDATE
#undef PROPSTRCPY

// set/get atom floating point property
int MolfileInterface::property(int propid, int idx, float *prop)
{
  if ((_info == NULL) || (prop == NULL) || (idx < 0) || (idx >= _natoms))
    return P_NONE;

  molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);

  if (_mode & M_WSTRUCT)
    _props |= write_atom_property(a[idx], propid, *prop);

  if (_mode & M_RSTRUCT)
    *prop = read_float_property(a[idx], propid);

  return _props;
}

// set/get per type floating point property
int MolfileInterface::property(int propid, int *types, float *prop)
{
  if ((_info == NULL) || (types == NULL) || (prop == NULL))
    return P_NONE;

  molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);

  if (_mode & M_WSTRUCT) {
    for (int i=0; i < _natoms; ++i)
      _props |= write_atom_property(a[i], propid, prop[types[i]]);
  }

  // useless for reading.
  if (_mode & M_RSTRUCT)
    return P_NONE;

  return _props;
}

// set/get per atom floating point property
int MolfileInterface::property(int propid, float *prop)
{
  if ((_info == NULL) || (prop == NULL))
    return P_NONE;

  molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);

  if (_mode & M_WSTRUCT) {
    for (int i=0; i < _natoms; ++i)
      _props |= write_atom_property(a[i], propid, prop[i]);
  }

  if (_mode & M_RSTRUCT) {
    for (int i=0; i < _natoms; ++i)
      prop[i] = read_float_property(a[i], propid);
  }

  return _props;
}

// set/get atom floating point property
int MolfileInterface::property(int propid, int idx, double *prop)
{
  if ((_info == NULL) || (prop == NULL) || (idx < 0) || (idx >= _natoms))
    return P_NONE;

  molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);

  if (_mode & M_WSTRUCT)
    return write_atom_property(a[idx], propid, *prop);

  if (_mode & M_RSTRUCT)
    *prop = static_cast<double>(read_float_property(a[idx], propid));

  return _props;
}

// set/get per type floating point property
int MolfileInterface::property(int propid, int *types, double *prop)
{
  if ((_info == NULL) || (types == NULL) || (prop == NULL))
    return P_NONE;

  molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);

  if (_mode & M_WSTRUCT) {
    for (int i=0; i < _natoms; ++i)
      _props |= write_atom_property(a[i], propid, prop[types[i]]);
  }

  // useless for reading
  if (_mode & M_RSTRUCT)
    return P_NONE;

  return _props;
}

// set/get per atom floating point property
int MolfileInterface::property(int propid, double *prop)
{
  if ((_info == NULL) || (prop == NULL))
    return P_NONE;

  molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);

  if (_mode & M_WSTRUCT) {
    for (int i=0; i < _natoms; ++i)
      _props |= write_atom_property(a[i], propid, prop[i]);
  }
  if (_mode & M_RSTRUCT) {
    for (int i=0; i < _natoms; ++i)
      prop[i] = static_cast<double>(read_float_property(a[i], propid));
  }

  return _props;
}

#define INT_TO_STRING_BODY(IDX)                         \
  buf[15] = 0;                                          \
  if (propid & P_NAME)                                  \
    _props |= write_atom_property(a[IDX],P_NAME,buf);   \
  if (propid & P_TYPE)                                  \
    _props |= write_atom_property(a[IDX],P_TYPE,buf);   \
  buf[7] = 0;                                           \
  if (propid & P_RESN)                                  \
    _props |= write_atom_property(a[IDX],P_RESN,buf);   \
  if (propid & P_SEGN)                                  \
    _props |= write_atom_property(a[IDX],P_SEGN,buf);   \
  buf[1] = 0;                                           \
  if (propid & P_CHAI)                                  \
    _props |= write_atom_property(a[IDX],P_CHAI,buf)

// set/get atom integer property
int MolfileInterface::property(int propid, int idx, int *prop)
{
  if ((_info == NULL) || (prop == NULL) || (idx < 0) || (idx >= _natoms))
    return P_NONE;

  molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);

  if (_mode & M_WSTRUCT) {
    char buf[64];

    _props |= write_atom_property(a[idx], propid, *prop);

    if (propid & (P_NAME|P_TYPE|P_RESN|P_SEGN|P_CHAI)) {
      sprintf(buf,"%d",*prop);
      INT_TO_STRING_BODY(idx);
    }
  }

  if (_mode & M_RSTRUCT)
    *prop = read_int_property(a[idx], propid);

  return _props;
}

// set/get per type integer property
int MolfileInterface::property(int propid, int *types, int *prop)
{
  if ((_info == NULL) || (types == NULL) || (prop == NULL))
    return P_NONE;

  molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);

  if (_mode & M_WSTRUCT) {
    char buf[64];

    for (int i=0; i < _natoms; ++i)
      _props |= write_atom_property(a[i], propid, prop[types[i]]);

    if (propid & (P_NAME|P_TYPE|P_RESN|P_SEGN|P_CHAI)) {
      for (int i=0; i < _natoms; ++i) {
        sprintf(buf,"%d",prop[types[i]]);
        INT_TO_STRING_BODY(i);
      }
    }
  }

  // useless when reading
  if (_mode & M_RSTRUCT)
    return P_NONE;

  return _props;
}

// set/get per atom integer property
int MolfileInterface::property(int propid, int *prop)
{
  if ((_info == NULL) || (prop == NULL))
    return P_NONE;

  molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);

  if (_mode & M_WSTRUCT) {
    char buf[64];

    for (int i=0; i < _natoms; ++i)
      _props |= write_atom_property(a[i],propid,prop[i]);

    if (propid & (P_NAME|P_TYPE|P_RESN|P_SEGN|P_CHAI)) {
      for (int i=0; i < _natoms; ++i) {
        sprintf(buf,"%d",prop[i]);
        INT_TO_STRING_BODY(i);
      }
    }
  }

  if (_mode & M_RSTRUCT) {
    for (int i=0; i < _natoms; ++i)
      prop[i] = read_int_property(a[i], propid);
  }

  return _props;
}
#undef INT_TO_STRING_BODY

// set/get atom string property
int MolfileInterface::property(int propid, int idx, char *prop)
{
  if ((_info == NULL) || (prop == NULL) || (idx < 0) || (idx >= _natoms))
    return P_NONE;

  molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);

  if (_mode & M_WSTRUCT) {
    _props |= write_atom_property(a[idx], propid, prop);
  }

  if (_mode & M_RSTRUCT)
    strcpy(prop,read_string_property(a[idx], propid));

  return _props;
}

// set/get per type string property
int MolfileInterface::property(int propid, int *types, char **prop)
{
  if ((_info == NULL) || (types == NULL) || (prop == NULL))
    return P_NONE;

  molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);

  if (_mode & M_WSTRUCT) {
    for (int i=0; i < _natoms; ++i) {
      _props |= write_atom_property(a[i], propid, prop[types[i]]);
    }
  }

  // useless when reading
  if (_mode & M_RSTRUCT)
    return P_NONE;

  return _props;
}

// set/get per atom string property
int MolfileInterface::property(int propid, char **prop)
{
  if ((_info == NULL) || (prop == NULL))
    return P_NONE;

  molfile_atom_t *a = static_cast<molfile_atom_t *>(_info);

  if (_mode & M_WSTRUCT) {
    for (int i=0; i < _natoms; ++i) {
      _props |= write_atom_property(a[i], propid, prop[i]);
    }
  }

  // not supported right now. XXX: should we use strdup() here?
  if (_mode & M_RSTRUCT)
    return P_NONE;

  return _props;
}
