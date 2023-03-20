// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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

#include "platform.h"
#include "tokenizer.h"
#include "utils.h"

#include "molfile_plugin.h"

#include <cstdio>
#include <cstring>
#include <cstdlib>

#if vmdplugin_ABIVERSION < 16
#error "unsupported VMD molfile plugin ABI version"
#endif

#define DEBUG 0

extern "C" {
  typedef int (*initfunc)();
  typedef int (*regfunc)(void *, vmdplugin_register_cb);
  typedef int (*finifunc)();

  typedef struct {
    void *p;
    const char *name;
  } plugin_reginfo_t;

  // callback function for plugin registration.
  static int plugin_register_cb(void *v, vmdplugin_t *p)
  {
    auto r = static_cast<plugin_reginfo_t *>(v);
    // make sure we have the proper plugin type (native reader)
    // for the desired file type (called "name" at this level)
    if ((strcmp(MOLFILE_PLUGIN_TYPE,p->type) == 0)
        && (strcmp(r->name, p->name) == 0)) {
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
    /* X  */ 0.00000f, 1.00794f, 4.00260f, 6.941f, 9.012182f, 10.811f,
    /* C  */ 12.0107f, 14.0067f, 15.9994f, 18.9984032f, 20.1797f,
    /* Na */ 22.989770f, 24.3050f, 26.981538f, 28.0855f, 30.973761f,
    /* S  */ 32.065f, 35.453f, 39.948f, 39.0983f, 40.078f, 44.955910f,
    /* Ti */ 47.867f, 50.9415f, 51.9961f, 54.938049f, 55.845f, 58.9332f,
    /* Ni */ 58.6934f, 63.546f, 65.409f, 69.723f, 72.64f, 74.92160f,
    /* Se */ 78.96f, 79.904f, 83.798f, 85.4678f, 87.62f, 88.90585f,
    /* Zr */ 91.224f, 92.90638f, 95.94f, 98.0f, 101.07f, 102.90550f,
    /* Pd */ 106.42f, 107.8682f, 112.411f, 114.818f, 118.710f, 121.760f,
    /* Te */ 127.60f, 126.90447f, 131.293f, 132.90545f, 137.327f,
    /* La */ 138.9055f, 140.116f, 140.90765f, 144.24f, 145.0f, 150.36f,
    /* Eu */ 151.964f, 157.25f, 158.92534f, 162.500f, 164.93032f,
    /* Er */ 167.259f, 168.93421f, 173.04f, 174.967f, 178.49f, 180.9479f,
    /* W  */ 183.84f, 186.207f, 190.23f, 192.217f, 195.078f, 196.96655f,
    /* Hg */ 200.59f, 204.3833f, 207.2f, 208.98038f, 209.0f, 210.0f, 222.0f,
    /* Fr */ 223.0f, 226.0f, 227.0f, 232.0381f, 231.03588f, 238.02891f,
    /* Np */ 237.0f, 244.0f, 243.0f, 247.0f, 247.0f, 251.0f, 252.0f, 257.0f,
    /* Md */ 258.0f, 259.0f, 262.0f, 261.0f, 262.0f, 266.0f, 264.0f, 269.0f,
    /* Mt */ 268.0f, 271.0f, 272.0f
  };

  /*
   * corresponding table of VDW radii.
   * van der Waals radii are taken from A. Bondi,
   * J. Phys. Chem., 68, 441 - 452, 1964,
   * except the value for H, which is taken from R.S. Rowland & R. Taylor,
   * J.Phys.Chem., 100, 7384 - 7391, 1996. Radii that are not available in
   * either of these publications have RvdW = 2.00 \AA.
   * The radii for Ions (Na, K, Cl, Ca, Mg, and Cs are based on the CHARMM27
   * Rmin/2 parameters for (SOD, POT, CLA, CAL, MG, CES) by default.
   */
  static const float pte_vdw_radius[] = {
    /* X  */ 1.5f, 1.2f, 1.4f, 1.82f, 2.0f, 2.0f,
    /* C  */ 1.7f, 1.55f, 1.52f, 1.47f, 1.54f,
    /* Na */ 1.36f, 1.18f, 2.0f, 2.1f, 1.8f,
    /* S  */ 1.8f, 2.27f, 1.88f, 1.76f, 1.37f, 2.0f,
    /* Ti */ 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
    /* Ni */ 1.63f, 1.4f, 1.39f, 1.07f, 2.0f, 1.85f,
    /* Se */ 1.9f, 1.85f, 2.02f, 2.0f, 2.0f, 2.0f,
    /* Zr */ 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
    /* Pd */ 1.63f, 1.72f, 1.58f, 1.93f, 2.17f, 2.0f,
    /* Te */ 2.06f, 1.98f, 2.16f, 2.1f, 2.0f,
    /* La */ 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
    /* Eu */ 2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
    /* Er */ 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
    /* W  */ 2.0f, 2.0f, 2.0f, 2.0f, 1.72f, 1.66f,
    /* Hg */ 1.55f, 1.96f, 2.02f, 2.0f, 2.0f, 2.0f, 2.0f,
    /* Fr */ 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 1.86f,
    /* Np */ 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
    /* Md */ 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
    /* Mt */ 2.0f, 2.0f, 2.0f
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

  static int get_pte_idx_from_string(const char *label) {
    int i, ind;
    char atom[3];

    if (label != nullptr) {
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

} // end of extern "C" region

using namespace LAMMPS_NS;

// constructor.
MolfileInterface::MolfileInterface(const char *type, const int mode)
  : _plugin(nullptr), _dso(nullptr), _ptr(nullptr), _info(nullptr), _natoms(0),
    _mode(mode), _caps(M_NONE), _props(0)
{
  _name = utils::strdup("none");
  _type = utils::strdup(type);
}

// destructor.
MolfileInterface::~MolfileInterface()
{
  forget_plugin();

  if (_info) {
    auto a = static_cast<molfile_atom_t *>(_info);
    delete[] a;
    _info = nullptr;
  }
  delete[] _name;
  delete[] _type;
}

// register the best matching plugin in a given directory
int MolfileInterface::find_plugin(const char *pluginpath)
{
  int retval = E_NONE;

  if (pluginpath == nullptr) return E_DIR;

  // search for suitable file names in provided path and try to inspect them
  // only look at .so files, since this is what VMD uses on all platforms

  for (const auto &dir : Tokenizer(pluginpath,":").as_vector()) {
    for (const auto &filename : platform::list_directory(dir)) {
      if (utils::strmatch(filename,"\\.so$")) {
        int rv = load_plugin(platform::path_join(dir,filename).c_str());
        if (rv > retval) retval = rv;
      }
    }
  }
  return retval;
}

// register the best matching plugin in a given directory
int MolfileInterface::load_plugin(const char *filename)
{
  void *dso;
  int retval = E_NONE;

  // access shared object
  dso = platform::dlopen(filename);
  if (dso == nullptr)
    return E_FILE;

  // check for required plugin symbols
  void *ifunc = platform::dlsym(dso,"vmdplugin_init");
  void *rfunc = platform::dlsym(dso,"vmdplugin_register");
  void *ffunc = platform::dlsym(dso,"vmdplugin_fini");
  if (ifunc == nullptr || rfunc == nullptr || ffunc == nullptr) {
    platform::dlclose(dso);
    return E_SYMBOL;
  }

  // initialize plugin. skip plugin if it fails.
  if (((initfunc)(ifunc))()) {
    platform::dlclose(dso);
    return E_SYMBOL;
  }

  // pre-register plugin.
  // the callback will be called for each plugin in the DSO and
  // check the file type. plugin->name will change if successful.
  plugin_reginfo_t reginfo;
  reginfo.p = nullptr;
  reginfo.name=_type;
  ((regfunc)rfunc)(&reginfo, plugin_register_cb);

  // make some checks to see if the plugin is suitable or not.
  auto plugin = static_cast<molfile_plugin_t *>(reginfo.p);

  // if the callback found a matching plugin and copied the struct,
  // its name element will point to a different location now.
  if (plugin == nullptr) {
    retval = E_TYPE;

    // check if the ABI matches the one used to compile this code
  } else if (plugin->abiversion != vmdplugin_ABIVERSION) {
    retval = E_ABI;

    // check if (basic) reading is supported
  } else if ((_mode & M_READ) &&
             ( (plugin->open_file_read == nullptr) ||
               (plugin->read_next_timestep  == nullptr) ||
               (plugin->close_file_read == nullptr) )) {
    retval = E_MODE;

    // check if (basic) writing is supported
  } else if ( (_mode & M_WRITE) &&
              ( (plugin->open_file_write == nullptr) ||
                (plugin->write_timestep  == nullptr) ||
                (plugin->close_file_write == nullptr) )) {
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
                (p->minorv >= plugin->minorv)) {
      retval = E_VERSION;
    }
  }

  // bingo! this one is a keeper.
  if (retval == E_NONE) {

    // make sure any existing plugin is wiped out
    forget_plugin();

    delete[] _name;
    _name = utils::strdup(fmt::format("{} v{}.{} by {}", plugin->prettyname,
                                      plugin->majorv, plugin->minorv, plugin->author));

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
  platform::dlclose(dso);
  return retval;
}

// deregister a plugin and close or reset all associated objects.
void MolfileInterface::forget_plugin()
{
  if (_ptr)
    close();

  if (_plugin)
    _plugin = nullptr;

  if (_dso) {
    void *ffunc = platform::dlsym(_dso,"vmdplugin_fini");
    if (ffunc)
      ((finifunc)ffunc)();
    platform::dlclose(_dso);
  }
  _dso = nullptr;

  delete[] _name;
  _name = utils::strdup("none");
  _caps = M_NONE;
}

// open file for reading or writing
int MolfileInterface::open(const char *name, int *natoms)
{
  if (!_plugin || !_dso || !natoms)
    return E_FILE;
  auto p = static_cast<molfile_plugin_t *>(_plugin);

  if (_mode & M_WRITE)
    _ptr = p->open_file_write(name,_type,*natoms);
  else if (_mode & M_READ)
    _ptr = p->open_file_read(name,_type,natoms);

  if (_ptr == nullptr)
    return E_FILE;

  _natoms = *natoms;
  // we need to deal with structure information,
  // so we allocate and initialize storage for it.
  if (_mode & (M_RSTRUCT|M_WSTRUCT)) {
    auto a = new molfile_atom_t[_natoms];
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
  auto p = static_cast<molfile_plugin_t *>(_plugin);

  int optflags = MOLFILE_NOOPTIONS;

  if (_mode & M_WSTRUCT) {
    optflags |= (_props & P_BFAC) ? MOLFILE_BFACTOR : 0;
    optflags |= (_props & P_OCCP) ? MOLFILE_OCCUPANCY : 0;
    optflags |= (_props & P_MASS) ? MOLFILE_MASS : 0;
    optflags |= (_props & P_CHRG) ? MOLFILE_CHARGE : 0;
    optflags |= (_props & P_RADS) ? MOLFILE_RADIUS : 0;
    optflags |= (_props & P_ATMN) ? MOLFILE_ATOMICNUMBER : 0;

    auto a = static_cast<molfile_atom_t *>(_info);
    p->write_structure(_ptr,optflags,a);
  } else if (_mode & M_RSTRUCT) {
    auto a = static_cast<molfile_atom_t *>(_info);
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

  auto p = static_cast<molfile_plugin_t *>(_plugin);

  if (_mode & M_WRITE) {
    p->close_file_write(_ptr);
  } else if (_mode & M_READ) {
    p->close_file_read(_ptr);
  }

  if (_info) {
    auto a = static_cast<molfile_atom_t *>(_info);
    delete[] a;
    _info = nullptr;
  }
  _ptr = nullptr;
  _natoms = 0;

  return E_NONE;
}


// read or write timestep
int MolfileInterface::timestep(float *coords, float *vels,
                               float *cell, double *simtime)
{
  if (!_plugin || !_dso || !_ptr)
    return 1;

  auto p = static_cast<molfile_plugin_t *>(_plugin);
  auto t = new molfile_timestep_t;
  int rv;

  if (_mode & M_WRITE) {
    t->coords = coords;
    t->velocities = vels;
    if (cell != nullptr) {
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
    if (coords == nullptr) {
      rv = p->read_next_timestep(_ptr, _natoms, nullptr);
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
      if (cell != nullptr) {
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

    if (rv == MOLFILE_EOF) {
      delete t;
      return 1;
    }
  }

  delete t;
  return 0;
}

// functions to read properties from molfile structure

#define PROPUPDATE(PROP,ENTRY,VAL)              \
  if (propid == (PROP)) { (VAL) = a.ENTRY; }

#define PROPSTRCPY(PROP,ENTRY,VAL)              \
  if (propid == (PROP)) { strcpy(VAL,a.ENTRY); }

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
  const char *prop = nullptr;
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
  if ((propid & (PROP)) == (PROP)) { a.ENTRY = VAL; plist |= (PROP); }

#define PROPSTRCPY(PROP,ENTRY,VAL)                                      \
  if ((propid & (PROP)) == (PROP)) { strcpy(a.ENTRY,VAL); plist |= (PROP); }

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
  if ((_info == nullptr) || (prop == nullptr) || (idx < 0) || (idx >= _natoms))
    return P_NONE;

  auto a = static_cast<molfile_atom_t *>(_info);

  if (_mode & M_WSTRUCT)
    _props |= write_atom_property(a[idx], propid, *prop);

  if (_mode & M_RSTRUCT)
    *prop = read_float_property(a[idx], propid);

  return _props;
}

// set/get per type floating point property
int MolfileInterface::property(int propid, int *types, float *prop)
{
  if ((_info == nullptr) || (types == nullptr) || (prop == nullptr))
    return P_NONE;

  auto a = static_cast<molfile_atom_t *>(_info);

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
  if ((_info == nullptr) || (prop == nullptr))
    return P_NONE;

  auto a = static_cast<molfile_atom_t *>(_info);

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
  if ((_info == nullptr) || (prop == nullptr) || (idx < 0) || (idx >= _natoms))
    return P_NONE;

  auto a = static_cast<molfile_atom_t *>(_info);

  if (_mode & M_WSTRUCT)
    return write_atom_property(a[idx], propid, *prop);

  if (_mode & M_RSTRUCT)
    *prop = static_cast<double>(read_float_property(a[idx], propid));

  return _props;
}

// set/get per type floating point property
int MolfileInterface::property(int propid, int *types, double *prop)
{
  if ((_info == nullptr) || (types == nullptr) || (prop == nullptr))
    return P_NONE;

  auto a = static_cast<molfile_atom_t *>(_info);

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
  if ((_info == nullptr) || (prop == nullptr))
    return P_NONE;

  auto a = static_cast<molfile_atom_t *>(_info);

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
  if ((_info == nullptr) || (prop == nullptr) || (idx < 0) || (idx >= _natoms))
    return P_NONE;

  auto a = static_cast<molfile_atom_t *>(_info);

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
  if ((_info == nullptr) || (types == nullptr) || (prop == nullptr))
    return P_NONE;

  auto a = static_cast<molfile_atom_t *>(_info);

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
  if ((_info == nullptr) || (prop == nullptr))
    return P_NONE;

  auto a = static_cast<molfile_atom_t *>(_info);

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
  if ((_info == nullptr) || (prop == nullptr) || (idx < 0) || (idx >= _natoms))
    return P_NONE;

  auto a = static_cast<molfile_atom_t *>(_info);

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
  if ((_info == nullptr) || (types == nullptr) || (prop == nullptr))
    return P_NONE;

  auto a = static_cast<molfile_atom_t *>(_info);

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
  if ((_info == nullptr) || (prop == nullptr))
    return P_NONE;

  auto a = static_cast<molfile_atom_t *>(_info);

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
