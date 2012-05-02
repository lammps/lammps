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

#ifndef LMP_MOLFILE_INTERFACE_H
#define LMP_MOLFILE_INTERFACE_H

namespace LAMMPS_NS {

// This class provides an abstract interface
// to the VMD plugin library.

class MolfileInterface 
{
 public:
  // plugin modes.
  enum {M_NONE=0,
	M_READ    = 1<<0, M_WRITE   = 1<<1,
	M_RSTRUCT = 1<<2, M_WSTRUCT = 1<<3,
	M_RBONDS  = 1<<4, M_WBONDS  = 1<<5,
	M_RANGLES = 1<<6, M_WANGLES = 1<<7,
	M_RVOL    = 1<<8, M_WVOL    = 1<<9};

  // plugin finder return values.
  enum {E_NONE=0,  //< nothing happened
	E_DIR,     //< path is not a directory or not readable
	E_FILE,    //< file not a DSO or not readable
	E_SYMBOL,  //< DSO is not a VMD plugin
	E_TYPE,    //< plugin is not of the correct type
	E_ABI,     //< plugin ABI does not match
	E_MODE,    //< plugin does not support desired mode
	E_VERSION, //< plugin is not newer as the current one
	E_MATCH    //< plugin matches
  };

  MolfileInterface(const char *type, const int mode);
  ~MolfileInterface();

  // disallowed default methods
 private:
  MolfileInterface() {};
  MolfileInterface(const MolfileInterface &) {};
  MolfileInterface &operator =(const MolfileInterface &) {return *this;};

 public:
  // search in the given directory path for a molfile plugin that
  // is of the right type and supports the desired mode.
  // if a plugin is already registered and a newer version is
  // found, this new version will override the old one.
  int find_plugin(const char *path);

  // try to register the plugin at given file name
  int load_plugin(const char *name);

  // deregister the current plugin/DSO and clean up.
  void forget_plugin();

  // return formatted string describing plugin
  char *get_plugin_name() const { return _name; };

  // return canonical plugin name
  char *get_plugin_type() const { return _type; };

  // open file through plugin
  int open(const char *name, int *natoms);

  // read/write timestep
  int timestep(float *coords, float *vels, float *cell, double *simtime);

  // close file managed by plugin
  int close();

  // inquire on interface status
  bool is_open() const { return (_ptr != 0); };

  // internal data
 protected:
  void *_plugin; // pointer to plugin struct
  void *_dso;    // handle to DSO
  void *_ptr;    // pointer to plugin data handle
  void *_info;   // pointer to atomic info data
  char *_type;   // canonical plugin name
  char *_name;   // plugin formatted name
  int   _natoms; // number of atoms
  int   _mode;   // plugin mode of operation
  int   _caps;   // plugin capabilities
};

}

#endif
