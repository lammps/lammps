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
   Contributing author: Axel Kohlmeyer (Temple U)
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
        M_READ    = 1<<0,  M_WRITE   = 1<<1,
        M_RSTRUCT = 1<<2,  M_WSTRUCT = 1<<3,
        M_RBONDS  = 1<<4,  M_WBONDS  = 1<<5,
        M_RANGLES = 1<<6,  M_WANGLES = 1<<7,
        M_RVOL    = 1<<8,  M_WVOL    = 1<<9,
        M_RVELS   = 1<<10, M_WVELS   = 1<<11,
        M_LAST    = 1<<12 };

  // plugin finder return values.
  enum {E_NONE=0,  //< nothing happened
        E_DIR,     //< path is not a directory or not readable
        E_FILE,    //< file not a DSO or not readable
        E_SYMBOL,  //< DSO is not a VMD plugin
        E_TYPE,    //< plugin is not of the correct type
        E_ABI,     //< plugin ABI does not match
        E_MODE,    //< plugin does not support desired mode
        E_VERSION, //< plugin is not newer as the current one
        E_MATCH,   //< plugin matches
        E_LAST     //< last entry
  };

  // atom structure properties. deliberately not complete.
  enum {P_NONE=0,     //< no structure information available
        P_NAME=1<<0,  //< atom name,    char[16]
        P_TYPE=1<<1,  //< atom type,    char[16]
        P_RESN=1<<2,  //< residue name, char[ 8]
        P_RESI=1<<3,  //< residue index, int
        P_SEGN=1<<4,  //< segment name, char[ 8]
        P_CHAI=1<<5,  //< chain id,     char[ 2]
        P_OCCP=1<<6,  //< occupancy,    float
        P_BFAC=1<<7,  //< B factor,     float
        P_MASS=1<<8,  //< atom mass,    float
        P_CHRG=1<<9,  //< atom charge,  float
        P_RADS=1<<10, //< atom radius,  float
        P_ATMN=1<<11, //< atomic number, int
        P_LAST=1<<12  //< last entry
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
  // return canonical plugin name (= file type)
  char *get_plugin_type() const { return _type; };

  // file operations

  // open file through plugin
  int open(const char *name, int *natoms);
  // read/write structure info
  int structure();
  // read/write timestep
  int timestep(float *coords, float *vels, float *cell, double *simtime);
  // close file managed by plugin
  int close();

  // inquire on interface status

  // true if file stream is active.
  bool is_open() const { return (_ptr != 0); };
  // true if file format requires or provides atom properties
  bool has_props() const { return (_mode & (M_RSTRUCT|M_WSTRUCT)) != 0; };
  // true if file format can read or write velocities
  bool has_vels() const { return (_mode & (M_RVELS|M_WVELS)) != 0; };

  // return number of atoms in current file. -1 if closed/invalid;
  bool get_natoms() const { return _natoms; };
  // return property bitmask
  bool get_props() const { return _props; };

  // atom property operations

  // set/get atom floating point property
  int property(int propid, int idx, float *prop);
  // set/get per type floating point property
  int property(int propid, int *types, float *prop);
  // set/get per atom floating point property
  int property(int propid, float *prop);
  // set/get atom floating point property
  int property(int propid, int idx, double *prop);
  // set/get per type floating point property
  int property(int propid, int *types, double *prop);
  // set/get per atom floating point property
  int property(int propid, double *prop);
  // set/get atom integer property
  int property(int propid, int idx, int *prop);
  // set/get per type integer property
  int property(int propid, int *types, int *prop);
  // set/get per atom integer property
  int property(int propid, int *prop);
  // set/get atom string property
  int property(int propid, int idx, char *prop);
  // set/get per type string property
  int property(int propid, int *types, char **prop);
  // set/get per atom string property
  int property(int propid, char **prop);

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
  int   _props;  // accumulated/available properties
};

}

#endif
