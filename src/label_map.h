/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_LABEL_MAP_H
#define LMP_LABEL_MAP_H

#include "pointers.h"  // IWYU pragma: export

#include <unordered_map>

namespace LAMMPS_NS {

class LabelMap : protected Pointers {
 public:
   int natomtypes,nbondtypes,nangletypes,ndihedraltypes,nimpropertypes;
   std::vector<std::string> typelabel,btypelabel,atypelabel;
   std::vector<std::string> dtypelabel,itypelabel;
   std::unordered_map<std::string, int> typelabel_map;
   std::unordered_map<std::string, int> btypelabel_map;
   std::unordered_map<std::string, int> atypelabel_map;
   std::unordered_map<std::string, int> dtypelabel_map;
   std::unordered_map<std::string, int> itypelabel_map;

   // per-type data struct mapping this label map to another

   struct Lmap2Lmap {
     int *atom;
     int *bond;
     int *angle;
     int *dihedral;
     int *improper;
   };

   Lmap2Lmap lmap2lmap;

   LabelMap(LAMMPS *lmp);
   ~LabelMap();

   void allocate_type_labels();
   void modify_lmap(int, char **);                // labelmap command in the input script
   void merge_lmap(class LabelMap *, int);        // copy another lmap into this one
   void create_lmap2lmap(class LabelMap *, int);  // index mapping between two lmaps
   int find(const std::string &, int);            // find numeric type of type label
   int is_complete(int);                          // check if all types are assigned

   // input/output for atom class label map

   void write_data(FILE *);
   void read_restart(FILE *fp);
   void write_restart(FILE *);

 private:
   int find_or_create(const std::string &, std::vector<std::string> &, std::unordered_map<std::string, int> &); // look up type or create new type
   int search(const std::string &, const std::unordered_map<std::string, int> &); // look up type index
   char *read_string(FILE *);
   void write_string(const std::string &, FILE *);
   int read_int(FILE *);
};

}

#endif
