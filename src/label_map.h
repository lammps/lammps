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

#ifndef LMP_LABEL_MAP_H
#define LMP_LABEL_MAP_H

#include "pointers.h"  // IWYU pragma: export

namespace LAMMPS_NS {

class LabelMap : protected Pointers {
 public:
   int natomtypes,nbondtypes,nangletypes;
   int ndihedraltypes,nimpropertypes;
   std::string id;
   std::vector<std::string> typelabel,btypelabel,atypelabel;
   std::vector<std::string> dtypelabel,itypelabel;

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
   void modify_lmap(int, char **);
   void merge_lmap(class LabelMap *, int);                           // copy another lmap into this one
   void create_lmap2lmap(class LabelMap *, int);                     // index mapping between two lmaps
   int find_or_create(std::string, std::vector<std::string> &, int); // look up type or create new type
   int find(std::string, int);                                       // find numeric type of type label
   int search(std::string, std::vector<std::string>, int);           // look up type index
   int is_complete(int);                                             // check if all types are assigned

   // input/output for atom class label map

   void write_data(FILE *);
   void read_restart(FILE *fp);
   void write_restart(FILE *);

 private:
   int me;

   char *read_string(FILE *);
   void write_string(std::string, FILE *);
   int read_int(FILE *);
};

}

#endif

/* ERROR/WARNING messages:

E: Topology type exceeds system topology type

The number of bond, angle, etc types exceeds the system setting. See
the create_box or read_data command for how to specify these values.

*/
