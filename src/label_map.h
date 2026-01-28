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

/*! \file label_map.h */

#ifndef LMP_LABEL_MAP_H
#define LMP_LABEL_MAP_H

#include "pointers.h"    // IWYU pragma: export

#include <unordered_map>

namespace LAMMPS_NS {

/*! \class LabelMap
 *  \brief Manage type labels for atoms, bonds, angles, dihedrals, and impropers
 *
 * The LabelMap class provides functionality to map between string labels and
 * numeric type indices for atoms, bonds, angles, dihedrals, and impropers in LAMMPS.
 * This enables users to reference types by symbolic names (e.g., "C", "H", "C-H")
 * instead of numeric indices, improving readability and maintainability of input scripts.
 *
 * Type labels for bonded interactions *may* (but are not required to) use a
 * hyphen-delimited format indicating the types of the constituent atoms.  Examples:
 * - Bond types: "atom1-atom2" (e.g., "C-H", "N-O")
 * - Angle types: "atom1-atom2-atom3" (e.g., "H-C-H", "C-N-C")
 * - Dihedral types: "atom1-atom2-atom3-atom4" (e.g., "C-C-N-H")
 * - Improper types: "atom1-atom2-atom3-atom4" (e.g., "C-N-C-C")
 *
 * The class supports bidirectional lookup (label <-> type) and can infer bonded
 * interaction types from constituent atom types when using a hyphen-delimited
 * format convention. */

class LabelMap : protected Pointers {
  friend class AtomVec;
  friend class DumpCustom;
  friend class DumpExtXYZ;
  friend class DumpXYZ;
  friend class ReadData;

 public:
  /*! Construct a LabelMap instance
   *
   * \param  lmp              Pointer to LAMMPS instance
   * \param  natomtypes       Number of atom types in map
   * \param  nbondtypes       Number of bond types in map
   * \param  nangletypes      Number of angle types in map
   * \param  ndihedraltypes   Number of dihedral types in map
   * \param  nimpropertypes   Number of improper types in map */
  LabelMap(LAMMPS *lmp, int, int, int, int, int);
  ~LabelMap() override;

  /*! Process labelmap command from input script
   *
   \verbatim embed:rst

Add or modify type label mappings from the LAMMPS
:doc:`labelmap <labelmap>` input command.

   \endverbatim
   *
   * \param  narg  Number of arguments
   * \param  arg   Array of argument strings */
  void modify_lmap(int, char **);

  /*! Copy another LabelMap into this one
   *
   \verbatim embed:rst

Merge type labels from another LabelMap instance into the current one.
Currently used when combining data from multiple sources with
:doc:`read_data add <read_data>` or when replicating the system with
:doc:`replicate <replicate>`.

   \endverbatim
   *
   * \param  lmap  Pointer to source LabelMap
   * \param  mode  Merge mode flag */
  void merge_lmap(LabelMap *, int);

  /*! Create index mapping between two LabelMaps
   *
   * Build a mapping structure (lmap2lmap) that translates type indices
   * from another LabelMap to the current one based on matching labels.
   *
   * \param  lmap  Pointer to source LabelMap
   * \param  mode  Mapping mode flag */
  void create_lmap2lmap(LabelMap *, int);

  /*! Find numeric type from type label
   *
   * Look up the numeric type index corresponding to a type label string.
   *
   * \param  mylabel  Type label string to search for
   * \param  mode     Type category: Atom::ATOM, Atom::BOND, Atom::ANGLE,
   *                  Atom::DIHEDRAL, or Atom::IMPROPER
   * \return          Numeric type index (1-based), or -1 if not found */
  int find_type(const std::string &, int) const;

  /*! Find type label from numeric type
   *
   * Reverse lookup: retrieve the type label string for a given numeric type.
   *
   * \param  i     Numeric type index (1-based)
   * \param  mode  Type category: Atom::ATOM, Atom::BOND, Atom::ANGLE,
   *               Atom::DIHEDRAL, or Atom::IMPROPER
   * \return       Reference to type label string, or empty string if not found */
  const std::string &find_label(int, int) const;

  /*! Check if all types have assigned labels
   *
   * Verify that every type in the specified category has a corresponding label.
   *
   * \param  mode  Type category: Atom::ATOM, Atom::BOND, Atom::ANGLE,
   *               Atom::DIHEDRAL, or Atom::IMPROPER
   * \return       True if all types have labels, false otherwise */
  bool is_complete(int) const;

  /*! \name Interaction type inference from hyphen-delimited labels
   *
   * These methods infer bonded interaction types (bonds, angles, dihedrals, impropers)
   * from constituent atom types using a hyphen-delimited format.  This requires that
   * type labels for bonded interactions were entered following this convention.
   * The inference functions consider the symmetry of the interaction and thus atom
   * types may be swapped accordingly and the bonded type will still be matched.
   * @{ */

  /*! Infer bond type from two numeric atom types
   *
   * Look up or create a bond type from two atom type indices by constructing
   * a hyphen-delimited label (e.g., "C-H") and searching the bond type labels.
   *
   * \param  atype1  First atom type index
   * \param  atype2  Second atom type index
   * \return         Bond type index, or -1 if not found */
  int infer_bondtype(int, int);

  /*! Infer bond type from atom type labels
   *
   * \overload
   *
   * Look up a bond type from two atom type labels.
   *
   * \param  labels  Vector of two atom type label strings
   * \return         Bond type index, or -1 if not found */
  int infer_bondtype(const std::vector<std::string> &);

  /*! Infer angle type from three numeric atom types
   *
   * Look up or create an angle type from three atom type indices by
   * constructing a hyphen-delimited label (e.g., "H1-C1-H2").
   *
   * \param  atype1  First atom type index
   * \param  atype2  Second atom type index (center atom)
   * \param  atype3  Third atom type index
   * \return         Angle type index, or -1 if not found */
  int infer_angletype(int, int, int);

  /*! Infer angle type from three atom type labels
   *
   * \overload
   *
   * Look up an angle type from three atom type labels.
   *
   * \param  labels  Vector of three atom type label strings
   * \return         Angle type index, or -1 if not found */
  int infer_angletype(const std::vector<std::string> &);

  /*! Infer dihedral type from four numeric atom types
   *
   * Look up a dihedral type from four atom type indices by
   * constructing a hyphen-delimited label (e.g., "C-C-N-H").
   *
   * \param  atype1  First atom type index
   * \param  atype2  Second atom type index
   * \param  atype3  Third atom type index
   * \param  atype4  Fourth atom type index
   * \return         Dihedral type index, or -1 if not found */
  int infer_dihedraltype(int, int, int, int);

  /*! Infer dihedral type from atom type labels
   *
   * \overload
   *
   * Look up a dihedral type from four atom type labels.
   *
   * \param  labels  Vector of four atom type label strings
   * \return         Dihedral type index, or -1 if not found */
  int infer_dihedraltype(const std::vector<std::string> &);

  /*! Infer improper type from four numeric atom types
   *
   * Look up an improper type from four atom type indices by
   * constructing a hyphen-delimited label (e.g., "C-N-C-C").
   *
   * \param  atype1  First atom type index (center atom)
   * \param  atype2  Second atom type index
   * \param  atype3  Third atom type index
   * \param  atype4  Fourth atom type index
   * \return         Improper type index, or -1 if not found */
  int infer_impropertype(int, int, int, int);

  /*! Infer improper type from atom type labels
   *
   * \overload
   *
   * Look up an improper type from four atom type labels.
   *
   * \param  labels  Vector of four atom type label strings
   * \return         Improper type index, or -1 if not found */
  int infer_impropertype(const std::vector<std::string> &);

  /*! @} */

  /*! Parse hyphen-delimited type label into components
   *
   * Split a hyphen-delimited label (e.g., "C-N-H") into individual type strings.
   * Validates that the number of components matches the expected count.
   *
   * \param  ntypes  Expected number of components
   * \param  label   Hyphen-delimited label string
   * \param  types   Output vector to store component strings
   * \return         0 on success, -1 if component count doesn't match ntypes */
  int parse_typelabel(int, const std::string &, std::vector<std::string> &);

  /*! \name I/O methods for label map persistence
   * @{ */

  /*! Write label map to data file
   *
   * Output all type labels as sections to a LAMMPS data file.
   *
   * \param  fp  File pointer for writing */
  void write_data(FILE *);

  /*! Read label map from restart file
   *
   * Restore label map data from a LAMMPS restart file.
   *
   * \param  fp  File pointer for reading */
  void read_restart(FILE *fp);

  /*! Write label map to restart file
   *
   * Save label map data to a LAMMPS restart file for later restoration.
   *
   * \param  fp  File pointer for writing */
  void write_restart(FILE *);

  /*! @} */

 protected:
  int natomtypes, nbondtypes, nangletypes, ndihedraltypes, nimpropertypes;    //!< Type counts
  std::vector<std::string> typelabel, btypelabel,
      atypelabel;                                     //!< Label storage (atoms, bonds, angles)
  std::vector<std::string> dtypelabel, itypelabel;    //!< Label storage (dihedrals, impropers)

  std::unordered_map<std::string, int> typelabel_map;     //!< Atom label → type mapping
  std::unordered_map<std::string, int> btypelabel_map;    //!< Bond label → type mapping
  std::unordered_map<std::string, int> atypelabel_map;    //!< Angle label → type mapping
  std::unordered_map<std::string, int> dtypelabel_map;    //!< Dihedral label → type mapping
  std::unordered_map<std::string, int> itypelabel_map;    //!< Improper label → type mapping

  /*! \struct Lmap2Lmap
   *  \brief Mapping structure between two LabelMaps
   *
   * Stores per-type index mappings from another LabelMap to this one,
   * enabling type translation when merging or comparing different label maps. */

  struct Lmap2Lmap {
    int *atom;        //!< Atom type mapping array
    int *bond;        //!< Bond type mapping array
    int *angle;       //!< Angle type mapping array
    int *dihedral;    //!< Dihedral type mapping array
    int *improper;    //!< Improper type mapping array
  };

  Lmap2Lmap lmap2lmap;    //!< Instance of inter-map translation data

  void reset_type_labels();    //!< Clear all type labels
  int find_or_create(
      const std::string &, std::vector<std::string> &,
      std::unordered_map<std::string, int> &);    //!< Look up type or create new type
  int search(const std::string &,
             const std::unordered_map<std::string, int> &) const;    //!< Look up type index
  char *read_string(FILE *);                         //!< Read string from binary file
  void write_string(const std::string &, FILE *);    //!< Write string to binary file
  int read_int(FILE *);                              //!< Read integer from binary file

  void write_map(const std::string &);    //!< Write label map to file for debugging
};

}    // namespace LAMMPS_NS

#endif
