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

#include <array>
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
  LabelMap(LAMMPS *lmp, int natomtypes, int nbondtypes, int nangletypes, int ndihedraltypes,
           int nimpropertypes);
  ~LabelMap() override;

  int checkflag;    //!< Flag to check for self-consistent type labels

  /*! Process labelmap command from input script
   *
   \verbatim embed:rst

Add or modify type label mappings from the LAMMPS
:doc:`labelmap <labelmap>` input command.

   \endverbatim
   *
   * \param  narg  Number of arguments
   * \param  arg   Array of argument strings */
  void modify_lmap(int narg, char **arg);

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
  void merge_lmap(LabelMap *lmap, int mode);

  /*! Create index mapping between two LabelMaps
   *
   * Build a mapping structure (lmap2lmap) that translates type indices
   * from another LabelMap to the current one based on matching labels.
   *
   * \param  lmap  Pointer to source LabelMap
   * \param  mode  Mapping mode flag */
  void create_lmap2lmap(LabelMap *lmap, int mode);

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
   * \return         Bond type index if types match in the specified order,
   *                 negative bond type index if types match in reverse order,
   *                 0 if there is no match found */
  int infer_bondtype(int atype1, int atype2);

  /*! Infer bond type from atom type labels
   *
   * \overload
   *
   * Look up a bond type from two atom type labels.
   *
   * \param  labels  Vector of two atom type label strings
   * \return         Bond type index if types match in the specified order,
   *                 negative bond type index if types match in reverse order,
   *                 0 if there is no match found */
  int infer_bondtype(const std::vector<std::string> &labels);

  /*! Infer angle type from three numeric atom types
   *
   * Look up or create an angle type from three atom type indices by
   * constructing a hyphen-delimited label (e.g., "H1-C1-H2").
   *
   * \param  atype1  First atom type index
   * \param  atype2  Second atom type index (center atom)
   * \param  atype3  Third atom type index
   * \return         Angle type index if types match in the specified order,
   *                 negative angle type index if types match in reverse order,
   *                 0 if there is no match found */
  int infer_angletype(int atype1, int atype2, int atype3);

  /*! Infer angle type from three atom type labels
   *
   * \overload
   *
   * Look up an angle type from three atom type labels.
   *
   * \param  labels  Vector of three atom type label strings
   * \return         Angle type index if types match in the specified order,
   *                 negative angle type index if types match in reverse order,
   *                 0 if there is no match found */
  int infer_angletype(const std::vector<std::string> &labels);

  /*! Infer dihedral type from four numeric atom types
   *
   * Look up a dihedral type from four atom type indices by
   * constructing a hyphen-delimited label (e.g., "C-C-N-H").
   *
   * \param  atype1  First atom type index
   * \param  atype2  Second atom type index
   * \param  atype3  Third atom type index
   * \param  atype4  Fourth atom type index
   * \return         Dihedral type index if types match in the specified order,
   *                 negative dihedral type index if types match in reverse order,
   *                 0 if there is no match found */
  int infer_dihedraltype(int atype1, int atype2, int atype3, int atype4);

  /*! Infer dihedral type from atom type labels
   *
   * \overload
   *
   * Look up a dihedral type from four atom type labels.
   *
   * \param  labels  Vector of four atom type label strings
   * \return         Dihedral type index if types match in the specified order,
   *                 negative dihedral type index if types match in reverse order,
   *                 0 if there is no match found */
  int infer_dihedraltype(const std::vector<std::string> &labels);

  /*! Infer improper type from four numeric atom types
   *
   * Look up an improper type from four atom type indices by
   * constructing a hyphen-delimited label (e.g., "C-N-C-C").
   *
   * \param  atype1  First atom type index (center atom)
   * \param  atype2  Second atom type index
   * \param  atype3  Third atom type index
   * \param  atype4  Fourth atom type index
   * \param  iorder  Order in which types were matched to improper type label
   * \return         Improper type index if types match in the specified order,
   *                 negative improper type index if types match but different order,
   *                 0 if there is no match found */
  int infer_impropertype(int atype1, int atype2, int atype3, int atype4,
                         std::array<int, 4> *iorder = nullptr);

  /*! Infer improper type from atom type labels
   *
   * \overload
   *
   * Look up an improper type from four atom type labels.
   *
   * \param  labels  Vector of four atom type label strings
   * \param  iorder  Order in which types were matched to improper type label
   * \return         Improper type index if types match in the specified order,
   *                 negative improper type index if types match but different order,
   *                 0 if there is no match found */
  int infer_impropertype(const std::vector<std::string> &labels,
                         std::array<int, 4> *iorder = nullptr);

  /*! @} */

  /*! Parse hyphen-delimited type label into components
   *
   * Split a hyphen-delimited label (e.g., "C-N-H") into individual type strings.
   * Validates that the number of components matches the expected count.
   *
   * \param       ntypes  Expected number of components
   * \param       label   Hyphen-delimited label string
   * \param[out]  types   Output vector to store component strings
   * \return      0 on success, -1 if component count doesn't match ntypes */
  int parse_typelabel(int ntypes, const std::string &label, std::vector<std::string> &types);

  /*! \name I/O methods for label map persistence
   * @{ */

  /*! Write label map to data file
   *
   * Output all type labels as sections to a LAMMPS data file.
   *
   * \param  fp  File pointer for writing */
  void write_data(FILE *fp);

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

  void check_labels();    //!< Check if type labels are self-consistent

  /*! @} */

 protected:
  int natomtypes, nbondtypes, nangletypes, ndihedraltypes, nimpropertypes;    //!< Type counts
  std::vector<std::string> typelabel, btypelabel,
      atypelabel;                                     //!< Label storage (atoms, bonds, angles)
  std::vector<std::string> dtypelabel, itypelabel;    //!< Label storage (dihedrals, impropers)

  std::unordered_map<std::string, int> typelabel_map;     //!< Atom label -> type mapping
  std::unordered_map<std::string, int> btypelabel_map;    //!< Bond label -> type mapping
  std::unordered_map<std::string, int> atypelabel_map;    //!< Angle label -> type mapping
  std::unordered_map<std::string, int> dtypelabel_map;    //!< Dihedral label -> type mapping
  std::unordered_map<std::string, int> itypelabel_map;    //!< Improper label -> type mapping

  int check_which_labels[4];    //!< Indicate check for bonds, angles, dihedrals, and/or impropers

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

  /*! Look up type label with given name or create new label if it doesn't exist

   * \param          mylabel     string with type label
   * \param[in,out]  labels      list of type labels
   * \param[in,out]  labels_map  label to numeric type hash table
   * \return numeric type
   */
  int find_or_create(const std::string &mylabel, std::vector<std::string> &labels,
                     std::unordered_map<std::string, int> &labels_map);

  /*! Look up numeric type of type label string in type map
   *
   * \param mylabel      type label
   * \param labels_map   label to type map
   * \return numeric type or -1 if not found */
  int search(const std::string &mylabel,
             const std::unordered_map<std::string, int> &labels_map) const;

  /*! Read a C-style string from a binary file and broadcast to world communicator
   *
   * the string buffer is allocated with new and must be freed by the calling code with delete[]
   *
   * \param fp  FILE pointer of the openend file
   * \return pointer to the allocated string buffer */
  char *read_string(FILE *fp);

  /*! Encode string to binary file.
   *
   * Must be only called from MPI rank 0.
   *
   * \param str string to write to the file
   * \param fp  FILE pointer of the opened file */
  void write_string(const std::string &str, FILE *fp);

  /*! Read integer from binary file and broadcast it to world communicator
   *
   * \param fp FILE pointer of the opened file
   * \return  integer value read from file */
  int read_int(FILE *fp);

  /*! Write out current label maps to a file for debugging
   *
   * \param filename  file name */
  void write_map(const std::string &filename);
};

}    // namespace LAMMPS_NS

#endif
