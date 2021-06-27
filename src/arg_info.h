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

#ifndef LMP_ARG_INFO_H
#define LMP_ARG_INFO_H

/*! \file arg_info.h */

#include <string>

namespace LAMMPS_NS {
class ArgInfo {
 public:
  // clang-format off
  /*! constants for argument types */
  enum ArgTypes {
    ERROR          =-2,
    UNKNOWN        =-1,
    NONE           = 0,
    X              = 1<<0,
    V              = 1<<1,
    F              = 1<<2,
    COMPUTE        = 1<<3,
    FIX            = 1<<4,
    VARIABLE       = 1<<5,
    KEYWORD        = 1<<6,
    TYPE           = 1<<7,
    MOLECULE       = 1<<8,
    DNAME          = 1<<9,
    INAME          = 1<<10,
    DENSITY_NUMBER = 1<<11,
    DENSITY_MASS   = 1<<12,
    MASS           = 1<<13,
    TEMPERATURE    = 1<<14,
    BIN1D          = 1<<15,
    BIN2D          = 1<<16,
    BIN3D          = 1<<17,
    BINSPHERE      = 1<<18,
    BINCYLINDER    = 1<<19
  };
  // clang-format on
  ArgInfo(const std::string &arg, int allowed = COMPUTE | FIX | VARIABLE);
  virtual ~ArgInfo() {}

 public:
  /*! get type of reference
     *
     * Return a type constant for the reference. This may be either
     * COMPUTE, FIX, VARIABLE (if not restricted to a subset of those
     * by the "allowed" argument of the constructor) or NONE, if it
     * if not a recognized or allowed reference, or UNKNOWN, in case
     * some error happened identifying or parsing the values of the indices
     *
     * \return  integer with a constant from ArgTypes enumerator */

  int get_type() const { return type; }

  /*! get dimension of reference
     *
     * This will return either 0, 1, 2 depending on whether the
     * reference has no, one or two "[{number}]" postfixes.
     *
     * \return  integer with the dimensionality of the reference */
  int get_dim() const { return dim; }

  /*! get index of first dimension
     *
     * This will return the number in the first "[{number}]"
     * postfix or 0 if there is no postfix.
     *
     * \return  integer with index or the postfix or 0 */
  int get_index1() const { return index1; }

  /*! get index of second dimension
     *
     * This will return the number in the second "[{number}]"
     * postfix or -1 if there is no second postfix.
     *
     * \return  integer with index of the postfix or -1 */
  int get_index2() const { return index2; }

  /*! return reference to the ID or name of the reference
     *
     * This string is pointing to an internal storage element and
     * is only valid to use while the ArgInfo class instance is
     * in scope. If you need a long-lived string make a copy
     * with copy_name().
     *
     * \return  C-style char * string */
  const char *get_name() const { return name.c_str(); }

  char *copy_name();

 private:
  std::string name;
  int type, dim, index1, index2;

  // disabled standard methods
  ArgInfo() {}
  ArgInfo(const ArgInfo &) {}
  void operator=(const ArgInfo &) {}
};
}    // namespace LAMMPS_NS
#endif
