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

#ifndef LMP_SAFE_POINTERS_H
#define LMP_SAFE_POINTERS_H

// collection of smart pointers for specific purposes

#include <cstdio>

namespace LAMMPS_NS {

/** Class to automatically close a FILE pointer when a class instance goes out of scope.

\verbatim embed:rst

Drop in replacement for ``FILE *``. Use as ``SafeFilePtr fp;`` instead of
``FILE *fp = nullptr;`` and there is no more need to explicitly call
``fclose(fp)``.

\endverbatim
*/
class SafeFilePtr {
 public:
  SafeFilePtr() : fp(nullptr) {};
  SafeFilePtr(FILE *_fp) : fp(_fp) {};

  SafeFilePtr(const SafeFilePtr &) = delete;
  SafeFilePtr(SafeFilePtr &&o) noexcept : fp(o.fp) { o.fp = nullptr; }
  SafeFilePtr &operator=(const SafeFilePtr &) = delete;

  ~SafeFilePtr()
  {
    if (fp) fclose(fp);
  }

  SafeFilePtr &operator=(FILE *_fp)
  {
    if (fp && (fp != _fp)) fclose(fp);
    fp = _fp;
    return *this;
  }
  operator FILE *() const { return fp; }

 private:
  FILE *fp;
};
}    // namespace LAMMPS_NS

#endif
