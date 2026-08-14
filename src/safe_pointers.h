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

/** Class to automatically close a FILE pointer when it goes out of scope

\verbatim embed:rst

This is a drop-in replacement for declaring a ``FILE *`` variable and can
be passed to functions or used in logical expressions the same way.  This
is particularly useful when a code block will throw an exception, or has
to close and re-open a file and similar.  It helps to simplify code and
reduces the risk of memory and file descriptor leaks.

Below are some usage examples:

.. code-block:: c++

  // Replace:
  FILE *fp = nullptr;
  // With:
  SafeFilePtr fp;

  // You can use "fp" as usual:
  fp = fopen("some.file","r");
  // a second assignment will automatically close the opened file
  fp = fopen("other.file", "r");
  // and assigning nullptr will just close it
  fp = nullptr;

  // There also is a custom constructor available as a shortcut
  SafeFilePtr fp(fopen("some.file", "r"));

  // You can indicate that a file was opened with popen() to call pclose() instead of fclose()
  SafeFilePtr fp;
  if (platform::has_compress_extension(filename)) {
    fp.set_pclose();
    fp = platform::compressed_write(filename);
  } else {
    fp = fopen(filename, "w");
  }
  if (!fp) error->one(FLERR, "Failed to open file {}: {}", filename, utils::getsyserror());

  // reading or writing works without needing to change the source code
  fputs("write text to file\n", fp);
  char buffer[100];
  utils::sfgets(FLERR, buffer, 100, fp, filename, error);

\endverbatim
*/
class SafeFilePtr {
 public:
  SafeFilePtr() : fp(nullptr), use_pclose(false) {};
  SafeFilePtr(FILE *_fp, bool _use_pclose = false) : fp(_fp), use_pclose(_use_pclose) {};

  SafeFilePtr(const SafeFilePtr &) = delete;
  SafeFilePtr(SafeFilePtr &&o) noexcept : fp(o.fp), use_pclose(o.use_pclose) { o.fp = nullptr; }
  SafeFilePtr &operator=(const SafeFilePtr &) = delete;

  ~SafeFilePtr();

  /** Assign new file pointer and close old one if still open.
   *
   * The value of use_pclose determines whether `pclose()` is called or `fclose()`.
   * Assigning `nullptr` closes the file and resets use_pclose
   *
   * \param _fp  new file pointer, may be `nullptr`
   * \return reference to updated class instance */
  SafeFilePtr &operator=(FILE *_fp);

  /** Flag that the file pointer needs to be closed with `pclose()` instead of `fclose()` */
  void set_pclose() { use_pclose = true; }

  /** Custom type cast operator so that SafeFilePtr can be used where FILE * was used
   *
   * \return currently stored/monitored file pointer */
  operator FILE *() const { return fp; }

 private:
  FILE *fp;
  bool use_pclose;
};
}    // namespace LAMMPS_NS

#endif
