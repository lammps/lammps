/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <https://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "reaxff_api.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "error.h"

namespace ReaxFF {

/* safe malloc */
void *smalloc(LAMMPS_NS::Error *error_ptr, rc_bigint n, const std::string &name)
{
  void *ptr;

  if (n <= 0) {
    auto errmsg = fmt::format("Invalid size {} for array {}. Returning NULL.", n, name);
    if (error_ptr)
      error_ptr->one(FLERR, errmsg);
    else
      fputs(errmsg.c_str(), stderr);

    return nullptr;
  }

  ptr = malloc(n);
  if (ptr == nullptr) {
    auto errmsg = fmt::format("Failed to allocate {} bytes for array {}", n, name);
    if (error_ptr)
      error_ptr->one(FLERR, errmsg);
    else
      fputs(errmsg.c_str(), stderr);
  }

  return ptr;
}

/* safe calloc */
void *scalloc(LAMMPS_NS::Error *error_ptr, rc_bigint n, rc_bigint size, const std::string &name)
{
  void *ptr;

  if (n <= 0) {
    auto errmsg = fmt::format("Invalid size {} for array {}. Returning NULL.\n", n, name);
    if (error_ptr)
      error_ptr->one(FLERR, errmsg);
    else
      fputs(errmsg.c_str(), stderr);
    return nullptr;
  }

  if (size <= 0) {
    auto errmsg = fmt::format("Elements size for array {} is {}. Returning NULL", name, size);
    if (error_ptr)
      error_ptr->one(FLERR, errmsg);
    else
      fputs(errmsg.c_str(), stderr);
    return nullptr;
  }

  ptr = calloc(n, size);
  if (ptr == nullptr) {
    auto errmsg = fmt::format("Failed to allocate {} bytes for array {}", n * size, name);
    if (error_ptr)
      error_ptr->one(FLERR, errmsg);
    else
      fputs(errmsg.c_str(), stderr);
  }

  return ptr;
}

/* safe free */
void sfree(LAMMPS_NS::Error *error_ptr, void *ptr, const std::string &name)
{
  if (ptr == nullptr) {
    auto errmsg = std::string("Trying to free the already free()'d pointer: ") + name;
    if (error_ptr)
      error_ptr->one(FLERR, errmsg);
    else
      fputs(errmsg.c_str(), stderr);
    return;
  }

  free(ptr);
  ptr = nullptr;
}
}    // namespace ReaxFF
