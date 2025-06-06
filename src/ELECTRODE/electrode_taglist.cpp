/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meissner (TUHH)
------------------------------------------------------------------------- */

#include "electrode_taglist.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "text_file_reader.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <unordered_map>

using namespace LAMMPS_NS;

ElectrodeTaglist::ElectrodeTaglist(LAMMPS *lmp, std::vector<int> group_bits) : Pointers(lmp)
{
  int const nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  int const nprocs = comm->nprocs;
  int *recvcounts = new int[nprocs];
  int *displs = new int[nprocs];

  // assign a tag to each matrix index sorted by group and by tag
  taglist_bygroup = std::vector<tagint>();
  for (int gbit : group_bits) {
    std::vector<tagint> taglist_local_group;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & gbit) { taglist_local_group.push_back(tag[i]); }
    }
    // gather from all cpus for this group
    int gnum_local = taglist_local_group.size();
    MPI_Allgather(&gnum_local, 1, MPI_INT, recvcounts, 1, MPI_INT, world);
    displs[0] = 0;
    for (int i = 1; i < nprocs; i++) { displs[i] = displs[i - 1] + recvcounts[i - 1]; }
    int const gnum = displs[nprocs - 1] + recvcounts[nprocs - 1];
    std::vector<tagint> taglist_all(gnum);
    MPI_Allgatherv(taglist_local_group.data(), gnum_local, MPI_LMP_TAGINT, taglist_all.data(),
                   recvcounts, displs, MPI_LMP_TAGINT, world);
    std::sort(taglist_all.begin(), taglist_all.end());
    for (tagint t : taglist_all) taglist_bygroup.push_back(t);
  }
  n = taglist_bygroup.size();

  // taglist only sorted by tag not group, same order as in computes
  taglist = taglist_bygroup;
  std::sort(taglist.begin(), taglist.end());

  tag_to_iele = std::unordered_map<tagint, int>();
  tag_to_iele.reserve(n);
  for (std::size_t i = 0; i < n; i++) { tag_to_iele.insert(std::pair<tagint, int>(taglist[i], i)); }

  // group_idx allows mapping a vector that is sorted by taglist to being
  // ordered by taglist_bygroup
  group_idx = std::vector<int>(n);
  for (std::size_t i{0}; i < n; i++) { group_idx[i] = tag_to_iele[taglist_bygroup[i]]; }

  delete[] recvcounts;
  delete[] displs;
}

/* ---------------------------------------------------------------------- */

ElectrodeTaglist::~ElectrodeTaglist() {}

/* ---------------------------------------------------------------------- */

std::unordered_map<tagint, int> ElectrodeTaglist::get_tag_to_iele()
{
  return tag_to_iele;
}

/* ---------------------------------------------------------------------- */

void ElectrodeTaglist::write_to_file(const std::string file_str, double *unsorted)
{
  if (comm->me == 0) write_to_file(file_str, sort_by_group(unsorted));
}
/* ---------------------------------------------------------------------- */

void ElectrodeTaglist::write_to_file(const std::string file_str, double **unsorted)
{
  if (comm->me == 0) write_to_file(file_str, sort_by_group(unsorted));
}

/* ---------------------------------------------------------------------- */

void ElectrodeTaglist::write_to_file(const std::string file_str,
                                     std::vector<std::vector<double>> mat)
{
  assert(comm->me == 0);
  FILE *file = fopen(file_str.c_str(), "w");
  if (file == nullptr) error->one(FLERR, "Cannot open file {}: {}", file_str, utils::getsyserror());
  for (const auto &t : taglist_bygroup) fmt::print(file, "{:20}", t);
  fputs("\n", file);
  for (const auto &vec : mat) {
    for (const auto &x : vec) fmt::print(file, "{:20.11e}", x);
    fputs("\n", file);
  }
  fclose(file);
}

/* ---------------------------------------------------------------------- */

std::vector<std::vector<double>> ElectrodeTaglist::sort_by_group(double **mat)
{
  std::vector<std::vector<double>> ordered_mat(n, std::vector<double>(n));
  for (std::size_t i = 0; i < n; i++) {
    bigint const gi = group_idx[i];
    assert(gi < n);
    for (std::size_t j = 0; j < n; j++) ordered_mat[gi][group_idx[j]] = mat[i][j];
  }
  return ordered_mat;
}

/* ---------------------------------------------------------------------- */

std::vector<std::vector<double>> ElectrodeTaglist::sort_by_group(double *vec)
{
  std::vector<std::vector<double>> ordered_vec(n, std::vector<double>(1));
  for (std::size_t i = 0; i < n; i++) ordered_vec[group_idx[i]][0] = vec[i];
  return ordered_vec;
}

/* ---------------------------------------------------------------------- */

void ElectrodeTaglist::read_from_file(const std::string &input_file, double **array,
                                      const std::string &filetype)
{
  if (comm->me == 0) {
    std::vector<std::vector<double>> matrix;
    std::vector<tagint> tags;
    try {
      TextFileReader reader(input_file, filetype);
      int bufsize = n * 20 + 4;
      reader.set_bufsize(bufsize > 100 ? bufsize : 100);

      // get line with tags
      auto values = reader.next_values(n);
      for (std::size_t i = 0; i < n; ++i) tags.push_back(values.next_tagint());

      std::vector<double> a_line;
      for (std::size_t i = 0; i < n; ++i) {
        a_line.clear();
        values = reader.next_values(n);
        for (std::size_t j = 0; j < n; ++j) a_line.push_back(values.next_double());
        matrix.push_back(a_line);
      }
    } catch (std::exception &e) {
      error->one(FLERR, "Error parsing {} file: {}", filetype, e.what());
    }

    std::vector<tagint> idx;
    for (const auto &t : taglist) {
      for (std::size_t i = 0; i < tags.size(); i++) {
        if (t == tags[i]) {
          idx.push_back(i);
          break;
        }
      }
    }
    if (idx.size() != n) error->all(FLERR, "Read tags do not match taglist of fix electrode");
    for (std::size_t i = 0; i < n; i++) {
      bigint const ii = idx[i];
      for (std::size_t j = 0; j < n; j++) array[i][j] = matrix[ii][idx[j]];
    }
  }
  MPI_Bcast(&array[0][0], n * n, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

double ElectrodeTaglist::memory_usage()
{
  double bytes = 0.;
  bytes += taglist.capacity() * sizeof(tagint);
  bytes += taglist_bygroup.capacity() * sizeof(tagint);
  bytes += group_idx.capacity() * sizeof(int);
  bytes += tag_to_iele.bucket_count() * (sizeof(tagint) + sizeof(int));    // TODO check
  return bytes;
}
