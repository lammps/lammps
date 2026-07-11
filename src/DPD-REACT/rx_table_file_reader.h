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
#ifndef LMP_RX_TABLE_FILE_READER_H
#define LMP_RX_TABLE_FILE_READER_H

#include "error.h"
#include "pointers.h"
#include "table_file_reader.h"
#include "tokenizer.h"

#include <queue>

#if __cplusplus >= 202002L
#include <concepts>
#endif

namespace LAMMPS_NS {

class RxTableFileReader : protected Pointers {
 public:
  using TableIndex_t = int;

  RxTableFileReader(LAMMPS *lmp, const std::string &keyword, const std::string &filename,
                    const std::string &type, bool check_for_unread_tokens = true);

  [[nodiscard]] bool has_next_param_token() const;

  std::string next_param_token_as_string();
  int next_param_token_as_int();
  double next_param_token_as_double();

  [[nodiscard]] TableIndex_t get_num_table_entries() const;

  template <typename TableLineInvocable_t>
#if __cplusplus >= 202002L
    requires std::invocable<TableLineInvocable_t, TableIndex_t, ValueTokenizer &>
#endif
  void read_in_table_data(TableLineInvocable_t table_line_invocable)
  {

    if (is_table_already_read)
      error->one(FLERR, Error::NOLASTLINE, "Attempted to read table more than once.");

    for (TableIndex_t i = 0; i < N; i++) {
      auto *line = reader.next_line();

      try {
        ValueTokenizer table_values(line);
        table_line_invocable(i, table_values);

        if (check_for_unread_tokens && table_values.has_next()) {
          error->one(FLERR,
                     "Misplaced characters at end of line in "
                     "file {}. Full line: {}",
                     filename, line);
        }

      } catch (TokenizerException &e) {
        error->one(FLERR, "{}. File name: {} Full line: {}", e.what(), filename, line);
      }
    }

    is_table_already_read = true;
  }

 private:
  std::string filename;
  bool check_for_unread_tokens;
  TableFileReader reader;
  bool is_table_already_read;
  TableIndex_t N;
  std::queue<std::string> param_tokens;
};

}    // namespace LAMMPS_NS

#endif
