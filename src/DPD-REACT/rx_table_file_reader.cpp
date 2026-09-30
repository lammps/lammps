#include "rx_table_file_reader.h"
#include "comm.h"
#include "utils.h"

using namespace LAMMPS_NS;

RxTableFileReader::RxTableFileReader(LAMMPS *lmp, const std::string &keyword,
                                     const std::string &filename_, const std::string &type,
                                     bool check_for_unread_tokens_) :
    Pointers(lmp),
    filename(filename_), check_for_unread_tokens(check_for_unread_tokens_),
    reader(lmp, filename_, type), is_table_already_read(false)
{

  if (comm->me != 0) { error->one(FLERR, "RxTableFileReader should only be called by proc 0!"); }

  auto *keyword_line = reader.find_section_start(keyword);

  if (!keyword_line) { error->one(FLERR, "Did not find keyword '{}' in '{}'", keyword, filename); }

  auto *line = reader.next_line();
  ValueTokenizer param_values(line);
  bool is_N_found = false;

  while (param_values.has_next()) {
    try {
      auto param_token = param_values.next_string();

      // The !is_N_found is necessary because "N" may appear more than
      // once on the parameter line, where it may be used, for
      // example, as a species name rather than an indicator of the
      // number of table entries.
      if (param_token == "N" && !is_N_found) {
        N = param_values.next_int();
        is_N_found = true;

      } else {
        param_tokens.emplace(std::move(param_token));
      }
    } catch (TokenizerException &e) {
      error->one(FLERR,
                 "Reading of table parameter line failed {}. File: {} "
                 "Full line: {}",
                 e.what(), filename, line);
    }
  }

  if (!is_N_found) {
    error->one(FLERR,
               "Failed to read value for N "
               "from line in table parameter file {}. "
               "Full line: {}",
               filename, line);
  }
}

bool RxTableFileReader::has_next_param_token() const
{
  return !param_tokens.empty();
}

std::string RxTableFileReader::next_param_token_as_string()
{
  std::string out_val(std::move(param_tokens.front()));
  param_tokens.pop();
  return out_val;
}

int RxTableFileReader::next_param_token_as_int()
{
  return utils::inumeric(FLERR, next_param_token_as_string(), true, lmp);
}

double RxTableFileReader::next_param_token_as_double()
{
  return utils::numeric(FLERR, next_param_token_as_string(), true, lmp);
}

RxTableFileReader::TableIndex_t RxTableFileReader::get_num_table_entries() const
{
  return N;
}
