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

#include "utils.h"

#include "arg_info.h"
#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "fix.h"
#include "fmt/chrono.h"
#include "input.h"
#include "label_map.h"
#include "memory.h"
#include "modify.h"
#include "text_file_reader.h"
#include "universe.h"
#include "update.h"
#include "variable.h"

#include <cctype>
#include <cerrno>
#include <cstring>
#include <ctime>

/*! \file utils.cpp */

/*
 * Mini regex-module adapted from https://github.com/kokke/tiny-regex-c
 * which is in the public domain.
 *
 * Supports:
 * ---------
 *   '.'        Dot, matches any character
 *   '^'        Start anchor, matches beginning of string
 *   '$'        End anchor, matches end of string
 *   '*'        Asterisk, match zero or more (greedy)
 *   '+'        Plus, match one or more (greedy)
 *   '?'        Question, match zero or one (non-greedy)
 *   '[abc]'    Character class, match if one of {'a', 'b', 'c'}
 *   '[a-zA-Z]' Character ranges, the character set of the ranges { a-z | A-Z }
 *   '\s'       Whitespace, \t \f \r \n \v and spaces
 *   '\S'       Non-whitespace
 *   '\w'       Alphanumeric, [a-zA-Z0-9_]
 *   '\W'       Non-alphanumeric
 *   '\d'       Digits, [0-9]
 *   '\D'       Non-digits
 *   '\i'       Integer chars, [0-9], '+' and '-'
 *   '\I'       Non-integers
 *   '\f'       Floating point number chars, [0-9], '.', 'e', 'E', '+' and '-'
 *   '\F'       Non-floats
 *
 * *NOT* supported:
 *   '[^abc]'   Inverted class
 *   'a|b'      Branches
 *   '(abc)+'   Groups
 */

extern "C" {
/** Match text against a (simplified) regular expression
   * (regexp will be compiled automatically). */
static int re_match(const char *text, const char *pattern);

/** Match find substring that matches a (simplified) regular expression
   * (regexp will be compiled automatically). */
static int re_find(const char *text, const char *pattern, int *matchlen);
}

////////////////////////////////////////////////////////////////////////
// Merge sort support functions

static void do_merge(int *idx, int *buf, int llo, int lhi, int rlo, int rhi, void *ptr,
                     int (*comp)(int, int, void *));
static void insertion_sort(int *index, int num, void *ptr, int (*comp)(int, int, void *));

////////////////////////////////////////////////////////////////////////

using namespace LAMMPS_NS;

/** More flexible and specific matching of a string against a pattern.
 *  This function is supposed to be a more safe, more specific and
 *  simple to use API to find pattern matches. The purpose is to replace
 *  uses of either strncmp() or strstr() in the code base to find
 *  sub-strings safely. With strncmp() finding prefixes, the number of
 *  characters to match must be counted, which can lead to errors,
 *  while using "^pattern" will do the same with less problems.
 *  Matching for suffixes using strstr() is not as specific as 'pattern$',
 *  and complex matches, e.g. "^rigid.*\/small.*", to match all small
 *  body optimized rigid fixes require only one test.
 *
 *  The use of std::string arguments allows for simple concatenation
 *  even with char * type variables.
 *  Example: utils::strmatch(text, std::string("^") + charptr)
 */
bool utils::strmatch(const std::string &text, const std::string &pattern)
{
  const int pos = re_match(text.c_str(), pattern.c_str());
  return (pos >= 0);
}

/** This function is a companion function to utils::strmatch(). Arguments
 *  and logic is the same, but instead of a boolean, it returns the
 *  sub-string that matches the regex pattern.  There can be only one match.
 *  This can be used as a more flexible alternative to strstr().
 */
std::string utils::strfind(const std::string &text, const std::string &pattern)
{
  int matchlen;
  const int pos = re_find(text.c_str(), pattern.c_str(), &matchlen);
  if ((pos >= 0) && (matchlen > 0))
    return text.substr(pos, matchlen);
  else
    return "";
}

void utils::missing_cmd_args(const std::string &file, int line, const std::string &cmd,
                             Error *error)
{
  if (error) error->all(file, line, "Illegal {} command: missing argument(s)", cmd);
}

/* specialization for the case of just a single string argument */

void utils::logmesg(LAMMPS *lmp, const std::string &mesg)
{
  if (lmp->screen) fputs(mesg.c_str(), lmp->screen);
  if (lmp->logfile) fputs(mesg.c_str(), lmp->logfile);
}

void utils::fmtargs_logmesg(LAMMPS *lmp, fmt::string_view format, fmt::format_args args)
{
  try {
    logmesg(lmp, fmt::vformat(format, args));
  } catch (fmt::format_error &e) {
    logmesg(lmp, std::string(e.what()) + "\n");
  }
}

std::string utils::errorurl(int errorcode)
{
  return fmt::format("\nFor more information see https://docs.lammps.org/err{:04d}", errorcode);
}

void utils::flush_buffers(LAMMPS *lmp)
{
  if (lmp->screen) fflush(lmp->screen);
  if (lmp->logfile) fflush(lmp->logfile);
  if (lmp->universe->uscreen) fflush(lmp->universe->uscreen);
  if (lmp->universe->ulogfile) fflush(lmp->universe->ulogfile);
}

/* define this here, so we won't have to include the headers
   everywhere and utils.h will more likely be included anyway. */

std::string utils::getsyserror()
{
  return {strerror(errno)};
}

// read line into buffer. if line is too long keep reading until EOL or EOF
// but return only the first part with a newline at the end.

char *utils::fgets_trunc(char *buf, int size, FILE *fp)
{
  constexpr int MAXDUMMY = 256;
  char dummy[MAXDUMMY];
  char *ptr = fgets(buf, size, fp);

  // EOF?
  if (!ptr) return nullptr;

  int n = strlen(buf);

  // check the string being read in:
  // - if string is shorter than the buffer make sure it has a final newline and return
  // - if string is exactly the size of the buffer and has a final newline return
  // - otherwise truncate with final newline and read into dummy buffer until EOF or newline is found
  if (n < size - 1) {
    if (buf[n - 1] != '\n') {
      buf[n] = '\n';
      buf[n + 1] = '\0';
    }
    return buf;
  } else if (buf[n - 1] == '\n') {
    return buf;
  } else
    buf[size - 2] = '\n';

  // continue reading into dummy buffer until end of line or file
  do {
    ptr = fgets(dummy, MAXDUMMY, fp);
    if (ptr)
      n = strlen(ptr);
    else
      n = 0;
  } while (n == MAXDUMMY - 1 && ptr[MAXDUMMY - 1] != '\n');

  // return truncated chunk
  return buf;
}

/* like fgets() but aborts with an error or EOF is encountered */
void utils::sfgets(const char *srcname, int srcline, char *s, int size, FILE *fp,
                   const char *filename, Error *error)
{
  constexpr int MAXPATHLENBUF = 1024;
  char *rv = fgets(s, size, fp);
  if (rv == nullptr) {    // something went wrong
    char buf[MAXPATHLENBUF];
    std::string errmsg;

    // try to figure out the file name from the file pointer
    if (!filename) filename = platform::guesspath(fp, buf, MAXPATHLENBUF);

    if (feof(fp)) {
      errmsg = "Unexpected end of file while reading file '";
    } else if (ferror(fp)) {
      errmsg = "Unexpected error while reading file '";
    } else {
      errmsg = "Unexpected short read while reading file '";
    }
    errmsg += filename;
    errmsg += "'";

    if (error) error->one(srcname, srcline, errmsg);
    if (s) *s = '\0';    // truncate string to empty in case error is null pointer
  }
}

/* like fread() but aborts with an error or EOF is encountered */
void utils::sfread(const char *srcname, int srcline, void *s, size_t size, size_t num, FILE *fp,
                   const char *filename, Error *error)
{
  constexpr int MAXPATHLENBUF = 1024;
  size_t rv = fread(s, size, num, fp);
  if (rv != num) {    // something went wrong
    char buf[MAXPATHLENBUF];
    std::string errmsg;

    // try to figure out the file name from the file pointer
    if (!filename) filename = platform::guesspath(fp, buf, MAXPATHLENBUF);

    if (feof(fp)) {
      errmsg = "Unexpected end of file while reading file '";
    } else if (ferror(fp)) {
      errmsg = "Unexpected error while reading file '";
    } else {
      errmsg = "Unexpected short read while reading file '";
    }
    errmsg += filename;
    errmsg += "'";

    if (error) error->one(srcname, srcline, errmsg);
  }
}

/* ------------------------------------------------------------------ */

/* read N lines and broadcast */
int utils::read_lines_from_file(FILE *fp, int nlines, int nmax, char *buffer, int me, MPI_Comm comm)
{
  char *ptr = buffer;
  *ptr = '\0';
  int mylines = 0;

  if (me == 0) {
    if (fp) {
      for (int i = 0; i < nlines; i++) {
        ptr = fgets_trunc(ptr, nmax, fp);
        if (!ptr) break;    // EOF?
        // advance ptr to end of string
        ptr += strlen(ptr);
        // ensure buffer is null terminated. null char is start of next line.
        *ptr = '\0';
        // count line
        ++mylines;
      }
    }
  }

  int n = strlen(buffer);
  if (nlines != mylines) n = 0;
  MPI_Bcast(&n, 1, MPI_INT, 0, comm);
  if (n == 0) return 1;
  MPI_Bcast(buffer, n + 1, MPI_CHAR, 0, comm);
  return 0;
}

/* ------------------------------------------------------------------ */

std::string utils::check_packages_for_style(const std::string &style, const std::string &name,
                                            LAMMPS *lmp)
{
  std::string errmsg = "Unrecognized " + style + " style '" + name + "'";
  const char *pkg = lmp->match_style(style.c_str(), name.c_str());

  if (pkg) {
    errmsg += fmt::format(" is part of the {} package", pkg);
    if (LAMMPS::is_installed_pkg(pkg))
      errmsg += ", but seems to be missing because of a dependency";
    else
      errmsg += " which is not enabled in this LAMMPS binary.";
  }
  return errmsg;
}

/* ----------------------------------------------------------------------
   read a boolean value from a string
   transform to lower case before checking
   generate an error if is not a legitimate boolean
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

int utils::logical(const char *file, int line, const std::string &str, bool do_abort, LAMMPS *lmp)
{
  if (str.empty()) {
    const char msg[] = "Expected boolean parameter instead of NULL or empty string "
                       "in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }

  // convert to ascii
  std::string buf(str);
  if (has_utf8(buf)) buf = utf8_subst(buf);

  int rv = 0;
  if ((buf == "yes") || (buf == "on") || (buf == "true") || (buf == "1")) {
    rv = 1;
  } else if ((buf == "no") || (buf == "off") || (buf == "false") || (buf == "0")) {
    rv = 0;
  } else {
    std::string msg("Expected boolean parameter instead of '");
    msg += buf + "' in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }
  return rv;
}

/* ----------------------------------------------------------------------
   wrapper for logical() that accepts a char pointer instead of a string
------------------------------------------------------------------------- */

int utils::logical(const char *file, int line, const char *str, bool do_abort, LAMMPS *lmp)
{
  if (str)
    return logical(file, line, std::string(str), do_abort, lmp);
  else
    return logical(file, line, std::string(""), do_abort, lmp);
}

/* ----------------------------------------------------------------------
   read a floating point value from a string
   generate an error if not a legitimate floating point value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

double utils::numeric(const char *file, int line, const std::string &str, bool do_abort,
                      LAMMPS *lmp)
{
  if (str.empty()) {
    const char msg[] = "Expected floating point parameter instead of"
                       " NULL or empty string in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }

  std::string buf(str);
  if (has_utf8(buf)) buf = utf8_subst(buf);

  if (!is_double(buf)) {
    std::string msg("Expected floating point parameter instead of '");
    msg += buf + "' in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }

  double rv = 0;
  try {
    rv = stod(buf);
  } catch (std::invalid_argument const &) {
    auto msg = fmt::format("Floating point number {} in input script or data file is invalid", buf);
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  } catch (std::out_of_range const &) {
    auto msg =
        fmt::format("Floating point number {} in input script or data file is out of range", buf);
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }
  return rv;
}

/* ----------------------------------------------------------------------
   wrapper for numeric() that accepts a char pointer instead of a string
------------------------------------------------------------------------- */

double utils::numeric(const char *file, int line, const char *str, bool do_abort, LAMMPS *lmp)
{
  if (str)
    return numeric(file, line, std::string(str), do_abort, lmp);
  else
    return numeric(file, line, std::string(""), do_abort, lmp);
}

/* ----------------------------------------------------------------------
   read an integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

int utils::inumeric(const char *file, int line, const std::string &str, bool do_abort, LAMMPS *lmp)
{
  if (str.empty()) {
    const char msg[] = "Expected integer parameter instead of"
                       " NULL or empty string in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }

  std::string buf(str);
  if (has_utf8(buf)) buf = utf8_subst(buf);

  if (!is_integer(buf)) {
    std::string msg("Expected integer parameter instead of '");
    msg += buf + "' in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }

  int rv = 0;
  try {
    rv = stoi(buf);
  } catch (std::out_of_range const &) {
    auto msg = fmt::format("Integer {} in input script or data file is out of range", buf);
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }
  return rv;
}

/* ----------------------------------------------------------------------
   wrapper for inumeric() that accepts a char pointer instead of a string
------------------------------------------------------------------------- */

int utils::inumeric(const char *file, int line, const char *str, bool do_abort, LAMMPS *lmp)
{
  if (str)
    return inumeric(file, line, std::string(str), do_abort, lmp);
  else
    return inumeric(file, line, std::string(""), do_abort, lmp);
}

/* ----------------------------------------------------------------------
   read a big integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

bigint utils::bnumeric(const char *file, int line, const std::string &str, bool do_abort,
                       LAMMPS *lmp)
{
  if (str.empty()) {
    const char msg[] = "Expected integer parameter instead of"
                       " NULL or empty string in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }

  std::string buf(str);
  if (has_utf8(buf)) buf = utf8_subst(buf);

  if (!is_integer(buf)) {
    std::string msg("Expected integer parameter instead of '");
    msg += buf + "' in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }

  long long rv = 0;
  try {
    rv = stoll(buf);
    if (rv > MAXBIGINT) throw std::out_of_range("64-bit");
  } catch (std::out_of_range const &) {
    auto msg = fmt::format("Integer {} in input script or data file is out of range", buf);
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }
  return static_cast<bigint>(rv);
}

/* ----------------------------------------------------------------------
   wrapper for bnumeric() that accepts a char pointer instead of a string
------------------------------------------------------------------------- */

bigint utils::bnumeric(const char *file, int line, const char *str, bool do_abort, LAMMPS *lmp)
{
  if (str)
    return bnumeric(file, line, std::string(str), do_abort, lmp);
  else
    return bnumeric(file, line, std::string(""), do_abort, lmp);
}

/* ----------------------------------------------------------------------
   read a tag integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

tagint utils::tnumeric(const char *file, int line, const std::string &str, bool do_abort,
                       LAMMPS *lmp)
{
  if (str.empty()) {
    const char msg[] = "Expected integer parameter instead of"
                       " NULL or empty string in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }

  std::string buf(str);
  if (has_utf8(buf)) buf = utf8_subst(buf);

  if (!is_integer(buf)) {
    std::string msg("Expected integer parameter instead of '");
    msg += buf + "' in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }

  long long rv = 0;
  try {
    rv = stoll(buf);
    if (rv > MAXTAGINT) throw std::out_of_range("64-bit");
  } catch (std::out_of_range const &) {
    auto msg = fmt::format("Integer {} in input script or data file is out of range", buf);
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }
  return static_cast<tagint>(rv);
}

/* ----------------------------------------------------------------------
   wrapper for tnumeric() that accepts a char pointer instead of a string
------------------------------------------------------------------------- */

tagint utils::tnumeric(const char *file, int line, const char *str, bool do_abort, LAMMPS *lmp)
{
  if (str)
    return tnumeric(file, line, std::string(str), do_abort, lmp);
  else
    return tnumeric(file, line, std::string(""), do_abort, lmp);
}

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterisk
------------------------------------------------------------------------- */
// clang-format off
template <typename TYPE>
void utils::bounds(const char *file, int line, const std::string &str,
                   bigint nmin, bigint nmax, TYPE &nlo, TYPE &nhi, Error *error)
{
  nlo = nhi = -1;

  // check for illegal charcters
  size_t found = str.find_first_not_of("*-0123456789");
  if (found != std::string::npos) {
    if (error) error->all(file, line, "Invalid range string: {}", str);
    return;
  }

  found = str.find_first_of('*');
  if (found == std::string::npos) {    // contains no '*'
    nlo = nhi = strtol(str.c_str(), nullptr, 10);
  } else if (str.size() == 1) {    // is only '*'
    nlo = nmin;
    nhi = nmax;
  } else if (found == 0) {    // is '*j'
    nlo = nmin;
    nhi = strtol(str.substr(1).c_str(), nullptr, 10);
  } else if (str.size() == found + 1) {    // is 'i*'
    nlo = strtol(str.c_str(), nullptr, 10);
    nhi = nmax;
  } else {    // is 'i*j'
    nlo = strtol(str.c_str(), nullptr, 10);
    nhi = strtol(str.substr(found + 1).c_str(), nullptr, 10);
  }

  if (error) {
    if ((nlo <= 0) || (nhi <= 0))
      error->all(file, line, "Invalid range string: {}", str);

    if (nlo < nmin)
      error->all(file, line, "Numeric index {} is out of bounds ({}-{})", nlo, nmin, nmax);
    else if (nhi > nmax)
      error->all(file, line, "Numeric index {} is out of bounds ({}-{})", nhi, nmin, nmax);
    else if (nlo > nhi)
      error->all(file, line, "Numeric index {} is out of bounds ({}-{})", nlo, nmin, nhi);
  }
}

template void utils::bounds<>(const char *, int, const std::string &,
                              bigint, bigint, int &, int &, Error *);
template void utils::bounds<>(const char *, int, const std::string &,
                              bigint, bigint, long &, long &, Error *);
template void utils::bounds<>(const char *, int, const std::string &,
                              bigint, bigint, long long &, long long &, Error *);
// clang-format on

/* -------------------------------------------------------------------------
   Expand list of arguments in arg to earg if arg contains wildcards
------------------------------------------------------------------------- */

int utils::expand_args(const char *file, int line, int narg, char **arg, int mode, char **&earg,
                       LAMMPS *lmp)
{
  int iarg;

  char *ptr = nullptr;
  for (iarg = 0; iarg < narg; iarg++) {
    ptr = strchr(arg[iarg], '*');
    if (ptr) break;
  }

  if (!ptr) {
    earg = arg;
    return narg;
  }

  // maxarg should always end up equal to newarg, so caller can free earg

  int maxarg = narg - iarg;
  earg = (char **) lmp->memory->smalloc(maxarg * sizeof(char *), "input:earg");

  int newarg = 0, expandflag, nlo, nhi, nmax;
  std::string id, wc, tail;

  for (iarg = 0; iarg < narg; iarg++) {
    std::string word(arg[iarg]);
    expandflag = 0;

    // match grids

    if (strmatch(word, "^[cf]_\\w+:\\w+:\\w+\\[\\d*\\*\\d*\\]")) {
      auto gridid = utils::parse_grid_id(FLERR, word, lmp->error);

      size_t first = gridid[2].find('[');
      size_t second = gridid[2].find(']', first + 1);
      id = gridid[2].substr(0, first);
      wc = gridid[2].substr(first + 1, second - first - 1);
      tail = gridid[2].substr(second + 1);

      // grids from compute

      if (gridid[0][0] == 'c') {

        auto compute = lmp->modify->get_compute_by_id(gridid[0].substr(2));
        if (compute && compute->pergrid_flag) {

          int dim = 0;
          int igrid = compute->get_grid_by_name(gridid[1], dim);

          if (igrid >= 0) {

            int ncol = 0;
            compute->get_griddata_by_name(igrid, id, ncol);
            nmax = ncol;

            expandflag = 1;
          }
        }
        // grids from fix

      } else if (gridid[0][0] == 'f') {

        auto fix = lmp->modify->get_fix_by_id(gridid[0].substr(2));
        if (fix && fix->pergrid_flag) {

          int dim = 0;
          int igrid = fix->get_grid_by_name(gridid[1], dim);

          if (igrid >= 0) {

            int ncol = 0;
            fix->get_griddata_by_name(igrid, id, ncol);
            nmax = ncol;

            expandflag = 1;
          }
        }
      }

      // expand wild card string to nlo/nhi numbers

      if (expandflag) {
        utils::bounds(file, line, wc, 1, nmax, nlo, nhi, lmp->error);

        if (newarg + nhi - nlo + 1 > maxarg) {
          maxarg += nhi - nlo + 1;
          earg = (char **) lmp->memory->srealloc(earg, maxarg * sizeof(char *), "input:earg");
        }

        for (int index = nlo; index <= nhi; index++) {
          earg[newarg] =
              utils::strdup(fmt::format("{}:{}:{}[{}]{}", gridid[0], gridid[1], id, index, tail));
          newarg++;
        }
      }

      // match compute, fix, or custom property array reference with a '*' wildcard
      // number range in the first pair of square brackets

    } else if (strmatch(word, "^[cfv]_\\w+\\[\\d*\\*\\d*\\]") ||
               strmatch(word, "^[id]2_\\w+\\[\\d*\\*\\d*\\]")) {

      // split off the compute/fix/property ID, the wildcard and trailing text

      size_t first = word.find('[');
      size_t second = word.find(']', first + 1);
      if (word[1] == '2')
        id = word.substr(3, first - 3);
      else
        id = word.substr(2, first - 2);

      wc = word.substr(first + 1, second - first - 1);
      tail = word.substr(second + 1);

      // compute

      if (word[0] == 'c') {
        auto compute = lmp->modify->get_compute_by_id(id);

        // check for global vector/array, peratom array, local array

        if (compute) {
          if (mode == 0 && compute->vector_flag) {
            nmax = compute->size_vector;
            expandflag = 1;
          } else if (mode == 1 && compute->array_flag) {
            nmax = compute->size_array_cols;
            expandflag = 1;
          } else if (compute->peratom_flag && compute->size_peratom_cols) {
            nmax = compute->size_peratom_cols;
            expandflag = 1;
          } else if (compute->local_flag && compute->size_local_cols) {
            nmax = compute->size_local_cols;
            expandflag = 1;
          }
        }

        // fix

      } else if (word[0] == 'f') {
        auto fix = lmp->modify->get_fix_by_id(id);

        // check for global vector/array, peratom array, local array

        if (fix) {
          if (mode == 0 && fix->vector_flag) {
            nmax = fix->size_vector;
            expandflag = 1;
          } else if (mode == 1 && fix->array_flag) {
            nmax = fix->size_array_cols;
            expandflag = 1;
          } else if (fix->peratom_flag && fix->size_peratom_cols) {
            nmax = fix->size_peratom_cols;
            expandflag = 1;
          } else if (fix->local_flag && fix->size_local_cols) {
            nmax = fix->size_local_cols;
            expandflag = 1;
          }
        }

        // vector variable

      } else if (word[0] == 'v') {
        int index = lmp->input->variable->find(id.c_str());

        // check for global vector/array, peratom array, local array

        if (index >= 0) {
          if (mode == 0 && lmp->input->variable->vectorstyle(index)) {
            utils::bounds(file, line, wc, 1, MAXSMALLINT, nlo, nhi, lmp->error);
            if (nhi < MAXSMALLINT) {
              nmax = nhi;
              expandflag = 1;
            }
          }
        }

        // only match custom array reference with a '*' wildcard
        // number range in the first pair of square brackets

      } else if ((word[0] == 'i') || (word[0] == 'd')) {
        int flag, cols;
        int icustom = lmp->atom->find_custom(id.c_str(), flag, cols);

        if ((icustom >= 0) && (mode == 1) && (cols > 0)) {

          // check for custom per-atom array

          if (((word[0] == 'i') && (flag == 0)) || ((word[0] == 'd') && (flag == 1))) {
            nmax = cols;
            expandflag = 1;
          }
        }
      }

      // expansion will take place

      if (expandflag) {

        // expand wild card string to nlo/nhi numbers

        utils::bounds(file, line, wc, 1, nmax, nlo, nhi, lmp->error);

        if (newarg + nhi - nlo + 1 > maxarg) {
          maxarg += nhi - nlo + 1;
          earg = (char **) lmp->memory->srealloc(earg, maxarg * sizeof(char *), "input:earg");
        }

        for (int index = nlo; index <= nhi; index++) {
          if (word[1] == '2')
            earg[newarg] = utils::strdup(fmt::format("{}2_{}[{}]{}", word[0], id, index, tail));
          else
            earg[newarg] = utils::strdup(fmt::format("{}_{}[{}]{}", word[0], id, index, tail));
          newarg++;
        }
      }
    }

    // no expansion: duplicate original string

    if (!expandflag) {
      if (newarg == maxarg) {
        maxarg++;
        earg = (char **) lmp->memory->srealloc(earg, maxarg * sizeof(char *), "input:earg");
      }
      earg[newarg] = utils::strdup(word);
      newarg++;
    }
  }

  // printf("NEWARG %d\n",newarg); for (int i = 0; i < newarg; i++) printf("  arg %d: %s\n",i,earg[i]);
  return newarg;
}

static const char *labeltypes[] = {"Atom", "Bond", "Angle", "Dihedral", "Improper"};

/* -------------------------------------------------------------------------
   Expand type string to numeric string from labelmap.
   Return copy of expanded type or null pointer.
------------------------------------------------------------------------- */

char *utils::expand_type(const char *file, int line, const std::string &str, int mode, LAMMPS *lmp)
{
  if (!lmp) return nullptr;
  if (!lmp->atom->labelmapflag) return nullptr;

  const std::string typestr = utils::utf8_subst(utils::trim(str));
  if (is_type(typestr) == 1) {
    if (!lmp->atom->labelmapflag)
      lmp->error->all(file, line, "{} type string {} cannot be used without a labelmap",
                      labeltypes[mode], typestr);

    int type = lmp->atom->lmap->find(typestr, mode);
    if (type == -1)
      lmp->error->all(file, line, "{} type string {} not found in labelmap", labeltypes[mode],
                      typestr);

    return utils::strdup(std::to_string(type));
  } else
    return nullptr;
}

/* ----------------------------------------------------------------------
   Check grid reference for valid Compute or Fix which produces per-grid data
   errstr = name of calling command used if error is generated
   ref = grid reference as it appears in an input script
     e.g. c_myCompute:grid:data[2], ditto for a fix
   nevery = frequency at which caller will access fix, not used if a compute
   return arguments:
     id = ID of compute or fix
     igrid = index of which grid in compute/fix (0 to N-1)
     idata = index of which data field in igrid (0 to N-1)
     index = index into data field (0 for vector, 1-N for column of array)
   method return = ArgInfo::COMPUTE or ArgInfo::FIX or -1 for neither
     caller decides what to do if arg is not a COMPUTE or FIX reference
------------------------------------------------------------------------- */

int utils::check_grid_reference(char *errstr, char *ref, int nevery, char *&id, int &igrid,
                                int &idata, int &index, LAMMPS *lmp)
{
  ArgInfo argi(ref, ArgInfo::COMPUTE | ArgInfo::FIX);
  index = argi.get_index1();
  auto name = argi.get_name();

  switch (argi.get_type()) {

    case ArgInfo::UNKNOWN: {
      lmp->error->all(FLERR, "%s grid reference %s is invalid", errstr, ref);
    } break;

      // compute value = c_ID

    case ArgInfo::COMPUTE: {

      // split name = idcompute:gname:dname into 3 strings

      auto words = parse_grid_id(FLERR, name, lmp->error);
      const auto &idcompute = words[0];
      const auto &gname = words[1];
      const auto &dname = words[2];

      auto icompute = lmp->modify->get_compute_by_id(idcompute);
      if (!icompute) lmp->error->all(FLERR, "{} compute ID {} not found", errstr, idcompute);
      if (icompute->pergrid_flag == 0)
        lmp->error->all(FLERR, "{} compute {} does not compute per-grid info", errstr, idcompute);

      int dim;
      igrid = icompute->get_grid_by_name(gname, dim);
      if (igrid < 0)
        lmp->error->all(FLERR, "{} compute {} does not recognize grid name {}", errstr, idcompute,
                        gname);

      int ncol;
      idata = icompute->get_griddata_by_name(igrid, dname, ncol);
      if (idata < 0)
        lmp->error->all(FLERR, "{} compute {} does not recognize data name {}", errstr, idcompute,
                        dname);

      if (argi.get_dim() == 0 && ncol)
        lmp->error->all(FLERR, "{} compute {} data {} is not per-grid vector", errstr, idcompute,
                        dname);
      if (argi.get_dim() && ncol == 0)
        lmp->error->all(FLERR, "{} compute {} data {} is not per-grid array", errstr, idcompute,
                        dname);
      if (argi.get_dim() && argi.get_index1() > ncol)
        lmp->error->all(FLERR, "{} compute {} array {} is accessed out-of-range", errstr, idcompute,
                        dname);

      id = utils::strdup(idcompute);
      return ArgInfo::COMPUTE;
    } break;

      // fix value = f_ID

    case ArgInfo::FIX: {

      // split name = idfix:gname:dname into 3 strings

      auto words = parse_grid_id(FLERR, name, lmp->error);
      const auto &idfix = words[0];
      const auto &gname = words[1];
      const auto &dname = words[2];

      auto ifix = lmp->modify->get_fix_by_id(idfix);
      if (!ifix) lmp->error->all(FLERR, "{} fix ID {} not found", errstr, idfix);
      if (ifix->pergrid_flag == 0)
        lmp->error->all(FLERR, "{} fix {} does not compute per-grid info", errstr, idfix);
      if (nevery % ifix->pergrid_freq)
        lmp->error->all(FLERR, "{} fix {} not computed at compatible time", errstr, idfix);

      int dim;
      igrid = ifix->get_grid_by_name(gname, dim);
      if (igrid < 0)
        lmp->error->all(FLERR, "{} fix {} does not recognize grid name {}", errstr, idfix, gname);

      int ncol;
      idata = ifix->get_griddata_by_name(igrid, dname, ncol);
      if (idata < 0)
        lmp->error->all(FLERR, "{} fix {} does not recognize data name {}", errstr, idfix, dname);

      if (argi.get_dim() == 0 && ncol)
        lmp->error->all(FLERR, "{} fix {} data {} is not per-grid vector", errstr, idfix, dname);
      if (argi.get_dim() > 0 && ncol == 0)
        lmp->error->all(FLERR, "{} fix {} data {} is not per-grid array", errstr, idfix, dname);
      if (argi.get_dim() > 0 && argi.get_index1() > ncol)
        lmp->error->all(FLERR, "{} fix {} array {} is accessed out-of-range", errstr, idfix, dname);

      id = utils::strdup(idfix);
      return ArgInfo::FIX;
    } break;
  }

  return -1;
}

/* ----------------------------------------------------------------------
   Parse grid reference into id:gridname:dataname
   return vector of 3 substrings
------------------------------------------------------------------------- */

std::vector<std::string> utils::parse_grid_id(const char *file, int line, const std::string &name,
                                              Error *error)
{
  auto words = Tokenizer(name, ":").as_vector();
  if (words.size() != 3) {
    if (error)
      error->all(file, line, "Grid ID {} must be 3 strings separated by 2 ':'characters", name);
    else
      return {"", "", ""};
  }

  return words;
}

/* ----------------------------------------------------------------------
   Make copy of string in new storage. Works like the (non-portable)
   C-style strdup() but also accepts a C++ string as argument.
------------------------------------------------------------------------- */

char *utils::strdup(const std::string &text)
{
  auto tmp = new char[text.size() + 1];
  strcpy(tmp, text.c_str());    // NOLINT
  return tmp;
}

/* ----------------------------------------------------------------------
   Return string converted to lowercase
------------------------------------------------------------------------- */

std::string utils::lowercase(const std::string &text)
{
  std::string converted(text);
  for (auto &c : converted) c = ::tolower(c);
  return converted;
}

/* ----------------------------------------------------------------------
   Return string converted to uppercase
------------------------------------------------------------------------- */

std::string utils::uppercase(const std::string &text)
{
  std::string converted(text);
  for (auto &c : converted) c = ::toupper(c);
  return converted;
}

/* ----------------------------------------------------------------------
   Return string without leading or trailing whitespace
------------------------------------------------------------------------- */

std::string utils::trim(const std::string &line)
{
  int beg = re_match(line.c_str(), "\\S+");
  int end = re_match(line.c_str(), "\\s+$");
  if (beg < 0) beg = 0;
  if (end < 0) end = line.size();

  return line.substr(beg, end - beg);
}

/* ----------------------------------------------------------------------
   Return string without trailing # comment
------------------------------------------------------------------------- */

std::string utils::trim_comment(const std::string &line)
{
  auto end = line.find('#');
  if (end != std::string::npos) { return line.substr(0, end); }
  return {line};
}

/* ----------------------------------------------------------------------
   Replace '*' with number and optional zero-padding
------------------------------------------------------------------------- */

std::string utils::star_subst(const std::string &name, bigint step, int pad)
{
  auto star = name.find('*');
  if (star == std::string::npos) return name;

  return fmt::format("{}{:0{}}{}", name.substr(0, star), step, pad, name.substr(star + 1));
}

/* ----------------------------------------------------------------------
   Remove accelerator style suffix from string
------------------------------------------------------------------------- */
std::string utils::strip_style_suffix(const std::string &style, LAMMPS *lmp)
{
  std::string newstyle = style;
  if (lmp->suffix_enable) {
    if (lmp->suffix) {
      if (utils::strmatch(style, fmt::format("/{}$", lmp->suffix))) {
        newstyle.resize(style.size() - strlen(lmp->suffix) - 1);
        return newstyle;
      }
    }
    if (lmp->suffix2) {
      if (utils::strmatch(style, fmt::format("/{}$", lmp->suffix2))) {
        newstyle.resize(style.size() - strlen(lmp->suffix2) - 1);
        return newstyle;
      }
    }
  }
  return newstyle;
}

/* ----------------------------------------------------------------------
   Replace UTF-8 encoded chars with known ASCII equivalents
------------------------------------------------------------------------- */

std::string utils::utf8_subst(const std::string &line)
{
  const auto *const in = (const unsigned char *) line.c_str();
  const int len = line.size();
  std::string out;

  for (int i = 0; i < len; ++i) {

    // UTF-8 2-byte character
    if ((in[i] & 0xe0U) == 0xc0U) {
      if ((i + 1) < len) {
        // NON-BREAKING SPACE (U+00A0)
        if ((in[i] == 0xc2U) && (in[i + 1] == 0xa0U)) out += ' ', ++i;
        // MODIFIER LETTER PLUS SIGN (U+02D6)
        if ((in[i] == 0xcbU) && (in[i + 1] == 0x96U)) out += '+', ++i;
        // MODIFIER LETTER MINUS SIGN (U+02D7)
        if ((in[i] == 0xcbU) && (in[i + 1] == 0x97U)) out += '-', ++i;
      }
      // UTF-8 3-byte character
    } else if ((in[i] & 0xf0U) == 0xe0U) {
      if ((i + 2) < len) {
        // EN QUAD (U+2000)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x80U)) out += ' ', i += 2;
        // EM QUAD (U+2001)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x81U)) out += ' ', i += 2;
        // EN SPACE (U+2002)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x82U)) out += ' ', i += 2;
        // EM SPACE (U+2003)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x83U)) out += ' ', i += 2;
        // THREE-PER-EM SPACE (U+2004)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x84U)) out += ' ', i += 2;
        // FOUR-PER-EM SPACE (U+2005)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x85U)) out += ' ', i += 2;
        // SIX-PER-EM SPACE (U+2006)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x86U)) out += ' ', i += 2;
        // FIGURE SPACE (U+2007)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x87U)) out += ' ', i += 2;
        // PUNCTUATION SPACE (U+2008)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x88U)) out += ' ', i += 2;
        // THIN SPACE (U+2009)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x89U)) out += ' ', i += 2;
        // HAIR SPACE (U+200A)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x8aU)) out += ' ', i += 2;
        // ZERO WIDTH SPACE (U+200B)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x8bU)) out += ' ', i += 2;
        // LEFT SINGLE QUOTATION MARK (U+2018)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x98U)) out += '\'', i += 2;
        // RIGHT SINGLE QUOTATION MARK (U+2019)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x99U)) out += '\'', i += 2;
        // LEFT DOUBLE QUOTATION MARK (U+201C)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x9cU)) out += '"', i += 2;
        // RIGHT DOUBLE QUOTATION MARK (U+201D)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x9dU)) out += '"', i += 2;
        // NARROW NO-BREAK SPACE (U+202F)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0xafU)) out += ' ', i += 2;
        // WORD JOINER (U+2060)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x81U) && (in[i + 2] == 0xa0U)) out += ' ', i += 2;
        // INVISIBLE SEPARATOR (U+2063)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x81U) && (in[i + 2] == 0xa3U)) out += ' ', i += 2;
        // INVISIBLE PLUS (U+2064)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x81U) && (in[i + 2] == 0xa4U)) out += '+', i += 2;
        // MINUS SIGN (U+2212)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x88U) && (in[i + 2] == 0x92U)) out += '-', i += 2;
        // ZERO WIDTH NO-BREAK SPACE (U+FEFF)
        if ((in[i] == 0xefU) && (in[i + 1] == 0xbbU) && (in[i + 2] == 0xbfU)) out += ' ', i += 2;
      }
      // UTF-8 4-byte character
    } else if ((in[i] & 0xf8U) == 0xf0U) {
      if ((i + 3) < len) { ; }
    } else
      out += in[i];
  }
  return out;
}

/* ----------------------------------------------------------------------
   return number of words
------------------------------------------------------------------------- */

size_t utils::count_words(const char *text)
{
  size_t count = 0;
  const char *buf = text;
  char c = *buf;

  while (c) {
    if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f') {
      c = *++buf;
      continue;
    };

    ++count;
    c = *++buf;

    while (c) {
      if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f') { break; }
      c = *++buf;
    }
  }

  return count;
}

/* ----------------------------------------------------------------------
   return number of words
------------------------------------------------------------------------- */

size_t utils::count_words(const std::string &text)
{
  return utils::count_words(text.c_str());
}

/* ----------------------------------------------------------------------
   Return number of words
------------------------------------------------------------------------- */

size_t utils::count_words(const std::string &text, const std::string &separators)
{
  size_t count = 0;
  size_t start = text.find_first_not_of(separators);

  while (start != std::string::npos) {
    size_t end = text.find_first_of(separators, start);
    ++count;

    if (end == std::string::npos) {
      return count;
    } else {
      start = text.find_first_not_of(separators, end + 1);
    }
  }
  return count;
}

/* ----------------------------------------------------------------------
   Trim comment from string and return number of words
------------------------------------------------------------------------- */

size_t utils::trim_and_count_words(const std::string &text, const std::string &separators)
{
  return utils::count_words(trim_comment(text), separators);
}

/* ----------------------------------------------------------------------
   combine words in vector to single string with separator added between words
------------------------------------------------------------------------- */
std::string utils::join_words(const std::vector<std::string> &words, const std::string &sep)
{
  std::string result;

  if (words.size() > 0) result = words[0];
  for (std::size_t i = 1; i < words.size(); ++i) result += sep + words[i];

  return result;
}

/* ----------------------------------------------------------------------
   Convert string into words on whitespace while handling single and
   double quotes.
------------------------------------------------------------------------- */
std::vector<std::string> utils::split_words(const std::string &text)
{
  std::vector<std::string> list;
  const char *buf = text.c_str();
  std::size_t beg = 0;
  std::size_t len = 0;
  std::size_t add = 0;
  char c = *buf;

  while (c) {
    // leading whitespace
    if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f') {
      c = *++buf;
      ++beg;
      continue;
    };
    len = 0;

  // handle escaped/quoted text.
  quoted:

    // handle single quote
    if (c == '\'') {
      ++beg;
      add = 1;
      c = *++buf;
      while (((c != '\'') && (c != '\0')) || ((c == '\\') && (buf[1] == '\''))) {
        if ((c == '\\') && (buf[1] == '\'')) {
          ++buf;
          ++len;
        }
        c = *++buf;
        ++len;
      }
      if (c != '\'') ++len;
      c = *++buf;

      // handle triple double quotation marks
    } else if ((c == '"') && (buf[1] == '"') && (buf[2] == '"') && (buf[3] != '"')) {
      len = 3;
      add = 1;
      buf += 3;
      c = *buf;

      // handle double quote
    } else if (c == '"') {
      ++beg;
      add = 1;
      c = *++buf;
      while (((c != '"') && (c != '\0')) || ((c == '\\') && (buf[1] == '"'))) {
        if ((c == '\\') && (buf[1] == '"')) {
          ++buf;
          ++len;
        }
        c = *++buf;
        ++len;
      }
      if (c != '"') ++len;
      c = *++buf;
    }

    // unquoted
    while (true) {
      if ((c == '\'') || (c == '"')) goto quoted;
      // skip escaped quote
      if ((c == '\\') && ((buf[1] == '\'') || (buf[1] == '"'))) {
        ++buf;
        ++len;
        c = *++buf;
        ++len;
      }
      if ((c == ' ') || (c == '\t') || (c == '\r') || (c == '\n') || (c == '\f') || (c == '\0')) {
        list.push_back(text.substr(beg, len));
        beg += len + add;
        break;
      }
      c = *++buf;
      ++len;
    }
  }
  return list;
}

/* ----------------------------------------------------------------------
   Convert multi-line string into lines
------------------------------------------------------------------------- */

std::vector<std::string> utils::split_lines(const std::string &text)
{
  return Tokenizer(text, "\r\n").as_vector();
}

/* ----------------------------------------------------------------------
   Return whether string is a valid integer number
------------------------------------------------------------------------- */

bool utils::is_integer(const std::string &str)
{
  if (str.empty()) return false;

  return strmatch(str, "^[+-]?\\d+$");
}

/* ----------------------------------------------------------------------
   Return whether string is a valid floating-point number
------------------------------------------------------------------------- */

bool utils::is_double(const std::string &str)
{
  if (str.empty()) return false;

  return strmatch(str, "^[+-]?\\d+\\.?\\d*$") ||
      strmatch(str, "^[+-]?\\d+\\.?\\d*[eE][+-]?\\d+$") || strmatch(str, "^[+-]?\\d*\\.?\\d+$") ||
      strmatch(str, "^[+-]?\\d*\\.?\\d+[eE][+-]?\\d+$");
}

/* ----------------------------------------------------------------------
   Return whether string is a valid ID string
------------------------------------------------------------------------- */

bool utils::is_id(const std::string &str)
{
  if (str.empty()) return false;

  for (const auto &c : str) {
    if (isalnum(c) || (c == '_')) continue;
    return false;
  }
  return true;
}

/* ----------------------------------------------------------------------
   Check whether string is a valid type or type label string
------------------------------------------------------------------------- */

int utils::is_type(const std::string &str)
{
  if (str.empty()) return -1;

  bool numeric = true;
  int nstar = 0;
  for (const auto &c : str) {
    if (isdigit(c)) continue;
    if (c == '*') {
      ++nstar;
      continue;
    }
    numeric = false;
  }
  if (numeric && (nstar < 2)) return 0;

  // TODO: the first two checks below are not really needed with this function.
  // If a type label has at least one character that is not a digit or '*'
  // it can be identified by this function as type label due to the check above.
  // Whitespace and multi-byte characters are not allowed.
  if (isdigit(str[0]) || (str[0] == '*') || (str[0] == '#')) return -1;
  if (str.find_first_of(" \t\r\n\f") != std::string::npos) return -1;
  if (has_utf8(utf8_subst(str))) return -1;
  return 1;
}

/* ----------------------------------------------------------------------
   try to find potential file as specified by name
   search current directory and the LAMMPS_POTENTIALS directory if
   specified
------------------------------------------------------------------------- */

std::string utils::get_potential_file_path(const std::string &path)
{
  if (platform::file_is_readable(path)) {
    return path;
  } else {
    for (const auto &dir : platform::list_pathenv("LAMMPS_POTENTIALS")) {
      auto pot = platform::path_basename(path);
      auto filepath = platform::path_join(dir, pot);
      if (platform::file_is_readable(filepath)) return filepath;
    }
  }
  return "";
}

/* ----------------------------------------------------------------------
   read first line of potential file
   if it has a DATE field, return the following word
------------------------------------------------------------------------- */

std::string utils::get_potential_date(const std::string &path, const std::string &potential_name)
{
  TextFileReader reader(path, potential_name);
  reader.ignore_comments = false;

  char *line = reader.next_line();
  if (line == nullptr) return "";
  Tokenizer words(line);
  while (words.has_next()) {
    if (words.next() == "DATE:") {
      if (words.has_next()) return words.next();
    }
  }
  return "";
}

/* ----------------------------------------------------------------------
   read first line of potential file
   if it has UNITS field, return following word
------------------------------------------------------------------------- */

std::string utils::get_potential_units(const std::string &path, const std::string &potential_name)
{
  TextFileReader reader(path, potential_name);
  reader.ignore_comments = false;

  char *line = reader.next_line();
  if (line == nullptr) return "";
  Tokenizer words(line);
  while (words.has_next()) {
    if (words.next() == "UNITS:") {
      if (words.has_next()) return words.next();
    }
  }
  return "";
}

/* ----------------------------------------------------------------------
   return bitmask of supported conversions for a given property
------------------------------------------------------------------------- */

int utils::get_supported_conversions(const int property)
{
  if (property == ENERGY)
    return METAL2REAL | REAL2METAL;
  else
    return NOCONVERT;
}

/* ----------------------------------------------------------------------
   return conversion factor for a given property and conversion setting
   return 0.0 if unknown.
------------------------------------------------------------------------- */

double utils::get_conversion_factor(const int property, const int conversion)
{
  if (property == ENERGY) {
    if (conversion == NOCONVERT) {
      return 1.0;
    } else if (conversion == METAL2REAL) {
      return 23.060549;
    } else if (conversion == REAL2METAL) {
      return 1.0 / 23.060549;
    }
  }
  return 0.0;
}

/* ----------------------------------------------------------------------
   open a potential file as specified by name
   if fails, search in dir specified by env variable LAMMPS_POTENTIALS
------------------------------------------------------------------------- */

FILE *utils::open_potential(const std::string &name, LAMMPS *lmp, int *auto_convert)
{
  auto error = lmp->error;
  auto me = lmp->comm->me;

  std::string filepath = get_potential_file_path(name);

  if (!filepath.empty()) {
    std::string unit_style = lmp->update->unit_style;
    std::string date = get_potential_date(filepath, "potential");
    std::string units = get_potential_units(filepath, "potential");

    if (!date.empty() && (me == 0))
      logmesg(lmp, "Reading potential file {} with DATE: {}\n", name, date);

    if (auto_convert == nullptr) {
      if (!units.empty() && (units != unit_style) && (me == 0)) {
        error->one(FLERR, "Potential file {} requires {} units but {} units are in use", name,
                   units, unit_style);
        return nullptr;
      }
    } else {
      if (units.empty() || units == unit_style) {
        *auto_convert = NOCONVERT;
      } else {
        if ((units == "metal") && (unit_style == "real") && (*auto_convert & METAL2REAL)) {
          *auto_convert = METAL2REAL;
        } else if ((units == "real") && (unit_style == "metal") && (*auto_convert & REAL2METAL)) {
          *auto_convert = REAL2METAL;
        } else {
          error->one(FLERR, "Potential file {} requires {} units but {} units are in use", name,
                     units, unit_style);
          return nullptr;
        }
      }
      if ((*auto_convert != NOCONVERT) && (me == 0))
        error->warning(FLERR, "Converting potential file in {} units to {} units", units,
                       unit_style);
    }
    return fopen(filepath.c_str(), "r");
  }
  return nullptr;
}

/* ----------------------------------------------------------------------
   convert a timespec ([[HH:]MM:]SS) to seconds
   the strings "off" and "unlimited" result in -1.0;
------------------------------------------------------------------------- */

double utils::timespec2seconds(const std::string &timespec)
{
  double vals[3];
  int i = 0;

  // first handle allowed textual inputs
  if (timespec == "off") return -1.0;
  if (timespec == "unlimited") return -1.0;

  vals[0] = vals[1] = vals[2] = 0;

  ValueTokenizer values(timespec, ":");

  try {
    for (i = 0; i < 3; i++) {
      if (!values.has_next()) break;
      vals[i] = values.next_int();
    }
  } catch (TokenizerException &) {
    return -1.0;
  }

  if (i == 3)
    return (vals[0] * 60 + vals[1]) * 60 + vals[2];
  else if (i == 2)
    return vals[0] * 60 + vals[1];
  return vals[0];
}

/* ----------------------------------------------------------------------
   convert a LAMMPS version date (1Jan01) to a number
------------------------------------------------------------------------- */

int utils::date2num(const std::string &date)
{
  std::size_t found = date.find_first_not_of("0123456789 ");
  int num = strtol(date.substr(0, found).c_str(), nullptr, 10);
  auto month = date.substr(found);
  found = month.find_first_of("0123456789 ");
  num += strtol(month.substr(found).c_str(), nullptr, 10) * 10000;
  if (num < 1000000) num += 20000000;

  if (strmatch(month, "^Jan"))
    num += 100;
  else if (strmatch(month, "^Feb"))
    num += 200;
  else if (strmatch(month, "^Mar"))
    num += 300;
  else if (strmatch(month, "^Apr"))
    num += 400;
  else if (strmatch(month, "^May"))
    num += 500;
  else if (strmatch(month, "^Jun"))
    num += 600;
  else if (strmatch(month, "^Jul"))
    num += 700;
  else if (strmatch(month, "^Aug"))
    num += 800;
  else if (strmatch(month, "^Sep"))
    num += 900;
  else if (strmatch(month, "^Oct"))
    num += 1000;
  else if (strmatch(month, "^Nov"))
    num += 1100;
  else if (strmatch(month, "^Dec"))
    num += 1200;
  return num;
}

/* ----------------------------------------------------------------------
   get formatted string of current date from fmtlib
------------------------------------------------------------------------- */

std::string utils::current_date()
{
  time_t tv = time(nullptr);
  std::tm today = fmt::localtime(tv);
  return fmt::format("{:%Y-%m-%d}", today);
}

/* ----------------------------------------------------------------------
   binary search in vector of ascending doubles
------------------------------------------------------------------------- */

int utils::binary_search(const double needle, const int n, const double *haystack)
{
  int lo = 0;
  int hi = n - 1;

  if (needle < haystack[lo]) return lo;
  if (needle >= haystack[hi]) return hi;

  // ensure haystack[lo] <= needle < haystack[hi] at every iteration
  // done when lo,hi are adjacent

  int index = (lo + hi) / 2;
  while (lo < hi - 1) {
    if (needle < haystack[index])
      hi = index;
    else if (needle >= haystack[index])
      lo = index;
    index = (lo + hi) / 2;
  }

  return index;
}

/* ----------------------------------------------------------------------
 * Merge sort part 1: Loop over sublists doubling in size with each iteration.
 * Pre-sort small sublists with insertion sort for better overall performance.
------------------------------------------------------------------------- */

void utils::merge_sort(int *index, int num, void *ptr, int (*comp)(int, int, void *))
{
  if (num < 2) return;

  int chunk, i, j;

  // do insertion sort on chunks of up to 64 elements

  chunk = 64;
  for (i = 0; i < num; i += chunk) {
    j = (i + chunk > num) ? num - i : chunk;
    insertion_sort(index + i, j, ptr, comp);
  }

  // already done?

  if (chunk >= num) return;

  // continue with merge sort on the pre-sorted chunks.
  // we need an extra buffer for temporary storage and two
  // pointers to operate on, so we can swap the pointers
  // rather than copying to the hold buffer in each pass

  int *buf = new int[num];
  int *dest = index;
  int *hold = buf;

  while (chunk < num) {
    int m;

    // swap hold and destination buffer

    int *tmp = dest;
    dest = hold;
    hold = tmp;

    // merge from hold array to destination array

    for (i = 0; i < num - 1; i += 2 * chunk) {
      j = i + 2 * chunk;
      if (j > num) j = num;
      m = i + chunk;
      if (m > num) m = num;
      do_merge(dest, hold, i, m, m, j, ptr, comp);
    }

    // copy all indices not handled by the chunked merge sort loop

    for (; i < num; i++) dest[i] = hold[i];
    chunk *= 2;
  }

  // if the final sorted data is in buf, copy back to index

  if (dest == buf) memcpy(index, buf, sizeof(int) * num);

  delete[] buf;
}

/* ------------------------------------------------------------------ */

/* ----------------------------------------------------------------------
 * Merge sort part 2: Insertion sort for pre-sorting of small chunks
------------------------------------------------------------------------- */

void insertion_sort(int *index, int num, void *ptr, int (*comp)(int, int, void *))
{
  if (num < 2) return;
  for (int i = 1; i < num; ++i) {
    int tmp = index[i];
    for (int j = i - 1; j >= 0; --j) {
      if ((*comp)(index[j], tmp, ptr) > 0) {
        index[j + 1] = index[j];
      } else {
        index[j + 1] = tmp;
        break;
      }
      if (j == 0) index[0] = tmp;
    }
  }
}

/* ----------------------------------------------------------------------
 * Merge sort part 3: Merge two sublists
------------------------------------------------------------------------- */

static void do_merge(int *idx, int *buf, int llo, int lhi, int rlo, int rhi, void *ptr,
                     int (*comp)(int, int, void *))
{
  int i = llo;
  int l = llo;
  int r = rlo;
  while ((l < lhi) && (r < rhi)) {
    if ((*comp)(buf[l], buf[r], ptr) < 0)
      idx[i++] = buf[l++];
    else
      idx[i++] = buf[r++];
  }

  while (l < lhi) idx[i++] = buf[l++];
  while (r < rhi) idx[i++] = buf[r++];
}

/* ------------------------------------------------------------------ */

extern "C" {

/* Typedef'd pointer to get abstract datatype. */
typedef struct regex_t *re_t;
typedef struct regex_context_t *re_ctx_t;

/* Compile regex string pattern to a regex_t-array. */
static re_t re_compile(re_ctx_t context, const char *pattern);

/* Find matches of the compiled pattern inside text. */
static int re_matchp(const char *text, re_t pattern, int *matchlen);

/* Definitions: */

#define MAX_REGEXP_OBJECTS 256 /* Max number of regex symbols in expression. */
#define MAX_CHAR_CLASS_LEN 256 /* Max length of character-class buffer in.   */

enum {
  RX_UNUSED,
  RX_DOT,
  RX_BEGIN,
  RX_END,
  RX_QUESTIONMARK,
  RX_STAR,
  RX_PLUS,
  RX_CHAR,
  RX_CHAR_CLASS,
  RX_INV_CHAR_CLASS,
  RX_DIGIT,
  RX_NOT_DIGIT,
  RX_INTEGER,
  RX_NOT_INTEGER,
  RX_FLOAT,
  RX_NOT_FLOAT,
  RX_ALPHA,
  RX_NOT_ALPHA,
  RX_WHITESPACE,
  RX_NOT_WHITESPACE /*, BRANCH */
};

typedef struct regex_t {
  unsigned char type; /* CHAR, STAR, etc.                      */
  union {
    unsigned char ch;   /*      the character itself             */
    unsigned char *ccl; /*  OR  a pointer to characters in class */
  } u;
} regex_t;

typedef struct regex_context_t {
  /* MAX_REGEXP_OBJECTS is the max number of symbols in the expression.
       MAX_CHAR_CLASS_LEN determines the size of buffer for chars in all char-classes in the expression. */
  regex_t re_compiled[MAX_REGEXP_OBJECTS];
  unsigned char ccl_buf[MAX_CHAR_CLASS_LEN];
} regex_context_t;

int re_match(const char *text, const char *pattern)
{
  regex_context_t context;
  int dummy;
  return re_matchp(text, re_compile(&context, pattern), &dummy);
}

int re_find(const char *text, const char *pattern, int *matchlen)
{
  regex_context_t context;
  return re_matchp(text, re_compile(&context, pattern), matchlen);
}

/* Private function declarations: */
static int matchpattern(regex_t *pattern, const char *text, int *matchlen);
static int matchcharclass(char c, const char *str);
static int matchstar(regex_t p, regex_t *pattern, const char *text, int *matchlen);
static int matchplus(regex_t p, regex_t *pattern, const char *text, int *matchlen);
static int matchone(regex_t p, char c);
static int matchdigit(char c);
static int matchint(char c);
static int matchfloat(char c);
static int matchalpha(char c);
static int matchwhitespace(char c);
static int matchmetachar(char c, const char *str);
static int matchrange(char c, const char *str);
static int matchdot(char c);
static int ismetachar(char c);

/* Semi-public functions: */
int re_matchp(const char *text, re_t pattern, int *matchlen)
{
  *matchlen = 0;
  if (pattern != nullptr) {
    if (pattern[0].type == RX_BEGIN) {
      return ((matchpattern(&pattern[1], text, matchlen)) ? 0 : -1);
    } else {
      int idx = -1;

      do {
        idx += 1;

        if (matchpattern(pattern, text, matchlen)) {
          if (text[0] == '\0') return -1;

          return idx;
        }
      } while (*text++ != '\0');
    }
  }
  return -1;
}

re_t re_compile(re_ctx_t context, const char *pattern)
{
  regex_t *const re_compiled = context->re_compiled;
  unsigned char *const ccl_buf = context->ccl_buf;
  int ccl_bufidx = 1;

  char c;    /* current char in pattern   */
  int i = 0; /* index into pattern        */
  int j = 0; /* index into re_compiled    */

  while (pattern[i] != '\0' && (j + 1 < MAX_REGEXP_OBJECTS)) {
    c = pattern[i];

    switch (c) {
        /* Meta-characters: */
      case '^': {
        re_compiled[j].type = RX_BEGIN;
      } break;
      case '$': {
        re_compiled[j].type = RX_END;
      } break;
      case '.': {
        re_compiled[j].type = RX_DOT;
      } break;
      case '*': {
        re_compiled[j].type = RX_STAR;
      } break;
      case '+': {
        re_compiled[j].type = RX_PLUS;
      } break;
      case '?': {
        re_compiled[j].type = RX_QUESTIONMARK;
      } break;

        /* Escaped character-classes (\s \w ...): */
      case '\\': {
        if (pattern[i + 1] != '\0') {
          /* Skip the escape-char '\\' */
          i += 1;
          /* ... and check the next */
          switch (pattern[i]) {
              /* Meta-character: */
            case 'd': {
              re_compiled[j].type = RX_DIGIT;
            } break;
            case 'D': {
              re_compiled[j].type = RX_NOT_DIGIT;
            } break;
            case 'i': {
              re_compiled[j].type = RX_INTEGER;
            } break;
            case 'I': {
              re_compiled[j].type = RX_NOT_INTEGER;
            } break;
            case 'f': {
              re_compiled[j].type = RX_FLOAT;
            } break;
            case 'F': {
              re_compiled[j].type = RX_NOT_FLOAT;
            } break;
            case 'w': {
              re_compiled[j].type = RX_ALPHA;
            } break;
            case 'W': {
              re_compiled[j].type = RX_NOT_ALPHA;
            } break;
            case 's': {
              re_compiled[j].type = RX_WHITESPACE;
            } break;
            case 'S': {
              re_compiled[j].type = RX_NOT_WHITESPACE;
            } break;

              /* Escaped character, e.g. '.' or '$' */
            default: {
              re_compiled[j].type = RX_CHAR;
              re_compiled[j].u.ch = pattern[i];
            } break;
          }
        }
        /* '\\' as last char in pattern -> invalid regular expression. */
      } break;

        /* Character class: */
      case '[': {
        /* Remember where the char-buffer starts. */
        int buf_begin = ccl_bufidx;

        /* Look-ahead to determine if negated */
        if (pattern[i + 1] == '^') {
          re_compiled[j].type = RX_INV_CHAR_CLASS;
          i += 1;                  /* Increment i to avoid including '^' in the char-buffer */
          if (pattern[i + 1] == 0) /* incomplete pattern, missing non-zero char after '^' */
          {
            return nullptr;
          }
        } else {
          re_compiled[j].type = RX_CHAR_CLASS;
        }

        /* Copy characters inside [..] to buffer */
        while ((pattern[++i] != ']') && (pattern[i] != '\0')) {
          /* Missing ] */
          if (pattern[i] == '\\') {
            if (ccl_bufidx >= MAX_CHAR_CLASS_LEN - 1) { return nullptr; }
            if (pattern[i + 1] == 0) /* incomplete pattern, missing non-zero char after '\\' */
            {
              return nullptr;
            }
            ccl_buf[ccl_bufidx++] = pattern[i++];
          } else if (ccl_bufidx >= MAX_CHAR_CLASS_LEN) {
            return nullptr;
          }
          ccl_buf[ccl_bufidx++] = pattern[i];
        }
        if (ccl_bufidx >= MAX_CHAR_CLASS_LEN) {
          /* Catches cases such as [00000000000000000000000000000000000000][ */
          return nullptr;
        }
        /* Null-terminate string end */
        ccl_buf[ccl_bufidx++] = 0;
        re_compiled[j].u.ccl = &ccl_buf[buf_begin];
      } break;

        /* Other characters: */
      default: {
        re_compiled[j].type = RX_CHAR;
        re_compiled[j].u.ch = c;
      } break;
    }
    /* no buffer-out-of-bounds access on invalid patterns -
     * see https://github.com/kokke/tiny-regex-c/commit/1a279e04014b70b0695fba559a7c05d55e6ee90b */
    if (pattern[i] == 0) { return nullptr; }

    i += 1;
    j += 1;
  }
  /* 'RX_UNUSED' is a sentinel used to indicate end-of-pattern */
  re_compiled[j].type = RX_UNUSED;

  return (re_t) re_compiled;
}

/* Private functions: */
static int matchdigit(char c)
{
  return isdigit(c);
}

static int matchint(char c)
{
  return (matchdigit(c) || (c == '-') || (c == '+'));
}

static int matchfloat(char c)
{
  return (matchint(c) || (c == '.') || (c == 'e') || (c == 'E'));
}

static int matchalpha(char c)
{
  return isalpha(c);
}

static int matchwhitespace(char c)
{
  return isspace(c);
}

static int matchalphanum(char c)
{
  return ((c == '_') || matchalpha(c) || matchdigit(c));
}

static int matchrange(char c, const char *str)
{
  return ((c != '-') && (str[0] != '\0') && (str[0] != '-') && (str[1] == '-') &&
          (str[1] != '\0') && (str[2] != '\0') && ((c >= str[0]) && (c <= str[2])));
}

static int matchdot(char c)
{
#if defined(RE_DOT_MATCHES_NEWLINE) && (RE_DOT_MATCHES_NEWLINE == 1)
  (void) c;
  return 1;
#else
  return c != '\n' && c != '\r';
#endif
}

static int ismetachar(char c)
{
  return ((c == 's') || (c == 'S') || (c == 'w') || (c == 'W') || (c == 'd') || (c == 'D'));
}

static int matchmetachar(char c, const char *str)
{
  switch (str[0]) {
    case 'd':
      return matchdigit(c);
    case 'D':
      return !matchdigit(c);
    case 'i':
      return matchint(c);
    case 'I':
      return !matchint(c);
    case 'f':
      return matchfloat(c);
    case 'F':
      return !matchfloat(c);
    case 'w':
      return matchalphanum(c);
    case 'W':
      return !matchalphanum(c);
    case 's':
      return matchwhitespace(c);
    case 'S':
      return !matchwhitespace(c);
    default:
      return (c == str[0]);
  }
}

static int matchcharclass(char c, const char *str)
{
  do {
    if (matchrange(c, str)) {
      return 1;
    } else if (str[0] == '\\') {
      /* Escape-char: increment str-ptr and match on next char */
      str += 1;
      if (matchmetachar(c, str)) {
        return 1;
      } else if ((c == str[0]) && !ismetachar(c)) {
        return 1;
      }
    } else if (c == str[0]) {
      if (c == '-') {
        return ((str[-1] == '\0') || (str[1] == '\0'));
      } else {
        return 1;
      }
    }
  } while (*str++ != '\0');

  return 0;
}

static int matchone(regex_t p, char c)
{
  switch (p.type) {
    case RX_DOT:
      return matchdot(c);
    case RX_CHAR_CLASS:
      return matchcharclass(c, (const char *) p.u.ccl);
    case RX_INV_CHAR_CLASS:
      return !matchcharclass(c, (const char *) p.u.ccl);
    case RX_DIGIT:
      return matchdigit(c);
    case RX_NOT_DIGIT:
      return !matchdigit(c);
    case RX_INTEGER:
      return matchint(c);
    case RX_NOT_INTEGER:
      return !matchint(c);
    case RX_FLOAT:
      return matchfloat(c);
    case RX_NOT_FLOAT:
      return !matchfloat(c);
    case RX_ALPHA:
      return matchalphanum(c);
    case RX_NOT_ALPHA:
      return !matchalphanum(c);
    case RX_WHITESPACE:
      return matchwhitespace(c);
    case RX_NOT_WHITESPACE:
      return !matchwhitespace(c);
    default:
      return (p.u.ch == c);
  }
}

static int matchstar(regex_t p, regex_t *pattern, const char *text, int *matchlen)
{
  int prelen = *matchlen;
  const char *prepos = text;
  while ((text[0] != '\0') && matchone(p, *text)) {
    text++;
    (*matchlen)++;
  }
  while (text >= prepos) {
    if (matchpattern(pattern, text--, matchlen)) return 1;
    (*matchlen)--;
  }

  *matchlen = prelen;
  return 0;
}

static int matchplus(regex_t p, regex_t *pattern, const char *text, int *matchlen)
{
  const char *prepos = text;
  while ((text[0] != '\0') && matchone(p, *text)) {
    text++;
    (*matchlen)++;
  }
  while (text > prepos) {
    if (matchpattern(pattern, text--, matchlen)) return 1;
    (*matchlen)--;
  }
  return 0;
}

static int matchquestion(regex_t p, regex_t *pattern, const char *text, int *matchlen)
{
  if (p.type == RX_UNUSED) return 1;
  if (matchpattern(pattern, text, matchlen)) return 1;
  if (*text && matchone(p, *text++)) {
    if (matchpattern(pattern, text, matchlen)) {
      (*matchlen)++;
      return 1;
    }
  }
  return 0;
}

/* Iterative matching */
static int matchpattern(regex_t *pattern, const char *text, int *matchlen)
{
  int pre = *matchlen;
  do {
    if ((pattern[0].type == RX_UNUSED) || (pattern[1].type == RX_QUESTIONMARK)) {
      return matchquestion(pattern[0], &pattern[2], text, matchlen);
    } else if (pattern[1].type == RX_STAR) {
      return matchstar(pattern[0], &pattern[2], text, matchlen);
    } else if (pattern[1].type == RX_PLUS) {
      return matchplus(pattern[0], &pattern[2], text, matchlen);
    } else if ((pattern[0].type == RX_END) && pattern[1].type == RX_UNUSED) {
      return (text[0] == '\0');
    }
    (*matchlen)++;
  } while ((text[0] != '\0') && matchone(*pattern++, *text++));

  *matchlen = pre;
  return 0;
}
}
