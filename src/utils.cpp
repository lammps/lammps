/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "utils.h"
#include <cstring>
#include <cstdlib>
#include <cerrno>
#include "lammps.h"
#include "error.h"
#include "tokenizer.h"
#include "text_file_reader.h"
#include "fmt/format.h"

#if defined(__linux__)
#include <unistd.h>  // for readlink
#endif

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

extern "C"
{
  /** Match text against a (simplified) regular expression
   * (regexp will be compiled automatically). */
  static int  re_match(const char *text, const char *pattern);
}

using namespace LAMMPS_NS;

/** More flexible and specific matching of a string against a pattern.
 *  This function is supposed to be a more safe, more specific and
 *  simple to use API to find pattern matches. The purpose is to replace
 *  uses of either strncmp() or strstr() in the code base to find
 *  substrings safely. With strncmp() finding prefixes, the number of
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
  const int pos = re_match(text.c_str(),pattern.c_str());
  return (pos >= 0);
}

/* This simplifies the repetitive task of outputting some
 * message to both the screen and/or the log file. In combination
 * with using fmt::format(), which returns the formatted text
 * in a std::string() instance, this can be used to reduce
 * operations previously requiring several lines of code to
 * a single statement. */

void utils::logmesg(LAMMPS *lmp, const std::string &mesg)
{
  if (lmp->screen)  fputs(mesg.c_str(), lmp->screen);
  if (lmp->logfile) fputs(mesg.c_str(), lmp->logfile);
}

/* define this here, so we won't have to include the headers
   everywhere and utils.h will more likely be included anyway. */

std::string utils::getsyserror()
{
  return std::string(strerror(errno));
}

/*
 * On Linux the folder /proc/self/fd holds symbolic links to the actual
 * pathnames associated with each open file descriptor of the current process.
 */
const char *utils::guesspath(char *buf, int len, FILE *fp)
{
  memset(buf,0,len);

#if defined(__linux__)
  int fd = fileno(fp);
  // get pathname from /proc or copy (unknown)
  if (readlink(fmt::format("/proc/self/fd/{}",fd).c_str(),buf,len-1) <= 0)
    strncpy(buf,"(unknown)",len-1);
#else
  strncpy(buf,"(unknown)",len-1);
#endif
  return buf;
}

#define MAXPATHLENBUF 1024
/* like fgets() but aborts with an error or EOF is encountered */
void utils::sfgets(const char *srcname, int srcline, char *s, int size,
                   FILE *fp, const char *filename, Error *error)
{
  char *rv = fgets(s,size,fp);
  if (rv == NULL) { // something went wrong
    char buf[MAXPATHLENBUF];
    std::string errmsg;

    // try to figure out the file name from the file pointer
    if (!filename)
      filename = guesspath(buf,MAXPATHLENBUF,fp);

    if (feof(fp)) {
      errmsg = "Unexpected end of file while reading file '";
    } else if (ferror(fp)) {
      errmsg = "Unexpected error while reading file '";
    } else {
      errmsg = "Unexpected short read while reading file '";
    }
    errmsg += filename;
    errmsg += "'";

    if (error) error->one(srcname,srcline,errmsg);
    if (s) *s = '\0'; // truncate string to empty in case error is NULL
  }
  return;
}

/* like fread() but aborts with an error or EOF is encountered */
void utils::sfread(const char *srcname, int srcline, void *s, size_t size,
                   size_t num, FILE *fp, const char *filename, Error *error)
{
  size_t rv = fread(s,size,num,fp);
  if (rv != num) { // something went wrong
    char buf[MAXPATHLENBUF];
    std::string errmsg;

    // try to figure out the file name from the file pointer
    if (!filename)
      filename = guesspath(buf,MAXPATHLENBUF,fp);

    if (feof(fp)) {
      errmsg = "Unexpected end of file while reading file '";
    } else if (ferror(fp)) {
      errmsg = "Unexpected error while reading file '";
    } else {
      errmsg = "Unexpected short read while reading file '";
    }
    errmsg += filename;
    errmsg += "'";

    if (error) error->one(srcname,srcline,errmsg);
  }
  return;
}

/* ------------------------------------------------------------------ */

std::string utils::check_packages_for_style(const std::string &style,
                                            const std::string &name,
                                            LAMMPS *lmp)
{
  std::string errmsg = "Unrecognized " + style + " style '" + name + "'";
  const char *pkg = lmp->match_style(style.c_str(),name.c_str());

  if (pkg) {
    errmsg += fmt::format(" is part of the {} package",pkg);
    if (lmp->is_installed_pkg(pkg))
      errmsg += ", but seems to be missing because of a dependency";
    else
      errmsg += " which is not enabled in this LAMMPS binary.";
  }
  return errmsg;
}


/* ----------------------------------------------------------------------
   read a floating point value from a string
   generate an error if not a legitimate floating point value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

double utils::numeric(const char *file, int line, const char *str,
                      bool do_abort, LAMMPS *lmp)
{
  int n = 0;

  if (str) n = strlen(str);
  if (n == 0) {
    if (do_abort)
      lmp->error->one(file,line,"Expected floating point parameter instead of"
                      " NULL or empty string in input script or data file");
    else
      lmp->error->all(file,line,"Expected floating point parameter instead of"
                      " NULL or empty string in input script or data file");
  }

  for (int i = 0; i < n; i++) {
    if (isdigit(str[i])) continue;
    if (str[i] == '-' || str[i] == '+' || str[i] == '.') continue;
    if (str[i] == 'e' || str[i] == 'E') continue;
    std::string msg("Expected floating point parameter instead of '");
    msg += str;
    msg += "' in input script or data file";
    if (do_abort)
      lmp->error->one(file,line,msg);
    else
      lmp->error->all(file,line,msg);
  }

  return atof(str);
}

/* ----------------------------------------------------------------------
   read an integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

int utils::inumeric(const char *file, int line, const char *str,
                    bool do_abort, LAMMPS *lmp)
{
  int n = 0;

  if (str) n = strlen(str);
  if (n == 0) {
    if (do_abort)
      lmp->error->one(file,line,"Expected integer parameter instead of "
                      "NULL or empty string in input script or data file");
    else
      lmp->error->all(file,line,"Expected integer parameter instead of "
                      "NULL or empty string in input script or data file");
  }

  for (int i = 0; i < n; i++) {
    if (isdigit(str[i]) || str[i] == '-' || str[i] == '+') continue;
    std::string msg("Expected integer parameter instead of '");
    msg += str;
    msg += "' in input script or data file";
    if (do_abort)
      lmp->error->one(file,line,msg);
    else
      lmp->error->all(file,line,msg);
  }

  return atoi(str);
}

/* ----------------------------------------------------------------------
   read a big integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

bigint utils::bnumeric(const char *file, int line, const char *str,
                       bool do_abort, LAMMPS *lmp)
{
  int n = 0;

  if (str) n = strlen(str);
  if (n == 0) {
    if (do_abort)
      lmp->error->one(file,line,"Expected integer parameter instead of "
                      "NULL or empty string in input script or data file");
    else
      lmp->error->all(file,line,"Expected integer parameter instead of "
                      "NULL or empty string in input script or data file");
  }

  for (int i = 0; i < n; i++) {
    if (isdigit(str[i]) || str[i] == '-' || str[i] == '+') continue;
    std::string msg("Expected integer parameter instead of '");
    msg += str;
    msg += "' in input script or data file";
    if (do_abort)
      lmp->error->one(file,line,msg);
    else
      lmp->error->all(file,line,msg);
  }

  return ATOBIGINT(str);
}

/* ----------------------------------------------------------------------
   read a tag integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

tagint utils::tnumeric(const char *file, int line, const char *str,
                       bool do_abort, LAMMPS *lmp)
{
  int n = 0;

  if (str) n = strlen(str);
  if (n == 0) {
    if (do_abort)
      lmp->error->one(file,line,"Expected integer parameter instead of "
                      "NULL or empty string in input script or data file");
    else
      lmp->error->all(file,line,"Expected integer parameter instead of "
                      "NULL or empty string in input script or data file");
  }

  for (int i = 0; i < n; i++) {
    if (isdigit(str[i]) || str[i] == '-' || str[i] == '+') continue;
    std::string msg("Expected integer parameter instead of '");
    msg += str;
    msg += "' in input script or data file";
    if (do_abort)
      lmp->error->one(file,line,msg);
    else
      lmp->error->all(file,line,msg);
  }

  return ATOTAGINT(str);
}

/* ----------------------------------------------------------------------
   Return string without leading or trailing whitespace
------------------------------------------------------------------------- */

std::string utils::trim(const std::string & line) {
  int beg = re_match(line.c_str(),"\\S+");
  int end = re_match(line.c_str(),"\\s+$");
  if (beg < 0) beg = 0;
  if (end < 0) end = line.size();

  return line.substr(beg,end-beg);
}

/* ----------------------------------------------------------------------
   Return string without trailing # comment
------------------------------------------------------------------------- */

std::string utils::trim_comment(const std::string & line) {
  auto end = line.find_first_of("#");
  if (end != std::string::npos) {
    return line.substr(0, end);
  }
  return std::string(line);
}

/* ----------------------------------------------------------------------
   return number of words
------------------------------------------------------------------------- */

size_t utils::count_words(const char * text) {
  size_t count = 0;
  const char * buf = text;
  char c = *buf;

  while (c) {
    if (c == ' ' || c == '\t' || c == '\r' ||  c == '\n' || c == '\f') {
      c = *++buf;
      continue;
    };

    ++count;
    c = *++buf;

    while (c) {
      if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f') {
        break;
      }
      c = *++buf;
    }
  }

  return count;
}

/* ----------------------------------------------------------------------
   return number of words
------------------------------------------------------------------------- */

size_t utils::count_words(const std::string & text) {
  return utils::count_words(text.c_str());
}

/* ----------------------------------------------------------------------
   Return number of words
------------------------------------------------------------------------- */

size_t utils::count_words(const std::string & text, const std::string & separators) {
  size_t count = 0;
  size_t start = text.find_first_not_of(separators);

  while (start != std::string::npos) {
    size_t end = text.find_first_of(separators, start);
    ++count;

    if(end == std::string::npos) {
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

size_t utils::trim_and_count_words(const std::string & text, const std::string & separators) {
  return utils::count_words(utils::trim_comment(text), separators);
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
    if (c == ' ' || c == '\t' || c == '\r' ||  c == '\n' || c == '\f') {
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
      while (((c != '\'') && (c != '\0'))
             || ((c == '\\') && (buf[1] == '\''))) {
        if ((c == '\\') && (buf[1] == '\'')) {
          ++buf;
          ++len;
        }
        c = *++buf;
        ++len;
      }
      if (c != '\'') ++len;
      c = *++buf;

      // handle double quote
    } else if (c == '"') {
      ++beg;
      add = 1;
      c = *++buf;
      while (((c != '"') && (c != '\0'))
             || ((c == '\\') && (buf[1] == '"'))) {
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
    while (1) {
      if ((c == '\'') || (c == '"')) goto quoted;
      // skip escaped quote
      if ((c == '\\') && ((buf[1] == '\'') || (buf[1] == '"'))) {
        ++buf;
        ++len;
        c = *++buf;
        ++len;
      }
      if ((c == ' ') || (c == '\t') || (c == '\r') || (c == '\n')
          || (c == '\f') || (c == '\0')) {
          list.push_back(text.substr(beg,len));
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
   Return whether string is a valid integer number
------------------------------------------------------------------------- */

bool utils::is_integer(const std::string & str) {
  if (str.size() == 0) {
    return false;
  }

  for (auto c : str) {
    if (isdigit(c) || c == '-' || c == '+') continue;
    return false;
  }
  return true;
}

/* ----------------------------------------------------------------------
   Return whether string is a valid floating-point number
------------------------------------------------------------------------- */

bool utils::is_double(const std::string & str) {
  if (str.size() == 0) {
    return false;
  }

  for (auto c : str) {
    if (isdigit(c)) continue;
    if (c == '-' || c == '+' || c == '.') continue;
    if (c == 'e' || c == 'E') continue;
    return false;
  }
  return true;
}

/* ----------------------------------------------------------------------
   strip off leading part of path, return just the filename
------------------------------------------------------------------------- */

std::string utils::path_basename(const std::string & path) {
#if defined(_WIN32)
  size_t start = path.find_last_of("/\\");
#else
  size_t start = path.find_last_of("/");
#endif

  if (start == std::string::npos) {
    start = 0;
  } else {
    start += 1;
  }

  return path.substr(start);
}

/* ----------------------------------------------------------------------
   join two paths
------------------------------------------------------------------------- */

std::string utils::path_join(const std::string & a, const std::string & b) {
  #if defined(_WIN32)
    return fmt::format("{}\\{}", a, b);
  #else
    return fmt::format("{}/{}", a, b);
  #endif
}

/* ----------------------------------------------------------------------
   try to open file for reading
------------------------------------------------------------------------- */

bool utils::file_is_readable(const std::string & path) {
  FILE * fp = fopen(path.c_str(), "r");
  if(fp) {
    fclose(fp);
    return true;
  }
  return false;
}

/* ----------------------------------------------------------------------
   try to find potential file as specified by name
   search current directory and the LAMMPS_POTENTIALS directory if
   specified
------------------------------------------------------------------------- */

std::string utils::get_potential_file_path(const std::string& path) {
  std::string filepath = path;
  std::string filename = utils::path_basename(path);

  if(utils::file_is_readable(filepath)) {
    return filepath;
  } else {
    // try the environment variable directory
    const char *path = getenv("LAMMPS_POTENTIALS");

    if (path != nullptr){
      std::string pot = utils::path_basename(filepath);
      filepath = utils::path_join(path, pot);

      if (utils::file_is_readable(filepath)) {
        return filepath;
      }
    }
  }
  return "";
}

/* ----------------------------------------------------------------------
   read first line of potential file
   if it has a DATE field, return the following word
------------------------------------------------------------------------- */

std::string utils::get_potential_date(const std::string & path, const std::string & potential_name) {
  TextFileReader reader(path, potential_name);
  reader.ignore_comments = false;

  char *line = reader.next_line();
  ValueTokenizer values(line);
  while (values.has_next()) {
    std::string word = values.next_string();
    if (word == "DATE:") {
      if (values.has_next()) {
        std::string date = values.next_string();
        return date;
      }
    }
  }
  return "";
}

/* ----------------------------------------------------------------------
   read first line of potential file
   if it has UNITS field, return following word
------------------------------------------------------------------------- */

std::string utils::get_potential_units(const std::string & path, const std::string & potential_name) {
  TextFileReader reader(path, potential_name);
  reader.ignore_comments = false;

  char *line = reader.next_line();
  ValueTokenizer values(line);
  while (values.has_next()) {
    std::string word = values.next_string();
    if (word == "UNITS:") {
      if (values.has_next()) {
        std::string units = values.next_string();
        return units;
      }
    }
  }
  return "";
}

/* ----------------------------------------------------------------------
   return bitmask of supported conversions for a given property
------------------------------------------------------------------------- */
int utils::get_supported_conversions(const int property)
{
  if (property == ENERGY) {
    return METAL2REAL | REAL2METAL;
  }
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
      return 1.0/23.060549;
    }
  }
  return 0.0;
}

/* ----------------------------------------------------------------------
   convert a timespec ([[HH:]MM:]SS) to seconds
   the strings "off" and "unlimited" result in -1.0;
------------------------------------------------------------------------- */

double utils::timespec2seconds(const std::string & timespec)
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
  } catch (TokenizerException & e) {
    return -1.0;
  }

  if (i == 3) return (vals[0]*60 + vals[1])*60 + vals[2];
  else if (i == 2) return vals[0]*60 + vals[1];
  return vals[0];
}

/* ------------------------------------------------------------------ */

extern "C" {
  /* Typedef'd pointer to get abstract datatype. */
  typedef struct regex_t *re_t;

  /* Compile regex string pattern to a regex_t-array. */
  static re_t re_compile(const char *pattern);


  /* Find matches of the compiled pattern inside text. */
  static int  re_matchp(const char *text, re_t pattern);


/* Definitions: */

#define MAX_REGEXP_OBJECTS 30 /* Max number of regex symbols in expression. */
#define MAX_CHAR_CLASS_LEN 40 /* Max length of character-class buffer in.   */


  enum { UNUSED, DOT, BEGIN, END, QUESTIONMARK, STAR, PLUS,
         CHAR, CHAR_CLASS, INV_CHAR_CLASS, DIGIT, NOT_DIGIT,
         INTEGER, NOT_INTEGER, FLOAT, NOT_FLOAT,
         ALPHA, NOT_ALPHA, WHITESPACE, NOT_WHITESPACE /*, BRANCH */ };

  typedef struct regex_t {
    unsigned char  type;   /* CHAR, STAR, etc.                      */
    union {
      unsigned char  ch;   /*      the character itself             */
      unsigned char *ccl;  /*  OR  a pointer to characters in class */
    };
  } regex_t;

/* Private function declarations: */
  static int matchpattern(regex_t *pattern, const char *text);
  static int matchcharclass(char c, const char *str);
  static int matchstar(regex_t p, regex_t *pattern, const char *text);
  static int matchplus(regex_t p, regex_t *pattern, const char *text);
  static int matchone(regex_t p, char c);
  static int matchdigit(char c);
  static int matchint(char c);
  static int matchfloat(char c);
  static int matchalpha(char c);
  static int matchwhitespace(char c);
  static int matchmetachar(char c, const char *str);
  static int matchrange(char c, const char *str);
  static int ismetachar(char c);

/* Semi-public functions: */
  int re_match(const char *text, const char *pattern)
  {
    return re_matchp(text, re_compile(pattern));
  }

  int re_matchp(const char *text, re_t pattern)
  {
    if (pattern != 0) {
      if (pattern[0].type == BEGIN) {
        return ((matchpattern(&pattern[1], text)) ? 0 : -1);
      } else {
        int idx = -1;

        do {
          idx += 1;

          if (matchpattern(pattern, text)) {
            if (text[0] == '\0')
              return -1;

            return idx;
          }
        }
        while (*text++ != '\0');
      }
    }
    return -1;
  }

  re_t re_compile(const char *pattern)
  {
    /* The sizes of the two static arrays below substantiates the static RAM usage of this module.
       MAX_REGEXP_OBJECTS is the max number of symbols in the expression.
       MAX_CHAR_CLASS_LEN determines the size of buffer for chars in all char-classes in the expression. */
    static regex_t re_compiled[MAX_REGEXP_OBJECTS];
    static unsigned char ccl_buf[MAX_CHAR_CLASS_LEN];
    int ccl_bufidx = 1;

    char c;     /* current char in pattern   */
    int i = 0;  /* index into pattern        */
    int j = 0;  /* index into re_compiled    */

    while (pattern[i] != '\0' && (j+1 < MAX_REGEXP_OBJECTS)) {
      c = pattern[i];

      switch (c) {
        /* Meta-characters: */
      case '^': {    re_compiled[j].type = BEGIN;           } break;
      case '$': {    re_compiled[j].type = END;             } break;
      case '.': {    re_compiled[j].type = DOT;             } break;
      case '*': {    re_compiled[j].type = STAR;            } break;
      case '+': {    re_compiled[j].type = PLUS;            } break;
      case '?': {    re_compiled[j].type = QUESTIONMARK;    } break;

        /* Escaped character-classes (\s \w ...): */
      case '\\': {
        if (pattern[i+1] != '\0') {
          /* Skip the escape-char '\\' */
          i += 1;
          /* ... and check the next */
          switch (pattern[i]) {
            /* Meta-character: */
          case 'd': {    re_compiled[j].type = DIGIT;            } break;
          case 'D': {    re_compiled[j].type = NOT_DIGIT;        } break;
          case 'i': {    re_compiled[j].type = INTEGER;          } break;
          case 'I': {    re_compiled[j].type = NOT_INTEGER;      } break;
          case 'f': {    re_compiled[j].type = FLOAT;            } break;
          case 'F': {    re_compiled[j].type = NOT_FLOAT;        } break;
          case 'w': {    re_compiled[j].type = ALPHA;            } break;
          case 'W': {    re_compiled[j].type = NOT_ALPHA;        } break;
          case 's': {    re_compiled[j].type = WHITESPACE;       } break;
          case 'S': {    re_compiled[j].type = NOT_WHITESPACE;   } break;

            /* Escaped character, e.g. '.' or '$' */
          default: {
            re_compiled[j].type = CHAR;
            re_compiled[j].ch = pattern[i];
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
        if (pattern[i+1] == '^') {
          re_compiled[j].type = INV_CHAR_CLASS;
          i += 1; /* Increment i to avoid including '^' in the char-buffer */
        } else {
          re_compiled[j].type = CHAR_CLASS;
        }

        /* Copy characters inside [..] to buffer */
        while ((pattern[++i] != ']') && (pattern[i] != '\0')) {
          /* Missing ] */
          if (pattern[i] == '\\') {
            if (ccl_bufidx >= MAX_CHAR_CLASS_LEN - 1) {
              return 0;
            }
            ccl_buf[ccl_bufidx++] = pattern[i++];
          } else if (ccl_bufidx >= MAX_CHAR_CLASS_LEN) {
            return 0;
          }
          ccl_buf[ccl_bufidx++] = pattern[i];
        }
        if (ccl_bufidx >= MAX_CHAR_CLASS_LEN) {
          /* Catches cases such as [00000000000000000000000000000000000000][ */
          return 0;
        }
        /* Null-terminate string end */
        ccl_buf[ccl_bufidx++] = 0;
        re_compiled[j].ccl = &ccl_buf[buf_begin];
      } break;

        /* Other characters: */
      default: {
        re_compiled[j].type = CHAR;
        re_compiled[j].ch = c;
      } break;
      }
      i += 1;
      j += 1;
    }
    /* 'UNUSED' is a sentinel used to indicate end-of-pattern */
    re_compiled[j].type = UNUSED;

    return (re_t) re_compiled;
  }


/* Private functions: */
  static int matchdigit(char c)
  {
    return ((c >= '0') && (c <= '9'));
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
    return ((c >= 'a') && (c <= 'z')) || ((c >= 'A') && (c <= 'Z'));
  }

  static int matchwhitespace(char c)
  {
    return ((c == ' ') || (c == '\t') || (c == '\n') || (c == '\r') || (c == '\f') || (c == '\v'));
  }

  static int matchalphanum(char c)
  {
    return ((c == '_') || matchalpha(c) || matchdigit(c));
  }

  static int matchrange(char c, const char *str)
  {
    return ((c != '-') && (str[0] != '\0')
            && (str[0] != '-') && (str[1] == '-')
            && (str[1] != '\0') && (str[2] != '\0')
            && ((c >= str[0]) && (c <= str[2])));
  }

  static int ismetachar(char c)
  {
    return ((c == 's') || (c == 'S')
            || (c == 'w') || (c == 'W')
            || (c == 'd') || (c == 'D'));
  }

  static int matchmetachar(char c, const char *str)
  {
    switch (str[0]) {
    case 'd': return  matchdigit(c);
    case 'D': return !matchdigit(c);
    case 'i': return  matchint(c);
    case 'I': return !matchint(c);
    case 'f': return  matchfloat(c);
    case 'F': return !matchfloat(c);
    case 'w': return  matchalphanum(c);
    case 'W': return !matchalphanum(c);
    case 's': return  matchwhitespace(c);
    case 'S': return !matchwhitespace(c);
    default:  return (c == str[0]);
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
    }
    while (*str++ != '\0');

    return 0;
  }

  static int matchone(regex_t p, char c)
  {
    switch (p.type) {
    case DOT:            return 1;
    case CHAR_CLASS:     return  matchcharclass(c, (const char *)p.ccl);
    case INV_CHAR_CLASS: return !matchcharclass(c, (const char *)p.ccl);
    case DIGIT:          return  matchdigit(c);
    case NOT_DIGIT:      return !matchdigit(c);
    case INTEGER:        return  matchint(c);
    case NOT_INTEGER:    return !matchint(c);
    case FLOAT:          return  matchfloat(c);
    case NOT_FLOAT:      return !matchfloat(c);
    case ALPHA:          return  matchalphanum(c);
    case NOT_ALPHA:      return !matchalphanum(c);
    case WHITESPACE:     return  matchwhitespace(c);
    case NOT_WHITESPACE: return !matchwhitespace(c);
    default:             return  (p.ch == c);
    }
  }

  static int matchstar(regex_t p, regex_t *pattern, const char *text)
  {
    do {
      if (matchpattern(pattern, text))
        return 1;
    }
    while ((text[0] != '\0') && matchone(p, *text++));

    return 0;
  }

  static int matchplus(regex_t p, regex_t *pattern, const char *text)
  {
    while ((text[0] != '\0') && matchone(p, *text++)) {
      if (matchpattern(pattern, text))
        return 1;
    }
    return 0;
  }

  static int matchquestion(regex_t p, regex_t *pattern, const char *text)
  {
    if (p.type == UNUSED)
      return 1;
    if (matchpattern(pattern, text))
      return 1;
    if (*text && matchone(p, *text++))
      return matchpattern(pattern, text);
    return 0;
  }

/* Iterative matching */
  static int matchpattern(regex_t *pattern, const char *text)
  {
    do {
      if ((pattern[0].type == UNUSED) || (pattern[1].type == QUESTIONMARK)) {
        return matchquestion(pattern[0], &pattern[2], text);
      } else if (pattern[1].type == STAR) {
        return matchstar(pattern[0], &pattern[2], text);
      } else if (pattern[1].type == PLUS) {
        return matchplus(pattern[0], &pattern[2], text);
      } else if ((pattern[0].type == END) && pattern[1].type == UNUSED) {
        return (text[0] == '\0');
      }
    }
    while ((text[0] != '\0') && matchone(*pattern++, *text++));

    return 0;
  }

}
