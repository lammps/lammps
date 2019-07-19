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

#include <mpi.h>
#include <cstring>
#include "utils.h"
#include "lammps.h"
#include "error.h"

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
bool utils::strmatch(std::string text, std::string pattern)
{
  const int pos = re_match(text.c_str(),pattern.c_str());
  return (pos >= 0);
}

/* utility function to avoid code repetition when parsing args */
int utils::cfvarg(std::string mode, const char *arg, char *&cfv_id)
{
  int rv = utils::NONE;
  cfv_id = NULL;

  if (!arg) return rv;

  if (utils::strmatch(arg,std::string("^[") + mode + "]_")) {
    if (*arg == 'c') rv = utils::COMPUTE;
    else if (*arg == 'f') rv = utils::FIX;
    else if (*arg == 'v') rv = utils::VARIABLE;
    else return rv;             // should not happen

    arg += 2;
    int n = strlen(arg)+1;
    cfv_id = new char[n];
    strcpy(cfv_id,arg);
  }

  return rv;
}

/* like fgets() but aborts with an error or EOF is encountered */
void utils::sfgets(const char *srcname, int srcline, char *s, int size,
                   FILE *fp, const char *filename, Error *error)
{
  char *rv = fgets(s,size,fp);
  if (rv == NULL) { // something went wrong
    std::string errmsg;

    if (feof(fp)) {
      errmsg = "Unexpected end of file while reading file '";
    } else if (ferror(fp)) {
      errmsg = "Unexpected error while reading file '";
    } else {
      errmsg = "Unexpected short read while reading file '";
    }
    errmsg += filename;
    errmsg += "'";

    if (error) error->one(srcname,srcline,errmsg.c_str());
    if (s) *s = '\0'; // truncate string to empty in case error is NULL
  }
  return;
}

/* ------------------------------------------------------------------ */

std::string utils::check_packages_for_style(std::string style,
                                            std::string name, LAMMPS *lmp)
{
  std::string errmsg = "Unrecognized " + style + " style '" + name + "'";
  const char *pkg = lmp->match_style(style.c_str(),name.c_str());

  if (pkg) {
    errmsg += " is part of the " + std::string(pkg) + " package";
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
      lmp->error->one(file,line,msg.c_str());
    else
      lmp->error->all(file,line,msg.c_str());
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
      lmp->error->one(file,line,msg.c_str());
    else
      lmp->error->all(file,line,msg.c_str());
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
      lmp->error->one(file,line,msg.c_str());
    else
      lmp->error->all(file,line,msg.c_str());
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
      lmp->error->one(file,line,msg.c_str());
    else
      lmp->error->all(file,line,msg.c_str());
  }

  return ATOTAGINT(str);
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
