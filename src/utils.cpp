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

#include <cstring>
#include <cstdio>
#include <cstddef>
#include "utils.h"
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

/**
 * JSON type identifier. Basic types are:
 * 	o Object
 * 	o Array
 * 	o String
 * 	o Other primitive: number, boolean (true/false) or null
 */
  typedef enum {
	JSMN_UNDEFINED = 0,
	JSMN_OBJECT = 1,
	JSMN_ARRAY = 2,
	JSMN_STRING = 3,
	JSMN_PRIMITIVE = 4
  } jsmntype_t;

  enum jsmnerr {
	/* Not enough tokens were provided */
	JSMN_ERROR_NOMEM = -1,
	/* Invalid character inside JSON string */
	JSMN_ERROR_INVAL = -2,
	/* The string is not a full JSON packet, more bytes expected */
	JSMN_ERROR_PART = -3
  };

/**
 * JSON token description.
 * type		type (object, array, string etc.)
 * start	start position in JSON data string
 * end		end position in JSON data string
 */
  typedef struct {
	jsmntype_t type;
	int start;
	int end;
	int size;
  } jsmntok_t;

/**
 * JSON parser. Contains an array of token blocks available. Also stores
 * the string being parsed now and current position in that string
 */
  typedef struct {
	unsigned int pos; /* offset in the JSON string */
	unsigned int toknext; /* next token to allocate */
	int toksuper; /* superior token node, e.g parent object or array */
  } jsmn_parser;

/**
 * Create JSON parser over an array of tokens
 */
  static void jsmn_init(jsmn_parser *parser);

/**
 * Run JSON parser. It parses a JSON data string into and array of tokens, each describing
 * a single JSON object.
 */
  static int jsmn_parse(jsmn_parser *parser, const char *js, size_t len,
                        jsmntok_t *tokens, unsigned int num_tokens);
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



int utils::kim_simulator_json_parse(int argc, char **argv)
{
    FILE *fp;
    jsmn_parser p;
    jsmntok_t *tok;
    char *buf;
    size_t nbytes;

    if (argc != 2) {
        printf("usage: %s <json-file>\n",argv[0]);
        return 1;
    }

    // open JSON file

    fp = fopen(argv[1],"rb");
    if (!fp) {
        perror("Error opening JSON file");
        return 2;
    }

    // determine file size and allocate suitable buffer

    fseek(fp,0,SEEK_END);
    long int flen = ftell(fp);
    rewind(fp);
    buf = new char[flen];
    nbytes = fread(buf,1,flen,fp);
    fclose(fp);

    // parse once to count number of tokens

    jsmn_init(&p);
    int ntok = jsmn_parse(&p,buf,nbytes,NULL,1);
    if (ntok < 0) {
        printf("failed to parse JSON: %d\n",ntok);
        return 3;
    }

    // allocate token storage and parse again

    jsmn_init(&p);
    tok = new jsmntok_t[ntok];
    int retval = jsmn_parse(&p,buf,nbytes,tok,ntok);
    if ((retval < 1) || (tok[0].type != JSMN_OBJECT)) {
        printf("failed to parse JSON: no root object\n");
        return 4;
    }

    for (int i=1; i < retval; ++i) {
        printf("key: %.*s\n",tok[i].end-tok[i].start,buf+tok[i].start);
        if (tok[i+1].type == JSMN_ARRAY) {
            printf("value is array of size %d\n",tok[i+1].size);
            i += tok[i+1].size + 1;
        } else {
            ++i;
        }
    }
    
    delete [] buf;
    delete [] tok;
    return 0;
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

  
/**
 * Allocates a fresh unused token from the token pool.
 */
  static jsmntok_t *jsmn_alloc_token(jsmn_parser *parser,
                                     jsmntok_t *tokens, size_t num_tokens) {
	jsmntok_t *tok;
	if (parser->toknext >= num_tokens) {
      return NULL;
	}
	tok = &tokens[parser->toknext++];
	tok->start = tok->end = -1;
	tok->size = 0;
	return tok;
  }

/**
 * Fills token type and boundaries.
 */
  static void jsmn_fill_token(jsmntok_t *token, jsmntype_t type,
                              int start, int end) {
	token->type = type;
	token->start = start;
	token->end = end;
	token->size = 0;
  }

/**
 * Fills next available token with JSON primitive.
 */
  static int jsmn_parse_primitive(jsmn_parser *parser, const char *js,
                                  size_t len, jsmntok_t *tokens, size_t num_tokens) {
	jsmntok_t *token;
	int start;

	start = parser->pos;

	for (; parser->pos < len && js[parser->pos] != '\0'; parser->pos++) {
      switch (js[parser->pos]) {
      case '\t' : case '\r' : case '\n' : case ' ' :
      case ','  : case ']'  : case '}' :
        goto found;
      }
      if (js[parser->pos] < 32 || js[parser->pos] >= 127) {
        parser->pos = start;
        return JSMN_ERROR_INVAL;
      }
	}

    found:
	if (tokens == NULL) {
      parser->pos--;
      return 0;
	}
	token = jsmn_alloc_token(parser, tokens, num_tokens);
	if (token == NULL) {
      parser->pos = start;
      return JSMN_ERROR_NOMEM;
	}
	jsmn_fill_token(token, JSMN_PRIMITIVE, start, parser->pos);
	parser->pos--;
	return 0;
  }

/**
 * Fills next token with JSON string.
 */
  static int jsmn_parse_string(jsmn_parser *parser, const char *js,
                               size_t len, jsmntok_t *tokens, size_t num_tokens) {
	jsmntok_t *token;

	int start = parser->pos;

	parser->pos++;

	/* Skip starting quote */
	for (; parser->pos < len && js[parser->pos] != '\0'; parser->pos++) {
      char c = js[parser->pos];

      /* Quote: end of string */
      if (c == '\"') {
        if (tokens == NULL) {
          return 0;
        }
        token = jsmn_alloc_token(parser, tokens, num_tokens);
        if (token == NULL) {
          parser->pos = start;
          return JSMN_ERROR_NOMEM;
        }
        jsmn_fill_token(token, JSMN_STRING, start+1, parser->pos);
        return 0;
      }

      /* Backslash: Quoted symbol expected */
      if (c == '\\' && parser->pos + 1 < len) {
        int i;
        parser->pos++;
        switch (js[parser->pos]) {
          /* Allowed escaped symbols */
        case '\"': case '/' : case '\\' : case 'b' :
        case 'f' : case 'r' : case 'n'  : case 't' :
          break;
          /* Allows escaped symbol \uXXXX */
        case 'u':
          parser->pos++;
          for(i = 0; i < 4 && parser->pos < len && js[parser->pos] != '\0'; i++) {
            /* If it isn't a hex character we have an error */
            if(!((js[parser->pos] >= 48 && js[parser->pos] <= 57) || /* 0-9 */
                 (js[parser->pos] >= 65 && js[parser->pos] <= 70) || /* A-F */
                 (js[parser->pos] >= 97 && js[parser->pos] <= 102))) { /* a-f */
              parser->pos = start;
              return JSMN_ERROR_INVAL;
            }
            parser->pos++;
          }
          parser->pos--;
          break;
          /* Unexpected symbol */
        default:
          parser->pos = start;
          return JSMN_ERROR_INVAL;
        }
      }
	}
	parser->pos = start;
	return JSMN_ERROR_PART;
  }

/**
 * Parse JSON string and fill tokens.
 */
  int jsmn_parse(jsmn_parser *parser, const char *js, size_t len,
                 jsmntok_t *tokens, unsigned int num_tokens) {
	int r;
	int i;
	jsmntok_t *token;
	int count = parser->toknext;

	for (; parser->pos < len && js[parser->pos] != '\0'; parser->pos++) {
      char c;
      jsmntype_t type;

      c = js[parser->pos];
      switch (c) {
      case '{': case '[':
        count++;
        if (tokens == NULL) {
          break;
        }
        token = jsmn_alloc_token(parser, tokens, num_tokens);
        if (token == NULL)
          return JSMN_ERROR_NOMEM;
        if (parser->toksuper != -1) {
          tokens[parser->toksuper].size++;
        }
        token->type = (c == '{' ? JSMN_OBJECT : JSMN_ARRAY);
        token->start = parser->pos;
        parser->toksuper = parser->toknext - 1;
        break;
      case '}': case ']':
        if (tokens == NULL)
          break;
        type = (c == '}' ? JSMN_OBJECT : JSMN_ARRAY);
        for (i = parser->toknext - 1; i >= 0; i--) {
          token = &tokens[i];
          if (token->start != -1 && token->end == -1) {
            if (token->type != type) {
              return JSMN_ERROR_INVAL;
            }
            parser->toksuper = -1;
            token->end = parser->pos + 1;
            break;
          }
        }
        /* Error if unmatched closing bracket */
        if (i == -1) return JSMN_ERROR_INVAL;
        for (; i >= 0; i--) {
          token = &tokens[i];
          if (token->start != -1 && token->end == -1) {
            parser->toksuper = i;
            break;
          }
        }
        break;
      case '\"':
        r = jsmn_parse_string(parser, js, len, tokens, num_tokens);
        if (r < 0) return r;
        count++;
        if (parser->toksuper != -1 && tokens != NULL)
          tokens[parser->toksuper].size++;
        break;
      case '\t' : case '\r' : case '\n' : case ' ':
        break;
      case ':':
        parser->toksuper = parser->toknext - 1;
        break;
      case ',':
        if (tokens != NULL && parser->toksuper != -1 &&
            tokens[parser->toksuper].type != JSMN_ARRAY &&
            tokens[parser->toksuper].type != JSMN_OBJECT) {
          for (i = parser->toknext - 1; i >= 0; i--) {
            if (tokens[i].type == JSMN_ARRAY || tokens[i].type == JSMN_OBJECT) {
              if (tokens[i].start != -1 && tokens[i].end == -1) {
                parser->toksuper = i;
                break;
              }
            }
          }
        }
        break;
        /* In non-strict mode every unquoted value is a primitive */
      default:
        r = jsmn_parse_primitive(parser, js, len, tokens, num_tokens);
        if (r < 0) return r;
        count++;
        if (parser->toksuper != -1 && tokens != NULL)
          tokens[parser->toksuper].size++;
        break;

      }
	}

	if (tokens != NULL) {
      for (i = parser->toknext - 1; i >= 0; i--) {
        /* Unmatched opened object or array */
        if (tokens[i].start != -1 && tokens[i].end == -1) {
          return JSMN_ERROR_PART;
        }
      }
	}

	return count;
  }

/**
 * Creates a new parser based over a given  buffer with an array of tokens
 * available.
 */
  void jsmn_init(jsmn_parser *parser) {
	parser->pos = 0;
	parser->toknext = 0;
	parser->toksuper = -1;
  }
}

