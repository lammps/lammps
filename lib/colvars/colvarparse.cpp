/// -*- c++ -*-

#include <sstream>
#include <iostream>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"


// space & tab
std::string const colvarparse::white_space = " \t";

std::string colvarparse::dummy_string = "";
size_t      colvarparse::dummy_pos = 0;


// definition of single-value keyword parsers

#define _get_keyval_scalar_string_(TYPE)                                \
                                                                        \
  bool colvarparse::get_keyval(std::string const &conf,                 \
                               char const *key,                         \
                               TYPE &value,                             \
                               TYPE const &def_value,                   \
                               Parse_Mode const parse_mode)             \
  {                                                                     \
    std::string data;                                                   \
    bool b_found = false, b_found_any = false;                          \
    size_t save_pos = 0, found_count = 0;                               \
                                                                        \
    do {                                                                \
      std::string data_this = "";                                       \
      b_found = key_lookup(conf, key, data_this, save_pos);             \
      if (b_found) {                                                    \
        if (!b_found_any)                                               \
          b_found_any = true;                                           \
        found_count++;                                                  \
        data = data_this;                                               \
      }                                                                 \
    } while (b_found);                                                  \
                                                                        \
    if (found_count > 1)                                                \
      cvm::log("Warning: found more than one instance of \""+           \
               std::string(key)+"\".\n");                               \
                                                                        \
    if (data.size()) {                                                  \
      std::istringstream is(data);                                      \
      TYPE x(def_value);                                                \
      if (is >> x)                                                      \
        value = x;                                                      \
      else                                                              \
        cvm::error("Error: in parsing \""+                              \
                   std::string(key)+"\".\n", INPUT_ERROR);              \
      if (parse_mode != parse_silent) {                                 \
        cvm::log("# "+std::string(key)+" = "+                           \
                 cvm::to_str(value)+"\n");                              \
      }                                                                 \
    } else {                                                            \
                                                                        \
      if (b_found_any)                                                  \
        cvm::error("Error: improper or missing value "                  \
                   "for \""+std::string(key)+"\".\n", INPUT_ERROR);     \
      value = def_value;                                                \
      if (parse_mode != parse_silent) {                                 \
        cvm::log("# "+std::string(key)+" = \""+                         \
                 cvm::to_str(def_value)+"\" [default]\n");              \
      }                                                                 \
    }                                                                   \
                                                                        \
    return b_found_any;                                                 \
  }


#define _get_keyval_scalar_(TYPE)                                       \
                                                                        \
  bool colvarparse::get_keyval(std::string const &conf,                 \
                               char const *key,                         \
                               TYPE &value,                             \
                               TYPE const &def_value,                   \
                               Parse_Mode const parse_mode)             \
  {                                                                     \
    std::string data;                                                   \
    bool b_found = false, b_found_any = false;                          \
    size_t save_pos = 0, found_count = 0;                               \
                                                                        \
    do {                                                                \
      std::string data_this = "";                                       \
      b_found = key_lookup(conf, key, data_this, save_pos);             \
      if (b_found) {                                                    \
        if (!b_found_any)                                               \
          b_found_any = true;                                           \
        found_count++;                                                  \
        data = data_this;                                               \
      }                                                                 \
    } while (b_found);                                                  \
                                                                        \
    if (found_count > 1)                                                \
      cvm::log("Warning: found more than one instance of \""+           \
               std::string(key)+"\".\n");                               \
                                                                        \
    if (data.size()) {                                                  \
      std::istringstream is(data);                                      \
      size_t data_count = 0;                                            \
      TYPE x(def_value);                                                \
      while (is >> x) {                                                 \
        value = x;                                                      \
        data_count++;                                                   \
      }                                                                 \
      if (data_count == 0)                                              \
        cvm::fatal_error("Error: in parsing \""+                        \
                         std::string(key)+"\".\n");                     \
      if (data_count > 1)                                               \
        cvm::error("Error: multiple values "                            \
                   "are not allowed for keyword \""+                    \
                   std::string(key)+"\".\n", INPUT_ERROR);              \
      if (parse_mode != parse_silent) {                                 \
        cvm::log("# "+std::string(key)+" = "+                           \
                 cvm::to_str(value)+"\n");                              \
      }                                                                 \
    } else {                                                            \
                                                                        \
      if (b_found_any)                                                  \
        cvm::error("Error: improper or missing value "                  \
                   "for \""+std::string(key)+"\".\n", INPUT_ERROR);     \
      value = def_value;                                                \
      if (parse_mode != parse_silent) {                                 \
        cvm::log("# "+std::string(key)+" = "+                           \
                 cvm::to_str(def_value)+" [default]\n");                \
      }                                                                 \
    }                                                                   \
                                                                        \
    return b_found_any;                                                 \
  }


// definition of multiple-value keyword parsers

#define _get_keyval_vector_(TYPE)                                       \
                                                                        \
  bool colvarparse::get_keyval(std::string const &conf,                 \
                               char const *key,                         \
                               std::vector<TYPE> &values,               \
                               std::vector<TYPE> const &def_values,     \
                               Parse_Mode const parse_mode)             \
  {                                                                     \
    std::string data;                                                   \
    bool b_found = false, b_found_any = false;                          \
    size_t save_pos = 0, found_count = 0;                               \
                                                                        \
    do {                                                                \
      std::string data_this = "";                                       \
      b_found = key_lookup(conf, key, data_this, save_pos);             \
      if (b_found) {                                                    \
        if (!b_found_any)                                               \
          b_found_any = true;                                           \
        found_count++;                                                  \
        data = data_this;                                               \
      }                                                                 \
    } while (b_found);                                                  \
                                                                        \
    if (found_count > 1)                                                \
      cvm::log("Warning: found more than one instance of \""+           \
               std::string(key)+"\".\n");                               \
                                                                        \
    if (data.size()) {                                                  \
      std::istringstream is(data);                                      \
                                                                        \
      if (values.size() == 0) {                                         \
                                                                        \
        std::vector<TYPE> x;                                            \
        if (def_values.size())                                          \
          x = def_values;                                               \
        else                                                            \
          x.assign(1, TYPE());                                          \
                                                                        \
        for (size_t i = 0;                                              \
             ( is >> x[ ((i<x.size()) ? i : x.size()-1) ] );            \
             i++) {                                                     \
          values.push_back(x[ ((i<x.size()) ? i : x.size()-1) ]);       \
        }                                                               \
                                                                        \
      } else {                                                          \
                                                                        \
        size_t i = 0;                                                   \
        for ( ; i < values.size(); i++) {                               \
          TYPE x(values[i]);                                            \
          if (is >> x)                                                  \
            values[i] = x;                                              \
          else                                                          \
            cvm::error("Error: in parsing \""+                          \
                       std::string(key)+"\".\n", INPUT_ERROR);          \
        }                                                               \
      }                                                                 \
                                                                        \
      if (parse_mode != parse_silent) {                                 \
        cvm::log("# "+std::string(key)+" = "+                           \
                 cvm::to_str(values)+"\n");                             \
      }                                                                 \
                                                                        \
    } else {                                                            \
                                                                        \
      if (b_found_any)                                                  \
        cvm::error("Error: improper or missing values for \""+          \
                   std::string(key)+"\".\n", INPUT_ERROR);              \
                                                                        \
      for (size_t i = 0; i < values.size(); i++)                        \
        values[i] = def_values[ (i > def_values.size()) ? 0 : i ];      \
                                                                        \
      if (parse_mode != parse_silent) {                                 \
        cvm::log("# "+std::string(key)+" = "+                           \
                 cvm::to_str(def_values)+" [default]\n");               \
      }                                                                 \
    }                                                                   \
                                                                        \
    return b_found_any;                                                 \
  }


// single-value keyword parsers

_get_keyval_scalar_(int);
_get_keyval_scalar_(size_t);
_get_keyval_scalar_string_(std::string);
_get_keyval_scalar_(cvm::real);
_get_keyval_scalar_(cvm::rvector);
_get_keyval_scalar_(cvm::quaternion);
_get_keyval_scalar_(colvarvalue);


// multiple-value keyword parsers

_get_keyval_vector_(int);
_get_keyval_vector_(size_t);
_get_keyval_vector_(std::string);
_get_keyval_vector_(cvm::real);
_get_keyval_vector_(cvm::rvector);
_get_keyval_vector_(cvm::quaternion);
_get_keyval_vector_(colvarvalue);



bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             bool &value,
                             bool const &def_value,
                             Parse_Mode const parse_mode)
{
  std::string data;
  bool b_found = false, b_found_any = false;
  size_t save_pos = 0, found_count = 0;

  do {
    std::string data_this = "";
    b_found = key_lookup(conf, key, data_this, save_pos);
    if (b_found) {
      if (!b_found_any)
        b_found_any = true;
      found_count++;
      data = data_this;
    }
  } while (b_found);

  if (found_count > 1)
    cvm::log("Warning: found more than one instance of \""+
             std::string(key)+"\".\n");

  if (data.size()) {
    if ( (data == std::string("on")) ||
         (data == std::string("yes")) ||
         (data == std::string("true")) ) {
      value = true;
    } else if ( (data == std::string("off")) ||
                (data == std::string("no")) ||
                (data == std::string("false")) ) {
      value = false;
    } else
      cvm::fatal_error("Error: boolean values only are allowed "
                       "for \""+std::string(key)+"\".\n");
    if (parse_mode != parse_silent) {
      cvm::log("# "+std::string(key)+" = "+
               (value ? "on" : "off")+"\n");
    }
  } else {

    if (b_found_any) {
      if (parse_mode != parse_silent) {
        cvm::log("# "+std::string(key)+" = on\n");
      }
      value = true;
    } else {
      value = def_value;
      if (parse_mode != parse_silent) {
        cvm::log("# "+std::string(key)+" = "+
                 (def_value ? "on" : "off")+" [default]\n");
      }
    }
  }

  return b_found_any;
}


void colvarparse::add_keyword(char const *key)
{
  for (std::list<std::string>::iterator ki = allowed_keywords.begin();
       ki != allowed_keywords.end(); ki++) {
    if (to_lower_cppstr(std::string(key)) == *ki)
      return;
  }
  // not found in the list
  //   if (cvm::debug())
  //     cvm::log("Registering a new keyword, \""+std::string (key)+"\".\n");
  allowed_keywords.push_back(to_lower_cppstr(std::string(key)));
}


void colvarparse::strip_values(std::string &conf)
{
  size_t offset = 0;
  data_begin_pos.sort();
  data_end_pos.sort();

  std::list<size_t>::iterator data_begin = data_begin_pos.begin();
  std::list<size_t>::iterator data_end   = data_end_pos.begin();

  for ( ; (data_begin != data_begin_pos.end()) &&
          (data_end   != data_end_pos.end()) ;
        data_begin++, data_end++) {

    //     std::cerr << "data_begin, data_end "
    //               << *data_begin << ", " << *data_end
    //               << "\n";

    size_t const nchars = *data_end-*data_begin;

    //     std::cerr << "conf[data_begin:data_end] = \""
    //               << std::string (conf, *data_begin - offset, nchars)
    //               << "\"\n";

    conf.erase(*data_begin - offset, nchars);
    offset += nchars;

    //     std::cerr << ("Stripped config = \"\n"+conf+"\"\n");

  }
}


int colvarparse::check_keywords(std::string &conf, char const *key)
{
  if (cvm::debug())
    cvm::log("Configuration string for \""+std::string(key)+
             "\": \"\n"+conf+"\".\n");

  strip_values(conf);
  // after stripping, the config string has either empty lines, or
  // lines beginning with a keyword

  std::string line;
  std::istringstream is(conf);
  while (std::getline(is, line)) {
    if (line.size() == 0)
      continue;
    if (line.find_first_not_of(white_space) ==
        std::string::npos)
      continue;

    std::string uk;
    std::istringstream line_is(line);
    line_is >> uk;
    // if (cvm::debug())
    //   cvm::log ("Checking the validity of \""+uk+"\" from line:\n" + line);
    uk = to_lower_cppstr(uk);

    bool found_keyword = false;
    for (std::list<std::string>::iterator ki = allowed_keywords.begin();
         ki != allowed_keywords.end(); ki++) {
      if (uk == *ki) {
        found_keyword = true;
        break;
      }
    }
    if (!found_keyword) {
      cvm::log("Error: keyword \""+uk+"\" is not supported, "
               "or not recognized in this context.\n");
      cvm::set_error_bits(INPUT_ERROR);
      return COLVARS_ERROR;
    }
  }
  allowed_keywords.clear();
  data_begin_pos.clear();
  data_end_pos.clear();
  return COLVARS_OK;
}


std::istream & colvarparse::getline_nocomments(std::istream &is,
                                               std::string &line,
                                               char const delim)
{
  std::getline(is, line, delim);
  size_t const comment = line.find('#');
  if (comment != std::string::npos) {
    line.erase(comment);
  }
  return is;
}


bool colvarparse::key_lookup(std::string const &conf,
                             char const *key_in,
                             std::string &data,
                             size_t &save_pos)
{
  // add this keyword to the register (in its camelCase version)
  add_keyword(key_in);

  // use the lowercase version from now on
  std::string const key(to_lower_cppstr(key_in));

  // "conf_lower" is only used to lookup the keyword, but its value
  // will be read from "conf", in order not to mess up file names
  std::string const conf_lower(to_lower_cppstr(conf));

  // by default, there is no value, unless we found one
  data = "";

  // when the function is invoked without save_pos, ensure that we
  // start from zero
  colvarparse::dummy_pos = 0;

  // start from the first occurrence of key
  size_t pos = conf_lower.find(key, save_pos);

  // iterate over all instances until it finds the isolated keyword
  while (true) {

    if (pos == std::string::npos) {
      // no valid instance of the keyword has been found
      return false;
    }

    bool b_isolated_left = true, b_isolated_right = true;

    if (pos > 0) {
      if ( std::string("\n"+white_space+
                       "}").find(conf[pos-1]) ==
           std::string::npos ) {
        // none of the valid delimiting characters is on the left of key
        b_isolated_left = false;
      }
    }

    if (pos < conf.size()-key.size()-1) {
      if ( std::string("\n"+white_space+
                       "{").find(conf[pos+key.size()]) ==
           std::string::npos ) {
        // none of the valid delimiting characters is on the right of key
        b_isolated_right = false;
      }
    }

    // check that there are matching braces between here and the end of conf
    bool const b_not_within_block = brace_check(conf, pos);

    bool const b_isolated = (b_isolated_left && b_isolated_right &&
                             b_not_within_block);

    if (b_isolated) {
      // found it
      break;
    } else {
      // try the next occurrence of key
      pos = conf_lower.find(key, pos+key.size());
    }
  }

  // check it is not between quotes
  //   if ( (conf.find_last_of  ("\"",
  //                             conf.find_last_of  (white_space, pos)) !=
  //         std::string::npos) &&
  //        (conf.find_first_of ("\"",
  //                             conf.find_first_of (white_space, pos)) !=
  //         std::string::npos) )
  //     return false;


  // save the pointer for a future call (when iterating over multiple
  // valid instances of the same keyword)
  save_pos = pos + key.size();

  // get the remainder of the line
  size_t pl = conf.rfind("\n", pos);
  size_t line_begin = (pl == std::string::npos) ? 0 : pos;
  size_t nl = conf.find  ("\n", pos);
  size_t line_end = (nl == std::string::npos) ? conf.size() : nl;
  std::string line(conf, line_begin, (line_end-line_begin));

  size_t data_begin = (to_lower_cppstr(line)).find(key) + key.size();
  data_begin = line.find_first_not_of(white_space, data_begin+1);

  //   size_t data_begin_absolute = data_begin + line_begin;
  //   size_t data_end_absolute   = data_begin;

  if (data_begin != std::string::npos) {

    size_t data_end = line.find_last_not_of(white_space) + 1;
    data_end = (data_end == std::string::npos) ? line.size() : data_end;
    //     data_end_absolute = data_end + line_begin;

    if (line.find('{', data_begin) != std::string::npos) {

      size_t brace_count = 1;
      size_t brace = line.find('{', data_begin);  // start from the first opening brace

      while (brace_count > 0) {

        // find the matching closing brace
        brace = line.find_first_of("{}", brace+1);
        while (brace == std::string::npos) {
          // add a new line
          if (line_end >= conf.size()) {
            cvm::fatal_error("Parse error: reached the end while "
                             "looking for closing brace; until now "
                             "the following was parsed: \"\n"+
                             line+"\".\n");
            return false;
          }
          size_t const old_end = line.size();
          //           data_end_absolute += old_end+1;

          line_begin = line_end;
          nl = conf.find('\n', line_begin+1);
          if (nl == std::string::npos)
            line_end = conf.size();
          else
            line_end = nl;
          line.append(conf, line_begin, (line_end-line_begin));

          brace = line.find_first_of("{}", old_end);
        }

        if (line[brace] == '{') brace_count++;
        if (line[brace] == '}') brace_count--;
      }

      // set data_begin after the opening brace
      data_begin = line.find_first_of('{', data_begin) + 1;
      data_begin = line.find_first_not_of(white_space,
                                          data_begin);
      // set data_end before the closing brace
      data_end = brace;
      data_end = line.find_last_not_of(white_space+"}",
                                       data_end) + 1;
      //       data_end_absolute = line_end;

      if (data_end > line.size())
        data_end = line.size();
    }

    data.append(line, data_begin, (data_end-data_begin));

    if (data.size() && save_delimiters) {
      data_begin_pos.push_back(conf.find(data, pos+key.size()));
      data_end_pos.push_back   (data_begin_pos.back()+data.size());
      //       std::cerr << "key = " << key << ", data = \""
      //                 << data << "\", data_begin, data_end = "
      //                 << data_begin_pos.back() << ", " << data_end_pos.back()
      //                 << "\n";
    }
  }

  save_pos = line_end;

  return true;
}


std::istream & operator>> (std::istream &is, colvarparse::read_block const &rb)
{
  size_t start_pos = is.tellg();
  std::string read_key, next;

  if ( !(is >> read_key) || !(read_key == rb.key) ||
       !(is >> next) ) {
    // the requested keyword has not been found, or it is not possible
    // to read data after it
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
    return is;
  }

  if (next != "{") {
    (*rb.data) = next;
    return is;
  }

  size_t brace_count = 1;
  std::string line;
  while (colvarparse::getline_nocomments(is, line)) {
    size_t br = 0, br_old = 0;
    while ( (br = line.find_first_of("{}", br)) != std::string::npos) {
      if (line[br] == '{') brace_count++;
      if (line[br] == '}') brace_count--;
      br_old = br;
      br++;
    }
    if (brace_count) (*rb.data).append(line + "\n");
    else {
      (*rb.data).append(line, 0, br_old);
      break;
    }
  }
  if (brace_count)  {
    // end-of-file reached
    // restore initial position
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
  }
  return is;
}


bool colvarparse::brace_check(std::string const &conf,
                              size_t const start_pos)
{
  size_t brace_count = 0;
  size_t brace = start_pos;
  while ( (brace = conf.find_first_of("{}", brace)) != std::string::npos) {
    if (conf[brace] == '{') brace_count++;
    if (conf[brace] == '}') brace_count--;
    brace++;
  }

  if (brace_count != 0)
    return false;
  else
    return true;
}
