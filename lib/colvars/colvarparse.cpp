// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.


#include <sstream>
#include <iostream>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"


// space & tab
char const * const colvarparse::white_space = " \t";


// definition of single-value keyword parsers

template<typename TYPE> bool colvarparse::_get_keyval_scalar_(std::string const &conf,
                                                              char const *key,
                                                              TYPE &value,
                                                              TYPE const &def_value,
                                                              Parse_Mode const parse_mode)
{
  std::string data;
  bool b_found = false, b_found_any = false;
  size_t save_pos = 0, found_count = 0;

  do {
    std::string data_this = "";
    b_found = key_lookup(conf, key, &data_this, &save_pos);
    if (b_found) {
      if (!b_found_any)
        b_found_any = true;
      found_count++;
      data = data_this;
    }
  } while (b_found);

  if (found_count > 1) {
    cvm::error("Error: found more than one instance of \""+
               std::string(key)+"\".\n", INPUT_ERROR);
  }

  if (data.size()) {
    std::istringstream is(data);
    TYPE x(def_value);
    if (is >> x) {
      value = x;
    } else {
      cvm::error("Error: in parsing \""+
                 std::string(key)+"\".\n", INPUT_ERROR);
    }
    if (parse_mode != parse_silent) {
      cvm::log("# "+std::string(key)+" = "+
               cvm::to_str(value)+"\n");
    }
  } else {

    if (b_found_any) {
      cvm::error("Error: improper or missing value "
                 "for \""+std::string(key)+"\".\n", INPUT_ERROR);
    }
    value = def_value;
    if (parse_mode != parse_silent) {
      cvm::log("# "+std::string(key)+" = "+
               cvm::to_str(def_value)+" [default]\n");
    }
  }

  return b_found_any;
}


bool colvarparse::_get_keyval_scalar_string_(std::string const &conf,
                                             char const *key,
                                             std::string &value,
                                             std::string const &def_value,
                                             Parse_Mode const parse_mode)
{
  std::string data;
  bool b_found = false, b_found_any = false;
  size_t save_pos = 0, found_count = 0;

  do {
    std::string data_this = "";
    b_found = key_lookup(conf, key, &data_this, &save_pos);
    if (b_found) {
      if (!b_found_any)
        b_found_any = true;
      found_count++;
      data = data_this;
    }
  } while (b_found);

  if (found_count > 1) {
    cvm::error("Error: found more than one instance of \""+
               std::string(key)+"\".\n", INPUT_ERROR);
  }

  if (data.size()) {
    std::istringstream is(data);
    size_t data_count = 0;
    std::string x(def_value);
    while (is >> x) {
      value = x;
      data_count++;
    }
    if (data_count == 0)
      cvm::error("Error: in parsing \""+
                 std::string(key)+"\".\n", INPUT_ERROR);
    if (data_count > 1) {
      cvm::error("Error: multiple values "
                 "are not allowed for keyword \""+
                 std::string(key)+"\".\n", INPUT_ERROR);
    }
    if (parse_mode != parse_silent) {
      cvm::log("# "+std::string(key)+" = \""+
               cvm::to_str(value)+"\"\n");
    }
  } else {

    if (b_found_any) {
      cvm::error("Error: improper or missing value "
                 "for \""+std::string(key)+"\".\n", INPUT_ERROR);
    }
    value = def_value;
    if (parse_mode != parse_silent) {
      cvm::log("# "+std::string(key)+" = \""+
               cvm::to_str(def_value)+"\" [default]\n");
    }
  }

  return b_found_any;
}


// multiple-value keyword parsers

template<typename TYPE> bool colvarparse::_get_keyval_vector_(std::string const &conf,
                                                              char const *key,
                                                              std::vector<TYPE> &values,
                                                              std::vector<TYPE> const &def_values,
                                                              Parse_Mode const parse_mode)
{
  std::string data;
  bool b_found = false, b_found_any = false;
  size_t save_pos = 0, found_count = 0;

  do {
    std::string data_this = "";
    b_found = key_lookup(conf, key, &data_this, &save_pos);
    if (b_found) {
      if (!b_found_any)
        b_found_any = true;
      found_count++;
      data = data_this;
    }
  } while (b_found);

  if (found_count > 1) {
    cvm::error("Error: found more than one instance of \""+
               std::string(key)+"\".\n", INPUT_ERROR);
  }

  if (data.size()) {
    std::istringstream is(data);

    if (values.size() == 0) {

      std::vector<TYPE> x;
      if (def_values.size())
        x = def_values;
      else
        x.assign(1, TYPE());

      for (size_t i = 0;
           ( is >> x[ ((i<x.size()) ? i : x.size()-1) ] );
           i++) {
        values.push_back(x[ ((i<x.size()) ? i : x.size()-1) ]);
      }

    } else {

      size_t i = 0;
      for ( ; i < values.size(); i++) {
        TYPE x(values[i]);
        if (is >> x) {
          values[i] = x;
        } else {
          cvm::error("Error: in parsing \""+
                     std::string(key)+"\".\n", INPUT_ERROR);
        }
      }
    }

    if (parse_mode != parse_silent) {
      cvm::log("# "+std::string(key)+" = "+
               cvm::to_str(values)+"\n");
    }

  } else {

    if (b_found_any) {
      cvm::error("Error: improper or missing values for \""+
                 std::string(key)+"\".\n", INPUT_ERROR);
    }

    for (size_t i = 0; i < values.size(); i++)
      values[i] = def_values[ (i > def_values.size()) ? 0 : i ];

    if (parse_mode != parse_silent) {
      cvm::log("# "+std::string(key)+" = "+
               cvm::to_str(def_values)+" [default]\n");
    }
  }

  return b_found_any;
}


// single-value keyword parsers


bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             int &value,
                             int const &def_value,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_scalar_<int>(conf, key, value, def_value, parse_mode);
}

bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             size_t &value,
                             size_t const &def_value,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_scalar_<size_t>(conf, key, value, def_value, parse_mode);
}

bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             long &value,
                             long const &def_value,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_scalar_<long>(conf, key, value, def_value, parse_mode);
}

bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             std::string &value,
                             std::string const &def_value,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_scalar_string_(conf, key, value, def_value, parse_mode);
}

bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             cvm::real &value,
                             cvm::real const &def_value,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_scalar_<cvm::real>(conf, key, value, def_value, parse_mode);
}

bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             cvm::rvector &value,
                             cvm::rvector const &def_value,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_scalar_<cvm::rvector>(conf, key, value, def_value, parse_mode);
}

bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             cvm::quaternion &value,
                             cvm::quaternion const &def_value,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_scalar_<cvm::quaternion>(conf, key, value, def_value, parse_mode);
}

bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             colvarvalue &value,
                             colvarvalue const &def_value,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_scalar_<colvarvalue>(conf, key, value, def_value, parse_mode);
}


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
    b_found = key_lookup(conf, key, &data_this, &save_pos);
    if (b_found) {
      if (!b_found_any)
        b_found_any = true;
      found_count++;
      data = data_this;
    }
  } while (b_found);

  if (found_count > 1) {
    cvm::error("Error: found more than one instance of \""+
               std::string(key)+"\".\n", INPUT_ERROR);
  }

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
      cvm::error("Error: boolean values only are allowed "
                 "for \""+std::string(key)+"\".\n", INPUT_ERROR);
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


// multiple-value keyword parsers

bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             std::vector<int> &values,
                             std::vector<int> const &def_values,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_vector_<int>(conf, key, values, def_values, parse_mode);
}

bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             std::vector<size_t> &values,
                             std::vector<size_t> const &def_values,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_vector_<size_t>(conf, key, values, def_values, parse_mode);
}

bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             std::vector<long> &values,
                             std::vector<long> const &def_values,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_vector_<long>(conf, key, values, def_values, parse_mode);
}

bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             std::vector<std::string> &values,
                             std::vector<std::string> const &def_values,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_vector_<std::string>(conf, key, values, def_values, parse_mode);
}

bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             std::vector<cvm::real> &values,
                             std::vector<cvm::real> const &def_values,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_vector_<cvm::real>(conf, key, values, def_values, parse_mode);
}

bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             std::vector<cvm::rvector> &values,
                             std::vector<cvm::rvector> const &def_values,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_vector_<cvm::rvector>(conf, key, values, def_values, parse_mode);
}

bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             std::vector<cvm::quaternion> &values,
                             std::vector<cvm::quaternion> const &def_values,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_vector_<cvm::quaternion>(conf, key, values, def_values, parse_mode);
}

bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             std::vector<colvarvalue> &values,
                             std::vector<colvarvalue> const &def_values,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_vector_<colvarvalue>(conf, key, values, def_values, parse_mode);
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


void colvarparse::clear_keyword_registry()
{
  allowed_keywords.clear();
  data_begin_pos.clear();
  data_end_pos.clear();
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
  while (cvm::getline(is, line)) {
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
      cvm::error("Error: keyword \""+uk+"\" is not supported, "
                 "or not recognized in this context.\n", INPUT_ERROR);
      return INPUT_ERROR;
    }
  }

  clear_keyword_registry();

  return COLVARS_OK;
}


std::istream & colvarparse::read_config_line(std::istream &is,
                                             std::string &line)
{
  cvm::getline(is, line);
  config_string += line+'\n';
  size_t const comment = line.find('#');
  if (comment != std::string::npos) {
    line.erase(comment);
  }
  return is;
}


std::istream & colvarparse::getline_nocomments(std::istream &is,
                                               std::string &line)
{
  cvm::getline(is, line);
  size_t const comment = line.find('#');
  if (comment != std::string::npos) {
    line.erase(comment);
  }
  return is;
}


bool colvarparse::key_lookup(std::string const &conf,
                             char const *key_in,
                             std::string *data,
                             size_t *save_pos)
{
  if (cvm::debug()) {
    cvm::log("Looking for the keyword \""+std::string(key_in)+
             "\" and its value.\n");
  }

  // add this keyword to the register (in its camelCase version)
  add_keyword(key_in);

  // use the lowercase version from now on
  std::string const key(to_lower_cppstr(key_in));

  // "conf_lower" is only used to lookup the keyword, but its value
  // will be read from "conf", in order not to mess up file names
  std::string const conf_lower(to_lower_cppstr(conf));

  // by default, there is no value, unless we found one
  if (data != NULL) {
    data->clear();
  }

  // start from the first occurrence of key
  size_t pos = conf_lower.find(key, (save_pos != NULL) ? *save_pos : 0);

  // iterate over all instances of the substring until it finds it as isolated keyword
  while (true) {

    if (pos == std::string::npos) {
      // no valid instance of the keyword has been found
      if (cvm::debug()) {
        cvm::log("Keyword \""+std::string(key_in)+"\" not found.\n");
      }
      return false;
    }

    bool b_isolated_left = true, b_isolated_right = true;

    if (pos > 0) {
      if ( std::string("\n"+std::string(white_space)+
                       "}").find(conf[pos-1]) ==
           std::string::npos ) {
        // none of the valid delimiting characters is on the left of key
        b_isolated_left = false;
      }
    }

    if (pos < conf.size()-key.size()-1) {
      if ( std::string("\n"+std::string(white_space)+
                       "{").find(conf[pos+key.size()]) ==
           std::string::npos ) {
        // none of the valid delimiting characters is on the right of key
        b_isolated_right = false;
      }
    }

    // check that there are matching braces between here and the end of conf
    bool const b_not_within_block = (check_braces(conf, pos) == COLVARS_OK);

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

  if (save_pos != NULL) {
  // save the pointer for a future call (when iterating over multiple
  // valid instances of the same keyword)
    *save_pos = pos + key.size();
  }

  // get the remainder of the line
  size_t pl = conf.rfind("\n", pos);
  size_t line_begin = (pl == std::string::npos) ? 0 : pos;
  size_t nl = conf.find("\n", pos);
  size_t line_end = (nl == std::string::npos) ? conf.size() : nl;
  std::string line(conf, line_begin, (line_end-line_begin));

  size_t data_begin = (to_lower_cppstr(line)).find(key) + key.size();
  data_begin = line.find_first_not_of(white_space, data_begin+1);

  if (data_begin != std::string::npos) {

    size_t data_end = line.find_last_not_of(white_space) + 1;
    data_end = (data_end == std::string::npos) ? line.size() : data_end;

    size_t brace = line.find('{', data_begin);  // look for an opening brace
    size_t brace_last = brace;

    if (brace != std::string::npos) {

      // find the matching closing brace

//       if (cvm::debug()) {
//         cvm::log("Multi-line value, config is now \""+line+"\".\n");
//       }

      int brace_count = 1;

      while (brace_count > 0) {

        brace = line.find_first_of("{}", brace_last+1);
        // find all braces within this line
        while (brace < std::string::npos) {
          brace_last = brace;
          if (line[brace] == '{') brace_count++;
          if (line[brace] == '}') brace_count--;
          if (brace_count == 0) {
            data_end = brace+1;
            break;
          }
          brace = line.find_first_of("{}", brace+1);
        }

        if (brace_count == 0) {
          data_end = brace+1;
          break;
        }

        if (brace == std::string::npos) {

          // add a new line
          if (line_end >= conf.size()) {
            cvm::error("Parse error: reached the end while "
                       "looking for closing brace; until now "
                       "the following was parsed: \"\n"+
                       line+"\".\n", INPUT_ERROR);
            return false;
          }

          line_begin = line_end;
          nl = conf.find('\n', line_begin+1);
          if (nl == std::string::npos)
            line_end = conf.size();
          else
            line_end = nl;
          line.append(conf, line_begin, (line_end-line_begin));

//           if (cvm::debug()) {
//             cvm::log("Added a new line, config is now \""+line+"\".\n");
//           }
        }

        if (brace_count < 0) {
          cvm::error("Error: found closing brace without opening brace.\n", INPUT_ERROR);
        }
      }

      // strip the leading and trailing braces
      data_begin = line.find_first_of('{') + 1;
      data_begin = line.find_first_not_of(white_space,
                                          data_begin);

      data_end = line.find_last_of('}', line.size()) - 1;
      data_end = line.find_last_not_of(white_space,
                                       data_end) + 1;
    }

    if (data != NULL) {
      data->append(line, data_begin, (data_end-data_begin));

      if (cvm::debug()) {
        cvm::log("Keyword value = \""+*data+"\".\n");
      }

      if (data->size()) {
        data_begin_pos.push_back(conf.find(*data, pos+key.size()));
        data_end_pos.push_back(data_begin_pos.back()+data->size());
      }
    }
  }

  if (save_pos != NULL) *save_pos = line_end;

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


int colvarparse::check_braces(std::string const &conf,
                              size_t const start_pos)
{
  int brace_count = 0;
  size_t brace = start_pos;
  while ((brace = conf.find_first_of("{}", brace)) != std::string::npos) {
    if (conf[brace] == '{') brace_count++;
    if (conf[brace] == '}') brace_count--;
    brace++;
  }
  return (brace_count != 0) ? INPUT_ERROR : COLVARS_OK;
}
