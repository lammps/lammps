// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <sstream>
#include <iostream>
#include <algorithm>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"


// space & tab
char const * const colvarparse::white_space = " \t";


namespace {

  // Avoid having to put the bool assignment in the template :-(
  void set_bool(void *p, bool val)
  {
    bool *v = reinterpret_cast<bool *>(p);
    *v = val;
  }

}


colvarparse::colvarparse()
{
  init();
}


void colvarparse::init()
{
  config_string.clear();
  clear_keyword_registry();
}


colvarparse::colvarparse(const std::string& conf)
{
  init(conf);
}


void colvarparse::init(std::string const &conf)
{
  if (! config_string.size()) {
    init();
    config_string = conf;
  }
}


colvarparse::~colvarparse()
{
  init();
}



bool colvarparse::get_key_string_value(std::string const &conf,
                                       char const *key, std::string &data)
{
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

  return b_found_any;
}

bool colvarparse::get_key_string_multi_value(std::string const &conf,
                                             char const *key, std::vector<std::string>& data)
{
  bool b_found = false, b_found_any = false;
  size_t save_pos = 0, found_count = 0;

  data.clear();

  do {
    std::string data_this = "";
    b_found = key_lookup(conf, key, &data_this, &save_pos);
    if (b_found) {
      if (!b_found_any)
        b_found_any = true;
      found_count++;
      data.push_back(data_this);
    }
  } while (b_found);

  return b_found_any;
}


template<typename TYPE>
void colvarparse::mark_key_set_user(std::string const &key_str,
                                    TYPE const &value,
                                    Parse_Mode const &parse_mode)
{
  key_set_modes[to_lower_cppstr(key_str)] = key_set_user;
  if (parse_mode & parse_echo) {
    cvm::log("# "+key_str+" = "+cvm::to_str(value)+"\n",
             cvm::log_user_params());
  }
}


template<typename TYPE>
void colvarparse::mark_key_set_default(std::string const &key_str,
                                       TYPE const &def_value,
                                       Parse_Mode const &parse_mode)
{
  key_set_modes[to_lower_cppstr(key_str)] = key_set_default;
  if (parse_mode & parse_echo_default) {
    cvm::log("# "+key_str+" = "+cvm::to_str(def_value)+
             " [default]\n", cvm::log_default_params());
  }
}


void colvarparse::error_key_required(std::string const &key_str,
                                     Parse_Mode const &parse_mode)
{
  if (key_already_set(key_str)) {
    return;
  }
  if (parse_mode & parse_restart) {
    cvm::error("Error: keyword \""+key_str+
               "\" is missing from the restart.\n", INPUT_ERROR);
  } else {
    cvm::error("Error: keyword \""+key_str+
               "\" is required.\n", INPUT_ERROR);
  }
}


template<typename TYPE>
int colvarparse::_get_keyval_scalar_value_(std::string const &key_str,
                                           std::string const &data,
                                           TYPE &value,
                                           TYPE const &def_value)
{
  std::istringstream is(data);
  size_t value_count = 0;
  TYPE x(def_value);

  while (is >> x) {
    value = x;
    value_count++;
  }

  if (value_count == 0) {
    return cvm::error("Error: in parsing \""+
                      key_str+"\".\n", INPUT_ERROR);
  }

  if (value_count > 1) {
    return cvm::error("Error: multiple values "
                      "are not allowed for keyword \""+
                      key_str+"\".\n", INPUT_ERROR);
  }

  return COLVARS_OK;
}


template<>
int colvarparse::_get_keyval_scalar_value_(std::string const &key_str,
                                           std::string const &data,
                                           bool &value,
                                           bool const & /* def_value */)
{
  if ( (data == std::string("on")) ||
       (data == std::string("yes")) ||
       (data == std::string("true")) ) {
    set_bool(reinterpret_cast<void *>(&value), true);
  } else if ( (data == std::string("off")) ||
              (data == std::string("no")) ||
              (data == std::string("false")) ) {
    set_bool(reinterpret_cast<void *>(&value), false);
  } else {
    return cvm::error("Error: boolean values only are allowed "
                      "for \""+key_str+"\".\n", INPUT_ERROR);
  }
  return COLVARS_OK;
}


template<typename TYPE>
int colvarparse::_get_keyval_scalar_novalue_(std::string const &key_str,
                                             TYPE & /* value */,
                                             Parse_Mode const & /* parse_mode */)
{
  return cvm::error("Error: improper or missing value "
                    "for \""+key_str+"\".\n", INPUT_ERROR);
}

template<>
int colvarparse::_get_keyval_scalar_novalue_(std::string const &key_str,
                                             bool &value,
                                             Parse_Mode const &parse_mode)
{
  set_bool(reinterpret_cast<void *>(&value), true);
  mark_key_set_user<bool>(key_str, value, parse_mode);
  return COLVARS_OK;
}


template<typename TYPE>
bool colvarparse::_get_keyval_scalar_(std::string const &conf,
                                      char const *key,
                                      TYPE &value,
                                      TYPE const &def_value,
                                      Parse_Mode const &parse_mode)
{
  std::string const key_str(key);

  std::string data;
  bool const b_found_any = get_key_string_value(conf, key, data);

  if (data.size()) {

    _get_keyval_scalar_value_<TYPE>(key_str, data, value, def_value);

    mark_key_set_user<TYPE>(key_str, value, parse_mode);

  } else { // No string value

    if (b_found_any) {

      _get_keyval_scalar_novalue_<TYPE>(key_str, value, parse_mode);

    } else {

      if (parse_mode & parse_required) {
        if (cvm::debug()) {
          cvm::log("get_keyval, parse_required = "+cvm::to_str(parse_mode & parse_required)+
                   "\n");
        }
        error_key_required(key_str, parse_mode);
        return false;
      }

      if ( (parse_mode & parse_override) || !(key_already_set(key)) ) {
        value = def_value;
        mark_key_set_default<TYPE>(key_str, value, parse_mode);
      }
    }
  }

  return b_found_any;
}


template<typename TYPE>
bool colvarparse::_get_keyval_vector_(std::string const &conf,
                                      char const *key,
                                      std::vector<TYPE> &values,
                                      std::vector<TYPE> const &def_values,
                                      Parse_Mode const &parse_mode)
{
  std::string const key_str(key);

  std::string data;
  bool const b_found_any = get_key_string_value(conf, key, data);

  if (data.size()) {
    std::istringstream is(data);

    if (values.size() == 0) {

      std::vector<TYPE> x;
      if (def_values.size()) {
        x = def_values;
      } else {
        x.assign(1, TYPE());
      }

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
                     key_str+"\".\n", INPUT_ERROR);
        }
      }
    }

    mark_key_set_user< std::vector<TYPE> >(key_str, values, parse_mode);

  } else {

    if (b_found_any) {
      cvm::error("Error: improper or missing values for \""+
                 key_str+"\".\n", INPUT_ERROR);
    } else {

      if ((values.size() > 0) && (values.size() != def_values.size())) {
        cvm::error("Error: the number of default values for \""+
                   key_str+"\" is different from the number of "
                   "current values.\n", BUG_ERROR);
      }

      if (parse_mode & parse_required) {
        error_key_required(key_str, parse_mode);
        return false;
      }

      if ( (parse_mode & parse_override) || !(key_already_set(key)) ) {
        for (size_t i = 0; i < values.size(); i++) {
          values[i] = def_values[i];
        }
        mark_key_set_default< std::vector<TYPE> >(key_str, def_values,
                                                  parse_mode);
      }

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
                             cvm::step_number &value,
                             cvm::step_number const &def_value,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_scalar_<cvm::step_number>(conf, key, value, def_value, parse_mode);
}

bool colvarparse::get_keyval(std::string const &conf,
                             char const *key,
                             std::string &value,
                             std::string const &def_value,
                             Parse_Mode const parse_mode)
{
  return _get_keyval_scalar_<std::string>(conf, key, value, def_value, parse_mode);
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
  return _get_keyval_scalar_<bool>(conf, key, value, def_value, parse_mode);
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
  std::string const key_str_lower(to_lower_cppstr(std::string(key)));

  if (key_set_modes.find(key_str_lower) != key_set_modes.end()) {
    return;
  }

  key_set_modes[key_str_lower] = key_not_set;

  allowed_keywords.push_back(key_str_lower);
}


bool colvarparse::key_already_set(std::string const &key_str)
{
  std::string const key_str_lower(to_lower_cppstr(key_str));

  if (key_set_modes.find(key_str_lower) == key_set_modes.end()) {
    return false;
  }

  return (key_set_modes[key_str_lower] > 0);
}


void colvarparse::strip_values(std::string &conf)
{
  size_t offset = 0;
  data_begin_pos.sort();
  data_end_pos.sort();
  std::list<size_t>::iterator data_begin_pos_last = std::unique(data_begin_pos.begin(), data_begin_pos.end());
  data_begin_pos.erase(data_begin_pos_last, data_begin_pos.end());
  std::list<size_t>::iterator data_end_pos_last = std::unique(data_end_pos.begin(), data_end_pos.end());
  data_end_pos.erase(data_end_pos_last, data_end_pos.end());

  std::list<size_t>::iterator data_begin = data_begin_pos.begin();
  std::list<size_t>::iterator data_end   = data_end_pos.begin();

  for ( ; (data_begin != data_begin_pos.end()) &&
          (data_end   != data_end_pos.end()) ;
        data_begin++, data_end++) {
    size_t const nchars = *data_end-*data_begin;
    conf.erase(*data_begin - offset, nchars);
    offset += nchars;
  }
}


void colvarparse::clear_keyword_registry()
{
  key_set_modes.clear();
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

void colvarparse::split_string(const std::string& data, const std::string& delim, std::vector<std::string>& dest) {
    size_t index = 0, new_index = 0;
    std::string tmpstr;
    while (index != data.length()) {
        new_index = data.find(delim, index);
        if (new_index != std::string::npos) tmpstr = data.substr(index, new_index - index);
        else tmpstr = data.substr(index, data.length());
        if (!tmpstr.empty()) {
            dest.push_back(tmpstr);
        }
        if (new_index == std::string::npos) break;
        index = new_index + 1;
    }
}
