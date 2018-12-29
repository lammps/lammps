// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPARSE_H
#define COLVARPARSE_H

#include <cstring>
#include <string>

#include "colvarmodule.h"
#include "colvarvalue.h"


/// \file colvarparse.h Parsing functions for collective variables


/// \brief Base class containing parsing functions; all objects which
/// need to parse input inherit from this
class colvarparse {

protected:

  /// \brief List of legal keywords for this object: this is updated
  /// by each call to colvarparse::get_keyval() or
  /// colvarparse::key_lookup()
  std::list<std::string> allowed_keywords;

  /// \brief List of delimiters for the values of each keyword in the
  /// configuration string; all keywords will be stripped of their
  /// values before the keyword check is performed
  std::list<size_t>      data_begin_pos;

  /// \brief List of delimiters for the values of each keyword in the
  /// configuration string; all keywords will be stripped of their
  /// values before the keyword check is performed
  std::list<size_t>      data_end_pos;

  /// \brief Add a new valid keyword to the list
  void add_keyword(char const *key);

  /// \brief Remove all the values from the config string
  void strip_values(std::string &conf);

  /// \brief Configuration string of the object (includes comments)
  std::string config_string;

public:


  inline colvarparse()
  {
    init();
  }

  /// Constructor that stores the object's config string
  inline colvarparse(const std::string& conf)
  {
    init(conf);
  }

  /// Set the object ready to parse a new configuration string
  inline void init()
  {
    config_string.clear();
    clear_keyword_registry();
  }

  /// Set a new config string for this object
  inline void init(std::string const &conf)
  {
    if (! config_string.size()) {
      init();
      config_string = conf;
    }
  }

  /// Get the configuration string (includes comments)
  inline std::string const & get_config() const
  {
    return config_string;
  }

  /// How a keyword is parsed in a string
  enum Parse_Mode {
    /// \brief(default) Read the first instance of a keyword (if
    /// any), report its value, and print a warning when there is more
    /// than one
    parse_normal,
    /// \brief Like parse_normal, but don't send any message to the log
    /// (useful e.g. in restart files when such messages are very
    /// numerous and redundant)
    parse_silent
  };

  /// \brief Check that all the keywords within "conf" are in the list
  /// of allowed keywords; this will invoke strip_values() first and
  /// then loop over all words
  int check_keywords(std::string &conf, char const *key);

  /// \brief Use this after parsing a config string (note that check_keywords() calls it already)
  void clear_keyword_registry();

  /// \fn get_keyval bool const get_keyval (std::string const &conf,
  /// char const *key, _type_ &value, _type_ const &def_value,
  /// Parse_Mode const parse_mode) \brief Helper function to parse
  /// keywords in the configuration and get their values
  ///
  /// In normal circumstances, you should use either version the
  /// get_keyval function.  Both of them look for the C string "key"
  /// in the C++ string "conf", and assign the corresponding value (if
  /// available) to the variable "value" (first version), or assign as
  /// many values as found to the vector "values" (second version).
  ///
  /// If "key" is found but no value is associated to it, the default
  /// value is provided (either one copy or as many copies as the
  /// current length of the vector "values" specifies).  A message
  /// will print, unless parse_mode is equal to parse_silent.  The
  /// return value of both forms of get_keyval is true if "key" is
  /// found (with or without value), and false when "key" is absent in
  /// the string "conf".  If there is more than one instance of the
  /// keyword, a warning will be raised; instead, to loop over
  /// multiple instances key_lookup() should be invoked directly.
  ///
  /// If you introduce a new data type, add two new instances of this
  /// functions, or insert this type in the \link colvarvalue \endlink
  /// wrapper class (colvarvalue.h).

  bool get_keyval(std::string const &conf,
                  char const *key,
                  int &value,
                  int const &def_value = (int)0,
                  Parse_Mode const parse_mode = parse_normal);
  bool get_keyval(std::string const &conf,
                  char const *key,
                  size_t &value,
                  size_t const &def_value = (size_t)0,
                  Parse_Mode const parse_mode = parse_normal);
  bool get_keyval(std::string const &conf,
                  char const *key,
                  long &value,
                  long const &def_value = 0,
                  Parse_Mode const parse_mode = parse_normal);
  bool get_keyval(std::string const &conf,
                  char const *key,
                  std::string &value,
                  std::string const &def_value = std::string(""),
                  Parse_Mode const parse_mode = parse_normal);
  bool get_keyval(std::string const &conf,
                  char const *key,
                  cvm::real &value,
                  cvm::real const &def_value = (cvm::real)0.0,
                  Parse_Mode const parse_mode = parse_normal);
  bool get_keyval(std::string const &conf,
                  char const *key,
                  cvm::rvector &value,
                  cvm::rvector const &def_value = cvm::rvector(),
                  Parse_Mode const parse_mode = parse_normal);
  bool get_keyval(std::string const &conf,
                  char const *key,
                  cvm::quaternion &value,
                  cvm::quaternion const &def_value = cvm::quaternion(),
                  Parse_Mode const parse_mode = parse_normal);
  bool get_keyval(std::string const &conf,
                  char const *key,
                  colvarvalue &value,
                  colvarvalue const &def_value = colvarvalue(colvarvalue::type_notset),
                  Parse_Mode const parse_mode = parse_normal);
  bool get_keyval(std::string const &conf,
                  char const *key,
                  bool &value,
                  bool const &def_value = false,
                  Parse_Mode const parse_mode = parse_normal);
  bool get_keyval(std::string const &conf,
                  char const *key,
                  std::vector<int> &values,
                  std::vector<int> const &def_values = std::vector<int>(0, (int)0),
                  Parse_Mode const parse_mode = parse_normal);
  bool get_keyval(std::string const &conf,
                  char const *key,
                  std::vector<size_t> &values,
                  std::vector<size_t> const &def_values = std::vector<size_t>(0, (size_t)0),
                  Parse_Mode const parse_mode = parse_normal);
  bool get_keyval(std::string const &conf,
                  char const *key,
                  std::vector<long> &values,
                  std::vector<long> const &def_values = std::vector<long>(0, (long)0),
                  Parse_Mode const parse_mode = parse_normal);
  bool get_keyval(std::string const &conf,
                  char const *key,
                  std::vector<std::string> &values,
                  std::vector<std::string> const &def_values = std::vector<std::string>(0, std::string("")),
                  Parse_Mode const parse_mode = parse_normal);
  bool get_keyval(std::string const &conf,
                  char const *key,
                  std::vector<cvm::real> &values,
                  std::vector<cvm::real> const &def_values = std::vector<cvm::real>(0, (cvm::real)0.0),
                  Parse_Mode const parse_mode = parse_normal);
  bool get_keyval(std::string const &conf,
                  char const *key,
                  std::vector<cvm::rvector> &values,
                  std::vector<cvm::rvector> const &def_values = std::vector<cvm::rvector>(0, cvm::rvector()),
                  Parse_Mode const parse_mode = parse_normal);
  bool get_keyval(std::string const &conf,
                  char const *key,
                  std::vector<cvm::quaternion> &values,
                  std::vector<cvm::quaternion> const &def_values = std::vector<cvm::quaternion>(0, cvm::quaternion()),
                  Parse_Mode const parse_mode = parse_normal);
  bool get_keyval(std::string const &conf,
                  char const *key,
                  std::vector<colvarvalue> &values,
                  std::vector<colvarvalue> const &def_values = std::vector<colvarvalue>(0, colvarvalue(colvarvalue::type_notset)),
                  Parse_Mode const parse_mode = parse_normal);

protected:

  // Templates
  template<typename TYPE> bool _get_keyval_scalar_(std::string const &conf,
                                                   char const *key,
                                                   TYPE &value,
                                                   TYPE const &def_value,
                                                   Parse_Mode const parse_mode);
  bool _get_keyval_scalar_string_(std::string const &conf,
                                  char const *key,
                                  std::string &value,
                                  std::string const &def_value,
                                  Parse_Mode const parse_mode);

  template<typename TYPE> bool _get_keyval_vector_(std::string const &conf,
                                                   char const *key,
                                                   std::vector<TYPE> &values,
                                                   std::vector<TYPE> const &def_values,
                                                   Parse_Mode const parse_mode);

public:

  /// \brief Return a lowercased copy of the string
  static inline std::string to_lower_cppstr(std::string const &in)
  {
    std::string out = "";
    for (size_t i = 0; i < in.size(); i++) {
      out.append(1, (char) ::tolower(in[i]) );
    }
    return out;
  }

  /// \brief Helper class to read a block of the type "key { ... }"
  /// from a stream and store it in a string
  ///
  /// Useful on restarts, where the file is too big to be loaded in a
  /// string by key_lookup; it can only check that the keyword is
  /// correct and the block is properly delimited by braces, not
  /// skipping other blocks
  class read_block {

    std::string const   key;
    std::string * const data;

  public:
    inline read_block(std::string const &key_in, std::string &data_in)
      : key(key_in), data(&data_in)
    {}
    inline ~read_block() {}
    friend std::istream & operator >> (std::istream &is, read_block const &rb);
  };


  /// Accepted white space delimiters, used in key_lookup()
  static const char * const white_space;

  /// \brief Low-level function for parsing configuration strings;
  /// automatically adds the requested keyword to the list of valid
  /// ones.  \param conf the content of the configuration file or one
  /// of its blocks \param key the keyword to search within "conf" \param
  /// data (optional) holds the string provided after "key", if any
  /// \param save_pos (optional) stores the position of the keyword
  /// within "conf", useful when doing multiple calls
  bool key_lookup(std::string const &conf,
                  char const *key,
                  std::string *data = NULL,
                  size_t *save_pos = NULL);

  /// \brief Reads a configuration line, adds it to config_string, and returns
  /// the stream \param is Input stream \param s String that will hold the
  /// configuration line, with comments stripped
  std::istream & read_config_line(std::istream &is, std::string &line);

  /// \brief Works as std::getline() but also removes everything
  /// between a comment character and the following newline
  static std::istream & getline_nocomments(std::istream &is, std::string &s);

  /// \brief Check if the content of a config string has matching braces
  /// \param conf The configuration string \param start_pos Start the count
  /// from this position
  static int check_braces(std::string const &conf, size_t const start_pos);

};


#endif
