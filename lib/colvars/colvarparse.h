// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPARSE_H
#define COLVARPARSE_H

#include <cstring>
#include <string>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparams.h"


/// \file colvarparse.h Parsing functions for collective variables


/// \brief Base class containing parsing functions; all objects which
/// need to parse input inherit from this
class colvarparse : public colvarparams {

public:

  /// Default constructor
  colvarparse();

  /// Constructor that stores the object's config string
  colvarparse(const std::string& conf);

  /// Set the object ready to parse a new configuration string
  void init();

  /// Set a new config string for this object
  void init(std::string const &conf);

  /// Default destructor
  virtual ~colvarparse();

  /// Get the configuration string (includes comments)
  inline std::string const & get_config() const
  {
    return config_string;
  }

  /// How a keyword is parsed in a string
  enum Parse_Mode {
    /// Zero for all flags
    parse_null = 0,
    /// Print the value of a keyword if it is given
    parse_echo = (1<<1),
    /// Print the default value of a keyword, if it is NOT given
    parse_echo_default = (1<<2),
    /// Do not print the keyword
    parse_silent = 0,
    /// Raise error if the keyword is not provided
    parse_required = (1<<16),
    /// Successive calls to get_keyval() will override the previous values
    /// when the keyword is not given any more
    parse_override = (1<<17),
    /// The call is being executed from a read_restart() function
    parse_restart = (1<<18),
    /// Alias for old default behavior (should be phased out)
    parse_normal = (1<<2) | (1<<1) | (1<<17)
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
                  cvm::step_number &value,
                  cvm::step_number const &def_value = 0,
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

  /// Get the string value of a keyword, and save it for later parsing
  bool get_key_string_value(std::string const &conf,
                            char const *key, std::string &data);

  /// Get multiple strings from repeated instances of a same keyword
  bool get_key_string_multi_value(std::string const &conf,
                                  char const *key, std::vector<std::string>& data);

  /// Template for single-value keyword parsers
  template<typename TYPE>
  bool _get_keyval_scalar_(std::string const &conf,
                           char const *key,
                           TYPE &value,
                           TYPE const &def_value,
                           Parse_Mode const &parse_mode);

  /// Template for multiple-value keyword parsers
  template<typename TYPE>
  bool _get_keyval_vector_(std::string const &conf,
                           char const *key,
                           std::vector<TYPE> &values,
                           std::vector<TYPE> const &def_values,
                           Parse_Mode const &parse_mode);

  /// Extract the value of a variable from a string
  template<typename TYPE>
  int _get_keyval_scalar_value_(std::string const &key_str,
                                std::string const &data,
                                TYPE &value,
                                TYPE const &def_value);

  /// Handle the case where the user provides a keyword without value
  template<typename TYPE>
  int _get_keyval_scalar_novalue_(std::string const &key_str,
                                  TYPE &value,
                                  Parse_Mode const &parse_mode);

  /// Record that the keyword has just been user-defined
  template<typename TYPE>
  void mark_key_set_user(std::string const &key_str,
                         TYPE const &value,
                         Parse_Mode const &parse_mode);

  /// Record that the keyword has just been set to its default value
  template<typename TYPE>
  void mark_key_set_default(std::string const &key_str,
                            TYPE const &def_value,
                            Parse_Mode const &parse_mode);

  /// Raise error condition due to the keyword being required!
  void error_key_required(std::string const &key_str,
                          Parse_Mode const &parse_mode);

  /// True if the keyword has been set already
  bool key_already_set(std::string const &key_str);

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
  /// the stream \param is Input stream \param line String that will hold the
  /// configuration line, with comments stripped
  std::istream & read_config_line(std::istream &is, std::string &line);

  /// \brief Works as std::getline() but also removes everything
  /// between a comment character and the following newline
  static std::istream & getline_nocomments(std::istream &is, std::string &s);

  /// \brief Check if the content of a config string has matching braces
  /// \param conf The configuration string \param start_pos Start the count
  /// from this position
  static int check_braces(std::string const &conf, size_t const start_pos);

  /// \brief Split a string with a specified delimiter into a vector
  /// \param data The string to be splitted
  /// \param delim A delimiter
  /// \param dest A destination vector to store the splitted results
  static void split_string(const std::string& data, const std::string& delim, std::vector<std::string>& dest);

protected:

  /// \brief List of legal keywords for this object: this is updated
  /// by each call to colvarparse::get_keyval() or
  /// colvarparse::key_lookup()
  std::list<std::string> allowed_keywords;

  /// How a keyword has been set
  enum key_set_mode {
    key_not_set = 0,
    key_set_user = 1,
    key_set_default = 2
  };

  /// Track which keywords have been already set, and how
  std::map<std::string, key_set_mode> key_set_modes;

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

};


/// Bitwise OR between two Parse_mode flags
inline colvarparse::Parse_Mode operator | (colvarparse::Parse_Mode const &mode1,
                                           colvarparse::Parse_Mode const &mode2)
{
  return static_cast<colvarparse::Parse_Mode>(static_cast<int>(mode1) |
                                              static_cast<int>(mode2));
}

#endif
