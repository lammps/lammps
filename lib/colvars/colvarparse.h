/// -*- c++ -*-

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

  /// \brief Whether or not to accumulate data_begin_pos and
  /// data_end_pos in key_lookup(); it may be useful to disable
  /// this after the constructor is called, because other files may be
  /// read (e.g. restart) that would mess up with the registry; in any
  /// case, nothing serious happens until check_keywords() is invoked
  /// (which should happen only right after construction)
  bool save_delimiters;

  /// \brief Add a new valid keyword to the list
  void add_keyword (char const *key);

  /// \brief Remove all the values from the config string
  void strip_values (std::string &conf);

public:

  inline colvarparse()
    : save_delimiters (true)
  {}

  /// How a keyword is parsed in a string
  enum Parse_Mode {
    /// \brief (default) Read the first instance of a keyword (if
    /// any), report its value, and print a warning when there is more
    /// than one
    parse_normal,
    /// \brief Like parse_normal, but don't send any message to the log
    /// (useful e.g. in restart files when such messages are very
    /// numerous and redundant)
    parse_silent
  };

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

#define _get_keyval_scalar_proto_(_type_,_def_value_)           \
  bool get_keyval (std::string const &conf,                     \
                   char const *key,                             \
                   _type_ &value,                               \
                   _type_ const &def_value = _def_value_,       \
                   Parse_Mode const parse_mode = parse_normal)

    _get_keyval_scalar_proto_ (int, (int)0);
    _get_keyval_scalar_proto_ (size_t, (size_t)0);
    _get_keyval_scalar_proto_ (std::string, std::string (""));
    _get_keyval_scalar_proto_ (cvm::real, (cvm::real)0.0);
    _get_keyval_scalar_proto_ (cvm::rvector, cvm::rvector());
    _get_keyval_scalar_proto_ (cvm::quaternion, cvm::quaternion());
    _get_keyval_scalar_proto_ (colvarvalue, colvarvalue (colvarvalue::type_notset));
    _get_keyval_scalar_proto_ (bool, false);

#define _get_keyval_vector_proto_(_type_,_def_value_)                   \
  bool get_keyval (std::string const &conf,                             \
                   char const *key,                                     \
                   std::vector<_type_> &values,                         \
                   std::vector<_type_> const &def_values =              \
                   std::vector<_type_> (0, static_cast<_type_>(_def_value_)),                \
                   Parse_Mode const parse_mode = parse_normal)

    _get_keyval_vector_proto_ (int, 0);
    _get_keyval_vector_proto_ (size_t, 0);
    _get_keyval_vector_proto_ (std::string, std::string (""));
    _get_keyval_vector_proto_ (cvm::real, 0.0);
    _get_keyval_vector_proto_ (cvm::rvector, cvm::rvector());
    _get_keyval_vector_proto_ (cvm::quaternion, cvm::quaternion());
    _get_keyval_vector_proto_ (colvarvalue, colvarvalue (colvarvalue::type_notset));


  /// \brief Check that all the keywords within "conf" are in the list
  /// of allowed keywords; this will invoke strip_values() first and
  /// then loop over all words
  int check_keywords (std::string &conf, char const *key);


  /// \brief Return a lowercased copy of the string
  static inline std::string to_lower_cppstr (std::string const &in)
  {
    std::string out = "";
    for (size_t i = 0; i < in.size(); i++) {
      out.append (1, (char) ::tolower (in[i]) );
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
    inline read_block (std::string const &key_in, std::string &data_in)
      : key (key_in), data (&data_in)
    {}
    inline ~read_block() {}
    friend std::istream & operator >> (std::istream &is, read_block const &rb);
  };


  /// Accepted white space delimiters, used in key_lookup()
  static std::string const white_space;

  /// \brief Low-level function for parsing configuration strings;
  /// automatically adds the requested keywords to the list of valid
  /// ones.  \param conf the content of the configuration file or one
  /// of its blocks \param key the keyword to search in "conf" \param
  /// data (optional) holds the string provided after "key", if any
  /// \param save_pos (optional) stores the position of the keyword
  /// within "conf", useful when doing multiple calls \param
  /// save_delimiters (optional)
  bool key_lookup (std::string const &conf,
                   char const *key,
                   std::string &data = dummy_string,
                   size_t &save_pos = dummy_pos);

  /// Used as a default argument by key_lookup
  static std::string dummy_string;
  /// Used as a default argument by key_lookup
  static size_t dummy_pos;

  /// \brief Works as std::getline() but also removes everything
  /// between a comment character and the following newline
  static std::istream & getline_nocomments (std::istream &is,
                                            std::string &s,
                                            char const delim = '\n');

  /// Check if the content of the file has matching braces
  bool brace_check (std::string const &conf,
                    size_t const start_pos = 0);

};


#endif
