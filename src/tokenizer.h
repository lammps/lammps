/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Richard Berger (Temple U)
------------------------------------------------------------------------- */

#ifndef LMP_TOKENIZER_H
#define LMP_TOKENIZER_H

#include "lmptype.h"

#include <exception>
#include <string>
#include <vector>

namespace LAMMPS_NS {

#define TOKENIZER_DEFAULT_SEPARATORS " \t\r\n\f"

class Tokenizer {
  std::string text;
  std::string separators;
  size_t start;
  size_t ntokens;

 public:
  Tokenizer(std::string str, std::string separators = TOKENIZER_DEFAULT_SEPARATORS);
  Tokenizer(Tokenizer &&);
  Tokenizer(const Tokenizer &);
  Tokenizer &operator=(const Tokenizer &);
  Tokenizer &operator=(Tokenizer &&);
  void swap(Tokenizer &);

  void reset();
  void skip(int n = 1);
  bool has_next() const;
  bool contains(const std::string &str) const;
  std::string next();

  size_t count();
  std::vector<std::string> as_vector();
};

/** General Tokenizer exception class */

class TokenizerException : public std::exception {
  std::string message;

 public:
  // remove unused default constructor
  TokenizerException() = delete;

  /** Thrown during retrieving or skipping tokens
   *
   * \param  msg    String with error message
   * \param  token  String of the token/word that caused the error */
  explicit TokenizerException(const std::string &msg, const std::string &token);

  /** Retrieve message describing the thrown exception
   * \return string with error message */
  const char *what() const noexcept override { return message.c_str(); }
};

/** Exception thrown by ValueTokenizer when trying to convert an invalid integer string */

class InvalidIntegerException : public TokenizerException {

 public:
  /** Thrown during converting string to integer number
   *
   * \param  token  String of the token/word that caused the error */
  explicit InvalidIntegerException(const std::string &token) :
      TokenizerException("Not a valid integer number", token)
  {
  }
};

/** Exception thrown by ValueTokenizer when trying to convert an floating point string */

class InvalidFloatException : public TokenizerException {
 public:
  /** Thrown during converting string to floating point number
   *
   * \param  token  String of the token/word that caused the error */
  explicit InvalidFloatException(const std::string &token) :
      TokenizerException("Not a valid floating-point number", token)
  {
  }
};

class ValueTokenizer {
  Tokenizer tokens;

 public:
  ValueTokenizer(const std::string &str,
                 const std::string &separators = TOKENIZER_DEFAULT_SEPARATORS);
  ValueTokenizer(const ValueTokenizer &) = default;
  ValueTokenizer(ValueTokenizer &&);
  ValueTokenizer &operator=(const ValueTokenizer &);
  ValueTokenizer &operator=(ValueTokenizer &&);
  void swap(ValueTokenizer &);

  std::string next_string();
  tagint next_tagint();
  bigint next_bigint();
  int next_int();
  double next_double();

  bool has_next() const;
  bool contains(const std::string &value) const;
  void skip(int ntokens = 1);

  size_t count();
};

}    // namespace LAMMPS_NS

#endif
