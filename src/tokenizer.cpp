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

#include "tokenizer.h"
#include "fmt/format.h"
#include "utils.h"

#include <utility>

using namespace LAMMPS_NS;

TokenizerException::TokenizerException(const std::string &msg, const std::string &token)
{
  if (token.empty()) {
    message = msg;
  } else {
    message = fmt::format("{}: '{}'", msg, token);
  }
}

/** Class for splitting text into words
 *
 * This tokenizer will break down a string into sub-strings (i.e words)
 * separated by the given separator characters. If the string contains
 * certain known UTF-8 characters they will be replaced by their ASCII
 * equivalents processing the string.
 *
\verbatim embed:rst

*See also*
   :cpp:class:`ValueTokenizer`, :cpp:func:`utils::split_words`, :cpp:func:`utils::utf8_subst`

\endverbatim
 *
 * \param  str          string to be processed
 * \param  _separators  string with separator characters (default: " \t\r\n\f") */

Tokenizer::Tokenizer(std::string str, std::string _separators) :
    text(std::move(str)), separators(std::move(_separators)), start(0), ntokens(std::string::npos)
{
  // replace known UTF-8 characters with ASCII equivalents
  if (utils::has_utf8(text)) text = utils::utf8_subst(text);
  reset();
}

Tokenizer::Tokenizer(const Tokenizer &rhs) :
    text(rhs.text), separators(rhs.separators), ntokens(rhs.ntokens)
{
  reset();
}

Tokenizer::Tokenizer(Tokenizer &&rhs) :
    text(std::move(rhs.text)), separators(std::move(rhs.separators)), ntokens(rhs.ntokens)
{
  reset();
}

Tokenizer &Tokenizer::operator=(const Tokenizer &other)
{
  Tokenizer tmp(other);
  swap(tmp);
  return *this;
}

Tokenizer &Tokenizer::operator=(Tokenizer &&other)
{
  Tokenizer tmp(std::move(other));
  swap(tmp);
  return *this;
}

void Tokenizer::swap(Tokenizer &other)
{
  std::swap(text, other.text);
  std::swap(separators, other.separators);
  std::swap(start, other.start);
  std::swap(ntokens, other.ntokens);
}

/*! Re-position the tokenizer state to the first word,
 * i.e. the first non-separator character */
void Tokenizer::reset()
{
  start = text.find_first_not_of(separators);
}

/*! Search the text to be processed for a sub-string.
 *
 * \param  str  string to be searched for
 * \return      true if string was found, false if not */
bool Tokenizer::contains(const std::string &str) const
{
  return text.find(str) != std::string::npos;
}

/*! Skip over a given number of tokens
 *
 * \param  n  number of tokens to skip over */
void Tokenizer::skip(int n)
{
  for (int i = 0; i < n; ++i) {
    if (!has_next()) throw TokenizerException("No more tokens", "");

    size_t end = text.find_first_of(separators, start);

    if (end == std::string::npos) {
      start = end;
    } else {
      start = text.find_first_not_of(separators, end + 1);
    }
  }
}

/*! Indicate whether more tokens are available
 *
 * \return   true if there are more tokens, false if not */
bool Tokenizer::has_next() const
{
  return start != std::string::npos;
}

/*! Retrieve next token.
 *
 * \return   string with the next token */
std::string Tokenizer::next()
{
  if (!has_next()) throw TokenizerException("No more tokens", "");

  size_t end = text.find_first_of(separators, start);

  if (end == std::string::npos) {
    std::string token = text.substr(start);
    start = end;
    return token;
  }

  std::string token = text.substr(start, end - start);
  start = text.find_first_not_of(separators, end + 1);
  return token;
}

/*! Count number of tokens in text.
 *
 * \return   number of counted tokens */
size_t Tokenizer::count()
{
  // lazy evaluation
  if (ntokens == std::string::npos) { ntokens = utils::count_words(text, separators); }
  return ntokens;
}

/*! Retrieve the entire text converted to an STL vector of tokens.
 *
 * \return   The STL vector */
std::vector<std::string> Tokenizer::as_vector()
{
  // store current state
  size_t current = start;

  reset();

  // generate vector
  std::vector<std::string> tokens;

  while (has_next()) { tokens.emplace_back(next()); }

  // restore state
  start = current;

  return tokens;
}

/*! Class for reading text with numbers
 *
\verbatim embed:rst

*See also*
   :cpp:class:`Tokenizer`

\endverbatim
 *
 * \param str         String to be processed
 * \param separators  String with separator characters (default: " \t\r\n\f")
 *
 * \see Tokenizer InvalidIntegerException InvalidFloatException */

ValueTokenizer::ValueTokenizer(const std::string &str, const std::string &separators) :
    tokens(str, separators)
{
}

ValueTokenizer::ValueTokenizer(ValueTokenizer &&rhs) : tokens(std::move(rhs.tokens)) {}

ValueTokenizer &ValueTokenizer::operator=(const ValueTokenizer &other)
{
  ValueTokenizer tmp(other);
  swap(tmp);
  return *this;
}

ValueTokenizer &ValueTokenizer::operator=(ValueTokenizer &&other)
{
  ValueTokenizer tmp(std::move(other));
  swap(tmp);
  return *this;
}

void ValueTokenizer::swap(ValueTokenizer &other)
{
  std::swap(tokens, other.tokens);
}

/*! Indicate whether more tokens are available
 *
 * \return   true if there are more tokens, false if not */
bool ValueTokenizer::has_next() const
{
  return tokens.has_next();
}

/*! Search the text to be processed for a sub-string.
 *
 * \param  value  string with value to be searched for
 * \return        true if string was found, false if not */
bool ValueTokenizer::contains(const std::string &value) const
{
  return tokens.contains(value);
}

/*! Retrieve next token
 *
 * \return   string with next token */
std::string ValueTokenizer::next_string()
{
  return tokens.next();
}

/*! Retrieve next token and convert to int
 *
 * \return   value of next token */
int ValueTokenizer::next_int()
{
  std::string current = tokens.next();
  if (!utils::is_integer(current)) { throw InvalidIntegerException(current); }
  return atoi(current.c_str());
}

/*! Retrieve next token and convert to bigint
 *
 * \return   value of next token */
bigint ValueTokenizer::next_bigint()
{
  std::string current = tokens.next();
  if (!utils::is_integer(current)) { throw InvalidIntegerException(current); }
  return ATOBIGINT(current.c_str());
}

/*! Retrieve next token and convert to tagint
 *
 * \return   value of next token */
tagint ValueTokenizer::next_tagint()
{
  std::string current = tokens.next();
  if (!utils::is_integer(current)) { throw InvalidIntegerException(current); }
  return ATOTAGINT(current.c_str());
}

/*! Retrieve next token and convert to double
 *
 * \return   value of next token */
double ValueTokenizer::next_double()
{
  std::string current = tokens.next();
  if (!utils::is_double(current)) { throw InvalidFloatException(current); }
  return atof(current.c_str());
}

/*! Skip over a given number of tokens
 *
 * \param  n  number of tokens to skip over */
void ValueTokenizer::skip(int n)
{
  tokens.skip(n);
}

/*! Count number of tokens in text.
 *
 * \return   number of counted tokens */
size_t ValueTokenizer::count()
{
  return tokens.count();
}
