// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.


#include "colvarmodule.h"
#include "colvartypes.h"
#include "colvarvalue.h"
#include "colvars_memstream.h"


bool cvm::memory_stream::expand_output_buffer(size_t add_bytes)
{
  auto &buffer = external_output_buffer_ ? *external_output_buffer_ : internal_buffer_;
  if ((buffer.size() + add_bytes) <= max_length_) {
    buffer.resize((buffer.size() + add_bytes));
  } else {
    setstate(std::ios::badbit);
  }
  return bool(*this);
}


template <> void cvm::memory_stream::write_object(std::string const &t)
{
  size_t const string_length = t.size();
  size_t const new_data_size = sizeof(size_t) + sizeof(char) * string_length;
  if (expand_output_buffer(new_data_size)) {
    std::memcpy(output_location(), &string_length, sizeof(size_t));
    incr_write_pos(sizeof(size_t));
    std::memcpy(output_location(), t.c_str(), t.size() * sizeof(char));
    incr_write_pos(t.size() * sizeof(char));
  }
}

template <> cvm::memory_stream &operator<<(cvm::memory_stream &os, std::string const &t)
{
  os.write_object<std::string>(t);
  return os;
}

template <> void cvm::memory_stream::write_object(colvarvalue const &t)
{
  *this << t;
}

template <> void cvm::memory_stream::write_object(cvm::vector1d<cvm::real> const &t)
{
  return write_vector<cvm::real>(t.data_array());
}

template <>
cvm::memory_stream &operator<<(cvm::memory_stream &os, cvm::vector1d<cvm::real> const &t)
{
  os.write_vector<cvm::real>(t.data_array());
  return os;
}


template <> void cvm::memory_stream::read_object(std::string &t)
{
  begin_reading();
  size_t string_length = 0;
  if (has_remaining(sizeof(size_t))) {
    std::memcpy(&string_length, input_location(), sizeof(size_t));
    incr_read_pos(sizeof(size_t));
    if (has_remaining(string_length * sizeof(char))) {
      t.assign(reinterpret_cast<char const *>(input_location()), string_length);
      incr_read_pos(string_length * sizeof(char));
      done_reading();
    } else {
      setstate(std::ios::failbit);
    }
  }
}

template <> cvm::memory_stream &operator>>(cvm::memory_stream &is, std::string &t)
{
  is.read_object<std::string>(t);
  return is;
}

template <> void cvm::memory_stream::read_object(colvarvalue &t)
{
  *this >> t;
}

template <> void cvm::memory_stream::read_object(cvm::vector1d<cvm::real> &t)
{
  return read_vector<cvm::real>(t.data_array());
}

template <> cvm::memory_stream &operator>>(cvm::memory_stream &is, cvm::vector1d<cvm::real> &t)
{
  is.read_vector<cvm::real>(t.data_array());
  return is;
}
