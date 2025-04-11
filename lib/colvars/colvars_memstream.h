// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef MEMORY_STREAM_H
#define MEMORY_STREAM_H

#include <cstring>
#include <string>
#include <typeinfo>
#include <vector>
#include <iomanip>


// Work around missing std::is_trivially_copyable in old GCC and Clang versions
// TODO remove this after CentOS 7 has been beyond EOL for a while
#if (defined(__GNUC__) && (__GNUC__ < 5) && !defined(__clang__)) || (defined(__clang__) && (__clang_major__ < 7))
// Clang needs an exception, because it defines __GNUC__ as well
#define IS_TRIVIALLY_COPYABLE(T) __has_trivial_copy(T)
#else
#define IS_TRIVIALLY_COPYABLE(T) std::is_trivially_copyable<T>::value
#endif


class cvm::memory_stream {

public:

  /// Set up an empty stream with an internal buffer, suitable for writing to
  /// \param max_length Maximum allowed capacity (default is 64 GiB)
  memory_stream(size_t max_length = (static_cast<size_t>(1L) << 36)) : max_length_(max_length) {}

  /// Set up a stream based on an external input buffer
  memory_stream(size_t n, unsigned char const *buf)
    : external_input_buffer_(buf), internal_buffer_(), data_length_(n), max_length_(data_length_)
  {
  }

  /// Set up a stream based on an external output buffer
  memory_stream(std::vector<unsigned char> &buf) : memory_stream()
  {
    external_output_buffer_ = &buf;
  }

  /// Length of the buffer
  inline size_t length() const { return data_length_; }

  /// Output buffer
  inline unsigned char *output_buffer()
  {
    return (external_output_buffer_ ? external_output_buffer_->data() : internal_buffer_.data());
  }

  /// Next location to write to
  inline unsigned char *output_location() { return output_buffer() + data_length_; }

  /// Input buffer
  inline unsigned char const *input_buffer() const
  {
    return (external_input_buffer_ ? external_input_buffer_ : internal_buffer_.data());
  }

  /// Next location to read from
  inline unsigned char const *input_location() const { return input_buffer() + read_pos_; }

  /// Cast operator to be used to test for errors
  inline explicit operator bool() const { return state_ == std::ios::goodbit; }

  /// Write a simple object to the output buffer
  template <typename T> void write_object(T const &t);

  /// Wrapper to write_object()
  template <typename T> friend memory_stream &operator<<(memory_stream &os, T const &t);

  /// Write a vector of simple objects to the output buffer
  template <typename T> void write_vector(std::vector<T> const &t);

  /// Wrapper to write_vector()
  template <typename T>
  friend memory_stream &operator<<(memory_stream &os, std::vector<T> const &t);

  /// Read a simple object from the buffer
  template <typename T> void read_object(T &t);

  /// Wrapper to read_object()
  template <typename T> friend memory_stream &operator>>(memory_stream &is, T &t);

  /// Read a vector of simple objects from the buffer
  template <typename T> void read_vector(std::vector<T> &t);

  /// Wrapper to read_vector()
  template <typename T> friend memory_stream &operator>>(memory_stream &is, std::vector<T> &t);


  // Compatibility with STL stream functions

  /// Report the current position in the buffer
  inline size_t tellg() const { return read_pos_; }

  /// Report the current position in the buffer
  inline memory_stream & seekg(size_t pos) { read_pos_ = pos; return *this; }

  /// Ignore formatting operators
  inline void setf(decltype(std::ios::fmtflags(0)), decltype(std::ios::floatfield)) {}

  /// Ignore formatting operators
  inline void flags(decltype(std::ios::fmtflags(0))) {}

  /// Get the current formatting flags (i.e. none because this stream is unformatted)
  inline decltype(std::ios::fmtflags(0)) flags() const { return std::ios::fmtflags(0); }

  /// Get the error code
  inline std::ios::iostate rdstate() const { return state_; }

  /// Set the error code
  inline void setstate(std::ios::iostate new_state) { state_ |= new_state; }

  /// Clear the error code
  inline void clear() { state_ = std::ios::goodbit; }

protected:

  /// External output buffer
  std::vector<unsigned char> *external_output_buffer_ = nullptr;

  /// External input buffer
  unsigned char const *external_input_buffer_ = nullptr;

  /// Internal buffer (may server for both input and output)
  std::vector<unsigned char> internal_buffer_;

  /// Length of the data buffer (either internal or external)
  size_t data_length_ = 0L;

  /// Largest allowed capacity of the data buffer
  size_t const max_length_;

  /// Error status
  std::ios::iostate state_ = std::ios::goodbit;

  /// Add the requester number of bytes to the array capacity; return false if buffer is external
  bool expand_output_buffer(size_t add_bytes);

  /// Move the buffer position past the data just written
  inline void incr_write_pos(size_t c) { data_length_ += c; }

  /// Current position when reading from the buffer
  size_t read_pos_ = 0L;

  /// Begin an attempt to read an object; assume EOF unless there is space remaining
  inline void begin_reading() { setstate(std::ios::eofbit); }

  /// Mark the reading attempt succesful
  inline void done_reading() { clear(); }

  /// Move the buffer position past the data just read
  inline void incr_read_pos(size_t c) { read_pos_ += c; }

  /// Check that the buffer contains enough bytes to read as the argument says, set error
  /// otherwise
  inline bool has_remaining(size_t c) { return c <= (data_length_ - read_pos_); }
  };

template <typename T> void cvm::memory_stream::write_object(T const &t)
{
  static_assert(IS_TRIVIALLY_COPYABLE(T), "Cannot use write_object() on complex type");
  size_t const new_data_size = sizeof(T);
  if (expand_output_buffer(new_data_size)) {
    std::memcpy(output_location(), &t, sizeof(T));
    incr_write_pos(new_data_size);
  }
}

template <typename T> cvm::memory_stream &operator<<(cvm::memory_stream &os, T const &t)
{
  os.write_object<T>(t);
  return os;
}

template <typename T> void cvm::memory_stream::write_vector(std::vector<T> const &t)
{
  static_assert(IS_TRIVIALLY_COPYABLE(T), "Cannot use write_vector() on complex type");
  size_t const vector_length = t.size();
  size_t const new_data_size = sizeof(size_t) + sizeof(T) * vector_length;
  if (expand_output_buffer(new_data_size)) {
    std::memcpy(output_location(), &vector_length, sizeof(size_t));
    incr_write_pos(sizeof(T));
    std::memcpy(output_location(), t.data(), t.size() * sizeof(T));
    incr_write_pos(t.size() * sizeof(T));
  }
}

template <typename T>
cvm::memory_stream &operator<<(cvm::memory_stream &os, std::vector<T> const &t)
{
  os.write_vector<T>(t);
  return os;
}

template <typename T> void cvm::memory_stream::read_object(T &t)
{
  static_assert(IS_TRIVIALLY_COPYABLE(T), "Cannot use read_object() on complex type");
  begin_reading();
  if (has_remaining(sizeof(T))) {
    std::memcpy(&t, input_location(), sizeof(T));
    incr_read_pos(sizeof(T));
    done_reading();
  }
}

template <typename T> cvm::memory_stream &operator>>(cvm::memory_stream &is, T &t)
{
  is.read_object<T>(t);
  return is;
}

template <typename T> void cvm::memory_stream::read_vector(std::vector<T> &t)
{
  static_assert(IS_TRIVIALLY_COPYABLE(T), "Cannot use read_vector() on complex type");
  begin_reading();
  size_t vector_length = 0;
  if (has_remaining(sizeof(size_t))) {
    std::memcpy(&vector_length, input_location(), sizeof(size_t));
    incr_read_pos(sizeof(size_t));
    if (has_remaining(vector_length * sizeof(T))) {
      t.resize(vector_length);
      std::memcpy(t.data(), input_location(), vector_length * sizeof(T));
      incr_read_pos(vector_length * sizeof(T));
      done_reading();
    } else {
      setstate(std::ios::failbit);
    }
  }
}

template <typename T> cvm::memory_stream &operator>>(cvm::memory_stream &is, std::vector<T> &t)
{
  is.read_vector<T>(t);
  return is;
}

template <typename T> cvm::memory_stream &operator<<(cvm::memory_stream &os,
                                                     decltype(std::setprecision(10)) const &)
{
  return os;
}

#if !defined(_MSC_VER) && !defined(__SUNPRO_CC)
// Visual Studio and MSVC use the same return type for both modifiers
template <typename T> cvm::memory_stream &operator<<(cvm::memory_stream &os,
                                                     decltype(std::setw(10)) const &)
{
  return os;
}
#endif

// Declare specializations

template <> void cvm::memory_stream::write_object(std::string const &t);

template <> cvm::memory_stream &operator<<(cvm::memory_stream &os, std::string const &t);

template <> void cvm::memory_stream::write_object(colvarvalue const &t);

template <> cvm::memory_stream &operator<<(cvm::memory_stream &os, colvarvalue const &x);

template <> void cvm::memory_stream::write_object(cvm::vector1d<cvm::real> const &t);

template <>
cvm::memory_stream &operator<<(cvm::memory_stream &os, cvm::vector1d<cvm::real> const &t);

template <> void cvm::memory_stream::read_object(std::string &t);

template <> cvm::memory_stream &operator>>(cvm::memory_stream &is, std::string &t);

template <> void cvm::memory_stream::read_object(colvarvalue &t);

template <> cvm::memory_stream &operator>>(cvm::memory_stream &is, colvarvalue &t);

template <> void cvm::memory_stream::read_object(cvm::vector1d<cvm::real> &t);

template <> cvm::memory_stream &operator>>(cvm::memory_stream &is, cvm::vector1d<cvm::real> &t);

#endif
