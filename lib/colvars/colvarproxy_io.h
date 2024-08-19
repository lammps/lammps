// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_IO_H
#define COLVARPROXY_IO_H

#include <iosfwd>
#include <list>
#include <map>
#include <string>


/// Methods for data input/output
class colvarproxy_io {

public:

  /// Constructor
  colvarproxy_io();

  /// Destructor
  virtual ~colvarproxy_io();

  /// Ensure that we're on the main thread (derived class will do actual check)
  virtual bool io_available();

  /// \brief Save the current frame number in the argument given
  // Returns error code
  virtual int get_frame(long int &);

  /// \brief Set the current frame number (as well as colvarmodule::it)
  // Returns error code
  virtual int set_frame(long int);

  /// \brief Rename the given file, before overwriting it
  virtual int backup_file(char const *filename);

  /// \brief Rename the given file, before overwriting it
  inline int backup_file(std::string const &filename)
  {
    return backup_file(filename.c_str());
  }

  /// Remove the given file (on Windows only, rename to filename.old)
  virtual int remove_file(char const *filename);

  /// Remove the given file (on Windows only, rename to filename.old)
  inline int remove_file(std::string const &filename)
  {
    return remove_file(filename.c_str());
  }

  /// Rename the given file
  virtual int rename_file(char const *filename, char const *newfilename);

  /// Rename the given file
  inline int rename_file(std::string const &filename,
                         std::string const &newfilename)
  {
    return rename_file(filename.c_str(), newfilename.c_str());
  }

  /// Prefix of the input state file to be read next
  inline std::string const & input_prefix() const
  {
    return input_prefix_str;
  }

  /// Initialize input_prefix (NOTE: it will be erased after state file is read)
  virtual int set_input_prefix(std::string const &prefix);

  /// Default prefix to be used for all output files (final configuration)
  inline std::string const & output_prefix() const
  {
    return output_prefix_str;
  }

  /// Set default output prefix
  virtual int set_output_prefix(std::string const &prefix);

  /// Prefix of the restart (checkpoint) file to be written next
  inline std::string const & restart_output_prefix() const
  {
    return restart_output_prefix_str;
  }

  /// Set default restart state file prefix
  virtual int set_restart_output_prefix(std::string const &prefix);

  /// Default restart frequency (as set by the simulation engine)
  inline int default_restart_frequency() const
  {
    return restart_frequency_engine;
  }

  /// Communicate/set the restart frequency of the simulation engine
  virtual int set_default_restart_frequency(int freq);

  // The input stream functions below are not virtual, because they currently
  // rely on the fact that either std::ifstream or std::istringstream is used

  /// Returns a reference to given input stream, creating it if needed
  /// \param input_name File name (later only a handle)
  /// \param description Purpose of the file
  /// \param error_on_fail Raise error when failing to open (allow testing)
  std::istream & input_stream(std::string const &input_name,
                              std::string const description = "file/channel",
                              bool error_on_fail = true);

  /// Returns a reference to given input stream, creating it if needed
  /// \param input_name Identifier of the input stream
  /// \param content Set this string as the stream's buffer
  /// \param description Purpose of the stream
  std::istream & input_stream_from_string(std::string const &input_name,
                                          std::string const &content,
                                          std::string const description = "string");

  /// Check if the file/channel is open (without opening it if not)
  bool input_stream_exists(std::string const &input_name);

  /// Closes the given input stream
  int close_input_stream(std::string const &input_name);

  /// Closes all input streams
  int close_input_streams();

  /// Same as close_input_stream(), but also removes the corresponding entry from memory
  int delete_input_stream(std::string const &input_name);

  /// List all input streams that were opened at some point
  std::list<std::string> list_input_stream_names() const;

  /// Returns a reference to the named output file/channel (open it if needed)
  /// \param output_name File name or identifier
  /// \param description Purpose of the file
  virtual std::ostream &output_stream(std::string const &output_name,
                                      std::string const description);

  /// Check if the file/channel is open (without opening it if not)
  virtual bool output_stream_exists(std::string const &output_name);

  /// Flushes the given output file/channel
  virtual int flush_output_stream(std::string const &output_name);

  /// Flushes all output files/channels
  virtual int flush_output_streams();

  /// Closes the given output file/channel
  virtual int close_output_stream(std::string const &output_name);

  /// Close all open files/channels to prevent data loss
  virtual int close_output_streams();

protected:

  /// Prefix of the input state file to be read next
  std::string input_prefix_str;

  /// Default prefix to be used for all output files (final configuration)
  std::string output_prefix_str;

  /// Prefix of the restart (checkpoint) file to be written next
  std::string restart_output_prefix_str;

  /// How often the simulation engine will write its own restart
  int restart_frequency_engine;

  /// Container of input files/channels indexed by path name
  std::map<std::string, std::istream *> input_streams_;

  /// Object whose reference is returned when read errors occur
  std::istream *input_stream_error_;

  /// Currently open output files/channels
  std::map<std::string, std::ostream *> output_streams_;

  /// Object whose reference is returned when write errors occur
  std::ostream *output_stream_error_;
};


#endif
