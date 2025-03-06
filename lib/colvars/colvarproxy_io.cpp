// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

// Using access() to check if a file exists (until we can assume C++14/17)
#if !defined(_WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#if defined(_WIN32)
#include <io.h>
#endif

#include <cerrno>
#include <cstdio>

#include <list>
#include <map>
#include <sstream>
#include <fstream>

#include "colvarmodule.h"
#include "colvarproxy_io.h"


colvarproxy_io::colvarproxy_io()
{
  restart_frequency_engine = 0;
  input_stream_error_ = new std::istringstream();
  input_stream_error_->setstate(std::ios::badbit);
  output_stream_error_ = new std::ostringstream();
  output_stream_error_->setstate(std::ios::badbit);
}


colvarproxy_io::~colvarproxy_io()
{
  delete input_stream_error_;
  close_input_streams();
  delete output_stream_error_;
  close_output_streams();
}


bool colvarproxy_io::io_available()
{
  return false;
}


int colvarproxy_io::get_frame(long int&)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_io::set_frame(long int)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_io::backup_file(char const *filename)
{
  // Simplified version of NAMD_file_exists()
  int exit_code;
  do {
#if defined(_WIN32) && !defined(__CYGWIN__)
    // We could use _access_s here, but it is probably too new
    exit_code = _access(filename, 00);
#else
    exit_code = access(filename, F_OK);
#endif
  } while ((exit_code != 0) && (errno == EINTR));
  if (exit_code != 0) {
    if (errno == ENOENT) {
      // File does not exist
      return COLVARS_OK;
    } else {
      return cvm::error("Unknown error while checking if file \""+
                        std::string(filename)+"\" exists.\n", COLVARS_ERROR);
    }
  }

  // The file exists, then rename it
  if (std::string(filename).rfind(std::string(".colvars.state")) !=
      std::string::npos) {
    return rename_file(filename, (std::string(filename)+".old").c_str());
  } else {
    return rename_file(filename, (std::string(filename)+".BAK").c_str());
  }
}


int colvarproxy_io::remove_file(char const *filename)
{
  int error_code = COLVARS_OK;
#if defined(_WIN32) && !defined(__CYGWIN__)
  // Because the file may be open by other processes, rename it to filename.old
  std::string const renamed_file(std::string(filename)+".old");
  // It may still be there from an interrupted run, so remove it to be safe
  std::remove(renamed_file.c_str());
  int rename_exit_code = 0;
  while ((rename_exit_code = std::rename(filename,
                                         renamed_file.c_str())) != 0) {
    if (errno == EINTR) continue;
    error_code |= COLVARS_FILE_ERROR;
    break;
  }
  // Ask to remove filename.old, but ignore any errors raised
  std::remove(renamed_file.c_str());
#else
  if (std::remove(filename)) {
    if (errno != ENOENT) {
      error_code |= COLVARS_FILE_ERROR;
    }
  }
#endif
  if (error_code != COLVARS_OK) {
    return cvm::error("Error: in removing file \""+std::string(filename)+
                      "\".\n.",
                      error_code);
  }
  return COLVARS_OK;
}


int colvarproxy_io::rename_file(char const *filename, char const *newfilename)
{
  int error_code = COLVARS_OK;
#if defined(_WIN32) && !defined(__CYGWIN__)
  // On straight Windows, must remove the destination before renaming it
  if (_access(newfilename, 00) == 0) {
    error_code |= remove_file(newfilename);
  }
#endif
  int rename_exit_code = 0;
  while ((rename_exit_code = std::rename(filename, newfilename)) != 0) {
    if (errno == EINTR) continue;
    // Call log() instead of error to allow the next try
    cvm::log("Error: in renaming file \""+std::string(filename)+"\" to \""+
             std::string(newfilename)+"\".\n.");
    error_code |= COLVARS_FILE_ERROR;
    if (errno == EXDEV) continue;
    break;
  }
  return rename_exit_code ? error_code : COLVARS_OK;
}


int colvarproxy_io::set_input_prefix(std::string const &prefix)
{
  // set input restart name and strip the extension, if present
  input_prefix_str = prefix;
  if (input_prefix_str.rfind(".colvars.state") != std::string::npos) {
    input_prefix_str.erase(input_prefix_str.rfind(".colvars.state"),
                           std::string(".colvars.state").size());
  }
  return COLVARS_OK;
}


int colvarproxy_io::set_output_prefix(std::string const &prefix)
{
  // set input restart name and strip the extension, if present
  output_prefix_str = prefix;
  if (output_prefix_str.rfind(".colvars.state") != std::string::npos) {
    output_prefix_str.erase(output_prefix_str.rfind(".colvars.state"),
                            std::string(".colvars.state").size());
  }
  return COLVARS_OK;
}


int colvarproxy_io::set_restart_output_prefix(std::string const &prefix)
{
  // set input restart name and strip the extension, if present
  restart_output_prefix_str = prefix;
  if (restart_output_prefix_str.rfind(".colvars.state") != std::string::npos) {
    restart_output_prefix_str.erase(restart_output_prefix_str.rfind(".colvars.state"),
                                    std::string(".colvars.state").size());
  }
  return COLVARS_OK;
}


int colvarproxy_io::set_default_restart_frequency(int freq)
{
  // TODO check for compatibility with colvarsRestartFrequency
  restart_frequency_engine = freq;
  return COLVARS_OK;
}



std::istream &colvarproxy_io::input_stream(std::string const &input_name,
                                           std::string const description,
                                           bool error_on_fail)
{
  if (!io_available()) {
    cvm::error("Error: trying to access an input file/channel "
               "from the wrong thread.\n", COLVARS_BUG_ERROR);
    return *input_stream_error_;
  }

  if (colvarproxy_io::input_stream_exists(input_name)) {
    std::ifstream *ifs =
      dynamic_cast<std::ifstream *>(input_streams_[input_name]);
    if (ifs && !ifs->is_open()) {
      // This file was opened before, re-open it.  Using std::ios::binary to
      // work around differences in line termination conventions
      // See https://github.com/Colvars/colvars/commit/8236879f7de4
      ifs->open(input_name.c_str(), std::ios::binary);
    }
  } else {
    input_streams_[input_name] = new std::ifstream(input_name.c_str(),
                                                   std::ios::binary);
  }

  if (input_streams_[input_name]->fail() && error_on_fail) {
    cvm::error("Error: cannot open "+description+" \""+input_name+"\".\n",
               COLVARS_FILE_ERROR);
  }

  return *(input_streams_[input_name]);
}


std::istream &
colvarproxy_io::input_stream_from_string(std::string const &input_name,
                                         std::string const &content,
                                         std::string const description)
{
  if (!io_available()) {
    cvm::error("Error: trying to access an input file/channel "
               "from the wrong thread.\n", COLVARS_BUG_ERROR);
    return *input_stream_error_;
  }

  if (colvarproxy_io::input_stream_exists(input_name)) {

    std::istringstream *iss =
      dynamic_cast<std::istringstream *>(input_streams_[input_name]);
    if (iss) {
      // If there is already a stringstream, replace it
      delete iss;
    } else {
      std::ifstream *ifs =
        dynamic_cast<std::ifstream *>(input_streams_[input_name]);
      if (ifs) {
        if (ifs->is_open()) {
          ifs->close();
        }
      }
    }
  }

  input_streams_[input_name] = new std::istringstream(content);

  return *(input_streams_[input_name]);
}


bool colvarproxy_io::input_stream_exists(std::string const &input_name)
{
  return (input_streams_.count(input_name) > 0);
}


int colvarproxy_io::close_input_stream(std::string const &input_name)
{
  if (colvarproxy_io::input_stream_exists(input_name)) {
    std::ifstream *ifs = dynamic_cast<std::ifstream *>(input_streams_[input_name]);
    if (ifs) {
      if (ifs->is_open()) {
        ifs->close();
      }
    } else {
      // From a string, just rewind to the begining
      std::istringstream * iss = dynamic_cast<std::istringstream *>(input_streams_[input_name]);
      if (iss) {
        iss->clear();
        iss->seekg(0);
      }
    }
    return COLVARS_OK;
  }
  return cvm::error("Error: input file/channel \""+input_name+
                    "\" does not exist.\n", COLVARS_FILE_ERROR);
}


int colvarproxy_io::delete_input_stream(std::string const &input_name)
{
  if (colvarproxy_io::close_input_stream(input_name) == COLVARS_OK) {
    std::ifstream *ifs = dynamic_cast<std::ifstream *>(input_streams_[input_name]);
    if (ifs) {
      delete ifs;
    } else {
      std::istringstream * iss = dynamic_cast<std::istringstream *>(input_streams_[input_name]);
      if (iss) {
        delete iss;
      }
    }
    input_streams_.erase(input_name);
    return COLVARS_OK;
  }
  return cvm::error("Error: input file/channel \""+input_name+
                    "\" does not exist.\n", COLVARS_FILE_ERROR);
}


int colvarproxy_io::close_input_streams()
{
  for (std::map<std::string,
         std::istream *>::iterator ii = input_streams_.begin();
       ii != input_streams_.end();
       ii++) {
    delete ii->second;
  }
  input_streams_.clear();
  return COLVARS_OK;
}


std::list<std::string> colvarproxy_io::list_input_stream_names() const
{
  std::list<std::string> result;
  for (std::map<std::string,
         std::istream *>::const_iterator ii = input_streams_.begin();
       ii != input_streams_.end();
       ii++) {
    result.push_back(ii->first);
  }
  return result;
}


std::ostream & colvarproxy_io::output_stream(std::string const &output_name,
                                             std::string const description)
{
  if (cvm::debug()) {
    cvm::log("Using colvarproxy_io::output_stream()\n");
  }

  if (!io_available()) {
    cvm::error("Error: trying to access an output file/channel "
               "from the wrong thread.\n", COLVARS_BUG_ERROR);
    return *output_stream_error_;
  }

  if (colvarproxy_io::output_stream_exists(output_name)) {
    return *(output_streams_[output_name]);
  }

  backup_file(output_name.c_str());

  output_streams_[output_name] = new std::ofstream(output_name.c_str(), std::ios::binary);
  if (!*(output_streams_[output_name])) {
    cvm::error("Error: cannot write to "+description+" \""+output_name+"\".\n",
               COLVARS_FILE_ERROR);
  }

  return *(output_streams_[output_name]);
}


bool colvarproxy_io::output_stream_exists(std::string const &output_name)
{
  return (output_streams_.count(output_name) > 0);
}


int colvarproxy_io::flush_output_stream(std::string const &output_name)
{
  if (!io_available()) {
    // No-op
    return COLVARS_OK;
  }

  if (colvarproxy_io::output_stream_exists(output_name)) {
    (dynamic_cast<std::ofstream *>(output_streams_[output_name]))->flush();
    return COLVARS_OK;
  }

  return COLVARS_OK;
}


int colvarproxy_io::flush_output_streams()
{
  if (!io_available()) {
    return COLVARS_OK;
  }

  for (std::map<std::string, std::ostream *>::iterator osi = output_streams_.begin();
       osi != output_streams_.end();
       osi++) {
    (dynamic_cast<std::ofstream *>(osi->second))->flush();
  }

  return COLVARS_OK;
}


int colvarproxy_io::close_output_stream(std::string const &output_name)
{
  if (!io_available()) {
    return cvm::error("Error: trying to access an output file/channel "
                      "from the wrong thread.\n", COLVARS_BUG_ERROR);
  }

  if (colvarproxy_io::output_stream_exists(output_name)) {
    (dynamic_cast<std::ofstream *>(output_streams_[output_name]))->close();
    delete output_streams_[output_name];
    output_streams_.erase(output_name);
  }

  return COLVARS_OK;
}


int colvarproxy_io::close_output_streams()
{
  if (! io_available()) {
    return COLVARS_OK;
  }

  for (std::map<std::string, std::ostream *>::iterator osi = output_streams_.begin();
       osi != output_streams_.end();
       osi++) {
    (dynamic_cast<std::ofstream *>(osi->second))->close();
    delete osi->second;
  }
  output_streams_.clear();

  return COLVARS_OK;
}
