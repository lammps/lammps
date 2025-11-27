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

#ifndef LMP_ERROR_H
#define LMP_ERROR_H

#include "exceptions.h"
#include "pointers.h"

namespace LAMMPS_NS {

/** \class LAMMPS_NS::Error
 * \brief Error handling class for LAMMPS
 *
 * The Error class provides a centralized mechanism for error handling and
 * warning messages throughout LAMMPS. It supports different error levels:
 *
 * - **Fatal errors** that terminate execution (all() and one() methods)
 * - **Warnings** that are reported but allow execution to continue
 * - **Universe-level errors** for multi-partition simulations
 *
 * Error messages include the source file name and line number where the error
 * occurred, and can optionally indicate which command-line argument caused the
 * error. Messages are formatted using the fmtlib library for flexible formatting.
 *
 * The class also tracks the number of warnings issued and supports a maximum
 * warning count to prevent excessive output. The last error message is stored
 * and can be retrieved, which is useful for the library interface and testing.
 *
 * \sa LAMMPS_NS::Pointers, LAMMPS_NS::ErrorType */
class Error : protected Pointers {
 public:
  /** Constructor for Error class
   * \param lmp Pointer to the main LAMMPS instance */
  Error(class LAMMPS *);

  /** Abort all processes in the universe after printing error message
   * \param file Source file name where error occurred
   * \param line Line number in source file
   * \param str Error message string */
  [[noreturn]] void universe_all(const std::string &, int, const std::string &);

  /** Abort this process in the universe after printing error message
   * \param file Source file name where error occurred
   * \param line Line number in source file
   * \param str Error message string */
  [[noreturn]] void universe_one(const std::string &, int, const std::string &);

  /** Print a warning message in the universe
   * \param file Source file name where warning occurred
   * \param line Line number in source file
   * \param str Warning message string */
  void universe_warn(const std::string &, int, const std::string &);

  static constexpr int NOPOINTER = -2;    /**< Sentinel value: no pointer to faulty argument */
  static constexpr int NOLASTLINE = -3;   /**< Sentinel value: no last input line to display */
  static constexpr int ARGZERO = -99;     /**< Sentinel value: error in argument zero (command) */

  // regular error calls

  /** Abort all MPI processes in the world communicator with an error message
   * \param file Source file name where error occurred
   * \param line Line number in source file
   * \param str Error message string */
  [[noreturn]] void all(const std::string &file, int line, const std::string &str)
  {
    all(file, line, NOPOINTER, str);
  }

  /** Abort all MPI processes with a formatted error message
   * \param file Source file name where error occurred
   * \param line Line number in source file
   * \param format Format string (fmtlib style)
   * \param args Arguments for format string */
  template <typename... Args>
  [[noreturn]] void all(const std::string &file, int line, const std::string &format,
                        Args &&...args)
  {
    _all(file, line, NOPOINTER, format, fmt::make_format_args(args...));
  }

  /** Abort this MPI process with an error message
   * \param file Source file name where error occurred
   * \param line Line number in source file
   * \param str Error message string */
  [[noreturn]] void one(const std::string &file, int line, const std::string &str)
  {
    one(file, line, NOPOINTER, str);
  }

  /** Abort this MPI process with a formatted error message
   * \param file Source file name where error occurred
   * \param line Line number in source file
   * \param format Format string (fmtlib style)
   * \param args Arguments for format string */
  template <typename... Args>
  [[noreturn]] void one(const std::string &file, int line, const std::string &format,
                        Args &&...args)
  {
    _one(file, line, NOPOINTER, format, fmt::make_format_args(args...));
  }

  // overloaded error calls indicating faulty argument in command line

  /** Abort all MPI processes with an error message indicating faulty argument
   * \param file Source file name where error occurred
   * \param line Line number in source file
   * \param failed Index of faulty argument (or NOPOINTER, NOLASTLINE, ARGZERO)
   * \param str Error message string */
  [[noreturn]] void all(const std::string &, int, int, const std::string &);

  /** Abort all MPI processes with a formatted error indicating faulty argument
   * \param file Source file name where error occurred
   * \param line Line number in source file
   * \param failed Index of faulty argument
   * \param format Format string (fmtlib style)
   * \param args Arguments for format string */
  template <typename... Args>
  [[noreturn]] void all(const std::string &file, int line, int failed, const std::string &format,
                        Args &&...args)
  {
    _all(file, line, failed, format, fmt::make_format_args(args...));
  }

  /** Abort this MPI process with an error message indicating faulty argument
   * \param file Source file name where error occurred
   * \param line Line number in source file
   * \param failed Index of faulty argument
   * \param str Error message string */
  [[noreturn]] void one(const std::string &, int, int, const std::string &);

  /** Abort this MPI process with a formatted error indicating faulty argument
   * \param file Source file name where error occurred
   * \param line Line number in source file
   * \param failed Index of faulty argument
   * \param format Format string (fmtlib style)
   * \param args Arguments for format string */
  template <typename... Args>
  [[noreturn]] void one(const std::string &file, int line, int failed, const std::string &format,
                        Args &&...args)
  {
    _one(file, line, failed, format, fmt::make_format_args(args...));
  }

  /** Print a warning message
   * \param file Source file name where warning occurred
   * \param line Line number in source file
   * \param str Warning message string */
  void warning(const std::string &, int, const std::string &);

  /** Print a formatted warning message
   * \param file Source file name where warning occurred
   * \param line Line number in source file
   * \param format Format string (fmtlib style)
   * \param args Arguments for format string */
  template <typename... Args>
  void warning(const std::string &file, int line, const std::string &format, Args &&...args)
  {
    _warning(file, line, format, fmt::make_format_args(args...));
  }

  /** Clean up and exit LAMMPS
   * \param status Exit status code (default 0) */
  [[noreturn]] void done(int = 0);    // 1 would be fully backwards compatible

  /** Get the current warning count
   * \return Number of warnings issued */
  int get_numwarn() const { return numwarn; }

  /** Get the maximum warning count
   * \return Maximum number of warnings before suppression */
  int get_maxwarn() const { return maxwarn; }

  /** Set the current warning count
   * \param val New warning count value */
  void set_numwarn(int val) { numwarn = val; }

  /** Set the maximum warning count
   * \param val Maximum number of warnings before suppression */
  void set_maxwarn(int val) { maxwarn = val; }

  /** Set whether all warnings are printed regardless of count
   * \param val 1 to print all warnings, 0 to suppress after maxwarn */
  void set_allwarn(int val) { allwarn = val; }

  /** Get the last error message
   * \return String containing the last error message */
  std::string get_last_error() const;

  /** Get the type of the last error
   * \return ErrorType enum value of the last error */
  ErrorType get_last_error_type() const;

  /** Store an error message for later retrieval
   * \param msg Error message string
   * \param type Type of error (default ERROR_NORMAL) */
  void set_last_error(const std::string &msg, ErrorType type = ERROR_NORMAL);

  /** Control whether error messages are displayed
   * \param flag 1 to show errors, 0 to suppress display
   * \return Previous value of showerror flag */
  int set_show_error(const int flag);

 private:
  std::string last_error_message;    /**< Storage for the last error message */
  ErrorType last_error_type;         /**< Type of the last error */

  int numwarn;      /**< Current count of warnings issued */
  int maxwarn;      /**< Maximum warnings before suppression */
  int allwarn;      /**< Flag to print all warnings regardless of count */
  int showerror;    /**< Flag controlling error message display */

  // internal versions that accept explicit fmtlib arguments
  [[noreturn]] void _all(const std::string &, int, int, fmt::string_view, fmt::format_args args);
  [[noreturn]] void _one(const std::string &, int, int, fmt::string_view, fmt::format_args args);
  void _warning(const std::string &, int, fmt::string_view, fmt::format_args args);
  void _message(const std::string &, int, fmt::string_view, fmt::format_args args);
};

}    // namespace LAMMPS_NS

#endif
