// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARSCRIPT_H
#define COLVARSCRIPT_H

#include <string>
#include <vector>
#include <map>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarproxy.h"


// Only these error values are part of the scripting interface
#define COLVARSCRIPT_ERROR -1
#define COLVARSCRIPT_OK 0


class colvardeps;

class colvarscript  {

private:

  colvarproxy *proxy_;
  colvarmodule *colvars;

  inline colvarscript() {} // no-argument construction forbidden

public:

  friend class colvarproxy;

  colvarscript(colvarproxy *p, colvarmodule *m);

  ~colvarscript();

  /// String representation of the result of a script call
  std::string str_result_;

  /// Run a script command with space-separated positional arguments (objects)
  int run(int objc, unsigned char *const objv[]);

  /// Get the string result of the current scripting call
  inline std::string const &str_result() const
  {
    return str_result_;
  }

  /// Modify the string result of the current scripting call
  inline std::string &modify_str_result()
  {
    return str_result_;
  }

  /// Set the return value to the given string
  int set_result_str(std::string const &s);

  /// Clear the string result
  int clear_str_result();

  /// Add the given string to the error message of the script interface
  void add_error_msg(std::string const &s);

  /// Commands available
  enum command {
#define CVSCRIPT_ENUM_COMM(COMM) COMM,
#undef CVSCRIPT
#define CVSCRIPT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)  \
  CVSCRIPT_ENUM_COMM(COMM)
#ifdef COLVARSCRIPT_COMMANDS_H
#undef COLVARSCRIPT_COMMANDS_H
#endif
#include "colvarscript_commands.h"
#undef COLVARSCRIPT_COMMANDS_H
#undef CVSCRIPT
#undef CVSCRIPT_ENUM_COMM
    cv_n_commands
  };

  /// Type of object handling a script command
  enum Object_type {
    use_module,
    use_colvar,
    use_bias
  };

  /// Return the prefix of the individual command for each object function
  std::string get_cmd_prefix(Object_type t);

  /// Get a pointer to the i-th argument of the command (NULL if not given)
  template<Object_type T>
  unsigned char *get_cmd_arg(int iarg, int objc, unsigned char *const objv[]);

  /// Instantiation of get_cmd_arg<> for module-level commands
  unsigned char *get_module_cmd_arg(int iarg, int objc,
                                    unsigned char *const objv[]);

  /// Instantiation of get_cmd_arg<> for colvar-level commands
  unsigned char *get_colvar_cmd_arg(int iarg, int objc,
                                    unsigned char *const objv[]);

  /// Instantiation of get_cmd_arg<> for bias-level commands
  unsigned char *get_bias_cmd_arg(int iarg, int objc,
                                  unsigned char *const objv[]);

  /// Check the argument count of the command
  template<Object_type T>
  int check_cmd_nargs(char const *cmd, int objc,
                      int n_args_min, int n_args_max);

  /// Instantiation of check_cmd_nargs<> for module-level commands
  int check_module_cmd_nargs(char const *cmd, int objc,
                             int n_args_min, int n_args_max);

  /// Instantiation of check_cmd_nargs<> for colvar-level commands
  int check_colvar_cmd_nargs(char const *cmd, int objc,
                             int n_args_min, int n_args_max);

  /// Instantiation of get_cmd_arg<> for bias-level commands
  int check_bias_cmd_nargs(char const *cmd, int objc,
                           int n_args_min, int n_args_max);

  /// Number of positional arguments to shift for each object type
  template<colvarscript::Object_type T>
  int cmd_arg_shift();

  /// Get names of all commands
  inline char const **get_command_names() const
  {
    return cmd_names;
  }

  /// Get one-line help summary for a command
  /// \param cmd Name of the command's function (e.g. "cv_units")
  char const *get_command_help(char const *cmd);

  /// Get description of the return value of a command
  /// \param cmd Name of the command's function (e.g. "cv_units")
  char const *get_command_rethelp(char const *cmd);

  /// Get description of the argument of a command (excluding prefix)
  /// \param cmd Name of the command's function (e.g. "cv_units")
  /// \param i Index of the argument; 0 is the first argument after the
  /// prefix, e.g. "value" has an index of 0 in the array of arguments:
  /// { "cv", "colvar", "xi", "value" }
  char const *get_command_arghelp(char const *cmd, int i);

  /// Get number of required arguments (excluding prefix)
  /// \param cmd Name of the command's function (e.g. "cv_units")
  int get_command_n_args_min(char const *cmd);

  /// Get number of total arguments (excluding prefix)
  /// \param cmd Name of the command's function (e.g. "cv_units")
  int get_command_n_args_max(char const *cmd);

  /// Set the main command for the CLI, when it is not "cv" (e.g. LAMMPS)
  inline void set_cmdline_main_cmd(std::string const &cmd) {
    cmdline_main_cmd_ = cmd;
  }

  /// Get help string for a command (does not specify how it is launched)
  /// \param cmd Name of the command's function (e.g. "cv_units")
  char const *get_command_full_help(char const *cmd);

  /// Get summary of command line syntax for all commands of a given context
  /// \param t One of use_module, use_colvar or use_bias
  std::string get_cmdline_help_summary(Object_type t);

  /// Get a description of how the command should be used in a command line
  /// \param t One of use_module, use_colvar or use_bias
  /// \param c Value of the \link command \endlink enum
  std::string get_command_cmdline_syntax(Object_type t, command c);

  /// Get the command line syntax following by the help string
  /// \param t One of use_module, use_colvar or use_bias
  /// \param cmd Name of the subcommand (e.g. "units")
  std::string get_command_cmdline_help(Object_type t, std::string const &cmd);

  /// Set error code for unsupported script operation
  int unsupported_op();

  /// Pointer to the Colvars main object
  inline colvarmodule *module()
  {
    return this->colvars;
  }

  /// Pointer to the colvarproxy object (interface with host engine)
  inline colvarproxy *proxy()
  {
    return this->proxy_;
  }

  // Input functions - get the string reps of script argument objects

  /// Get the string representation of an object (by default, a simple cast)
  char *obj_to_str(unsigned char *obj);

  /// Get a list of strings from an object (does not work with a simple cast)
  std::vector<std::string> obj_to_str_vector(unsigned char *obj);


  // Output functions - convert internal objects to representations suitable
  // for use in the scripting language.  At the moment only conversion to C
  // strings is supported, and obj is assumed to be a char * pointer.

  /// Copy x into obj if not NULL, or into the script object's result otherwise
  int set_result_int(int const &x, unsigned char *obj = NULL);

  /// Copy x into obj if not NULL, or into the script object's result otherwise
  int set_result_int_vec(std::vector<int> const &x, unsigned char *obj = NULL);

  /// Copy x into obj if not NULL, or into the script object's result otherwise
  int set_result_long_int(long int const &x, unsigned char *obj = NULL);

  /// Copy x into obj if not NULL, or into the script object's result otherwise
  int set_result_long_int_vec(std::vector<long int> const &x,
                              unsigned char *obj = NULL);

  /// Copy x into obj if not NULL, or into the script object's result otherwise
  int set_result_real(cvm::real const &x, unsigned char *obj = NULL);

  /// Copy x into obj if not NULL, or into the script object's result otherwise
  int set_result_real_vec(std::vector<cvm::real> const &x,
                          unsigned char *obj = NULL);

  /// Copy x into obj if not NULL, or into the script object's result otherwise
  int set_result_rvector(cvm::rvector const &x, unsigned char *obj = NULL);

  /// Copy x into obj if not NULL, or into the script object's result otherwise
  int set_result_rvector_vec(std::vector<cvm::rvector> const &x,
                             unsigned char *obj = NULL);

  /// Copy x into obj if not NULL, or into the script object's result otherwise
  int set_result_colvarvalue(colvarvalue const &x, unsigned char *obj = NULL);

  /// Copy x into obj if not NULL, or into the script object's result otherwise
  int set_result_colvarvalue_vec(std::vector<colvarvalue> const &x,
                                 unsigned char *obj = NULL);

private:

  /// Set up all script API functions
  int init_commands();

  /// Set up a single script API function
  int init_command(colvarscript::command const &comm,
                   char const *name, char const *help,
                   int n_args_min, int n_args_max, char const *arghelp,
                   int (*fn)(void *, int, unsigned char * const *));

public: // TODO this function will be removed soon

  /// Run subcommands on base colvardeps object (colvar, bias, ...)
  int proc_features(colvardeps *obj,
                    int argc, unsigned char *const argv[]);

private: // TODO

  /// Internal identifiers of command strings
  std::map<std::string, command> cmd_str_map;

  /// Main command used in command line ("cv" by default)
  std::string cmdline_main_cmd_;

  /// Inverse of cmd_str_map (to be exported outside this class)
  char const **cmd_names;

  /// Help strings for each command
  std::vector<std::string> cmd_help;

  /// Description of the return values of each command (may be empty)
  std::vector<std::string> cmd_rethelp;

  /// Minimum number of arguments for each command
  std::vector<size_t> cmd_n_args_min;

  /// Maximum number of arguments for each command
  std::vector<size_t> cmd_n_args_max;

  /// Help strings for each command argument
  std::vector< std::vector<std::string> > cmd_arghelp;

  /// Full help strings for each command
  std::vector<std::string> cmd_full_help;

  /// Implementations of each command
  std::vector<int (*)(void *, int, unsigned char * const *)> cmd_fns;

  /// Get a pointer to the implementation of the given command
  inline int (*get_cmd_fn(std::string const &cmd_key))(void *,
                                                       int,
                                                       unsigned char * const *)
  {
    if (cmd_str_map.count(cmd_key) > 0) {
      return cmd_fns[cmd_str_map[cmd_key]];
    }
    return NULL;
  }

  /// Set obj equal to x, using its string representation
  template <typename T>
  int set_result_text(T const &x, unsigned char *obj);

  /// Code reused by instances of set_result_text()
  template <typename T>
  int pack_vector_elements_text(std::vector<T> const &x, std::string &x_str);

  /// Code reused by all instances of set_result_text()
  int set_result_text_from_str(std::string const &x_str, unsigned char *obj);


};


/// Get a pointer to the main colvarscript object
inline static colvarscript *colvarscript_obj()
{
  return cvm::main()->proxy->script;
}


/// Get a pointer to the colvar object pointed to by pobj
inline static colvar *colvar_obj(void *pobj)
{
  return reinterpret_cast<colvar *>(pobj);
}


/// Get a pointer to the colvarbias object pointed to by pobj
inline static colvarbias *colvarbias_obj(void *pobj)
{
  return reinterpret_cast<colvarbias *>(pobj);
}



template<colvarscript::Object_type T>
unsigned char *colvarscript::get_cmd_arg(int iarg,
                                         int objc,
                                         unsigned char *const objv[])
{
  int const shift = cmd_arg_shift<T>();
  return (shift+iarg < objc) ? objv[shift+iarg] : NULL;
}


inline unsigned char *colvarscript::get_module_cmd_arg(int iarg, int objc,
                                                       unsigned char *const objv[])
{
  return get_cmd_arg<use_module>(iarg, objc, objv);
}


inline unsigned char *colvarscript::get_colvar_cmd_arg(int iarg, int objc,
                                                       unsigned char *const objv[])
{
  return get_cmd_arg<use_colvar>(iarg, objc, objv);
}


inline unsigned char *colvarscript::get_bias_cmd_arg(int iarg, int objc,
                                                     unsigned char *const objv[])
{
  return get_cmd_arg<use_bias>(iarg, objc, objv);
}


template<colvarscript::Object_type T>
int colvarscript::check_cmd_nargs(char const *cmd,
                                  int objc,
                                  int n_args_min,
                                  int n_args_max)
{
  int const shift = cmd_arg_shift<T>();
  if (objc < shift+n_args_min) {
    add_error_msg("Insufficient number of arguments ("+cvm::to_str(objc)+
                  ") for script function \""+std::string(cmd)+
                  "\":\n"+get_command_full_help(cmd));
    return COLVARSCRIPT_ERROR;
  }
  if (objc > shift+n_args_max) {
    add_error_msg("Too many arguments ("+cvm::to_str(objc)+
                  ") for script function \""+std::string(cmd)+
                  "\":\n"+get_command_full_help(cmd));
    return COLVARSCRIPT_ERROR;
  }
  return COLVARSCRIPT_OK;
}


inline int colvarscript::check_module_cmd_nargs(char const *cmd,
                                                int objc,
                                                int n_args_min,
                                                int n_args_max)
{
  return check_cmd_nargs<use_module>(cmd, objc, n_args_min, n_args_max);
}


inline int colvarscript::check_colvar_cmd_nargs(char const *cmd,
                                                int objc,
                                                int n_args_min,
                                                int n_args_max)
{
  return check_cmd_nargs<use_colvar>(cmd, objc, n_args_min, n_args_max);
}


inline int colvarscript::check_bias_cmd_nargs(char const *cmd,
                                              int objc,
                                              int n_args_min,
                                              int n_args_max)
{
  return check_cmd_nargs<use_bias>(cmd, objc, n_args_min, n_args_max);
}


template<colvarscript::Object_type T>
int colvarscript::cmd_arg_shift()
{
  int shift = 0;
  if (T == use_module) {
    // "cv" and "COMMAND" are 1st and 2nd argument, and shift is equal to 2
    shift = 2;
  } else if (T == use_colvar) {
    // Same as above with additional arguments "colvar" and "NAME"
    shift = 4;
  } else if (T == use_bias) {
    shift = 4;
  }
  return shift;
}


extern "C" {

  /// Generic wrapper for string-based scripting
  int run_colvarscript_command(int objc, unsigned char *const objv[]);

  /// Get the string result of a script call
  const char * get_colvarscript_result();

}


#endif // #ifndef COLVARSCRIPT_H
