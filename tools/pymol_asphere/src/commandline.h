/***************************************************************************
                                commandline.h
                               W. Michael Brown
                             -------------------

  Command line parsing stuff..
                             
    begin                : Sun Jun 11 2003
    copyright            : (C) 2003 by W. Michael Brown
    email                : mbrown@nirvana.unm.edu
 ***************************************************************************/

#ifndef COMMANDLINE_H
#define COMMANDLINE_H

#include "error.h"
#include "misc.h"
#include <map>
#include <vector>
#include <string>
#include <stdlib.h>
#include <iostream>
using namespace std;
   
/// Parsing of command-line parameters and automatic man page generation
/** Allows manditory and optional command line arguments to be specified along
  * with parameters set using flags. For arguments that do not require flags,
  * a space char is used to set the flag. To specify the command line:
  *
  *   foo input_file output_file -t signed_number [-o]
  *
  * \verbatim CommandLine cl;
  cl.addmanditory(' ',2);
  cl.addargname(' ',"input_file"); cl.addargname(' ',"output_file");
  cl.addsigned('t',1,1);
  cl.addargname('t',"arg_name");
  cl.adddescription('t',"Manditory parameter allowing signed_numbers");
  cl.add('o',0);
  \endverbatim
  * To instead specify that the outputfile is an optional argument:
  * \verbatim cl.addmanditory(' ',2,1) \endverbatim
  *
  * A help flag can also be set which does not require other command line
  * arguments to be set
  *
  * When the commandline is parsed, it is asserted that all manditory arguments
  * and flags are set, that each flag is set with the correct number of
  * parameters, and that no unknown flags have been set.
  *
  * One can check whether or not optional flags have been set using set() or
  * the [] operator:
  * \verbatim cl.set('o')  <-->  cl['o']  \endverbatim
  *
  * Man pages can be generated automatically by adding chapters. The SYNOPSIS
  * and description of commandline arguments and flags is generated
  * automatically. Examples of chapters that are consistently used are
  * NAME, VERSION, DESCRIPTION, PARAMETERS, USAGE, EXAMPLES, AUTHORS, BUGS,
  * SEE ALSO.
  *
  * The following characters can be used for formatting in parameter
  * descriptions and chapters:
  * \verbatim
  \t		tab
  \n		forced newline
  .TP		new paragraph (margin indented)
  \\fB  all following text is bold
  \\fI	all following text is italic
  \\fR	all following text is regular
  \endverbatim
  *
  * The arguments, flags and parameter_names are automatically formatted
  * throughout the man page with bold and italic typeface.
  **/
class CommandLine {
 public:
  CommandLine();
  ~CommandLine();

  /// Add an optional argument which takes num parameters
  /** \note Use ' ' for n to specify an argument with no flag
    * \param n flag used to pass arguments (i.e. -f)
    * \param num the number of arguments that must be specified after flag */
  void add(char n, unsigned num);
  /// Add a manditory argument which takes num parameters
  /** \sa add() **/
  void addmanditory(char n, unsigned num);
  /// Add a flag which can take signed numbers as parameters
  /** \sa add() **/
  void addsigned(char n, unsigned num);

 	/// Add an optional flag with minimum and maximum number of args
  /** \sa add() **/
	void add(char n, unsigned num, unsigned man_num);
 	/// Add a manditory flag with minimum and maximum number of args
  /** \sa add() **/
  void addmanditory(char n, unsigned num, unsigned man_num);
 	/// Add a flag with minimum and maximum number of args that takes signed nums
  /** \sa add() **/
  void addsigned(char n, unsigned num, unsigned man_num);

  /// Add a help argument (does not require other manditory arguments be set)
  /** \sa add() **/
  void addhelp(char n, unsigned num);
  
  /// Specify the names for arguments in the synopsis (for man page)
  /** The names can be added one by one in the order of the flag parameters **/
  void addargname(char n, const string &an);
  /// Specify the names for arguments in the synopsis (for man page)
  void addargnames(char n, unsigned num, const string args[]);
  /// Specify the description for an argument (for the man page SYNOPSIS)
  /** See the class description for formating characters **/
  void adddescription(char n, const string &d);
  /// Specify the description for an argument (for the man page SYNOPSIS)
  /** See the class description for formating characters **/
  void adddescription(char n, unsigned num, const string d[]);
  /// Add a man page chapter with title 'name'
  void addtoman_chapter(const string &name,const string &body);
  /// Add a man page chapter with title 'name'
  void addtoman_chapter(const string &name,unsigned line_count,
  											const string body[]);
  
  /// Returns the number of arguments that are not flags (char n=' ')
  unsigned argsize();
  /// Returns the number of parameters passed for a given flag
  unsigned argsize(char n);
  /// Returns true if the optional flag was set on the commandline
  bool set(char n);
  /// Returns true if the optional flag was set on the commandline
  bool operator [](char n);

  /// Force a parameter to be unset
  void unset(char n);
  
  /// Return flag parameter or argument as a cstring (0-based index)
  char *arg(char n, unsigned num);         
  /// Return flag parameter or argument as a integer (0-based index)
  int argint(char n, unsigned num);      
  /// Return flag parameter or argument as a double (0-based index)
  double argdouble(char n, unsigned num);
  /// Return flag parameter or argument as a string (0-based index)
  string argstring(char n, unsigned num);

  /// Parse the command line arguments
  bool parse(int argc, char* argv[], Error *error);
  
  /// Return the program name
  string program_name();
  /// Return a string with the entire commandline as entered
  string full_command_line();
  /// Write out a man_page
	void write_man_page(ostream & out, const string &version,
 											const string &header);
        
  /// Advanced writing of man page
  void writeman_chapter(ostream &out, const string &name, const string &bold,
												const string &italic, const string &regular);

  /// Return a string with all of the options in a man page synopsis format
  string format_synopsis(const string &bold, const string &italic,
  											 const string &regular);
  /// Write a synopsis in plain text format fitted to a given column width
  void write_text_synopsis(ostream &out, unsigned column_width);
  
 private:
  string full_line;
  void check(char n);
  // Returns true if parameters have been set with optional arguments
  bool optargparams();

  class Parameter;
  map<char,Parameter> parameters;

  bool help_set;		// True if there is a parameter for getting help
  char help_param;  // Parameter for help
  
  // Stuff for man page
  string progname; // The name for the program (argv[0])
  unsigned man_numchapters; // Number of chapters in man page
  map<string,vector<string> > man_chapters;
  map<unsigned, string> manchapter_order;

  // Format a parameter as a string with argument names
  string format_parameter(char param, const string &bold, const string &italic,
													const string &regular);
  // Formats the strings in input for a man page by replacing newlines and
  // adding bold and italic fonts to parameters and arguments respectively
	vector<string> man_format(const vector<string> &input, const string &bold,
														const string &italic, const string &regular);

};

// This class is for internal use within CommandLine
class CommandLine::Parameter {
 public:
  friend class CommandLine;
  Parameter();
  ~Parameter();

 private:
  unsigned num_args;       // Maximum number of arguments for the parameter
  unsigned manditory_args; // Minimum number of arguments for the parameter
  bool manditory;			     // This parameter MUST be set

  bool dash;  				     // Allow for signed numbers (dashes in args)

  bool set;                // This parameter has been set
  vector<char *> args;

  // Strings for man page description
  vector<string> argnames; // Names for each argument to the parameter
  vector<string> description; // Description for how to set this parameter
};

#endif
