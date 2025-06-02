/***************************************************************************
                                  error.h 
                             -------------------

  Class for error handling

  __________________________________________________________________________

    begin                : Thu Oct 9 2003
    copyright            : (C) 2003 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

#ifndef ERRORCLASS
#define ERRORCLASS

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>

#ifdef MUSE_MPI
#include <mpi.h>
#endif

using namespace std;
// forward declarations
namespace a {
	string itoa(unsigned int);
	void format_fit(unsigned, const string &, vector<string> &);
}

/// Notice Class for Handling Object Output
/** A notice object stores an ostream for output and a notice_level.
  * All output is sent along with a level.  Any output whose level is
  * greater than notice_level is sent to a null stream. The C++ output
  * operator '<<' can be used with the Notice operator '[]' which passes
  * the level:
  * \verbatim notice_object[29] << "This notice has level 29" << endl;
  * \endverbatim
  *
  * The guidelines for output notice levels are:
  * - \e     0: Information that must be output
  * - \e 1 - 9: Normal program output
  * - \e 10-19: Parameters useful for storing how the program was run
  * - \e 20-29: Extraneous information useful that may be useful to user
  * - \e 30-  : Debugging information */
class Notice {
 public:
  /// Standard output (cout) is the default notice output ostream
  /** The default maximum notice level is 9 \e (notice_level=10) */
  Notice();
  ~Notice();
  
  /// Set the output stream for notice output
  void setostream(ostream &out);
  
  /// Returns a null stream if level is two high, else returns notice stream
  ostream & operator[] (const unsigned level);
  
	/// Set the degree of program output
  void set_notice_level(unsigned l);
  /// Get the degree of program output
  unsigned get_notice_level();
  /// Generate a notice with a specified calling class
  void notice(unsigned level, const string calling_class, const string note);
  /// Generate a notice with a specified calling class
	void notice(unsigned level, const string calling_class,
 							vector<string> &notes);
	/// Generate a notice
	void notice(unsigned level, const string note);
	/// Generate a notice
	void notice(unsigned level, vector<string> &note);

 private:
  unsigned notice_level;

	ostream *nullout; // Null stream for redirecting output to nowhere
	ostream *noteout; // Output for notices
};

/// Error and Notice Handling
/** This class is intended to handle all output to the user. Output is
  * divided into notices and warnings/errors. Output of any message is
  * associated with a level. For notices, if level is greater than or equal to
  * max_notice_level, no output is generated. For warnings, if level is less
  * than min_warning_level, it is dismissed with no output. If the level is
  * greater than max_warning_level, the program is terminated and all warnings
  * and errors are output.
  *
  * \note By default, on destruction of an Error object, all unhandled
  * warnings and errors are output
  *
  * A log file can be specified for each object. In this case, all notices
  * are output to the log file only and errors are output to both stderr
  * and the log file
  *
  * Errors can be generated with a string or using the internal message buffer:
  \verbatim
  Error error;
  error.buffer() << "Incorrect file format for file: " << filename;
  error.addbuf(512,19,"FooClass";
  // --- OR
  string message = "Incorrect file format for file: "+filename;
  error.addwarning(512,19,"FooClass",message);
  \endverbatim
  *
  * Newlines will be inserted into the error message automatically in order
  * to format it for the string. Forced newlines can be specified with \n
  *
  * Programs can check whether or not errors have been generated using the []
  * operator and can 'handle' them by outputting the message or dismissing
  * them without any output
  *
  * Notices are generated using the public Notice class (see Notice())
  *
  * \b Error \b  Levels:
  *  - \e 0 - 1:  Errors expected to happen during normal execution
  *  - \e 2 - 9:	Errors that a non-interactive program can handle and continue
  *  - \e 10-19:  Errors that interactive program can handle (file not found,etc.)
  *  - \e 20-  :  Serious errors that should terminate execution 
  **/
class Error {
 	public:
 		/// Default constructor (use cerr for output and no log file)
 		/** Default max notice level is 9, min warning level is 2, and max warning
      * level is 9 */
    Error();
		~Error();	

    /// Set a log file for error AND notice output
    void set_logfile(ostream &out);
			
		/// Returns the number of errors (if any) generated with id
	  unsigned operator[](unsigned id);
	  
	  /// Add warning, terminate if level is greater than max level
	  /** Newlines will be inserted into the message automatically when the
      * message is formatted for output. However, forced newlines can also
      * be inserted. **/
		void addwarning(unsigned ID, unsigned level, const string calling_class,
										const string warning);
		/// Add serious error (terminates execution)
		void generate_error(unsigned ID, const string calling_class,
												const string error);

	  /// Add an message to the error buffer. Warning generated with addbuf()
	  /** Newlines will be inserted into the message automatically when the
      * message is formatted for output. However, forced newlines can also
      * be inserted.
      *
      \verbatim
      Error error;
      error.buffer() << "Choice not supported: " << choice;
      error.addbuf(512,9,"FooClass");
      \endverbatim **/
    ostringstream & buffer();
    /// Generate warning with message in buffer
    /** \sa buffer() **/
    void addbuf(unsigned ID, unsigned level, const string calling_class);
    /// Generate serious error with message in buffer
    /** \sa buffer() **/
    void addbuf(unsigned ID, const string calling_class);

		/// Number of Unhandled Warnings
		unsigned warnings();
		/// Total number of warnings
		unsigned total_warnings();
		
		/// Handle all warnings with this ID by writing them to out
		void writewarning(unsigned ID);
		/// Handle LAST warning with this ID WITHOUT writing it
		void dismiss_warning(unsigned ID);
		/// Handle ALL warnings with this ID WITHOUT writing it
		void dismiss_all_warnings(unsigned ID);
		/// Handle all warnings by writing them out
		void writewarnings();
		/// Handle all warnings without writing
		void dismiss_warnings();

		/// Write out the total warnings (write errorcount errors)
		void writetotals(unsigned errorcount);

    /// For generating notices
    /** \sa Notice **/
    Notice note;
	private:
		struct ErrCom;
		map<unsigned,vector<ErrCom> > warning_list;
		typedef multimap<unsigned, vector<ErrCom> >::iterator warning_iter;
		unsigned handled_warnings;
		unsigned unhandled_warnings;
    bool handleatend;             // Write any unhandled errors on destruct
    bool writetotalatend;         // Write totals on destruct if not 0
		
    unsigned min_level;           // Don't output warnings less than min_level
		unsigned max_level; 					// if a warning has a level>max_level error!
		ostream *errout, *logout;			// Output for errors and warnings!
		ostream *nullout;             // No output

		ostringstream buffer_stream;  // For creating messages for warnings

		unsigned column_width;
		void write_err(unsigned ID, ErrCom &err);
		void writeline();
};

struct Error::ErrCom {
	unsigned level;
	string calling_class;
	string message;
};

#endif

