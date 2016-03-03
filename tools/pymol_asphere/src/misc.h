/***************************************************************************
                                    misc.h
                             -------------------
                               W. Michael Brown

  Miscellaneous functions that do not deserve their own class

 __________________________________________________________________________
    This file is part of the "All" Library
 __________________________________________________________________________

    begin                : May 30 2003
    copyright            : (C) 2003 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

#ifndef MISC_H
#define MISC_H

#include "error.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>

using namespace std;

/// Miscellaneous functions that do not deserve their own class
/** \e a contains functions for \n
  * - fileio
  * - string manipulation
  * - conversion from numbers to strings */
namespace a {
  /// Output the SNL Copyright info for publication \e year.
  void copyright(ostream &out, unsigned year);
  /// Get the SNL Copyright info as a string
  string copyrightstring(unsigned year);

	/// Open a file for input. Generates error ID \b 1-15 on fail.
	void fileopen(ifstream &in, const char *filename, Error &error);
	/// Open a file for input. Generates error ID \b 1-15 on fail.
	void fileopen(ifstream &in, const string &filename, Error &error);
	/// Open a binary file for input. Generates error ID \b 1-15 on fail.
	void fileopenbinary(ifstream &in, const string &filename, Error &error);

	/// Open a file for output. Generates error ID \b 2-15 on fail.
	void fileopen(ofstream &out, const string &filename, Error &error);
	/// Open a file for output. Generates error ID \b 2-15 on fail.
	void fileopenbinary(ofstream &out, const string &filename, Error &error);
	/// Open a file for append. Generates error ID \b 2-15 on fail.
	void fileopenapp(ofstream &out, const string &filename, Error &error);
	/// Open a file for output. Generates error ID \b 2-15 on fail.
	void fileopen(ofstream &out, const char *filename, Error &error);

	/// Close an input file. Generates error ID \b 10-15 on fail.
	void fileclose(ifstream &in, Error &error);
	/// Close an output file. Generates error ID \b 11-15 on fail.
	void fileclose(ofstream &out, Error &error);

  /// Put a string back into an istream
  void putback(istream &in, const string &s);
  
  /// Get the current date in the format "January 1, 2003"
  string date();

	/// Returns the filename without the extension
	string namewoext(const string &filename);
	/// Returns the filename without extension or directory
	string filenameonly(const string &filename);
	/// Return the extension of a filename
	string extension(const string &filename);
	/// Centers a string over the specified length
	string strcenter(const string &s, unsigned length);
  /// True if a character is whitespace
  bool whitespace(char c);
	/// True if a string is only whitespace
	bool whitespace(const string &s);
	/// Remove all whitespace from a string
	string remove_whitespace(const string &s);
	/// Replace any instance of \e source with \e target within the string \e s
	void str_replace(const string &source, const string &target, string &s);
  /// Convert all alpha characters to lower case
  string tolower(const string &s);

  /// The tokens parsed from cstring are \e added to the input vector
  /** \param line The line with 'white space' delimeted tokens
    * \param tokens Each token parsed is added to the vector */
  void get_tokens(const char *line, vector<string> &tokens);
  /// The tokens parsed from string are \e added to the input vector
  /** \param line The line with 'white space' delimeted tokens
    * \param tokens Each token parsed is added to the vector */
  void get_tokens(const string &line, vector<string> &tokens);
  /// Parse a string into tokens based on delimiter
  /** \param line The line with 'white space' delimeted tokens
    * \param tokens Each token parsed is added to the vector */
  void get_tokens(char delimiter, const string &line, vector<string> &tokens);
  /// Return the first token in a string
  string get_first_token(const char *line);
  /// Format a string to fit within a specified column width
  /** Newlines are inserted between whitespace if possible, otherwise line is
    * wrapped. Each string in the vector represents one line of text
    * \note The output vector is not emptied! **/
  void format_fit(unsigned column_width, const string &input,
  								vector<string> &output);
  
  /// Return a string of num underscores
  string underline(unsigned num);
	
	/// Returns string representation of unsigned number
	string itoa(unsigned i);
	/// Returns string representation of int number
	string itoa(int i);
	/// Returns string representation of double number
	string ftoa(double i);

	/// Seed the random number generator
	void seedrandom(unsigned seed);
	/// Seed the random number generator with the current time
	void seedrandom_time();
	/// Returns a random integer between 0 and max
	unsigned irandom(unsigned max);
	/// Returns a random integer between 0 and max
	long lrandom(long max); 
	/// Returns a random double between 0 and max
 	double frandom(double max);
}

/// Iterate through file names to give each unique numbers
class FileIterator {
 public:
  /// Empty constructer with no header, no extension, no lead zeros
  FileIterator();
  /// Specify the filename format using leading zeros
  /** Files are generated according to the following format:
    * \verbatim header+%0'digits'file_number+extension \endverbatim */
	FileIterator(const string &header, const string &extension, unsigned digits);
  /// Specify the filename format without leading zeros
  /** Files are generated according to the following format:
    * \verbatim header+file_number+extension \endverbatim */
	FileIterator(const string &header, const string &extension);

  /// Set the current file number
  void set_file_num(unsigned fnum);
  /// Set the file header
  void set_file_header(const string &head);
  /// Set the file extension
  void set_file_extensions(const string &ext);
  /// Set the number of leading zeros
  void set_lead_zeros(unsigned digs);
  /// Set the current file number, header, extension, leading zeros
  void set(unsigned fnum, const string &head, const string &ext, unsigned digs);
  
  /// Returns the next filename.
  string nextfilename();
 private:
  string header,extension;
  unsigned digits;
  unsigned file_num;	
};

#endif
