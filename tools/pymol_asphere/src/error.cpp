/***************************************************************************
                                  error.cpp 
                             -------------------

  Class for error handling

  __________________________________________________________________________

    begin                : Thu Oct 9 2003
    copyright            : (C) 2003 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

#include "error.h"
#include <cstring>
 
Notice::Notice() {
	nullout=new ostream(NULL);
	noteout=&cout;
	notice_level=9;
}

Notice::~Notice() {
	if (nullout!=NULL)
		delete nullout;
}

// Returns a null stream if level is two high, else returns notice stream
ostream & Notice::operator[] (const unsigned level) {
	if (level<=notice_level)
		return *noteout;
	else
		return *nullout;
}

void Notice::setostream(ostream &out) {
	noteout=&out;
}

void Notice::set_notice_level(unsigned l) {
	notice_level=l;
}

unsigned Notice::get_notice_level() {
  return notice_level;
}

void Notice::notice(unsigned level, const string calling_class,
									  const string note) {
  if (level<notice_level)
    *noteout << calling_class << "::" << note;
}

void Notice::notice(unsigned level,const string calling_class,
										vector<string> &notes) {
  if (level<=notice_level) {
    *noteout << calling_class << "::";
      for (unsigned i=0; i<notes.size(); i++)
        *noteout << notes[i];
  }
}

void Notice::notice(unsigned level, const string note) {
  if (level<=notice_level)
    *noteout << note;
}

void Notice::notice(unsigned level, vector<string> &notes) {
  if (level<=notice_level)
    for (unsigned i=0; i<notes.size(); i++)
      *noteout << notes[i];
}

Error::Error() {
	nullout=new ostream(NULL);
	unhandled_warnings=0;
	handled_warnings=0;
	max_level=9;
	min_level=2;

	handleatend=true;
  writetotalatend=true;

	errout=&cerr;
	logout=nullout;
	column_width=70;
}
 
Error::~Error() {
	if (handleatend && unhandled_warnings!=0)
		writewarnings();
	if (writetotalatend && total_warnings()!=0)
		writetotals(0);
	if (nullout!=NULL)
		delete nullout;
}
	
// Set a log file for error AND notice output
void Error::set_logfile(ostream &out) {
	logout=&out;
	note.setostream(out);
}

// Total number of warnings
unsigned Error::total_warnings() {
	return unhandled_warnings+handled_warnings;
}

// Returns the number of errors generated with ID
unsigned Error::operator[](unsigned id) {
	warning_iter m;
	m=warning_list.find(id);
	if (m==warning_list.end())
		return 0;
  return m->second.size();
}

void Error::addwarning(unsigned ID, unsigned level, const string calling_class,
											 const string warning) {
  if (level<min_level)
  	return;
  if (level>max_level)
  	generate_error(ID,calling_class,warning);
  vector<ErrCom> *e=&(warning_list[ID]);
  e->push_back(ErrCom());
  e->back().level=level;
  e->back().calling_class=calling_class;
  e->back().message=warning;
  unhandled_warnings++;
}

void Error::generate_error(unsigned ID, const string calling_class,
													 const string error) {
	ErrCom err;
	err.level=max_level+1;
	err.calling_class=calling_class;
	err.message=error;

	if (warnings()!=0)
		writewarnings();
	write_err(ID,err);
  writetotals(1);
  #ifdef MUSE_MPI
  MPI_Abort ( MPI_COMM_WORLD, 1 );
  #else
  exit(1);
  #endif
}

// Add an error/warning (Program termination if level >= max level
ostringstream & Error::buffer() {
  return buffer_stream;
}

// Generate warning with message in buffer
void Error::addbuf(unsigned ID, unsigned level, const string calling_class) {
	addwarning(ID,level,calling_class,buffer_stream.str());
	buffer_stream.str("");
}

// Generate serious error with message in buffer
void Error::addbuf(unsigned ID, const string calling_class) {
	generate_error(ID,calling_class,buffer_stream.str());
}

unsigned Error::warnings() {
	return unhandled_warnings;
}

void Error::writeline() {
  *errout << "+";
  *logout << "+";
  for (unsigned i=0; i<column_width-2; i++) {
  	*errout << "-";
  	*logout << "-";
  }
  *errout << "+\n";
  *logout << "+\n";
}

void Error::write_err(unsigned ID, ErrCom &err) {
  if (err.level<min_level)
		return;
		
  *errout << "\n";
  *logout << "\n";
  writeline();
  (*errout).setf(ios::left);
  (*errout).unsetf(ios::right);

  // Output the IDs
  unsigned width12=(unsigned)floor((double)(column_width-10)/3.0);
  unsigned width3=column_width-10-width12-width12;
  string et;
  unsigned width1;
  if (err.level>max_level) {
  	et="| Error: "; width1=width12-7;
  } else {
  	et="| Warning: "; width1=width12-9;
  }
  *errout << et << setw(width1) << ID << " | Level: "
  				<< setw(width12-7) << err.level << " | " << setw(width3)
      		<< err.calling_class << " |\n";
  *logout << et << setw(width1) << ID << " | Level: "
  				<< setw(width12-7) << err.level << " | " << setw(width3)
      		<< err.calling_class << " |\n";
  writeline();

  // Output the message
  vector <string> messages;
  a::format_fit(column_width-3,err.message,messages);
  for (unsigned i=0; i<messages.size(); i++) {
  	*errout << "| " << setw(column_width-3) << messages[i] << "|\n";
  	*logout << "| " << setw(column_width-3) << messages[i] << "|\n";
  }
  writeline();

	return;
}

void Error::writewarning(unsigned ID) {
 	for (unsigned i=0; i<warning_list[ID].size(); i++)
 		write_err(ID,warning_list[ID][i]);
  handled_warnings+=warning_list[ID].size();
	dismiss_all_warnings(ID);
	return;
}

void Error::dismiss_warning(unsigned ID) {
	warning_iter m;
	m=warning_list.find(ID);
	if (m==warning_list.end())
		return;
	unhandled_warnings--;
 	if (m->second.size()==1)
  	warning_list.erase(m);
  else
    m->second.erase(m->second.end()--);
}

void Error::dismiss_all_warnings(unsigned ID) {
	warning_iter m;
	m=warning_list.find(ID);
	if (m==warning_list.end())
		return;
	unhandled_warnings-=m->second.size();
  warning_list.erase(m);
}

void Error::writewarnings() {
	while (warning_list.size()>0)
		writewarning(warning_list.begin()->first);
	return;
}

void Error::dismiss_warnings() {
	while (warning_list.size()>0)
		dismiss_warning(warning_list.begin()->first);
	return;
}

// Write out the total warnings and errors
void Error::writetotals(unsigned errorcount) {
  *errout << "\n";
  *logout << "\n";
  writeline();
  (*errout).setf(ios::left);
  (*errout).unsetf(ios::right);
  unsigned width1=(unsigned)floor((double)(column_width-7)/2.0);
  unsigned width2=column_width-7-width1;
  string swarnings="Total Warnings: "+a::itoa(handled_warnings+
  																					 unhandled_warnings);
  string serrors="Total Errors: "+a::itoa(errorcount);
  *errout << "| " << setw(width1) << swarnings << " | " << setw(width2)
      		<< serrors << " |\n";
  *logout << "| " << setw(width1) << swarnings << " | " << setw(width2)
      		<< serrors << " |\n";
  writeline();
}

