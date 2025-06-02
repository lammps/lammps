/***************************************************************************
															commandline.cpp
															W. Michael Brown
														-------------------

	Command line parsing stuff..
	*  Add argument with ' ' for arguments that do not use - (i.e. filenames)
	*  Manditory for the ' ' requires only 1 argument regardless of the number
		- use argsize;

		begin                : Sun Jun 11 2003
		copyright            : (C) 2003 by W. Michael Brown
		email                : mbrown@nirvana.unm.edu
***************************************************************************/

#include "commandline.h"

CommandLine::Parameter::Parameter() {
  num_args=0;
  manditory_args=0;
  manditory=false;
  dash=false;
  set=false;
}

CommandLine::Parameter::~Parameter() {}

CommandLine::CommandLine() {
	help_set=false;
	progname="";

	// Set up the default man chapters
	manchapter_order[0]="NAME";
	man_chapters["NAME"].clear();
	manchapter_order[1]="VERSION";
	man_chapters["VERSION"].clear();
	manchapter_order[2]="SYNOPSIS";
	man_chapters["SYNOPSIS"].clear();
	manchapter_order[3]="DESCRIPTION";
	man_chapters["DESCRIPTION"].clear();
	manchapter_order[4]="PARAMETERS";
	man_chapters["PARAMETERS"].clear();
	man_numchapters=5;
}

CommandLine::~CommandLine() {
}

void CommandLine::addargname(char n, const string &an) {
	parameters[n].argnames.push_back(an);
}

void CommandLine::adddescription(char n, const string &d) {
	parameters[n].description.push_back(d);
}

// Returns the argument size for a parameter
unsigned CommandLine::argsize(char n) {
	return parameters[n].args.size();
}

bool CommandLine::operator [](char n) {
	return set(n);
}

string CommandLine::program_name() {
	return progname;
}

string CommandLine::full_command_line() {
	return full_line;
}

// ----------------- Add allowed arguments

void CommandLine::add(char n, unsigned num) {
  check(n);
  parameters[n].num_args=num;
  parameters[n].manditory_args=num;
}

// Where num represents maximum number of arguments and man_num represents
// manditory number of arguments for a parameter
void CommandLine::add(char n, unsigned num, unsigned man_num) {
  check(n);
  parameters[n].num_args=num;
  parameters[n].manditory_args=man_num;
}

void CommandLine::addmanditory(char n, unsigned num) {
  check(n);
  parameters[n].manditory_args=num;
  parameters[n].num_args=num;
  parameters[n].manditory=true;
}

// Where num represents maximum number of arguments and man_num represents
// manditory number of arguments for a parameter
void CommandLine::addmanditory(char n, unsigned num, unsigned man_num) {
  check(n);
  parameters[n].num_args=num;
  parameters[n].manditory_args=man_num;
  parameters[n].manditory=true;
}

// Allow parameter that takes signed numbers as input
void CommandLine::addsigned(char n, unsigned num) {
  check(n);
  parameters[n].num_args=num;
  parameters[n].manditory_args=num;
  parameters[n].dash=true;  // Allow dashes in arguments
}    

// Allow parameter that takes signed numbers as input
void CommandLine::addsigned(char n, unsigned num, unsigned man_num) {
  check(n);
  parameters[n].num_args=num;
  parameters[n].manditory_args=man_num;
  parameters[n].dash=true;  // Allow dashes in arguments
}

void CommandLine::addhelp(char n, unsigned num) {
  check(n);
  parameters[n].num_args=num;
  parameters[n].manditory_args=num;

  help_set=true;
  help_param=n;
}

// Add descriptions for man pages
void CommandLine::addargnames(char n, unsigned num, const string args[]) {
	for (unsigned i=0; i<num; i++)
		parameters[n].argnames.push_back(args[i]);
}

void CommandLine::adddescription(char n, unsigned num, const string d[]) {
	for (unsigned i=0; i<num; i++)
  	parameters[n].description.push_back(d[i]);
}

void CommandLine::addtoman_chapter(const string &name, const string &body) {
  if (man_chapters.find(name)==man_chapters.end()) {
		manchapter_order[man_numchapters]=name;
		man_numchapters++;
  }
	man_chapters[name].push_back(body);
}

void CommandLine::addtoman_chapter(const string &name, unsigned line_count,
																	 const string body[]) {
  if (man_chapters.find(name)==man_chapters.end()) {
		manchapter_order[man_numchapters]=name;
		man_numchapters++;
  }
	for (unsigned i=0; i<line_count; i++)
		man_chapters[name].push_back(body[i]);
}

void CommandLine::check(char n) {
  map<char,Parameter>::iterator m;
  m=parameters.find(n);
  if (m!=parameters.end()) {
    cerr << "DEVELOPER ERROR: Parameter: " << n << " set twice!\n";
    cerr << "commandline.h: Exiting...\n\n";
    exit(1);
  }
}

bool CommandLine::optargparams() {
  map<char,Parameter>::iterator m;
  for (m=parameters.begin(); m!=parameters.end(); m++)
    if (m->second.manditory_args<m->second.num_args)
      return true;

  return false;
}  

// Returns the number of arguments that are not parameters (' ' name)
unsigned CommandLine::argsize() {
  return parameters[' '].args.size();
}

// Whether or not this parameter was set
bool CommandLine::set(char n) {
  return parameters[n].set;
}

// Force a parameter to be unset
void CommandLine::unset(char n) {
  parameters[n].set=false;
}

char* CommandLine::arg(char n, unsigned num) {
  return parameters[n].args[num];
}

int CommandLine::argint(char n, unsigned num) {
  return atoi(parameters[n].args[num]);
}

double CommandLine::argdouble(char n, unsigned num) {
  return atof(parameters[n].args[num]);
}

string CommandLine::argstring(char n, unsigned num) {
  string s;
  s=parameters[n].args[num];
  return s;
}

bool CommandLine::parse(int argc, char* argv[], Error *error) {
  unsigned i;
  map<char,Parameter>::iterator m;
  char flag;
  int parameter; // Set to the parameter that arguments are being parsed

  progname=a::filenameonly(argv[0]);
  full_line=string(argv[0]);
  for (unsigned i=1; i<unsigned(argc); i++)
  	full_line+=' '+string(argv[i]);

  int argnum=1;
  while (argnum<argc) {
    // Set an argument
    if (argv[argnum][0]!='-') {
      if (parameters[' '].args.size()==parameters[' '].num_args) {
				error->addwarning(3,9,"CommandLine","Invalid Command Line Argument: "+
															              string(argv[argnum]));
				return false;
      }
      parameters[' '].set=true;
      parameters[' '].args.push_back(argv[argnum]);
      argnum++;
      continue;
    }

    // Set a parameter
    flag=argv[argnum][1];
    m=parameters.find(flag);
    if (m==parameters.end()) {
			error->addwarning(4,9,"CommandLine","Invalid Command Line Parameter: "+
			                                    string(argv[argnum]));
      return false;
    }
    // Make sure all required arguments are set before parameters with
    // optional arguments
    if (m->second.manditory_args<m->second.num_args)
      if (parameters[' '].args.size()<parameters[' '].num_args) {
        error->buffer() << "Parameters with optional arguments must be "
                        << "placed after manditory commandline arguments";
        error->addbuf(5,9,"CommandLine");
        return false;
      }
    
    parameter=argnum;
    argnum++;
    m->second.set=true;

    // Get the parameter arguments
    for (i=0; i<m->second.num_args; i++) {
      // Make sure we are not at the end of the parameter list
      if (argnum>=argc)
				if (m->second.args.size()<m->second.manditory_args) {
					error->buffer() << "Invalid Number of Arguments: -" << m->first
													<< " requires " << m->second.num_args
             							<< " arguments!";
					error->addbuf(5,9,"CommandLine");
					return false;
        } else
          break;
      
      // Assure the right number of arguments for signed numbers
			if (argv[argnum][0]=='-' && m->second.dash==true)
				if (!isdigit(argv[argnum][1]))
					if (m->second.args.size()<m->second.manditory_args) {
					error->buffer() << "Invalid Number of Arguments: -" << m->first
													<< " requires " << m->second.num_args
             							<< " arguments!";
					error->addbuf(5,9,"CommandLine");
					return false;
				} else
				  break;
							
      // Assure the right number of arguments for other numbers
      if (argv[argnum][0]=='-' && m->second.dash==false)
        if (m->second.args.size()<m->second.manditory_args) {
					error->buffer() << "Invalid Number of Arguments: -" << m->first
													<< " requires " << m->second.num_args
             							<< " arguments!";
					error->addbuf(5,9,"CommandLine");
					return false;
			  } else
          break;

      m->second.args.push_back(argv[argnum]);
      argnum++;
    }
  }


  // If help was set, we do not need to check for manditory args
  if (help_set)
		if (parameters[help_param].set)
			return true;
  
  // Assure that the manditory arguments were set for commandline
  if (parameters[' '].manditory) 
		if (parameters[' '].args.size()<parameters[' '].manditory_args) {
			error->addwarning(6,9,"CommandLine","Missing Required Argument!");
			return false;
		}

  // Assure that manditory arguments were set for parameters!
  for (m=parameters.begin(); m!=parameters.end(); m++) {
    if (m->second.manditory)
      if (m->second.set==false) {
      	error->buffer() << "Missing Required Argument: \n\n";
				if (m->first!=' ')
					error->buffer() << "-" << m->first << " must be set!";
				error->addbuf(7,9,"CommandLine");
				return false;
      }
  }
  return true;
}

void CommandLine::write_man_page(ostream & out, const string &version,
																 const string &header) {
	unsigned i;
	map<char,Parameter>::iterator m;
	string bold="\\fB";
	string italic="\\fI";
	string regular="\\fR";
																			
	out << ".TH " << program_name() << " 1 \"" << a::date()
			<< "\"" << " \"" << program_name() << " (" << header << ") "
			<< version << "\" \"" << header << '\"' << endl;

  // Go through the chapters
  map<unsigned,string>::iterator mc;
  for (mc=manchapter_order.begin(); mc!=manchapter_order.end(); mc++) {
    // NAME Section
		if (mc->second=="NAME") {
			out << ".SH NAME\n" << bold << program_name() << regular;
			if (man_chapters["NAME"].size()!=0) {
				out << " - ";
				for (i=0; i<man_chapters["NAME"].size(); i++)
					out << man_chapters["NAME"][i];
			}
			out << endl;
			continue;
		} // end NAME

		// SYNOPSIS section
		if (mc->second=="SYNOPSIS" && man_chapters["SYNOPSIS"].empty()) {
			out << ".PD 2\n.SH SYNOPSIS\n.PD 1\n.TP\n";
      if (program_name().length()>6)
				out << bold << program_name() << regular << " ";
			else
				out << ".B " << program_name() << endl;
			out << format_synopsis(bold,italic,regular) << endl;
      if (program_name().length()>6)
				out << ".br\n";
      if (parameters[' '].num_args!=0 && optargparams())
        out << ".PD 2\n.PP\nParameters with optional arguments should be placed "
            << "after manditory commandline arguments.\n";
			continue;
		} // end SYNOPSIS
	
		// PARAMETERS
		if (mc->second=="PARAMETERS") {
			if (!(parameters.size()==0 ||
				   (parameters.size()==1 && (parameters.begin())->first==' '))) {
			  out << ".PD 2\n.SH PARAMETERS\n.PD 1\n";
			  for (m=parameters.begin(); m!=parameters.end(); m++) {
					if (m->first==' ')
						continue;
					out << ".TP\n";
					out << format_parameter(m->first,bold,italic,regular) << endl;
					if (format_parameter(m->first,"","","").length()>6)
						out << ".PD 0\n.TP\n.PP\n.PD 1\n";
					// Output the description
					vector<string> fdesc=man_format(m->second.description,bold,
																					italic,regular);
					for (i=0; i<fdesc.size(); i++)
						out << fdesc[i];
					out << endl;
			  } // end for m=
			}
			continue;
		} // end PARAMETERS
		writeman_chapter(out,mc->second,bold,italic,regular);
	}
}

void CommandLine::writeman_chapter(ostream &out,const string &name,
																	 const string &bold, const string &italic,
                  								 const string &regular) {
	if (man_chapters[name].empty())
		return;

	out << ".PD 2\n.SH " << name << "\n.PD 1\n";

	vector<string> formatted=man_format(man_chapters[name],bold,italic,regular);
	for (unsigned i=0; i<formatted.size(); i++)
		out << formatted[i];
	out << endl;
}
	
string CommandLine::format_synopsis(const string &bold, const string &italic,
																		const string &regular) {
  string s;
	map<char,Parameter>::iterator m;
  bool space=false;
	
	s=format_parameter(' ',bold,italic,regular)+" ";
	for (m=parameters.begin(); m!=parameters.end(); m++) {
		if (m->first==' ')
			continue;
		if (space)
			s+=' ';
		space=true;
		if (m->second.manditory==false)
			s+='[';
		s+=format_parameter(m->first,bold,italic,regular);
		if (m->second.manditory==false)
			s+="]";
	} // end for
	return s;
}
	
// Write a synopsis in plain text format fitted to a given column width
void CommandLine::write_text_synopsis(ostream &out, unsigned column_width) {
  string s=format_synopsis("","","");
	vector<string> lines;
	a::format_fit(column_width-program_name().length()-1,s,lines);
	out << program_name() << " ";
	for (unsigned i=0; i<lines.size(); i++) {
		out << lines[i] << endl;
		for (unsigned j=0; j<(program_name().length()+1); j++)
			out << " ";
  }
  out << endl;
}

string CommandLine::format_parameter(char param, const string &bold,
																		 const string &italic,
                   									 const string &regular) {
  Parameter &p=parameters[param];
												
  string strp="";
	if (param!=' ')
		strp=bold+"-"+param+regular;
		
	if (p.num_args!=0)
 		for (unsigned i=0; i<p.argnames.size(); i++) {
      if (!(param==' ' && i==0)) // Ignore first space on command args
				strp+=' ';
		  if (i>=p.manditory_args)
 				strp+='[';
      if (param!=' ')
				strp+=italic+p.argnames[i]+regular;
	    else
				strp+=p.argnames[i];
 			if (i>=p.manditory_args)
 				strp+=']';
    }

  if (p.argnames.size()<p.num_args)
   	strp+=" [...]";

	return strp;
}

vector<string> CommandLine::man_format(const vector<string> &input,
																			 const string &bold, const string &italic,
																			 const string &regular) {
	vector<string> output;
	map<char,Parameter>::iterator m;

	for (unsigned i=0; i<input.size(); i++) {
		string temp=input[i];
		a::str_replace("\n","\n.PD 0\n.PP\n.PD 1\n",temp);
		a::str_replace(".EN","\n.EN\n",temp);
		a::str_replace(".EQ","\n.EQ\n",temp);
		a::str_replace(".TP","\n.PD 0\n.TP\n.PP\n.PD 1\n",temp);
		a::str_replace(program_name(),bold+program_name()+regular,temp);
		for (m=parameters.begin(); m!=parameters.end(); m++) {
      if (m->first!=' ') {
				string source="-";
				source+=m->first;
				a::str_replace(source,bold+source+regular,temp);
		  }
			for (unsigned j=0; j<m->second.argnames.size(); j++)
				a::str_replace(m->second.argnames[j],
											 italic+m->second.argnames[j]+regular,temp);
		}
		output.push_back(temp);
	}

	return output;
}

