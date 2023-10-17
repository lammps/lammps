// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "input.h"

#include "accelerator_kokkos.h"
#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "comm.h"
#include "comm_brick.h"
#include "comm_tiled.h"
#include "command.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "integrate.h"
#include "kspace.h"
#include "label_map.h"
#include "memory.h"
#include "min.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "pair.h"
#include "special.h"
#include "style_command.h"      // IWYU pragma: keep
#include "thermo.h"
#include "timer.h"
#include "universe.h"
#include "update.h"
#include "variable.h"

#include <cstring>
#include <cerrno>
#include <cctype>

using namespace LAMMPS_NS;

#define DELTALINE 256
#define DELTA 4

// maximum nesting level of input files
static constexpr int LMP_MAXFILE = 16;

/* ----------------------------------------------------------------------
   one instance per command in style_command.h
------------------------------------------------------------------------- */

template <typename T> static Command *command_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

/* ---------------------------------------------------------------------- */

/** \class LAMMPS_NS::Input
 *  \brief Class for processing commands and input files
 *
\verbatim embed:rst

The Input class contains methods for reading, pre-processing and
parsing LAMMPS commands and input files and will dispatch commands
to the respective class instances or contains the code to execute
the commands directly.  It also contains the instance of the
Variable class which performs computations and text substitutions.

\endverbatim */

/** Input class constructor
 *
\verbatim embed:rst

This sets up the input processing, processes the *-var* and *-echo*
command line flags, holds the factory of commands and creates and
initializes an instance of the Variable class.

To execute a command, a specific class instance, derived from
:cpp:class:`Command`, is created, then its ``command()`` member
function executed, and finally the class instance is deleted.

\endverbatim
 *
 * \param  lmp   pointer to the base LAMMPS class
 * \param  argc  number of entries in *argv*
 * \param  argv  argument vector  */

Input::Input(LAMMPS *lmp, int argc, char **argv) :
    Pointers(lmp), variable(nullptr), labelstr(nullptr), infiles(nullptr), inlines(nullptr),
    command_map(nullptr)
{
  MPI_Comm_rank(world, &me);

  maxline = maxcopy = maxwork = 0;
  line = copy = work = nullptr;
  narg = maxarg = 0;
  arg = nullptr;

  echo_screen = 0;
  echo_log = 1;

  label_active = 0;
  jump_skip = 0;
  utf8_warn = true;

  if (me == 0) {
    nfile = 1;
    infiles = new FILE *[LMP_MAXFILE];
    infiles[0] = infile;
    inlines = new int[LMP_MAXFILE];
  }

  variable = new Variable(lmp);

  // fill map with commands listed in style_command.h

  command_map = new CommandCreatorMap();

#define COMMAND_CLASS
#define CommandStyle(key,Class) \
  (*command_map)[#key] = &command_creator<Class>;
#include "style_command.h"      // IWYU pragma: keep
#undef CommandStyle
#undef COMMAND_CLASS

  // process command-line args
  // check for args "-var" and "-echo"
  // caller has already checked that sufficient arguments exist

  int iarg = 1;
  while (iarg < argc) {
    if (strcmp(argv[iarg],"-var") == 0 || strcmp(argv[iarg],"-v") == 0) {
      int jarg = iarg+3;
      while (jarg < argc && argv[jarg][0] != '-') jarg++;
      variable->set(argv[iarg+1],jarg-iarg-2,&argv[iarg+2]);
      iarg = jarg;
    } else if (strcmp(argv[iarg],"-echo") == 0 ||
               strcmp(argv[iarg],"-e") == 0) {
      narg = 1;
      char **tmp = arg;        // trick echo() into using argv instead of arg
      arg = &argv[iarg+1];
      echo();
      arg = tmp;
      iarg += 2;
     } else iarg++;
  }
}

/* ---------------------------------------------------------------------- */

Input::~Input()
{
  // don't free command and arg strings
  // they just point to other allocated memory

  memory->sfree(line);
  memory->sfree(copy);
  memory->sfree(work);
  delete[] labelstr;
  memory->sfree(arg);
  delete[] infiles;
  delete[] inlines;
  delete variable;

  delete command_map;
}

/** Process all input from the ``FILE *`` pointer *infile*
 *
\verbatim embed:rst

This will read lines from *infile*, parse and execute them until the end
of the file is reached.  The *infile* pointer will usually point to
``stdin`` or the input file given with the ``-in`` command line flag.

\endverbatim */

void Input::file()
{
  int m,n,mstart,ntriple,endfile;
  int nline = *output->thermo->get_line();

  while (true) {

    // read a line from input script
    // when done, n = length of line including str terminator, 0 if end of file
    // if line ends in continuation char '&', concatenate next line
    // if triple quotes are used, read until closing triple quotes

    if (me == 0) {
      ntriple = 0;
      endfile = 0;
      m = 0;

      while (true) {

        if (infile == nullptr) {
          n = 0;
          break;
        }

        mstart = m;

        while (true) {
          if (maxline-m < 2) reallocate(line,maxline,0);

          // end of file reached, so break
          // n == 0 if nothing read, else n = line with str terminator

          if (fgets(&line[m],maxline-m,infile) == nullptr) {
            endfile = 1;
            if (m) n = strlen(line) + 1;
            else n = 0;
            break;
          }

          // continue if last char read was not a newline
          // can happen if line is very long

          m += strlen(&line[m]);
          if (line[m-1] != '\n') continue;
          break;
        }

        if (endfile) break;

        // add # of triple quotes in just-read line to ntriple

        ntriple += numtriple(&line[mstart]);

        // trim whitespace from end of line
        // line[m] = last printable char

        m--;
        while (m >= 0 && isspace(line[m])) m--;

        // continue reading if final printable char is "&", count line

        if (m >= 0 && line[m] == '&') {
          ++nline;
          continue;
        }

        // continue reading if odd number of triple quotes

        if (ntriple % 2) {
          line[m+1] = '\n';
          m += 2;
          ++nline;
          continue;
        }

        // done, break with n = length of line with str terminator

        line[m+1] = '\0';
        n = m+2;
        break;
      }
    }
    output->thermo->set_line(++nline);

    // bcast the line
    // if n = 0, end-of-file
    // error if label_active is set, since label wasn't encountered
    // if original input file, code is done
    // else go back to previous input file

    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n == 0) {
      if (label_active) error->all(FLERR,"Label wasn't found in input script");
      break;
    }

    if (n > maxline) reallocate(line,maxline,n);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // echo the command unless scanning for label

    if (me == 0 && label_active == 0) {
      if (echo_screen && screen) fprintf(screen,"%s\n",line);
      if (echo_log && logfile) fprintf(logfile,"%s\n",line);
    }

    // parse the line
    // if no command, skip to next line in input script

    parse();
    if (command == nullptr) continue;

    // if scanning for label, skip command unless it's a label command

    if (label_active && strcmp(command,"label") != 0) continue;

    // execute the command

    if (execute_command() && line)
      error->all(FLERR,"Unknown command: {}",line);
    nline = *output->thermo->get_line();
  }
}

/** Process all input from the file *filename*
 *
\verbatim embed:rst

This function opens the file at the path *filename*, put the current
file pointer stored in *infile* on a stack and instead assign *infile*
with the newly opened file pointer.  Then it will call the
:cpp:func:`Input::file() <LAMMPS_NS::Input::file()>` function to read,
parse and execute the contents of that file.  When the end of the file
is reached, it is closed and the previous file pointer from the infile
file pointer stack restored to *infile*.

\endverbatim
 *
 * \param  filename  name of file with LAMMPS commands */

void Input::file(const char *filename)
{
  // error if another nested file still open, should not be possible
  // open new filename and set infile, infiles[0], nfile
  // call to file() will close filename and decrement nfile

  if (me == 0) {
    if (nfile == LMP_MAXFILE) error->one(FLERR,"Too many nested levels of input scripts");

    if (filename) {
      infile = fopen(filename,"r");
      if (infile == nullptr)
        error->one(FLERR,"Cannot open input script {}: {}", filename, utils::getsyserror());
      if (nfile > 0) inlines[nfile - 1] = *output->thermo->get_line();
      inlines[nfile] = -1;
      infiles[nfile++] = infile;
    }
  }

  // process contents of file

  file();

  if (me == 0) {
    if (filename) {
      fclose(infile);
      nfile--;
      infile = infiles[nfile-1];
      output->thermo->set_line(inlines[nfile-1]);
    }
  }
}

/** Process a single command from a string in *single*
 *
\verbatim embed:rst

This function takes the text in *single*, makes a copy, parses that,
executes the command and returns the name of the command (without the
arguments).  If there was no command in *single* it will return
``nullptr``.

\endverbatim
 *
 * \param  single  string with LAMMPS command
 * \return         string with name of the parsed command w/o arguments */

char *Input::one(const std::string &single)
{
  int n = single.size() + 1;
  if (n > maxline) reallocate(line,maxline,n);
  strcpy(line,single.c_str());

  // echo the command unless scanning for label

  if (me == 0 && label_active == 0) {
    if (echo_screen && screen) fprintf(screen,"%s\n",line);
    if (echo_log && logfile) fprintf(logfile,"%s\n",line);
  }

  // parse the line
  // if no command, just return a null pointer

  parse();
  if (command == nullptr) return nullptr;

  // if scanning for label, skip command unless it's a label command

  if (label_active && strcmp(command,"label") != 0) return nullptr;

  // execute the command and return its name

  if (execute_command())
    error->all(FLERR,"Unknown command: {}",line);

  return command;
}

/* ----------------------------------------------------------------------
   send text to active echo file pointers
------------------------------------------------------------------------- */

void Input::write_echo(const std::string &txt)
{
  if (me == 0) {
    if (echo_screen && screen) fputs(txt.c_str(),screen);
    if (echo_log && logfile) fputs(txt.c_str(),logfile);
  }
}

/* ----------------------------------------------------------------------
   parse copy of command line by inserting string terminators
   strip comment = all chars from # on
   replace all $ via variable substitution except within quotes
   command = first word
   narg = # of args
   arg[] = individual args
   treat text between single/double/triple quotes as one arg via nextword()
------------------------------------------------------------------------- */

void Input::parse()
{
  // duplicate line into copy string to break into words

  int n = strlen(line) + 1;
  if (n > maxcopy) reallocate(copy,maxcopy,n);
  strcpy(copy,line);

  // strip a # comment by replacing it with 0
  // do not treat a # inside single/double/triple quotes as a comment

  char *ptrmatch;
  char *ptr = copy;

  while (*ptr) {
    if (*ptr == '#') {
      *ptr = '\0';
      break;
    }
    if (*ptr == '\'') {
      ptrmatch = strchr(ptr+1,'\'');
      if (ptrmatch == nullptr)
        error->all(FLERR,"Unmatched single quote in command");
      ptr = ptrmatch + 1;
    } else if (*ptr == '"') {
      if (strstr(ptr,"\"\"\"") == ptr) {
        ptrmatch = strstr(ptr+3,"\"\"\"");
        if (ptrmatch == nullptr)
          error->all(FLERR,"Unmatched triple quote in command");
        ptr = ptrmatch + 3;
      } else {
        ptrmatch = strchr(ptr+1,'"');
        if (ptrmatch == nullptr)
          error->all(FLERR,"Unmatched double quote in command");
        ptr = ptrmatch + 1;
      }
    } else ptr++;
  }

  if (utils::has_utf8(copy)) {
    std::string buf = utils::utf8_subst(copy);
    strcpy(copy,buf.c_str());
    if (utf8_warn && (comm->me == 0))
      error->warning(FLERR,"Detected non-ASCII characters in input. "
                     "Will try to continue by replacing with ASCII "
                     "equivalents where known.");
    utf8_warn = false;
  }

  // perform $ variable substitution (print changes)
  // except if searching for a label since earlier variable may not be defined

  if (!label_active) substitute(copy,work,maxcopy,maxwork,1);

  // command = 1st arg in copy string

  char *next;
  command = nextword(copy,&next);
  if (command == nullptr) return;

  // point arg[] at each subsequent arg in copy string
  // nextword() inserts string terminators into copy string to delimit args
  // nextword() treats text between single/double/triple quotes as one arg

  narg = 0;
  ptr = next;
  while (ptr) {
    if (narg == maxarg) {
      maxarg += DELTA;
      arg = (char **) memory->srealloc(arg,maxarg*sizeof(char *),"input:arg");
    }
    arg[narg] = nextword(ptr,&next);
    if (!arg[narg]) break;
    narg++;
    ptr = next;
  }
}

/* ----------------------------------------------------------------------
   find next word in str
   insert 0 at end of word
   ignore leading whitespace
   treat text between single/double/triple quotes as one arg
   matching quote must be followed by whitespace char if not end of string
   strip quotes from returned word
   return ptr to start of word or null pointer if no word in string
   also return next = ptr after word
------------------------------------------------------------------------- */

char *Input::nextword(char *str, char **next)
{
  char *start,*stop;

  // start = first non-whitespace char

  start = &str[strspn(str," \t\n\v\f\r")];
  if (*start == '\0') return nullptr;

  // if start is single/double/triple quote:
  //   start = first char beyond quote
  //   stop = first char of matching quote
  //   next = first char beyond matching quote
  //   next must be null char or whitespace
  // if start is not single/double/triple quote:
  //   stop = first whitespace char after start
  //   next = char after stop, or stop itself if stop is null char

  if (strstr(start,"\"\"\"") == start) {
    stop = strstr(&start[3],"\"\"\"");
    if (!stop) error->all(FLERR,"Unbalanced quotes in input line");
    start += 3;
    *next = stop+3;
    if (**next && !isspace(**next))
      error->all(FLERR,"Input line quote not followed by white-space");
  } else if (*start == '"' || *start == '\'') {
    stop = strchr(&start[1],*start);
    if (!stop) error->all(FLERR,"Unbalanced quotes in input line");
    start++;
    *next = stop+1;
    if (**next && !isspace(**next))
      error->all(FLERR,"Input line quote not followed by white-space");
  } else {
    stop = &start[strcspn(start," \t\n\v\f\r")];
    if (*stop == '\0') *next = stop;
    else *next = stop+1;
  }

  // set stop to null char to terminate word

  *stop = '\0';
  return start;
}

/* ----------------------------------------------------------------------
   substitute for $ variables in str using work str2 and return it
   reallocate str/str2 to hold expanded version if necessary & reset max/max2
   print updated string if flag is set and not searching for label
   label_active will be 0 if called from external class
------------------------------------------------------------------------- */

void Input::substitute(char *&str, char *&str2, int &max, int &max2, int flag)
{
  // use str2 as scratch space to expand str, then copy back to str
  // reallocate str and str2 as necessary
  // do not replace variables inside single/double/triple quotes
  // var = pts at variable name, ended by null char
  //   if $ is followed by '{', trailing '}' becomes null char
  //   else $x becomes x followed by null char
  // beyond = points to text following variable

  int i,n,paren_count,nchars;
  char immediate[256];
  char *var,*value,*beyond;
  int quoteflag = 0;
  char *ptrmatch;

  char *ptr = str;

  n = strlen(str) + 1;
  if (n > max2) reallocate(str2,max2,n);
  *str2 = '\0';
  char *ptr2 = str2;

  while (*ptr) {

    // variable substitution

    if (*ptr == '$' && !quoteflag) {

      // value = ptr to expanded variable
      // variable name between curly braces, e.g. ${a}

      if (*(ptr+1) == '{') {
        var = ptr+2;
        i = 0;

        while (var[i] != '\0' && var[i] != '}') i++;

        if (var[i] == '\0') error->one(FLERR,"Invalid variable name");
        var[i] = '\0';
        beyond = ptr + strlen(var) + 3;
        value = variable->retrieve(var);

      // immediate variable between parenthesis, e.g. $(1/3) or $(1/3:%.6g)

      } else if (*(ptr+1) == '(') {
        var = ptr+2;
        paren_count = 0;
        i = 0;

        while (var[i] != '\0' && (var[i] != ')' || paren_count != 0)) {
          switch (var[i]) {
          case '(': paren_count++; break;
          case ')': paren_count--; break;
          default: ;
          }
          i++;
        }

        if (var[i] == '\0') error->one(FLERR,"Invalid immediate variable");
        var[i] = '\0';
        beyond = ptr + strlen(var) + 3;

        // check if an inline format specifier was appended with a colon

        char fmtstr[64] = "%.20g";
        char *fmtflag;
        if ((fmtflag=strrchr(var, ':')) && (fmtflag[1]=='%')) {
          strncpy(fmtstr,&fmtflag[1],sizeof(fmtstr)-1);
          *fmtflag='\0';
        }

        // quick check for proper format string

        if (!utils::strmatch(fmtstr,"%[0-9 ]*\\.[0-9]+[efgEFG]"))
          error->all(FLERR,"Incorrect conversion in format string");

        snprintf(immediate,256,fmtstr,variable->compute_equal(var));
        value = immediate;

      // single character variable name, e.g. $a

      } else {
        var = ptr;
        var[0] = var[1];
        var[1] = '\0';
        beyond = ptr + 2;
        value = variable->retrieve(var);
      }

      if (value == nullptr)
        error->one(FLERR,"Substitution for illegal variable {}",var);

      // check if storage in str2 needs to be expanded
      // re-initialize ptr and ptr2 to the point beyond the variable.

      n = strlen(str2) + strlen(value) + strlen(beyond) + 1;
      if (n > max2) reallocate(str2,max2,n);
      strcat(str2,value);
      ptr2 = str2 + strlen(str2);
      ptr = beyond;

      // output substitution progress if requested

      if (flag && me == 0 && label_active == 0) {
        if (echo_screen && screen) fprintf(screen,"%s%s\n",str2,beyond);
        if (echo_log && logfile) fprintf(logfile,"%s%s\n",str2,beyond);
      }

    // check for single/double/triple quotes and skip past them

    } else if (*ptr == '\'') {
      ptrmatch = strchr(ptr+1,'\'');
      if (ptrmatch == nullptr)
        error->all(FLERR,"Unmatched single quote in command");
      nchars = ptrmatch+1 - ptr;
      strncpy(ptr2,ptr,nchars);
      ptr += nchars;
      ptr2 += nchars;
    } else if (*ptr == '"') {
      if (strstr(ptr,"\"\"\"") == ptr) {
        ptrmatch = strstr(ptr+3,"\"\"\"");
        if (ptrmatch == nullptr)
          error->all(FLERR,"Unmatched triple quote in command");
        nchars = ptrmatch+3 - ptr;
        strncpy(ptr2,ptr,nchars);
        ptr += nchars;
        ptr2 += nchars;
      } else {
        ptrmatch = strchr(ptr+1,'"');
        if (ptrmatch == nullptr)
          error->all(FLERR,"Unmatched double quote in command");
        nchars = ptrmatch+1 - ptr;
        strncpy(ptr2,ptr,nchars);
        ptr += nchars;
        ptr2 += nchars;
      }

    // else copy current single character into str2

    } else *ptr2++ = *ptr++;

    // terminate current str2 so variable sub can perform strlen()

    *ptr2 = '\0';
  }

  // set length of input str to length of work str2
  // copy work string back to input str

  if (max2 > max) reallocate(str,max,max2);
  strcpy(str,str2);
}

/* ----------------------------------------------------------------------
   return number of triple quotes in line
------------------------------------------------------------------------- */

int Input::numtriple(char *line)
{
  int count = 0;
  char *ptr = line;
  while ((ptr = strstr(ptr,"\"\"\""))) {
    ptr += 3;
    count++;
  }
  return count;
}

/* ----------------------------------------------------------------------
   rellocate a string
   if n > 0: set max >= n in increments of DELTALINE
   if n = 0: just increment max by DELTALINE
------------------------------------------------------------------------- */

void Input::reallocate(char *&str, int &max, int n)
{
  if (n) {
    while (n > max) max += DELTALINE;
  } else max += DELTALINE;

  str = (char *) memory->srealloc(str,max*sizeof(char),"input:str");
}

/* ----------------------------------------------------------------------
   process a single parsed command
   return 0 if successful, -1 if did not recognize command
------------------------------------------------------------------------- */

int Input::execute_command()
{
  int flag = 1;

  std::string mycmd = command;
  if (mycmd == "clear") clear();
  else if (mycmd == "echo") echo();
  else if (mycmd == "if") ifthenelse();
  else if (mycmd == "include") include();
  else if (mycmd == "jump") jump();
  else if (mycmd == "label") label();
  else if (mycmd == "log") log();
  else if (mycmd == "next") next_command();
  else if (mycmd == "partition") partition();
  else if (mycmd == "print") print();
  else if (mycmd == "python") python();
  else if (mycmd == "quit") quit();
  else if (mycmd == "shell") shell();
  else if (mycmd == "variable") variable_command();

  else if (mycmd == "angle_coeff") angle_coeff();
  else if (mycmd == "angle_style") angle_style();
  else if (mycmd == "atom_modify") atom_modify();
  else if (mycmd == "atom_style") atom_style();
  else if (mycmd == "bond_coeff") bond_coeff();
  else if (mycmd == "bond_style") bond_style();
  else if (mycmd == "bond_write") bond_write();
  else if (mycmd == "boundary") boundary();
  else if (mycmd == "comm_modify") comm_modify();
  else if (mycmd == "comm_style") comm_style();
  else if (mycmd == "compute") compute();
  else if (mycmd == "compute_modify") compute_modify();
  else if (mycmd == "dielectric") dielectric();
  else if (mycmd == "dihedral_coeff") dihedral_coeff();
  else if (mycmd == "dihedral_style") dihedral_style();
  else if (mycmd == "dimension") dimension();
  else if (mycmd == "dump") dump();
  else if (mycmd == "dump_modify") dump_modify();
  else if (mycmd == "fix") fix();
  else if (mycmd == "fix_modify") fix_modify();
  else if (mycmd == "group") group_command();
  else if (mycmd == "improper_coeff") improper_coeff();
  else if (mycmd == "improper_style") improper_style();
  else if (mycmd == "kspace_modify") kspace_modify();
  else if (mycmd == "kspace_style") kspace_style();
  else if (mycmd == "labelmap") labelmap();
  else if (mycmd == "lattice") lattice();
  else if (mycmd == "mass") mass();
  else if (mycmd == "min_modify") min_modify();
  else if (mycmd == "min_style") min_style();
  else if (mycmd == "molecule") molecule();
  else if (mycmd == "neigh_modify") neigh_modify();
  else if (mycmd == "neighbor") neighbor_command();
  else if (mycmd == "newton") newton();
  else if (mycmd == "package") package();
  else if (mycmd == "pair_coeff") pair_coeff();
  else if (mycmd == "pair_modify") pair_modify();
  else if (mycmd == "pair_style") pair_style();
  else if (mycmd == "pair_write") pair_write();
  else if (mycmd == "processors") processors();
  else if (mycmd == "region") region();
  else if (mycmd == "reset_timestep") reset_timestep();
  else if (mycmd == "restart") restart();
  else if (mycmd == "run_style") run_style();
  else if (mycmd == "special_bonds") special_bonds();
  else if (mycmd == "suffix") suffix();
  else if (mycmd == "thermo") thermo();
  else if (mycmd == "thermo_modify") thermo_modify();
  else if (mycmd == "thermo_style") thermo_style();
  else if (mycmd == "timestep") timestep();
  else if (mycmd == "timer") timer_command();
  else if (mycmd == "uncompute") uncompute();
  else if (mycmd == "undump") undump();
  else if (mycmd == "unfix") unfix();
  else if (mycmd == "units") units();

  else flag = 0;

  // return if command was listed above

  if (flag) return 0;

  // process "meta-commands", i.e. commands that may have sub-commands
  // they return 1 if there was a match and 0 if not

  if (mycmd == "reset_atoms") flag = meta(mycmd);
  if (flag) return 0;

  // invoke commands added via style_command.h
  // try suffixed version first

  if (lmp->suffix_enable && lmp->non_pair_suffix()) {
    mycmd = command + std::string("/") + lmp->non_pair_suffix();
    if (command_map->find(mycmd) == command_map->end()) {
      if (lmp->suffix2) {
        mycmd = command + std::string("/") + lmp->suffix2;
        if (command_map->find(mycmd) == command_map->end())
          mycmd = command;
      } else mycmd = command;
    }
  }
  if (command_map->find(mycmd) != command_map->end()) {
    CommandCreator &command_creator = (*command_map)[mycmd];
    Command *cmd = command_creator(lmp);
    cmd->command(narg,arg);
    delete cmd;
    return 0;
  }

  // unrecognized command

  return -1;
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void Input::clear()
{
  if (narg > 0) error->all(FLERR,"Illegal clear command: unexpected arguments but found {}", narg);
  output->thermo->set_line(-1);
  lmp->destroy();
  lmp->create();
  lmp->post_create();
}

/* ---------------------------------------------------------------------- */

void Input::echo()
{
  if (narg != 1) error->all(FLERR,"Illegal echo command: expected 1 argument but found {}", narg);

  if (strcmp(arg[0],"none") == 0) {
    echo_screen = 0;
    echo_log = 0;
  } else if (strcmp(arg[0],"screen") == 0) {
    echo_screen = 1;
    echo_log = 0;
  } else if (strcmp(arg[0],"log") == 0) {
    echo_screen = 0;
    echo_log = 1;
  } else if (strcmp(arg[0],"both") == 0) {
    echo_screen = 1;
    echo_log = 1;
  } else error->all(FLERR,"Unknown echo keyword: {}", arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::ifthenelse()
{
  if (narg < 3) utils::missing_cmd_args(FLERR, "if", error);

  // substitute for variables in Boolean expression for "if"
  // in case expression was enclosed in quotes
  // must substitute on copy of arg else will step on subsequent args

  int n = strlen(arg[0]) + 1;
  if (n > maxline) reallocate(line,maxline,n);
  strcpy(line,arg[0]);
  substitute(line,work,maxline,maxwork,0);

  // evaluate Boolean expression for "if"

  double btest = variable->evaluate_boolean(line);

  // bound "then" commands

  if (strcmp(arg[1],"then") != 0) error->all(FLERR,"Illegal if command: expected \"then\" but found \"{}\"", arg[1]);

  int first = 2;
  int iarg = first;
  while (iarg < narg &&
         (strcmp(arg[iarg],"elif") != 0 && strcmp(arg[iarg],"else") != 0))
    iarg++;
  int last = iarg-1;

  // execute "then" commands
  // make copies of all arg string commands
  // required because re-parsing a command via one() will wipe out args

  if (btest != 0.0) {
    int ncommands = last-first + 1;
    if (ncommands <= 0) utils::missing_cmd_args(FLERR, "if then", error);

    auto commands = new char*[ncommands];
    ncommands = 0;
    for (int i = first; i <= last; i++) {
      n = strlen(arg[i]) + 1;
      if (n == 1) error->all(FLERR,"Illegal if then command: execute command is empty");
      commands[ncommands] = new char[n];
      strcpy(commands[ncommands],arg[i]);
      ncommands++;
    }

    for (int i = 0; i < ncommands; i++) {
      one(commands[i]);
      delete[] commands[i];
    }
    delete[] commands;

    return;
  }

  // done if no "elif" or "else"

  if (iarg == narg) return;

  // check "elif" or "else" until find commands to execute
  // substitute for variables and evaluate Boolean expression for "elif"
  // must substitute on copy of arg else will step on subsequent args
  // bound and execute "elif" or "else" commands

  while (iarg != narg) {
    if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "if then", error);
    if (strcmp(arg[iarg],"elif") == 0) {
      n = strlen(arg[iarg+1]) + 1;
      if (n > maxline) reallocate(line,maxline,n);
      strcpy(line,arg[iarg+1]);
      substitute(line,work,maxline,maxwork,0);
      btest = variable->evaluate_boolean(line);
      first = iarg+2;
    } else {
      btest = 1.0;
      first = iarg+1;
    }

    iarg = first;
    while (iarg < narg &&
           (strcmp(arg[iarg],"elif") != 0 && strcmp(arg[iarg],"else") != 0))
      iarg++;
    last = iarg-1;

    if (btest == 0.0) continue;

    int ncommands = last-first + 1;
    if (ncommands <= 0) utils::missing_cmd_args(FLERR, "if elif/else", error);

    auto commands = new char*[ncommands];
    ncommands = 0;
    for (int i = first; i <= last; i++) {
      n = strlen(arg[i]) + 1;
      if (n == 1) error->all(FLERR,"Illegal if elif/else command: execute command is empty");
      commands[ncommands] = new char[n];
      strcpy(commands[ncommands],arg[i]);
      ncommands++;
    }

    // execute the list of commands

    for (int i = 0; i < ncommands; i++) {
      one(commands[i]);
      delete[] commands[i];
    }
    delete[] commands;

    return;
  }
}

/* ---------------------------------------------------------------------- */

void Input::include()
{
  if (narg != 1) error->all(FLERR,"Illegal include command");

  if (me == 0) {
    if (nfile == LMP_MAXFILE)
      error->one(FLERR,"Too many nested levels of input scripts");

    // expand variables
    int n = strlen(arg[0]) + 1;
    if (n > maxline) reallocate(line,maxline,n);
    strcpy(line,arg[0]);
    substitute(line,work,maxline,maxwork,0);

    infile = fopen(line,"r");
    if (infile == nullptr)
      error->one(FLERR,"Cannot open input script {}: {}", line, utils::getsyserror());

    infiles[nfile++] = infile;
  }

  // process contents of file

  file();

  if (me == 0) {
    fclose(infile);
    nfile--;
    infile = infiles[nfile-1];
  }
}

/* ---------------------------------------------------------------------- */

void Input::jump()
{
  if (narg < 1 || narg > 2) error->all(FLERR,"Illegal jump command: expected 1 or 2 argument(s) but found {}", narg);

  if (jump_skip) {
    jump_skip = 0;
    return;
  }

  if (me == 0) {
    output->thermo->set_line(-1);
    if (strcmp(arg[0],"SELF") == 0) {
      rewind(infile);
    } else {
      if (infile && infile != stdin) fclose(infile);
      infile = fopen(arg[0],"r");
      if (infile == nullptr)
        error->one(FLERR,"Cannot open input script {}: {}", arg[0], utils::getsyserror());
      inlines[nfile-1] = -1;
      infiles[nfile-1] = infile;
    }
  }

  if (narg == 2) {
    label_active = 1;
    delete[] labelstr;
    labelstr = utils::strdup(arg[1]);
  }
}

/* ---------------------------------------------------------------------- */

void Input::label()
{
  if (narg != 1) error->all(FLERR,"Illegal label command: expected 1 argument but found {}", narg);
  if (label_active && strcmp(labelstr,arg[0]) == 0) label_active = 0;
}

/* ---------------------------------------------------------------------- */

void Input::log()
{
  if ((narg < 1) || (narg > 2)) error->all(FLERR,"Illegal log command: expected 1 or 2 argument(s) but found {}", narg);

  int appendflag = 0;
  if (narg == 2) {
    if (strcmp(arg[1],"append") == 0) appendflag = 1;
    else error->all(FLERR,"Unknown log keyword: {}", arg[1]);
  }

  if (me == 0) {
    if (logfile) fclose(logfile);
    if (strcmp(arg[0],"none") == 0) logfile = nullptr;
    else {
      if (appendflag) logfile = fopen(arg[0],"a");
      else logfile = fopen(arg[0],"w");

      if (logfile == nullptr)
        error->one(FLERR,"Cannot open logfile {}: {}",
                                     arg[0], utils::getsyserror());

    }
    if (universe->nworlds == 1) universe->ulogfile = logfile;
  }
}

/* ---------------------------------------------------------------------- */

void Input::next_command()
{
  if (variable->next(narg,arg)) jump_skip = 1;
}

/* ---------------------------------------------------------------------- */

void Input::partition()
{
  if (narg < 3) utils::missing_cmd_args(FLERR, "partition", error);

  int ilo,ihi;
  int yesflag = utils::logical(FLERR,arg[0],false,lmp);
  utils::bounds(FLERR,arg[1],1,universe->nworlds,ilo,ihi,error);

  // new command starts at the 3rd argument,
  // which must not be another partition command

  if (strcmp(arg[2],"partition") == 0) error->all(FLERR,"Illegal partition command");

  char *cmd = strstr(line,arg[2]);

  // execute the remaining command line on requested partitions

  if (yesflag) {
    if (universe->iworld+1 >= ilo && universe->iworld+1 <= ihi) one(cmd);
  } else {
    if (universe->iworld+1 < ilo || universe->iworld+1 > ihi) one(cmd);
  }
}

/* ---------------------------------------------------------------------- */

void Input::print()
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "print", error);

  // copy 1st arg back into line (copy is being used)
  // check maxline since arg[0] could have been expanded by variables
  // substitute for $ variables (no printing) and print arg

  int n = strlen(arg[0]) + 1;
  if (n > maxline) reallocate(line,maxline,n);
  strcpy(line,arg[0]);
  substitute(line,work,maxline,maxwork,0);

  // parse optional args

  FILE *fp = nullptr;
  int screenflag = 1;
  int universeflag = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0 || strcmp(arg[iarg],"append") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal print {} command: missing argument(s)", arg[iarg]);
      if (me == 0) {
        if (fp != nullptr) fclose(fp);
        if (strcmp(arg[iarg],"file") == 0) fp = fopen(arg[iarg+1],"w");
        else fp = fopen(arg[iarg+1],"a");
        if (fp == nullptr)
          error->one(FLERR,"Cannot open print file {}: {}", arg[iarg+1], utils::getsyserror());
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"screen") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "print screen", error);
      screenflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"universe") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "print universe", error);
      universeflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Unknown print keyword: {}", arg[iarg]);
  }

  if (me == 0) {
    if (screenflag && screen) fprintf(screen,"%s\n",line);
    if (screenflag && logfile) fprintf(logfile,"%s\n",line);
    if (fp) {
      fprintf(fp,"%s\n",line);
      fclose(fp);
    }
  }
  if (universeflag && (universe->me == 0)) {
    if (universe->uscreen)  fprintf(universe->uscreen, "%s\n",line);
    if (universe->ulogfile) fprintf(universe->ulogfile,"%s\n",line);
  }
}

/* ---------------------------------------------------------------------- */

void Input::python()
{
  variable->python_command(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::quit()
{
  if (narg == 0) error->done(0); // 1 would be fully backwards compatible
  if (narg == 1) error->done(utils::inumeric(FLERR,arg[0],false,lmp));
  error->all(FLERR,"Illegal quit command: expected 0 or 1 argument but found {}", narg);
}

/* ---------------------------------------------------------------------- */

void Input::shell()
{
  int rv,err;

  if (narg < 1) utils::missing_cmd_args(FLERR, "shell", error);

  if (strcmp(arg[0],"cd") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal shell command: expected 2 argument but found {}", narg);
    rv = (platform::chdir(arg[1]) < 0) ? errno : 0;
    MPI_Reduce(&rv,&err,1,MPI_INT,MPI_MAX,0,world);
    errno = err;
    if (me == 0 && err != 0) {
      error->warning(FLERR, "Shell command 'cd {}' failed with error '{}'", arg[1], utils::getsyserror());
    }
  } else if (strcmp(arg[0],"mkdir") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "shell mkdir", error);
    if (me == 0) {
      for (int i = 1; i < narg; i++) {
        rv = (platform::mkdir(arg[i]) < 0) ? errno : 0;
        if (rv != 0)
          error->warning(FLERR, "Shell command 'mkdir {}' failed with error '{}'", arg[i],
                         utils::getsyserror());
      }
    }
  } else if (strcmp(arg[0],"mv") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal shell command: expected 3 argument but found {}", narg);
    if (me == 0) {
      if (platform::path_is_directory(arg[2])) {
        if (system(fmt::format("mv {} {}", arg[1], arg[2]).c_str()))
          error->warning(FLERR,"Shell command 'mv {} {}' returned with non-zero status", arg[1], arg[2]);
      } else {
        if (rename(arg[1],arg[2]) < 0) {
          error->warning(FLERR, "Shell command 'mv {} {}' failed with error '{}'",
                         arg[1],arg[2],utils::getsyserror());
        }
      }
    }
  } else if (strcmp(arg[0],"rm") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "shell rm", error);
    if (me == 0) {
      int i = 1;
      bool warn = true;
      if (strcmp(arg[i], "-f") == 0) {
        warn = false;
        ++i;
      }
      for (;i < narg; i++) {
        if (platform::unlink(arg[i]) < 0)
          if (warn)
            error->warning(FLERR, "Shell command 'rm {}' failed with error '{}'",
                           arg[i], utils::getsyserror());
      }
    }
  } else if (strcmp(arg[0],"rmdir") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "shell rmdir", error);
    if (me == 0) {
      for (int i = 1; i < narg; i++) {
        if (platform::rmdir(arg[i]) < 0)
          error->warning(FLERR, "Shell command 'rmdir {}' failed with error '{}'",
                         arg[i], utils::getsyserror());
      }
    }
  } else if (strcmp(arg[0],"putenv") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "shell putenv", error);
    for (int i = 1; i < narg; i++) {
      rv = 0;
      if (arg[i]) rv = platform::putenv(arg[i]);
      rv = (rv < 0) ? errno : 0;
      MPI_Reduce(&rv,&err,1,MPI_INT,MPI_MAX,0,world);
      errno = err;
      if (me == 0 && err != 0)
        error->warning(FLERR, "Shell command 'putenv {}' failed with error '{}'",
                       arg[i], utils::getsyserror());
    }

  // concat arguments and invoke string in shell via system()

  } else {
    if (me == 0) {
      std::string cmd = arg[0];
      for (int i = 1; i < narg; i++) {
        cmd += " ";
        cmd += arg[i];
      }

      if (system(cmd.c_str()) != 0)
        error->warning(FLERR,"Shell command {} returned with non-zero status", cmd);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Input::variable_command()
{
  variable->set(narg,arg);
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   one function for each LAMMPS-specific input script command
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void Input::angle_coeff()
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Angle_coeff command before simulation box is defined");
  if (force->angle == nullptr)
    error->all(FLERR,"Angle_coeff command before angle_style is defined");
  if (atom->avec->angles_allow == 0)
    error->all(FLERR,"Angle_coeff command when no angles allowed");
  char *newarg = utils::expand_type(FLERR, arg[0], Atom::ANGLE, lmp);
  if (newarg) arg[0] = newarg;
  force->angle->coeff(narg,arg);
  delete[] newarg;
}

/* ---------------------------------------------------------------------- */

void Input::angle_style()
{
  if (narg < 1) error->all(FLERR,"Illegal angle_style command");
  if (atom->avec->angles_allow == 0)
    error->all(FLERR,"Angle_style command when no angles allowed");
  force->create_angle(arg[0],1);
  if (force->angle) force->angle->settings(narg-1,&arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::atom_modify()
{
  atom->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::atom_style()
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "atom_style", error);
  if (domain->box_exist)
    error->all(FLERR,"Atom_style command after simulation box is defined");
  atom->create_avec(arg[0],narg-1,&arg[1],1);
}

/* ---------------------------------------------------------------------- */

void Input::bond_coeff()
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Bond_coeff command before simulation box is defined");
  if (force->bond == nullptr)
    error->all(FLERR,"Bond_coeff command before bond_style is defined");
  if (atom->avec->bonds_allow == 0)
    error->all(FLERR,"Bond_coeff command when no bonds allowed");
  char *newarg = utils::expand_type(FLERR, arg[0], Atom::BOND, lmp);
  if (newarg) arg[0] = newarg;
  force->bond->coeff(narg,arg);
  delete[] newarg;
}

/* ---------------------------------------------------------------------- */

void Input::bond_style()
{
  if (narg < 1) error->all(FLERR,"Illegal bond_style command");
  if (atom->avec->bonds_allow == 0)
    error->all(FLERR,"Bond_style command when no bonds allowed");
  force->create_bond(arg[0],1);
  if (force->bond) force->bond->settings(narg-1,&arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::bond_write()
{
  if (atom->avec->bonds_allow == 0)
    error->all(FLERR,"Bond_write command when no bonds allowed");
  if (force->bond == nullptr)
    error->all(FLERR,"Bond_write command before bond_style is defined");
  else force->bond->write_file(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::boundary()
{
  if (domain->box_exist)
    error->all(FLERR,"Boundary command after simulation box is defined");
  domain->set_boundary(narg,arg,0);
}

/* ---------------------------------------------------------------------- */

void Input::comm_modify()
{
  comm->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::comm_style()
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "comm_style", error);
  if (strcmp(arg[0],"brick") == 0) {
    if (comm->style == Comm::BRICK) return;
    Comm *oldcomm = comm;
    comm = new CommBrick(lmp,oldcomm);
    delete oldcomm;
  } else if (strcmp(arg[0],"tiled") == 0) {
    if (comm->style == Comm::TILED) return;
    Comm *oldcomm = comm;
    if (lmp->kokkos) comm = new CommTiledKokkos(lmp,oldcomm);
    else comm = new CommTiled(lmp,oldcomm);
    delete oldcomm;
  } else error->all(FLERR,"Unknown comm_style argument: {}", arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::compute()
{
  modify->add_compute(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::compute_modify()
{
  modify->modify_compute(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::dielectric()
{
  if (narg != 1) error->all(FLERR,"Illegal dielectric command");
  force->dielectric = utils::numeric(FLERR,arg[0],false,lmp);
}

/* ---------------------------------------------------------------------- */

void Input::dihedral_coeff()
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Dihedral_coeff command before simulation box is defined");
  if (force->dihedral == nullptr)
    error->all(FLERR,"Dihedral_coeff command before dihedral_style is defined");
  if (atom->avec->dihedrals_allow == 0)
    error->all(FLERR,"Dihedral_coeff command when no dihedrals allowed");
  char *newarg = utils::expand_type(FLERR, arg[0], Atom::DIHEDRAL, lmp);
  if (newarg) arg[0] = newarg;
  force->dihedral->coeff(narg,arg);
  delete[] newarg;
}

/* ---------------------------------------------------------------------- */

void Input::dihedral_style()
{
  if (narg < 1) error->all(FLERR,"Illegal dihedral_style command");
  if (atom->avec->dihedrals_allow == 0)
    error->all(FLERR,"Dihedral_style command when no dihedrals allowed");
  force->create_dihedral(arg[0],1);
  if (force->dihedral) force->dihedral->settings(narg-1,&arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::dimension()
{
  if (narg != 1) error->all(FLERR, "Dimension command expects exactly 1 argument");
  if (domain->box_exist)
    error->all(FLERR,"Dimension command after simulation box is defined");
  domain->dimension = utils::inumeric(FLERR,arg[0],false,lmp);
  if (domain->dimension != 2 && domain->dimension != 3)
    error->all(FLERR, "Invalid dimension argument: {}", arg[0]);

  // must reset default extra_dof of all computes
  // since some were created before dimension command is encountered

  for (auto &c : modify->get_compute_list()) c->reset_extra_dof();
}

/* ---------------------------------------------------------------------- */

void Input::dump()
{
  output->add_dump(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::dump_modify()
{
  output->modify_dump(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::fix()
{
  modify->add_fix(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::fix_modify()
{
  modify->modify_fix(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::group_command()
{
  group->assign(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::improper_coeff()
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Improper_coeff command before simulation box is defined");
  if (force->improper == nullptr)
    error->all(FLERR,"Improper_coeff command before improper_style is defined");
  if (atom->avec->impropers_allow == 0)
    error->all(FLERR,"Improper_coeff command when no impropers allowed");
  char *newarg = utils::expand_type(FLERR, arg[0], Atom::IMPROPER, lmp);
  if (newarg) arg[0] = newarg;
  force->improper->coeff(narg,arg);
  delete[] newarg;
}

/* ---------------------------------------------------------------------- */

void Input::improper_style()
{
  if (narg < 1) error->all(FLERR,"Illegal improper_style command");
  if (atom->avec->impropers_allow == 0)
    error->all(FLERR,"Improper_style command when no impropers allowed");
  force->create_improper(arg[0],1);
  if (force->improper) force->improper->settings(narg-1,&arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::kspace_modify()
{
  if (force->kspace == nullptr)
    error->all(FLERR,"KSpace style has not yet been set");
  force->kspace->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::kspace_style()
{
  force->create_kspace(arg[0],1);
  if (force->kspace) force->kspace->settings(narg-1,&arg[1]);
}
/* ---------------------------------------------------------------------- */

void Input::labelmap()
{
  if (domain->box_exist == 0) error->all(FLERR,"Labelmap command before simulation box is defined");
  if (!atom->labelmapflag) atom->add_label_map();
  atom->lmap->modify_lmap(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::lattice()
{
  domain->set_lattice(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::mass()
{
  if (narg != 2) error->all(FLERR,"Illegal mass command: expected 2 arguments but found {}", narg);
  if (domain->box_exist == 0)
    error->all(FLERR,"Mass command before simulation box is defined");
  atom->set_mass(FLERR,narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::min_modify()
{
  update->minimize->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::min_style()
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Min_style command before simulation box is defined");
  update->create_minimize(narg,arg,1);
}

/* ---------------------------------------------------------------------- */

void Input::molecule()
{
  atom->add_molecule(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::neigh_modify()
{
  neighbor->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::neighbor_command()
{
  neighbor->set(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::newton()
{
  int newton_pair=1,newton_bond=1;

  if (narg == 1) {
    newton_pair = newton_bond = utils::logical(FLERR,arg[0],false,lmp);
  } else if (narg == 2) {
    newton_pair = utils::logical(FLERR,arg[0],false,lmp);
    newton_bond = utils::logical(FLERR,arg[1],false,lmp);
  } else error->all(FLERR,"Illegal newton command");

  force->newton_pair = newton_pair;

  if (domain->box_exist && (newton_bond != force->newton_bond))
    error->all(FLERR,"Newton bond change after simulation box is defined");
  force->newton_bond = newton_bond;

  if (newton_pair || newton_bond) force->newton = 1;
  else force->newton = 0;
}

/* ---------------------------------------------------------------------- */

void Input::package()
{
  if (domain->box_exist)
    error->all(FLERR,"Package command after simulation box is defined");
  if (narg < 1) error->all(FLERR,"Illegal package command");

  // same checks for packages existing as in LAMMPS::post_create()
  // since can be invoked here by package command in input script

  if (strcmp(arg[0],"gpu") == 0) {
    if (!modify->check_package("GPU"))
      error->all(FLERR,"Package gpu command without GPU package installed");

    std::string fixcmd = "package_gpu all GPU";
    for (int i = 1; i < narg; i++) fixcmd += std::string(" ") + arg[i];
    modify->add_fix(fixcmd);

  } else if (strcmp(arg[0],"kokkos") == 0) {
    if (lmp->kokkos == nullptr || lmp->kokkos->kokkos_exists == 0)
      error->all(FLERR, "Package kokkos command without KOKKOS package enabled");
    lmp->kokkos->accelerator(narg-1,&arg[1]);

  } else if (strcmp(arg[0],"omp") == 0) {
    if (!modify->check_package("OMP"))
      error->all(FLERR, "Package omp command without OPENMP package installed");

    std::string fixcmd = "package_omp all OMP";
    for (int i = 1; i < narg; i++) fixcmd += std::string(" ") + arg[i];
    modify->add_fix(fixcmd);

 } else if (strcmp(arg[0],"intel") == 0) {
    if (!modify->check_package("INTEL"))
      error->all(FLERR, "Package intel command without INTEL package installed");

    std::string fixcmd = "package_intel all INTEL";
    for (int i = 1; i < narg; i++) fixcmd += std::string(" ") + arg[i];
    modify->add_fix(fixcmd);

  } else error->all(FLERR,"Unknown package keyword: {}", arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::pair_coeff()
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Pair_coeff command before simulation box is defined");
  if (force->pair == nullptr) error->all(FLERR,"Pair_coeff command without a pair style");
  if (narg < 2) utils::missing_cmd_args(FLERR,"pair_coeff", error);
  if (force->pair->one_coeff && ((strcmp(arg[0],"*") != 0) || (strcmp(arg[1],"*") != 0)))
    error->all(FLERR,"Pair_coeff must start with * * for pair style {}", force->pair_style);

  char *newarg0 = utils::expand_type(FLERR, arg[0], Atom::ATOM, lmp);
  if (newarg0) arg[0] = newarg0;
  char *newarg1 = utils::expand_type(FLERR, arg[1], Atom::ATOM, lmp);
  if (newarg1) arg[1] = newarg1;

  // if arg[1] < arg[0], and neither contain a wildcard, reorder

  int itype,jtype;
  if (utils::strmatch(arg[0],"^\\d+$") && utils::strmatch(arg[1],"^\\d+$")) {
    itype = utils::inumeric(FLERR,arg[0],false,lmp);
    jtype = utils::inumeric(FLERR,arg[1],false,lmp);
    if (jtype < itype) {
      char *str = arg[0];
      arg[0] = arg[1];
      arg[1] = str;
    }
  }

  force->pair->coeff(narg,arg);
  delete[] newarg0;
  delete[] newarg1;
}

/* ---------------------------------------------------------------------- */

void Input::pair_modify()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Pair_modify command before pair_style is defined");
  force->pair->modify_params(narg,arg);
}

/* ----------------------------------------------------------------------
   if old pair style exists and new style is same, just change settings
   else create new pair class
------------------------------------------------------------------------- */

void Input::pair_style()
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "pair_style", error);
  if (force->pair) {
    std::string style = arg[0];
    int match = 0;
    if (style == force->pair_style) match = 1;
    if (!match && lmp->suffix_enable) {
      if (lmp->suffix)
        if (style + "/" + lmp->suffix == force->pair_style) match = 1;

      if (lmp->suffix2)
        if (style + "/" + lmp->suffix2 == force->pair_style) match = 1;
    }
    if (match) {
      force->pair->settings(narg-1,&arg[1]);
      return;
    }
  }

  force->create_pair(arg[0],1);
  if (force->pair) force->pair->settings(narg-1,&arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::pair_write()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Pair_write command before pair_style is defined");
  force->pair->write_file(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::processors()
{
  if (domain->box_exist)
    error->all(FLERR,"Processors command after simulation box is defined");
  comm->set_processors(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::region()
{
  domain->add_region(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::reset_timestep()
{
  update->reset_timestep(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::restart()
{
  output->create_restart(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::run_style()
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Run_style command before simulation box is defined");
  update->create_integrate(narg,arg,1);
}

/* ---------------------------------------------------------------------- */

void Input::special_bonds()
{
  // store 1-3,1-4 and dihedral/extra flag values before change
  // change in 1-2 coeffs will not change the special list

  double lj2 = force->special_lj[2];
  double lj3 = force->special_lj[3];
  double coul2 = force->special_coul[2];
  double coul3 = force->special_coul[3];
  int onefive = force->special_onefive;
  int angle = force->special_angle;
  int dihedral = force->special_dihedral;

  force->set_special(narg,arg);

  // if simulation box defined and saved values changed, redo special list

  if (domain->box_exist && atom->molecular == Atom::MOLECULAR) {
    if (lj2 != force->special_lj[2] || lj3 != force->special_lj[3] ||
        coul2 != force->special_coul[2] || coul3 != force->special_coul[3] ||
        onefive != force->special_onefive ||
        angle != force->special_angle ||
        dihedral != force->special_dihedral) {
      Special special(lmp);
      special.build();
    }
  }
}

/* ---------------------------------------------------------------------- */

void Input::suffix()
{
  if (narg < 1) error->all(FLERR,"Illegal suffix command");

  const std::string firstarg = arg[0];

  if ((firstarg == "off") || (firstarg == "no") || (firstarg == "false")) {
    lmp->suffix_enable = 0;
  } else if ((firstarg == "on") || (firstarg == "yes") || (firstarg == "true")) {
    lmp->suffix_enable = 1;
    if (!lmp->suffix) error->all(FLERR,"May only enable suffixes after defining one");
  } else {
    lmp->suffix_enable = 1;

    delete[] lmp->suffix;
    delete[] lmp->suffix2;
    lmp->suffix = lmp->suffix2 = nullptr;

    if (firstarg == "hybrid") {
      if (narg != 3) error->all(FLERR,"Illegal suffix command");
      lmp->suffix = utils::strdup(arg[1]);
      lmp->suffix2 = utils::strdup(arg[2]);
    } else {
      if (narg != 1) error->all(FLERR,"Illegal suffix command");
      lmp->suffix = utils::strdup(arg[0]);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Input::thermo()
{
  output->set_thermo(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::thermo_modify()
{
  output->thermo->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::thermo_style()
{
  int nline = *output->thermo->get_line();
  output->create_thermo(narg,arg);
  output->thermo->set_line(nline);
}

/* ---------------------------------------------------------------------- */

void Input::timer_command()
{
  timer->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::timestep()
{
  if (narg != 1) error->all(FLERR,"Illegal timestep command");

  update->update_time();
  update->dt = utils::numeric(FLERR,arg[0],false,lmp);
  update->dt_default = 0;

  // timestep command can be invoked between runs or by run every
  // calls to other classes that need to know timestep size changed
  // similar logic is in FixDtReset::end_of_step()
  // only need to do this if a run has already occurred

  if (update->first_update == 0) return;

  int respaflag = 0;
  if (utils::strmatch(update->integrate_style, "^respa")) respaflag = 1;
  if (respaflag) update->integrate->reset_dt();

  if (force->pair) force->pair->reset_dt();
  for (auto &ifix : modify->get_fix_list()) ifix->reset_dt();
  output->reset_dt();
}

/* ---------------------------------------------------------------------- */

void Input::uncompute()
{
  if (narg != 1) error->all(FLERR,"Illegal uncompute command");
  modify->delete_compute(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::undump()
{
  if (narg != 1) error->all(FLERR,"Illegal undump command");
  output->delete_dump(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::unfix()
{
  if (narg != 1) error->all(FLERR,"Illegal unfix command");
  modify->delete_fix(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::units()
{
  if (narg != 1) error->all(FLERR,"Illegal units command: expected 1 argument but found {}", narg);
  if (domain->box_exist)
    error->all(FLERR,"Units command after simulation box is defined");
  update->set_units(arg[0]);
}

/* ---------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   function for meta commands
------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

int Input::meta(const std::string &prefix)
{
  auto mycmd = fmt::format("{}_{}", utils::uppercase(prefix), utils::uppercase(arg[0]));
  if (command_map->find(mycmd) != command_map->end()) {
    CommandCreator &command_creator = (*command_map)[mycmd];
    Command *cmd = command_creator(lmp);
    cmd->command(narg-1,arg+1);
    delete cmd;
    return 1;
  } else return 0;
}
