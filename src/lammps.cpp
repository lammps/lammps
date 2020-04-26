/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "lammps.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <map>
#include <string>
#include "style_angle.h"     // IWYU pragma: keep
#include "style_atom.h"      // IWYU pragma: keep
#include "style_bond.h"      // IWYU pragma: keep
#include "style_command.h"   // IWYU pragma: keep
#include "style_compute.h"   // IWYU pragma: keep
#include "style_dihedral.h"  // IWYU pragma: keep
#include "style_dump.h"      // IWYU pragma: keep
#include "style_fix.h"       // IWYU pragma: keep
#include "style_improper.h"  // IWYU pragma: keep
#include "style_integrate.h" // IWYU pragma: keep
#include "style_kspace.h"    // IWYU pragma: keep
#include "style_minimize.h"  // IWYU pragma: keep
#include "style_pair.h"      // IWYU pragma: keep
#include "style_region.h"    // IWYU pragma: keep
#include "universe.h"
#include "input.h"
#include "info.h"
#include "atom.h"
#include "update.h"
#include "neighbor.h"
#include "comm.h"
#include "comm_brick.h"
#include "domain.h"
#include "force.h"
#include "modify.h"
#include "group.h"
#include "output.h"
#include "citeme.h"
#include "accelerator_kokkos.h"
#include "accelerator_omp.h"    // IWYU pragma: keep
#include "timer.h"
#include "lmppython.h"
#include "version.h"
#include "memory.h"
#include "error.h"

#include "lmpinstalledpkgs.h"
#include "lmpgitversion.h"

static void print_style(FILE *fp, const char *str, int &pos);

struct LAMMPS_NS::package_styles_lists {
  std::map<std::string,std::string> angle_styles;
  std::map<std::string,std::string> atom_styles;
  std::map<std::string,std::string> body_styles;
  std::map<std::string,std::string> bond_styles;
  std::map<std::string,std::string> command_styles;
  std::map<std::string,std::string> compute_styles;
  std::map<std::string,std::string> dihedral_styles;
  std::map<std::string,std::string> dump_styles;
  std::map<std::string,std::string> fix_styles;
  std::map<std::string,std::string> improper_styles;
  std::map<std::string,std::string> integrate_styles;
  std::map<std::string,std::string> kspace_styles;
  std::map<std::string,std::string> minimize_styles;
  std::map<std::string,std::string> pair_styles;
  std::map<std::string,std::string> reader_styles;
  std::map<std::string,std::string> region_styles;
};

using namespace LAMMPS_NS;

/** \class LAMMPS_NS::LAMMPS
 * \brief Composite class of a LAMMPS simulation instance
 *
 * The LAMMPS class contains pointers of all constituent
 * class instances and global variables that are used by
 * a LAMMPS simulation instance and thus represent the
 * state of the simulation.
 */

/** \brief LAMMPS constructor
 *
 * The LAMMPS constructor starts up a simulation by allocating all
 * fundamental classes in the necessary order, parses input switches
 * and their arguments, initializes communicators, screen and logfile
 * output FILE pointers.
 *
 * \param narg number of arguments
 * \param arg list of arguments
 * \param communicator MPI communicator used by this LAMMPS instance
 */
LAMMPS::LAMMPS(int narg, char **arg, MPI_Comm communicator) :
  memory(NULL), error(NULL), universe(NULL), input(NULL), atom(NULL),
  update(NULL), neighbor(NULL), comm(NULL), domain(NULL), force(NULL),
  modify(NULL), group(NULL), output(NULL), timer(NULL), kokkos(NULL),
  atomKK(NULL), memoryKK(NULL), python(NULL), citeme(NULL)
{
  memory = new Memory(this);
  error = new Error(this);
  universe = new Universe(this,communicator);

  clientserver = 0;
  cslib = NULL;
  cscomm = 0;

  screen = NULL;
  logfile = NULL;
  infile = NULL;

  initclock = MPI_Wtime();

  init_pkg_lists();

  // check if -mpi is first arg
  // if so, then 2 apps were launched with one mpirun command
  //   this means passed communicator (e.g. MPI_COMM_WORLD) is bigger than LAMMPS
  //     e.g. for client/server coupling with another code
  //     in the future LAMMPS might leverage this in other ways
  //   universe communicator needs to shrink to be just LAMMPS
  // syntax: -mpi color
  //   color = integer for this app, different than other app(s)
  // do the following:
  //   perform an MPI_Comm_split() to create a new LAMMPS-only subcomm
  //   NOTE: this assumes other app(s) does same thing, else will hang!
  //   re-create universe with subcomm
  //   store full multi-app comm in cscomm
  //   cscomm is used by CSLIB package to exchange messages w/ other app

  int iarg = 1;
  if (narg-iarg >= 2 && (strcmp(arg[iarg],"-mpi") == 0 ||
                         strcmp(arg[iarg],"-m") == 0)) {
    int me,nprocs;
    MPI_Comm_rank(communicator,&me);
    MPI_Comm_size(communicator,&nprocs);
    int color = atoi(arg[iarg+1]);
    MPI_Comm subcomm;
    MPI_Comm_split(communicator,color,me,&subcomm);
    cscomm = communicator;
    communicator = subcomm;
    delete universe;
    universe = new Universe(this,communicator);
  }

  // parse input switches

  int inflag = 0;
  int screenflag = 0;
  int logflag = 0;
  int partscreenflag = 0;
  int partlogflag = 0;
  int kokkosflag = 0;
  int restart2data = 0;
  int restart2dump = 0;
  int restartremap = 0;
  int citeflag = 1;
  int helpflag = 0;

  suffix = suffix2 = NULL;
  suffix_enable = 0;
  if (arg) exename = arg[0];
  else exename = NULL;
  packargs = NULL;
  num_package = 0;
  char *restartfile = NULL;
  int wfirst,wlast;
  int kkfirst,kklast;

  int npack = 0;
  int *pfirst = NULL;
  int *plast = NULL;

  iarg = 1;
  while (iarg < narg) {

    if (strcmp(arg[iarg],"-echo") == 0 ||
               strcmp(arg[iarg],"-e") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 2;

    } else if (strcmp(arg[iarg],"-help") == 0 ||
               strcmp(arg[iarg],"-h") == 0) {
      if (iarg+1 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      helpflag = 1;
      citeflag = 0;
      iarg += 1;

    } else if (strcmp(arg[iarg],"-in") == 0 ||
               strcmp(arg[iarg],"-i") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      inflag = iarg + 1;
      iarg += 2;

    } else if (strcmp(arg[iarg],"-kokkos") == 0 ||
               strcmp(arg[iarg],"-k") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      if (strcmp(arg[iarg+1],"on") == 0) kokkosflag = 1;
      else if (strcmp(arg[iarg+1],"off") == 0) kokkosflag = 0;
      else error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 2;
      // delimit any extra args for the Kokkos instantiation
      kkfirst = iarg;
      while (iarg < narg && arg[iarg][0] != '-') iarg++;
      kklast = iarg;

    } else if (strcmp(arg[iarg],"-log") == 0 ||
               strcmp(arg[iarg],"-l") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      logflag = iarg + 1;
      iarg += 2;

    } else if (strcmp(arg[iarg],"-mpi") == 0 ||
               strcmp(arg[iarg],"-m") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      if (iarg != 1) error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 2;

    } else if (strcmp(arg[iarg],"-nocite") == 0 ||
               strcmp(arg[iarg],"-nc") == 0) {
      citeflag = 0;
      iarg++;

    } else if (strcmp(arg[iarg],"-package") == 0 ||
               strcmp(arg[iarg],"-pk") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      memory->grow(pfirst,npack+1,"lammps:pfirst");
      memory->grow(plast,npack+1,"lammps:plast");
      // delimit args for package command invocation
      // any package arg with leading "-" will be followed by numeric digit
      iarg++;
      pfirst[npack] = iarg;
      while (iarg < narg) {
        if (arg[iarg][0] != '-') iarg++;
        else if (isdigit(arg[iarg][1])) iarg++;
        else break;
      }
      plast[npack++] = iarg;

    } else if (strcmp(arg[iarg],"-partition") == 0 ||
        strcmp(arg[iarg],"-p") == 0) {
      universe->existflag = 1;
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      iarg++;
      while (iarg < narg && arg[iarg][0] != '-') {
        universe->add_world(arg[iarg]);
        iarg++;
      }

    } else if (strcmp(arg[iarg],"-plog") == 0 ||
               strcmp(arg[iarg],"-pl") == 0) {
      if (iarg+2 > narg)
       error->universe_all(FLERR,"Invalid command-line argument");
      partlogflag = iarg + 1;
      iarg += 2;

    } else if (strcmp(arg[iarg],"-pscreen") == 0 ||
               strcmp(arg[iarg],"-ps") == 0) {
      if (iarg+2 > narg)
       error->universe_all(FLERR,"Invalid command-line argument");
      partscreenflag = iarg + 1;
      iarg += 2;

    } else if (strcmp(arg[iarg],"-reorder") == 0 ||
               strcmp(arg[iarg],"-ro") == 0) {
      if (iarg+3 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      if (universe->existflag)
        error->universe_all(FLERR,"Cannot use -reorder after -partition");
      universe->reorder(arg[iarg+1],arg[iarg+2]);
      iarg += 3;

    } else if (strcmp(arg[iarg],"-restart2data") == 0 ||
               strcmp(arg[iarg],"-r2data") == 0) {
      if (iarg+3 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      if (restart2dump)
        error->universe_all(FLERR,
                            "Cannot use both -restart2data and -restart2dump");
      restart2data = 1;
      restartfile = arg[iarg+1];
      // check for restart remap flag
      if (strcmp(arg[iarg+2],"remap") == 0) {
        if (iarg+4 > narg)
          error->universe_all(FLERR,"Invalid command-line argument");
        restartremap = 1;
        iarg++;
      }
      iarg += 2;
      // delimit args for the write_data command
      wfirst = iarg;
      while (iarg < narg && arg[iarg][0] != '-') iarg++;
      wlast = iarg;

    } else if (strcmp(arg[iarg],"-restart2dump") == 0 ||
               strcmp(arg[iarg],"-r2dump") == 0) {
      if (iarg+3 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      if (restart2data)
        error->universe_all(FLERR,
                            "Cannot use both -restart2data and -restart2dump");
      restart2dump = 1;
      restartfile = arg[iarg+1];
      // check for restart remap flag
      if (strcmp(arg[iarg+2],"remap") == 0) {
        if (iarg+4 > narg)
          error->universe_all(FLERR,"Invalid command-line argument");
        restartremap = 1;
        iarg++;
      }
      iarg += 2;
      // delimit args for the write_dump command
      wfirst = iarg;
      while (iarg < narg && arg[iarg][0] != '-') iarg++;
      wlast = iarg;

    } else if (strcmp(arg[iarg],"-screen") == 0 ||
               strcmp(arg[iarg],"-sc") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      screenflag = iarg + 1;
      iarg += 2;

    } else if (strcmp(arg[iarg],"-suffix") == 0 ||
               strcmp(arg[iarg],"-sf") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      delete [] suffix;
      delete [] suffix2;
      suffix2 = NULL;
      suffix_enable = 1;
      // hybrid option to set fall-back for suffix2
      if (strcmp(arg[iarg+1],"hybrid") == 0) {
        if (iarg+4 > narg)
          error->universe_all(FLERR,"Invalid command-line argument");
        int n = strlen(arg[iarg+2]) + 1;
        suffix = new char[n];
        strcpy(suffix,arg[iarg+2]);
        n = strlen(arg[iarg+3]) + 1;
        suffix2 = new char[n];
        strcpy(suffix2,arg[iarg+3]);
        iarg += 4;
      } else {
        int n = strlen(arg[iarg+1]) + 1;
        suffix = new char[n];
        strcpy(suffix,arg[iarg+1]);
        iarg += 2;
      }

    } else if (strcmp(arg[iarg],"-var") == 0 ||
               strcmp(arg[iarg],"-v") == 0) {
      if (iarg+3 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 3;
      while (iarg < narg && arg[iarg][0] != '-') iarg++;

    } else error->universe_all(FLERR,"Invalid command-line argument");
  }

  // if no partition command-line switch, universe is one world with all procs

  if (universe->existflag == 0) universe->add_world(NULL);

  // sum of procs in all worlds must equal total # of procs

  if (!universe->consistent())
    error->universe_all(FLERR,"Processor partitions do not match "
                        "number of allocated processors");

  // universe cannot use stdin for input file

  if (universe->existflag && inflag == 0)
    error->universe_all(FLERR,"Must use -in switch with multiple partitions");

  // if no partition command-line switch, cannot use -pscreen option

  if (universe->existflag == 0 && partscreenflag)
    error->universe_all(FLERR,"Can only use -pscreen with multiple partitions");

  // if no partition command-line switch, cannot use -plog option

  if (universe->existflag == 0 && partlogflag)
    error->universe_all(FLERR,"Can only use -plog with multiple partitions");

  // set universe screen and logfile

  if (universe->me == 0) {
    if (screenflag == 0)
      universe->uscreen = stdout;
    else if (strcmp(arg[screenflag],"none") == 0)
      universe->uscreen = NULL;
    else {
      universe->uscreen = fopen(arg[screenflag],"w");
      if (universe->uscreen == NULL)
        error->universe_one(FLERR,"Cannot open universe screen file");
    }
    if (logflag == 0) {
      if (helpflag == 0) {
        universe->ulogfile = fopen("log.lammps","w");
        if (universe->ulogfile == NULL)
          error->universe_warn(FLERR,"Cannot open log.lammps for writing");
      }
    } else if (strcmp(arg[logflag],"none") == 0)
      universe->ulogfile = NULL;
    else {
      universe->ulogfile = fopen(arg[logflag],"w");
      if (universe->ulogfile == NULL)
        error->universe_one(FLERR,"Cannot open universe log file");
    }
  }

  if (universe->me > 0) {
    if (screenflag == 0) universe->uscreen = stdout;
    else universe->uscreen = NULL;
    universe->ulogfile = NULL;
  }

  // make universe and single world the same, since no partition switch
  // world inherits settings from universe
  // set world screen, logfile, communicator, infile
  // open input script if from file

  if (universe->existflag == 0) {
    screen = universe->uscreen;
    logfile = universe->ulogfile;
    world = universe->uworld;

    if (universe->me == 0) {
      if (inflag == 0) infile = stdin;
      else infile = fopen(arg[inflag],"r");
      if (infile == NULL) {
        char str[128];
        snprintf(str,128,"Cannot open input script %s",arg[inflag]);
        error->one(FLERR,str);
      }
    }

    if ((universe->me == 0) && !helpflag) {
      if (screen) fprintf(screen,"LAMMPS (%s)\n",universe->version);
      if (logfile) fprintf(logfile,"LAMMPS (%s)\n",universe->version);
    }

  // universe is one or more worlds, as setup by partition switch
  // split universe communicator into separate world communicators
  // set world screen, logfile, communicator, infile
  // open input script

  } else {
    int me;
    MPI_Comm_split(universe->uworld,universe->iworld,0,&world);
    MPI_Comm_rank(world,&me);

    if (me == 0)
      if (partscreenflag == 0)
       if (screenflag == 0) {
         char str[32];
         sprintf(str,"screen.%d",universe->iworld);
         screen = fopen(str,"w");
         if (screen == NULL) error->one(FLERR,"Cannot open screen file");
       } else if (strcmp(arg[screenflag],"none") == 0)
         screen = NULL;
       else {
         char str[128];
         snprintf(str,128,"%s.%d",arg[screenflag],universe->iworld);
         screen = fopen(str,"w");
         if (screen == NULL) error->one(FLERR,"Cannot open screen file");
       }
      else if (strcmp(arg[partscreenflag],"none") == 0)
        screen = NULL;
      else {
        char str[128];
        snprintf(str,128,"%s.%d",arg[partscreenflag],universe->iworld);
        screen = fopen(str,"w");
        if (screen == NULL) error->one(FLERR,"Cannot open screen file");
      } else screen = NULL;

    if (me == 0)
      if (partlogflag == 0)
       if (logflag == 0) {
         char str[32];
         sprintf(str,"log.lammps.%d",universe->iworld);
         logfile = fopen(str,"w");
         if (logfile == NULL) error->one(FLERR,"Cannot open logfile");
       } else if (strcmp(arg[logflag],"none") == 0)
         logfile = NULL;
       else {
         char str[128];
         snprintf(str,128,"%s.%d",arg[logflag],universe->iworld);
         logfile = fopen(str,"w");
         if (logfile == NULL) error->one(FLERR,"Cannot open logfile");
       }
      else if (strcmp(arg[partlogflag],"none") == 0)
        logfile = NULL;
      else {
        char str[128];
        snprintf(str,128,"%s.%d",arg[partlogflag],universe->iworld);
        logfile = fopen(str,"w");
        if (logfile == NULL) error->one(FLERR,"Cannot open logfile");
      } else logfile = NULL;

    if (me == 0) {
      infile = fopen(arg[inflag],"r");
      if (infile == NULL) {
        char str[128];
        snprintf(str,128,"Cannot open input script %s",arg[inflag]);
        error->one(FLERR,str);
      }
    } else infile = NULL;

    // screen and logfile messages for universe and world

    if ((universe->me == 0) && (!helpflag)) {
      if (universe->uscreen) {
        fprintf(universe->uscreen,"LAMMPS (%s)\n",universe->version);
        fprintf(universe->uscreen,"Running on %d partitions of processors\n",
                universe->nworlds);
      }
      if (universe->ulogfile) {
        fprintf(universe->ulogfile,"LAMMPS (%s)\n",universe->version);
        fprintf(universe->ulogfile,"Running on %d partitions of processors\n",
                universe->nworlds);
      }
    }

    if ((me == 0) && (!helpflag)) {
      if (screen) {
        fprintf(screen,"LAMMPS (%s)\n",universe->version);
        fprintf(screen,"Processor partition = %d\n",universe->iworld);
      }
      if (logfile) {
        fprintf(logfile,"LAMMPS (%s)\n",universe->version);
        fprintf(logfile,"Processor partition = %d\n",universe->iworld);
      }
    }
  }

  // check consistency of datatype settings in lmptype.h

  if (sizeof(smallint) != sizeof(int))
    error->all(FLERR,"Smallint setting in lmptype.h is invalid");
  if (sizeof(imageint) < sizeof(smallint))
    error->all(FLERR,"Imageint setting in lmptype.h is invalid");
  if (sizeof(tagint) < sizeof(smallint))
    error->all(FLERR,"Tagint setting in lmptype.h is invalid");
  if (sizeof(bigint) < sizeof(imageint) || sizeof(bigint) < sizeof(tagint))
    error->all(FLERR,"Bigint setting in lmptype.h is invalid");

  int mpisize;
  MPI_Type_size(MPI_LMP_TAGINT,&mpisize);
  if (mpisize != sizeof(tagint))
      error->all(FLERR,"MPI_LMP_TAGINT and tagint in "
                 "lmptype.h are not compatible");
  MPI_Type_size(MPI_LMP_BIGINT,&mpisize);
  if (mpisize != sizeof(bigint))
      error->all(FLERR,"MPI_LMP_BIGINT and bigint in "
                 "lmptype.h are not compatible");

#ifdef LAMMPS_SMALLBIG
  if (sizeof(smallint) != 4 || sizeof(imageint) != 4 ||
      sizeof(tagint) != 4 || sizeof(bigint) != 8)
    error->all(FLERR,"Small to big integers are not sized correctly");
#endif
#ifdef LAMMPS_BIGBIG
  if (sizeof(smallint) != 4 || sizeof(imageint) != 8 ||
      sizeof(tagint) != 8 || sizeof(bigint) != 8)
    error->all(FLERR,"Small to big integers are not sized correctly");
#endif
#ifdef LAMMPS_SMALLSMALL
  if (sizeof(smallint) != 4 || sizeof(imageint) != 4 ||
      sizeof(tagint) != 4 || sizeof(bigint) != 4)
    error->all(FLERR,"Small to big integers are not sized correctly");
#endif

  // create Kokkos class if KOKKOS installed, unless explicitly switched off
  // instantiation creates dummy Kokkos class if KOKKOS is not installed
  // add args between kkfirst and kklast to Kokkos instantiation

  kokkos = NULL;
  if (kokkosflag == 1) {
    kokkos = new KokkosLMP(this,kklast-kkfirst,&arg[kkfirst]);
    if (!kokkos->kokkos_exists)
      error->all(FLERR,"Cannot use -kokkos on without KOKKOS installed");
  }

  // allocate CiteMe class if enabled

  if (citeflag) citeme = new CiteMe(this);
  else citeme = NULL;

  // allocate input class now that MPI is fully setup

  input = new Input(this,narg,arg);

  // copy package cmdline arguments

  if (npack > 0) {
    num_package = npack;
    packargs = new char**[npack];
    for (int i=0; i < npack; ++i) {
      int n = plast[i] - pfirst[i];
      packargs[i] = new char*[n+1];
      for (int j=0; j < n; ++j)
        packargs[i][j] = strdup(arg[pfirst[i]+j]);
      packargs[i][n] = NULL;
    }
    memory->destroy(pfirst);
    memory->destroy(plast);
  }

  // if helpflag set, print help and quit with "success" status
  // otherwise allocate top level classes.

  if (helpflag) {
    if (universe->me == 0 && screen) help();
    error->done(0);
  } else {
    create();
    post_create();
  }

  // if either restart conversion option was used, invoke 2 commands and quit
  // add args between wfirst and wlast to write_data or write_data command
  // add "noinit" to write_data to prevent a system init
  // write_dump will just give a warning message about no init

  if (restart2data || restart2dump) {
    char cmd[256];
    snprintf(cmd,248,"read_restart %s\n",restartfile);
    if (restartremap) strcat(cmd," remap\n");
    input->one(cmd);
    if (restart2data) strcpy(cmd,"write_data");
    else strcpy(cmd,"write_dump");
    for (iarg = wfirst; iarg < wlast; iarg++)
      snprintf(&cmd[strlen(cmd)],246-strlen(cmd)," %s",arg[iarg]);
    if (restart2data) strcat(cmd," noinit\n");
    else strcat(cmd,"\n");
    input->one(cmd);
    error->done(0);
  }
}

/* ----------------------------------------------------------------------
   shutdown LAMMPS
   delete top-level classes
   close screen and log files in world and universe
   output files were already closed in destroy()
   delete fundamental classes
------------------------------------------------------------------------- */

LAMMPS::~LAMMPS()
{
  const int me = comm->me;

  destroy();
  delete citeme;

  if (num_package) {
    for (int i = 0; i < num_package; i++) {
      for (char **ptr = packargs[i]; *ptr != NULL; ++ptr)
        free(*ptr);
      delete[] packargs[i];
    }
    delete[] packargs;
  }
  num_package = 0;
  packargs = NULL;

  double totalclock = MPI_Wtime() - initclock;
  if ((me == 0) && (screen || logfile)) {
    char outtime[128];
    int seconds = fmod(totalclock,60.0);
    totalclock  = (totalclock - seconds) / 60.0;
    int minutes = fmod(totalclock,60.0);
    int hours = (totalclock - minutes) / 60.0;
    sprintf(outtime,"Total wall time: "
            "%d:%02d:%02d\n", hours, minutes, seconds);
    if (screen) fputs(outtime,screen);
    if (logfile) fputs(outtime,logfile);
  }

  if (universe->nworlds == 1) {
    if (screen && screen != stdout) fclose(screen);
    if (logfile) fclose(logfile);
    logfile = NULL;
    if (screen != stdout) screen = NULL;
  } else {
    if (screen && screen != stdout) fclose(screen);
    if (logfile) fclose(logfile);
    if (universe->ulogfile) fclose(universe->ulogfile);
    logfile = NULL;
    if (screen != stdout) screen = NULL;
  }

  if (infile && infile != stdin) fclose(infile);

  if (world != universe->uworld) MPI_Comm_free(&world);

  delete python;
  delete kokkos;
  delete [] suffix;
  delete [] suffix2;

  // free the MPI comm created by -mpi command-line arg processed in constructor
  // it was passed to universe as if original universe world
  // may have been split later by partitions, universe will free the splits
  // free a copy of uorig here, so check in universe destructor will still work

  MPI_Comm copy = universe->uorig;
  if (cscomm) MPI_Comm_free(&copy);

  delete input;
  delete universe;
  delete error;
  delete memory;

  delete pkg_lists;
}

/* ----------------------------------------------------------------------
   allocate single instance of top-level classes
   fundamental classes are allocated in constructor
   some classes have package variants
------------------------------------------------------------------------- */

void LAMMPS::create()
{
  force = NULL;         // Domain->Lattice checks if Force exists

  // Comm class must be created before Atom class
  // so that nthreads is defined when create_avec invokes grow()

  if (kokkos) comm = new CommKokkos(this);
  else comm = new CommBrick(this);

  if (kokkos) neighbor = new NeighborKokkos(this);
  else neighbor = new Neighbor(this);

  if (kokkos) domain = new DomainKokkos(this);
#ifdef LMP_USER_OMP
  else domain = new DomainOMP(this);
#else
  else domain = new Domain(this);
#endif

  if (kokkos) atom = new AtomKokkos(this);
  else atom = new Atom(this);

  if (kokkos)
    atom->create_avec("atomic/kk",0,NULL,1);
  else
    atom->create_avec("atomic",0,NULL,1);

  group = new Group(this);
  force = new Force(this);    // must be after group, to create temperature

  if (kokkos) modify = new ModifyKokkos(this);
  else modify = new Modify(this);

  output = new Output(this);  // must be after group, so "all" exists
                              // must be after modify so can create Computes
  update = new Update(this);  // must be after output, force, neighbor
  timer = new Timer(this);

  python = new Python(this);
}

/* ----------------------------------------------------------------------
   check suffix consistency with installed packages
   invoke package-specific default package commands
     only invoke if suffix is set and enabled
     also check if suffix2 is set
   called from LAMMPS constructor and after clear() command
     so that package-specific core classes have been instantiated
------------------------------------------------------------------------- */

void LAMMPS::post_create()
{
  // default package command triggered by "-k on"

  if (kokkos && kokkos->kokkos_exists) input->one("package kokkos");

  // suffix will always be set if suffix_enable = 1
  // check that KOKKOS package classes were instantiated
  // check that GPU, INTEL, USER-OMP fixes were compiled with LAMMPS

  if (suffix_enable) {

    if (strcmp(suffix,"gpu") == 0 && !modify->check_package("GPU"))
      error->all(FLERR,"Using suffix gpu without GPU package installed");
    if (strcmp(suffix,"intel") == 0 && !modify->check_package("INTEL"))
      error->all(FLERR,"Using suffix intel without USER-INTEL package installed");
    if (strcmp(suffix,"kk") == 0 &&
        (kokkos == NULL || kokkos->kokkos_exists == 0))
      error->all(FLERR,"Using suffix kk without KOKKOS package enabled");
    if (strcmp(suffix,"omp") == 0 && !modify->check_package("OMP"))
      error->all(FLERR,"Using suffix omp without USER-OMP package installed");

    if (strcmp(suffix,"gpu") == 0) input->one("package gpu 1");
    if (strcmp(suffix,"intel") == 0) input->one("package intel 1");
    if (strcmp(suffix,"omp") == 0) input->one("package omp 0");

    if (suffix2) {
      if (strcmp(suffix2,"gpu") == 0) input->one("package gpu 1");
      if (strcmp(suffix2,"intel") == 0) input->one("package intel 1");
      if (strcmp(suffix2,"omp") == 0) input->one("package omp 0");
    }
  }

  // invoke any command-line package commands

  if (num_package) {
    char str[256];
    for (int i = 0; i < num_package; i++) {
      strcpy(str,"package");
      for (char **ptr = packargs[i]; *ptr != NULL; ++ptr) {
        if (strlen(str) + strlen(*ptr) + 2 > 256)
          error->all(FLERR,"Too many -pk arguments in command line");
        strcat(str," ");
        strcat(str,*ptr);
      }
      input->one(str);
    }
  }
}

/* ----------------------------------------------------------------------
   initialize top-level classes
   do not initialize Timer class, other classes like Run() do that explicitly
------------------------------------------------------------------------- */

void LAMMPS::init()
{
  update->init();
  force->init();         // pair must come after update due to minimizer
  domain->init();
  atom->init();          // atom must come after force and domain
                         //   atom deletes extra array
                         //   used by fix shear_history::unpack_restart()
                         //     when force->pair->gran_history creates fix
                         //   atom_vec init uses deform_vremap
  modify->init();        // modify must come after update, force, atom, domain
  neighbor->init();      // neighbor must come after force, modify
  comm->init();          // comm must come after force, modify, neighbor, atom
  output->init();        // output must come after domain, force, modify
}

/* ----------------------------------------------------------------------
   delete single instance of top-level classes
   fundamental classes are deleted in destructor
------------------------------------------------------------------------- */

void LAMMPS::destroy()
{
  delete update;
  update = NULL;

  delete neighbor;
  neighbor = NULL;

  delete comm;
  comm = NULL;

  delete force;
  force = NULL;

  delete group;
  group = NULL;

  delete output;
  output = NULL;

  delete modify;          // modify must come after output, force, update
                          //   since they delete fixes
  modify = NULL;

  delete domain;          // domain must come after modify
                          //   since fix destructors access domain
  domain = NULL;

  delete atom;            // atom must come after modify, neighbor
                          //   since fixes delete callbacks in atom
  atom = NULL;

  delete timer;
  timer = NULL;

  delete python;
  python = NULL;
}

/* ----------------------------------------------------------------------
   initialize lists of styles in packages
------------------------------------------------------------------------- */

void _noopt LAMMPS::init_pkg_lists()
{
  pkg_lists = new package_styles_lists;
#define PACKAGE "UNKNOWN"
#define ANGLE_CLASS
#define AngleStyle(key,Class)                   \
  pkg_lists->angle_styles[#key] = PACKAGE;
#include "packages_angle.h"
#undef AngleStyle
#undef ANGLE_CLASS
#define ATOM_CLASS
#define AtomStyle(key,Class)                    \
  pkg_lists->atom_styles[#key] = PACKAGE;
#include "packages_atom.h"
#undef AtomStyle
#undef ATOM_CLASS
#define BODY_CLASS
#define BodyStyle(key,Class)                    \
  pkg_lists->body_styles[#key] = PACKAGE;
#include "packages_body.h"
#undef BodyStyle
#undef BODY_CLASS
#define BOND_CLASS
#define BondStyle(key,Class)                    \
  pkg_lists->bond_styles[#key] = PACKAGE;
#include "packages_bond.h"
#undef BondStyle
#undef BOND_CLASS
#define COMMAND_CLASS
#define CommandStyle(key,Class)                 \
  pkg_lists->command_styles[#key] = PACKAGE;
#include "packages_command.h"
#undef CommandStyle
#undef COMMAND_CLASS
#define COMPUTE_CLASS
#define ComputeStyle(key,Class)                 \
  pkg_lists->compute_styles[#key] = PACKAGE;
#include "packages_compute.h"
#undef ComputeStyle
#undef COMPUTE_CLASS
#define DIHEDRAL_CLASS
#define DihedralStyle(key,Class)                \
  pkg_lists->dihedral_styles[#key] = PACKAGE;
#include "packages_dihedral.h"
#undef DihedralStyle
#undef DIHEDRAL_CLASS
#define DUMP_CLASS
#define DumpStyle(key,Class)                    \
  pkg_lists->dump_styles[#key] = PACKAGE;
#include "packages_dump.h"
#undef DumpStyle
#undef DUMP_CLASS
#define FIX_CLASS
#define FixStyle(key,Class)                     \
  pkg_lists->fix_styles[#key] = PACKAGE;
#include "packages_fix.h"
#undef FixStyle
#undef FIX_CLASS
#define IMPROPER_CLASS
#define ImproperStyle(key,Class)                \
  pkg_lists->improper_styles[#key] = PACKAGE;
#include "packages_improper.h"
#undef ImproperStyle
#undef IMPROPER_CLASS
#define INTEGRATE_CLASS
#define IntegrateStyle(key,Class)               \
  pkg_lists->integrate_styles[#key] = PACKAGE;
#include "packages_integrate.h"
#undef IntegrateStyle
#undef INTEGRATE_CLASS
#define KSPACE_CLASS
#define KSpaceStyle(key,Class)                  \
  pkg_lists->kspace_styles[#key] = PACKAGE;
#include "packages_kspace.h"
#undef KSpaceStyle
#undef KSPACE_CLASS
#define MINIMIZE_CLASS
#define MinimizeStyle(key,Class)                \
  pkg_lists->minimize_styles[#key] = PACKAGE;
#include "packages_minimize.h"
#undef MinimizeStyle
#undef MINIMIZE_CLASS
#define PAIR_CLASS
#define PairStyle(key,Class)                    \
  pkg_lists->pair_styles[#key] = PACKAGE;
#include "packages_pair.h"
#undef PairStyle
#undef PAIR_CLASS
#define READER_CLASS
#define ReaderStyle(key,Class)                  \
  pkg_lists->reader_styles[#key] = PACKAGE;
#include "packages_reader.h"
#undef ReaderStyle
#undef READER_CLASS
#define REGION_CLASS
#define RegionStyle(key,Class)                  \
  pkg_lists->region_styles[#key] = PACKAGE;
#include "packages_region.h"
#undef RegionStyle
#undef REGION_CLASS
}

bool LAMMPS::is_installed_pkg(const char *pkg)
{
  for (int i=0; installed_packages[i] != NULL; ++i)
    if (strcmp(installed_packages[i],pkg) == 0) return true;

  return false;
}

#define check_for_match(style,list,name)                                \
  if (strcmp(list,#style) == 0) {                                       \
    std::map<std::string,std::string> &styles(pkg_lists-> style ## _styles); \
    if (styles.find(name) != styles.end()) {                            \
      return styles[name].c_str();                                      \
    }                                                                   \
  }

const char *LAMMPS::match_style(const char *style, const char *name)
{
  check_for_match(angle,style,name);
  check_for_match(atom,style,name);
  check_for_match(body,style,name);
  check_for_match(bond,style,name);
  check_for_match(command,style,name);
  check_for_match(compute,style,name);
  check_for_match(dump,style,name);
  check_for_match(fix,style,name);
  check_for_match(compute,style,name);
  check_for_match(improper,style,name);
  check_for_match(integrate,style,name);
  check_for_match(kspace,style,name);
  check_for_match(minimize,style,name);
  check_for_match(pair,style,name);
  check_for_match(reader,style,name);
  check_for_match(region,style,name);
  return NULL;
}

/* ----------------------------------------------------------------------
   help message for command line options and styles present in executable
------------------------------------------------------------------------- */

void _noopt LAMMPS::help()
{
  FILE *fp = screen;
  const char *pager = NULL;

  // if output is "stdout", use a pipe to a pager for paged output.
  // this will avoid the most important help text to rush past the
  // user. scrollback buffers are often not large enough. this is most
  // beneficial to windows users, who are not used to command line.

  if (fp == stdout) {
    pager = getenv("PAGER");
    if (pager == NULL) pager = "more";
#if defined(_WIN32)
    fp = _popen(pager,"w");
#else
    fp = popen(pager,"w");
#endif

    // reset to original state, if pipe command failed
    if (fp == NULL) {
      fp = stdout;
      pager = NULL;
    }
  }

  // general help message about command line and flags

  if (has_git_info) {
    fprintf(fp,"\nLarge-scale Atomic/Molecular Massively Parallel Simulator - "
            LAMMPS_VERSION "\nGit info (%s / %s)\n\n",git_branch, git_descriptor);
  } else {
    fprintf(fp,"\nLarge-scale Atomic/Molecular Massively Parallel Simulator - "
            LAMMPS_VERSION "\n\n");
  }
  fprintf(fp,
          "Usage example: %s -var t 300 -echo screen -in in.alloy\n\n"
          "List of command line options supported by this LAMMPS executable:\n\n"
          "-echo none/screen/log/both  : echoing of input script (-e)\n"
          "-help                       : print this help message (-h)\n"
          "-in filename                : read input from file, not stdin (-i)\n"
          "-kokkos on/off ...          : turn KOKKOS mode on or off (-k)\n"
          "-log none/filename          : where to send log output (-l)\n"
          "-mpicolor color             : which exe in a multi-exe mpirun cmd (-m)\n"
          "-nocite                     : disable writing log.cite file (-nc)\n"
          "-package style ...          : invoke package command (-pk)\n"
          "-partition size1 size2 ...  : assign partition sizes (-p)\n"
          "-plog basename              : basename for partition logs (-pl)\n"
          "-pscreen basename           : basename for partition screens (-ps)\n"
          "-restart2data rfile dfile ... : convert restart to data file (-r2data)\n"
          "-restart2dump rfile dgroup dstyle dfile ... \n"
          "                            : convert restart to dump file (-r2dump)\n"
          "-reorder topology-specs     : processor reordering (-r)\n"
          "-screen none/filename       : where to send screen output (-sc)\n"
          "-suffix gpu/intel/opt/omp   : style suffix to apply (-sf)\n"
          "-var varname value          : set index style variable (-v)\n\n",
          exename);


  print_config(fp);
  fprintf(fp,"List of individual style options included in this LAMMPS executable\n\n");

  int pos = 80;
  fprintf(fp,"* Atom styles:\n");
#define ATOM_CLASS
#define AtomStyle(key,Class) print_style(fp,#key,pos);
#include "style_atom.h"
#undef ATOM_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Integrate styles:\n");
#define INTEGRATE_CLASS
#define IntegrateStyle(key,Class) print_style(fp,#key,pos);
#include "style_integrate.h"
#undef INTEGRATE_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Minimize styles:\n");
#define MINIMIZE_CLASS
#define MinimizeStyle(key,Class) print_style(fp,#key,pos);
#include "style_minimize.h"
#undef MINIMIZE_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Pair styles:\n");
#define PAIR_CLASS
#define PairStyle(key,Class) print_style(fp,#key,pos);
#include "style_pair.h"
#undef PAIR_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Bond styles:\n");
#define BOND_CLASS
#define BondStyle(key,Class) print_style(fp,#key,pos);
#include "style_bond.h"
#undef BOND_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Angle styles:\n");
#define ANGLE_CLASS
#define AngleStyle(key,Class) print_style(fp,#key,pos);
#include "style_angle.h"
#undef ANGLE_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Dihedral styles:\n");
#define DIHEDRAL_CLASS
#define DihedralStyle(key,Class) print_style(fp,#key,pos);
#include "style_dihedral.h"
#undef DIHEDRAL_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Improper styles:\n");
#define IMPROPER_CLASS
#define ImproperStyle(key,Class) print_style(fp,#key,pos);
#include "style_improper.h"
#undef IMPROPER_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* KSpace styles:\n");
#define KSPACE_CLASS
#define KSpaceStyle(key,Class) print_style(fp,#key,pos);
#include "style_kspace.h"
#undef KSPACE_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Fix styles\n");
#define FIX_CLASS
#define FixStyle(key,Class) print_style(fp,#key,pos);
#include "style_fix.h"
#undef FIX_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Compute styles:\n");
#define COMPUTE_CLASS
#define ComputeStyle(key,Class) print_style(fp,#key,pos);
#include "style_compute.h"
#undef COMPUTE_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Region styles:\n");
#define REGION_CLASS
#define RegionStyle(key,Class) print_style(fp,#key,pos);
#include "style_region.h"
#undef REGION_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Dump styles:\n");
#define DUMP_CLASS
#define DumpStyle(key,Class) print_style(fp,#key,pos);
#include "style_dump.h"
#undef DUMP_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Command styles\n");
#define COMMAND_CLASS
#define CommandStyle(key,Class) print_style(fp,#key,pos);
#include "style_command.h"
#undef COMMAND_CLASS
  fprintf(fp,"\n\n");

  // close pipe to pager, if active

  if (pager != NULL) pclose(fp);
}

/* ----------------------------------------------------------------------
   print style names in columns
   skip any style that starts with upper-case letter, since internal
------------------------------------------------------------------------- */

void print_style(FILE *fp, const char *str, int &pos)
{
  if (isupper(str[0])) return;

  int len = strlen(str);
  if (pos+len > 80) {
    fprintf(fp,"\n");
    pos = 0;
  }

  if (len < 16) {
    fprintf(fp,"%-16s",str);
    pos += 16;
  } else if (len < 32) {
    fprintf(fp,"%-32s",str);
    pos += 32;
  } else if (len < 48) {
    fprintf(fp,"%-48s",str);
    pos += 48;
  } else if (len < 64) {
    fprintf(fp,"%-64s",str);
    pos += 64;
  } else {
    fprintf(fp,"%-80s",str);
    pos += 80;
  }
}

void LAMMPS::print_config(FILE *fp)
{
  const char *pkg;
  int ncword, ncline = 0;

  const char *infobuf = Info::get_os_info();
  fprintf(fp,"OS: %s\n\n",infobuf);
  delete[] infobuf;

  infobuf = Info::get_compiler_info();
  fprintf(fp,"Compiler: %s with %s\n",infobuf,Info::get_openmp_info());
  delete[] infobuf;
  fprintf(fp,"C++ standard: %s\n",Info::get_cxx_info());

  int major,minor;
  infobuf = Info::get_mpi_info(major,minor);
  fprintf(fp,"MPI v%d.%d: %s\n\n",major,minor,infobuf);

  fputs("Active compile time flags:\n\n",fp);
  if (Info::has_gzip_support()) fputs("-DLAMMPS_GZIP\n",fp);
  if (Info::has_png_support()) fputs("-DLAMMPS_PNG\n",fp);
  if (Info::has_jpeg_support()) fputs("-DLAMMPS_JPEG\n",fp);
  if (Info::has_ffmpeg_support()) fputs("-DLAMMPS_FFMPEG\n",fp);
  if (Info::has_exceptions()) fputs("-DLAMMPS_EXCEPTIONS\n",fp);
#if defined(LAMMPS_BIGBIG)
  fputs("-DLAMMPS_BIGBIG\n",fp);
#elif defined(LAMMPS_SMALLBIG)
  fputs("-DLAMMPS_SMALLBIG\n",fp);
#else // defined(LAMMPS_SMALLSMALL)
  fputs("-DLAMMPS_SMALLSMALL\n",fp);
#endif
  fprintf(fp,"\nsizeof(smallint): %3d-bit\n",(int)sizeof(smallint)*8);
  fprintf(fp,"sizeof(imageint): %3d-bit\n",(int)sizeof(imageint)*8);
  fprintf(fp,"sizeof(tagint):   %3d-bit\n",(int)sizeof(tagint)*8);
  fprintf(fp,"sizeof(bigint):   %3d-bit\n",(int)sizeof(bigint)*8);


  fputs("\nInstalled packages:\n\n",fp);
  for (int i = 0; NULL != (pkg = installed_packages[i]); ++i) {
    ncword = strlen(pkg);
    if (ncline + ncword > 78) {
      ncline = 0;
      fputs("\n",fp);
    }
    fprintf(fp,"%s ",pkg);
    ncline += ncword + 1;
  }
  fputs("\n\n",fp);
}
