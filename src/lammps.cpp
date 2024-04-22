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

#include "lammps.h"

#include "style_angle.h"     // IWYU pragma: keep
#include "style_atom.h"      // IWYU pragma: keep
#include "style_body.h"      // IWYU pragma: keep
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
#include "style_reader.h"    // IWYU pragma: keep
#include "style_region.h"    // IWYU pragma: keep

#include "accelerator_kokkos.h"
#include "accelerator_omp.h"
#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "comm_brick.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "info.h"
#include "input.h"
#include "lmppython.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "suffix.h"
#include "timer.h"
#include "universe.h"
#include "update.h"
#include "variable.h"
#include "version.h"

#if defined(LMP_PLUGIN)
#include "plugin.h"
#endif

#include <cctype>
#include <cmath>
#include <cstring>
#include <map>

#include "lmpinstalledpkgs.h"
#include "lmpgitversion.h"

#if defined(LAMMPS_UPDATE)
#define UPDATE_STRING " - " LAMMPS_UPDATE
#else
#define UPDATE_STRING ""
#endif

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
 * \brief LAMMPS simulation instance
 *
 * The LAMMPS class contains pointers of all constituent class instances
 * and global variables that are used by a LAMMPS simulation. Its contents
 * represent the entire state of the simulation.
 *
 * The LAMMPS class manages the components of an MD simulation by creating,
 * deleting, and initializing instances of the classes it is composed of,
 * processing command line flags, and providing access to some global properties.
 * The specifics of setting up and running a simulation are handled by the
 * individual component class instances. */

/** Create a LAMMPS simulation instance
 *
 * \param args list of arguments
 * \param communicator MPI communicator used by this LAMMPS instance
 */
LAMMPS::LAMMPS(argv & args, MPI_Comm communicator) :
  LAMMPS(args.size(), argv_pointers(args).data(), communicator) {
}

/** Create a LAMMPS simulation instance
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
  memory(nullptr), error(nullptr), universe(nullptr), input(nullptr), atom(nullptr),
  update(nullptr), neighbor(nullptr), comm(nullptr), domain(nullptr), force(nullptr),
  modify(nullptr), group(nullptr), output(nullptr), timer(nullptr), kokkos(nullptr),
  atomKK(nullptr), memoryKK(nullptr), python(nullptr), citeme(nullptr)
{
  memory = new Memory(this);
  error = new Error(this);
  universe = new Universe(this,communicator);

  version = (const char *) LAMMPS_VERSION;
  num_ver = utils::date2num(version);
  restart_ver = -1;

  // append git descriptor info to update string when compiling development or maintenance version

  std::string update_string = UPDATE_STRING;
  if (has_git_info() && ((update_string == " - Development") || (update_string == " - Maintenance")))
    update_string += fmt::format(" - {}", git_descriptor());

  external_comm = MPI_COMM_NULL;
  mdicomm = nullptr;

  skiprunflag = 0;

  screen = nullptr;
  logfile = nullptr;
  infile = nullptr;

  initclock = platform::walltime();

  init_pkg_lists();

#if defined(LMP_PYTHON) && defined(_WIN32)
  // If the LAMMPSHOME environment variable is set, it should point
  // to the location of the LAMMPS installation tree where we bundle
  // the matching Python installation for use with the PYTHON package.
  // This is currently only used on Windows with the Windows installer packages
  const char *lmpenv = getenv("LAMMPSHOME");
  if (lmpenv) {
    platform::putenv(fmt::format("PYTHONHOME={}",lmpenv));
  }
#endif

  // check if -mpicolor is first arg
  // if so, then 2 or more apps were launched with one mpirun command
  //   this means passed communicator (e.g. MPI_COMM_WORLD) is bigger than LAMMPS
  //   universe communicator needs to shrink to be just LAMMPS
  // syntax: -mpicolor color
  //   color = integer for this app, different than any other app(s)
  // do the following:
  //   perform an MPI_Comm_split() to create a new LAMMPS-only subcomm
  //   NOTE: this assumes other app(s) make same call, else will hang!
  //   re-create universe with subcomm
  //   store comm that all apps belong to in external_comm

  int iarg = 1;
  if (narg-iarg >= 2 && (strcmp(arg[iarg],"-mpicolor") == 0 ||
                         strcmp(arg[iarg],"-m") == 0)) {
    int me,nprocs;
    MPI_Comm_rank(communicator,&me);
    MPI_Comm_size(communicator,&nprocs);
    int color = atoi(arg[iarg+1]);
    MPI_Comm subcomm;
    MPI_Comm_split(communicator,color,me,&subcomm);
    external_comm = communicator;
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
  int citescreen = CiteMe::TERSE;
  int citelogfile = CiteMe::VERBOSE;
  char *citefile = nullptr;
  int helpflag = 0;
  int nonbufflag = 0;

  suffix = suffix2 = nullptr;
  suffix_enable = 0;
  pair_only_flag = 0;
  if (arg) exename = utils::strdup(arg[0]);
  else exename = nullptr;
  packargs = nullptr;
  num_package = 0;
  char *restartfile = nullptr;
  int wfirst,wlast;
  int kkfirst,kklast;

  int npack = 0;
  int *pfirst = nullptr;
  int *plast = nullptr;

  iarg = 1;
  while (iarg < narg) {

    if (strcmp(arg[iarg],"-cite") == 0 ||
               strcmp(arg[iarg],"-c") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");

      if (strcmp(arg[iarg+1],"both") == 0) {
        citescreen = CiteMe::VERBOSE;
        citelogfile = CiteMe::VERBOSE;
        citefile = nullptr;
      } else if (strcmp(arg[iarg+1],"none") == 0) {
        citescreen = CiteMe::TERSE;
        citelogfile = CiteMe::TERSE;
        citefile = nullptr;
      } else if (strcmp(arg[iarg+1],"screen") == 0) {
        citescreen = CiteMe::VERBOSE;
        citelogfile = CiteMe::TERSE;
        citefile = nullptr;
      } else if (strcmp(arg[iarg+1],"log") == 0) {
        citescreen = CiteMe::TERSE;
        citelogfile = CiteMe::VERBOSE;
        citefile = nullptr;
      } else {
        citescreen = CiteMe::TERSE;
        citelogfile = CiteMe::TERSE;
        citefile = arg[iarg+1];
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"-echo") == 0 ||
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
      inflag = -1;              // skip inflag check
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
      const std::string kokkosarg = arg[iarg+1];
      if ((kokkosarg == "on") || (kokkosarg == "yes") || (kokkosarg == "true"))
        kokkosflag = 1;
      else if ((kokkosarg == "off") || (kokkosarg == "no") || (kokkosarg == "false"))
        kokkosflag = 0;
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

    } else if (strcmp(arg[iarg],"-mpicolor") == 0 ||
               strcmp(arg[iarg],"-m") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      if (iarg != 1) error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 2;

    } else if (strcmp(arg[iarg],"-nocite") == 0 ||
               strcmp(arg[iarg],"-nc") == 0) {
      citeflag = 0;
      iarg++;

    } else if (strcmp(arg[iarg],"-nonbuf") == 0 ||
               strcmp(arg[iarg],"-nb") == 0) {
      nonbufflag = 1;
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
      inflag = -1;               // skip inflag check
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
      inflag = -1;               // skip inflag check
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

    } else if (strcmp(arg[iarg],"-skiprun") == 0 ||
               strcmp(arg[iarg],"-sr") == 0) {
      skiprunflag = 1;
      ++iarg;

    } else if (strcmp(arg[iarg],"-suffix") == 0 ||
               strcmp(arg[iarg],"-sf") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      delete[] suffix;
      delete[] suffix2;
      suffix = suffix2 = nullptr;
      suffix_enable = 1;
      // hybrid option to set fall-back for suffix2
      if (strcmp(arg[iarg+1],"hybrid") == 0) {
        if (iarg+4 > narg)
          error->universe_all(FLERR,"Invalid command-line argument");
        suffix = utils::strdup(arg[iarg+2]);
        suffix2 = utils::strdup(arg[iarg+3]);
        iarg += 4;
      } else {
        suffix = utils::strdup(arg[iarg+1]);
        iarg += 2;
      }

    } else if (strcmp(arg[iarg],"-var") == 0 ||
               strcmp(arg[iarg],"-v") == 0) {
      if (iarg+3 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 3;
      while (iarg < narg && arg[iarg][0] != '-') iarg++;

    } else {
      error->universe_all(FLERR, fmt::format("Invalid command-line argument: {}", arg[iarg]) );
    }
  }

  // if no partition command-line switch, universe is one world with all procs

  if (universe->existflag == 0) universe->add_world(nullptr);

  // sum of procs in all worlds must equal total # of procs

  if (!universe->consistent())
    error->universe_all(FLERR,"Processor partitions do not match number of allocated processors");

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
      universe->uscreen = nullptr;
    else {
      universe->uscreen = fopen(arg[screenflag],"w");
      if (universe->uscreen == nullptr)
        error->universe_one(FLERR,fmt::format("Cannot open universe screen file {}: {}",
                                              arg[screenflag],utils::getsyserror()));
    }
    if (logflag == 0) {
      if (helpflag == 0) {
        universe->ulogfile = fopen("log.lammps","w");
        if (universe->ulogfile == nullptr)
          error->universe_warn(FLERR,"Cannot open log.lammps for writing: "
                               + utils::getsyserror());
      }
    } else if (strcmp(arg[logflag],"none") == 0)
      universe->ulogfile = nullptr;
    else {
      universe->ulogfile = fopen(arg[logflag],"w");
      if (universe->ulogfile == nullptr)
        error->universe_one(FLERR,fmt::format("Cannot open universe log file {}: {}",
                                              arg[logflag],utils::getsyserror()));
    }
  }

  if (universe->me > 0) {
    if (screenflag == 0) universe->uscreen = stdout;
    else universe->uscreen = nullptr;
    universe->ulogfile = nullptr;
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
      if (inflag <= 0) infile = stdin;
      else if (strcmp(arg[inflag], "none") == 0) infile = stdin;
      else infile = fopen(arg[inflag],"r");
      if (infile == nullptr)
        error->one(FLERR,"Cannot open input script {}: {}", arg[inflag], utils::getsyserror());
      if (!helpflag)
        utils::logmesg(this,"LAMMPS ({}{})\n", version, update_string);

     // warn against using I/O redirection in parallel runs
      if ((inflag == 0) && (universe->nprocs > 1))
        error->warning(FLERR, "Using I/O redirection is unreliable with parallel runs. "
                       "Better to use the -in switch to read input files.");
      utils::flush_buffers(this);
    }

  // universe is one or more worlds, as setup by partition switch
  // split universe communicator into separate world communicators
  // set world screen, logfile, communicator, infile
  // open input script

  } else {
    int me;
    MPI_Comm_split(universe->uworld,universe->iworld,0,&world);
    MPI_Comm_rank(world,&me);

    screen = logfile = infile = nullptr;
    if (me == 0) {
      std::string str;
      if (partscreenflag == 0) {
        if (screenflag == 0) {
          str = fmt::format("screen.{}",universe->iworld);
          screen = fopen(str.c_str(),"w");
          if (screen == nullptr)
            error->one(FLERR,"Cannot open screen file {}: {}",str,utils::getsyserror());
        } else if (strcmp(arg[screenflag],"none") == 0) {
          screen = nullptr;
        } else {
          str = fmt::format("{}.{}",arg[screenflag],universe->iworld);
          screen = fopen(str.c_str(),"w");
          if (screen == nullptr)
            error->one(FLERR,"Cannot open screen file {}: {}",arg[screenflag],utils::getsyserror());
        }
      } else if (strcmp(arg[partscreenflag],"none") == 0) {
        screen = nullptr;
      } else {
        str = fmt::format("{}.{}",arg[partscreenflag],universe->iworld);
        screen = fopen(str.c_str(),"w");
        if (screen == nullptr)
          error->one(FLERR,"Cannot open screen file {}: {}",str,utils::getsyserror());
      }

      if (partlogflag == 0) {
        if (logflag == 0) {
          str = fmt::format("log.lammps.{}",universe->iworld);
          logfile = fopen(str.c_str(),"w");
          if (logfile == nullptr)
            error->one(FLERR,"Cannot open logfile {}: {}",str, utils::getsyserror());
        } else if (strcmp(arg[logflag],"none") == 0) {
          logfile = nullptr;
        } else {
          str = fmt::format("{}.{}",arg[logflag],universe->iworld);
          logfile = fopen(str.c_str(),"w");
          if (logfile == nullptr)
            error->one(FLERR,"Cannot open logfile {}: {}",str, utils::getsyserror());
        }
      } else if (strcmp(arg[partlogflag],"none") == 0) {
        logfile = nullptr;
      } else {
        str = fmt::format("{}.{}",arg[partlogflag],universe->iworld);
        logfile = fopen(str.c_str(),"w");
        if (logfile == nullptr)
          error->one(FLERR,"Cannot open logfile {}: {}",str, utils::getsyserror());
      }

      if (strcmp(arg[inflag], "none") != 0) {
        infile = fopen(arg[inflag],"r");
        if (infile == nullptr)
          error->one(FLERR,"Cannot open input script {}: {}",arg[inflag], utils::getsyserror());
      }
    }

    // make all screen and logfile output unbuffered for debugging crashes

    if (nonbufflag) {
      if (universe->uscreen) setbuf(universe->uscreen, nullptr);
      if (universe->ulogfile) setbuf(universe->ulogfile, nullptr);
      if (screen) setbuf(screen, nullptr);
      if (logfile) setbuf(logfile, nullptr);
    }

    // screen and logfile messages for universe and world

    if ((universe->me == 0) && (!helpflag)) {
      const char fmt[] = "LAMMPS ({})\nRunning on {} partitions of processors\n";
      if (universe->uscreen)
        fmt::print(universe->uscreen,fmt,version,universe->nworlds);

      if (universe->ulogfile)
        fmt::print(universe->ulogfile,fmt,version,universe->nworlds);
    }

    if ((me == 0) && (!helpflag))
      utils::logmesg(this,"LAMMPS ({})\nProcessor partition = {}\n", version, universe->iworld);
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
      error->all(FLERR,"MPI_LMP_TAGINT and tagint in lmptype.h are not compatible");
  MPI_Type_size(MPI_LMP_BIGINT,&mpisize);
  if (mpisize != sizeof(bigint))
      error->all(FLERR,"MPI_LMP_BIGINT and bigint in lmptype.h are not compatible");

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

  kokkos = nullptr;
  if (kokkosflag == 1) {
    kokkos = new KokkosLMP(this,kklast-kkfirst,&arg[kkfirst]);
    if (!kokkos->kokkos_exists)
      error->all(FLERR,"Cannot use -kokkos on without KOKKOS installed");
  }

  // allocate CiteMe class if enabled

  if (citeflag) citeme = new CiteMe(this,citescreen,citelogfile,citefile);
  else citeme = nullptr;

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
        packargs[i][j] = utils::strdup(arg[pfirst[i]+j]);
      packargs[i][n] = nullptr;
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
    std::string cmd = fmt::format("read_restart {}",restartfile);
    if (restartremap) cmd += " remap\n";
    input->one(cmd);
    if (restart2data) cmd = "write_data ";
    else cmd = "write_dump";
    for (iarg = wfirst; iarg < wlast; iarg++)
       cmd += fmt::format(" {}", arg[iarg]);
    if (restart2data) cmd += " noinit";
    input->one(cmd);
    error->done(0);
  }
}

/** Shut down a LAMMPS simulation instance
 *
 * The LAMMPS destructor shuts down the simulation by deleting top-level class
 * instances, closing screen and log files for the global instance (aka "world")
 * and files and MPI communicators in sub-partitions ("universes"). Then it
 * deletes the fundamental class instances and copies of data inside the class.
 */
LAMMPS::~LAMMPS() noexcept(false)
{
  const int me = comm->me;

  delete citeme;
  destroy();

  if (num_package) {
    for (int i = 0; i < num_package; i++) {
      for (char **ptr = packargs[i]; *ptr != nullptr; ++ptr)
        delete[] *ptr;
      delete[] packargs[i];
    }
    delete[] packargs;
  }
  num_package = 0;
  packargs = nullptr;

  double totalclock = platform::walltime() - initclock;
  if ((me == 0) && (screen || logfile)) {
    int seconds = fmod(totalclock,60.0);
    totalclock  = (totalclock - seconds) / 60.0;
    int minutes = fmod(totalclock,60.0);
    int hours = (totalclock - minutes) / 60.0;
    utils::logmesg(this, "Total wall time: {}:{:02d}:{:02d}\n", hours, minutes, seconds);
  }

  if (universe->nworlds == 1) {
    if (screen && screen != stdout) fclose(screen);
    if (logfile) fclose(logfile);
    logfile = nullptr;
    if (screen != stdout) screen = nullptr;
  } else {
    if (screen && screen != stdout) fclose(screen);
    if (logfile) fclose(logfile);
    if (universe->ulogfile) fclose(universe->ulogfile);
    logfile = nullptr;
    if (screen != stdout) screen = nullptr;
  }

  if (infile && infile != stdin) fclose(infile);

  if (world != universe->uworld) MPI_Comm_free(&world);

  delete python;
  delete kokkos;
  delete[] suffix;
  delete[] suffix2;

  // free the MPI comm created by -mpicolor cmdline arg processed in constructor
  // it was passed to universe as if original universe world
  // may have been split later by partitions, universe will free the splits
  // free a copy of uorig here, so check in universe destructor will still work

  MPI_Comm copy = universe->uorig;
  if (external_comm != MPI_COMM_NULL) MPI_Comm_free(&copy);

  delete input;
  delete universe;
  delete error;
  delete memory;

  delete pkg_lists;
  delete[] exename;
}

/* ----------------------------------------------------------------------
   allocate single instance of top-level classes
   fundamental classes are allocated in constructor
   some classes have package variants
------------------------------------------------------------------------- */

void LAMMPS::create()
{
  force = nullptr;         // Domain->Lattice checks if Force exists

  // Comm class must be created before Atom class
  // so that nthreads is defined when create_avec invokes grow()

  if (kokkos) comm = new CommBrickKokkos(this);
  else comm = new CommBrick(this);

  if (kokkos) neighbor = new NeighborKokkos(this);
  else neighbor = new Neighbor(this);

  if (kokkos) domain = new DomainKokkos(this);
#ifdef LMP_OPENMP
  else domain = new DomainOMP(this);
#else
  else domain = new Domain(this);
#endif

  if (kokkos) atom = new AtomKokkos(this);
  else atom = new Atom(this);

  if (kokkos)
    atom->create_avec("atomic/kk",0,nullptr,1);
  else
    atom->create_avec("atomic",0,nullptr,1);

  group = new Group(this);
  force = new Force(this);    // must be after group, to create temperature

  if (kokkos) modify = new ModifyKokkos(this);
  else modify = new Modify(this);

  output = new Output(this);  // must be after group, so "all" exists
                              // must be after modify so can create Computes
  update = new Update(this);  // must be after output, force, neighbor
  timer = new Timer(this);

  python = new Python(this);

  // auto-load plugins
#if defined(LMP_PLUGIN)
  plugin_auto_load(this);
#endif
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
  if (skiprunflag) input->one("timer timeout 0 every 1");

  // Don't unnecessarily reissue a package command via suffix
  int package_issued = Suffix::NONE;

  // default package command triggered by "-k on"

  if (kokkos && kokkos->kokkos_exists) input->one("package kokkos");

  // invoke any command-line package commands

  if (num_package) {
    std::string str;
    for (int i = 0; i < num_package; i++) {
      str = "package";
      char *pkg_name = *(packargs[i]);
      if (pkg_name != nullptr) {
        if (strcmp("gpu", pkg_name) == 0) package_issued |= Suffix::GPU;
        if (strcmp("omp", pkg_name) == 0) package_issued |= Suffix::OMP;
        if (strcmp("intel", pkg_name) == 0) package_issued |= Suffix::INTEL;
      }
      for (char **ptr = packargs[i]; *ptr != nullptr; ++ptr) {
        str += " ";
        str += *ptr;
      }
      input->one(str);
    }
  }

  // check that KOKKOS package classes were instantiated
  // check that GPU, INTEL, OPENMP fixes were compiled with LAMMPS
  // do not re-issue package command if already issued

  if (suffix_enable) {

    if (strcmp(suffix,"gpu") == 0 && !modify->check_package("GPU"))
      error->all(FLERR,"Using suffix gpu without GPU package installed");
    if (strcmp(suffix,"intel") == 0 && !modify->check_package("INTEL"))
      error->all(FLERR,"Using suffix intel without INTEL package installed");
    if (strcmp(suffix,"kk") == 0 &&
        (kokkos == nullptr || kokkos->kokkos_exists == 0))
      error->all(FLERR,"Using suffix kk without KOKKOS package enabled");
    if (strcmp(suffix,"omp") == 0 && !modify->check_package("OMP"))
      error->all(FLERR,"Using suffix omp without OPENMP package installed");

    if (strcmp(suffix,"gpu") == 0 && !(package_issued & Suffix::GPU))
      input->one("package gpu 0");
    if (strcmp(suffix,"intel") == 0 && !(package_issued & Suffix::INTEL))
      input->one("package intel 1");
    if (strcmp(suffix,"omp") == 0 && !(package_issued & Suffix::OMP))
      input->one("package omp 0");

    if (suffix2) {
      if (strcmp(suffix2,"gpu") == 0 && !(package_issued & Suffix::GPU))
        input->one("package gpu 0");
      if (strcmp(suffix2,"intel") == 0 && !(package_issued & Suffix::INTEL))
        input->one("package intel 1");
      if (strcmp(suffix2,"omp") == 0 && !(package_issued & Suffix::OMP))
        input->one("package omp 0");
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
  // must wipe out all plugins first, if configured
#if defined(LMP_PLUGIN)
  plugin_clear(this);
#endif

  delete update;
  update = nullptr;

  delete neighbor;
  neighbor = nullptr;

  delete force;
  force = nullptr;

  delete group;
  group = nullptr;

  delete output;
  output = nullptr;

  // undefine atomfile variables because they use a fix for backing storage
  input->variable->purge_atomfile();

  delete modify;          // modify must come after output, force, update
                          //   since they delete fixes
  modify = nullptr;

  delete comm;            // comm must come after modify
                          //   since fix destructors may access comm
  comm = nullptr;

  delete domain;          // domain must come after modify
                          //   since fix destructors access domain
  domain = nullptr;

  delete atom;            // atom must come after modify, neighbor
                          //   since fixes delete callbacks in atom
  atom = nullptr;

  delete timer;
  timer = nullptr;

  delete python;
  python = nullptr;

  restart_ver = -1;       // reset last restart version id
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

/** Return true if a LAMMPS package is enabled in this binary
 *
 * \param pkg name of package
 * \return true if yes, else false
 */
bool LAMMPS::is_installed_pkg(const char *pkg)
{
  for (int i=0; installed_packages[i] != nullptr; ++i)
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

/** \brief Return name of package that a specific style belongs to
 *
 * This function checks the given name against all list of styles
 * for all type of styles and if the name and the style match, it
 * returns which package this style belongs to.
 *
 * \param style Type of style (e.g. atom, pair, fix, etc.)
 * \param name Name of style
 * \return Name of the package this style is part of
 */
const char *LAMMPS::match_style(const char *style, const char *name)
{
  check_for_match(angle,style,name);
  check_for_match(atom,style,name);
  check_for_match(body,style,name);
  check_for_match(bond,style,name);
  check_for_match(command,style,name);
  check_for_match(compute,style,name);
  check_for_match(dihedral,style,name);
  check_for_match(dump,style,name);
  check_for_match(fix,style,name);
  check_for_match(improper,style,name);
  check_for_match(integrate,style,name);
  check_for_match(kspace,style,name);
  check_for_match(minimize,style,name);
  check_for_match(pair,style,name);
  check_for_match(reader,style,name);
  check_for_match(region,style,name);
  return nullptr;
}

#undef check_for_match

/** \brief  Return suffix for non-pair styles depending on pair_only_flag
 *
 * \return  suffix or null pointer
 */
const char *LAMMPS::non_pair_suffix() const
{
  const char *mysuffix;
  if (pair_only_flag) {
#ifdef LMP_KOKKOS_GPU
    if (utils::strmatch(suffix,"^kk")) mysuffix = "kk/host";
    else mysuffix = nullptr;
#else
    mysuffix = nullptr;
#endif
  } else {
    mysuffix = suffix;
  }
  return mysuffix;
}

/* ----------------------------------------------------------------------
   help message for command line options and styles present in executable
------------------------------------------------------------------------- */

void _noopt LAMMPS::help()
{
  FILE *fp = screen;
  const char *pager = nullptr;

  // if output is a console, use a pipe to a pager for paged output.
  // this will avoid the most important help text to rush past the
  // user. scrollback buffers are often not large enough. this is most
  // beneficial to windows users, who are not used to command line.

  int use_pager = platform::is_console(fp);

  // cannot use this with OpenMPI since its console is non-functional

#if defined(OPEN_MPI)
  use_pager = 0;
#endif

  if (use_pager) {
    pager = getenv("PAGER");
    if (pager == nullptr) pager = "more";
    fp = platform::popen(pager,"w");

    // reset to original state, if pipe command failed
    if (fp == nullptr) {
      fp = screen;
      pager = nullptr;
    }
  }

  // general help message about command line and flags

  if (has_git_info()) {
    fprintf(fp,"\nLarge-scale Atomic/Molecular Massively Parallel Simulator - "
            LAMMPS_VERSION UPDATE_STRING "\nGit info (%s / %s)\n\n",git_branch(), git_descriptor());
  } else {
    fprintf(fp,"\nLarge-scale Atomic/Molecular Massively Parallel Simulator - "
            LAMMPS_VERSION UPDATE_STRING "\n\n");
  }
  fprintf(fp,
          "Usage example: %s -var t 300 -echo screen -in in.alloy\n\n"
          "List of command line options supported by this LAMMPS executable:\n\n"
          "-echo none/screen/log/both  : echoing of input script (-e)\n"
          "-help                       : print this help message (-h)\n"
          "-in none/filename           : read input from file or stdin (default) (-i)\n"
          "-kokkos on/off ...          : turn KOKKOS mode on or off (-k)\n"
          "-log none/filename          : where to send log output (-l)\n"
          "-mdi '<mdi flags>'          : pass flags to the MolSSI Driver Interface\n"
          "-mpicolor color             : which exe in a multi-exe mpirun cmd (-m)\n"
          "-cite                       : select citation reminder style (-c)\n"
          "-nocite                     : disable citation reminder (-nc)\n"
          "-nonbuf                     : disable screen/logfile buffering (-nb)\n"
          "-package style ...          : invoke package command (-pk)\n"
          "-partition size1 size2 ...  : assign partition sizes (-p)\n"
          "-plog basename              : basename for partition logs (-pl)\n"
          "-pscreen basename           : basename for partition screens (-ps)\n"
          "-restart2data rfile dfile ... : convert restart to data file (-r2data)\n"
          "-restart2dump rfile dgroup dstyle dfile ... \n"
          "                            : convert restart to dump file (-r2dump)\n"
          "-reorder topology-specs     : processor reordering (-r)\n"
          "-screen none/filename       : where to send screen output (-sc)\n"
          "-skiprun                    : skip loops in run and minimize (-sr)\n"
          "-suffix gpu/intel/kk/opt/omp: style suffix to apply (-sf)\n"
          "-var varname value          : set index style variable (-v)\n\n",
          exename);

  print_config(fp);
  fprintf(fp,"List of individual style options included in this LAMMPS executable\n\n");

  int pos = 80;
  fprintf(fp,"* Atom styles:\n");
#define ATOM_CLASS
#define AtomStyle(key,Class) print_style(fp,#key,pos);
#include "style_atom.h"  // IWYU pragma: keep
#undef ATOM_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Integrate styles:\n");
#define INTEGRATE_CLASS
#define IntegrateStyle(key,Class) print_style(fp,#key,pos);
#include "style_integrate.h"  // IWYU pragma: keep
#undef INTEGRATE_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Minimize styles:\n");
#define MINIMIZE_CLASS
#define MinimizeStyle(key,Class) print_style(fp,#key,pos);
#include "style_minimize.h"  // IWYU pragma: keep
#undef MINIMIZE_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Pair styles:\n");
#define PAIR_CLASS
#define PairStyle(key,Class) print_style(fp,#key,pos);
#include "style_pair.h"  // IWYU pragma: keep
#undef PAIR_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Bond styles:\n");
#define BOND_CLASS
#define BondStyle(key,Class) print_style(fp,#key,pos);
#include "style_bond.h"  // IWYU pragma: keep
#undef BOND_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Angle styles:\n");
#define ANGLE_CLASS
#define AngleStyle(key,Class) print_style(fp,#key,pos);
#include "style_angle.h"  // IWYU pragma: keep
#undef ANGLE_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Dihedral styles:\n");
#define DIHEDRAL_CLASS
#define DihedralStyle(key,Class) print_style(fp,#key,pos);
#include "style_dihedral.h"  // IWYU pragma: keep
#undef DIHEDRAL_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Improper styles:\n");
#define IMPROPER_CLASS
#define ImproperStyle(key,Class) print_style(fp,#key,pos);
#include "style_improper.h"  // IWYU pragma: keep
#undef IMPROPER_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* KSpace styles:\n");
#define KSPACE_CLASS
#define KSpaceStyle(key,Class) print_style(fp,#key,pos);
#include "style_kspace.h"  // IWYU pragma: keep
#undef KSPACE_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Fix styles\n");
#define FIX_CLASS
#define FixStyle(key,Class) print_style(fp,#key,pos);
#include "style_fix.h"  // IWYU pragma: keep
#undef FIX_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Compute styles:\n");
#define COMPUTE_CLASS
#define ComputeStyle(key,Class) print_style(fp,#key,pos);
#include "style_compute.h"  // IWYU pragma: keep
#undef COMPUTE_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Region styles:\n");
#define REGION_CLASS
#define RegionStyle(key,Class) print_style(fp,#key,pos);
#include "style_region.h"  // IWYU pragma: keep
#undef REGION_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Dump styles:\n");
#define DUMP_CLASS
#define DumpStyle(key,Class) print_style(fp,#key,pos);
#include "style_dump.h"  // IWYU pragma: keep
#undef DUMP_CLASS
  fprintf(fp,"\n\n");

  pos = 80;
  fprintf(fp,"* Command styles\n");
#define COMMAND_CLASS
#define CommandStyle(key,Class) print_style(fp,#key,pos);
#include "style_command.h"  // IWYU pragma: keep
#undef COMMAND_CLASS
  fprintf(fp,"\n\n");

  // close pipe to pager, if active

  if (pager != nullptr) platform::pclose(fp);
}

/* ----------------------------------------------------------------------
   print style names in columns
   skip any internal style that starts with an upper-case letter
   also skip "redundant" KOKKOS styles ending in kk/host or kk/device
------------------------------------------------------------------------- */

void print_style(FILE *fp, const char *str, int &pos)
{
  if (isupper(str[0]) || utils::strmatch(str,"/kk/host$")
      || utils::strmatch(str,"/kk/device$")) return;

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

  fmt::print(fp,"OS: {}\n\n",platform::os_info());

  fmt::print(fp,"Compiler: {} with {}\nC++ standard: {}\n",
             platform::compiler_info(),platform::openmp_standard(),
             platform::cxx_standard());

  int major,minor;
  std::string infobuf = platform::mpi_info(major,minor);
  fmt::print(fp,"MPI v{}.{}: {}\n\n",major,minor,infobuf);

  fmt::print(fp,"Accelerator configuration:\n\n{}\n",
             Info::get_accelerator_info());
#if defined(LMP_GPU)
  fmt::print(fp,"Compatible GPU present: {}\n\n",Info::has_gpu_device() ? "yes" : "no");
#endif

  fputs("FFT information:\n\n",fp);
  fputs(Info::get_fft_info().c_str(),fp);

  fputs("\nActive compile time flags:\n\n",fp);
  if (Info::has_gzip_support()) fputs("-DLAMMPS_GZIP\n",fp);
  if (Info::has_png_support()) fputs("-DLAMMPS_PNG\n",fp);
  if (Info::has_jpeg_support()) fputs("-DLAMMPS_JPEG\n",fp);
  if (Info::has_ffmpeg_support()) fputs("-DLAMMPS_FFMPEG\n",fp);
  if (Info::has_fft_single_support()) fputs("-DFFT_SINGLE\n",fp);
#if defined(LAMMPS_BIGBIG)
  fputs("-DLAMMPS_BIGBIG\n",fp);
#elif defined(LAMMPS_SMALLBIG)
  fputs("-DLAMMPS_SMALLBIG\n",fp);
#else // defined(LAMMPS_SMALLSMALL)
  fputs("-DLAMMPS_SMALLSMALL\n",fp);
#endif

  fmt::print(fp,"sizeof(smallint): {}-bit\n"
             "sizeof(imageint): {}-bit\n"
             "sizeof(tagint):   {}-bit\n"
             "sizeof(bigint):   {}-bit\n",
             sizeof(smallint)*8, sizeof(imageint)*8,
             sizeof(tagint)*8, sizeof(bigint)*8);

  if (Info::has_gzip_support()) fmt::print(fp,"\n{}\n",platform::compress_info());

  fputs("\nInstalled packages:\n\n",fp);
  for (int i = 0; nullptr != (pkg = installed_packages[i]); ++i) {
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

/** Create vector of argv string pointers including terminating nullptr element
 *
 * \param args list of arguments
 */
std::vector<char*> LAMMPS::argv_pointers(argv & args){
  std::vector<char*> r;
  r.reserve(args.size()+1);
  for(auto & a : args) {
    r.push_back((char*)a.data());
  }
  r.push_back(nullptr);
  return r;
}
