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

#include "mpi.h"
#include "string.h"
#include "ctype.h"
#include "lammps.h"
#include "style_angle.h"
#include "style_atom.h"
#include "style_bond.h"
#include "style_command.h"
#include "style_compute.h"
#include "style_dihedral.h"
#include "style_dump.h"
#include "style_fix.h"
#include "style_improper.h"
#include "style_integrate.h"
#include "style_kspace.h"
#include "style_minimize.h"
#include "style_pair.h"
#include "style_region.h"
#include "universe.h"
#include "input.h"
#include "atom.h"
#include "update.h"
#include "neighbor.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "modify.h"
#include "group.h"
#include "output.h"
#include "citeme.h"
#include "accelerator_cuda.h"
#include "accelerator_omp.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   start up LAMMPS
   allocate fundamental classes (memory, error, universe, input)
   parse input switches
   initialize communicators, screen & logfile output
   input is allocated at end after MPI info is setup
------------------------------------------------------------------------- */

LAMMPS::LAMMPS(int narg, char **arg, MPI_Comm communicator)
{
  memory = new Memory(this);
  error = new Error(this);
  universe = new Universe(this,communicator);
  output = NULL;

  screen = NULL;
  logfile = NULL;

  // parse input switches

  int inflag = 0;
  int screenflag = 0;
  int logflag = 0;
  int partscreenflag = 0;
  int partlogflag = 0;
  int cudaflag = -1;
  int restartflag = 0;
  int citeflag = 1;
  int helpflag = 0;

  suffix = NULL;
  suffix_enable = 0;
  char *rfile = NULL;
  char *dfile = NULL;
  int wdfirst,wdlast;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"-partition") == 0 ||
        strcmp(arg[iarg],"-p") == 0) {
      universe->existflag = 1;
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      iarg++;
      while (iarg < narg && arg[iarg][0] != '-') {
        universe->add_world(arg[iarg]);
        iarg++;
      }
    } else if (strcmp(arg[iarg],"-in") == 0 ||
               strcmp(arg[iarg],"-i") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      inflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-screen") == 0 ||
               strcmp(arg[iarg],"-sc") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      screenflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-log") == 0 ||
               strcmp(arg[iarg],"-l") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      logflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-var") == 0 ||
               strcmp(arg[iarg],"-v") == 0) {
      if (iarg+3 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 3;
      while (iarg < narg && arg[iarg][0] != '-') iarg++;
    } else if (strcmp(arg[iarg],"-echo") == 0 ||
               strcmp(arg[iarg],"-e") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 2;
    } else if (strcmp(arg[iarg],"-pscreen") == 0 ||
               strcmp(arg[iarg],"-ps") == 0) {
      if (iarg+2 > narg)
       error->universe_all(FLERR,"Invalid command-line argument");
      partscreenflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-plog") == 0 ||
               strcmp(arg[iarg],"-pl") == 0) {
      if (iarg+2 > narg)
       error->universe_all(FLERR,"Invalid command-line argument");
      partlogflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-cuda") == 0 ||
               strcmp(arg[iarg],"-c") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      if (strcmp(arg[iarg+1],"on") == 0) cudaflag = 1;
      else if (strcmp(arg[iarg+1],"off") == 0) cudaflag = 0;
      else error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 2;
    } else if (strcmp(arg[iarg],"-suffix") == 0 ||
               strcmp(arg[iarg],"-sf") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      delete [] suffix;
      int n = strlen(arg[iarg+1]) + 1;
      suffix = new char[n];
      strcpy(suffix,arg[iarg+1]);
      suffix_enable = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-reorder") == 0 ||
               strcmp(arg[iarg],"-ro") == 0) {
      if (iarg+3 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      if (universe->existflag)
        error->universe_all(FLERR,"Cannot use -reorder after -partition");
      universe->reorder(arg[iarg+1],arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"-restart") == 0 ||
               strcmp(arg[iarg],"-r") == 0) {
      if (iarg+3 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      restartflag = 1;
      rfile = arg[iarg+1];
      dfile = arg[iarg+2];
      iarg += 3;
      // delimit any extra args for the write_data command
      wdfirst = iarg;
      while (iarg < narg && arg[iarg][0] != '-') iarg++;
      wdlast = iarg;
    } else if (strcmp(arg[iarg],"-nocite") == 0 ||
               strcmp(arg[iarg],"-nc") == 0) {
      citeflag = 0;
      iarg++;
    } else if (strcmp(arg[iarg],"-help") == 0 ||
               strcmp(arg[iarg],"-h") == 0) {
      if (iarg+1 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      helpflag = 1;
      iarg += 1;
    } else error->universe_all(FLERR,"Invalid command-line argument");
  }

  // if no partition command-line switch, universe is one world with all procs

  if (universe->existflag == 0) universe->add_world(NULL);

  // sum of procs in all worlds must equal total # of procs

  if (!universe->consistent())
    error->universe_all(FLERR,"Processor partitions are inconsistent");

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
    infile = NULL;

    if (universe->me == 0) {
      if (inflag == 0) infile = stdin;
      else infile = fopen(arg[inflag],"r");
      if (infile == NULL) {
        char str[128];
        sprintf(str,"Cannot open input script %s",arg[inflag]);
        error->one(FLERR,str);
      }
    }

    if (universe->me == 0) {
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
         sprintf(str,"%s.%d",arg[screenflag],universe->iworld);
         screen = fopen(str,"w");
         if (screen == NULL) error->one(FLERR,"Cannot open screen file");
       }
      else if (strcmp(arg[partscreenflag],"none") == 0)
        screen = NULL;
      else {
        char str[128];
        sprintf(str,"%s.%d",arg[partscreenflag],universe->iworld);
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
         sprintf(str,"%s.%d",arg[logflag],universe->iworld);
         logfile = fopen(str,"w");
         if (logfile == NULL) error->one(FLERR,"Cannot open logfile");
       }
      else if (strcmp(arg[partlogflag],"none") == 0)
        logfile = NULL;
      else {
        char str[128];
        sprintf(str,"%s.%d",arg[partlogflag],universe->iworld);
        logfile = fopen(str,"w");
        if (logfile == NULL) error->one(FLERR,"Cannot open logfile");
      } else logfile = NULL;

    if (me == 0) {
      infile = fopen(arg[inflag],"r");
      if (infile == NULL) {
        char str[128];
        sprintf(str,"Cannot open input script %s",arg[inflag]);
        error->one(FLERR,str);
      }
    } else infile = NULL;

    // screen and logfile messages for universe and world

    if (universe->me == 0) {
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

    if (me == 0) {
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

  // check datatype settings in lmptype.h

  if (sizeof(smallint) != sizeof(int))
    error->all(FLERR,"Smallint setting in lmptype.h is invalid");
  if (sizeof(imageint) < sizeof(smallint))
    error->all(FLERR,"Imageint setting in lmptype.h is invalid");
  if (sizeof(tagint) < sizeof(smallint))
    error->all(FLERR,"Tagint setting in lmptype.h is invalid");
  if (sizeof(bigint) < sizeof(imageint) || sizeof(bigint) < sizeof(tagint))
    error->all(FLERR,"Bigint setting in lmptype.h is invalid");

  int mpisize;
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

  // create CUDA class if USER-CUDA installed, unless explicitly switched off
  // instantiation creates dummy CUDA class if USER-CUDA is not installed

  if (cudaflag == 0) {
    cuda = NULL;
  } else if (cudaflag == 1) {
    cuda = new Cuda(this);
    if (!cuda->cuda_exists)
      error->all(FLERR,"Cannot use -cuda on without USER-CUDA installed");
  } else {
    cuda = new Cuda(this);
    if (!cuda->cuda_exists) {
      delete cuda;
      cuda = NULL;
    }
  }

  int me;
  MPI_Comm_rank(world,&me);
  if (cuda && me == 0) error->message(FLERR,"USER-CUDA mode is enabled");

  // allocate CiteMe class if enabled

  if (citeflag) citeme = new CiteMe(this);
  else citeme = NULL;

  // allocate input class now that MPI is fully setup

  input = new Input(this,narg,arg);

  // allocate top-level classes

  create();
  post_create();

  // if helpflag set, print help and quit

  if (helpflag) {
    if (universe->me == 0 && screen) help();
    error->done();
  }

  // if restartflag set, process 2 command and quit
  // add args between wdfirst and wdlast to write_data
  // also add "noinit" to prevent write_data from doing system init

  if (restartflag) {
    char cmd[128];
    sprintf(cmd,"read_restart %s\n",rfile);
    input->one(cmd);
    sprintf(cmd,"write_data %s",dfile);
    for (iarg = wdfirst; iarg < wdlast; iarg++)
      sprintf(&cmd[strlen(cmd)]," %s",arg[iarg]);
    strcat(cmd," noinit\n");
    input->one(cmd);
    error->done();
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
  destroy();

  delete citeme;

  if (universe->nworlds == 1) {
    if (logfile) fclose(logfile);
  } else {
    if (screen && screen != stdout) fclose(screen);
    if (logfile) fclose(logfile);
    if (universe->ulogfile) fclose(universe->ulogfile);
  }

  if (world != universe->uworld) MPI_Comm_free(&world);

  delete cuda;
  delete [] suffix;

  delete input;
  delete universe;
  delete error;
  delete memory;
}

/* ----------------------------------------------------------------------
   allocate single instance of top-level classes
   fundamental classes are allocated in constructor
   some classes have package variants
------------------------------------------------------------------------- */

void LAMMPS::create()
{
  // Comm class must be created before Atom class
  // so that nthreads is defined when create_avec invokes grow()

  if (cuda) comm = new CommCuda(this);
  else comm = new Comm(this);

  if (cuda) neighbor = new NeighborCuda(this);
  else neighbor = new Neighbor(this);

  if (cuda) domain = new DomainCuda(this);
#ifdef LMP_USER_OMP
  else domain = new DomainOMP(this);
#else
  else domain = new Domain(this);
#endif

  atom = new Atom(this);
  atom->create_avec("atomic",0,NULL,suffix);

  group = new Group(this);
  force = new Force(this);    // must be after group, to create temperature

  if (cuda) modify = new ModifyCuda(this);
  else modify = new Modify(this);

  output = new Output(this);  // must be after group, so "all" exists
                              // must be after modify so can create Computes
  update = new Update(this);  // must be after output, force, neighbor
  timer = new Timer(this);
}

/* ----------------------------------------------------------------------
   invoke package-specific setup commands
   called from LAMMPS constructor and after clear() command
   only invoke if suffix is set and enabled
------------------------------------------------------------------------- */

void LAMMPS::post_create()
{
  if (suffix && suffix_enable) {
    if (strcmp(suffix,"gpu") == 0) input->one("package gpu force/neigh 0 0 1");
    if (strcmp(suffix,"omp") == 0) input->one("package omp *");
  }
}

/* ----------------------------------------------------------------------
   initialize top-level classes
   do not initialize Timer class, other classes like Run() do that explicitly
------------------------------------------------------------------------- */

void LAMMPS::init()
{
  if (cuda) cuda->accelerator(0,NULL);

  update->init();
  force->init();         // pair must come after update due to minimizer
  domain->init();
  atom->init();          // atom must come after force and domain
                         //   atom deletes extra array
                         //   used by fix shear_history::unpack_restart()
                         //   when force->pair->gran_history creates fix ??
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
  delete neighbor;
  delete comm;
  delete force;
  delete group;
  delete output;
  delete modify;          // modify must come after output, force, update
                          //   since they delete fixes
  delete domain;          // domain must come after modify
                          //   since fix destructors access domain
  delete atom;            // atom must come after modify, neighbor
                          //   since fixes delete callbacks in atom
  delete timer;

  modify = NULL;          // necessary since input->variable->varreader
                          // will be destructed later
}

/* ----------------------------------------------------------------------
   help message for command line options and styles present in executable
------------------------------------------------------------------------- */

void LAMMPS::help()
{
  fprintf(screen,
          "\nCommand line options:\n\n"
          "-cuda on/off                : turn CUDA mode on or off (-c)\n"
          "-echo none/screen/log/both  : echoing of input script (-e)\n"
          "-in filename                : read input from file, not stdin (-i)\n"
          "-help                       : print this help message (-h)\n"
          "-log none/filename          : where to send log output (-l)\n"
          "-nocite                     : disable writing log.cite file (-nc)\n"
          "-partition size1 size2 ...  : assign partition sizes (-p)\n"
          "-plog basename              : basename for partition logs (-pl)\n"
          "-pscreen basename           : basename for partition screens (-ps)\n"
          "-reorder topology-specs     : processor reordering (-r)\n"
          "-screen none/filename       : where to send screen output (-sc)\n"
          "-suffix cuda/gpu/opt/omp    : style suffix to apply (-sf)\n"
          "-var varname value          : set index style variable (-v)\n\n");
  
  fprintf(screen,"Style options compiled with this executable\n\n");

  int pos = 80;
  fprintf(screen,"* Atom styles:\n");
#define ATOM_CLASS
#define AtomStyle(key,Class) print_style(#key,pos);
#include "style_atom.h"
#undef ATOM_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Integrate styles:\n");
#define INTEGRATE_CLASS
#define IntegrateStyle(key,Class) print_style(#key,pos);
#include "style_integrate.h"
#undef INTEGRATE_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Minimize styles:\n");
#define MINIMIZE_CLASS
#define MinimizeStyle(key,Class) print_style(#key,pos);
#include "style_minimize.h"
#undef MINIMIZE_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Pair styles:\n");
#define PAIR_CLASS
#define PairStyle(key,Class) print_style(#key,pos);
#include "style_pair.h"
#undef PAIR_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Bond styles:\n");
#define BOND_CLASS
#define BondStyle(key,Class) print_style(#key,pos);
#include "style_bond.h"
#undef BOND_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Angle styles:\n");
#define ANGLE_CLASS
#define AngleStyle(key,Class) print_style(#key,pos);
#include "style_angle.h"
#undef ANGLE_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Dihedral styles:\n");
#define DIHEDRAL_CLASS
#define DihedralStyle(key,Class) print_style(#key,pos);
#include "style_dihedral.h"
#undef DIHEDRAL_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Improper styles:\n");
#define IMPROPER_CLASS
#define ImproperStyle(key,Class) print_style(#key,pos);
#include "style_improper.h"
#undef IMPROPER_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* KSpace styles:\n");
#define KSPACE_CLASS
#define KSpaceStyle(key,Class) print_style(#key,pos);
#include "style_kspace.h"
#undef KSPACE_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Fix styles\n");
#define FIX_CLASS
#define FixStyle(key,Class) print_style(#key,pos);
#include "style_fix.h"
#undef FIX_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Compute styles:\n");
#define COMPUTE_CLASS
#define ComputeStyle(key,Class) print_style(#key,pos);
#include "style_compute.h"
#undef COMPUTE_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Region styles:\n");
#define REGION_CLASS
#define RegionStyle(key,Class) print_style(#key,pos);
#include "style_region.h"
#undef REGION_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Dump styles:\n");
#define DUMP_CLASS
#define DumpStyle(key,Class) print_style(#key,pos);
#include "style_dump.h"
#undef DUMP_CLASS
  fprintf(screen,"\n\n");

  pos = 80;
  fprintf(screen,"* Command styles\n");
#define COMMAND_CLASS
#define CommandStyle(key,Class) print_style(#key,pos);
#include "style_command.h"
#undef COMMAND_CLASS
  fprintf(screen,"\n");
}

/* ----------------------------------------------------------------------
   print style names in columns
   skip any style that starts with upper-case letter, since internal
------------------------------------------------------------------------- */

void LAMMPS::print_style(const char *str, int &pos)
{
  if (isupper(str[0])) return;

  int len = strlen(str);
  if (pos+len > 80) { 
    fprintf(screen,"\n");
    pos = 0;
  }

  if (len < 16) {
    fprintf(screen,"%-16s",str);
    pos += 16;
  } else if (len < 32) {
    fprintf(screen,"%-32s",str);
    pos += 32;
  } else if (len < 48) {
    fprintf(screen,"%-48s",str);
    pos += 48;
  } else if (len < 64) {
    fprintf(screen,"%-64s",str);
    pos += 64;
  } else {
    fprintf(screen,"%-80s",str);
    pos += 80;
  }
}
