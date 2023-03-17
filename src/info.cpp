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

/* ----------------------------------------------------------------------
   Contributing authors:  Axel Kohlmeyer (Temple U),
                          Richard Berger (Temple U)
------------------------------------------------------------------------- */

#include "info.h"

#include "accelerator_kokkos.h"
#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "dump.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "input.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "region.h"
#include "update.h"
#include "variable.h"
#include "fmt/chrono.h"

#include <cctype>
#include <cmath>
#include <cstring>
#include <ctime>
#include <map>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#define PSAPI_VERSION 1
#include <windows.h>
#include <cstdint>
#include <psapi.h>
#else
#include <sys/resource.h>
#endif

#if defined(__linux__)
#include <features.h>
#include <malloc.h>
#include <sys/types.h>
#endif

namespace LAMMPS_NS {

enum {COMPUTES=1<<0,
      DUMPS=1<<1,
      FIXES=1<<2,
      GROUPS=1<<3,
      REGIONS=1<<4,
      CONFIG=1<<5,
      TIME=1<<6,
      MEMORY=1<<7,
      VARIABLES=1<<8,
      SYSTEM=1<<9,
      COMM=1<<10,
      COEFFS=1<<11,
      ATOM_STYLES=1<<12,
      INTEGRATE_STYLES=1<<13,
      MINIMIZE_STYLES=1<<14,
      PAIR_STYLES=1<<15,
      BOND_STYLES=1<<16,
      ANGLE_STYLES=1<<17,
      DIHEDRAL_STYLES=1<<18,
      IMPROPER_STYLES=1<<19,
      KSPACE_STYLES=1<<20,
      FIX_STYLES=1<<21,
      COMPUTE_STYLES=1<<22,
      REGION_STYLES=1<<23,
      DUMP_STYLES=1<<24,
      COMMAND_STYLES=1<<25,
      ACCELERATOR=1<<26,
      ALL=~0};

static const int STYLES = ATOM_STYLES | INTEGRATE_STYLES | MINIMIZE_STYLES
                        | PAIR_STYLES | BOND_STYLES | ANGLE_STYLES
                        | DIHEDRAL_STYLES | IMPROPER_STYLES | KSPACE_STYLES
                        | FIX_STYLES | COMPUTE_STYLES | REGION_STYLES
                        | DUMP_STYLES | COMMAND_STYLES;
}

using namespace LAMMPS_NS;

// must match enumerator in variable.h
static const char *varstyles[] = {
  "index", "loop", "world", "universe", "uloop", "string", "getenv",
  "file", "atomfile", "format", "equal", "atom", "vector", "python",
  "timer", "internal", "(unknown)"};

static const char *mapstyles[] = { "none", "array", "hash", "yes" };

static const char *commstyles[] = { "brick", "tiled" };
static const char *commlayout[] = { "uniform", "nonuniform", "irregular" };

static const char bstyles[] = "pfsm";

template<typename ValueType>
static void print_columns(FILE *fp, std::map<std::string, ValueType> *styles);

template<typename ValueType>
static bool find_style(const LAMMPS *lmp, std::map<std::string, ValueType> *styles,
                       const std::string &name, bool suffix_check);

template<typename ValueType>
static std::vector<std::string> get_style_names(std::map<std::string, ValueType> *styles);

/* ---------------------------------------------------------------------- */

void Info::command(int narg, char **arg)
{
  FILE *out=screen;
  int flags=0;

  if (comm->me != 0) return;

  // parse arguments
  int idx = 0;

  while (idx < narg) {
    if (strncmp(arg[idx],"all",3) == 0) {
      flags |= ALL;
      ++idx;
    } else if ((idx+1 < narg) && (strncmp(arg[idx],"out",3) == 0)
               && (strncmp(arg[idx+1],"screen",3) == 0)) {
      if ((out != screen) && (out != logfile)) fclose(out);
      out = screen;
      idx += 2;
    } else if ((idx+1 < narg) && (strncmp(arg[idx],"out",3) == 0)
               && (strncmp(arg[idx+1],"log",3) == 0)) {
      if ((out != screen) && (out != logfile)) fclose(out);
      out = logfile;
      idx += 2;
    } else if ((idx+2 < narg) && (strncmp(arg[idx],"out",3) == 0)
               && (strncmp(arg[idx+1],"append",3) == 0)) {
      if ((out != screen) && (out != logfile)) fclose(out);
      out = fopen(arg[idx+2],"a");
      idx += 3;
    } else if ((idx+2 < narg) && (strncmp(arg[idx],"out",3) == 0)
               && (strncmp(arg[idx+1],"overwrite",3) == 0)) {
      if ((out != screen) && (out != logfile)) fclose(out);
      out = fopen(arg[idx+2],"w");
      idx += 3;
    } else if (strncmp(arg[idx],"communication",5) == 0) {
      flags |= COMM;
      ++idx;
    } else if (strncmp(arg[idx],"computes",5) == 0) {
      flags |= COMPUTES;
      ++idx;
    } else if (strncmp(arg[idx],"dumps",5) == 0) {
      flags |= DUMPS;
      ++idx;
    } else if (strncmp(arg[idx],"fixes",5) == 0) {
      flags |= FIXES;
      ++idx;
    } else if (strncmp(arg[idx],"groups",3) == 0) {
      flags |= GROUPS;
      ++idx;
    } else if (strncmp(arg[idx],"regions",3) == 0) {
      flags |= REGIONS;
      ++idx;
    } else if (strncmp(arg[idx],"config",3) == 0) {
      flags |= CONFIG;
      ++idx;
    } else if (strncmp(arg[idx],"time",3) == 0) {
      flags |= TIME;
      ++idx;
    } else if (strncmp(arg[idx],"memory",3) == 0) {
      flags |= MEMORY;
      ++idx;
    } else if (strncmp(arg[idx],"variables",3) == 0) {
      flags |= VARIABLES;
      ++idx;
    } else if (strncmp(arg[idx],"system",3) == 0) {
      flags |= SYSTEM;
      ++idx;
    } else if (strncmp(arg[idx],"coeffs",3) == 0) {
      flags |= COEFFS;
      ++idx;
    } else if (strncmp(arg[idx],"accelerator",3) == 0) {
      flags |= ACCELERATOR;
      ++idx;
    } else if (strncmp(arg[idx],"styles",3) == 0) {
      if (idx+1 < narg) {
        ++idx;
        if (strncmp(arg[idx],"all",3) == 0) {
          flags |= STYLES;
          ++idx;
        } else if (strncmp(arg[idx],"atom",3) == 0) {
          flags |= ATOM_STYLES;
          ++idx;
        } else if (strncmp(arg[idx],"integrate",3) == 0) {
          flags |= INTEGRATE_STYLES;
          ++idx;
        } else if (strncmp(arg[idx],"minimize",3) == 0) {
          flags |= MINIMIZE_STYLES;
          ++idx;
        } else if (strncmp(arg[idx],"pair",3) == 0) {
          flags |= PAIR_STYLES;
          ++idx;
        } else if (strncmp(arg[idx],"bond",3) == 0) {
          flags |= BOND_STYLES;
          ++idx;
        } else if (strncmp(arg[idx],"angle",3) == 0) {
          flags |= ANGLE_STYLES;
          ++idx;
        } else if (strncmp(arg[idx],"dihedral",3) == 0) {
          flags |= DIHEDRAL_STYLES;
          ++idx;
        } else if (strncmp(arg[idx],"improper",3) == 0) {
          flags |= IMPROPER_STYLES;
          ++idx;
        } else if (strncmp(arg[idx],"kspace",3) == 0) {
          flags |= KSPACE_STYLES;
          ++idx;
        } else if (strncmp(arg[idx],"fix",3) == 0) {
          flags |= FIX_STYLES;
          ++idx;
        } else if (strncmp(arg[idx],"compute",4) == 0) {
          flags |= COMPUTE_STYLES;
          ++idx;
        } else if (strncmp(arg[idx],"region",3) == 0) {
          flags |= REGION_STYLES;
          ++idx;
        } else if (strncmp(arg[idx],"dump",3) == 0) {
          flags |= DUMP_STYLES;
          ++idx;
        } else if (strncmp(arg[idx],"command",4) == 0) {
          flags |= COMMAND_STYLES;
          ++idx;
        } else {
          flags |= STYLES;
        }
      } else {
        flags |= STYLES;
        ++idx;
      }
    } else {
      error->warning(FLERR,"Ignoring unknown or incorrect info command flag");
      ++idx;
    }
  }

  if (out == nullptr) return;

  fputs("\nInfo-Info-Info-Info-Info-Info-Info-Info-Info-Info-Info\n",out);
  std::time_t now = std::time(nullptr);
  fmt::print(out,"Printed on {:%a %b %d %H:%M:%S %Y}\n", fmt::localtime(now));

  if (flags & CONFIG) {
    fmt::print(out,"\nLAMMPS version: {} / {}\n", lmp->version, lmp->num_ver);

    if (LAMMPS::has_git_info())
      fmt::print(out,"Git info: {} / {} / {}\n",
                 LAMMPS::git_branch(), LAMMPS::git_descriptor(),LAMMPS::git_commit());

    fmt::print(out,"\nOS information: {}\n\n",platform::os_info());

    fmt::print(out,"sizeof(smallint): {}-bit\n"
               "sizeof(imageint): {}-bit\n"
               "sizeof(tagint):   {}-bit\n"
               "sizeof(bigint):   {}-bit\n",
               sizeof(smallint)*8, sizeof(imageint)*8,
               sizeof(tagint)*8, sizeof(bigint)*8);

    fmt::print(out,"\nCompiler: {} with {}\nC++ standard: {}\n",
               platform::compiler_info(),platform::openmp_standard(),platform::cxx_standard());

    fputs("\nActive compile time flags:\n\n",out);
    if (has_gzip_support()) fputs("-DLAMMPS_GZIP\n",out);
    if (has_png_support()) fputs("-DLAMMPS_PNG\n",out);
    if (has_jpeg_support()) fputs("-DLAMMPS_JPEG\n",out);
    if (has_ffmpeg_support()) fputs("-DLAMMPS_FFMPEG\n",out);
    if (has_fft_single_support()) fputs("-DFFT_SINGLE\n",out);
    if (has_exceptions()) fputs("-DLAMMPS_EXCEPTIONS\n",out);

#if defined(LAMMPS_BIGBIG)
    fputs("-DLAMMPS_BIGBIG\n",out);
#elif defined(LAMMPS_SMALLBIG)
    fputs("-DLAMMPS_SMALLBIG\n",out);
#else // defined(LAMMPS_SMALLSMALL)
    fputs("-DLAMMPS_SMALLSMALL\n",out);
#endif
    if (has_gzip_support()) fmt::print(out,"\n{}\n",platform::compress_info());

    int ncword, ncline = 0;
    fputs("\nInstalled packages:\n\n",out);
    for (const char **pkg = LAMMPS::installed_packages; *pkg != nullptr; ++pkg) {
      ncword = strlen(*pkg);
      if (ncline + ncword > 78) {
        ncline = 0;
        fputs("\n",out);
      }
      fmt::print(out,"{} ",*pkg);
      ncline += ncword + 1;
    }
    fputs("\n",out);
  }

  if (flags & ACCELERATOR) {
    fmt::print(out,"\nAccelerator configuration:\n\n{}",
               get_accelerator_info());
    if (Info::has_gpu_device())
      fmt::print(out,"\nAvailable GPU devices:\n{}\n",get_gpu_device_info());
  }

  if (flags & MEMORY) {
    double meminfo[3];

    get_memory_info(meminfo);

    fputs("\nMemory allocation information (MPI rank 0):\n\n",out);
    fmt::print(out,"Total dynamically allocated memory: {:.4} Mbyte\n",
               meminfo[0]);

#if defined(_WIN32)
    fmt::print(out,"Non-shared memory use: {:.4} Mbyte\n",meminfo[1]);
    fmt::print(out,"Maximum working set size: {:.4} Mbyte\n",meminfo[2]);
#else
#if defined(__linux__)
    fmt::print(out,"Current reserved memory pool size: {:.4} Mbyte\n",
               meminfo[1]);
#endif
    fmt::print(out,"Maximum resident set size: {:.4} Mbyte\n",meminfo[2]);
#endif
  }

  if (flags & COMM) {
    int major,minor;
    std::string version = platform::mpi_info(major,minor);

    fmt::print(out,"\nCommunication information:\n"
               "MPI library level: MPI v{}.{}\n"
               "MPI version: {}\n",major,minor,version);

    fmt::print(out,"Comm style = {},  Comm layout = {}\n"
               "Communicate velocities for ghost atoms = {}\n",
               commstyles[comm->style], commlayout[comm->layout],
               comm->ghost_velocity ? "yes" : "no");

    if (domain->box_exist) {
      if (comm->mode == 0)
        fmt::print(out,"Communication mode = single\n"
                   "Communication cutoff = {}\n",
                   comm->get_comm_cutoff());

      if (comm->mode == 1) {
        fputs("Communication mode = multi\n",out);
        double cut;
        for (int i=0; i < neighbor->ncollections; ++i) {
          if (comm->cutusermulti) cut = comm->cutusermulti[i];
          else cut = 0.0;
          for (int j=0; j < neighbor->ncollections; ++j) {
            cut = MAX(cut,sqrt(neighbor->cutcollectionsq[i][j]));
          }

          if (comm->cutusermulti) cut = MAX(cut,comm->cutusermulti[i]);
          fmt::print(out,"Communication cutoff for collection {} = {:.8}\n", i, cut);
        }
      }

      if (comm->mode == 2) {
        fputs("Communication mode = multi/old\n",out);
        double cut;
        for (int i=1; i <= atom->ntypes && neighbor->cuttype; ++i) {
          cut = neighbor->cuttype[i];
          if (comm->cutusermultiold) cut = MAX(cut,comm->cutusermultiold[i]);
          fmt::print(out,"Communication cutoff for type {} = {:.8}\n", i, cut);
        }
      }
    }
    fmt::print(out,"Nprocs = {},   Nthreads = {}\n",comm->nprocs,comm->nthreads);
    if (domain->box_exist)
      fmt::print(out,"Processor grid = {} x {} x {}\n",comm->procgrid[0],
                 comm->procgrid[1], comm->procgrid[2]);
  }

  if (flags & SYSTEM) {
    fputs("\nSystem information:\n",out);
    fmt::print(out,"Units         = {}\n", update->unit_style);
    fmt::print(out,"Atom style    = {}\n", atom->get_style());
    fmt::print(out,"Atom map      = {}\n", mapstyles[atom->map_style]);
    if (atom->molecular != Atom::ATOMIC) {
      const char *msg;
      msg = (atom->molecular == Atom::TEMPLATE) ? "template" : "standard";
      fmt::print(out,"Molecule type = {}\n",msg);
    }
    fmt::print(out,"Atoms     = {:12},  types = {:8d},  style = {}\n",
               atom->natoms, atom->ntypes, force->pair_style);

    if (force->pair && utils::strmatch(force->pair_style,"^hybrid")) {
      auto hybrid = dynamic_cast<PairHybrid *>(force->pair);
      fmt::print(out,"Hybrid sub-styles:");
      for (int i=0; i < hybrid->nstyles; ++i)
        fmt::print(out," {}", hybrid->keywords[i]);
      fputc('\n',out);
    }
    if (atom->molecular != Atom::ATOMIC) {
      const char *msg;
      msg = force->bond_style ? force->bond_style : "none";
      fmt::print(out,"Bonds     = {:12},  types = {:8},  style = {}\n",
                 atom->nbonds, atom->nbondtypes, msg);

      msg = force->angle_style ? force->angle_style : "none";
      fmt::print(out,"Angles    = {:12},  types = {:8},  style = {}\n",
                 atom->nangles, atom->nangletypes, msg);

      msg = force->dihedral_style ? force->dihedral_style : "none";
      fmt::print(out,"Dihedrals = {:12},  types = {:8},  style = {}\n",
                 atom->ndihedrals, atom->ndihedraltypes, msg);

      msg = force->improper_style ? force->improper_style : "none";
      fmt::print(out,"Impropers = {:12},  types = {:8},  style = {}\n",
                 atom->nimpropers, atom->nimpropertypes, msg);

      const double * const special_lj   = force->special_lj;
      const double * const special_coul = force->special_coul;

      fmt::print(out,"Special bond factors lj =    {:<8} {:<8} {:<8}\n"
                 "Special bond factors coul =  {:<8} {:<8} {:<8}\n",
                 special_lj[1],special_lj[2],special_lj[3],
                 special_coul[1],special_coul[2],special_coul[3]);
    }

    fmt::print(out,"Kspace style = {}\n",
               force->kspace ? force->kspace_style : "none");

    if (domain->box_exist) {
      fmt::print(out,"\nDimensions = {}\n",domain->dimension);
      fmt::print(out,"{} box = {:.8} x {:.8} x {:.8}\n",
                 domain->triclinic ? "Triclinic" : "Orthogonal",
                 domain->xprd, domain->yprd, domain->zprd);
      fmt::print(out,"Boundaries = {},{} {},{} {},{}\n",
                 bstyles[domain->boundary[0][0]],bstyles[domain->boundary[0][1]],
                 bstyles[domain->boundary[1][0]],bstyles[domain->boundary[1][1]],
                 bstyles[domain->boundary[2][0]],bstyles[domain->boundary[2][1]]);
      fmt::print(out,"xlo, xhi = {:.8}, {:.8}\n", domain->boxlo[0], domain->boxhi[0]);
      fmt::print(out,"ylo, yhi = {:.8}, {:.8}\n", domain->boxlo[1], domain->boxhi[1]);
      fmt::print(out,"zlo, zhi = {:.8}, {:.8}\n", domain->boxlo[2], domain->boxhi[2]);
      if (domain->triclinic)
        fmt::print(out,"Xy, xz, yz = {:.8}, {:.8}, {:.8}\n",
                   domain->xy, domain->xz, domain->yz);
    } else {
      fputs("\nBox has not yet been created\n",out);
    }
  }

  if (domain->box_exist && (flags & COEFFS)) {
    Pair *pair=force->pair;

    fputs("\nCoeff status information:\n",out);
    if (pair) {
      fputs("\nPair Coeffs:\n",out);
      for (int i=1; i <= atom->ntypes; ++i)
        for (int j=i; j <= atom->ntypes; ++j) {
          fmt::print(out,"{:6d} {:6d}:",i,j);
          if (pair->allocated && pair->setflag[i][j]) fputs(" is set\n",out);
          else fputs(" is not set\n",out);
        }
    }
    if (force->bond) {
      Bond *bond=force->bond;

      if (bond) {
        fputs("\nBond Coeffs:\n",out);
        for (int i=1; i <= atom->nbondtypes; ++i) {
          fmt::print(out,"{:6d}:",i);
          if (bond->allocated && bond->setflag[i]) fputs(" is set\n",out);
          else fputs (" is not set\n",out);
        }
      }
    }
    if (force->angle) {
      Angle *angle=force->angle;

      if (angle) {
        fputs("\nAngle Coeffs:\n",out);
        for (int i=1; i <= atom->nangletypes; ++i) {
          fmt::print(out,"{:6d}:",i);
          if (angle->allocated && angle->setflag[i]) fputs(" is set\n",out);
          else fputs (" is not set\n",out);
        }
      }
    }
    if (force->dihedral) {
      Dihedral *dihedral=force->dihedral;

      if (dihedral) {
        fputs("\nDihedral Coeffs:\n",out);
        for (int i=1; i <= atom->ndihedraltypes; ++i) {
          fmt::print(out,"{:6d}:",i);
          if (dihedral->allocated && dihedral->setflag[i]) fputs(" is set\n",out);
          else fputs (" is not set\n",out);
        }
      }
    }
    if (force->improper) {
      Improper *b=force->improper;

      if (b) {
        fputs("\nImproper Coeffs:\n",out);
        for (int i=1; i <= atom->nimpropertypes; ++i) {
          fmt::print(out,"{:6d}:",i);
          if (b->allocated && b->setflag[i]) fputs(" is set\n",out);
          else fputs (" is not set\n",out);
        }
      }
    }
  }

  if (flags & GROUPS) {
    int ngroup = group->ngroup;
    char **names = group->names;
    int *dynamic = group->dynamic;
    fputs("\nGroup information:\n",out);
    for (int i=0; i < ngroup; ++i) {
      if (names[i])
        fmt::print(out,"Group[{:2d}]:     {:16} ({})\n",
                   i, names[i], dynamic[i] ? "dynamic" : "static");
    }
  }

  if (flags & REGIONS) {
    fputs("\nRegion information:\n",out);
    int i=0;
    for (auto &reg : domain->get_region_list()) {
      fmt::print(out,"Region[{:3d}]:  {:16}  style = {:16}  side = {}\n",
                 i, std::string(reg->id)+',', std::string(reg->style)+',',
                 reg->interior ? "in" : "out");
      if (reg->bboxflag)
        fmt::print(out,"   Boundary:  lo {:.8} {:.8} {:.8}  hi {:.8} {:.8} {:.8}\n",
                   reg->extent_xlo, reg->extent_ylo,
                   reg->extent_zlo, reg->extent_xhi,
                   reg->extent_yhi, reg->extent_zhi);
      else fputs("   No Boundary\n",out);
    }
  }

  if (flags & COMPUTES) {
    int i = 0;
    char **names = group->names;
    fputs("\nCompute information:\n",out);
    for (const auto &compute : modify->get_compute_list())
      fmt::print(out,"Compute[{:3d}]:  {:16}  style = {:16}  group = {}\n", i++,
                 std::string(compute->id)+',',std::string(compute->style)+',',
                 names[compute->igroup]);
  }

  if (flags & DUMPS) {
    int ndump = output->ndump;
    Dump **dump = output->dump;
    int *nevery = output->every_dump;           \
    char **vnames = output->var_dump;
    char **names = group->names;
    fputs("\nDump information:\n",out);
    for (int i=0; i < ndump; ++i) {
      fmt::print(out,"Dump[{:3d}]:     {:16}  file = {:16}  style = {:16}  group = {:16}  ",
                 i, std::string(dump[i]->id)+',',std::string(dump[i]->filename)+',',
                 std::string(dump[i]->style)+',',std::string(names[dump[i]->igroup])+',');
      if (nevery[i]) {
        fmt::print(out,"every = {}\n", nevery[i]);
      } else {
        fmt::print(out,"every = {}\n", vnames[i]);
      }
    }
  }

  if (flags & FIXES) {
    int i = 0;
    char **names = group->names;
    fputs("\nFix information:\n",out);
    for (const auto &fix : modify->get_fix_list())
      fmt::print(out, "Fix[{:3d}]:      {:16}  style = {:16}  group = {}\n",i++,
                 std::string(fix->id)+',',std::string(fix->style)+',',names[fix->igroup]);
  }

  if (flags & VARIABLES) {
    int nvar = input->variable->nvar;
    int *style = input->variable->style;
    char **names = input->variable->names;
    char ***data = input->variable->data;
    fputs("\nVariable information:\n",out);
    for (int i=0; i < nvar; ++i) {
      int ndata = 1;
      fmt::print(out,"Variable[{:3d}]: {:16}  style = {:16}  def =",
                 i,std::string(names[i])+',',std::string(varstyles[style[i]])+',');
      if (style[i] == Variable::INTERNAL) {
        fmt::print(out,"{:.8}\n",input->variable->dvalue[i]);
        continue;
      }
      if ((style[i] != Variable::LOOP) && (style[i] != Variable::ULOOP))
        ndata = input->variable->num[i];
      for (int j=0; j < ndata; ++j)
        if (data[i][j]) fmt::print(out," {}",data[i][j]);
      fputs("\n",out);
    }
  }

  if (flags & TIME) {
    double wallclock = platform::walltime() - lmp->initclock;
    double cpuclock = platform::cputime();

    int cpuh,cpum,cpus,wallh,wallm,walls;
    cpus = fmod(cpuclock,60.0);
    cpuclock = (cpuclock - cpus) / 60.0;
    cpum = fmod(cpuclock,60.0);
    cpuh = (cpuclock - cpum) / 60.0;
    walls = fmod(wallclock,60.0);
    wallclock = (wallclock - walls) / 60.0;
    wallm = fmod(wallclock,60.0);
    wallh = (wallclock - wallm) / 60.0;
    fmt::print(out,"\nTotal time information (MPI rank 0):\n"
               "  CPU time: {:4d}:{:02d}:{:02d}\n"
               " Wall time: {:4d}:{:02d}:{:02d}\n",
               cpuh,cpum,cpus,wallh,wallm,walls);
  }

  if (flags & STYLES) {
    available_styles(out, flags);
  }

  fputs("\nInfo-Info-Info-Info-Info-Info-Info-Info-Info-Info-Info\n\n",out);

  // close output file pointer if opened locally thus forcing a hard sync.
  if ((out != screen) && (out != logfile))
    fclose(out);
}

void Info::available_styles(FILE * out, int flags)
{

  fputs("\nStyles information:\n",out);

  if (flags & ATOM_STYLES)      atom_styles(out);
  if (flags & INTEGRATE_STYLES) integrate_styles(out);
  if (flags & MINIMIZE_STYLES)  minimize_styles(out);
  if (flags & PAIR_STYLES)      pair_styles(out);
  if (flags & BOND_STYLES)      bond_styles(out);
  if (flags & ANGLE_STYLES)     angle_styles(out);
  if (flags & DIHEDRAL_STYLES)  dihedral_styles(out);
  if (flags & IMPROPER_STYLES)  improper_styles(out);
  if (flags & KSPACE_STYLES)    kspace_styles(out);
  if (flags & FIX_STYLES)       fix_styles(out);
  if (flags & COMPUTE_STYLES)   compute_styles(out);
  if (flags & REGION_STYLES)    region_styles(out);
  if (flags & DUMP_STYLES)      dump_styles(out);
  if (flags & COMMAND_STYLES)   command_styles(out);
}

void Info::atom_styles(FILE *out)
{
  fputs("\nAtom styles:\n",out);
  print_columns(out, atom->avec_map);
  fputs("\n\n\n",out);
}

void Info::integrate_styles(FILE *out)
{
  fputs("\nIntegrate styles:\n",out);
  print_columns(out, update->integrate_map);
  fputs("\n\n\n",out);
}

void Info::minimize_styles(FILE *out)
{
  fputs("\nMinimize styles:\n",out);
  print_columns(out, update->minimize_map);
  fputs("\n\n\n",out);
}

void Info::pair_styles(FILE *out)
{
  fputs("\nPair styles:\n",out);
  print_columns(out, force->pair_map);
  fputs("\n\n\n",out);
}

void Info::bond_styles(FILE *out)
{
  fputs("\nBond styles:\n",out);
  print_columns(out, force->bond_map);
  fputs("\n\n\n",out);
}

void Info::angle_styles(FILE *out)
{
  fputs("\nAngle styles:\n",out);
  print_columns(out, force->angle_map);
  fputs("\n\n\n",out);
}

void Info::dihedral_styles(FILE *out)
{
  fputs("\nDihedral styles:\n",out);
  print_columns(out, force->dihedral_map);
  fputs("\n\n\n",out);
}

void Info::improper_styles(FILE *out)
{
  fputs("\nImproper styles:\n",out);
  print_columns(out, force->improper_map);
  fputs("\n\n\n",out);
}

void Info::kspace_styles(FILE *out)
{
  fputs("\nKSpace styles:\n",out);
  print_columns(out, force->kspace_map);
  fputs("\n\n\n",out);
}

void Info::fix_styles(FILE *out)
{
  fputs("\nFix styles:\n",out);
  print_columns(out, modify->fix_map);
  fputs("\n\n\n",out);
}

void Info::compute_styles(FILE *out)
{
  fputs("\nCompute styles:\n",out);
  print_columns(out, modify->compute_map);
  fputs("\n\n\n",out);
}

void Info::region_styles(FILE *out)
{
  fputs("\nRegion styles:\n",out);
  print_columns(out, domain->region_map);
  fputs("\n\n\n",out);
}

void Info::dump_styles(FILE *out)
{
  fputs("\nDump styles:\n",out);
  print_columns(out,output->dump_map);
  fputs("\n\n\n",out);
}

void Info::command_styles(FILE *out)
{
  fputs("\nCommand styles (add-on input script commands):\n",out);
  print_columns(out, input->command_map);
  fputs("\n\n\n",out);
}


/* ---------------------------------------------------------------------- */

// the is_active() function returns true if the selected style or name
// in the selected category is currently in use.

bool Info::is_active(const char *category, const char *name)
{
  if ((category == nullptr) || (name == nullptr)) return false;
  const char *style = "none";

  if (strcmp(category,"package") == 0) {
    if (strcmp(name,"gpu") == 0) {
      return modify->get_fix_by_id("package_gpu") != nullptr;
    } else if (strcmp(name,"intel") == 0) {
      return modify->get_fix_by_id("package_intel") != nullptr;
    } else if (strcmp(name,"kokkos") == 0) {
      return lmp->kokkos && lmp->kokkos->kokkos_exists;
    } else if (strcmp(name,"omp") == 0) {
      return modify->get_fix_by_id("package_omp") != nullptr;
    } else error->all(FLERR,"Unknown name for info package category: {}", name);

  } else if (strcmp(category,"newton") == 0) {
    if (strcmp(name,"pair") == 0) return (force->newton_pair != 0);
    else if (strcmp(name,"bond") == 0) return (force->newton_bond != 0);
    else if (strcmp(name,"any") == 0) return (force->newton != 0);
    else error->all(FLERR,"Unknown name for info newton category: {}", name);

  } else if (strcmp(category,"pair") == 0) {
    if (force->pair == nullptr) return false;
    if (strcmp(name,"single") == 0) return (force->pair->single_enable != 0);
    else if (strcmp(name,"respa") == 0) return (force->pair->respa_enable != 0);
    else if (strcmp(name,"manybody") == 0) return (force->pair->manybody_flag != 0);
    else if (strcmp(name,"tail") == 0) return (force->pair->tail_flag != 0);
    else if (strcmp(name,"shift") == 0) return (force->pair->offset_flag != 0);
    else error->all(FLERR,"Unknown name for info pair category: {}", name);

  } else if (strcmp(category,"comm_style") == 0) {
    style = commstyles[comm->style];
  } else if (strcmp(category,"min_style") == 0) {
    style = update->minimize_style;
  } else if (strcmp(category,"run_style") == 0) {
    style = update->integrate_style;
  } else if (strcmp(category,"atom_style") == 0) {
    style = atom->atom_style;
  } else if (strcmp(category,"pair_style") == 0) {
    style = force->pair_style;
  } else if (strcmp(category,"bond_style") == 0) {
    style = force->bond_style;
  } else if (strcmp(category,"angle_style") == 0) {
    style = force->angle_style;
  } else if (strcmp(category,"dihedral_style") == 0) {
    style = force->dihedral_style;
  } else if (strcmp(category,"improper_style") == 0) {
    style = force->improper_style;
  } else if (strcmp(category,"kspace_style") == 0) {
    style = force->kspace_style;
  } else error->all(FLERR,"Unknown category for info is_active(): {}", category);

  int match = 0;
  if (strcmp(style,name) == 0) match = 1;

  if (!match && lmp->suffix_enable) {
    if (lmp->suffix) {
      std::string name_w_suffix = name + std::string("/") + lmp->suffix;
      if (name_w_suffix == style) match = 1;
    }
    if (!match && lmp->suffix2) {
      std::string name_w_suffix = name + std::string("/") + lmp->suffix2;
      if (name_w_suffix == style) match = 1;
    }
  }
  return match != 0;
}

/* ---------------------------------------------------------------------- */

// the is_available() function returns true if the selected style
// or name in the selected category is available for use (but need
// not be currently active).

bool Info::is_available(const char *category, const char *name)
{
  if ((category == nullptr) || (name == nullptr)) return false;

  if (has_style(category, name)) {
    return true;
  } else if (strcmp(category,"feature") == 0) {
    if (strcmp(name,"gzip") == 0) {
      return has_gzip_support();
    } else if (strcmp(name,"png") == 0) {
      return has_png_support();
    } else if (strcmp(name,"jpeg") == 0) {
      return has_jpeg_support();
    } else if (strcmp(name,"ffmpeg") == 0) {
      return has_ffmpeg_support();
    } else if (strcmp(name,"fft_single") == 0) {
      return has_fft_single_support();
    } else if (strcmp(name,"exceptions") == 0) {
      return has_exceptions();
    }
  } else error->all(FLERR,"Unknown category for info is_available(): {}", category);

  return false;
}

/* ---------------------------------------------------------------------- */

// the is_defined() function returns true if a particular ID of the
// selected category (e.g. fix ID, group ID, region ID etc.) has been
// defined and thus can be accessed. It does *NOT* check whether a
// particular ID has a particular style.

bool Info::is_defined(const char *category, const char *name)
{
  if ((category == nullptr) || (name == nullptr)) return false;

  if (strcmp(category,"compute") == 0) {
    if (modify->get_compute_by_id(name)) return true;
  } else if (strcmp(category,"dump") == 0) {
    if (output->get_dump_by_id(name)) return true;
  } else if (strcmp(category,"fix") == 0) {
    if (modify->get_fix_by_id(name)) return true;
  } else if (strcmp(category,"group") == 0) {
    if (group->find(name) >= 0) return true;
  } else if (strcmp(category,"region") == 0) {
    if (domain->get_region_by_id(name)) return true;
  } else if (strcmp(category,"variable") == 0) {
    if (input->variable->find(name) >= 0) return true;
  } else error->all(FLERR,"Unknown category for info is_defined(): {}", category);

  return false;
}

bool Info::has_style(const std::string &category, const std::string &name)
{
  if (category == "atom") {
    return find_style(lmp, atom->avec_map, name, false);
  } else if (category == "integrate") {
    return find_style(lmp, update->integrate_map, name, true);
  } else if (category == "minimize") {
    return find_style(lmp, update->minimize_map, name, true);
  } else if (category == "pair") {
    return find_style(lmp, force->pair_map, name, true);
  } else if (category == "bond") {
    return find_style(lmp, force->bond_map, name, true);
  } else if (category == "angle") {
    return find_style(lmp, force->angle_map, name, true);
  } else if (category == "dihedral") {
    return find_style(lmp, force->dihedral_map, name, true);
  } else if (category == "improper") {
    return find_style(lmp, force->improper_map, name, true);
  } else if (category == "kspace") {
    return find_style(lmp, force->kspace_map, name, true);
  } else if (category == "fix") {
    return find_style(lmp, modify->fix_map, name, true);
  } else if (category == "compute") {
    return find_style(lmp, modify->compute_map, name, true);
  } else if (category == "region") {
    return find_style(lmp, domain->region_map, name, false);
  } else if (category == "dump") {
    return find_style(lmp, output->dump_map, name, false);
  } else if (category == "command") {
    return find_style(lmp, input->command_map, name, false);
  }
  return false;
}

std::vector<std::string> Info::get_available_styles(const std::string &category)
{
  if (category == "atom") {
    return get_style_names(atom->avec_map);
  } else if (category == "integrate") {
    return get_style_names(update->integrate_map);
  } else if (category == "minimize") {
    return get_style_names(update->minimize_map);
  } else if (category == "pair") {
    return get_style_names(force->pair_map);
  } else if (category == "bond") {
    return get_style_names(force->bond_map);
  } else if (category == "angle") {
    return get_style_names(force->angle_map);
  } else if (category == "dihedral") {
    return get_style_names(force->dihedral_map);
  } else if (category == "improper") {
    return get_style_names(force->improper_map);
  } else if (category == "kspace") {
    return get_style_names(force->kspace_map);
  } else if (category == "fix") {
    return get_style_names(modify->fix_map);
  } else if (category == "compute") {
    return get_style_names(modify->compute_map);
  } else if (category == "region") {
    return get_style_names(domain->region_map);
  } else if (category == "dump") {
    return get_style_names(output->dump_map);
  } else if (category == "command") {
    return get_style_names(input->command_map);
  }
  return {};
}

template<typename ValueType>
static std::vector<std::string> get_style_names(std::map<std::string, ValueType> *styles)
{
  std::vector<std::string> names;

  names.reserve(styles->size());
  for (auto const &kv : *styles) {
    // skip "secret" styles
    if (isupper(kv.first[0])) continue;
    names.push_back(kv.first);
  }

  return names;
}

template<typename ValueType>
static bool find_style(const LAMMPS *lmp, std::map<std::string, ValueType> *styles,
                       const std::string &name, bool suffix_check)
{
  if (styles->find(name) != styles->end()) {
    return true;
  }

  if (suffix_check && lmp->suffix_enable) {
    if (lmp->suffix) {
      std::string name_w_suffix = name + "/" + lmp->suffix;
      if (find_style(lmp, styles, name_w_suffix, false)) {
        return true;
      }
    }
    if (lmp->suffix2) {
      std::string name_w_suffix = name + "/" + lmp->suffix2;
      if (find_style(lmp, styles, name_w_suffix, false)) {
        return true;
      }
    }
  }
  return false;
}

template<typename ValueType>
static void print_columns(FILE *fp, std::map<std::string, ValueType> *styles)
{
  if (styles->empty()) {
    fprintf(fp, "\nNone");
    return;
  }

  // std::map keys are already sorted
  int pos = 80;
  for (auto it = styles->begin(); it != styles->end(); ++it) {
    const std::string &style_name = it->first;

    // skip "internal" styles
    if (isupper(style_name[0]) || utils::strmatch(style_name,"/kk/host$")
        || utils::strmatch(style_name,"/kk/device$")) continue;

    int len = style_name.length();
    if (pos + len > 80) {
      fprintf(fp,"\n");
      pos = 0;
    }

    if (len < 16) {
      fprintf(fp,"%-16s", style_name.c_str());
      pos += 16;
    } else if (len < 32) {
      fprintf(fp,"%-32s", style_name.c_str());
      pos += 32;
    } else if (len < 48) {
      fprintf(fp,"%-48s", style_name.c_str());
      pos += 48;
    } else if (len < 64) {
      fprintf(fp,"%-64s", style_name.c_str());
      pos += 64;
    } else {
      fprintf(fp,"%-80s", style_name.c_str());
      pos += 80;
    }
  }
}

bool Info::has_gzip_support() {
#ifdef LAMMPS_GZIP
  return true;
#else
  return false;
#endif
}

bool Info::has_png_support() {
#ifdef LAMMPS_PNG
  return true;
#else
  return false;
#endif
}

bool Info::has_jpeg_support() {
#ifdef LAMMPS_JPEG
  return true;
#else
  return false;
#endif
}

bool Info::has_ffmpeg_support() {
#ifdef LAMMPS_FFMPEG
  return true;
#else
  return false;
#endif
}

bool Info::has_fft_single_support() {
#ifdef FFT_SINGLE
  return true;
#else
  return false;
#endif
}

bool Info::has_exceptions() {
#ifdef LAMMPS_EXCEPTIONS
  return true;
#else
  return false;
#endif
}

bool Info::has_package(const std::string &package_name) {
  for (int i = 0; LAMMPS::installed_packages[i] != nullptr; ++i) {
    if (package_name == LAMMPS::installed_packages[i]) {
      return true;
    }
  }
  return false;
}

#if defined(LMP_GPU)
extern bool lmp_gpu_config(const std::string &, const std::string &);
extern bool lmp_has_compatible_gpu_device();
extern std::string lmp_gpu_device_info();

// we will only report compatible GPUs, i.e. when a GPU device is
// available *and* supports the required floating point precision
bool Info::has_gpu_device()
{
  return lmp_has_compatible_gpu_device();
}

std::string Info::get_gpu_device_info()
{
  return lmp_gpu_device_info();
}
#else
bool Info::has_gpu_device()
{
  return false;
}
std::string Info::get_gpu_device_info()
{
  return "";
}
#endif

#if defined(LMP_KOKKOS)
#include "Kokkos_Macros.hpp"
#endif

bool Info::has_accelerator_feature(const std::string &package,
                                   const std::string &category,
                                   const std::string &setting)
{
#if defined(LMP_KOKKOS)
  if (package == "KOKKOS") {
    if (category == "precision") {
      return setting == "double";
    }
    if (category == "api") {
#if defined(KOKKOS_ENABLE_OPENMP)
      if (setting == "openmp") return true;
#endif
#if defined(KOKKOS_ENABLE_SERIAL)
      if (setting == "serial") return true;
#endif
#if defined(KOKKOS_ENABLE_THREADS)
      if (setting == "pthreads") return true;
#endif
#if defined(KOKKOS_ENABLE_CUDA)
      if (setting == "cuda") return true;
#endif
#if defined(KOKKOS_ENABLE_HIP)
      if (setting == "hip") return true;
#endif
#if defined(KOKKOS_ENABLE_SYCL)
      if (setting == "sycl") return true;
#endif
      return false;
    }
  }
#endif
#if defined(LMP_GPU)
  if (package == "GPU") {
    return lmp_gpu_config(category,setting);
  }
#endif
#if defined(LMP_OPENMP)
  if (package == "OPENMP") {
    if (category == "precision") {
      return setting == "double";
    }
    if (category == "api") {
#if defined(_OPENMP)
      return setting == "openmp";
#else
      return setting == "serial";
#endif
    }
  }
#endif
#if defined(LMP_INTEL)
  if (package == "INTEL") {
    if (category == "precision") {
      if (setting == "double") return true;
      else if (setting == "mixed") return true;
      else if (setting == "single") return true;
      else return false;
    }
    if (category == "api") {
#if defined(LMP_INTEL_OFFLOAD)
      if (setting == "phi") return true;
#elif defined(_OPENMP)
      if (setting == "openmp") return true;
#else
      if (setting == "serial") return true;
#endif
      return false;
    }
  }
#endif
  return false;
}

std::string Info::get_accelerator_info(const std::string &package)
{
  std::string mesg;
  if ((package.empty() || (package == "GPU")) && has_package("GPU")) {
    mesg += "GPU package API:";
    if (has_accelerator_feature("GPU","api","cuda"))   mesg += " CUDA";
    if (has_accelerator_feature("GPU","api","hip"))    mesg += " HIP";
    if (has_accelerator_feature("GPU","api","opencl")) mesg += " OpenCL";
    mesg +=  "\nGPU package precision:";
    if (has_accelerator_feature("GPU","precision","single")) mesg += " single";
    if (has_accelerator_feature("GPU","precision","mixed"))  mesg += " mixed";
    if (has_accelerator_feature("GPU","precision","double")) mesg += " double";
    mesg += "\n";
  }
  if ((package.empty() || (package == "KOKKOS")) && has_package("KOKKOS")) {
    mesg += "KOKKOS package API:";
    if (has_accelerator_feature("KOKKOS","api","cuda"))     mesg += " CUDA";
    if (has_accelerator_feature("KOKKOS","api","hip"))      mesg += " HIP";
    if (has_accelerator_feature("KOKKOS","api","sycl"))     mesg += " SYCL";
    if (has_accelerator_feature("KOKKOS","api","openmp"))   mesg += " OpenMP";
    if (has_accelerator_feature("KOKKOS","api","serial"))   mesg += " Serial";
    if (has_accelerator_feature("KOKKOS","api","pthreads")) mesg += " Pthreads";
    mesg +=  "\nKOKKOS package precision:";
    if (has_accelerator_feature("KOKKOS","precision","single")) mesg += " single";
    if (has_accelerator_feature("KOKKOS","precision","mixed"))  mesg += " mixed";
    if (has_accelerator_feature("KOKKOS","precision","double")) mesg += " double";
    mesg += "\n";
  }
  if ((package.empty() || (package == "OPENMP")) && has_package("OPENMP")) {
    mesg += "OPENMP package API:";
    if (has_accelerator_feature("OPENMP","api","openmp"))   mesg += " OpenMP";
    if (has_accelerator_feature("OPENMP","api","serial"))   mesg += " Serial";
    mesg +=  "\nOPENMP package precision:";
    if (has_accelerator_feature("OPENMP","precision","single")) mesg += " single";
    if (has_accelerator_feature("OPENMP","precision","mixed"))  mesg += " mixed";
    if (has_accelerator_feature("OPENMP","precision","double")) mesg += " double";
    mesg += "\n";
  }
  if ((package.empty() || (package == "INTEL")) && has_package("INTEL")) {
    mesg += "INTEL package API:";
    if (has_accelerator_feature("INTEL","api","phi"))      mesg += " Phi";
    if (has_accelerator_feature("INTEL","api","openmp"))   mesg += " OpenMP";
    mesg +=  "\nINTEL package precision:";
    if (has_accelerator_feature("INTEL","precision","single")) mesg += " single";
    if (has_accelerator_feature("INTEL","precision","mixed"))  mesg += " mixed";
    if (has_accelerator_feature("INTEL","precision","double")) mesg += " double";
    mesg += "\n";
  }
  return mesg;
}

/* ---------------------------------------------------------------------- */

void Info::get_memory_info(double *meminfo)
{
  double bytes = 0;
  bytes += atom->memory_usage();
  bytes += neighbor->memory_usage();
  bytes += comm->memory_usage();
  bytes += update->memory_usage();
  bytes += force->memory_usage();
  bytes += modify->memory_usage();
  for (int i = 0; i < output->ndump; i++)
    bytes += output->dump[i]->memory_usage();
  meminfo[0] = bytes/1024.0/1024.0;
  meminfo[1] = 0;
  meminfo[2] = 0;

#if defined(_WIN32)
  HANDLE phandle = GetCurrentProcess();
  PROCESS_MEMORY_COUNTERS_EX pmc;
  GetProcessMemoryInfo(phandle,(PROCESS_MEMORY_COUNTERS *)&pmc,sizeof(pmc));
  meminfo[1] = (double)pmc.PrivateUsage/1048576.0;
  meminfo[2] = (double)pmc.PeakWorkingSetSize/1048576.0;
#else
#if defined(__linux__)
// __GLIBC_MINOR__ is only defined on real glibc (i.e. not on musl)
#if defined(__GLIBC_MINOR__)
// newer glibc versions have mallinfo2
#if defined(__GLIBC__) && __GLIBC_PREREQ(2, 33)
  struct mallinfo2 mi;
  mi = mallinfo2();
#else
  struct mallinfo mi;
  mi = mallinfo();
#endif
  meminfo[1] = (double)mi.uordblks/1048576.0+(double)mi.hblkhd/1048576.0;
#endif
// not glibc => may not have mallinfo/mallinfo2
#else
  meminfo[1] = 0.0;
#endif
  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) == 0)
    meminfo[2] = (double)ru.ru_maxrss/1024.0;
#endif
}

/* ---------------------------------------------------------------------- */

char **Info::get_variable_names(int &num) {
  num = input->variable->nvar;
  return input->variable->names;
}
