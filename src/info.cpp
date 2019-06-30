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

/* ----------------------------------------------------------------------
   Contributing authors:  Axel Kohlmeyer (Temple U),
                          Richard Berger (Temple U)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <cstring>
#include <cctype>
#include <ctime>
#include <map>
#include <string>
#include "info.h"
#include "accelerator_kokkos.h"
#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "dump.h"
#include "fix.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "group.h"
#include "input.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "region.h"
#include "universe.h"
#include "variable.h"
#include "update.h"
#include "error.h"
#include "utils.h"

#ifdef _WIN32
#define PSAPI_VERSION 1
#include <windows.h>
#include <stdint.h> // <cstdint> requires C++-11
#include <psapi.h>
#else
#include <sys/resource.h>
#include <sys/utsname.h>
#endif

#if defined __linux
#include <malloc.h>
#endif

namespace LAMMPS_NS {
// same as in variable.cpp
enum {INDEX,LOOP,WORLD,UNIVERSE,ULOOP,STRING,GETENV,
      SCALARFILE,ATOMFILE,FORMAT,EQUAL,ATOM,VECTOR,PYTHON,INTERNAL};

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
      ALL=~0};

static const int STYLES = ATOM_STYLES | INTEGRATE_STYLES | MINIMIZE_STYLES
                        | PAIR_STYLES | BOND_STYLES | ANGLE_STYLES
                        | DIHEDRAL_STYLES | IMPROPER_STYLES | KSPACE_STYLES
                        | FIX_STYLES | COMPUTE_STYLES | REGION_STYLES
                        | DUMP_STYLES | COMMAND_STYLES;
}

static const char *varstyles[] = {
  "index", "loop", "world", "universe", "uloop", "string", "getenv",
  "file", "atomfile", "format", "equal", "atom", "vector", "python", "internal", "(unknown)"};

static const char *mapstyles[] = { "none", "array", "hash" };

static const char *commstyles[] = { "brick", "tiled" };
static const char *commlayout[] = { "uniform", "nonuniform", "irregular" };

static const char bstyles[] = "pfsm";

using namespace LAMMPS_NS;
using namespace std;

template<typename ValueType>
static void print_columns(FILE* fp, map<string, ValueType> * styles);

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

  if (out == NULL) return;

  fputs("\nInfo-Info-Info-Info-Info-Info-Info-Info-Info-Info-Info\n",out);
  time_t now = time(NULL);
  fprintf(out,"Printed on %s\n",ctime(&now));

  if (flags & CONFIG) {
    if (lmp->has_git_info) {
      fprintf(out,"\nLAMMPS version: %s / %s\nGit info: %s / %s / %s\n\n",
              universe->version, universe->num_ver,lmp->git_branch,
              lmp->git_descriptor,lmp->git_commit);
    } else {
      fprintf(out,"\nLAMMPS version: %s / %s\n\n",
              universe->version, universe->num_ver);
    }
    const char *infobuf = get_os_info();
    fprintf(out,"OS information: %s\n\n",infobuf);
    delete[] infobuf;

    fprintf(out,"sizeof(smallint): %3d-bit\n",(int)sizeof(smallint)*8);
    fprintf(out,"sizeof(imageint): %3d-bit\n",(int)sizeof(imageint)*8);
    fprintf(out,"sizeof(tagint):   %3d-bit\n",(int)sizeof(tagint)*8);
    fprintf(out,"sizeof(bigint):   %3d-bit\n",(int)sizeof(bigint)*8);

    infobuf = get_compiler_info();
    fprintf(out,"\nCompiler: %s with %s\n",infobuf,get_openmp_info());
    delete[] infobuf;

    fputs("\nActive compile time flags:\n\n",out);
    if (has_gzip_support()) fputs("-DLAMMPS_GZIP\n",out);
    if (has_png_support()) fputs("-DLAMMPS_PNG\n",out);
    if (has_jpeg_support()) fputs("-DLAMMPS_JPEG\n",out);
    if (has_ffmpeg_support()) fputs("-DLAMMPS_FFMPEG\n",out);
    if (has_exceptions()) fputs("-DLAMMPS_EXCEPTIONS\n",out);

#if defined(LAMMPS_BIGBIG)
    fputs("-DLAMMPS_BIGBIG\n",out);
#elif defined(LAMMPS_SMALLBIG)
    fputs("-DLAMMPS_SMALLBIG\n",out);
#else // defined(LAMMPS_SMALLSMALL)
    fputs("-DLAMMPS_SMALLSMALL\n",out);
#endif

    const char *pkg;
    int ncword, ncline = 0;

    fputs("\nInstalled packages:\n\n",out);
    for (int i = 0; NULL != (pkg = lmp->installed_packages[i]); ++i) {
      ncword = strlen(pkg);
      if (ncline + ncword > 78) {
        ncline = 0;
        fputs("\n",out);
      }
      fprintf(out,"%s ",pkg);
      ncline += ncword + 1;
    }
    fputs("\n",out);
  }

  if (flags & MEMORY) {

    fprintf(out,"\nMemory allocation information (MPI rank 0):\n\n");

    bigint bytes = 0;
    bytes += atom->memory_usage();
    bytes += neighbor->memory_usage();
    bytes += comm->memory_usage();
    bytes += update->memory_usage();
    bytes += force->memory_usage();
    bytes += modify->memory_usage();
    for (int i = 0; i < output->ndump; i++)
      bytes += output->dump[i]->memory_usage();
    double mbytes = bytes/1024.0/1024.0;
    fprintf(out,"Total dynamically allocated memory: %.4g Mbyte\n",mbytes);

#if defined(_WIN32)
    HANDLE phandle = GetCurrentProcess();
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(phandle,(PROCESS_MEMORY_COUNTERS *)&pmc,sizeof(pmc));
    fprintf(out,"Non-shared memory use: %.4g Mbyte\n",
            (double)pmc.PrivateUsage/1048576.0);
    fprintf(out,"Maximum working set size: %.4g Mbyte\n",
            (double)pmc.PeakWorkingSetSize/1048576.0);
#else
#if defined(__linux)
    struct mallinfo mi;
    mi = mallinfo();
    fprintf(out,"Current reserved memory pool size: %.4g Mbyte\n",
            (double)mi.uordblks/1048576.0+(double)mi.hblkhd/1048576.0);
#endif
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) == 0) {
      fprintf(out,"Maximum resident set size: %.4g Mbyte\n",
              (double)ru.ru_maxrss/1024.0);
    }
#endif
  }

  if (flags & COMM) {
    int major,minor;
    MPI_Get_version(&major,&minor);

    fprintf(out,"\nCommunication information:\n");
    fprintf(out,"MPI library level: MPI v%d.%d\n",major,minor);
    fprintf(out,"Comm style = %s,  Comm layout = %s\n",
            commstyles[comm->style], commlayout[comm->layout]);
    fprintf(out,"Communicate velocities for ghost atoms = %s\n",
            comm->ghost_velocity ? "yes" : "no");

    if (comm->mode == 0) {
      fprintf(out,"Communication mode = single\n");
      fprintf(out,"Communication cutoff = %g\n",
              MAX(comm->cutghostuser,neighbor->cutneighmax));
    }

    if (comm->mode == 1) {
      fprintf(out,"Communication mode = multi\n");
      double cut;
      for (int i=1; i <= atom->ntypes && neighbor->cuttype; ++i) {
        cut = neighbor->cuttype[i];
        if (comm->cutusermulti) cut = MAX(cut,comm->cutusermulti[i]);
        fprintf(out,"Communication cutoff for type %d = %g\n", i, cut);
      }
    }
    fprintf(out,"Nprocs = %d,   Nthreads = %d\n",
            comm->nprocs, comm->nthreads);
    if (domain->box_exist)
      fprintf(out,"Processor grid = %d x %d x %d\n",comm->procgrid[0],
            comm->procgrid[1], comm->procgrid[2]);
  }

  if (flags & SYSTEM) {
    fprintf(out,"\nSystem information:\n");
    fprintf(out,"Units      = %s\n",update->unit_style);
    fprintf(out,"Atom style = %s\n", atom->atom_style);
    fprintf(out,"Atom map   = %s\n", mapstyles[atom->map_style]);
    if (atom->molecular > 0) {
      const char *msg;
      msg = (atom->molecular == 2) ? "template" : "standard";
      fprintf(out,"Molecule type = %s\n",msg);
    }
    fprintf(out,"Atoms = " BIGINT_FORMAT ",  types = %d,  style = %s\n",
            atom->natoms, atom->ntypes, force->pair_style);

    if (force->pair && utils::strmatch(force->pair_style,"^hybrid")) {
      PairHybrid *hybrid = (PairHybrid *)force->pair;
      fprintf(out,"Hybrid sub-styles:");
      for (int i=0; i < hybrid->nstyles; ++i)
        fprintf(out," %s", hybrid->keywords[i]);
      fputc('\n',out);
    }
    if (atom->molecular > 0) {
      const char *msg;
      msg = force->bond_style ? force->bond_style : "none";
      fprintf(out,"Bonds = " BIGINT_FORMAT ",  types = %d,  style = %s\n",
              atom->nbonds, atom->nbondtypes, msg);

      msg = force->angle_style ? force->angle_style : "none";
      fprintf(out,"Angles = " BIGINT_FORMAT ",  types = %d,  style = %s\n",
              atom->nangles, atom->nangletypes, msg);

      msg = force->dihedral_style ? force->dihedral_style : "none";
      fprintf(out,"Dihedrals = " BIGINT_FORMAT ",  types = %d,  style = %s\n",
              atom->ndihedrals, atom->ndihedraltypes, msg);

      msg = force->improper_style ? force->improper_style : "none";
      fprintf(out,"Impropers = " BIGINT_FORMAT ",  types = %d,  style = %s\n",
              atom->nimpropers, atom->nimpropertypes, msg);

      const double * const special_lj   = force->special_lj;
      const double * const special_coul = force->special_coul;

      fprintf(out,"Special bond factors lj =   %-10g %-10g %-10g\n"
              "Special bond factors coul = %-10g %-10g %-10g\n",
              special_lj[1],special_lj[2],special_lj[3],
              special_coul[1],special_coul[2],special_coul[3]);
    }

    fprintf(out,"Kspace style = %s\n",
            force->kspace ? force->kspace_style : "none");

    if (domain->box_exist) {
      fprintf(out,"\nDimensions = %d\n",domain->dimension);
      fprintf(out,"%s box = %g x %g x %g\n",
              domain->triclinic ? "Triclinic" : "Orthogonal",
              domain->xprd, domain->yprd, domain->zprd);
      fprintf(out,"Boundaries = %c,%c %c,%c %c,%c\n",
              bstyles[domain->boundary[0][0]],bstyles[domain->boundary[0][1]],
              bstyles[domain->boundary[1][0]],bstyles[domain->boundary[1][1]],
              bstyles[domain->boundary[2][0]],bstyles[domain->boundary[2][1]]);
      fprintf(out,"xlo, xhi = %g, %g\n", domain->boxlo[0], domain->boxhi[0]);
      fprintf(out,"ylo, yhi = %g, %g\n", domain->boxlo[1], domain->boxhi[1]);
      fprintf(out,"zlo, zhi = %g, %g\n", domain->boxlo[2], domain->boxhi[2]);
      if (domain->triclinic)
          fprintf(out,"Xy, xz, yz = %g, %g, %g\n",
                  domain->xy, domain->xz, domain->yz);
    } else {
      fputs("\nBox has not yet been created\n",out);
    }
  }

  if (domain->box_exist && (flags & COEFFS)) {
    Pair *pair=force->pair;

    fprintf(out,"\nCoeff information:\n");
    if (pair) {
      fprintf(out,"Pair Coeffs:\n");
      for (int i=1; i <= atom->ntypes; ++i)
        for (int j=i; j <= atom->ntypes; ++j) {
          fprintf(out,"%3d %3d :",i,j);
          if (pair->allocated && pair->setflag[i][j]) fputs(" is set\n",out);
          else fputs (" is not set\n",out);
        }
    }
    if (force->bond) {
      Bond *bond=force->bond;

      if (bond) {
        fprintf(out,"Bond Coeffs:\n");
        for (int i=1; i <= atom->nbondtypes; ++i) {
          fprintf(out,"%3d :",i);
          if (bond->allocated && bond->setflag[i]) fputs(" is set\n",out);
          else fputs (" is not set\n",out);
        }
      }
    }
    if (force->angle) {
      Angle *angle=force->angle;

      if (angle) {
        fprintf(out,"Angle Coeffs:\n");
        for (int i=1; i <= atom->nangletypes; ++i) {
          fprintf(out,"%3d :",i);
          if (angle->allocated && angle->setflag[i]) fputs(" is set\n",out);
          else fputs (" is not set\n",out);
        }
      }
    }
    if (force->dihedral) {
      Dihedral *dihedral=force->dihedral;

      if (dihedral) {
        fprintf(out,"Dihedral Coeffs:\n");
        for (int i=1; i <= atom->ndihedraltypes; ++i) {
          fprintf(out,"%3d :",i);
          if (dihedral->allocated && dihedral->setflag[i]) fputs(" is set\n",out);
          else fputs (" is not set\n",out);
        }
      }
    }
    if (force->improper) {
      Improper *b=force->improper;

      if (b) {
        fprintf(out,"Improper Coeffs:\n");
        for (int i=1; i <= atom->nimpropertypes; ++i) {
          fprintf(out,"%3d :",i);
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
    fprintf(out,"\nGroup information:\n");
    for (int i=0; i < ngroup; ++i) {
      if (names[i])
        fprintf(out,"Group[%2d]: %s (%s)\n",
                i, names[i], dynamic[i] ? "dynamic" : "static");
    }
  }

  if (flags & REGIONS) {
    int nreg = domain->nregion;
    Region **regs = domain->regions;
    fprintf(out,"\nRegion information:\n");
    for (int i=0; i < nreg; ++i) {
      fprintf(out,"Region[%3d]: %s,  style = %s,  side = %s\n",
              i, regs[i]->id, regs[i]->style,
              regs[i]->interior ? "in" : "out");
      if (regs[i]->bboxflag)
        fprintf(out,"     Boundary: lo %g %g %g  hi %g %g %g\n",
                regs[i]->extent_xlo, regs[i]->extent_ylo,
                regs[i]->extent_zlo, regs[i]->extent_xhi,
                regs[i]->extent_yhi, regs[i]->extent_zhi);
      else fprintf(out,"     No Boundary\n");
    }
  }

  if (flags & COMPUTES) {
    int ncompute = modify->ncompute;
    Compute **compute = modify->compute;
    char **names = group->names;
    fprintf(out,"\nCompute information:\n");
    for (int i=0; i < ncompute; ++i) {
      fprintf(out,"Compute[%3d]: %s,  style = %s,  group = %s\n",
              i, compute[i]->id, compute[i]->style,
              names[compute[i]->igroup]);
    }
  }

  if (flags & DUMPS) {
    int ndump = output->ndump;
    Dump **dump = output->dump;
    int *nevery = output->every_dump;           \
    char **vnames = output->var_dump;
    char **names = group->names;
    fprintf(out,"\nDump information:\n");
    for (int i=0; i < ndump; ++i) {
      fprintf(out,"Dump[%3d]: %s,  file = %s,  style = %s,  group = %s,  ",
              i, dump[i]->id, dump[i]->filename,
              dump[i]->style, names[dump[i]->igroup]);
      if (nevery[i]) {
        fprintf(out,"every = %d\n", nevery[i]);
      } else {
        fprintf(out,"every = %s\n", vnames[i]);
      }
    }
  }

  if (flags & FIXES) {
    int nfix = modify->nfix;
    Fix **fix = modify->fix;
    char **names = group->names;
    fprintf(out,"\nFix information:\n");
    for (int i=0; i < nfix; ++i) {
      fprintf(out,"Fix[%3d]: %s,  style = %s,  group = %s\n",
              i, fix[i]->id, fix[i]->style, names[fix[i]->igroup]);
    }
  }

  if (flags & VARIABLES) {
    int nvar = input->variable->nvar;
    int *style = input->variable->style;
    char **names = input->variable->names;
    char ***data = input->variable->data;
    fprintf(out,"\nVariable information:\n");
    for (int i=0; i < nvar; ++i) {
      int ndata = 1;
      fprintf(out,"Variable[%3d]: %-10s,  style = %-10s,  def =",
              i,names[i],varstyles[style[i]]);
      if (style[i] == INTERNAL) {
        fprintf(out,"%g\n",input->variable->dvalue[i]);
        continue;
      }
      if ((style[i] != LOOP) && (style[i] != ULOOP))
        ndata = input->variable->num[i];
      for (int j=0; j < ndata; ++j)
        fprintf(out," %s",data[i][j]);
      fputs("\n",out);
    }
  }

  if (flags & TIME) {
    double wallclock = MPI_Wtime() - lmp->initclock;
    double cpuclock = 0.0;

#if defined(_WIN32)
    // from MSD docs.
    FILETIME ct,et,kt,ut;
    union { FILETIME ft; uint64_t ui; } cpu;
    if (GetProcessTimes(GetCurrentProcess(),&ct,&et,&kt,&ut)) {
      cpu.ft = ut;
      cpuclock = cpu.ui * 0.0000001;
    }
#else /* POSIX */
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) == 0) {
      cpuclock  = (double) ru.ru_utime.tv_sec;
      cpuclock += (double) ru.ru_utime.tv_usec * 0.000001;
    }
#endif /* ! _WIN32 */

    int cpuh,cpum,cpus,wallh,wallm,walls;
    cpus = fmod(cpuclock,60.0);
    cpuclock = (cpuclock - cpus) / 60.0;
    cpum = fmod(cpuclock,60.0);
    cpuh = (cpuclock - cpum) / 60.0;
    walls = fmod(wallclock,60.0);
    wallclock = (wallclock - walls) / 60.0;
    wallm = fmod(wallclock,60.0);
    wallh = (wallclock - wallm) / 60.0;
    fprintf(out,"\nTotal time information (MPI rank 0):\n"
            "  CPU time: %4d:%02d:%02d\n"
            " Wall time: %4d:%02d:%02d\n",
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

  fprintf(out,"\nStyles information:\n");

  if(flags & ATOM_STYLES)      atom_styles(out);
  if(flags & INTEGRATE_STYLES) integrate_styles(out);
  if(flags & MINIMIZE_STYLES)  minimize_styles(out);
  if(flags & PAIR_STYLES)      pair_styles(out);
  if(flags & BOND_STYLES)      bond_styles(out);
  if(flags & ANGLE_STYLES)     angle_styles(out);
  if(flags & DIHEDRAL_STYLES)  dihedral_styles(out);
  if(flags & IMPROPER_STYLES)  improper_styles(out);
  if(flags & KSPACE_STYLES)    kspace_styles(out);
  if(flags & FIX_STYLES)       fix_styles(out);
  if(flags & COMPUTE_STYLES)   compute_styles(out);
  if(flags & REGION_STYLES)    region_styles(out);
  if(flags & DUMP_STYLES)      dump_styles(out);
  if(flags & COMMAND_STYLES)   command_styles(out);
}

void Info::atom_styles(FILE * out)
{
  fprintf(out, "\nAtom styles:\n");
  print_columns(out, atom->avec_map);
  fprintf(out, "\n\n\n");
}

void Info::integrate_styles(FILE * out)
{
  fprintf(out, "\nIntegrate styles:\n");
  print_columns(out, update->integrate_map);
  fprintf(out, "\n\n\n");
}

void Info::minimize_styles(FILE * out)
{
  fprintf(out, "\nMinimize styles:\n");
  print_columns(out, update->minimize_map);
  fprintf(out, "\n\n\n");
}

void Info::pair_styles(FILE * out)
{
  fprintf(out, "\nPair styles:\n");
  print_columns(out, force->pair_map);
  fprintf(out, "\n\n\n");
}

void Info::bond_styles(FILE * out)
{
  fprintf(out, "\nBond styles:\n");
  print_columns(out, force->bond_map);
  fprintf(out, "\n\n\n");
}

void Info::angle_styles(FILE * out)
{
  fprintf(out, "\nAngle styles:\n");
  print_columns(out, force->angle_map);
  fprintf(out, "\n\n\n");
}

void Info::dihedral_styles(FILE * out)
{
  fprintf(out, "\nDihedral styles:\n");
  print_columns(out, force->dihedral_map);
  fprintf(out, "\n\n\n");
}

void Info::improper_styles(FILE * out)
{
  fprintf(out, "\nImproper styles:\n");
  print_columns(out, force->improper_map);
  fprintf(out, "\n\n\n");
}

void Info::kspace_styles(FILE * out)
{
  fprintf(out, "\nKSpace styles:\n");
  print_columns(out, force->kspace_map);
  fprintf(out, "\n\n\n");
}

void Info::fix_styles(FILE * out)
{
  fprintf(out, "\nFix styles:\n");
  print_columns(out, modify->fix_map);
  fprintf(out, "\n\n\n");
}

void Info::compute_styles(FILE * out)
{
  fprintf(out, "\nCompute styles:\n");
  print_columns(out, modify->compute_map);
  fprintf(out, "\n\n\n");
}

void Info::region_styles(FILE * out)
{
  fprintf(out, "\nRegion styles:\n");
  print_columns(out, domain->region_map);
  fprintf(out, "\n\n\n");
}

void Info::dump_styles(FILE * out)
{
  fprintf(out, "\nDump styles:\n");
  print_columns(out, output->dump_map);
  fprintf(out, "\n\n\n");
}

void Info::command_styles(FILE * out)
{
  fprintf(out, "\nCommand styles (add-on input script commands):\n");
  print_columns(out, input->command_map);
  fprintf(out, "\n\n\n");
}


/* ---------------------------------------------------------------------- */

// the is_active() function returns true if the selected style or name
// in the selected category is currently in use.

bool Info::is_active(const char *category, const char *name)
{
  if ((category == NULL) || (name == NULL)) return false;
  const char *style = "none";
  const int len = strlen(name);

  if (strcmp(category,"package") == 0) {
    if (strcmp(name,"gpu") == 0) {
      return (modify->find_fix("package_gpu") >= 0) ? true : false;
    } else if (strcmp(name,"intel") == 0) {
      return (modify->find_fix("package_intel") >= 0) ? true : false;
    } else if (strcmp(name,"kokkos") == 0) {
      return (lmp->kokkos && lmp->kokkos->kokkos_exists) ? true : false;
    } else if (strcmp(name,"omp") == 0) {
      return (modify->find_fix("package_omp") >= 0) ? true : false;
    } else error->all(FLERR,"Unknown name for info package category");

  } else if (strcmp(category,"newton") == 0) {
    if (strcmp(name,"pair") == 0) return (force->newton_pair != 0);
    else if (strcmp(name,"bond") == 0) return (force->newton_bond != 0);
    else if (strcmp(name,"any") == 0) return (force->newton != 0);
    else error->all(FLERR,"Unknown name for info newton category");

  } else if (strcmp(category,"pair") == 0) {
    if (force->pair == NULL) return false;
    if (strcmp(name,"single") == 0) return (force->pair->single_enable != 0);
    else if (strcmp(name,"respa") == 0) return (force->pair->respa_enable != 0);
    else if (strcmp(name,"manybody") == 0) return (force->pair->manybody_flag != 0);
    else if (strcmp(name,"tail") == 0) return (force->pair->tail_flag != 0);
    else if (strcmp(name,"shift") == 0) return (force->pair->offset_flag != 0);
    else error->all(FLERR,"Unknown name for info pair category");

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
  } else error->all(FLERR,"Unknown category for info is_active()");

  int match = 0;
  if (strcmp(style,name) == 0) match = 1;

  if (!match && lmp->suffix_enable) {
    if (lmp->suffix) {
      char *name_w_suffix = new char [len + 2 + strlen(lmp->suffix)];
      sprintf(name_w_suffix,"%s/%s",name,lmp->suffix);
      if (strcmp(style,name_w_suffix) == 0) match = 1;
      delete[] name_w_suffix;
    }
    if (!match && lmp->suffix2) {
      char *name_w_suffix = new char [len + 2 + strlen(lmp->suffix2)];
      sprintf(name_w_suffix,"%s/%s",name,lmp->suffix2);
      if (strcmp(style,name_w_suffix) == 0) match = 1;
      delete[] name_w_suffix;
    }
  }
  return match ? true : false;
}

/* ---------------------------------------------------------------------- */

// the is_available() function returns true if the selected style
// or name in the selected category is available for use (but need
// not be currently active).

bool Info::is_available(const char *category, const char *name)
{
  if ((category == NULL) || (name == NULL)) return false;
  const int len = strlen(name);
  int match = 0;

  if (strcmp(category,"command") == 0) {
    if (input->command_map->find(name) != input->command_map->end())
      match = 1;

  } else if (strcmp(category,"compute") == 0) {
    if (modify->compute_map->find(name) != modify->compute_map->end())
      match = 1;

    if (!match && lmp->suffix_enable) {
      if (lmp->suffix) {
        char *name_w_suffix = new char [len + 2 + strlen(lmp->suffix)];
        sprintf(name_w_suffix,"%s/%s",name,lmp->suffix);
        if (modify->compute_map->find(name_w_suffix) != modify->compute_map->end())
          match = 1;
        delete[] name_w_suffix;
      }
      if (!match && lmp->suffix2) {
        char *name_w_suffix = new char [len + 2 + strlen(lmp->suffix2)];
        sprintf(name_w_suffix,"%s/%s",name,lmp->suffix2);
        if (modify->compute_map->find(name_w_suffix) != modify->compute_map->end())
          match = 1;
        delete[] name_w_suffix;
      }
    }
  } else if (strcmp(category,"fix") == 0) {
    if (modify->fix_map->find(name) != modify->fix_map->end())
      match = 1;

    if (!match && lmp->suffix_enable) {
      if (lmp->suffix) {
        char *name_w_suffix = new char [len + 2 + strlen(lmp->suffix)];
        sprintf(name_w_suffix,"%s/%s",name,lmp->suffix);
        if (modify->fix_map->find(name_w_suffix) != modify->fix_map->end())
          match = 1;
        delete[] name_w_suffix;
      }
      if (!match && lmp->suffix2) {
        char *name_w_suffix = new char [len + 2 + strlen(lmp->suffix2)];
        sprintf(name_w_suffix,"%s/%s",name,lmp->suffix2);
        if (modify->fix_map->find(name_w_suffix) != modify->fix_map->end())
          match = 1;
        delete[] name_w_suffix;
      }
    }
  } else if (strcmp(category,"pair_style") == 0) {
    if (force->pair_map->find(name) != force->pair_map->end())
      match = 1;

    if (!match && lmp->suffix_enable) {
      if (lmp->suffix) {
        char *name_w_suffix = new char [len + 2 + strlen(lmp->suffix)];
        sprintf(name_w_suffix,"%s/%s",name,lmp->suffix);
        if (force->pair_map->find(name_w_suffix) != force->pair_map->end())
          match = 1;
        delete[] name_w_suffix;
      }
      if (!match && lmp->suffix2) {
        char *name_w_suffix = new char [len + 2 + strlen(lmp->suffix2)];
        sprintf(name_w_suffix,"%s/%s",name,lmp->suffix2);
        if (force->pair_map->find(name_w_suffix) != force->pair_map->end())
          match = 1;
        delete[] name_w_suffix;
      }
    }
  } else if (strcmp(category,"feature") == 0) {
    if (strcmp(name,"gzip") == 0) {
      return has_gzip_support();
    } else if (strcmp(name,"png") == 0) {
      return has_png_support();
    } else if (strcmp(name,"jpeg") == 0) {
      return has_jpeg_support();
    } else if (strcmp(name,"ffmpeg") == 0) {
      return has_ffmpeg_support();
    } else if (strcmp(name,"exceptions") == 0) {
      return has_exceptions();
    }
  } else error->all(FLERR,"Unknown category for info is_available()");

  return match ? true : false;
}

/* ---------------------------------------------------------------------- */

// the is_defined() function returns true if a particular ID of the
// selected category (e.g. fix ID, group ID, region ID etc.) has been
// defined and thus can be accessed. It does *NOT* check whether a
// particular ID has a particular style.

bool Info::is_defined(const char *category, const char *name)
{
  if ((category == NULL) || (name == NULL)) return false;

  if (strcmp(category,"compute") == 0) {
    int ncompute = modify->ncompute;
    Compute **compute = modify->compute;
    for (int i=0; i < ncompute; ++i) {
      if (strcmp(compute[i]->id,name) == 0)
        return true;
    }
  } else if (strcmp(category,"dump") == 0) {
    int ndump = output->ndump;
    Dump **dump = output->dump;
    for (int i=0; i < ndump; ++i) {
      if (strcmp(dump[i]->id,name) == 0)
        return true;
    }
  } else if (strcmp(category,"fix") == 0) {
    int nfix = modify->nfix;
    Fix **fix = modify->fix;
    for (int i=0; i < nfix; ++i) {
      if (strcmp(fix[i]->id,name) == 0)
        return true;
    }
  } else if (strcmp(category,"group") == 0) {
    int ngroup = group->ngroup;
    char **names = group->names;
    for (int i=0; i < ngroup; ++i) {
      if (strcmp(names[i],name) == 0)
        return true;
    }
  } else if (strcmp(category,"region") == 0) {
    int nreg = domain->nregion;
    Region **regs = domain->regions;
    for (int i=0; i < nreg; ++i) {
      if (strcmp(regs[i]->id,name) == 0)
        return true;
    }
  } else if (strcmp(category,"variable") == 0) {
    int nvar = input->variable->nvar;
    char **names = input->variable->names;
    for (int i=0; i < nvar; ++i) {
      if (strcmp(names[i],name) == 0)
        return true;
    }
  } else error->all(FLERR,"Unknown category for info is_defined()");

  return false;
}

template<typename ValueType>
static void print_columns(FILE* fp, map<string, ValueType> * styles)
{
  if (styles->empty()) {
    fprintf(fp, "\nNone");
    return;
  }

  // std::map keys are already sorted
  int pos = 80;
  for(typename map<string, ValueType>::iterator it = styles->begin(); it != styles->end(); ++it) {
    const string & style_name = it->first;

    // skip "secret" styles
    if (isupper(style_name[0])) continue;

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

bool Info::has_exceptions() {
#ifdef LAMMPS_EXCEPTIONS
  return true;
#else
  return false;
#endif
}

bool Info::has_package(const char * package_name) {
  for(int i = 0; LAMMPS::installed_packages[i] != NULL; ++i) {
    if(strcmp(package_name, LAMMPS::installed_packages[i]) == 0) {
      return true;
    }
  }
  return false;
}

/* ---------------------------------------------------------------------- */
#define _INFOBUF_SIZE 256

char *Info::get_os_info()
{
  char *buf = new char[_INFOBUF_SIZE];

#if defined(_WIN32)
  DWORD fullversion,majorv,minorv,buildv=0;

  fullversion = GetVersion();
  majorv = (DWORD) (LOBYTE(LOWORD(fullversion)));
  minorv = (DWORD) (HIBYTE(LOWORD(fullversion)));
  if (fullversion < 0x80000000)
    buildv = (DWORD) (HIWORD(fullversion));

  SYSTEM_INFO si;
  GetSystemInfo(&si);

  const char *machine;
  switch (si.wProcessorArchitecture) {
  case PROCESSOR_ARCHITECTURE_AMD64:
    machine = (const char *) "x86_64";
    break;
  case PROCESSOR_ARCHITECTURE_ARM:
    machine = (const char *) "arm";
    break;
  case PROCESSOR_ARCHITECTURE_IA64:
    machine = (const char *) "ia64";
    break;
  case PROCESSOR_ARCHITECTURE_INTEL:
    machine = (const char *) "i386";
    break;
  default:
    machine = (const char *) "(unknown)";
  }
  snprintf(buf,_INFOBUF_SIZE,"Windows %d.%d (%d) on %s",
           majorv,minorv,buildv,machine);
#else
  struct utsname ut;
  uname(&ut);
  snprintf(buf,_INFOBUF_SIZE,"%s %s on %s",
           ut.sysname, ut.release, ut.machine);
#endif
  return buf;
}

char *Info::get_compiler_info()
{
  char *buf = new char[_INFOBUF_SIZE];
#if __clang__
  snprintf(buf,_INFOBUF_SIZE,"Clang C++ %s", __VERSION__);
#elif __INTEL_COMPILER
  snprintf(buf,_INFOBUF_SIZE,"Intel C++ %s", __VERSION__);
#elif __GNUC__
  snprintf(buf,_INFOBUF_SIZE,"GNU C++ %s",   __VERSION__);
#else
  snprintf(buf,_INFOBUF_SIZE,"(Unknown)");
#endif
  return buf;
}

const char *Info::get_openmp_info()
{

#if !defined(_OPENMP)
  return (const char *)"OpenMP not enabled";
#else

// Supported OpenMP version corresponds to the release date of the
// specifications as posted at https://www.openmp.org/specifications/

#if _OPENMP > 201811
  return (const char *)"OpenMP newer than version 5.0";
#elif _OPENMP == 201811
  return (const char *)"OpenMP 5.0";
#elif _OPENMP == 201611
  return (const char *)"OpenMP 5.0 preview 1";
#elif _OPENMP == 201511
  return (const char *)"OpenMP 4.5";
#elif _OPENMP == 201307
  return (const char *)"OpenMP 4.0";
#elif _OPENMP == 201107
  return (const char *)"OpenMP 3.1";
#elif _OPENMP == 200805
  return (const char *)"OpenMP 3.0";
#elif _OPENMP == 200505
  return (const char *)"OpenMP 2.5";
#elif _OPENMP == 200203
  return (const char *)"OpenMP 2.0";
#else
  return (const char *)"unknown OpenMP version";
#endif

#endif
}

/* ---------------------------------------------------------------------- */

char **Info::get_variable_names(int &num) {
    num = input->variable->nvar;
    return input->variable->names;
}
