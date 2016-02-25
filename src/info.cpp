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
   Contributing author:  Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include <string.h>
#include "info.h"
#include "accelerator_cuda.h"
#include "accelerator_kokkos.h"
#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "dump.h"
#include "fix.h"
#include "force.h"
#include "pair.h"
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

#include <time.h>

#ifdef _WIN32
#define PSAPI_VERSION=1
#include <windows.h>
#include <stdint.h>
#include <psapi.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/utsname.h>
#endif

#if defined __linux
#include <malloc.h>
#endif

namespace LAMMPS_NS {
// same as in variable.cpp
enum {INDEX,LOOP,WORLD,UNIVERSE,ULOOP,STRING,GETENV,
      SCALARFILE,ATOMFILE,FORMAT,EQUAL,ATOM,PYTHON};

enum {COMPUTES=1<<0,
      DUMPS=1<<1,
      FIXES=1<<2,
      GROUPS=1<<3,
      REGIONS=1<<4,
      CONFIG=1<<5,
      TIME=1<<6,
      VARIABLES=1<<7,
      SYSTEM=1<<8,
      COMM=1<<9,
      ALL=~0};
}

static const char *varstyles[] = {
  "index", "loop", "world", "universe", "uloop", "string", "getenv",
  "file", "atomfile", "format", "equal", "atom", "python", "(unknown)"};

static const char *mapstyles[] = { "none", "array", "hash" };

static const char *commstyles[] = { "brick", "tiled" };
static const char *commlayout[] = { "uniform", "nonuniform", "irregular" };

static const char bstyles[] = "pfsm";

using namespace LAMMPS_NS;

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
    } else if (strncmp(arg[idx],"communication",4) == 0) {
      flags |= COMM;
      ++idx;
    } else if (strncmp(arg[idx],"computes",4) == 0) {
      flags |= COMPUTES;
      ++idx;
    } else if (strncmp(arg[idx],"dumps",3) == 0) {
      flags |= DUMPS;
      ++idx;
    } else if (strncmp(arg[idx],"fixes",3) == 0) {
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
    } else if (strncmp(arg[idx],"variables",3) == 0) {
      flags |= VARIABLES;
      ++idx;
    } else if (strncmp(arg[idx],"system",3) == 0) {
      flags |= SYSTEM;
      ++idx;
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

    fprintf(out,"\nLAMMPS version: %s / %s\n",
            universe->version, universe->num_ver);
    fprintf(out,"sizeof(smallint): %3d-bit\n",(int)sizeof(smallint)*8);
    fprintf(out,"sizeof(imageint): %3d-bit\n",(int)sizeof(imageint)*8);
    fprintf(out,"sizeof(tagint):   %3d-bit\n",(int)sizeof(tagint)*8);
    fprintf(out,"sizeof(bigint):   %3d-bit\n",(int)sizeof(bigint)*8);

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
    fprintf(out,"\nOS information: Windows %d.%d (%d) on %s\n",
            majorv,minorv,buildv,machine);
#else
    struct utsname ut;
    uname(&ut);
    fprintf(out,"\nOS information: %s %s on %s\n",
            ut.sysname, ut.release, ut.machine);
#endif

    fprintf(out,"\nMemory allocation information (MPI rank 0)\n");

#if defined(_WIN32)
    HANDLE phandle = GetCurrentProcess();
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(phandle,(PROCESS_MEMORY_COUNTERS *)&pmc,sizeof(pmc));
    fprintf(out,"Non-shared memory use: %.3g Mbyte\n",
            (double)pmc.PrivateUsage/1048576.0);
    fprintf(out,"Maximum working set size: %.3g Mbyte\n",
            (double)pmc.PeakWorkingSetSize/1048576.0);
#else
#if defined(__linux)
    struct mallinfo mi;
    mi = mallinfo();
    fprintf(out,"Total dynamically allocated memory: %.3g Mbyte\n",
            (double)mi.uordblks/1048576.0);
#endif
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) == 0) {
      fprintf(out,"Maximum resident set size: %.3g Mbyte\n",
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
    fprintf(out,"Nprocs = %d    Nthreads = %d\n",
            comm->nprocs, comm->nthreads);
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
      fprintf(out,"Xlo, zhi = %g, %g\n", domain->boxlo[0], domain->boxhi[0]);
      fprintf(out,"Ylo, zhi = %g, %g\n", domain->boxlo[1], domain->boxhi[1]);
      fprintf(out,"Zlo, zhi = %g, %g\n", domain->boxlo[2], domain->boxhi[2]);
      if (domain->triclinic)
          fprintf(out,"Xy, xz, yz = %g, %g, %g\n",
                  domain->xy, domain->xz, domain->yz);
    } else {
      fputs("\nBox has not yet been created\n",out);
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
      fprintf(out,"Variable[%3d]: %-10s  style = %-10s  def =",
              i,names[i],varstyles[style[i]]);
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

  fputs("\nInfo-Info-Info-Info-Info-Info-Info-Info-Info-Info-Info\n\n",out);

  // close output file pointer if opened locally thus forcing a hard sync.
  if ((out != screen) && (out != logfile))
    fclose(out);
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
    if (strcmp(name,"cuda") == 0) {
      return (lmp->cuda && lmp->cuda->cuda_exists) ? true : false;
    } else if (strcmp(name,"gpu") == 0) {
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
