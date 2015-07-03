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

#include "string.h"
#include "info.h"
#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "dump.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "modify.h"
#include "output.h"
#include "region.h"
#include "variable.h"
#include "update.h"
#include "error.h"

#ifdef _WIN32
#include <windows.h>
#include <stdint.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif

using namespace LAMMPS_NS;

// same as in variable.cpp
enum{INDEX,LOOP,WORLD,UNIVERSE,ULOOP,STRING,GETENV,
     SCALARFILE,ATOMFILE,FORMAT,EQUAL,ATOM,PYTHON};
static const char *varstyles[] = {
  "index", "loop", "world", "universe", "uloop", "string", "getenv",
  "file", "atomfile", "format", "equal", "atom", "python", "(unknown)"};

static const char *mapstyles[] = { "none", "array", "hash" };

static const char bstyles[] = "pfsm";

/* ---------------------------------------------------------------------- */

void Info::command(int narg, char **arg)
{
  if (!screen || comm->me) return;

  for (int idx = 0; idx < narg; ++idx) {
    fputs("\n",screen);

    if (strcmp(arg[idx],"computes") == 0) {
      int ncompute = modify->ncompute;
      Compute **compute = modify->compute;
      char **names = group->names;
      fprintf(screen,"Compute information:\n");
      for (int i=0; i < ncompute; ++i) {
        fprintf(screen,"Compute[%3d]: %s,  style: %s,  group: %s\n",
                i, compute[i]->id, compute[i]->style, names[compute[i]->igroup]);
      }

    } else if (strcmp(arg[idx],"dumps") == 0) {
      int ndump = output->ndump;
      Dump **dump = output->dump;
      int *nevery = output->every_dump;\
      char **vnames = output->var_dump;
      char **names = group->names;
      fprintf(screen,"Dump information:\n");
      for (int i=0; i < ndump; ++i) {
        fprintf(screen,"Dump[%3d]: %s,  file: %s,  style: %s,  group: %s,  ",
                i, dump[i]->id, dump[i]->filename,
                dump[i]->style, names[dump[i]->igroup]);
        if (nevery[i]) {
          fprintf(screen,"every: %d\n", nevery[i]);
        } else {
          fprintf(screen,"every: %s\n", vnames[i]);
        }
      }

    } else if (strcmp(arg[idx],"fixes") == 0) {
      int nfix = modify->nfix;
      Fix **fix = modify->fix;
      char **names = group->names;
      fprintf(screen,"Fix information:\n");
      for (int i=0; i < nfix; ++i) {
        fprintf(screen,"Fix[%3d]: %s,  style: %s,  group: %s\n",
                i, fix[i]->id, fix[i]->style, names[fix[i]->igroup]);
      }

    } else if (strcmp(arg[idx],"groups") == 0) {
      int ngroup = group->ngroup;
      char **names = group->names;
      int *dynamic = group->dynamic;
      fprintf(screen,"Group information:\n");
      for (int i=0; i < ngroup; ++i) {
        fprintf(screen,"Group[%2d]: %s (%s)\n",
                i, names[i], dynamic[i] ? "dynamic" : "static");
      }

    } else if (strcmp(arg[idx],"regions") == 0) {
      int nreg = domain->nregion;
      Region **regs = domain->regions;
      fprintf(screen,"Region information:\n");
      for (int i=0; i < nreg; ++i) {
        fprintf(screen,"Region[%3d]: %s,  style: %s,  side: %s\n",
                i, regs[i]->id, regs[i]->style, regs[i]->interior ? "in" : "out");
      }

    } else if (strcmp(arg[idx],"memory") == 0) {

    } else if (strcmp(arg[idx],"time") == 0) {

      double wallclock = MPI_Wtime() - lmp->initclock;
      double cpuclock = 0.0;
      
#ifdef _WIN32

      // from MSD docs.
      FILETIME ct,et,kt,ut;
      union { FILETIME ft; uint64_t ui; } cpu;
      if (GetProcessTimes(GetCurrentProcess(),&ct,&et,&kt,&ut)) {
        cpu.ft = ut;
        cpuclock = cpu.ui * 0.0000001;
      }

#else /* ! _WIN32 */

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
      fprintf(screen,"Total time information (MPI rank 0):\n"
              "  CPU time: %4d:%02d:%02d\n"
              " Wall time: %4d:%02d:%02d\n",
              cpuh,cpum,cpus,wallh,wallm,walls);

    } else if (strcmp(arg[idx],"variables") == 0) {
      int nvar = input->variable->nvar;
      int *style = input->variable->style;
      char **names = input->variable->names;
      char ***data = input->variable->data;
      fprintf(screen,"Variable information:\n");
      for (int i=0; i < nvar; ++i) {
        int ndata = 1;
        fprintf(screen,"Variable[%3d]: %-10s  style = %-10s  def =",
                i,names[i],varstyles[style[i]]);
        if ((style[i] != LOOP) && (style[i] != ULOOP))
          ndata = input->variable->num[i];
        for (int j=0; j < ndata; ++j)
          fprintf(screen," %s",data[i][j]);
        fputs("\n",screen);
      }

    } else if (strcmp(arg[idx],"system") == 0) {
      fprintf(screen,"System information:\n");
      fprintf(screen,"Units      = %s\n",update->unit_style);
      fprintf(screen,"Atom style = %s\n", atom->atom_style);
      fprintf(screen,"Atom map   = %s\n", mapstyles[atom->map_style]);
      if (atom->molecular > 0) {
        const char *msg;
        msg = (atom->molecular == 2) ? "template" : "standard";
        fprintf(screen,"Molecule type = %s\n",msg);
      }
      fprintf(screen,"Atoms: " BIGINT_FORMAT ",  types: %d,  style: %s\n",
              atom->natoms, atom->ntypes, force->pair_style);

      if (atom->molecular > 0) {
        const char *msg;
        msg = force->bond_style ? force->bond_style : "none";
        fprintf(screen,"Bonds: " BIGINT_FORMAT ",  types: %d,  style: %s\n",
                atom->nbonds, atom->nbondtypes, msg);

        msg = force->angle_style ? force->angle_style : "none";
        fprintf(screen,"Angles: " BIGINT_FORMAT ",  types: %d,  style: %s\n",
                atom->nangles, atom->nangletypes, msg);

        msg = force->dihedral_style ? force->dihedral_style : "none";
        fprintf(screen,"Dihedrals: " BIGINT_FORMAT ",  types: %d,  style: %s\n",
                atom->ndihedrals, atom->ndihedraltypes, msg);

        msg = force->improper_style ? force->improper_style : "none";
        fprintf(screen,"Impropers: " BIGINT_FORMAT ",  types: %d,  style: %s\n",
                atom->nimpropers, atom->nimpropertypes, msg);

        const double * const special_lj   = force->special_lj;
        const double * const special_coul = force->special_coul;

        fprintf(screen,"Special bond factors lj:   %-10g %-10g %-10g\n"
                "Special bond factors coul: %-10g %-10g %-10g\n",
                special_lj[1],special_lj[2],special_lj[3],
                special_coul[1],special_coul[2],special_coul[3]);
      }

      fprintf(screen,"Kspace style = %s\n",
              force->kspace ? force->kspace_style : "none");

      if (domain->box_exist) {
        fprintf(screen,"\n%s box: (%g x %g %g)\n",
                domain->triclinic ? "Triclinic" : "Orthogonal",
                domain->xprd, domain->yprd, domain->zprd);
        fprintf(screen,"Boundaries = %c,%c %c,%c %c,%c\n",
                bstyles[domain->boundary[0][0]],bstyles[domain->boundary[0][1]],
                bstyles[domain->boundary[1][0]],bstyles[domain->boundary[1][1]],
                bstyles[domain->boundary[2][0]],bstyles[domain->boundary[2][1]]);
        fprintf(screen,"Xlo, zhi = %g, %g\n", domain->boxlo[0], domain->boxhi[0]);
        fprintf(screen,"Ylo, zhi = %g, %g\n", domain->boxlo[1], domain->boxhi[1]);
        fprintf(screen,"Zlo, zhi = %g, %g\n", domain->boxlo[2], domain->boxhi[2]);
        if (domain->triclinic)
          fprintf(screen,"Xy, xz, yz = %g, %g, %g\n",
                  domain->xy, domain->xz, domain->yz);

      } else {
        fputs("\nBox has not yet been created\n",screen);
      }
    } else {
      error->all(FLERR,"Unknown info command style");
    }
  }
}
