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

#include "fix_print.h"

#include "error.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPrint::FixPrint(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  fp(nullptr), text(nullptr), copy(nullptr), work(nullptr), var_print(nullptr)
{
  if (narg < 5) error->all(FLERR,"Illegal fix print command");
  if (utils::strmatch(arg[3],"^v_")) {
    var_print = utils::strdup(arg[3]+2);
    nevery = 1;
  } else {
    nevery = utils::inumeric(FLERR,arg[3],false,lmp);
    if (nevery <= 0) error->all(FLERR,"Illegal fix print command");
  }

  MPI_Comm_rank(world,&me);

  text = utils::strdup(arg[4]);
  int n = strlen(text)+1;
  copy = (char *) memory->smalloc(n*sizeof(char),"fix/print:copy");
  work = (char *) memory->smalloc(n*sizeof(char),"fix/print:work");
  maxcopy = maxwork = n;

  // parse optional args

  fp = nullptr;
  screenflag = 1;
  char *title = nullptr;

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0 || strcmp(arg[iarg],"append") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix print command");
      if (me == 0) {
        if (strcmp(arg[iarg],"file") == 0) fp = fopen(arg[iarg+1],"w");
        else fp = fopen(arg[iarg+1],"a");
        if (fp == nullptr)
          error->one(FLERR,"Cannot open fix print file {}: {}",
                     arg[iarg+1], utils::getsyserror());
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"screen") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix print command");
      screenflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix print command");
      delete [] title;
      title = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix print command");
  }

  // print file comment line

  if (fp && me == 0) {
    if (title) fprintf(fp,"%s\n",title);
    else fprintf(fp,"# Fix print output for fix %s\n",id);
  }

  delete [] title;
}

/* ---------------------------------------------------------------------- */

FixPrint::~FixPrint()
{
  delete [] text;
  delete [] var_print;
  memory->sfree(copy);
  memory->sfree(work);

  if (fp && me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixPrint::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPrint::init()
{
  if (var_print) {
    ivar_print = input->variable->find(var_print);
    if (ivar_print < 0)
      error->all(FLERR,"Variable name for fix print timestep does not exist");
    if (!input->variable->equalstyle(ivar_print))
      error->all(FLERR,"Variable for fix print timestep is invalid style");
    next_print = static_cast<bigint>
      (input->variable->compute_equal(ivar_print));
    if (next_print <= update->ntimestep)
      error->all(FLERR,"Fix print timestep variable returned a bad timestep");
  } else {
    if (update->ntimestep % nevery)
      next_print = (update->ntimestep/nevery)*nevery + nevery;
    else
      next_print = update->ntimestep;
  }

  // add next_print to all computes that store invocation times
  // since don't know a priori which are invoked via variables by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  modify->addstep_compute_all(next_print);
}

/* ---------------------------------------------------------------------- */

void FixPrint::setup(int /* vflag */)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixPrint::end_of_step()
{
  if (update->ntimestep != next_print) return;

  // make a copy of text to work on
  // substitute for $ variables (no printing)
  // append a newline and print final copy
  // variable evaluation may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  strncpy(copy,text,maxcopy);
  input->substitute(copy,work,maxcopy,maxwork,0);

  if (var_print) {
    next_print = static_cast<bigint>
      (input->variable->compute_equal(ivar_print));
    if (next_print <= update->ntimestep)
      error->all(FLERR,"Fix print timestep variable returned a bad timestep");
  } else {
    next_print = (update->ntimestep/nevery)*nevery + nevery;
  }

  modify->addstep_compute(next_print);

  if (me == 0) {
    if (screenflag) utils::logmesg(lmp,std::string(copy) + "\n");
    if (fp) {
      fmt::print(fp,"{}\n",copy);
      fflush(fp);
    }
  }
}
