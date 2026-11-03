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

#include "universe.h"

#include "error.h"
#include "memory.h"
#include "safe_pointers.h"
#include "text_file_reader.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   create & initialize the universe of processors in communicator
------------------------------------------------------------------------- */

Universe::Universe(LAMMPS *lmp, MPI_Comm communicator) :
    Pointers(lmp), uscreen(stdout), ulogfile(nullptr), procs_per_world(nullptr), root_proc(nullptr),
    uni2orig(nullptr)
{
  uworld = uorig = communicator;
  MPI_Comm_rank(uworld, &me);
  MPI_Comm_size(uworld, &nprocs);

  existflag = 0;
  nworlds = 0;

  memory->create(uni2orig, nprocs, "universe:uni2orig");
  for (int i = 0; i < nprocs; i++) uni2orig[i] = i;
}

/* ---------------------------------------------------------------------- */

Universe::~Universe()
{
  if (uworld != uorig) MPI_Comm_free(&uworld);
  memory->destroy(procs_per_world);
  memory->destroy(root_proc);
  memory->destroy(uni2orig);
}

/* ----------------------------------------------------------------------
   reorder universe processors
   create uni2orig as inverse mapping
   re-create uworld communicator with new ordering via Comm_split()
   style = "nth", arg = N
   move every Nth proc to end of rankings
   style = "custom", arg = filename
   file has nprocs lines with I J
   I = universe proc ID in original communicator uorig
   J = universe proc ID in reordered communicator uworld
------------------------------------------------------------------------- */

void Universe::reorder(char *style, char *arg)
{
  if (uworld != uorig) MPI_Comm_free(&uworld);

  if (strcmp(style,"nth") == 0) {
    int n = utils::inumeric(FLERR,arg,false,lmp);
    if (n <= 0)
      error->universe_all(FLERR,"Invalid -reorder N value");
    if (nprocs % n)
      error->universe_all(FLERR,"Nprocs not a multiple of N for -reorder");
    for (int i = 0; i < nprocs; i++) {
      if (i < (n-1)*nprocs/n) uni2orig[i] = i/(n-1) * n + (i % (n-1));
      else uni2orig[i] = (i - (n-1)*nprocs/n) * n + n-1;
    }

  } else if (strcmp(style,"custom") == 0) {

    if (me == 0) {
      SafeFilePtr fp = fopen(arg,"r");
      if (fp == nullptr)
        error->universe_one(FLERR,fmt::format("Cannot open -reorder file {}: {}", arg,
                                              utils::getsyserror()));

      // read nprocs lines skipping over blank and comment lines
      // uni2orig = inverse mapping
      // each original and each new rank must appear exactly once

      std::vector<int> seen_orig(nprocs,0), seen_new(nprocs,0);

      try {
        TextFileReader reader(fp, "-reorder");
        for (int i = 0; i < nprocs; i++) {
          auto values = reader.next_values(2);
          int me_orig = values.next_int();
          int me_new = values.next_int();
          if ((me_orig < 0) || (me_orig >= nprocs) || (me_new < 0) || (me_new >= nprocs))
            error->one(FLERR,"Invalid entry '{} {}' in -reorder file", me_orig, me_new);
          if (seen_orig[me_orig] || seen_new[me_new])
            error->one(FLERR,"Duplicate entry '{} {}' in -reorder file", me_orig, me_new);
          seen_orig[me_orig] = seen_new[me_new] = 1;
          uni2orig[me_new] = me_orig;
        }
      } catch (std::exception &e) {
        error->one(FLERR,"Error reading -reorder file: {}", e.what());
      }
    }

    // bcast uni2org from proc 0 to all other universe procs

    MPI_Bcast(uni2orig,nprocs,MPI_INT,0,uorig);

  } else error->universe_all(FLERR,"Invalid command-line argument");

  // create new uworld communicator

  int ome,key = -1;
  MPI_Comm_rank(uorig,&ome);
  for (int i = 0; i < nprocs; i++)
    if (uni2orig[i] == ome) key = i;
  if (key < 0) error->universe_one(FLERR,"Invalid -reorder file: not a permutation of all ranks");

  MPI_Comm_split(uorig,0,key,&uworld);
  MPI_Comm_rank(uworld,&me);
  MPI_Comm_size(uworld,&nprocs);
}

/* ----------------------------------------------------------------------
   add 1 or more worlds to universe
   str == nullptr -> add 1 world with all procs in universe
   str = NxM -> add N worlds, each with M procs
   str = P -> add 1 world with P procs
------------------------------------------------------------------------- */

void Universe::add_world(char *str)
{
  int n,nper;

  n = 1;
  nper = 0;

  if (str != nullptr) {

    // check for valid partition argument

    bool valid = true;

    // str may not be empty and may only consist of digits or 'x'

    std::string part(str);
    if (part.empty()) valid = false;
    if (part.find_first_not_of("0123456789x") != std::string::npos) valid = false;

    if (valid) {
      std::size_t found = part.find_first_of('x');

      // 'x' may not be the first or last character

      if ((found == 0) || (found == (part.size() - 1))) {
        valid = false;
      } else if (found == std::string::npos) {
        nper = std::stoi(part);
      } else {
        n = std::stoi(part.substr(0,found));
        nper = std::stoi(part.substr(found+1));
      }
    }

    // require minimum of 1 partition with 1 processor

    if (n < 1 || nper < 1) valid = false;

    if (!valid)
      error->universe_all(FLERR, fmt::format("Invalid partition string '{}'", str));
  } else nper = nprocs;

  memory->grow(procs_per_world,nworlds+n,"universe:procs_per_world");
  memory->grow(root_proc,(nworlds+n),"universe:root_proc");

  for (int i = 0; i < n; i++) {
    procs_per_world[nworlds] = nper;
    if (nworlds == 0) root_proc[nworlds] = 0;
    else
      root_proc[nworlds] = root_proc[nworlds-1] + procs_per_world[nworlds-1];
    if (me >= root_proc[nworlds]) iworld = nworlds;
    nworlds++;
  }
}

/* ----------------------------------------------------------------------
   check if total procs in all worlds = procs in universe
------------------------------------------------------------------------- */

int Universe::consistent()
{
  int n = 0;
  for (int i = 0; i < nworlds; i++) n += procs_per_world[i];
  if (n == nprocs) return 1;
  else return 0;
}
