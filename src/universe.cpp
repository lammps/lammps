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

#include "universe.h"
#include <mpi.h>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include "version.h"
#include "error.h"
#include "force.h"
#include "memory.h"

using namespace LAMMPS_NS;

#define MAXLINE 256

static char *date2num(const char *version);

/* ----------------------------------------------------------------------
   create & initialize the universe of processors in communicator
------------------------------------------------------------------------- */

Universe::Universe(LAMMPS *lmp, MPI_Comm communicator) : Pointers(lmp)
{
  version = (const char *) LAMMPS_VERSION;
  num_ver = date2num(version);

  uworld = uorig = communicator;
  MPI_Comm_rank(uworld,&me);
  MPI_Comm_size(uworld,&nprocs);

  uscreen = stdout;
  ulogfile = NULL;

  existflag = 0;
  nworlds = 0;
  procs_per_world = NULL;
  root_proc = NULL;

  memory->create(uni2orig,nprocs,"universe:uni2orig");
  for (int i = 0; i < nprocs; i++) uni2orig[i] = i;
}

/* ---------------------------------------------------------------------- */

Universe::~Universe()
{
  if (uworld != uorig) MPI_Comm_free(&uworld);
  memory->destroy(procs_per_world);
  memory->destroy(root_proc);
  memory->destroy(uni2orig);
  delete [] num_ver;
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
  char line[MAXLINE];

  if (uworld != uorig) MPI_Comm_free(&uworld);

  if (strcmp(style,"nth") == 0) {
    int n = force->inumeric(FLERR,arg);
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
      FILE *fp = fopen(arg,"r");
      if (fp == NULL) error->universe_one(FLERR,"Cannot open -reorder file");

      // skip header = blank and comment lines

      char *ptr;
      if (!fgets(line,MAXLINE,fp))
        error->one(FLERR,"Unexpected end of -reorder file");
      while (1) {
        if ((ptr = strchr(line,'#'))) *ptr = '\0';
        if (strspn(line," \t\n\r") != strlen(line)) break;
        if (!fgets(line,MAXLINE,fp))
          error->one(FLERR,"Unexpected end of -reorder file");
      }

      // read nprocs lines
      // uni2orig = inverse mapping

      int me_orig,me_new,rv;
      rv = sscanf(line,"%d %d",&me_orig,&me_new);
      if (me_orig < 0 || me_orig >= nprocs ||
          me_new < 0 || me_new >= nprocs || rv != 2)
        error->one(FLERR,"Invalid entry in -reorder file");
      uni2orig[me_new] = me_orig;

      for (int i = 1; i < nprocs; i++) {
        if (!fgets(line,MAXLINE,fp))
          error->one(FLERR,"Unexpected end of -reorder file");
        rv = sscanf(line,"%d %d",&me_orig,&me_new);
        if (me_orig < 0 || me_orig >= nprocs ||
            me_new < 0 || me_new >= nprocs || rv != 2)
          error->one(FLERR,"Invalid entry in -reorder file");
        uni2orig[me_new] = me_orig;
      }
      fclose(fp);
    }

    // bcast uni2org from proc 0 to all other universe procs

    MPI_Bcast(uni2orig,nprocs,MPI_INT,0,uorig);

  } else error->universe_all(FLERR,"Invalid command-line argument");

  // create new uworld communicator

  int ome,key;
  MPI_Comm_rank(uorig,&ome);
  for (int i = 0; i < nprocs; i++)
    if (uni2orig[i] == ome) key = i;

  MPI_Comm_split(uorig,0,key,&uworld);
  MPI_Comm_rank(uworld,&me);
  MPI_Comm_size(uworld,&nprocs);
}

/* ----------------------------------------------------------------------
   add 1 or more worlds to universe
   str == NULL -> add 1 world with all procs in universe
   str = NxM -> add N worlds, each with M procs
   str = P -> add 1 world with P procs
------------------------------------------------------------------------- */

void Universe::add_world(char *str)
{
  int n,nper;
  char *ptr;

  n = 1;
  nper = 0;

  if (str != NULL) {

    // check for valid partition argument

    bool valid = true;

    // str may not be empty and may only consist of digits or 'x'

    size_t len = strlen(str);
    if (len < 1) valid = false;
    for (size_t i=0; i < len; ++i)
      if (isdigit(str[i]) || str[i] == 'x') continue;
      else valid = false;

    if (valid) {
      if ((ptr = strchr(str,'x')) != NULL) {

        // 'x' may not be the first or last character

        if (ptr == str) {
          valid = false;
        } else if (strlen(str) == len-1) {
          valid = false;
        } else {
          *ptr = '\0';
          n = atoi(str);
          nper = atoi(ptr+1);
          *ptr = 'x';
        }
      } else nper = atoi(str);
    }

    // require minimum of 1 partition with 1 processor

    if (n < 1 || nper < 1) valid = false;

    if (!valid) {
      char msg[128];
      snprintf(msg,128,"Invalid partition string '%s'",str);
      error->universe_all(FLERR,msg);
    }
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

// helper function to convert the LAMMPS date string to a version id
// that can be used for both string and numerical comparisons
// where newer versions are larger than older ones.

char *date2num(const char *version)
{
  int day,month,year;
  day = month = year = 0;

  if (version) {

    day = atoi(version);

    while (*version != '\0' && (isdigit(*version) || *version == ' '))
      ++version;

    if (strncmp(version,"Jan",3) == 0) month = 1;
    if (strncmp(version,"Feb",3) == 0) month = 2;
    if (strncmp(version,"Mar",3) == 0) month = 3;
    if (strncmp(version,"Apr",3) == 0) month = 4;
    if (strncmp(version,"May",3) == 0) month = 5;
    if (strncmp(version,"Jun",3) == 0) month = 6;
    if (strncmp(version,"Jul",3) == 0) month = 7;
    if (strncmp(version,"Aug",3) == 0) month = 8;
    if (strncmp(version,"Sep",3) == 0) month = 9;
    if (strncmp(version,"Oct",3) == 0) month = 10;
    if (strncmp(version,"Nov",3) == 0) month = 11;
    if (strncmp(version,"Dec",3) == 0) month = 12;

    while (*version != '\0' && !isdigit(*version))
      ++version;

    year = atoi(version);
  }

  char *ver = new char[64];
  sprintf(ver,"%04d%02d%02d", year % 10000, month, day % 100);

  return ver;
}
