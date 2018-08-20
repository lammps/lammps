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
   Contributing authors: Axel Kohlmeyer (Temple U),
                         Ryan S. Elliott (UMN)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the Free
   Software Foundation; either version 2 of the License, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along with
   this program; if not, see <https://www.gnu.org/licenses>.

   Linking LAMMPS statically or dynamically with other modules is making a
   combined work based on LAMMPS. Thus, the terms and conditions of the GNU
   General Public License cover the whole combination.

   In addition, as a special exception, the copyright holders of LAMMPS give
   you permission to combine LAMMPS with free software programs or libraries
   that are released under the GNU LGPL and with code included in the standard
   release of the "kim-api" under the CDDL (or modified versions of such code,
   with unchanged license). You may copy and distribute such a system following
   the terms of the GNU GPL for LAMMPS and the licenses of the other code
   concerned, provided that you include the source code of that other code
   when and as the GNU GPL requires distribution of source code.

   Note that people who make modified versions of LAMMPS are not obligated to
   grant this special exception for their modified versions; it is their choice
   whether to do so. The GNU General Public License gives permission to release
   a modified version without this exception; this exception also makes it
   possible to release a modified version which carries forward this exception.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Designed for use with the kim-api-v2.0.0-beta.1 (and newer) package
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cstring>
#include "kim_query.h"
#include "comm.h"
#include "error.h"
#include "input.h"
#include "variable.h"

using namespace LAMMPS_NS;

static char *do_query(char *, char *, int, MPI_Comm);

/* ---------------------------------------------------------------------- */

void KimQuery::command(int narg, char **arg)
{
  char *model, *property, *varname, *value;

  if (narg != 3) error->all(FLERR,"Illegal kim_query command");

  model = arg[0];
  property = arg[1];
  varname = arg[2];

  value = do_query(model, property, comm->me, world);

  if (comm->me == 0)
    printf("property %s for model %s is %s\n",property,model,value);

  char **varcmd = new char*[3];
  varcmd[0] = varname;
  varcmd[1] = (char *) "string";
  varcmd[2] = value;

  input->variable->set(3,varcmd);
  
  delete[] varcmd;
  delete[] value;
}


char *do_query(char *model, char *property, int rank, MPI_Comm comm)
{
  char val[512], *retval;
  int len;

  // only run query from rank 0
  if (rank == 0) {

    // fake query
    strcpy(val,(const char*)"4.25");
  }
  MPI_Bcast(val, 512, MPI_CHAR, 0, comm);
  len = strlen(val) + 1;
  retval = new char[len];
  strcpy(retval,val);
    
  return retval;
}
