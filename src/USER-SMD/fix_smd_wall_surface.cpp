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
   Contributing authors: Mike Parks (SNL), Ezwanur Rahman, J.T. Foster (UTSA)
------------------------------------------------------------------------- */

#include "fix_smd_wall_surface.h"
#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include <Eigen/Eigen>
#include "atom.h"
#include "domain.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "atom_vec.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace Eigen;
using namespace std;
#define DELTA 16384
#define EPSILON 1.0e-6

/* ---------------------------------------------------------------------- */

FixSMDWallSurface::FixSMDWallSurface(LAMMPS *lmp, int narg, char **arg) :
                Fix(lmp, narg, arg) {

        restart_global = 0;
        restart_peratom = 0;
        first = 1;

        //atom->add_callback(0);
        //atom->add_callback(1);

        if (narg != 6)
                error->all(FLERR, "Illegal number of arguments for fix smd/wall_surface");

        filename = strdup(arg[3]);
        wall_particle_type = force->inumeric(FLERR, arg[4]);
        wall_molecule_id = force->inumeric(FLERR, arg[5]);
        if (wall_molecule_id < 65535) {
                error->one(FLERR, "wall molcule id must be >= 65535\n");
        }

        if (comm->me == 0) {
                printf("\n>>========>>========>>========>>========>>========>>========>>========>>========\n");
                printf("fix smd/wall_surface reads trianglulated surface from file: %s\n", filename);
                printf("fix smd/wall_surface has particle type %d \n", wall_particle_type);
                printf("fix smd/wall_surface has molecule id %d \n", wall_molecule_id);
                printf(">>========>>========>>========>>========>>========>>========>>========>>========\n");
        }
}

/* ---------------------------------------------------------------------- */

FixSMDWallSurface::~FixSMDWallSurface() {
        free(filename);
        filename = NULL;
        // unregister this fix so atom class doesn't invoke it any more

        //atom->delete_callback(id, 0);
        //atom->delete_callback(id, 1);
}

/* ---------------------------------------------------------------------- */

int FixSMDWallSurface::setmask() {
        int mask = 0;
        return mask;
}

/* ---------------------------------------------------------------------- */

void FixSMDWallSurface::init() {
        if (!first)
                return;
}

/* ----------------------------------------------------------------------
 For minimization: setup as with dynamics
 ------------------------------------------------------------------------- */

void FixSMDWallSurface::min_setup(int vflag) {
        setup(vflag);
}

/* ----------------------------------------------------------------------
 create initial list of neighbor partners via call to neighbor->build()
 must be done in setup (not init) since fix init comes before neigh init
 ------------------------------------------------------------------------- */

void FixSMDWallSurface::setup(int /*vflag*/) {

        if (!first)
                return;
        first = 0;

        // set bounds for my proc
        // if periodic and I am lo/hi proc, adjust bounds by EPSILON
        // insures all data atoms will be owned even with round-off

        int triclinic = domain->triclinic;

        double epsilon[3];
        if (triclinic)
                epsilon[0] = epsilon[1] = epsilon[2] = EPSILON;
        else {
                epsilon[0] = domain->prd[0] * EPSILON;
                epsilon[1] = domain->prd[1] * EPSILON;
                epsilon[2] = domain->prd[2] * EPSILON;
        }

        if (triclinic == 0) {
                sublo[0] = domain->sublo[0];
                subhi[0] = domain->subhi[0];
                sublo[1] = domain->sublo[1];
                subhi[1] = domain->subhi[1];
                sublo[2] = domain->sublo[2];
                subhi[2] = domain->subhi[2];
        } else {
                sublo[0] = domain->sublo_lamda[0];
                subhi[0] = domain->subhi_lamda[0];
                sublo[1] = domain->sublo_lamda[1];
                subhi[1] = domain->subhi_lamda[1];
                sublo[2] = domain->sublo_lamda[2];
                subhi[2] = domain->subhi_lamda[2];
        }

        if (comm->layout != Comm::LAYOUT_TILED) {
                if (domain->xperiodic) {
                        if (comm->myloc[0] == 0)
                                sublo[0] -= epsilon[0];
                        if (comm->myloc[0] == comm->procgrid[0] - 1)
                                subhi[0] += epsilon[0];
                }
                if (domain->yperiodic) {
                        if (comm->myloc[1] == 0)
                                sublo[1] -= epsilon[1];
                        if (comm->myloc[1] == comm->procgrid[1] - 1)
                                subhi[1] += epsilon[1];
                }
                if (domain->zperiodic) {
                        if (comm->myloc[2] == 0)
                                sublo[2] -= epsilon[2];
                        if (comm->myloc[2] == comm->procgrid[2] - 1)
                                subhi[2] += epsilon[2];
                }

        } else {
                if (domain->xperiodic) {
                        if (comm->mysplit[0][0] == 0.0)
                                sublo[0] -= epsilon[0];
                        if (comm->mysplit[0][1] == 1.0)
                                subhi[0] += epsilon[0];
                }
                if (domain->yperiodic) {
                        if (comm->mysplit[1][0] == 0.0)
                                sublo[1] -= epsilon[1];
                        if (comm->mysplit[1][1] == 1.0)
                                subhi[1] += epsilon[1];
                }
                if (domain->zperiodic) {
                        if (comm->mysplit[2][0] == 0.0)
                                sublo[2] -= epsilon[2];
                        if (comm->mysplit[2][1] == 1.0)
                                subhi[2] += epsilon[2];
                }
        }

        read_triangles(0);
}

/* ----------------------------------------------------------------------
 function to determine number of values in a text line
 ------------------------------------------------------------------------- */

int FixSMDWallSurface::count_words(const char *line) {
        int n = strlen(line) + 1;
        char *copy;
        memory->create(copy, n, "atom:copy");
        strcpy(copy, line);

        char *ptr;
        if ((ptr = strchr(copy, '#')))
                *ptr = '\0';

        if (strtok(copy, " \t\n\r\f") == NULL) {
                memory->destroy(copy);
                return 0;
        }
        n = 1;
        while (strtok(NULL, " \t\n\r\f"))
                n++;

        memory->destroy(copy);
        return n;
}

/* ----------------------------------------------------------------------
 size of atom nlocal's restart data
 ------------------------------------------------------------------------- */

void FixSMDWallSurface::read_triangles(int pass) {

  double coord[3];

  int nlocal_previous = atom->nlocal;
  int ilocal = nlocal_previous;
  int m;
  int me;

  bigint natoms_previous = atom->natoms;
  Vector3d *vert;
  vert = new Vector3d[3];
  Vector3d normal, center;

  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    char str[128];
    snprintf(str,128, "Cannot open file %s", filename);
    error->one(FLERR, str);
  }

  MPI_Comm_rank(world, &me);
  if (me == 0) {
    if (screen) {
      if (pass == 0) {
        printf("\n>>========>>========>>========>>========>>========>>========>>========>>========\n");
        fprintf(screen, "  scanning triangle pairs ...\n");
      } else {
        fprintf(screen, "  reading triangle pairs ...\n");
      }
    }
    if (logfile) {
      if (pass == 0) {
        fprintf(logfile, "  scanning triangle pairs ...\n");
      } else {
        fprintf(logfile, "  reading triangle pairs ...\n");
      }
    }
  }

  char line[256];
  char *retpointer;
  char **values;
  int nwords;

  // read STL solid name
  retpointer = fgets(line, sizeof(line), fp);
  if (retpointer == NULL) {
    error->one(FLERR,"error reading number of triangle pairs");
  }

  nwords = count_words(line);
  if (nwords < 1) {
    error->one(FLERR,"first line of file is incorrect");
  }

//      values = new char*[nwords];
//      values[0] = strtok(line, " \t\n\r\f");
//      if (values[0] == NULL)
//              error->all(FLERR, "Incorrect atom format in data file");
//      for (m = 1; m < nwords; m++) {
//              values[m] = strtok(NULL, " \t\n\r\f");
//              if (values[m] == NULL)
//                      error->all(FLERR, "Incorrect atom format in data file");
//      }
//      delete[] values;
//
//      if (comm->me == 0) {
//              cout << "STL file contains solid body with name: " << values[1] << endl;
//      }

  // iterate over STL facets util end of body is reached

  while (fgets(line, sizeof(line), fp)) { // read a line, should be the facet line

    // evaluate facet line
    nwords = count_words(line);
    if (nwords != 5) {
      //sprintf(str, "found end solid line");
      //error->message(FLERR, str);
      break;
    } else {
      // should be facet line
    }

    values = new char*[nwords];
    values[0] = strtok(line, " \t\n\r\f");
    if (values[0] == NULL)
      error->all(FLERR, "Incorrect atom format in data file");
    for (m = 1; m < nwords; m++) {
      values[m] = strtok(NULL, " \t\n\r\f");
      if (values[m] == NULL)
        error->all(FLERR, "Incorrect atom format in data file");
    }

    normal << force->numeric(FLERR, values[2]), force->numeric(FLERR, values[3]), force->numeric(FLERR, values[4]);
    //cout << "normal is " << normal << endl;

    delete[] values;

    // read outer loop line
    retpointer = fgets(line, sizeof(line), fp);
    if (retpointer == NULL) {
      error->one(FLERR, "error reading outer loop");
    }

    nwords = count_words(line);
    if (nwords != 2) {
      error->one(FLERR,"error reading outer loop");
    }

    // read vertex lines

    for (int k = 0; k < 3; k++) {
      retpointer = fgets(line, sizeof(line), fp);
      if (retpointer == NULL) {
        error->one(FLERR,"error reading vertex line");
      }

      nwords = count_words(line);
      if (nwords != 4) {
        error->one(FLERR,"error reading vertex line");
      }

      values = new char*[nwords];
      values[0] = strtok(line, " \t\n\r\f");
      if (values[0] == NULL)
        error->all(FLERR,"Incorrect vertex line");
      for (m = 1; m < nwords; m++) {
        values[m] = strtok(NULL, " \t\n\r\f");
        if (values[m] == NULL)
          error->all(FLERR, "Incorrect vertex line");
      }

      vert[k] << force->numeric(FLERR, values[1]), force->numeric(FLERR, values[2]), force->numeric(FLERR, values[3]);
      //cout << "vertex is " << vert[k] << endl;
      //printf("%s %s %s\n", values[1], values[2], values[3]);
      delete[] values;
      //exit(1);

    }

    // read end loop line
    retpointer = fgets(line, sizeof(line), fp);
    if (retpointer == NULL) {
      error->one(FLERR, "error reading endloop");
    }

    nwords = count_words(line);
    if (nwords != 1) {
      error->one(FLERR,"error reading endloop");
    }

    // read end facet line
    retpointer = fgets(line, sizeof(line), fp);
    if (retpointer == NULL) {
      error->one(FLERR,"error reading endfacet");
    }

    nwords = count_words(line);
    if (nwords != 1) {
      error->one(FLERR,"error reading endfacet");
    }

    // now we have a normal and three vertices ... proceed with adding triangle

    center = (vert[0] + vert[1] + vert[2]) / 3.0;

    //      cout << "center is " << center << endl;

    double r1 = (center - vert[0]).norm();
    double r2 = (center - vert[1]).norm();
    double r3 = (center - vert[2]).norm();
    double r = MAX(r1, r2);
    r = MAX(r, r3);

    /*
     * if atom/molecule is in my subbox, create it
     * ... use x0 to hold triangle normal.
     * ... use smd_data_9 to hold the three vertices
     * ... use x to hold triangle center
     * ... radius is the mmaximal distance from triangle center to all vertices
     */

    //      printf("coord: %f %f %f\n", coord[0], coord[1], coord[2]);
    //      printf("sublo: %f %f %f\n", sublo[0], sublo[1], sublo[2]);
    //      printf("subhi: %f %f %f\n", subhi[0], subhi[1], subhi[2]);
    //printf("ilocal = %d\n", ilocal);
    if (center(0) >= sublo[0] && center(0) < subhi[0] && center(1) >= sublo[1] && center(1) < subhi[1] && center(2) >= sublo[2]
        && center(2) < subhi[2]) {
      //printf("******* KERATIN nlocal=%d ***\n", nlocal);
      coord[0] = center(0);
      coord[1] = center(1);
      coord[2] = center(2);
      atom->avec->create_atom(wall_particle_type, coord);

      /*
       * need to initialize pointers to atom vec arrays here, because they could have changed
       * due to calling grow() in create_atoms() above;
       */

      tagint *mol = atom->molecule;
      int *type = atom->type;
      double *radius = atom->radius;
      double *contact_radius = atom->contact_radius;
      double **smd_data_9 = atom->smd_data_9;
      double **x0 = atom->x0;

      radius[ilocal] = r; //ilocal;
      contact_radius[ilocal] = r; //ilocal;
      mol[ilocal] = wall_molecule_id;
      type[ilocal] = wall_particle_type;
      x0[ilocal][0] = normal(0);
      x0[ilocal][1] = normal(1);
      x0[ilocal][2] = normal(2);
      smd_data_9[ilocal][0] = vert[0](0);
      smd_data_9[ilocal][1] = vert[0](1);
      smd_data_9[ilocal][2] = vert[0](2);
      smd_data_9[ilocal][3] = vert[1](0);
      smd_data_9[ilocal][4] = vert[1](1);
      smd_data_9[ilocal][5] = vert[1](2);
      smd_data_9[ilocal][6] = vert[2](0);
      smd_data_9[ilocal][7] = vert[2](1);
      smd_data_9[ilocal][8] = vert[2](2);

      ilocal++;
    }

  }

// set new total # of atoms and error check

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT)
    error->all(FLERR, "Too many total atoms");

// add IDs for newly created atoms
// check that atom IDs are valid

  if (atom->tag_enable)
    atom->tag_extend();
  atom->tag_check();

// create global mapping of atoms
// zero nghost in case are adding new atoms to existing atoms

  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

// print status
  if (comm->me == 0) {
    if (screen) {
      printf("... fix smd/wall_surface finished reading triangulated surface\n");
      fprintf(screen, "fix smd/wall_surface created " BIGINT_FORMAT " atoms\n", atom->natoms - natoms_previous);
      printf(">>========>>========>>========>>========>>========>>========>>========>>========\n");
    }
    if (logfile) {
      fprintf(logfile, "... fix smd/wall_surface finished reading triangulated surface\n");
      fprintf(logfile, "fix smd/wall_surface created " BIGINT_FORMAT " atoms\n", atom->natoms - natoms_previous);
      fprintf(logfile, ">>========>>========>>========>>========>>========>>========>>========>>========\n");
    }
  }

  delete[] vert;
  fclose(fp);
}

