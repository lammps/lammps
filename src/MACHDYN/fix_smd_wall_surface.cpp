// clang-format off
/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 https://www.lammps.org/, Sandia National Laboratories
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

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "text_file_reader.h"

#include <cstring>
#include <Eigen/Eigen>

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

        //atom->add_callback(Atom::GROW);
        //atom->add_callback(Atom::RESTART);

        if (narg != 6)
                error->all(FLERR, "Illegal number of arguments for fix smd/wall_surface");

        filename = strdup(arg[3]);
        wall_particle_type = utils::inumeric(FLERR,arg[4],false,lmp);
        wall_molecule_id = utils::inumeric(FLERR,arg[5],false,lmp);
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
        filename = nullptr;
        // unregister this fix so atom class doesn't invoke it any more

        //atom->delete_callback(id,Atom::GROW);
        //atom->delete_callback(id,Atom::RESTART);
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
 size of atom nlocal's restart data
 ------------------------------------------------------------------------- */

void FixSMDWallSurface::read_triangles(int pass) {

  double coord[3];
  int nlocal_previous = atom->nlocal;
  int ilocal = nlocal_previous;
  bigint natoms_previous = atom->natoms;
  Vector3d *vert;
  vert = new Vector3d[3];
  Vector3d normal, center;

  FILE *fp = fopen(filename, "r");
  if (fp == nullptr)
    error->one(FLERR, "Cannot open file {}: {}", filename, utils::getsyserror());

  if (comm->me == 0) {
    utils::logmesg(lmp, "\n>>========>>========>>========>>========>>========>>========\n");
    if (pass == 0)
      utils::logmesg(lmp, "  scanning triangle pairs ...\n");
    else
      utils::logmesg(lmp, "  reading triangle pairs ...\n");
  }

  TextFileReader reader(fp, "triangles");
  try {
    char *line = reader.next_line();
    if (!line || !utils::strmatch(line, "^solid"))
      throw TokenizerException("Invalid triangles file format","");

    if (comm->me == 0)
      utils::logmesg(lmp, "  reading STL object '{}' from {}\n", utils::trim(line+6), filename);

    while((line = reader.next_line())) {

      // next line is facet line with 5 words
      auto values = utils::split_words(line);
      // otherwise stop reading
      if ((values.size() != 5) || !utils::strmatch(values[0],"^facet")) break;

      normal << utils::numeric(FLERR, values[2], false, lmp),
        utils::numeric(FLERR, values[3], false, lmp),
        utils::numeric(FLERR, values[4], false, lmp);

      line = reader.next_line(2);
      if (!line || !utils::strmatch(line, "^ *outer *loop"))
        throw TokenizerException("Error reading outer loop","");

      for (int k = 0; k < 3; ++k) {
        line = reader.next_line(4);
        values = utils::split_words(line);
        if ((values.size() != 4) || !utils::strmatch(values[0],"^vertex"))
          throw TokenizerException("Error reading vertex","");

        vert[k] << utils::numeric(FLERR, values[1], false, lmp),
          utils::numeric(FLERR, values[2], false, lmp),
          utils::numeric(FLERR, values[3], false, lmp);
      }

      line = reader.next_line(1);
      if (!line || !utils::strmatch(line, "^ *endloop"))
        throw TokenizerException("Error reading endloop","");
      line = reader.next_line(1);
      if (!line || !utils::strmatch(line, "^ *endfacet"))
        throw TokenizerException("Error reading endfacet","");

      // now we have a normal and three vertices ... proceed with adding triangle

      center = (vert[0] + vert[1] + vert[2]) / 3.0;

      double r1 = (center - vert[0]).norm();
      double r2 = (center - vert[1]).norm();
      double r3 = (center - vert[2]).norm();
      double r = MAX(MAX(r1, r2), r3);

      /*
       * if atom/molecule is in my subbox, create it
       * ... use x0 to hold triangle normal.
       * ... use smd_data_9 to hold the three vertices
       * ... use x to hold triangle center
       * ... radius is the mmaximal distance from triangle center to all vertices
       */

      if (center(0) >= sublo[0] && center(0) < subhi[0] && center(1) >= sublo[1] && center(1) < subhi[1] && center(2) >= sublo[2]
          && center(2) < subhi[2]) {

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
  } catch (std::exception &e) {
    error->all(FLERR, "Error reading triangles from file {}: {}", filename, e.what());
  }

  // set new total # of atoms and error check

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT)
    error->all(FLERR, "Too many total atoms");

  // add IDs for newly created atoms
  // check that atom IDs are valid

  if (atom->tag_enable) atom->tag_extend();
  atom->tag_check();

  // create global mapping of atoms
  // zero nghost in case are adding new atoms to existing atoms

  if (atom->map_style != Atom::MAP_NONE) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // print status
  if (comm->me == 0) {
    utils::logmesg(lmp, "... fix smd/wall_surface finished reading triangulated surface\n");
    utils::logmesg(lmp, "fix smd/wall_surface created {} atoms\n", atom->natoms - natoms_previous);
    utils::logmesg(lmp, ">>========>>========>>========>>========>>========>>========\n");
  }

  delete[] vert;
  fclose(fp);
}
