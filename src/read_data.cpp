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

// lmptype.h must be first b/c this file uses MAXBIGINT and includes mpi.h
// due to OpenMPI bug which sets INT64_MAX via its mpi.h
//   before lmptype.h can set flags to insure it is done correctly

#include "read_data.h"
#include "lmptype.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "force.h"
#include "molecule.h"
#include "group.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "special.h"
#include "irregular.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

#define MAXLINE 256
#define LB_FACTOR 1.1
#define CHUNK 1024
#define DELTA 4            // must be 2 or larger
#define MAXBODY 32         // max # of lines in one body

                           // customize for new sections
#define NSECTIONS 25       // change when add to header::section_keywords

enum{NONE,APPEND,VALUE,MERGE};

// pair style suffixes to ignore
// when matching Pair Coeffs comment to currently-defined pair style

const char *suffixes[] = {"/cuda","/gpu","/opt","/omp","/kk",
                          "/coul/cut","/coul/long","/coul/msm",
                          "/coul/dsf","/coul/debye","/coul/charmm",
                          NULL};

/* ---------------------------------------------------------------------- */

ReadData::ReadData(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  line = new char[MAXLINE];
  copy = new char[MAXLINE];
  keyword = new char[MAXLINE];
  style = new char[MAXLINE];
  buffer = new char[CHUNK*MAXLINE];
  narg = maxarg = 0;
  arg = NULL;
  fp = NULL;

  // customize for new sections
  // pointers to atom styles that store extra info

  nellipsoids = 0;
  avec_ellipsoid = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  nlines = 0;
  avec_line = (AtomVecLine *) atom->style_match("line");
  ntris = 0;
  avec_tri = (AtomVecTri *) atom->style_match("tri");
  nbodies = 0;
  avec_body = (AtomVecBody *) atom->style_match("body");
}

/* ---------------------------------------------------------------------- */

ReadData::~ReadData()
{
  delete [] line;
  delete [] copy;
  delete [] keyword;
  delete [] style;
  delete [] buffer;
  memory->sfree(arg);

  for (int i = 0; i < nfix; i++) {
    delete [] fix_header[i];
    delete [] fix_section[i];
  }
  memory->destroy(fix_index);
  memory->sfree(fix_header);
  memory->sfree(fix_section);
}

/* ---------------------------------------------------------------------- */

void ReadData::command(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal read_data command");

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  // optional args

  addflag = NONE;
  coeffflag = 1;
  id_offset = mol_offset = 0;
  offsetflag = shiftflag = 0;
  toffset = boffset = aoffset = doffset = ioffset = 0;
  shift[0] = shift[1] = shift[2] = 0.0;
  extra_atom_types = extra_bond_types = extra_angle_types =
    extra_dihedral_types = extra_improper_types = 0;

  groupbit = 0;

  nfix = 0;
  fix_index = NULL;
  fix_header = NULL;
  fix_section = NULL;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"add") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_data command");
      if (strcmp(arg[iarg+1],"append") == 0) addflag = APPEND;
      else if (strcmp(arg[iarg+1],"merge") == 0) addflag = MERGE;
      else {
        if (atom->molecule_flag && (iarg+3 > narg))
          error->all(FLERR,"Illegal read_data command");
        addflag = VALUE;
        bigint offset = force->bnumeric(FLERR,arg[iarg+1]);
        if (offset > MAXTAGINT)
          error->all(FLERR,"Read data add atomID offset is too big");
        id_offset = offset;

        if (atom->molecule_flag) {
          offset = force->bnumeric(FLERR,arg[iarg+2]);
          if (offset > MAXTAGINT)
            error->all(FLERR,"Read data add molID offset is too big");
          mol_offset = offset;
          iarg++;
        }
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"offset") == 0) {
      if (iarg+6 > narg) error->all(FLERR,"Illegal read_data command");
      offsetflag = 1;
      toffset = force->inumeric(FLERR,arg[iarg+1]);
      boffset = force->inumeric(FLERR,arg[iarg+2]);
      aoffset = force->inumeric(FLERR,arg[iarg+3]);
      doffset = force->inumeric(FLERR,arg[iarg+4]);
      ioffset = force->inumeric(FLERR,arg[iarg+5]);
      if (toffset < 0 || boffset < 0 || aoffset < 0 ||
          doffset < 0 || ioffset < 0)
        error->all(FLERR,"Illegal read_data command");
      iarg += 6;
    } else if (strcmp(arg[iarg],"shift") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal read_data command");
      shiftflag = 1;
      shift[0] = force->numeric(FLERR,arg[iarg+1]);
      shift[1] = force->numeric(FLERR,arg[iarg+2]);
      shift[2] = force->numeric(FLERR,arg[iarg+3]);
      if (domain->dimension == 2 && shift[2] != 0.0)
        error->all(FLERR,"Non-zero read_data shift z value for 2d simulation");
      iarg += 4;
    } else if (strcmp(arg[iarg],"nocoeff") == 0) {
      coeffflag = 0;
      iarg ++;
    } else if (strcmp(arg[iarg],"extra/atom/types") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_data command");
      extra_atom_types = force->inumeric(FLERR,arg[iarg+1]);
      if (extra_atom_types < 0) error->all(FLERR,"Illegal read_data command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra/bond/types") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_data command");
      if (!atom->avec->bonds_allow)
        error->all(FLERR,"No bonds allowed with this atom style");
      extra_bond_types = force->inumeric(FLERR,arg[iarg+1]);
      if (extra_bond_types < 0) error->all(FLERR,"Illegal read_data command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra/angle/types") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_data command");
      if (!atom->avec->angles_allow)
        error->all(FLERR,"No angles allowed with this atom style");
      extra_angle_types = force->inumeric(FLERR,arg[iarg+1]);
      if (extra_angle_types < 0) error->all(FLERR,"Illegal read_data command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra/dihedral/types") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_data command");
      if (!atom->avec->dihedrals_allow)
        error->all(FLERR,"No dihedrals allowed with this atom style");
      extra_dihedral_types = force->inumeric(FLERR,arg[iarg+1]);
      if (extra_dihedral_types < 0)
        error->all(FLERR,"Illegal read_data command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra/improper/types") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_data command");
      if (!atom->avec->impropers_allow)
        error->all(FLERR,"No impropers allowed with this atom style");
      extra_improper_types = force->inumeric(FLERR,arg[iarg+1]);
      if (extra_improper_types < 0)
        error->all(FLERR,"Illegal read_data command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra/bond/per/atom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_data command");
      if (! atom->molecular)
        error->all(FLERR,"No bonds allowed with this atom style");
      atom->extra_bond_per_atom = force->inumeric(FLERR,arg[iarg+1]);
      if (atom->extra_bond_per_atom < 0)
        error->all(FLERR,"Illegal read_data command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra/angle/per/atom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_data command");
      if (! atom->molecular)
        error->all(FLERR,"No angles allowed with this atom style");
      atom->extra_angle_per_atom = force->inumeric(FLERR,arg[iarg+1]);
      if (atom->extra_angle_per_atom < 0)
        error->all(FLERR,"Illegal read_data command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra/dihedral/per/atom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_data command");
      if (! atom->molecular)
        error->all(FLERR,"No dihedrals allowed with this atom style");
      atom->extra_dihedral_per_atom = force->inumeric(FLERR,arg[iarg+1]);
      if (atom->extra_dihedral_per_atom < 0)
        error->all(FLERR,"Illegal read_data command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra/improper/per/atom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_data command");
      if (! atom->molecular)
        error->all(FLERR,"No impropers allowed with this atom style");
      atom->extra_improper_per_atom = force->inumeric(FLERR,arg[iarg+1]);
      if (atom->extra_improper_per_atom < 0)
        error->all(FLERR,"Illegal read_data command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra/special/per/atom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_data command");
      if (! atom->molecular)
        error->all(FLERR,"No bonded interactions allowed with this atom style");
      force->special_extra = force->inumeric(FLERR,arg[iarg+1]);
      if (force->special_extra < 0)
        error->all(FLERR,"Illegal read_data command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_data command");
      int igroup = group->find_or_create(arg[iarg+1]);
      groupbit = group->bitmask[igroup];
      iarg += 2;
    } else if (strcmp(arg[iarg],"fix") == 0) {
      if (iarg+4 > narg)
        error->all(FLERR,"Illegal read_data command");
      memory->grow(fix_index,nfix+1,"read_data:fix_index");
      fix_header = (char **)
        memory->srealloc(fix_header,(nfix+1)*sizeof(char *),
                         "read_data:fix_header");
      fix_section = (char **)
        memory->srealloc(fix_section,(nfix+1)*sizeof(char *),
                         "read_data:fix_section");
      fix_index[nfix] = modify->find_fix(arg[iarg+1]);
      if (fix_index[nfix] < 0)
        error->all(FLERR,"Fix ID for read_data does not exist");
      if (strcmp(arg[iarg+2],"NULL") == 0) fix_header[nfix] = NULL;
      else {
        int n = strlen(arg[iarg+2]) + 1;
        fix_header[nfix] = new char[n];
        strcpy(fix_header[nfix],arg[iarg+2]);
      }
      int n = strlen(arg[iarg+3]) + 1;
      fix_section[nfix] = new char[n];
      strcpy(fix_section[nfix],arg[iarg+3]);
      nfix++;
      iarg += 4;

    } else error->all(FLERR,"Illegal read_data command");
  }

  // error checks

  if (domain->dimension == 2 && domain->zperiodic == 0)
    error->all(FLERR,"Cannot run 2d simulation with nonperiodic Z dimension");
  if (domain->box_exist && !addflag)
    error->all(FLERR,"Cannot read_data without add keyword "
               "after simulation box is defined");
  if (!domain->box_exist && addflag)
    error->all(FLERR,"Cannot use read_data add before "
               "simulation box is defined");
  if (offsetflag && addflag == NONE)
    error->all(FLERR,"Cannot use read_data offset without add flag");
  if (shiftflag && addflag == NONE)
    error->all(FLERR,"Cannot use read_data shift without add flag");
  if (addflag != NONE &&
      (extra_atom_types || extra_bond_types || extra_angle_types ||
       extra_dihedral_types || extra_improper_types))
    error->all(FLERR,"Cannot use read_data extra with add flag");

  // first time system initialization

  if (addflag == NONE) {
    domain->box_exist = 1;
    update->ntimestep = 0;
  }

  // compute atomID and optionally moleculeID offset for addflag = APPEND

  if (addflag == APPEND) {
    tagint *tag = atom->tag;
    tagint *molecule = atom->molecule;
    int nlocal = atom->nlocal;
    tagint maxid = 0, maxmol = 0;
    for (int i = 0; i < nlocal; i++) maxid = MAX(maxid,tag[i]);
    if (atom->molecule_flag)
      for (int i = 0; i < nlocal; i++) maxmol = MAX(maxmol,molecule[i]);
    MPI_Allreduce(&maxid,&id_offset,1,MPI_LMP_TAGINT,MPI_MAX,world);
    MPI_Allreduce(&maxmol,&mol_offset,1,MPI_LMP_TAGINT,MPI_MAX,world);
  }

  // set up pointer to hold original styles while we replace them with "zero"

  Pair *saved_pair = NULL;
  Bond *saved_bond = NULL;
  Angle *saved_angle = NULL;
  Dihedral *saved_dihedral = NULL;
  Improper *saved_improper = NULL;
  KSpace *saved_kspace = NULL;

  if (coeffflag == 0) {
    char *coeffs[2];
    coeffs[0] = (char *) "10.0";
    coeffs[1] = (char *) "nocoeff";

    saved_pair = force->pair;
    force->pair = NULL;
    force->create_pair("zero",0);
    if (force->pair) force->pair->settings(2,coeffs);

    coeffs[0] = coeffs[1];
    saved_bond = force->bond;
    force->bond = NULL;
    force->create_bond("zero",0);
    if (force->bond) force->bond->settings(1,coeffs);

    saved_angle = force->angle;
    force->angle = NULL;
    force->create_angle("zero",0);
    if (force->angle) force->angle->settings(1,coeffs);

    saved_dihedral = force->dihedral;
    force->dihedral = NULL;
    force->create_dihedral("zero",0);
    if (force->dihedral) force->dihedral->settings(1,coeffs);

    saved_improper = force->improper;
    force->improper = NULL;
    force->create_improper("zero",0);
    if (force->improper) force->improper->settings(1,coeffs);

    saved_kspace = force->kspace;
    force->kspace = NULL;
  }

  // -----------------------------------------------------------------

  // perform 1-pass read if no molecular topology in file
  // perform 2-pass read if molecular topology,
  //   first pass calculates max topology/atom

  // flags for this data file

  int atomflag,topoflag;
  int bondflag,angleflag,dihedralflag,improperflag;
  int ellipsoidflag,lineflag,triflag,bodyflag;

  atomflag = topoflag = 0;
  bondflag = angleflag = dihedralflag = improperflag = 0;
  ellipsoidflag = lineflag = triflag = bodyflag = 0;

  // values in this data file

  natoms = 0;
  ntypes = 0;
  nbonds = nangles = ndihedrals = nimpropers = 0;
  nbondtypes = nangletypes = ndihedraltypes = nimpropertypes = 0;

  boxlo[0] = boxlo[1] = boxlo[2] = -0.5;
  boxhi[0] = boxhi[1] = boxhi[2] = 0.5;
  triclinic = 0;
  keyword[0] = '\0';

  nlocal_previous = atom->nlocal;
  int firstpass = 1;

  while (1) {

    // open file on proc 0

    if (me == 0) {
      if (firstpass && screen) fprintf(screen,"Reading data file ...\n");
      open(arg[0]);
    } else fp = NULL;

    // read header info

    header(firstpass);

    // problem setup using info from header
    // only done once, if firstpass and first data file
    // apply extra settings before grow(), even if no topology in file
    // deallocate() insures new settings are used for topology arrays
    // if per-atom topology is in file, another grow() is done below

    if (firstpass && addflag == NONE) {
      atom->bond_per_atom = atom->extra_bond_per_atom;
      atom->angle_per_atom = atom->extra_angle_per_atom;
      atom->dihedral_per_atom = atom->extra_dihedral_per_atom;
      atom->improper_per_atom = atom->extra_improper_per_atom;

      int n;
      if (comm->nprocs == 1) n = static_cast<int> (atom->natoms);
      else n = static_cast<int> (LB_FACTOR * atom->natoms / comm->nprocs);

      atom->allocate_type_arrays();
      atom->deallocate_topology();
      atom->avec->grow(n);

      domain->boxlo[0] = boxlo[0]; domain->boxhi[0] = boxhi[0];
      domain->boxlo[1] = boxlo[1]; domain->boxhi[1] = boxhi[1];
      domain->boxlo[2] = boxlo[2]; domain->boxhi[2] = boxhi[2];

      if (triclinic) {
        domain->triclinic = 1;
        domain->xy = xy; domain->xz = xz; domain->yz = yz;
      }

      domain->print_box("  ");
      domain->set_initial_box();
      domain->set_global_box();
      comm->set_proc_grid();
      domain->set_local_box();
    }

    // change simulation box to be union of existing box and new box + shift
    // only done if firstpass and not first data file

    if (firstpass && addflag != NONE) {
      domain->boxlo[0] = MIN(domain->boxlo[0],boxlo[0]+shift[0]);
      domain->boxhi[0] = MAX(domain->boxhi[0],boxhi[0]+shift[0]);
      domain->boxlo[1] = MIN(domain->boxlo[1],boxlo[1]+shift[1]);
      domain->boxhi[1] = MAX(domain->boxhi[1],boxhi[1]+shift[1]);
      domain->boxlo[2] = MIN(domain->boxlo[2],boxlo[2]+shift[2]);
      domain->boxhi[2] = MAX(domain->boxhi[2],boxhi[2]+shift[2]);

      // NOTE: not sure what to do about tilt value in subsequent data files
      //if (triclinic) {
      //  domain->xy = xy; domain->xz = xz; domain->yz = yz;
      // }

      domain->print_box("  ");
      domain->set_initial_box();
      domain->set_global_box();
      comm->set_proc_grid();
      domain->set_local_box();
    }

    // customize for new sections
    // read rest of file in free format

    while (strlen(keyword)) {

      // if special fix matches, it processes section

      if (nfix) {
        int i;
        for (i = 0; i < nfix; i++)
          if (strcmp(keyword,fix_section[i]) == 0) {
            if (firstpass) fix(fix_index[i],keyword);
            else skip_lines(modify->fix[fix_index[i]]->
                            read_data_skip_lines(keyword));
            parse_keyword(0);
            break;
          }
        if (i < nfix) continue;
      }

      if (strcmp(keyword,"Atoms") == 0) {
        atomflag = 1;
        if (firstpass) {
          if (me == 0 && !style_match(style,atom->atom_style))
            error->warning(FLERR,"Atom style in data file differs "
                           "from currently defined atom style");
          atoms();
        } else skip_lines(natoms);
      } else if (strcmp(keyword,"Velocities") == 0) {
        if (atomflag == 0)
          error->all(FLERR,"Must read Atoms before Velocities");
        if (firstpass) velocities();
        else skip_lines(natoms);

      } else if (strcmp(keyword,"Bonds") == 0) {
        topoflag = bondflag = 1;
        if (nbonds == 0)
          error->all(FLERR,"Invalid data file section: Bonds");
        if (atomflag == 0) error->all(FLERR,"Must read Atoms before Bonds");
        bonds(firstpass);
      } else if (strcmp(keyword,"Angles") == 0) {
        topoflag = angleflag = 1;
        if (nangles == 0)
          error->all(FLERR,"Invalid data file section: Angles");
        if (atomflag == 0) error->all(FLERR,"Must read Atoms before Angles");
        angles(firstpass);
      } else if (strcmp(keyword,"Dihedrals") == 0) {
        topoflag = dihedralflag = 1;
        if (ndihedrals == 0)
          error->all(FLERR,"Invalid data file section: Dihedrals");
        if (atomflag == 0) error->all(FLERR,"Must read Atoms before Dihedrals");
        dihedrals(firstpass);
      } else if (strcmp(keyword,"Impropers") == 0) {
        topoflag = improperflag = 1;
        if (nimpropers == 0)
          error->all(FLERR,"Invalid data file section: Impropers");
        if (atomflag == 0) error->all(FLERR,"Must read Atoms before Impropers");
        impropers(firstpass);

      } else if (strcmp(keyword,"Ellipsoids") == 0) {
        ellipsoidflag = 1;
        if (!avec_ellipsoid)
          error->all(FLERR,"Invalid data file section: Ellipsoids");
        if (atomflag == 0)
          error->all(FLERR,"Must read Atoms before Ellipsoids");
        if (firstpass)
          bonus(nellipsoids,(AtomVec *) avec_ellipsoid,"ellipsoids");
        else skip_lines(nellipsoids);
      } else if (strcmp(keyword,"Lines") == 0) {
        lineflag = 1;
        if (!avec_line)
          error->all(FLERR,"Invalid data file section: Lines");
        if (atomflag == 0) error->all(FLERR,"Must read Atoms before Lines");
        if (firstpass) bonus(nlines,(AtomVec *) avec_line,"lines");
        else skip_lines(nlines);
      } else if (strcmp(keyword,"Triangles") == 0) {
        triflag = 1;
        if (!avec_tri)
          error->all(FLERR,"Invalid data file section: Triangles");
        if (atomflag == 0) error->all(FLERR,"Must read Atoms before Triangles");
        if (firstpass) bonus(ntris,(AtomVec *) avec_tri,"triangles");
        else skip_lines(ntris);
      } else if (strcmp(keyword,"Bodies") == 0) {
        bodyflag = 1;
        if (!avec_body)
          error->all(FLERR,"Invalid data file section: Bodies");
        if (atomflag == 0) error->all(FLERR,"Must read Atoms before Bodies");
        bodies(firstpass);

      } else if (strcmp(keyword,"Masses") == 0) {
        if (firstpass) mass();
        else skip_lines(ntypes);
      } else if (strcmp(keyword,"Pair Coeffs") == 0) {
        if (force->pair == NULL)
          error->all(FLERR,"Must define pair_style before Pair Coeffs");
        if (firstpass) {
          if (me == 0 && !style_match(style,force->pair_style))
            error->warning(FLERR,"Pair style in data file differs "
                           "from currently defined pair style");
          paircoeffs();
        } else skip_lines(ntypes);
      } else if (strcmp(keyword,"PairIJ Coeffs") == 0) {
        if (force->pair == NULL)
          error->all(FLERR,"Must define pair_style before PairIJ Coeffs");
        if (firstpass) {
          if (me == 0 && !style_match(style,force->pair_style))
            error->warning(FLERR,"Pair style in data file differs "
                           "from currently defined pair style");
          pairIJcoeffs();
        } else skip_lines(ntypes*(ntypes+1)/2);
      } else if (strcmp(keyword,"Bond Coeffs") == 0) {
        if (atom->avec->bonds_allow == 0)
          error->all(FLERR,"Invalid data file section: Bond Coeffs");
        if (force->bond == NULL)
          error->all(FLERR,"Must define bond_style before Bond Coeffs");
        if (firstpass) {
          if (me == 0 && !style_match(style,force->bond_style))
            error->warning(FLERR,"Bond style in data file differs "
                           "from currently defined bond style");
          bondcoeffs();
        } else skip_lines(nbondtypes);
      } else if (strcmp(keyword,"Angle Coeffs") == 0) {
        if (atom->avec->angles_allow == 0)
          error->all(FLERR,"Invalid data file section: Angle Coeffs");
        if (force->angle == NULL)
          error->all(FLERR,"Must define angle_style before Angle Coeffs");
        if (firstpass) {
          if (me == 0 && !style_match(style,force->angle_style))
            error->warning(FLERR,"Angle style in data file differs "
                           "from currently defined angle style");
          anglecoeffs(0);
        } else skip_lines(nangletypes);
      } else if (strcmp(keyword,"Dihedral Coeffs") == 0) {
        if (atom->avec->dihedrals_allow == 0)
          error->all(FLERR,"Invalid data file section: Dihedral Coeffs");
        if (force->dihedral == NULL)
          error->all(FLERR,"Must define dihedral_style before Dihedral Coeffs");
        if (firstpass) {
          if (me == 0 && !style_match(style,force->dihedral_style))
            error->warning(FLERR,"Dihedral style in data file differs "
                           "from currently defined dihedral style");
          dihedralcoeffs(0);
        } else skip_lines(ndihedraltypes);
      } else if (strcmp(keyword,"Improper Coeffs") == 0) {
        if (atom->avec->impropers_allow == 0)
          error->all(FLERR,"Invalid data file section: Improper Coeffs");
        if (force->improper == NULL)
          error->all(FLERR,"Must define improper_style before Improper Coeffs");
        if (firstpass) {
          if (me == 0 && !style_match(style,force->improper_style))
            error->warning(FLERR,"Improper style in data file differs "
                           "from currently defined improper style");
          impropercoeffs(0);
        } else skip_lines(nimpropertypes);

      } else if (strcmp(keyword,"BondBond Coeffs") == 0) {
        if (atom->avec->angles_allow == 0)
          error->all(FLERR,"Invalid data file section: BondBond Coeffs");
        if (force->angle == NULL)
          error->all(FLERR,"Must define angle_style before BondBond Coeffs");
        if (firstpass) anglecoeffs(1);
        else skip_lines(nangletypes);
      } else if (strcmp(keyword,"BondAngle Coeffs") == 0) {
        if (atom->avec->angles_allow == 0)
          error->all(FLERR,"Invalid data file section: BondAngle Coeffs");
        if (force->angle == NULL)
          error->all(FLERR,"Must define angle_style before BondAngle Coeffs");
        if (firstpass) anglecoeffs(2);
        else skip_lines(nangletypes);

      } else if (strcmp(keyword,"MiddleBondTorsion Coeffs") == 0) {
        if (atom->avec->dihedrals_allow == 0)
          error->all(FLERR,
                     "Invalid data file section: MiddleBondTorsion Coeffs");
        if (force->dihedral == NULL)
          error->all(FLERR,
                     "Must define dihedral_style before "
                     "MiddleBondTorsion Coeffs");
        if (firstpass) dihedralcoeffs(1);
        else skip_lines(ndihedraltypes);
      } else if (strcmp(keyword,"EndBondTorsion Coeffs") == 0) {
        if (atom->avec->dihedrals_allow == 0)
          error->all(FLERR,"Invalid data file section: EndBondTorsion Coeffs");
        if (force->dihedral == NULL)
          error->all(FLERR,
                     "Must define dihedral_style before EndBondTorsion Coeffs");
        if (firstpass) dihedralcoeffs(2);
        else skip_lines(ndihedraltypes);
      } else if (strcmp(keyword,"AngleTorsion Coeffs") == 0) {
        if (atom->avec->dihedrals_allow == 0)
          error->all(FLERR,"Invalid data file section: AngleTorsion Coeffs");
        if (force->dihedral == NULL)
          error->all(FLERR,
                     "Must define dihedral_style before AngleTorsion Coeffs");
        if (firstpass) dihedralcoeffs(3);
        else skip_lines(ndihedraltypes);
      } else if (strcmp(keyword,"AngleAngleTorsion Coeffs") == 0) {
        if (atom->avec->dihedrals_allow == 0)
          error->all(FLERR,
                     "Invalid data file section: AngleAngleTorsion Coeffs");
        if (force->dihedral == NULL)
          error->all(FLERR,
                     "Must define dihedral_style before "
                     "AngleAngleTorsion Coeffs");
        if (firstpass) dihedralcoeffs(4);
        else skip_lines(ndihedraltypes);
      } else if (strcmp(keyword,"BondBond13 Coeffs") == 0) {
        if (atom->avec->dihedrals_allow == 0)
          error->all(FLERR,"Invalid data file section: BondBond13 Coeffs");
        if (force->dihedral == NULL)
          error->all(FLERR,
                     "Must define dihedral_style before BondBond13 Coeffs");
        if (firstpass) dihedralcoeffs(5);
        else skip_lines(ndihedraltypes);

      } else if (strcmp(keyword,"AngleAngle Coeffs") == 0) {
        if (atom->avec->impropers_allow == 0)
          error->all(FLERR,"Invalid data file section: AngleAngle Coeffs");
        if (force->improper == NULL)
          error->all(FLERR,
                     "Must define improper_style before AngleAngle Coeffs");
        if (firstpass) impropercoeffs(1);
        else skip_lines(nimpropertypes);

      } else {
        char str[128];
        snprintf(str,128,"Unknown identifier in data file: %s",keyword);
        error->all(FLERR,str);
      }

      parse_keyword(0);
    }

    // error if natoms > 0 yet no atoms were read

    if (natoms > 0 && atomflag == 0)
      error->all(FLERR,"No atoms in data file");

    // close file

    if (me == 0) {
      if (compressed) pclose(fp);
      else fclose(fp);
      fp = NULL;
    }

    // done if this was 2nd pass

    if (!firstpass) break;

    // at end of 1st pass, error check for required sections
    // customize for new sections

    if ((nbonds && !bondflag) || (nangles && !angleflag) ||
        (ndihedrals && !dihedralflag) || (nimpropers && !improperflag))
      error->one(FLERR,"Needed molecular topology not in data file");

    if ((nellipsoids && !ellipsoidflag) || (nlines && !lineflag) ||
        (ntris && !triflag) || (nbodies && !bodyflag))
      error->one(FLERR,"Needed bonus data not in data file");

    // break out of loop if no molecular topology in file
    // else make 2nd pass

    if (!topoflag) break;
    firstpass = 0;

    // reallocate bond,angle,diehdral,improper arrays via grow()
    // will use new bond,angle,dihedral,improper per-atom values from 1st pass
    // will also observe extra settings even if bond/etc topology not in file
    // leaves other atom arrays unchanged, since already nmax in length

    if (addflag == NONE) atom->deallocate_topology();
    atom->avec->grow(atom->nmax);
  }

  // init per-atom fix/compute/variable values for created atoms

  atom->data_fix_compute_variable(nlocal_previous,atom->nlocal);

  // assign atoms added by this data file to specified group

  if (groupbit) {
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    for (int i = nlocal_previous; i < nlocal; i++)
      mask[i] |= groupbit;
  }

  // create special bond lists for molecular systems

  if (atom->molecular == 1) {
    Special special(lmp);
    special.build();
  }

  // for atom style template systems, count total bonds,angles,etc

  if (atom->molecular == 2) {
    Molecule **onemols = atom->avec->onemols;
    int *molindex = atom->molindex;
    int *molatom = atom->molatom;
    int nlocal = atom->nlocal;

    int imol,iatom;
    bigint nbonds,nangles,ndihedrals,nimpropers;
    nbonds = nangles = ndihedrals = nimpropers = 0;

    for (int i = 0; i < nlocal; i++) {
      imol = molindex[i];
      iatom = molatom[i];
      nbonds += onemols[imol]->num_bond[iatom];
      nangles += onemols[imol]->num_angle[iatom];
      ndihedrals += onemols[imol]->num_dihedral[iatom];
      nimpropers += onemols[imol]->num_improper[iatom];
    }

    MPI_Allreduce(&nbonds,&atom->nbonds,1,MPI_LMP_BIGINT,MPI_SUM,world);
    MPI_Allreduce(&nangles,&atom->nangles,1,MPI_LMP_BIGINT,MPI_SUM,world);
    MPI_Allreduce(&ndihedrals,&atom->ndihedrals,1,MPI_LMP_BIGINT,MPI_SUM,world);
    MPI_Allreduce(&nimpropers,&atom->nimpropers,1,MPI_LMP_BIGINT,MPI_SUM,world);

    if (!force->newton_bond) {
      atom->nbonds /= 2;
      atom->nangles /= 3;
      atom->ndihedrals /= 4;
      atom->nimpropers /= 4;
    }

    if (me == 0) {
      if (atom->nbonds) {
        if (screen)
          fprintf(screen,"  " BIGINT_FORMAT " template bonds\n",atom->nbonds);
        if (logfile)
          fprintf(logfile,"  " BIGINT_FORMAT " template bonds\n",atom->nbonds);
      }
      if (atom->nangles) {
        if (screen)
          fprintf(screen,"  " BIGINT_FORMAT " template angles\n",
                  atom->nangles);
        if (logfile)
          fprintf(logfile,"  " BIGINT_FORMAT " template angles\n",
                  atom->nangles);
      }
      if (atom->ndihedrals) {
        if (screen)
          fprintf(screen,"  " BIGINT_FORMAT " template dihedrals\n",
                  atom->nbonds);
        if (logfile)
          fprintf(logfile,"  " BIGINT_FORMAT " template bonds\n",
                  atom->ndihedrals);
      }
      if (atom->nimpropers) {
        if (screen)
          fprintf(screen,"  " BIGINT_FORMAT " template impropers\n",
                  atom->nimpropers);
        if (logfile)
          fprintf(logfile,"  " BIGINT_FORMAT " template impropers\n",
                  atom->nimpropers);
      }
    }
  }

  // for atom style template systems
  // insure nbondtypes,etc are still consistent with template molecules,
  //   in case data file re-defined them

  if (atom->molecular == 2) atom->avec->onemols[0]->check_attributes(1);

  // if adding atoms, migrate atoms to new processors
  // use irregular() b/c box size could have changed dramaticaly
  // resulting in procs now owning very different subboxes
  // with their previously owned atoms now far outside the subbox

  if (addflag != NONE) {
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    Irregular *irregular = new Irregular(lmp);
    irregular->migrate_atoms(1);
    delete irregular;
    if (domain->triclinic) domain->lamda2x(atom->nlocal);
  }

  // shrink-wrap the box if necessary and move atoms to new procs
  // if atoms are lost is b/c data file box was far from shrink-wrapped
  // do not use irregular() comm, which would not lose atoms,
  //   b/c then user could specify data file box as far too big and empty
  // do comm->init() but not comm->setup() b/c pair/neigh cutoffs not yet set
  // need call to map_set() b/c comm->exchange clears atom map

  if (domain->nonperiodic == 2) {
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    domain->reset_box();
    Irregular *irregular = new Irregular(lmp);
    irregular->migrate_atoms(1);
    delete irregular;
    if (domain->triclinic) domain->lamda2x(atom->nlocal);

    bigint natoms;
    bigint nblocal = atom->nlocal;
    MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
    if (natoms != atom->natoms)
      error->all(FLERR,
                 "Read_data shrink wrap did not assign all atoms correctly");
  }

  // restore old styles, when reading with nocoeff flag given

  if (coeffflag == 0) {
    if (force->pair) delete force->pair;
    force->pair = saved_pair;

    if (force->bond) delete force->bond;
    force->bond = saved_bond;

    if (force->angle) delete force->angle;
    force->angle = saved_angle;

    if (force->dihedral) delete force->dihedral;
    force->dihedral = saved_dihedral;

    if (force->improper) delete force->improper;
    force->improper = saved_improper;

    force->kspace = saved_kspace;
  }

  // total time

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"  read_data CPU = %g secs\n",time2-time1);
    if (logfile)
      fprintf(logfile,"  read_data CPU = %g secs\n",time2-time1);
  }
}

/* ----------------------------------------------------------------------
   read free-format header of data file
   1st line and blank lines are skipped
   non-blank lines are checked for header keywords and leading value is read
   header ends with EOF or non-blank line containing no header keyword
     if EOF, line is set to blank line
     else line has first keyword line for rest of file
   some logic differs if adding atoms
------------------------------------------------------------------------- */

void ReadData::header(int firstpass)
{
  int n;
  char *ptr;

  // customize for new sections

  const char *section_keywords[NSECTIONS] =
    {"Atoms","Velocities","Ellipsoids","Lines","Triangles","Bodies",
     "Bonds","Angles","Dihedrals","Impropers",
     "Masses","Pair Coeffs","PairIJ Coeffs","Bond Coeffs","Angle Coeffs",
     "Dihedral Coeffs","Improper Coeffs",
     "BondBond Coeffs","BondAngle Coeffs","MiddleBondTorsion Coeffs",
     "EndBondTorsion Coeffs","AngleTorsion Coeffs",
     "AngleAngleTorsion Coeffs","BondBond13 Coeffs","AngleAngle Coeffs"};

  // skip 1st line of file

  if (me == 0) {
    char *eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
  }

  while (1) {

    // read a line and bcast length

    if (me == 0) {
      if (fgets(line,MAXLINE,fp) == NULL) n = 0;
      else n = strlen(line) + 1;
    }
    MPI_Bcast(&n,1,MPI_INT,0,world);

    // if n = 0 then end-of-file so return with blank line

    if (n == 0) {
      line[0] = '\0';
      return;
    }

    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // trim anything from '#' onward
    // if line is blank, continue

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    if (strspn(line," \t\n\r") == strlen(line)) continue;

    // allow special fixes first chance to match and process the line
    // if fix matches, continue to next header line

    if (nfix) {
      for (n = 0; n < nfix; n++) {
        if (!fix_header[n]) continue;
        if (strstr(line,fix_header[n])) {
          modify->fix[fix_index[n]]->read_data_header(line);
          break;
        }
      }
      if (n < nfix) continue;
    }

    // search line for header keyword and set corresponding variable
    // customize for new header lines
    // check for triangles before angles so "triangles" not matched as "angles"
    int extra_flag_value = 0;

    if (strstr(line,"atoms")) {
      sscanf(line,BIGINT_FORMAT,&natoms);
      if (addflag == NONE) atom->natoms = natoms;
      else if (firstpass) atom->natoms += natoms;

    } else if (strstr(line,"ellipsoids")) {
      if (!avec_ellipsoid)
        error->all(FLERR,"No ellipsoids allowed with this atom style");
      sscanf(line,BIGINT_FORMAT,&nellipsoids);
      if (addflag == NONE) atom->nellipsoids = nellipsoids;
      else if (firstpass) atom->nellipsoids += nellipsoids;

    } else if (strstr(line,"lines")) {
      if (!avec_line)
        error->all(FLERR,"No lines allowed with this atom style");
      sscanf(line,BIGINT_FORMAT,&nlines);
      if (addflag == NONE) atom->nlines = nlines;
      else if (firstpass) atom->nlines += nlines;

    } else if (strstr(line,"triangles")) {
      if (!avec_tri)
        error->all(FLERR,"No triangles allowed with this atom style");
      sscanf(line,BIGINT_FORMAT,&ntris);
      if (addflag == NONE) atom->ntris = ntris;
      else if (firstpass) atom->ntris += ntris;

    } else if (strstr(line,"bodies")) {
      if (!avec_body)
        error->all(FLERR,"No bodies allowed with this atom style");
      sscanf(line,BIGINT_FORMAT,&nbodies);
      if (addflag == NONE) atom->nbodies = nbodies;
      else if (firstpass) atom->nbodies += nbodies;

    } else if (strstr(line,"bonds")) {
      sscanf(line,BIGINT_FORMAT,&nbonds);
      if (addflag == NONE) atom->nbonds = nbonds;
      else if (firstpass) atom->nbonds += nbonds;
    } else if (strstr(line,"angles")) {
      sscanf(line,BIGINT_FORMAT,&nangles);
      if (addflag == NONE) atom->nangles = nangles;
      else if (firstpass) atom->nangles += nangles;
    } else if (strstr(line,"dihedrals")) {
      sscanf(line,BIGINT_FORMAT,&ndihedrals);
      if (addflag == NONE) atom->ndihedrals = ndihedrals;
      else if (firstpass) atom->ndihedrals += ndihedrals;
    } else if (strstr(line,"impropers")) {
      sscanf(line,BIGINT_FORMAT,&nimpropers);
      if (addflag == NONE) atom->nimpropers = nimpropers;
      else if (firstpass) atom->nimpropers += nimpropers;

    // Atom class type settings are only set by first data file

    } else if (strstr(line,"atom types")) {
      sscanf(line,"%d",&ntypes);
      if (addflag == NONE) atom->ntypes = ntypes + extra_atom_types;
    } else if (strstr(line,"bond types")) {
      sscanf(line,"%d",&nbondtypes);
      if (addflag == NONE) atom->nbondtypes = nbondtypes + extra_bond_types;
    } else if (strstr(line,"angle types")) {
      sscanf(line,"%d",&nangletypes);
      if (addflag == NONE) atom->nangletypes = nangletypes + extra_angle_types;
    } else if (strstr(line,"dihedral types")) {
      sscanf(line,"%d",&ndihedraltypes);
      if (addflag == NONE)
        atom->ndihedraltypes = ndihedraltypes + extra_dihedral_types;
    } else if (strstr(line,"improper types")) {
      sscanf(line,"%d",&nimpropertypes);
      if (addflag == NONE)
        atom->nimpropertypes = nimpropertypes + extra_improper_types;

    // these settings only used by first data file
    // also, these are obsolescent. we parse them to maintain backward
    // compatibility, but the recommended way is to set them via keywords
    // in the LAMMPS input file. In case these flags are set in both,
    // the input and the data file, we use the larger of the two.

    } else if (strstr(line,"extra bond per atom")) {
      if (addflag == NONE) sscanf(line,"%d",&extra_flag_value);
      atom->extra_bond_per_atom = MAX(atom->extra_bond_per_atom,extra_flag_value);
    } else if (strstr(line,"extra angle per atom")) {
      if (addflag == NONE) sscanf(line,"%d",&extra_flag_value);
      atom->extra_angle_per_atom = MAX(atom->extra_angle_per_atom,extra_flag_value);
    } else if (strstr(line,"extra dihedral per atom")) {
      if (addflag == NONE) sscanf(line,"%d",&extra_flag_value);
      atom->extra_dihedral_per_atom = MAX(atom->extra_dihedral_per_atom,extra_flag_value);
    } else if (strstr(line,"extra improper per atom")) {
      if (addflag == NONE) sscanf(line,"%d",&extra_flag_value);
      atom->extra_improper_per_atom = MAX(atom->extra_improper_per_atom,extra_flag_value);
    } else if (strstr(line,"extra special per atom")) {
      if (addflag == NONE) sscanf(line,"%d",&extra_flag_value);
      force->special_extra = MAX(force->special_extra,extra_flag_value);

    // local copy of box info
    // so can treat differently for first vs subsequent data files

    } else if (strstr(line,"xlo xhi")) {
      sscanf(line,"%lg %lg",&boxlo[0],&boxhi[0]);
    } else if (strstr(line,"ylo yhi")) {
      sscanf(line,"%lg %lg",&boxlo[1],&boxhi[1]);
    } else if (strstr(line,"zlo zhi")) {
      sscanf(line,"%lg %lg",&boxlo[2],&boxhi[2]);
    } else if (strstr(line,"xy xz yz")) {
      triclinic = 1;
      sscanf(line,"%lg %lg %lg",&xy,&xz,&yz);

    } else break;
  }

  // error check on total system size

  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT ||
      atom->nellipsoids < 0 || atom->nellipsoids >= MAXBIGINT ||
      atom->nlines < 0 || atom->nlines >= MAXBIGINT ||
      atom->ntris < 0 || atom->ntris >= MAXBIGINT ||
      atom->nbodies < 0 || atom->nbodies >= MAXBIGINT ||
      atom->nbonds < 0 || atom->nbonds >= MAXBIGINT ||
      atom->nangles < 0 || atom->nangles >= MAXBIGINT ||
      atom->ndihedrals < 0 || atom->ndihedrals >= MAXBIGINT ||
      atom->nimpropers < 0 || atom->nimpropers >= MAXBIGINT)
    error->all(FLERR,"System in data file is too big");

  // check that exiting string is a valid section keyword

  parse_keyword(1);
  for (n = 0; n < NSECTIONS; n++)
    if (strcmp(keyword,section_keywords[n]) == 0) break;
  if (n == NSECTIONS) {
    char str[128];
    sprintf(str,"Unknown identifier in data file: %s",keyword);
    error->all(FLERR,str);
  }

  // error checks on header values
  // must be consistent with atom style and other header values

  if ((atom->nbonds || atom->nbondtypes) &&
      atom->avec->bonds_allow == 0)
    error->all(FLERR,"No bonds allowed with this atom style");
  if ((atom->nangles || atom->nangletypes) &&
      atom->avec->angles_allow == 0)
    error->all(FLERR,"No angles allowed with this atom style");
  if ((atom->ndihedrals || atom->ndihedraltypes) &&
      atom->avec->dihedrals_allow == 0)
    error->all(FLERR,"No dihedrals allowed with this atom style");
  if ((atom->nimpropers || atom->nimpropertypes) &&
      atom->avec->impropers_allow == 0)
    error->all(FLERR,"No impropers allowed with this atom style");

  if (atom->nbonds > 0 && atom->nbondtypes <= 0)
    error->all(FLERR,"Bonds defined but no bond types");
  if (atom->nangles > 0 && atom->nangletypes <= 0)
    error->all(FLERR,"Angles defined but no angle types");
  if (atom->ndihedrals > 0 && atom->ndihedraltypes <= 0)
    error->all(FLERR,"Dihedrals defined but no dihedral types");
  if (atom->nimpropers > 0 && atom->nimpropertypes <= 0)
    error->all(FLERR,"Impropers defined but no improper types");

  if (atom->molecular == 2) {
    if (atom->nbonds || atom->nangles || atom->ndihedrals || atom->nimpropers)
      error->all(FLERR,"No molecule topology allowed with atom style template");
  }
}

/* ----------------------------------------------------------------------
   read all atoms
------------------------------------------------------------------------- */

void ReadData::atoms()
{
  int nchunk,eof;

  if (me == 0) {
    if (screen) fprintf(screen,"  reading atoms ...\n");
    if (logfile) fprintf(logfile,"  reading atoms ...\n");
  }

  bigint nread = 0;

  while (nread < natoms) {
    nchunk = MIN(natoms-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_atoms(nchunk,buffer,id_offset,mol_offset,toffset,shiftflag,shift);
    nread += nchunk;
  }

  // check that all atoms were assigned correctly

  bigint n = atom->nlocal;
  bigint sum;
  MPI_Allreduce(&n,&sum,1,MPI_LMP_BIGINT,MPI_SUM,world);
  bigint nassign = sum - (atom->natoms - natoms);

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " atoms\n",nassign);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " atoms\n",nassign);
  }

  if (sum != atom->natoms)
    error->all(FLERR,"Did not assign all atoms correctly");

  // check that atom IDs are valid

  atom->tag_check();

  // check that bonus data has been reserved as needed

  atom->bonus_check();

  // create global mapping of atoms

  if (atom->map_style) {
    atom->map_init();
    atom->map_set();
  }
}

/* ----------------------------------------------------------------------
   read all velocities
   to find atoms, must build atom map if not a molecular system
------------------------------------------------------------------------- */

void ReadData::velocities()
{
  int nchunk,eof;

  if (me == 0) {
    if (screen) fprintf(screen,"  reading velocities ...\n");
    if (logfile) fprintf(logfile,"  reading velocities ...\n");
  }

  int mapflag = 0;
  if (atom->map_style == 0) {
    mapflag = 1;
    atom->map_init();
    atom->map_set();
  }

  bigint nread = 0;

  while (nread < natoms) {
    nchunk = MIN(natoms-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_vels(nchunk,buffer,id_offset);
    nread += nchunk;
  }

  if (mapflag) {
    atom->map_delete();
    atom->map_style = 0;
  }

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " velocities\n",natoms);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " velocities\n",natoms);
  }
}

/* ----------------------------------------------------------------------
   scan or read all bonds
------------------------------------------------------------------------- */

void ReadData::bonds(int firstpass)
{
  int nchunk,eof;

  if (me == 0) {
    if (firstpass) {
      if (screen) fprintf(screen,"  scanning bonds ...\n");
      if (logfile) fprintf(logfile,"  scanning bonds ...\n");
    } else {
      if (screen) fprintf(screen,"  reading bonds ...\n");
      if (logfile) fprintf(logfile,"  reading bonds ...\n");
    }
  }

  // allocate count if firstpass

  int nlocal = atom->nlocal;
  int *count = NULL;
  if (firstpass) {
    memory->create(count,nlocal,"read_data:count");
    memset(count,0,nlocal*sizeof(int));
  }

  // read and process bonds

  bigint nread = 0;

  while (nread < nbonds) {
    nchunk = MIN(nbonds-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_bonds(nchunk,buffer,count,id_offset,boffset);
    nread += nchunk;
  }

  // if firstpass: tally max bond/atom and return
  // if addflag = NONE, store max bond/atom with extra
  // else just check actual max does not exceed existing max

  if (firstpass) {
    int max = 0;
    for (int i = nlocal_previous; i < nlocal; i++) max = MAX(max,count[i]);
    int maxall;
    MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);
    if (addflag == NONE) maxall += atom->extra_bond_per_atom;

    if (me == 0) {
      if (screen) fprintf(screen,"  %d = max bonds/atom\n",maxall);
      if (logfile) fprintf(logfile,"  %d = max bonds/atom\n",maxall);
    }

    if (addflag != NONE) {
      if (maxall > atom->bond_per_atom)
        error->all(FLERR,"Subsequent read data induced "
                   "too many bonds per atom");
    } else atom->bond_per_atom = maxall;

    memory->destroy(count);
    return;
  }

  // if 2nd pass: check that bonds were assigned correctly

  bigint n = 0;
  for (int i = nlocal_previous; i < nlocal; i++) n += atom->num_bond[i];
  bigint sum;
  MPI_Allreduce(&n,&sum,1,MPI_LMP_BIGINT,MPI_SUM,world);
  int factor = 1;
  if (!force->newton_bond) factor = 2;

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " bonds\n",sum/factor);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " bonds\n",sum/factor);
  }

  if (sum != factor*nbonds)
    error->all(FLERR,"Bonds assigned incorrectly");
}

/* ----------------------------------------------------------------------
   scan or read all angles
------------------------------------------------------------------------- */

void ReadData::angles(int firstpass)
{
  int nchunk,eof;

  if (me == 0) {
    if (firstpass) {
      if (screen) fprintf(screen,"  scanning angles ...\n");
      if (logfile) fprintf(logfile,"  scanning angles ...\n");
    } else {
      if (screen) fprintf(screen,"  reading angles ...\n");
      if (logfile) fprintf(logfile,"  reading angles ...\n");
    }
  }

  // allocate count if firstpass

  int nlocal = atom->nlocal;
  int *count = NULL;
  if (firstpass) {
    memory->create(count,nlocal,"read_data:count");
    memset(count,0,nlocal*sizeof(int));
  }

  // read and process angles

  bigint nread = 0;

  while (nread < nangles) {
    nchunk = MIN(nangles-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_angles(nchunk,buffer,count,id_offset,aoffset);
    nread += nchunk;
  }

  // if firstpass: tally max angle/atom and return
  // if addflag = NONE, store max angle/atom with extra
  // else just check actual max does not exceed existing max

  if (firstpass) {
    int max = 0;
    for (int i = nlocal_previous; i < nlocal; i++) max = MAX(max,count[i]);
    int maxall;
    MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);
    if (addflag == NONE) maxall += atom->extra_angle_per_atom;

    if (me == 0) {
      if (screen) fprintf(screen,"  %d = max angles/atom\n",maxall);
      if (logfile) fprintf(logfile,"  %d = max angles/atom\n",maxall);
    }

    if (addflag != NONE) {
      if (maxall > atom->angle_per_atom)
        error->all(FLERR,"Subsequent read data induced "
                   "too many angles per atom");
    } else atom->angle_per_atom = maxall;

    memory->destroy(count);
    return;
  }

  // if 2nd pass: check that angles were assigned correctly

  bigint n = 0;
  for (int i = nlocal_previous; i < nlocal; i++) n += atom->num_angle[i];
  bigint sum;
  MPI_Allreduce(&n,&sum,1,MPI_LMP_BIGINT,MPI_SUM,world);
  int factor = 1;
  if (!force->newton_bond) factor = 3;

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " angles\n",sum/factor);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " angles\n",sum/factor);
  }

  if (sum != factor*nangles)
    error->all(FLERR,"Angles assigned incorrectly");
}

/* ----------------------------------------------------------------------
   scan or read all dihedrals
------------------------------------------------------------------------- */

void ReadData::dihedrals(int firstpass)
{
  int nchunk,eof;

  if (me == 0) {
    if (firstpass) {
      if (screen) fprintf(screen,"  scanning dihedrals ...\n");
      if (logfile) fprintf(logfile,"  scanning dihedrals ...\n");
    } else {
      if (screen) fprintf(screen,"  reading dihedrals ...\n");
      if (logfile) fprintf(logfile,"  reading dihedrals ...\n");
    }
  }

  // allocate count if firstpass

  int nlocal = atom->nlocal;
  int *count = NULL;
  if (firstpass) {
    memory->create(count,nlocal,"read_data:count");
    memset(count,0,nlocal*sizeof(int));
  }

  // read and process dihedrals

  bigint nread = 0;

  while (nread < ndihedrals) {
    nchunk = MIN(ndihedrals-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_dihedrals(nchunk,buffer,count,id_offset,doffset);
    nread += nchunk;
  }

  // if firstpass: tally max dihedral/atom and return
  // if addflag = NONE, store max dihedral/atom with extra
  // else just check actual max does not exceed existing max

  if (firstpass) {
    int max = 0;
    for (int i = nlocal_previous; i < nlocal; i++) max = MAX(max,count[i]);
    int maxall;
    MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);
    if (addflag == NONE) maxall += atom->extra_dihedral_per_atom;

    if (me == 0) {
      if (screen) fprintf(screen,"  %d = max dihedrals/atom\n",maxall);
      if (logfile) fprintf(logfile,"  %d = max dihedrals/atom\n",maxall);
    }

    if (addflag != NONE) {
      if (maxall > atom->dihedral_per_atom)
        error->all(FLERR,"Subsequent read data induced "
                   "too many dihedrals per atom");
    } else atom->dihedral_per_atom = maxall;

    memory->destroy(count);
    return;
  }

  // if 2nd pass: check that dihedrals were assigned correctly

  bigint n = 0;
  for (int i = nlocal_previous; i < nlocal; i++) n += atom->num_dihedral[i];
  bigint sum;
  MPI_Allreduce(&n,&sum,1,MPI_LMP_BIGINT,MPI_SUM,world);
  int factor = 1;
  if (!force->newton_bond) factor = 4;

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " dihedrals\n",sum/factor);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " dihedrals\n",sum/factor);
  }

  if (sum != factor*ndihedrals)
    error->all(FLERR,"Dihedrals assigned incorrectly");
}

/* ----------------------------------------------------------------------
   scan or read all impropers
------------------------------------------------------------------------- */

void ReadData::impropers(int firstpass)
{
  int nchunk,eof;

  if (me == 0) {
    if (firstpass) {
      if (screen) fprintf(screen,"  scanning impropers ...\n");
      if (logfile) fprintf(logfile,"  scanning impropers ...\n");
    } else {
      if (screen) fprintf(screen,"  reading impropers ...\n");
      if (logfile) fprintf(logfile,"  reading impropers ...\n");
    }
  }

  // allocate count if firstpass

  int nlocal = atom->nlocal;
  int *count = NULL;
  if (firstpass) {
    memory->create(count,nlocal,"read_data:count");
    memset(count,0,nlocal*sizeof(int));
  }

  // read and process impropers

  bigint nread = 0;

  while (nread < nimpropers) {
    nchunk = MIN(nimpropers-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_impropers(nchunk,buffer,count,id_offset,ioffset);
    nread += nchunk;
  }

  // if firstpass: tally max improper/atom and return
  // if addflag = NONE, store max improper/atom
  // else just check it does not exceed existing max

  if (firstpass) {
    int max = 0;
    for (int i = nlocal_previous; i < nlocal; i++) max = MAX(max,count[i]);
    int maxall;
    MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);
    if (addflag == NONE) maxall += atom->extra_improper_per_atom;

    if (me == 0) {
      if (screen) fprintf(screen,"  %d = max impropers/atom\n",maxall);
      if (logfile) fprintf(logfile,"  %d = max impropers/atom\n",maxall);
    }

    if (addflag != NONE) {
      if (maxall > atom->improper_per_atom)
        error->all(FLERR,"Subsequent read data induced "
                   "too many impropers per atom");
    } else atom->improper_per_atom = maxall;

    memory->destroy(count);
    return;
  }

  // if 2nd pass: check that impropers were assigned correctly

  bigint n = 0;
  for (int i = nlocal_previous; i < nlocal; i++) n += atom->num_improper[i];
  bigint sum;
  MPI_Allreduce(&n,&sum,1,MPI_LMP_BIGINT,MPI_SUM,world);
  int factor = 1;
  if (!force->newton_bond) factor = 4;

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " impropers\n",sum/factor);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " impropers\n",sum/factor);
  }

  if (sum != factor*nimpropers)
    error->all(FLERR,"Impropers assigned incorrectly");
}

/* ----------------------------------------------------------------------
   read all bonus data
   to find atoms, must build atom map if not a molecular system
------------------------------------------------------------------------- */

void ReadData::bonus(bigint nbonus, AtomVec *ptr, const char *type)
{
  int nchunk,eof;

  int mapflag = 0;
  if (atom->map_style == 0) {
    mapflag = 1;
    atom->map_init();
    atom->map_set();
  }

  bigint nread = 0;
  bigint natoms = nbonus;

  while (nread < natoms) {
    nchunk = MIN(natoms-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_bonus(nchunk,buffer,ptr,id_offset);
    nread += nchunk;
  }

  if (mapflag) {
    atom->map_delete();
    atom->map_style = 0;
  }

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " %s\n",natoms,type);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " %s\n",natoms,type);
  }
}

/* ----------------------------------------------------------------------
   read all body data
   variable amount of info per body, described by ninteger and ndouble
   to find atoms, must build atom map if not a molecular system
   if not firstpass, just read past data, but no processing of data
------------------------------------------------------------------------- */

void ReadData::bodies(int firstpass)
{
  int m,nchunk,nline,nmax,ninteger,ndouble,nword,ncount,onebody,tmp;
  char *eof;

  int mapflag = 0;
  if (atom->map_style == 0 && firstpass) {
    mapflag = 1;
    atom->map_init();
    atom->map_set();
  }

  // nmax = max # of bodies to read in this chunk
  // nchunk = actual # read

  bigint nread = 0;
  bigint natoms = nbodies;

  while (nread < natoms) {
    if (natoms-nread > CHUNK) nmax = CHUNK;
    else nmax = natoms-nread;

    if (me == 0) {
      nchunk = 0;
      nline = 0;
      m = 0;

      while (nchunk < nmax && nline <= CHUNK-MAXBODY) {
        eof = fgets(&buffer[m],MAXLINE,fp);
        if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
        sscanf(&buffer[m],"%d %d %d",&tmp,&ninteger,&ndouble);
        m += strlen(&buffer[m]);

        // read lines one at a time into buffer and count words
        // count to ninteger and ndouble until have enough lines

        onebody = 0;

        nword = 0;
        while (nword < ninteger) {
          eof = fgets(&buffer[m],MAXLINE,fp);
          if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
          ncount = atom->count_words(&buffer[m],copy);
          if (ncount == 0)
            error->one(FLERR,"Too few values in body lines in data file");
          nword += ncount;
          m += strlen(&buffer[m]);
          onebody++;
        }
        if (nword > ninteger)
          error->one(FLERR,"Too many values in body lines in data file");

        nword = 0;
        while (nword < ndouble) {
          eof = fgets(&buffer[m],MAXLINE,fp);
          if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
          ncount = atom->count_words(&buffer[m],copy);
          if (ncount == 0)
            error->one(FLERR,"Too few values in body lines in data file");
          nword += ncount;
          m += strlen(&buffer[m]);
          onebody++;
        }
        if (nword > ndouble)
          error->one(FLERR,"Too many values in body lines in data file");

        if (onebody+1 > MAXBODY)
          error->one(FLERR,
                     "Too many lines in one body in data file - boost MAXBODY");

        nchunk++;
        nline += onebody+1;
      }

      if (buffer[m-1] != '\n') strcpy(&buffer[m++],"\n");
      m++;
    }

    MPI_Bcast(&nchunk,1,MPI_INT,0,world);
    MPI_Bcast(&m,1,MPI_INT,0,world);
    MPI_Bcast(buffer,m,MPI_CHAR,0,world);

    if (firstpass) atom->data_bodies(nchunk,buffer,avec_body,id_offset);
    nread += nchunk;
  }

  if (mapflag && firstpass) {
    atom->map_delete();
    atom->map_style = 0;
  }

  if (me == 0 && firstpass) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " bodies\n",natoms);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " bodies\n",natoms);
  }
}

/* ---------------------------------------------------------------------- */

void ReadData::mass()
{
  char *next;
  char *buf = new char[ntypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,ntypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < ntypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    atom->set_mass(FLERR,buf,toffset);
    buf = next + 1;
  }
  delete [] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::paircoeffs()
{
  char *next;
  char *buf = new char[ntypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,ntypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < ntypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    parse_coeffs(buf,NULL,1,2,toffset);
    if (narg == 0)
      error->all(FLERR,"Unexpected empty line in PairCoeffs section");
    force->pair->coeff(narg,arg);
    buf = next + 1;
  }
  delete [] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::pairIJcoeffs()
{
  int i,j;
  char *next;

  int nsq = ntypes * (ntypes+1) / 2;
  char *buf = new char[nsq * MAXLINE];

  int eof = comm->read_lines_from_file(fp,nsq,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (i = 0; i < ntypes; i++)
    for (j = i; j < ntypes; j++) {
      next = strchr(buf,'\n');
      *next = '\0';
      parse_coeffs(buf,NULL,0,2,toffset);
      if (narg == 0)
        error->all(FLERR,"Unexpected empty line in PairCoeffs section");
      force->pair->coeff(narg,arg);
      buf = next + 1;
    }
  delete [] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::bondcoeffs()
{
  if (!nbondtypes) return;

  char *next;
  char *buf = new char[nbondtypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,nbondtypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < nbondtypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    parse_coeffs(buf,NULL,0,1,boffset);
    if (narg == 0)
      error->all(FLERR,"Unexpected empty line in BondCoeffs section");
    force->bond->coeff(narg,arg);
    buf = next + 1;
  }
  delete [] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::anglecoeffs(int which)
{
  if (!nangletypes) return;

  char *next;
  char *buf = new char[nangletypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,nangletypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < nangletypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    if (which == 0) parse_coeffs(buf,NULL,0,1,aoffset);
    else if (which == 1) parse_coeffs(buf,"bb",0,1,aoffset);
    else if (which == 2) parse_coeffs(buf,"ba",0,1,aoffset);
    if (narg == 0) error->all(FLERR,"Unexpected empty line in AngleCoeffs section");
    force->angle->coeff(narg,arg);
    buf = next + 1;
  }
  delete [] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::dihedralcoeffs(int which)
{
  if (!ndihedraltypes) return;

  char *next;
  char *buf = new char[ndihedraltypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,ndihedraltypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < ndihedraltypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    if (which == 0) parse_coeffs(buf,NULL,0,1,doffset);
    else if (which == 1) parse_coeffs(buf,"mbt",0,1,doffset);
    else if (which == 2) parse_coeffs(buf,"ebt",0,1,doffset);
    else if (which == 3) parse_coeffs(buf,"at",0,1,doffset);
    else if (which == 4) parse_coeffs(buf,"aat",0,1,doffset);
    else if (which == 5) parse_coeffs(buf,"bb13",0,1,doffset);
    if (narg == 0)
      error->all(FLERR,"Unexpected empty line in DihedralCoeffs section");
    force->dihedral->coeff(narg,arg);
    buf = next + 1;
  }
  delete [] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::impropercoeffs(int which)
{
  if (!nimpropertypes) return;

  char *next;
  char *buf = new char[nimpropertypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,nimpropertypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < nimpropertypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    if (which == 0) parse_coeffs(buf,NULL,0,1,ioffset);
    else if (which == 1) parse_coeffs(buf,"aa",0,1,ioffset);
    if (narg == 0) error->all(FLERR,"Unexpected empty line in ImproperCoeffs section");
    force->improper->coeff(narg,arg);
    buf = next + 1;
  }
  delete [] original;
}

/* ----------------------------------------------------------------------
   read fix section, pass lines to fix to process
   n = index of fix
------------------------------------------------------------------------- */

void ReadData::fix(int ifix, char *keyword)
{
  int nchunk,eof;

  bigint nline = modify->fix[ifix]->read_data_skip_lines(keyword);

  bigint nread = 0;
  while (nread < nline) {
    nchunk = MIN(nline-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    modify->fix[ifix]->read_data_section(keyword,nchunk,buffer,id_offset);
    nread += nchunk;
  }
}

/* ----------------------------------------------------------------------
   reallocate the count vector from cmax to amax+1 and return new length
   zero new locations
------------------------------------------------------------------------- */

int ReadData::reallocate(int **pcount, int cmax, int amax)
{
  int *count = *pcount;
  memory->grow(count,amax+1,"read_data:count");
  for (int i = cmax; i <= amax; i++) count[i] = 0;
  *pcount = count;
  return amax+1;
}

/* ----------------------------------------------------------------------
   proc 0 opens data file
   test if gzipped
------------------------------------------------------------------------- */

void ReadData::open(char *file)
{
  compressed = 0;
  char *suffix = file + strlen(file) - 3;
  if (suffix > file && strcmp(suffix,".gz") == 0) compressed = 1;
  if (!compressed) fp = fopen(file,"r");
  else {
#ifdef LAMMPS_GZIP
    char gunzip[128];
    snprintf(gunzip,128,"gzip -c -d %s",file);

#ifdef _WIN32
    fp = _popen(gunzip,"rb");
#else
    fp = popen(gunzip,"r");
#endif

#else
    error->one(FLERR,"Cannot open gzipped file");
#endif
  }

  if (fp == NULL) {
    char str[128];
    snprintf(str,128,"Cannot open file %s",file);
    error->one(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   grab next keyword
   read lines until one is non-blank
   keyword is all text on line w/out leading & trailing white space
   optional style can be appended after comment char '#'
   read one additional line (assumed blank)
   if any read hits EOF, set keyword to empty
   if first = 1, line variable holds non-blank line that ended header
------------------------------------------------------------------------- */

void ReadData::parse_keyword(int first)
{
  int eof = 0;
  int done = 0;

  // proc 0 reads upto non-blank line plus 1 following line
  // eof is set to 1 if any read hits end-of-file

  if (me == 0) {
    if (!first) {
      if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
    }
    while (eof == 0 && done == 0) {
      int blank = strspn(line," \t\n\r");
      if ((blank == (int)strlen(line)) || (line[blank] == '#')) {
        if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
      } else done = 1;
    }
    if (fgets(buffer,MAXLINE,fp) == NULL) {
      eof = 1;
      buffer[0] = '\0';
    }
  }

  // if eof, set keyword empty and return

  MPI_Bcast(&eof,1,MPI_INT,0,world);
  if (eof) {
    keyword[0] = '\0';
    return;
  }

  // bcast keyword line to all procs

  int n;
  if (me == 0) n = strlen(line) + 1;
  MPI_Bcast(&n,1,MPI_INT,0,world);
  MPI_Bcast(line,n,MPI_CHAR,0,world);

  // store optional "style" following comment char '#' after keyword

  char *ptr;
  if ((ptr = strchr(line,'#'))) {
    *ptr++ = '\0';
    while (*ptr == ' ' || *ptr == '\t') ptr++;
    int stop = strlen(ptr) - 1;
    while (ptr[stop] == ' ' || ptr[stop] == '\t'
           || ptr[stop] == '\n' || ptr[stop] == '\r') stop--;
    ptr[stop+1] = '\0';
    strcpy(style,ptr);
  } else style[0] = '\0';

  // copy non-whitespace portion of line into keyword

  int start = strspn(line," \t\n\r");
  int stop = strlen(line) - 1;
  while (line[stop] == ' ' || line[stop] == '\t'
         || line[stop] == '\n' || line[stop] == '\r') stop--;
  line[stop+1] = '\0';
  strcpy(keyword,&line[start]);
}

/* ----------------------------------------------------------------------
   proc 0 reads N lines from file
   could be skipping Natoms lines, so use bigints
------------------------------------------------------------------------- */

void ReadData::skip_lines(bigint n)
{
  if (me) return;
  if (n <= 0) return;
  char *eof = NULL;
  for (bigint i = 0; i < n; i++) eof = fgets(line,MAXLINE,fp);
  if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
}

/* ----------------------------------------------------------------------
   parse a line of coeffs into words, storing them in narg,arg
   trim anything from '#' onward
   word strings remain in line, are not copied
   if addstr != NULL, add addstr as extra arg for class2 angle/dihedral/improper
     if 2nd word starts with letter, then is hybrid style, add addstr after it
     else add addstr before 2nd word
   if dupflag, duplicate 1st word, so pair_coeff "2" becomes "2 2"
   if noffset, add offset to first noffset args, which are atom/bond/etc types
------------------------------------------------------------------------- */

void ReadData::parse_coeffs(char *line, const char *addstr,
                            int dupflag, int noffset, int offset)
{
  char *ptr;
  if ((ptr = strchr(line,'#'))) *ptr = '\0';

  narg = 0;
  char *word = strtok(line," \t\n\r\f");
  while (word) {
    if (narg == maxarg) {
      maxarg += DELTA;
      arg = (char **)
        memory->srealloc(arg,maxarg*sizeof(char *),"read_data:arg");
    }
    if (addstr && narg == 1 && !islower(word[0])) arg[narg++] = (char *) addstr;
    arg[narg++] = word;
    if (addstr && narg == 2 && islower(word[0])) arg[narg++] = (char *) addstr;
    if (dupflag && narg == 1) arg[narg++] = word;
    word = strtok(NULL," \t\n\r\f");
  }

  // to avoid segfaults on empty lines

  if (narg == 0) return;
  
  if (noffset) {
    int value = force->inumeric(FLERR,arg[0]);
    sprintf(argoffset1,"%d",value+offset);
    arg[0] = argoffset1;
    if (noffset == 2) {
      value = force->inumeric(FLERR,arg[1]);
      sprintf(argoffset2,"%d",value+offset);
      arg[1] = argoffset2;
    }
  }
}

/* ----------------------------------------------------------------------
   compare two style strings if they both exist
   one = comment in data file section, two = currently-defined style
   ignore suffixes listed in suffixes array at top of file
------------------------------------------------------------------------- */

int ReadData::style_match(const char *one, const char *two)
{
  int i,delta,len,len1,len2;

  if ((one == NULL) || (two == NULL)) return 1;

  len1 = strlen(one);
  len2 = strlen(two);

  for (i = 0; suffixes[i] != NULL; i++) {
    len = strlen(suffixes[i]);
    if ((delta = len1 - len) > 0)
      if (strcmp(one+delta,suffixes[i]) == 0) len1 = delta;
    if ((delta = len2 - len) > 0)
      if (strcmp(two+delta,suffixes[i]) == 0) len2 = delta;
  }

  if ((len1 == 0) || (len1 == len2) || (strncmp(one,two,len1) == 0)) return 1;
  return 0;
}
