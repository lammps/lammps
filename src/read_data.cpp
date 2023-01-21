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

#include "read_data.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_body.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "bond.h"
#include "comm.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "irregular.h"
#include "label_map.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "pair.h"
#include "special.h"
#include "tokenizer.h"
#include "update.h"

#include <cctype>
#include <cstring>
#include <unordered_map>
#include <unordered_set>

using namespace LAMMPS_NS;

static constexpr int MAXLINE = 256;
static constexpr double LB_FACTOR = 1.1;
static constexpr int CHUNK = 1024;
static constexpr int DELTA = 4;       // must be 2 or larger
static constexpr int MAXBODY = 32;    // max # of lines in one body

// customize for new sections

// clang-format off
static std::unordered_set<std::string> section_keywords = {
  "Atoms", "Velocities", "Ellipsoids", "Lines", "Triangles", "Bodies",
  "Bonds", "Angles", "Dihedrals", "Impropers",
  "Masses", "Pair Coeffs", "PairIJ Coeffs", "Bond Coeffs", "Angle Coeffs",
  "Dihedral Coeffs", "Improper Coeffs",
  "BondBond Coeffs", "BondAngle Coeffs", "MiddleBondTorsion Coeffs",
  "EndBondTorsion Coeffs", "AngleTorsion Coeffs",
  "AngleAngleTorsion Coeffs", "BondBond13 Coeffs", "AngleAngle Coeffs",
  "Atom Type Labels", "Bond Type Labels", "Angle Type Labels",
  "Dihedral Type Labels", "Improper Type Labels"
};

// function to check whether a string is a known data section name
// made a static class member, so it can be called from other classes

bool ReadData::is_data_section(const std::string &keyword)
{
  return section_keywords.count(keyword) > 0;
}

enum{NONE, APPEND, VALUE, MERGE};

// pair style suffixes to ignore
// when matching Pair Coeffs comment to currently-defined pair style

static const char *suffixes[] = {"/cuda", "/gpu", "/opt", "/omp", "/kk", "/coul/cut", "/coul/long",
                                 "/coul/msm", "/coul/dsf", "/coul/debye", "/coul/charmm", nullptr};

static const char *labeltypes[] = {"Atom", "Bond", "Angle", "Dihedral", "Improper" };

// clang-format on
/* ---------------------------------------------------------------------- */
ReadData::ReadData(LAMMPS *_lmp) : Command(_lmp), fp(nullptr), coeffarg(nullptr), lmap(nullptr)
{
  MPI_Comm_rank(world, &me);
  line = new char[MAXLINE];
  keyword = new char[MAXLINE];
  style = new char[MAXLINE];
  buffer = new char[CHUNK * MAXLINE];
  ncoeffarg = maxcoeffarg = 0;

  // customize for new sections
  // pointers to atom styles that store bonus info

  nellipsoids = 0;
  avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  nlines = 0;
  avec_line = dynamic_cast<AtomVecLine *>(atom->style_match("line"));
  ntris = 0;
  avec_tri = dynamic_cast<AtomVecTri *>(atom->style_match("tri"));
  nbodies = 0;
  avec_body = dynamic_cast<AtomVecBody *>(atom->style_match("body"));
}

/* ---------------------------------------------------------------------- */

ReadData::~ReadData()
{
  delete[] line;
  delete[] keyword;
  delete[] style;
  delete[] buffer;
  memory->sfree(coeffarg);

  for (int i = 0; i < nfix; i++) {
    delete[] fix_header[i];
    delete[] fix_section[i];
  }
  memory->sfree(fix_index);
  memory->sfree(fix_header);
  memory->sfree(fix_section);

  delete lmap;
}

// clang format off

/* ---------------------------------------------------------------------- */

void ReadData::command(int narg, char **arg)
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "read_data", error);

  MPI_Barrier(world);
  double time1 = platform::walltime();

  // optional args

  addflag = NONE;
  coeffflag = 1;
  id_offset = mol_offset = 0;
  offsetflag = shiftflag = settypeflag = 0;
  tlabelflag = blabelflag = alabelflag = dlabelflag = ilabelflag = 0;
  toffset = boffset = aoffset = doffset = ioffset = 0;
  shift[0] = shift[1] = shift[2] = 0.0;
  extra_atom_types = extra_bond_types = extra_angle_types = extra_dihedral_types =
      extra_improper_types = 0;

  groupbit = 0;

  nfix = 0;
  fix_index = nullptr;
  fix_header = nullptr;
  fix_section = nullptr;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "add") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "read_data add", error);
      if (strcmp(arg[iarg + 1], "append") == 0)
        addflag = APPEND;
      else if (strcmp(arg[iarg + 1], "merge") == 0)
        addflag = MERGE;
      else {
        if (atom->molecule_flag && (iarg + 3 > narg))
          utils::missing_cmd_args(FLERR, "read_data add", error);
        addflag = VALUE;
        bigint offset = utils::bnumeric(FLERR, arg[iarg + 1], false, lmp);
        if (offset > MAXTAGINT)
          error->all(FLERR, "Read data add atomID offset {} is too big", offset);
        id_offset = offset;

        if (atom->molecule_flag) {
          offset = utils::bnumeric(FLERR, arg[iarg + 2], false, lmp);
          if (offset > MAXTAGINT)
            error->all(FLERR, "Read data add molID offset {} is too big", offset);
          mol_offset = offset;
          iarg++;
        }
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "offset") == 0) {
      if (iarg + 6 > narg) utils::missing_cmd_args(FLERR, "read_data offset", error);
      offsetflag = 1;
      toffset = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      boffset = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
      aoffset = utils::inumeric(FLERR, arg[iarg + 3], false, lmp);
      doffset = utils::inumeric(FLERR, arg[iarg + 4], false, lmp);
      ioffset = utils::inumeric(FLERR, arg[iarg + 5], false, lmp);
      if (toffset < 0 || boffset < 0 || aoffset < 0 || doffset < 0 || ioffset < 0)
        error->all(FLERR, "Illegal read_data offset value(s)");
      iarg += 6;
    } else if (strcmp(arg[iarg], "shift") == 0) {
      if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "read_data shift", error);
      shiftflag = 1;
      shift[0] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      shift[1] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      shift[2] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      if (domain->dimension == 2 && shift[2] != 0.0)
        error->all(FLERR, "Non-zero read_data shift z value for 2d simulation not allowed");
      iarg += 4;
    } else if (strcmp(arg[iarg], "nocoeff") == 0) {
      coeffflag = 0;
      iarg++;
    } else if (strcmp(arg[iarg], "extra/atom/types") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "read_data extra/atom/types", error);
      extra_atom_types = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (extra_atom_types < 0)
        error->all(FLERR, "Illegal read_data extra/atom/types value {}", extra_atom_types);
      iarg += 2;
    } else if (strcmp(arg[iarg], "extra/bond/types") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "read_data extra/bond/types", error);
      if (!atom->avec->bonds_allow)
        error->all(FLERR, "No bonds allowed with atom style {}", atom->get_style());
      extra_bond_types = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (extra_bond_types < 0)
        error->all(FLERR, "Illegal read_data extra/bond/types value {}", extra_bond_types);
      iarg += 2;
    } else if (strcmp(arg[iarg], "extra/angle/types") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "read_data extra/angle/types", error);
      if (!atom->avec->angles_allow)
        error->all(FLERR, "No angles allowed with atom style {}", atom->get_style());
      extra_angle_types = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (extra_angle_types < 0)
        error->all(FLERR, "Illegal read_data extra/angle/types value {}", extra_angle_types);
      iarg += 2;
    } else if (strcmp(arg[iarg], "extra/dihedral/types") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "read_data extra/dihedral/types", error);
      if (!atom->avec->dihedrals_allow)
        error->all(FLERR, "No dihedrals allowed with atom style {}", atom->get_style());
      extra_dihedral_types = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (extra_dihedral_types < 0)
        error->all(FLERR, "Illegal read_data extra/dihedral/types value {}", extra_dihedral_types);
      iarg += 2;
    } else if (strcmp(arg[iarg], "extra/improper/types") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "read_data extra/improper/types", error);
      if (!atom->avec->impropers_allow)
        error->all(FLERR, "No impropers allowed with atom style {}", atom->get_style());
      extra_improper_types = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (extra_improper_types < 0)
        error->all(FLERR, "Illegal read_data extra/improper/types value {}", extra_improper_types);
      iarg += 2;
    } else if (strcmp(arg[iarg], "extra/bond/per/atom") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "read_data extra/bond/per/atom", error);
      if (atom->molecular == Atom::ATOMIC)
        error->all(FLERR, "No bonds allowed with atom style {}", atom->get_style());
      atom->extra_bond_per_atom = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (atom->extra_bond_per_atom < 0)
        error->all(FLERR, "Illegal read_data extra/bond/per/atom value {}",
                   atom->extra_bond_per_atom);
      iarg += 2;
    } else if (strcmp(arg[iarg], "extra/angle/per/atom") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "read_data extra/angle/per/atom", error);
      if (atom->molecular == Atom::ATOMIC)
        error->all(FLERR, "No angles allowed with atom style {}", atom->get_style());
      atom->extra_angle_per_atom = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (atom->extra_angle_per_atom < 0)
        error->all(FLERR, "Illegal read_data extra/angle/per/atom value {}",
                   atom->extra_angle_per_atom);
      iarg += 2;
    } else if (strcmp(arg[iarg], "extra/dihedral/per/atom") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, "read_data extra/dihedral/per/atom", error);
      if (atom->molecular == Atom::ATOMIC)
        error->all(FLERR, "No dihedrals allowed with atom style {}", atom->get_style());
      atom->extra_dihedral_per_atom = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (atom->extra_dihedral_per_atom < 0)
        error->all(FLERR, "Illegal read_data extra/dihedral/per/atom value {}",
                   atom->extra_dihedral_per_atom);
      iarg += 2;
    } else if (strcmp(arg[iarg], "extra/improper/per/atom") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, "read_data extra/improper/per/atom", error);
      if (atom->molecular == Atom::ATOMIC)
        error->all(FLERR, "No impropers allowed with atom style {}", atom->get_style());
      atom->extra_improper_per_atom = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (atom->extra_improper_per_atom < 0)
        error->all(FLERR, "Illegal read_data extra/improper/per/atom value {}",
                   atom->extra_improper_per_atom);
      iarg += 2;
    } else if (strcmp(arg[iarg], "extra/special/per/atom") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, "read_data extra/special/per/atom", error);
      if (atom->molecular == Atom::ATOMIC)
        error->all(FLERR, "No bonded interactions allowed with atom style {}", atom->get_style());
      force->special_extra = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (force->special_extra < 0)
        error->all(FLERR, "Illegal read_data extra/special/per/atom value {}",
                   force->special_extra);
      iarg += 2;
    } else if (strcmp(arg[iarg], "group") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "read_data group", error);
      int igroup = group->find_or_create(arg[iarg + 1]);
      groupbit = group->bitmask[igroup];
      iarg += 2;
    } else if (strcmp(arg[iarg], "fix") == 0) {
      if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "read_data fix", error);
      fix_index =
          (Fix **) memory->srealloc(fix_index, (nfix + 1) * sizeof(Fix *), "read_data:fix_index");
      fix_header = (char **) memory->srealloc(fix_header, (nfix + 1) * sizeof(char *),
                                              "read_data:fix_header");
      fix_section = (char **) memory->srealloc(fix_section, (nfix + 1) * sizeof(char *),
                                               "read_data:fix_section");
      if (is_data_section(arg[iarg + 3]))
        error->all(FLERR,
                   "Custom data section name {} for fix {} collides with existing data section",
                   arg[iarg + 3], arg[iarg + 1]);
      fix_index[nfix] = modify->get_fix_by_id(arg[iarg + 1]);
      if (!fix_index[nfix])
        error->all(FLERR, "Fix ID {} for read_data does not exist", arg[iarg + 1]);
      if (strcmp(arg[iarg + 2], "NULL") == 0)
        fix_header[nfix] = nullptr;
      else
        fix_header[nfix] = utils::strdup(arg[iarg + 2]);
      if (strcmp(arg[iarg + 3], "NULL") == 0)
        fix_section[nfix] = utils::strdup(arg[iarg + 1]);
      else
        fix_section[nfix] = utils::strdup(arg[iarg + 3]);
      nfix++;
      iarg += 4;

    } else
      error->all(FLERR, "Unknown read_data keyword {}", arg[iarg]);
  }

  // error checks

  if ((domain->dimension == 2) && (domain->zperiodic == 0))
    error->all(FLERR, "Cannot run 2d simulation with nonperiodic Z dimension");
  if ((domain->nonperiodic == 2) && utils::strmatch(force->kspace_style, "^msm"))
    error->all(FLERR,
               "Reading a data file with shrinkwrap boundaries is "
               "not compatible with a MSM KSpace style");
  if (domain->box_exist && !addflag)
    error->all(FLERR, "Cannot use read_data without add keyword after simulation box is defined");
  if (!domain->box_exist && addflag)
    error->all(FLERR, "Cannot use read_data add before simulation box is defined");
  if (offsetflag) {
    if (addflag == NONE) {
      error->all(FLERR, "Cannot use read_data offset without add keyword");
    } else {
      if (atom->labelmapflag) {
        if (comm->me == 0)
          error->warning(FLERR,
                         "Using read_data offset with a labelmap. Offsets will be only "
                         "applied to numeric types and not to type labels");
      }
    }
  }
  if (shiftflag && addflag == NONE)
    error->all(FLERR, "Cannot use read_data shift without add keyword");
  if (addflag != NONE &&
      (extra_atom_types || extra_bond_types || extra_angle_types || extra_dihedral_types ||
       extra_improper_types))
    error->all(FLERR, "Cannot use any read_data extra/*/types keyword with add keyword");

  // check if data file is available and readable

  if (!platform::file_is_readable(arg[0]))
    error->all(FLERR, fmt::format("Cannot open file {}: {}", arg[0], utils::getsyserror()));

  // reset so we can warn about reset image flags exactly once per data file

  atom->reset_image_flag[0] = atom->reset_image_flag[1] = atom->reset_image_flag[2] = false;

  // first time system initialization

  if (addflag == NONE) {
    domain->box_exist = 1;
    update->ntimestep = 0;
  } else {

    // clear global->local map for owned and ghost atoms
    // clear ghost count and any ghost bonus data internal to AtomVec
    // same logic as beginning of Comm::exchange()
    // do it now b/c adding atoms will overwrite ghost atoms

    if (atom->map_style != Atom::MAP_NONE) atom->map_clear();
    atom->nghost = 0;
    atom->avec->clear_bonus();
  }

  // compute atomID and optionally moleculeID offset for addflag = APPEND

  if (addflag == APPEND) {
    tagint *tag = atom->tag;
    tagint *molecule = atom->molecule;
    int nlocal = atom->nlocal;
    tagint maxid = 0, maxmol = 0;
    for (int i = 0; i < nlocal; i++) maxid = MAX(maxid, tag[i]);
    if (atom->molecule_flag)
      for (int i = 0; i < nlocal; i++) maxmol = MAX(maxmol, molecule[i]);
    MPI_Allreduce(&maxid, &id_offset, 1, MPI_LMP_TAGINT, MPI_MAX, world);
    MPI_Allreduce(&maxmol, &mol_offset, 1, MPI_LMP_TAGINT, MPI_MAX, world);
  }

  // set up pointer to hold original styles while we replace them with "zero"

  Pair *saved_pair = nullptr;
  Bond *saved_bond = nullptr;
  Angle *saved_angle = nullptr;
  Dihedral *saved_dihedral = nullptr;
  Improper *saved_improper = nullptr;
  KSpace *saved_kspace = nullptr;
  char *saved_pair_style = nullptr;
  char *saved_bond_style = nullptr;
  char *saved_angle_style = nullptr;
  char *saved_dihedral_style = nullptr;
  char *saved_improper_style = nullptr;
  char *saved_kspace_style = nullptr;

  if (coeffflag == 0) {
    char *coeffs[2];
    coeffs[0] = (char *) "10.0";
    coeffs[1] = (char *) "nocoeff";

    saved_pair = force->pair;
    saved_pair_style = force->pair_style;
    force->pair = nullptr;
    force->pair_style = nullptr;
    force->create_pair("zero", 0);
    if (force->pair) force->pair->settings(2, coeffs);

    coeffs[0] = coeffs[1];
    saved_bond = force->bond;
    saved_bond_style = force->bond_style;
    force->bond = nullptr;
    force->bond_style = nullptr;
    force->create_bond("zero", 0);
    if (force->bond) force->bond->settings(1, coeffs);

    saved_angle = force->angle;
    saved_angle_style = force->angle_style;
    force->angle = nullptr;
    force->angle_style = nullptr;
    force->create_angle("zero", 0);
    if (force->angle) force->angle->settings(1, coeffs);

    saved_dihedral = force->dihedral;
    saved_dihedral_style = force->dihedral_style;
    force->dihedral = nullptr;
    force->dihedral_style = nullptr;
    force->create_dihedral("zero", 0);
    if (force->dihedral) force->dihedral->settings(1, coeffs);

    saved_improper = force->improper;
    saved_improper_style = force->improper_style;
    force->improper = nullptr;
    force->improper_style = nullptr;
    force->create_improper("zero", 0);
    if (force->improper) force->improper->settings(1, coeffs);

    saved_kspace = force->kspace;
    saved_kspace_style = force->kspace_style;
    force->kspace = nullptr;
    force->kspace_style = nullptr;
  }

  // -----------------------------------------------------------------

  // perform 1-pass read if no molecular topology in file
  // perform 2-pass read if molecular topology,
  //   first pass calculates max topology/atom

  // flags for this data file

  int atomflag, topoflag;
  int bondflag, angleflag, dihedralflag, improperflag;
  int ellipsoidflag, lineflag, triflag, bodyflag;

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

  while (true) {

    // open file on proc 0

    if (me == 0) {
      if (firstpass) utils::logmesg(lmp, "Reading data file ...\n");
      open(arg[0]);
    } else
      fp = nullptr;

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
      if (comm->nprocs == 1)
        n = static_cast<int>(atom->natoms);
      else
        n = static_cast<int>(LB_FACTOR * atom->natoms / comm->nprocs);

      atom->allocate_type_arrays();
      atom->deallocate_topology();

      // allocate atom arrays to N, rounded up by AtomVec->DELTA

      bigint nbig = n;
      nbig = atom->avec->roundup(nbig);
      n = static_cast<int>(nbig);
      atom->avec->grow(n);

      domain->boxlo[0] = boxlo[0];
      domain->boxhi[0] = boxhi[0];
      domain->boxlo[1] = boxlo[1];
      domain->boxhi[1] = boxhi[1];
      domain->boxlo[2] = boxlo[2];
      domain->boxhi[2] = boxhi[2];

      if (triclinic) {
        domain->triclinic = 1;
        domain->xy = xy;
        domain->xz = xz;
        domain->yz = yz;
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
      double oldboxlo[3] = { domain->boxlo[0],  domain->boxlo[1] , domain->boxlo[2]};
      double oldboxhi[3] = { domain->boxhi[0],  domain->boxhi[1] , domain->boxhi[2]};
      domain->boxlo[0] = MIN(domain->boxlo[0], boxlo[0] + shift[0]);
      domain->boxhi[0] = MAX(domain->boxhi[0], boxhi[0] + shift[0]);
      domain->boxlo[1] = MIN(domain->boxlo[1], boxlo[1] + shift[1]);
      domain->boxhi[1] = MAX(domain->boxhi[1], boxhi[1] + shift[1]);
      domain->boxlo[2] = MIN(domain->boxlo[2], boxlo[2] + shift[2]);
      domain->boxhi[2] = MAX(domain->boxhi[2], boxhi[2] + shift[2]);

      // check of box has changed. If yes, warn about non-zero image flags
      if ((oldboxlo[0] != domain->boxlo[0]) || (oldboxlo[1] != domain->boxlo[1]) ||
          (oldboxlo[2] != domain->boxlo[2]) || (oldboxhi[0] != domain->boxhi[0]) ||
          (oldboxhi[1] != domain->boxhi[1]) || (oldboxhi[2] != domain->boxhi[2])) {
        int iflag = 1;
        for (int i=0; i < atom->nlocal; ++i) {
          int xbox = (atom->image[i] & IMGMASK) - IMGMAX;
          int ybox = (atom->image[i] >> IMGBITS & IMGMASK) - IMGMAX;
          int zbox = (atom->image[i] >> IMG2BITS) - IMGMAX;
          if (xbox != 0) iflag = 1;
          if (ybox != 0) iflag = 1;
          if (zbox != 0) iflag = 1;
        }
        int flag_all;
        MPI_Allreduce(&iflag,&flag_all, 1, MPI_INT, MPI_SUM, world);
        if ((flag_all > 0) && (comm->me == 0))
          error->warning(FLERR,"Non-zero image flags with growing box leads to bad coordinates");
      }

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

    // allocate space for type label map

    if (firstpass) {
      delete lmap;
      lmap = new LabelMap(lmp, ntypes, nbondtypes, nangletypes, ndihedraltypes, nimpropertypes);
    }

    // customize for new sections
    // read rest of file in free format

    while (strlen(keyword)) {

      if (strcmp(keyword, "Atoms") == 0) {
        atomflag = 1;
        if (firstpass) {
          if (me == 0 && !style_match(style, atom->atom_style))
            error->warning(FLERR,
                           "Atom style in data file differs from currently defined atom style");
          atoms();
        } else
          skip_lines(natoms);

      } else if (strcmp(keyword, "Velocities") == 0) {
        if (atomflag == 0) error->all(FLERR, "Must read Atoms before Velocities");
        if (firstpass)
          velocities();
        else
          skip_lines(natoms);

      } else if (strcmp(keyword, "Bonds") == 0) {
        topoflag = bondflag = 1;
        if (nbonds == 0) error->all(FLERR, "Invalid data file section: Bonds");
        if (atomflag == 0) error->all(FLERR, "Must read Atoms before Bonds");
        bonds(firstpass);

      } else if (strcmp(keyword, "Angles") == 0) {
        topoflag = angleflag = 1;
        if (nangles == 0) error->all(FLERR, "Invalid data file section: Angles");
        if (atomflag == 0) error->all(FLERR, "Must read Atoms before Angles");
        angles(firstpass);

      } else if (strcmp(keyword, "Dihedrals") == 0) {
        topoflag = dihedralflag = 1;
        if (ndihedrals == 0) error->all(FLERR, "Invalid data file section: Dihedrals");
        if (atomflag == 0) error->all(FLERR, "Must read Atoms before Dihedrals");
        dihedrals(firstpass);

      } else if (strcmp(keyword, "Impropers") == 0) {
        topoflag = improperflag = 1;
        if (nimpropers == 0) error->all(FLERR, "Invalid data file section: Impropers");
        if (atomflag == 0) error->all(FLERR, "Must read Atoms before Impropers");
        impropers(firstpass);

      } else if (strcmp(keyword, "Ellipsoids") == 0) {
        ellipsoidflag = 1;
        if (!avec_ellipsoid) error->all(FLERR, "Invalid data file section: Ellipsoids");
        if (atomflag == 0) error->all(FLERR, "Must read Atoms before Ellipsoids");
        if (firstpass)
          bonus(nellipsoids, (AtomVec *) avec_ellipsoid, "ellipsoids");
        else
          skip_lines(nellipsoids);

      } else if (strcmp(keyword, "Lines") == 0) {
        lineflag = 1;
        if (!avec_line) error->all(FLERR, "Invalid data file section: Lines");
        if (atomflag == 0) error->all(FLERR, "Must read Atoms before Lines");
        if (firstpass)
          bonus(nlines, (AtomVec *) avec_line, "lines");
        else
          skip_lines(nlines);

      } else if (strcmp(keyword, "Triangles") == 0) {
        triflag = 1;
        if (!avec_tri) error->all(FLERR, "Invalid data file section: Triangles");
        if (atomflag == 0) error->all(FLERR, "Must read Atoms before Triangles");
        if (firstpass)
          bonus(ntris, (AtomVec *) avec_tri, "triangles");
        else
          skip_lines(ntris);

      } else if (strcmp(keyword, "Bodies") == 0) {
        bodyflag = 1;
        if (!avec_body) error->all(FLERR, "Invalid data file section: Bodies");
        if (atomflag == 0) error->all(FLERR, "Must read Atoms before Bodies");
        bodies(firstpass, (AtomVec *) avec_body);

      } else if (strcmp(keyword, "Masses") == 0) {
        if (firstpass)
          mass();
        else
          skip_lines(ntypes);
      } else if (strcmp(keyword, "Pair Coeffs") == 0) {
        if (force->pair == nullptr) error->all(FLERR, "Must define pair_style before Pair Coeffs");
        if (firstpass) {
          if (me == 0 && !style_match(style, force->pair_style))
            error->warning(FLERR,
                           "Pair style in data file differs from currently defined pair style");
          paircoeffs();
        } else
          skip_lines(ntypes);
      } else if (strcmp(keyword, "PairIJ Coeffs") == 0) {
        if (force->pair == nullptr)
          error->all(FLERR, "Must define pair_style before PairIJ Coeffs");
        if (firstpass) {
          if (me == 0 && !style_match(style, force->pair_style))
            error->warning(FLERR,
                           "Pair style in data file differs "
                           "from currently defined pair style");
          pairIJcoeffs();
        } else
          skip_lines(ntypes * (ntypes + 1) / 2);
      } else if (strcmp(keyword, "Bond Coeffs") == 0) {
        if (atom->avec->bonds_allow == 0)
          error->all(FLERR, "Invalid data file section: Bond Coeffs");
        if (force->bond == nullptr) error->all(FLERR, "Must define bond_style before Bond Coeffs");
        if (firstpass) {
          if (me == 0 && !style_match(style, force->bond_style))
            error->warning(FLERR,
                           "Bond style in data file differs from currently defined bond style");
          bondcoeffs();
        } else
          skip_lines(nbondtypes);
      } else if (strcmp(keyword, "Angle Coeffs") == 0) {
        if (atom->avec->angles_allow == 0)
          error->all(FLERR, "Invalid data file section: Angle Coeffs");
        if (force->angle == nullptr)
          error->all(FLERR, "Must define angle_style before Angle Coeffs");
        if (firstpass) {
          if (me == 0 && !style_match(style, force->angle_style))
            error->warning(FLERR,
                           "Angle style in data file differs from currently defined angle style");
          anglecoeffs(0);
        } else
          skip_lines(nangletypes);
      } else if (strcmp(keyword, "Dihedral Coeffs") == 0) {
        if (atom->avec->dihedrals_allow == 0)
          error->all(FLERR, "Invalid data file section: Dihedral Coeffs");
        if (force->dihedral == nullptr)
          error->all(FLERR, "Must define dihedral_style before Dihedral Coeffs");
        if (firstpass) {
          if (me == 0 && !style_match(style, force->dihedral_style))
            error->warning(FLERR,
                           "Dihedral style in data file differs "
                           "from currently defined dihedral style");
          dihedralcoeffs(0);
        } else
          skip_lines(ndihedraltypes);
      } else if (strcmp(keyword, "Improper Coeffs") == 0) {
        if (atom->avec->impropers_allow == 0)
          error->all(FLERR, "Invalid data file section: Improper Coeffs");
        if (force->improper == nullptr)
          error->all(FLERR, "Must define improper_style before Improper Coeffs");
        if (firstpass) {
          if (me == 0 && !style_match(style, force->improper_style))
            error->warning(FLERR,
                           "Improper style in data file differs "
                           "from currently defined improper style");
          impropercoeffs(0);
        } else
          skip_lines(nimpropertypes);

      } else if (strcmp(keyword, "BondBond Coeffs") == 0) {
        if (atom->avec->angles_allow == 0)
          error->all(FLERR, "Invalid data file section: BondBond Coeffs");
        if (force->angle == nullptr)
          error->all(FLERR, "Must define angle_style before BondBond Coeffs");
        if (firstpass)
          anglecoeffs(1);
        else
          skip_lines(nangletypes);

      } else if (strcmp(keyword, "BondAngle Coeffs") == 0) {
        if (atom->avec->angles_allow == 0)
          error->all(FLERR, "Invalid data file section: BondAngle Coeffs");
        if (force->angle == nullptr)
          error->all(FLERR, "Must define angle_style before BondAngle Coeffs");
        if (firstpass)
          anglecoeffs(2);
        else
          skip_lines(nangletypes);
      } else if (strcmp(keyword, "UreyBradley Coeffs") == 0) {
        if (atom->avec->angles_allow == 0)
          error->all(FLERR, "Invalid data file section: UreyBradley Coeffs");
        if (force->angle == nullptr)
          error->all(FLERR, "Must define angle_style before UreyBradley Coeffs");
        if (firstpass)
          anglecoeffs(3);
        else
          skip_lines(nangletypes);

      } else if (strcmp(keyword, "MiddleBondTorsion Coeffs") == 0) {
        if (atom->avec->dihedrals_allow == 0)
          error->all(FLERR, "Invalid data file section: MiddleBondTorsion Coeffs");
        if (force->dihedral == nullptr)
          error->all(FLERR, "Must define dihedral_style before MiddleBondTorsion Coeffs");
        if (firstpass)
          dihedralcoeffs(1);
        else
          skip_lines(ndihedraltypes);

      } else if (strcmp(keyword, "EndBondTorsion Coeffs") == 0) {
        if (atom->avec->dihedrals_allow == 0)
          error->all(FLERR, "Invalid data file section: EndBondTorsion Coeffs");
        if (force->dihedral == nullptr)
          error->all(FLERR, "Must define dihedral_style before EndBondTorsion Coeffs");
        if (firstpass)
          dihedralcoeffs(2);
        else
          skip_lines(ndihedraltypes);

      } else if (strcmp(keyword, "AngleTorsion Coeffs") == 0) {
        if (atom->avec->dihedrals_allow == 0)
          error->all(FLERR, "Invalid data file section: AngleTorsion Coeffs");
        if (force->dihedral == nullptr)
          error->all(FLERR, "Must define dihedral_style before AngleTorsion Coeffs");
        if (firstpass)
          dihedralcoeffs(3);
        else
          skip_lines(ndihedraltypes);

      } else if (strcmp(keyword, "AngleAngleTorsion Coeffs") == 0) {
        if (atom->avec->dihedrals_allow == 0)
          error->all(FLERR, "Invalid data file section: AngleAngleTorsion Coeffs");
        if (force->dihedral == nullptr)
          error->all(FLERR, "Must define dihedral_style before AngleAngleTorsion Coeffs");
        if (firstpass)
          dihedralcoeffs(4);
        else
          skip_lines(ndihedraltypes);

      } else if (strcmp(keyword, "BondBond13 Coeffs") == 0) {
        if (atom->avec->dihedrals_allow == 0)
          error->all(FLERR, "Invalid data file section: BondBond13 Coeffs");
        if (force->dihedral == nullptr)
          error->all(FLERR, "Must define dihedral_style before BondBond13 Coeffs");
        if (firstpass)
          dihedralcoeffs(5);
        else
          skip_lines(ndihedraltypes);

      } else if (strcmp(keyword, "AngleAngle Coeffs") == 0) {
        if (atom->avec->impropers_allow == 0)
          error->all(FLERR, "Invalid data file section: AngleAngle Coeffs");
        if (force->improper == nullptr)
          error->all(FLERR, "Must define improper_style before AngleAngle Coeffs");
        if (firstpass)
          impropercoeffs(1);
        else
          skip_lines(nimpropertypes);

      } else if (strcmp(keyword, "Atom Type Labels") == 0) {
        if (firstpass) {
          if (atomflag == 1) error->all(FLERR, "Must read Atom Type Labels before Atoms");
          tlabelflag = 1;
          typelabels(Atom::ATOM);
        } else
          skip_lines(ntypes);

      } else if (strcmp(keyword, "Bond Type Labels") == 0) {
        if (nbondtypes) {
          if (firstpass) {
            if (bondflag == 1) error->all(FLERR, "Must read Bond Type Labels before Bonds");
            blabelflag = 1;
            typelabels(Atom::BOND);
          } else
            skip_lines(nbondtypes);
        }

      } else if (strcmp(keyword, "Angle Type Labels") == 0) {
        if (nangletypes) {
          if (firstpass) {
            if (angleflag == 1) error->all(FLERR, "Must read Angle Type Labels before Angles");
            alabelflag = 1;
            typelabels(Atom::ANGLE);
          } else
            skip_lines(nangletypes);
        }

      } else if (strcmp(keyword, "Dihedral Type Labels") == 0) {
        if (ndihedraltypes) {
          if (firstpass) {
            if (dihedralflag == 1)
              error->all(FLERR, "Must read Dihedral Type Labels before Dihedrals");
            dlabelflag = 1;
            typelabels(Atom::DIHEDRAL);
          } else
            skip_lines(ndihedraltypes);
        }

      } else if (strcmp(keyword, "Improper Type Labels") == 0) {
        if (nimpropertypes) {
          if (firstpass) {
            if (improperflag == 1)
              error->all(FLERR, "Must read Improper Type Labels before Impropers");
            ilabelflag = 1;
            typelabels(Atom::IMPROPER);
          } else
            skip_lines(nimpropertypes);
        }

        // if specified fix matches, it processes section

      } else if (nfix) {
        int i;
        for (i = 0; i < nfix; i++)
          if (strcmp(keyword, fix_section[i]) == 0) {
            if (firstpass) {
              fix(fix_index[i], keyword);
            } else {
              auto nskip = fix_index[i]->read_data_skip_lines(keyword);
              if (nskip < 0) nskip = natoms;
              skip_lines(nskip);
            }
            break;
          }
        if (i == nfix)
          error->all(FLERR, "Unknown identifier in data file: {}{}", keyword, utils::errorurl(1));

      } else
        error->all(FLERR, "Unknown identifier in data file: {}{}", keyword, utils::errorurl(1));

      parse_keyword(0);
    }

    // error if natoms > 0 yet no atoms were read

    if (natoms > 0 && atomflag == 0) error->all(FLERR, "No valid atoms found in data file");

    // close file

    if (me == 0) {
      if (compressed)
        platform::pclose(fp);
      else
        fclose(fp);
      fp = nullptr;
    }

    // done if this was 2nd pass

    if (!firstpass) break;

    // at end of 1st pass, error check for required sections
    // customize for new sections

    if (nbonds && !bondflag) error->one(FLERR, "Bonds section for {} bonds not found", nbonds);
    if (nangles && !angleflag) error->one(FLERR, "Angles section for {} angles not found", nangles);
    if (ndihedrals && !dihedralflag)
      error->one(FLERR, "Dihedrals section for {} dihedrals not found", ndihedrals);
    if (nimpropers && !improperflag)
      error->one(FLERR, "Impropers section for {} impropers not found", nimpropers);

    if (nellipsoids && !ellipsoidflag)
      error->one(FLERR, "Ellipsoids section for {} ellipsoids not found", nellipsoids);
    if (nlines && !lineflag) error->one(FLERR, "Lines section for {} lines not found", nlines);
    if (ntris && !triflag) error->one(FLERR, "Triangles section for {} triangles not found", ntris);
    if (nbodies && !bodyflag) error->one(FLERR, "Bodies section for {} bodies not found", nbodies);

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

  atom->data_fix_compute_variable(nlocal_previous, atom->nlocal);

  // assign atoms added by this data file to specified group

  if (groupbit) {
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    for (int i = nlocal_previous; i < nlocal; i++) mask[i] |= groupbit;
  }

  // create special bond lists for molecular systems

  if (atom->molecular == Atom::MOLECULAR) {
    Special special(lmp);
    special.build();
  }

  // for atom style template just count total bonds, etc. from template(s)

  if (atom->molecular == Atom::TEMPLATE) {
    Molecule **onemols = atom->avec->onemols;
    int *molindex = atom->molindex;
    int *molatom = atom->molatom;
    int nlocal = atom->nlocal;

    int imol, iatom;
    bigint nbonds, nangles, ndihedrals, nimpropers;
    nbonds = nangles = ndihedrals = nimpropers = 0;

    for (int i = 0; i < nlocal; i++) {
      imol = molindex[i];
      iatom = molatom[i];
      if (imol >= 0) {
        nbonds += onemols[imol]->num_bond[iatom];
        nangles += onemols[imol]->num_angle[iatom];
        ndihedrals += onemols[imol]->num_dihedral[iatom];
        nimpropers += onemols[imol]->num_improper[iatom];
      }
    }

    MPI_Allreduce(&nbonds, &atom->nbonds, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    MPI_Allreduce(&nangles, &atom->nangles, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    MPI_Allreduce(&ndihedrals, &atom->ndihedrals, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    MPI_Allreduce(&nimpropers, &atom->nimpropers, 1, MPI_LMP_BIGINT, MPI_SUM, world);

    if (me == 0) {
      std::string mesg;

      if (atom->nbonds) mesg += fmt::format("  {} template bonds\n", atom->nbonds);
      if (atom->nangles) mesg += fmt::format("  {} template angles\n", atom->nangles);
      if (atom->ndihedrals) mesg += fmt::format("  {} template dihedrals\n", atom->ndihedrals);
      if (atom->nimpropers) mesg += fmt::format("  {} template impropers\n", atom->nimpropers);

      utils::logmesg(lmp, mesg);
    }
  }

  // for atom style template systems
  // insure nbondtypes,etc are still consistent with template molecules,
  //   in case data file re-defined them

  if (atom->molecular == Atom::TEMPLATE) {
    int nset = MAX(1, atom->avec->onemols[0]->nset);
    for (int i = 0; i < nset; ++i) atom->avec->onemols[i]->check_attributes();
  }

  // if adding atoms, migrate atoms to new processors
  // use irregular() b/c box size could have changed dramaticaly
  // resulting in procs now owning very different subboxes
  // with their previously owned atoms now far outside the subbox

  if (addflag != NONE) {
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    auto irregular = new Irregular(lmp);
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
    auto irregular = new Irregular(lmp);
    irregular->migrate_atoms(1);
    delete irregular;
    if (domain->triclinic) domain->lamda2x(atom->nlocal);

    bigint natoms;
    bigint nblocal = atom->nlocal;
    MPI_Allreduce(&nblocal, &natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    if (natoms != atom->natoms)
      error->all(FLERR, "Read_data shrink wrap did not assign all atoms correctly");
  }

  // restore old styles, when reading with nocoeff flag given

  if (coeffflag == 0) {
    delete force->pair;
    delete[] force->pair_style;
    force->pair = saved_pair;
    force->pair_style = saved_pair_style;

    delete force->bond;
    delete[] force->bond_style;
    force->bond = saved_bond;
    force->bond_style = saved_bond_style;

    delete force->angle;
    delete[] force->angle_style;
    force->angle = saved_angle;
    force->angle_style = saved_angle_style;

    delete force->dihedral;
    delete[] force->dihedral_style;
    force->dihedral = saved_dihedral;
    force->dihedral_style = saved_dihedral_style;

    delete force->improper;
    delete[] force->improper_style;
    force->improper = saved_improper;
    force->improper_style = saved_improper_style;

    force->kspace = saved_kspace;
    force->kspace_style = saved_kspace_style;
  }

  // total time

  MPI_Barrier(world);

  if (comm->me == 0)
    utils::logmesg(lmp, "  read_data CPU = {:.3f} seconds\n", platform::walltime() - time1);
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

  // initialize type counts by the "extra" numbers so they get counted
  // in case the corresponding "types" line is missing and thus the extra
  // value will not be processed.
  if (addflag == NONE) {
    atom->ntypes = extra_atom_types;
    atom->nbondtypes = extra_bond_types;
    atom->nangletypes = extra_angle_types;
    atom->ndihedraltypes = extra_dihedral_types;
    atom->nimpropertypes = extra_improper_types;
  }

  // skip 1st line of file

  if (me == 0) {
    char *eof = utils::fgets_trunc(line, MAXLINE, fp);
    if (eof == nullptr) error->one(FLERR, "Unexpected end of data file");
  }

  while (true) {

    // read a line and bcast length

    if (me == 0) {
      if (utils::fgets_trunc(line, MAXLINE, fp) == nullptr)
        n = 0;
      else
        n = strlen(line) + 1;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, world);

    // if n = 0 then end-of-file so return with blank line

    if (n == 0) {
      line[0] = '\0';
      return;
    }

    MPI_Bcast(line, n, MPI_CHAR, 0, world);

    // trim anything from '#' onward
    // if line is blank, continue

    if ((ptr = strchr(line, '#'))) *ptr = '\0';
    if (strspn(line, " \t\n\r") == strlen(line)) continue;

    // allow special fixes first chance to match and process the line
    // if fix matches, continue to next header line

    if (nfix) {
      for (n = 0; n < nfix; n++) {
        if (!fix_header[n]) continue;
        if (strstr(line, fix_header[n])) {
          fix_index[n]->read_data_header(line);
          break;
        }
      }
      if (n < nfix) continue;
    }

    // search line for header keyword and set corresponding variable
    // customize for new header lines

    int extra_flag_value = 0;
    auto words = utils::split_words(line);

    if (utils::strmatch(line, "^\\s*\\d+\\s+atoms\\s")) {
      natoms = utils::bnumeric(FLERR, words[0], false, lmp);
      if (addflag == NONE)
        atom->natoms = natoms;
      else if (firstpass)
        atom->natoms += natoms;

    } else if (utils::strmatch(line, "^\\s*\\d+\\s+ellipsoids\\s")) {
      if (!avec_ellipsoid) error->all(FLERR, "No ellipsoids allowed with this atom style");
      nellipsoids = utils::bnumeric(FLERR, words[0], false, lmp);
      if (addflag == NONE)
        atom->nellipsoids = nellipsoids;
      else if (firstpass)
        atom->nellipsoids += nellipsoids;

    } else if (utils::strmatch(line, "^\\s*\\d+\\s+lines\\s")) {
      if (!avec_line) error->all(FLERR, "No lines allowed with this atom style");
      nlines = utils::bnumeric(FLERR, words[0], false, lmp);
      if (addflag == NONE)
        atom->nlines = nlines;
      else if (firstpass)
        atom->nlines += nlines;

    } else if (utils::strmatch(line, "^\\s*\\d+\\s+triangles\\s")) {
      if (!avec_tri) error->all(FLERR, "No triangles allowed with this atom style");
      ntris = utils::bnumeric(FLERR, words[0], false, lmp);
      if (addflag == NONE)
        atom->ntris = ntris;
      else if (firstpass)
        atom->ntris += ntris;

    } else if (utils::strmatch(line, "^\\s*\\d+\\s+bodies\\s")) {
      if (!avec_body) error->all(FLERR, "No bodies allowed with this atom style");
      nbodies = utils::bnumeric(FLERR, words[0], false, lmp);
      if (addflag == NONE)
        atom->nbodies = nbodies;
      else if (firstpass)
        atom->nbodies += nbodies;

    } else if (utils::strmatch(line, "^\\s*\\d+\\s+bonds\\s")) {
      nbonds = utils::bnumeric(FLERR, words[0], false, lmp);
      if (addflag == NONE)
        atom->nbonds = nbonds;
      else if (firstpass)
        atom->nbonds += nbonds;

    } else if (utils::strmatch(line, "^\\s*\\d+\\s+angles\\s")) {
      nangles = utils::bnumeric(FLERR, words[0], false, lmp);
      if (addflag == NONE)
        atom->nangles = nangles;
      else if (firstpass)
        atom->nangles += nangles;

    } else if (utils::strmatch(line, "^\\s*\\d+\\s+dihedrals\\s")) {
      ndihedrals = utils::bnumeric(FLERR, words[0], false, lmp);
      if (addflag == NONE)
        atom->ndihedrals = ndihedrals;
      else if (firstpass)
        atom->ndihedrals += ndihedrals;

    } else if (utils::strmatch(line, "^\\s*\\d+\\s+impropers\\s")) {
      nimpropers = utils::bnumeric(FLERR, words[0], false, lmp);
      if (addflag == NONE)
        atom->nimpropers = nimpropers;
      else if (firstpass)
        atom->nimpropers += nimpropers;

      // Atom class type settings are only set by first data file

    } else if (utils::strmatch(line, "^\\s*\\d+\\s+atom\\s+types\\s")) {
      ntypes = utils::inumeric(FLERR, words[0], false, lmp);
      if (addflag == NONE) atom->ntypes = ntypes + extra_atom_types;

    } else if (utils::strmatch(line, "\\s*\\d+\\s+bond\\s+types\\s")) {
      nbondtypes = utils::inumeric(FLERR, words[0], false, lmp);
      if (addflag == NONE) atom->nbondtypes = nbondtypes + extra_bond_types;

    } else if (utils::strmatch(line, "^\\s*\\d+\\s+angle\\s+types\\s")) {
      nangletypes = utils::inumeric(FLERR, words[0], false, lmp);
      if (addflag == NONE) atom->nangletypes = nangletypes + extra_angle_types;

    } else if (utils::strmatch(line, "^\\s*\\d+\\s+dihedral\\s+types\\s")) {
      ndihedraltypes = utils::inumeric(FLERR, words[0], false, lmp);
      if (addflag == NONE) atom->ndihedraltypes = ndihedraltypes + extra_dihedral_types;

    } else if (utils::strmatch(line, "^\\s*\\d+\\s+improper\\s+types\\s")) {
      nimpropertypes = utils::inumeric(FLERR, words[0], false, lmp);
      if (addflag == NONE) atom->nimpropertypes = nimpropertypes + extra_improper_types;

      // these settings only used by first data file
      // also, these are obsolescent. we parse them to maintain backward
      // compatibility, but the recommended way is to set them via keywords
      // in the LAMMPS input file. In case these flags are set in both,
      // the input and the data file, we use the larger of the two.

    } else if (strstr(line, "extra bond per atom")) {
      if (addflag == NONE) extra_flag_value = utils::inumeric(FLERR, words[0], false, lmp);
      atom->extra_bond_per_atom = MAX(atom->extra_bond_per_atom, extra_flag_value);
    } else if (strstr(line, "extra angle per atom")) {
      if (addflag == NONE) extra_flag_value = utils::inumeric(FLERR, words[0], false, lmp);
      atom->extra_angle_per_atom = MAX(atom->extra_angle_per_atom, extra_flag_value);
    } else if (strstr(line, "extra dihedral per atom")) {
      if (addflag == NONE) extra_flag_value = utils::inumeric(FLERR, words[0], false, lmp);
      atom->extra_dihedral_per_atom = MAX(atom->extra_dihedral_per_atom, extra_flag_value);
    } else if (strstr(line, "extra improper per atom")) {
      if (addflag == NONE) extra_flag_value = utils::inumeric(FLERR, words[0], false, lmp);
      atom->extra_improper_per_atom = MAX(atom->extra_improper_per_atom, extra_flag_value);
    } else if (strstr(line, "extra special per atom")) {
      if (addflag == NONE) extra_flag_value = utils::inumeric(FLERR, words[0], false, lmp);
      force->special_extra = MAX(force->special_extra, extra_flag_value);

      // local copy of box info
      // so can treat differently for first vs subsequent data files

    } else if (utils::strmatch(line, "^\\s*\\f+\\s+\\f+\\s+xlo\\s+xhi\\s")) {
      boxlo[0] = utils::numeric(FLERR, words[0], false, lmp);
      boxhi[0] = utils::numeric(FLERR, words[1], false, lmp);

    } else if (utils::strmatch(line, "^\\s*\\f+\\s+\\f+\\s+ylo\\s+yhi\\s")) {
      boxlo[1] = utils::numeric(FLERR, words[0], false, lmp);
      boxhi[1] = utils::numeric(FLERR, words[1], false, lmp);

    } else if (utils::strmatch(line, "^\\s*\\f+\\s+\\f+\\s+zlo\\s+zhi\\s")) {
      boxlo[2] = utils::numeric(FLERR, words[0], false, lmp);
      boxhi[2] = utils::numeric(FLERR, words[1], false, lmp);

    } else if (utils::strmatch(line, "^\\s*\\f+\\s+\\f+\\s+\\f+\\s+xy\\s+xz\\s+yz\\s")) {
      triclinic = 1;
      xy = utils::numeric(FLERR, words[0], false, lmp);
      xz = utils::numeric(FLERR, words[1], false, lmp);
      yz = utils::numeric(FLERR, words[2], false, lmp);

    } else
      break;
  }

  // error check on total system size

  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT || atom->nellipsoids < 0 ||
      atom->nellipsoids >= MAXBIGINT || atom->nlines < 0 || atom->nlines >= MAXBIGINT ||
      atom->ntris < 0 || atom->ntris >= MAXBIGINT || atom->nbodies < 0 ||
      atom->nbodies >= MAXBIGINT || atom->nbonds < 0 || atom->nbonds >= MAXBIGINT ||
      atom->nangles < 0 || atom->nangles >= MAXBIGINT || atom->ndihedrals < 0 ||
      atom->ndihedrals >= MAXBIGINT || atom->nimpropers < 0 || atom->nimpropers >= MAXBIGINT)
    error->all(FLERR, "System in data file is too big");

  // check that exiting string is a valid section keyword

  parse_keyword(1);
  if (!is_data_section(keyword)) error->all(FLERR, "Unknown identifier in data file: {}", keyword);

  // error checks on header values
  // must be consistent with atom style and other header values

  if ((atom->nbonds || atom->nbondtypes) && atom->avec->bonds_allow == 0)
    error->all(FLERR, "No bonds allowed with this atom style");
  if ((atom->nangles || atom->nangletypes) && atom->avec->angles_allow == 0)
    error->all(FLERR, "No angles allowed with this atom style");
  if ((atom->ndihedrals || atom->ndihedraltypes) && atom->avec->dihedrals_allow == 0)
    error->all(FLERR, "No dihedrals allowed with this atom style");
  if ((atom->nimpropers || atom->nimpropertypes) && atom->avec->impropers_allow == 0)
    error->all(FLERR, "No impropers allowed with this atom style");

  if (atom->nbonds > 0 && atom->nbondtypes <= 0)
    error->all(FLERR, "Bonds defined but no bond types");
  if (atom->nangles > 0 && atom->nangletypes <= 0)
    error->all(FLERR, "Angles defined but no angle types");
  if (atom->ndihedrals > 0 && atom->ndihedraltypes <= 0)
    error->all(FLERR, "Dihedrals defined but no dihedral types");
  if (atom->nimpropers > 0 && atom->nimpropertypes <= 0)
    error->all(FLERR, "Impropers defined but no improper types");

  if (atom->molecular == Atom::TEMPLATE) {
    if (atom->nbonds || atom->nangles || atom->ndihedrals || atom->nimpropers)
      error->all(FLERR, "No molecule topology allowed with atom style template");
  }
}

/* ----------------------------------------------------------------------
   read all atoms
------------------------------------------------------------------------- */

void ReadData::atoms()
{
  int nchunk, eof;

  if (me == 0) utils::logmesg(lmp, "  reading atoms ...\n");

  bigint nread = 0;

  while (nread < natoms) {
    nchunk = MIN(natoms - nread, CHUNK);
    eof = utils::read_lines_from_file(fp, nchunk, MAXLINE, buffer, me, world);
    if (eof) error->all(FLERR, "Unexpected end of data file");
    if (tlabelflag && !lmap->is_complete(Atom::ATOM))
      error->all(FLERR, "Label map is incomplete: all types must be assigned a unique type label");
    atom->data_atoms(nchunk, buffer, id_offset, mol_offset, toffset, shiftflag, shift, tlabelflag,
                     lmap->lmap2lmap.atom);
    nread += nchunk;
  }

  // warn if we have read data with non-zero image flags for non-periodic boundaries.
  // we may want to turn this into an error at some point, since this essentially
  // creates invalid position information that works by accident most of the time.

  if (comm->me == 0) {
    if (atom->reset_image_flag[0])
      error->warning(FLERR,
                     "Non-zero imageflag(s) in x direction for "
                     "non-periodic boundary reset to zero");
    if (atom->reset_image_flag[1])
      error->warning(FLERR,
                     "Non-zero imageflag(s) in y direction for "
                     "non-periodic boundary reset to zero");
    if (atom->reset_image_flag[2])
      error->warning(FLERR,
                     "Non-zero imageflag(s) in z direction for "
                     "non-periodic boundary reset to zero");
  }

  // check that all atoms were assigned correctly

  bigint n = atom->nlocal;
  bigint sum;
  MPI_Allreduce(&n, &sum, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  bigint nassign = sum - (atom->natoms - natoms);

  if (me == 0) utils::logmesg(lmp, "  {} atoms\n", nassign);

  if (sum != atom->natoms) error->all(FLERR, "Did not assign all atoms correctly");

  // check that atom IDs are valid

  atom->tag_check();

  // check that bonus data has been reserved as needed

  atom->bonus_check();

  // create global mapping of atoms

  if (atom->map_style != Atom::MAP_NONE) {
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
  int nchunk, eof;

  if (me == 0) utils::logmesg(lmp, "  reading velocities ...\n");

  int mapflag = 0;
  if (atom->map_style == Atom::MAP_NONE) {
    mapflag = 1;
    atom->map_init();
    atom->map_set();
  }

  bigint nread = 0;

  while (nread < natoms) {
    nchunk = MIN(natoms - nread, CHUNK);
    eof = utils::read_lines_from_file(fp, nchunk, MAXLINE, buffer, me, world);
    if (eof) error->all(FLERR, "Unexpected end of data file");
    atom->data_vels(nchunk, buffer, id_offset);
    nread += nchunk;
  }

  if (mapflag) {
    atom->map_delete();
    atom->map_style = Atom::MAP_NONE;
  }

  if (me == 0) utils::logmesg(lmp, "  {} velocities\n", natoms);
}

/* ----------------------------------------------------------------------
   scan or read all bonds
------------------------------------------------------------------------- */

void ReadData::bonds(int firstpass)
{
  int nchunk, eof;

  if (me == 0) {
    if (firstpass)
      utils::logmesg(lmp, "  scanning bonds ...\n");
    else
      utils::logmesg(lmp, "  reading bonds ...\n");
  }

  // allocate count if firstpass

  int nlocal = atom->nlocal;
  int *count = nullptr;
  if (firstpass) {
    memory->create(count, nlocal, "read_data:count");
    if (count) memset(count, 0, nlocal * sizeof(int));
  }

  // read and process bonds

  bigint nread = 0;

  while (nread < nbonds) {
    nchunk = MIN(nbonds - nread, CHUNK);
    eof = utils::read_lines_from_file(fp, nchunk, MAXLINE, buffer, me, world);
    if (eof) error->all(FLERR, "Unexpected end of data file");
    if (blabelflag && !lmap->is_complete(Atom::BOND))
      error->all(FLERR,
                 "Label map is incomplete: "
                 "all types must be assigned a unique type label");
    atom->data_bonds(nchunk, buffer, count, id_offset, boffset, blabelflag, lmap->lmap2lmap.bond);
    nread += nchunk;
  }

  // if firstpass: tally max bond/atom and return
  // if addflag = NONE, store max bond/atom with extra
  // else just check actual max does not exceed existing max

  if (firstpass) {
    int max = 0;
    for (int i = nlocal_previous; i < nlocal; i++) max = MAX(max, count[i]);
    int maxall;
    MPI_Allreduce(&max, &maxall, 1, MPI_INT, MPI_MAX, world);
    if (addflag == NONE) maxall += atom->extra_bond_per_atom;

    if (me == 0) utils::logmesg(lmp, "  {} = max bonds/atom\n", maxall);

    if (addflag != NONE) {
      if (maxall > atom->bond_per_atom)
        error->all(FLERR,
                   "Subsequent read data induced "
                   "too many bonds per atom");
    } else
      atom->bond_per_atom = maxall;

    memory->destroy(count);
    return;
  }

  // if 2nd pass: check that bonds were assigned correctly

  bigint n = 0;
  for (int i = nlocal_previous; i < nlocal; i++) n += atom->num_bond[i];
  bigint sum;
  MPI_Allreduce(&n, &sum, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  int factor = 1;
  if (!force->newton_bond) factor = 2;

  if (me == 0) utils::logmesg(lmp, "  {} bonds\n", sum / factor);

  if (sum != factor * nbonds) error->all(FLERR, "Bonds assigned incorrectly");
}

/* ----------------------------------------------------------------------
   scan or read all angles
------------------------------------------------------------------------- */

void ReadData::angles(int firstpass)
{
  int nchunk, eof;

  if (me == 0) {
    if (firstpass)
      utils::logmesg(lmp, "  scanning angles ...\n");
    else
      utils::logmesg(lmp, "  reading angles ...\n");
  }

  // allocate count if firstpass

  int nlocal = atom->nlocal;
  int *count = nullptr;
  if (firstpass) {
    memory->create(count, nlocal, "read_data:count");
    if (count) memset(count, 0, nlocal * sizeof(int));
  }

  // read and process angles

  bigint nread = 0;

  while (nread < nangles) {
    nchunk = MIN(nangles - nread, CHUNK);
    eof = utils::read_lines_from_file(fp, nchunk, MAXLINE, buffer, me, world);
    if (eof) error->all(FLERR, "Unexpected end of data file");
    if (alabelflag && !lmap->is_complete(Atom::ANGLE))
      error->all(FLERR,
                 "Label map is incomplete: "
                 "all types must be assigned a unique type label");
    atom->data_angles(nchunk, buffer, count, id_offset, aoffset, alabelflag, lmap->lmap2lmap.angle);
    nread += nchunk;
  }

  // if firstpass: tally max angle/atom and return
  // if addflag = NONE, store max angle/atom with extra
  // else just check actual max does not exceed existing max

  if (firstpass) {
    int max = 0;
    for (int i = nlocal_previous; i < nlocal; i++) max = MAX(max, count[i]);
    int maxall;
    MPI_Allreduce(&max, &maxall, 1, MPI_INT, MPI_MAX, world);
    if (addflag == NONE) maxall += atom->extra_angle_per_atom;

    if (me == 0) utils::logmesg(lmp, "  {} = max angles/atom\n", maxall);

    if (addflag != NONE) {
      if (maxall > atom->angle_per_atom)
        error->all(FLERR,
                   "Subsequent read data induced "
                   "too many angles per atom");
    } else
      atom->angle_per_atom = maxall;

    memory->destroy(count);
    return;
  }

  // if 2nd pass: check that angles were assigned correctly

  bigint n = 0;
  for (int i = nlocal_previous; i < nlocal; i++) n += atom->num_angle[i];
  bigint sum;
  MPI_Allreduce(&n, &sum, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  int factor = 1;
  if (!force->newton_bond) factor = 3;

  if (me == 0) utils::logmesg(lmp, "  {} angles\n", sum / factor);

  if (sum != factor * nangles) error->all(FLERR, "Angles assigned incorrectly");
}

/* ----------------------------------------------------------------------
   scan or read all dihedrals
------------------------------------------------------------------------- */

void ReadData::dihedrals(int firstpass)
{
  int nchunk, eof;

  if (me == 0) {
    if (firstpass)
      utils::logmesg(lmp, "  scanning dihedrals ...\n");
    else
      utils::logmesg(lmp, "  reading dihedrals ...\n");
  }

  // allocate count if firstpass

  int nlocal = atom->nlocal;
  int *count = nullptr;
  if (firstpass) {
    memory->create(count, nlocal, "read_data:count");
    if (count) memset(count, 0, nlocal * sizeof(int));
  }

  // read and process dihedrals

  bigint nread = 0;

  while (nread < ndihedrals) {
    nchunk = MIN(ndihedrals - nread, CHUNK);
    eof = utils::read_lines_from_file(fp, nchunk, MAXLINE, buffer, me, world);
    if (eof) error->all(FLERR, "Unexpected end of data file");
    if (dlabelflag && !lmap->is_complete(Atom::DIHEDRAL))
      error->all(FLERR,
                 "Label map is incomplete: "
                 "all types must be assigned a unique type label");
    atom->data_dihedrals(nchunk, buffer, count, id_offset, doffset, dlabelflag,
                         lmap->lmap2lmap.dihedral);
    nread += nchunk;
  }

  // if firstpass: tally max dihedral/atom and return
  // if addflag = NONE, store max dihedral/atom with extra
  // else just check actual max does not exceed existing max

  if (firstpass) {
    int max = 0;
    for (int i = nlocal_previous; i < nlocal; i++) max = MAX(max, count[i]);
    int maxall;
    MPI_Allreduce(&max, &maxall, 1, MPI_INT, MPI_MAX, world);
    if (addflag == NONE) maxall += atom->extra_dihedral_per_atom;

    if (me == 0) utils::logmesg(lmp, "  {} = max dihedrals/atom\n", maxall);

    if (addflag != NONE) {
      if (maxall > atom->dihedral_per_atom)
        error->all(FLERR,
                   "Subsequent read data induced "
                   "too many dihedrals per atom");
    } else
      atom->dihedral_per_atom = maxall;

    memory->destroy(count);
    return;
  }

  // if 2nd pass: check that dihedrals were assigned correctly

  bigint n = 0;
  for (int i = nlocal_previous; i < nlocal; i++) n += atom->num_dihedral[i];
  bigint sum;
  MPI_Allreduce(&n, &sum, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  int factor = 1;
  if (!force->newton_bond) factor = 4;

  if (me == 0) utils::logmesg(lmp, "  {} dihedrals\n", sum / factor);

  if (sum != factor * ndihedrals) error->all(FLERR, "Dihedrals assigned incorrectly");
}

/* ----------------------------------------------------------------------
   scan or read all impropers
------------------------------------------------------------------------- */

void ReadData::impropers(int firstpass)
{
  int nchunk, eof;

  if (me == 0) {
    if (firstpass)
      utils::logmesg(lmp, "  scanning impropers ...\n");
    else
      utils::logmesg(lmp, "  reading impropers ...\n");
  }

  // allocate count if firstpass

  int nlocal = atom->nlocal;
  int *count = nullptr;
  if (firstpass) {
    memory->create(count, nlocal, "read_data:count");
    if (count) memset(count, 0, nlocal * sizeof(int));
  }

  // read and process impropers

  bigint nread = 0;

  while (nread < nimpropers) {
    nchunk = MIN(nimpropers - nread, CHUNK);
    eof = utils::read_lines_from_file(fp, nchunk, MAXLINE, buffer, me, world);
    if (eof) error->all(FLERR, "Unexpected end of data file");
    if (ilabelflag && !lmap->is_complete(Atom::IMPROPER))
      error->all(FLERR,
                 "Label map is incomplete: "
                 "all types must be assigned a unique type label");
    atom->data_impropers(nchunk, buffer, count, id_offset, ioffset, ilabelflag,
                         lmap->lmap2lmap.improper);
    nread += nchunk;
  }

  // if firstpass: tally max improper/atom and return
  // if addflag = NONE, store max improper/atom
  // else just check it does not exceed existing max

  if (firstpass) {
    int max = 0;
    for (int i = nlocal_previous; i < nlocal; i++) max = MAX(max, count[i]);
    int maxall;
    MPI_Allreduce(&max, &maxall, 1, MPI_INT, MPI_MAX, world);
    if (addflag == NONE) maxall += atom->extra_improper_per_atom;

    if (me == 0) utils::logmesg(lmp, "  {} = max impropers/atom\n", maxall);

    if (addflag != NONE) {
      if (maxall > atom->improper_per_atom)
        error->all(FLERR, "Subsequent read data induced too many impropers per atom");
    } else
      atom->improper_per_atom = maxall;

    memory->destroy(count);
    return;
  }

  // if 2nd pass: check that impropers were assigned correctly

  bigint n = 0;
  for (int i = nlocal_previous; i < nlocal; i++) n += atom->num_improper[i];
  bigint sum;
  MPI_Allreduce(&n, &sum, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  int factor = 1;
  if (!force->newton_bond) factor = 4;

  if (me == 0) utils::logmesg(lmp, "  {} impropers\n", sum / factor);

  if (sum != factor * nimpropers) error->all(FLERR, "Impropers assigned incorrectly");
}

/* ----------------------------------------------------------------------
   read all bonus data
   to find atoms, must build atom map if not a molecular system
------------------------------------------------------------------------- */

void ReadData::bonus(bigint nbonus, AtomVec *ptr, const char *type)
{
  int nchunk, eof;

  int mapflag = 0;
  if (atom->map_style == Atom::MAP_NONE) {
    mapflag = 1;
    atom->map_init();
    atom->map_set();
  }

  bigint nread = 0;
  bigint natoms = nbonus;

  while (nread < natoms) {
    nchunk = MIN(natoms - nread, CHUNK);
    eof = utils::read_lines_from_file(fp, nchunk, MAXLINE, buffer, me, world);
    if (eof) error->all(FLERR, "Unexpected end of data file");
    atom->data_bonus(nchunk, buffer, ptr, id_offset);
    nread += nchunk;
  }

  if (mapflag) {
    atom->map_delete();
    atom->map_style = Atom::MAP_NONE;
  }

  if (me == 0) utils::logmesg(lmp, "  {} {}\n", natoms, type);
}

/* ----------------------------------------------------------------------
   read all body data
   variable amount of info per body, described by ninteger and ndouble
   to find atoms, must build atom map if not a molecular system
   if not firstpass, just read past body data and only process body header
------------------------------------------------------------------------- */

void ReadData::bodies(int firstpass, AtomVec *ptr)
{
  int m, nchunk, nline, nmax, ninteger, ndouble, nword, ncount, onebody;
  char *eof;

  int mapflag = 0;
  if (atom->map_style == Atom::MAP_NONE && firstpass) {
    mapflag = 1;
    atom->map_init();
    atom->map_set();
  }

  // nmax = max # of bodies to read in this chunk
  // nchunk = actual # read

  bigint nread = 0;
  bigint nblocks = nbodies;

  while (nread < nblocks) {
    if (nblocks - nread > CHUNK)
      nmax = CHUNK;
    else
      nmax = nblocks - nread;

    if (me == 0) {
      nchunk = 0;
      nline = 0;
      m = 0;

      while (nchunk < nmax && nline <= CHUNK - MAXBODY) {
        eof = utils::fgets_trunc(&buffer[m], MAXLINE, fp);
        const char *buf = &buffer[m];
        if (eof == nullptr) error->one(FLERR, "Unexpected end of data file");
        try {
          auto values = ValueTokenizer(utils::trim_comment(buf));
          tagint tagdata = values.next_tagint() + id_offset;
          ninteger = values.next_int();
          ndouble = values.next_double();
          if (tagdata <= 0 || tagdata > atom->map_tag_max)
            throw TokenizerException("Invalid atom ID in body header", utils::trim(buf));
          if (ninteger < 0)
            throw TokenizerException("Invalid number of integers", utils::trim(buf));
          if (ndouble < 0) throw TokenizerException("Invalid number of doubles", utils::trim(buf));
          if (values.has_next())
            throw TokenizerException("Too many tokens in body header", utils::trim(buf));
        } catch (TokenizerException &e) {
          error->one(FLERR, std::string(e.what()) + " while reading Bodies section of data file");
        }
        m += strlen(buf);

        // read lines one at a time into buffer and count words
        // count to ninteger and ndouble until have enough lines

        onebody = 0;

        nword = 0;
        while (nword < ninteger) {
          eof = utils::fgets_trunc(&buffer[m], MAXLINE, fp);
          if (eof == nullptr) error->one(FLERR, "Unexpected end of data file");
          ncount = utils::trim_and_count_words(&buffer[m]);
          if (ncount == 0) error->one(FLERR, "Too few values in body lines in data file");
          nword += ncount;
          m += strlen(&buffer[m]);
          onebody++;
        }
        if (nword > ninteger) error->one(FLERR, "Too many values in body lines in data file");

        nword = 0;
        while (nword < ndouble) {
          eof = utils::fgets_trunc(&buffer[m], MAXLINE, fp);
          if (eof == nullptr) error->one(FLERR, "Unexpected end of data file");
          ncount = utils::trim_and_count_words(&buffer[m]);
          if (ncount == 0) error->one(FLERR, "Too few values in body lines in data file");
          nword += ncount;
          m += strlen(&buffer[m]);
          onebody++;
        }
        if (nword > ndouble) error->one(FLERR, "Too many values in body lines in data file");

        if (onebody + 1 > MAXBODY)
          error->one(FLERR, "Too many lines in one body in data file - boost MAXBODY");

        nchunk++;
        nline += onebody + 1;
      }

      if (buffer[m - 1] != '\n') strcpy(&buffer[m++], "\n");
      m++;
    }

    MPI_Bcast(&nchunk, 1, MPI_INT, 0, world);
    MPI_Bcast(&m, 1, MPI_INT, 0, world);
    MPI_Bcast(buffer, m, MPI_CHAR, 0, world);

    if (firstpass) atom->data_bodies(nchunk, buffer, ptr, id_offset);
    nread += nchunk;
  }

  if (mapflag && firstpass) {
    atom->map_delete();
    atom->map_style = Atom::MAP_NONE;
  }

  if (me == 0 && firstpass) utils::logmesg(lmp, "  {} bodies\n", nblocks);
}

/* ---------------------------------------------------------------------- */

void ReadData::mass()
{
  settypeflag = 1;
  char *next;
  auto buf = new char[ntypes * MAXLINE];

  int eof = utils::read_lines_from_file(fp, ntypes, MAXLINE, buf, me, world);
  if (eof) error->all(FLERR, "Unexpected end of data file");
  if (tlabelflag && !lmap->is_complete(Atom::ATOM))
    error->all(FLERR, "Label map is incomplete: all types must be assigned a unique type label");

  char *original = buf;
  for (int i = 0; i < ntypes; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    atom->set_mass(FLERR, buf, toffset, tlabelflag, lmap->lmap2lmap.atom);
    buf = next + 1;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::paircoeffs()
{
  char *next;
  auto buf = new char[ntypes * MAXLINE];

  int eof = utils::read_lines_from_file(fp, ntypes, MAXLINE, buf, me, world);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  if (tlabelflag && !lmap->is_complete(Atom::ATOM))
    error->all(FLERR, "Label map is incomplete: all types must be assigned a unique type label");

  char *original = buf;
  for (int i = 0; i < ntypes; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    parse_coeffs(buf, nullptr, 1, 2, toffset, tlabelflag, lmap->lmap2lmap.atom);
    if (ncoeffarg == 0) error->all(FLERR, "Unexpected empty line in PairCoeffs section");
    force->pair->coeff(ncoeffarg, coeffarg);
    buf = next + 1;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::pairIJcoeffs()
{
  int i, j;
  char *next;

  int nsq = ntypes * (ntypes + 1) / 2;
  auto buf = new char[nsq * MAXLINE];

  int eof = utils::read_lines_from_file(fp, nsq, MAXLINE, buf, me, world);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  if (tlabelflag && !lmap->is_complete(Atom::ATOM))
    error->all(FLERR,
               "Label map is incomplete: "
               "all types must be assigned a unique type label");

  char *original = buf;
  for (i = 0; i < ntypes; i++)
    for (j = i; j < ntypes; j++) {
      next = strchr(buf, '\n');
      *next = '\0';
      parse_coeffs(buf, nullptr, 0, 2, toffset, tlabelflag, lmap->lmap2lmap.atom);
      if (ncoeffarg == 0) error->all(FLERR, "Unexpected empty line in PairCoeffs section");
      force->pair->coeff(ncoeffarg, coeffarg);
      buf = next + 1;
    }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::bondcoeffs()
{
  if (!nbondtypes) return;

  char *next;
  auto buf = new char[nbondtypes * MAXLINE];

  int eof = utils::read_lines_from_file(fp, nbondtypes, MAXLINE, buf, me, world);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  if (blabelflag && !lmap->is_complete(Atom::BOND))
    error->all(FLERR,
               "Label map is incomplete: "
               "all types must be assigned a unique type label");

  char *original = buf;
  for (int i = 0; i < nbondtypes; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    parse_coeffs(buf, nullptr, 0, 1, boffset, blabelflag, lmap->lmap2lmap.bond);
    if (ncoeffarg == 0) error->all(FLERR, "Unexpected empty line in BondCoeffs section");
    force->bond->coeff(ncoeffarg, coeffarg);
    buf = next + 1;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::anglecoeffs(int which)
{
  if (!nangletypes) return;

  char *next;
  auto buf = new char[nangletypes * MAXLINE];

  int eof = utils::read_lines_from_file(fp, nangletypes, MAXLINE, buf, me, world);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  if (alabelflag && !lmap->is_complete(Atom::ANGLE))
    error->all(FLERR,
               "Label map is incomplete: "
               "all types must be assigned a unique type label");

  char *original = buf;
  for (int i = 0; i < nangletypes; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    if (which == 0)
      parse_coeffs(buf, nullptr, 0, 1, aoffset, alabelflag, lmap->lmap2lmap.angle);
    else if (which == 1)
      parse_coeffs(buf, "bb", 0, 1, aoffset, alabelflag, lmap->lmap2lmap.angle);
    else if (which == 2)
      parse_coeffs(buf, "ba", 0, 1, aoffset, alabelflag, lmap->lmap2lmap.angle);
    else if (which == 3)
      parse_coeffs(buf, "ub", 0, 1, aoffset, alabelflag, lmap->lmap2lmap.angle);
    if (ncoeffarg == 0) error->all(FLERR, "Unexpected empty line in AngleCoeffs section");
    force->angle->coeff(ncoeffarg, coeffarg);
    buf = next + 1;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::dihedralcoeffs(int which)
{
  if (!ndihedraltypes) return;

  char *next;
  auto buf = new char[ndihedraltypes * MAXLINE];

  int eof = utils::read_lines_from_file(fp, ndihedraltypes, MAXLINE, buf, me, world);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  if (dlabelflag && !lmap->is_complete(Atom::DIHEDRAL))
    error->all(FLERR,
               "Label map is incomplete: "
               "all types must be assigned a unique type label");

  char *original = buf;
  for (int i = 0; i < ndihedraltypes; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    if (which == 0)
      parse_coeffs(buf, nullptr, 0, 1, doffset, dlabelflag, lmap->lmap2lmap.dihedral);
    else if (which == 1)
      parse_coeffs(buf, "mbt", 0, 1, doffset, dlabelflag, lmap->lmap2lmap.dihedral);
    else if (which == 2)
      parse_coeffs(buf, "ebt", 0, 1, doffset, dlabelflag, lmap->lmap2lmap.dihedral);
    else if (which == 3)
      parse_coeffs(buf, "at", 0, 1, doffset, dlabelflag, lmap->lmap2lmap.dihedral);
    else if (which == 4)
      parse_coeffs(buf, "aat", 0, 1, doffset, dlabelflag, lmap->lmap2lmap.dihedral);
    else if (which == 5)
      parse_coeffs(buf, "bb13", 0, 1, doffset, dlabelflag, lmap->lmap2lmap.dihedral);
    if (ncoeffarg == 0) error->all(FLERR, "Unexpected empty line in DihedralCoeffs section");
    force->dihedral->coeff(ncoeffarg, coeffarg);
    buf = next + 1;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::impropercoeffs(int which)
{
  if (!nimpropertypes) return;

  char *next;
  auto buf = new char[nimpropertypes * MAXLINE];

  int eof = utils::read_lines_from_file(fp, nimpropertypes, MAXLINE, buf, me, world);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  if (ilabelflag && !lmap->is_complete(Atom::IMPROPER))
    error->all(FLERR,
               "Label map is incomplete: "
               "all types must be assigned a unique type label");

  char *original = buf;
  for (int i = 0; i < nimpropertypes; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    if (which == 0)
      parse_coeffs(buf, nullptr, 0, 1, ioffset, ilabelflag, lmap->lmap2lmap.improper);
    else if (which == 1)
      parse_coeffs(buf, "aa", 0, 1, ioffset, ilabelflag, lmap->lmap2lmap.improper);
    if (ncoeffarg == 0) error->all(FLERR, "Unexpected empty line in ImproperCoeffs section");
    force->improper->coeff(ncoeffarg, coeffarg);
    buf = next + 1;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::typelabels(int mode)
{
  if (settypeflag) error->all(FLERR, "Must read Type Labels before any section involving types");

  int lntypes = 0;
  std::vector<std::string> *labels;
  std::unordered_map<std::string, int> *labels_map;
  if (me == 0)
    utils::logmesg(lmp, "  reading {} labelmap ...\n", utils::lowercase(labeltypes[mode]));

  switch (mode) {
    case Atom::ATOM:
      lntypes = lmap->natomtypes;
      labels = &lmap->typelabel;
      labels_map = &lmap->typelabel_map;
      break;
    case Atom::BOND:
      lntypes = lmap->nbondtypes;
      labels = &lmap->btypelabel;
      labels_map = &lmap->btypelabel_map;
      break;
    case Atom::ANGLE:
      lntypes = lmap->nangletypes;
      labels = &lmap->atypelabel;
      labels_map = &lmap->atypelabel_map;
      break;
    case Atom::DIHEDRAL:
      lntypes = lmap->ndihedraltypes;
      labels = &lmap->dtypelabel;
      labels_map = &lmap->dtypelabel_map;
      break;
    case Atom::IMPROPER:
      lntypes = lmap->nimpropertypes;
      labels = &lmap->itypelabel;
      labels_map = &lmap->itypelabel_map;
      break;
  }

  if (lntypes == 0) return;

  if (!atom->labelmapflag) atom->add_label_map();

  char *buf = new char[lntypes * MAXLINE];
  int eof = utils::read_lines_from_file(fp, lntypes, MAXLINE, buf, me, world);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  char *original = buf;
  char *next;
  for (int i = 0; i < lntypes; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    auto values = Tokenizer(buf).as_vector();
    int nwords = values.size();
    for (std::size_t ii = 0; ii < values.size(); ++ii) {
      if (utils::strmatch(values[ii], "^#")) {
        nwords = ii;
        break;
      }
    }
    if (nwords != 2)
      error->all(FLERR, "Invalid format in section: {} Type Labels: {}", labeltypes[mode], buf);
    values[1] = utils::utf8_subst(values[1]);
    if (utils::is_type(values[1]) != 1) error->all(FLERR, "Invalid type label {}", values[1]);
    int itype = utils::inumeric(FLERR, values[0], false, lmp);
    if ((itype < 1) || (itype > lntypes))
      error->all(FLERR, "Invalid type {} in section: {} Type Labels: {}", labeltypes[mode], buf);

    (*labels)[itype - 1] = values[1];
    (*labels_map)[values[1]] = itype;
    buf = next + 1;
  }
  delete[] original;

  for (int i = 0; i < lntypes; ++i) {
    if ((*labels)[i].empty())
      error->all(FLERR, "{} Type Labels map is incomplete", labeltypes[mode]);
  }
  if ((int) (*labels_map).size() != lntypes)
    error->all(FLERR, "{} Type Labels map is incomplete", labeltypes[mode]);

  // merge this read_data label map to atom class
  // determine mapping to let labels override numeric types
  // valid operations for first or subsequent data files

  atom->lmap->merge_lmap(lmap, mode);
  lmap->create_lmap2lmap(atom->lmap, mode);
}

/* ----------------------------------------------------------------------
   read fix section, pass lines to fix to process
   n = index of fix
------------------------------------------------------------------------- */

void ReadData::fix(Fix *ifix, char *keyword)
{
  int nchunk, eof;

  bigint nline = ifix->read_data_skip_lines(keyword);
  if (nline < 0) nline = natoms;

  bigint nread = 0;
  while (nread < nline) {
    nchunk = MIN(nline - nread, CHUNK);
    eof = utils::read_lines_from_file(fp, nchunk, MAXLINE, buffer, me, world);
    if (eof) error->all(FLERR, "Unexpected end of data file while reading section {}", keyword);
    ifix->read_data_section(keyword, nchunk, buffer, id_offset);
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
  memory->grow(count, amax + 1, "read_data:count");
  for (int i = cmax; i <= amax; i++) count[i] = 0;
  *pcount = count;
  return amax + 1;
}

/* ----------------------------------------------------------------------
   proc 0 opens data file
   test if compressed
------------------------------------------------------------------------- */

void ReadData::open(const std::string &file)
{
  if (platform::has_compress_extension(file)) {
    compressed = 1;
    fp = platform::compressed_read(file);
    if (!fp) error->one(FLERR, "Cannot open compressed file {}", file);
  } else {
    compressed = 0;
    fp = fopen(file.c_str(), "r");
    if (!fp) error->one(FLERR, "Cannot open file {}: {}", file, utils::getsyserror());
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
      if (utils::fgets_trunc(line, MAXLINE, fp) == nullptr) eof = 1;
    }
    while (eof == 0 && done == 0) {
      int blank = strspn(line, " \t\n\r");
      if ((blank == (int) strlen(line)) || (line[blank] == '#')) {
        if (utils::fgets_trunc(line, MAXLINE, fp) == nullptr) eof = 1;
      } else
        done = 1;
    }
    if (utils::fgets_trunc(buffer, MAXLINE, fp) == nullptr) {
      eof = 1;
      buffer[0] = '\0';
    }
  }

  // if eof, set keyword empty and return

  MPI_Bcast(&eof, 1, MPI_INT, 0, world);
  if (eof) {
    keyword[0] = '\0';
    return;
  }

  // bcast keyword line to all procs

  int n;
  if (me == 0) n = strlen(line) + 1;
  MPI_Bcast(&n, 1, MPI_INT, 0, world);
  MPI_Bcast(line, n, MPI_CHAR, 0, world);

  // store optional "style" following comment char '#' after keyword

  char *ptr;
  if ((ptr = strchr(line, '#'))) {
    *ptr++ = '\0';
    while (*ptr == ' ' || *ptr == '\t') ptr++;
    int stop = strlen(ptr) - 1;
    while (ptr[stop] == ' ' || ptr[stop] == '\t' || ptr[stop] == '\n' || ptr[stop] == '\r') stop--;
    ptr[stop + 1] = '\0';
    strcpy(style, ptr);
  } else
    style[0] = '\0';

  // copy non-whitespace portion of line into keyword

  int start = strspn(line, " \t\n\r");
  int stop = strlen(line) - 1;
  while (line[stop] == ' ' || line[stop] == '\t' || line[stop] == '\n' || line[stop] == '\r')
    stop--;
  line[stop + 1] = '\0';
  strcpy(keyword, &line[start]);
}

/* ----------------------------------------------------------------------
   proc 0 reads N lines from file
   could be skipping Natoms lines, so use bigints
------------------------------------------------------------------------- */

void ReadData::skip_lines(bigint n)
{
  if (me) return;
  if (n <= 0) return;
  char *eof = nullptr;
  for (bigint i = 0; i < n; i++) eof = utils::fgets_trunc(line, MAXLINE, fp);
  if (eof == nullptr) error->one(FLERR, "Unexpected end of data file");
}

/* ----------------------------------------------------------------------
   parse a line of coeffs into words, storing them in ncoeffarg,coeffarg
   trim anything from '#' onward
   word strings remain in line, are not copied
   if addstr != nullptr, add addstr as extra arg for class2 angle/dihedral/improper
     if 2nd word starts with letter, then is hybrid style, add addstr after it
     else add addstr before 2nd word
   if dupflag, duplicate 1st word, so pair_coeff "2" becomes "2 2"
   if noffset, add offset to first noffset args, which are atom/bond/etc types
   if labelflag, use ilabel to find the correct remapping of numeric type
------------------------------------------------------------------------- */

void ReadData::parse_coeffs(char *line, const char *addstr, int dupflag, int noffset, int offset,
                            int labelmode, int *ilabel)
{
  settypeflag = 1;
  char *ptr;
  if ((ptr = strchr(line, '#'))) *ptr = '\0';

  ncoeffarg = 0;
  char *word = line;
  char *end = line + strlen(line) + 1;

  while (word < end) {
    word += strspn(word, " \t\r\n\f");
    word[strcspn(word, " \t\r\n\f")] = '\0';
    if (strlen(word) == 0) break;
    if (ncoeffarg == maxcoeffarg) {
      maxcoeffarg += DELTA;
      coeffarg =
          (char **) memory->srealloc(coeffarg, maxcoeffarg * sizeof(char *), "read_data:coeffarg");
    }
    if (addstr && ncoeffarg == 1 && !islower(word[0])) coeffarg[ncoeffarg++] = (char *) addstr;
    coeffarg[ncoeffarg++] = word;
    if (addstr && ncoeffarg == 2 && islower(word[0])) coeffarg[ncoeffarg++] = (char *) addstr;
    if (dupflag && ncoeffarg == 1) coeffarg[ncoeffarg++] = word;
    word += strlen(word) + 1;
  }

  // to avoid segfaults on empty lines

  if (ncoeffarg == 0) return;

  if (noffset) {
    int value = utils::inumeric(FLERR, coeffarg[0], false, lmp);
    if (labelmode) value = ilabel[value - 1];
    argoffset1 = std::to_string(value + offset);
    coeffarg[0] = (char *)argoffset1.c_str();
    if (noffset == 2) {
      value = utils::inumeric(FLERR, coeffarg[1], false, lmp);
      if (labelmode) value = ilabel[value - 1];
      argoffset2 = std::to_string(value + offset);
      coeffarg[1] = (char *)argoffset2.c_str();
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
  int i, delta, len, len1, len2;

  if ((one == nullptr) || (two == nullptr)) return 1;

  len1 = strlen(one);
  len2 = strlen(two);

  for (i = 0; suffixes[i] != nullptr; i++) {
    len = strlen(suffixes[i]);
    if ((delta = len1 - len) > 0)
      if (strcmp(one + delta, suffixes[i]) == 0) len1 = delta;
    if ((delta = len2 - len) > 0)
      if (strcmp(two + delta, suffixes[i]) == 0) len2 = delta;
  }

  if ((len1 == 0) || (len1 == len2) || (strncmp(one, two, len1) == 0)) return 1;
  return 0;
}
