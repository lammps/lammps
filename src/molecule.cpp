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

#include "molecule.h"

#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_body.h"
#include "body.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "json.h"
#include "label_map.h"
#include "math_eigen.h"
#include "math_extra.h"
#include "math_special.h"
#include "memory.h"
#include "tokenizer.h"
#include "update.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>

using namespace LAMMPS_NS;
using MathSpecial::powint;

static constexpr int MAXLINE = 1024;
static constexpr double EPSILON = 1.0e-7;
static constexpr double BIG = 1.0e20;

static constexpr double SINERTIA = 0.4;    // moment of inertia prefactor for sphere

/* ---------------------------------------------------------------------- */

Molecule::Molecule(LAMMPS *lmp) :
    Pointers(lmp), id(nullptr), x(nullptr), type(nullptr), molecule(nullptr), q(nullptr),
    radius(nullptr), rmass(nullptr), mu(nullptr), num_bond(nullptr), bond_type(nullptr),
    bond_atom(nullptr), num_angle(nullptr), angle_type(nullptr), angle_atom1(nullptr),
    angle_atom2(nullptr), angle_atom3(nullptr), num_dihedral(nullptr), dihedral_type(nullptr),
    dihedral_atom1(nullptr), dihedral_atom2(nullptr), dihedral_atom3(nullptr),
    dihedral_atom4(nullptr), num_improper(nullptr), improper_type(nullptr), improper_atom1(nullptr),
    improper_atom2(nullptr), improper_atom3(nullptr), improper_atom4(nullptr), nspecial(nullptr),
    special(nullptr), shake_flag(nullptr), shake_atom(nullptr), shake_type(nullptr),
    avec_body(nullptr), ibodyparams(nullptr), dbodyparams(nullptr), fragmentmask(nullptr),
    dx(nullptr), dxcom(nullptr), dxbody(nullptr), quat_external(nullptr), fp(nullptr),
    count(nullptr)
{
  // parse args until reach unknown arg (next file)

  toffset = 0;
  boffset = aoffset = doffset = ioffset = 0;
  sizescale = 1.0;
  json_format = 0;

  // initialize all fields to empty

  Molecule::initialize();
}

// ------------------------------------------------------------------------------
//   process arguments from "molecule" command
// ------------------------------------------------------------------------------

void Molecule::command(int narg, char **arg, int &index)
{
  if (index >= narg) utils::missing_cmd_args(FLERR, "molecule", error);

  id = utils::strdup(arg[0]);
  if (!utils::is_id(id))
    error->all(FLERR, Error::ARGZERO,
               "Molecule template ID {} must have only alphanumeric or underscore characters", id);

  // parse args until reach unknown arg (next file)

  fileiarg = index;

  int iarg = fileiarg + 1;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "offset") == 0) {
      if (iarg + 6 > narg) utils::missing_cmd_args(FLERR, "molecule offset", error);
      toffset = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      boffset = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
      aoffset = utils::inumeric(FLERR, arg[iarg + 3], false, lmp);
      doffset = utils::inumeric(FLERR, arg[iarg + 4], false, lmp);
      ioffset = utils::inumeric(FLERR, arg[iarg + 5], false, lmp);
      if (toffset < 0) error->all(FLERR, iarg + 1, "Illegal atom type offset {}", toffset);
      if (boffset < 0) error->all(FLERR, iarg + 2, "Illegal bond type offset {}", boffset);
      if (aoffset < 0) error->all(FLERR, iarg + 3, "Illegal angle type offset {}", aoffset);
      if (doffset < 0) error->all(FLERR, iarg + 4, "Illegal dihedral type offset {}", doffset);
      if (ioffset < 0) error->all(FLERR, iarg + 5, "Illegal improper type offset {}", ioffset);
      iarg += 6;
    } else if (strcmp(arg[iarg], "toff") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "molecule toff", error);
      toffset = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (toffset < 0) error->all(FLERR, iarg + 1, "Illegal atom type offset {}", toffset);
      iarg += 2;
    } else if (strcmp(arg[iarg], "boff") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "molecule boff", error);
      boffset = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (boffset < 0) error->all(FLERR, iarg + 1, "Illegal bond type offset {}", boffset);
      iarg += 2;
    } else if (strcmp(arg[iarg], "aoff") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "molecule aoff", error);
      aoffset = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (aoffset < 0) error->all(FLERR, iarg + 1, "Illegal angle type offset {}", aoffset);
      iarg += 2;
    } else if (strcmp(arg[iarg], "doff") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "molecule doff", error);
      doffset = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (doffset < 0) error->all(FLERR, iarg + 1, "Illegal dihedral type offset {}", doffset);
      iarg += 2;
    } else if (strcmp(arg[iarg], "ioff") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "molecule ioff", error);
      ioffset = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (ioffset < 0) error->all(FLERR, iarg + 1, "Illegal improper type offset {}", ioffset);
      iarg += 2;
    } else if (strcmp(arg[iarg], "scale") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "molecule scale", error);
      sizescale = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (sizescale <= 0.0) error->all(FLERR, iarg + 1, "Illegal scale factor {}", sizescale);
      iarg += 2;
    } else
      break;
  }
  // clang-format on
  index = iarg;

  if (atom->labelmapflag &&
      ((toffset > 0) || (boffset > 0) || (aoffset > 0) || (doffset > 0) || (ioffset > 0))) {
    if (comm->me == 0)
      error->warning(FLERR,
                     "Using molecule command with type offsets and a labelmap. "
                     "Offsets will be only applied to numeric types and not to type labels");
  }

  // last molecule if have scanned all args

  if (iarg == narg)
    last = 1;
  else
    last = 0;

  // JSON files must have the extension .json

  std::string filename = arg[fileiarg];
  if (utils::strmatch(filename, "\\.json$")) {

    json moldata;
    std::vector<std::uint8_t> jsondata;
    int jsondata_size = 0;

    if (comm->me == 0) {
      fp = fopen(filename.c_str(), "r");
      if (fp == nullptr)
        error->one(FLERR, fileiarg, "Cannot open molecule file {}: {}", filename,
                   utils::getsyserror());
      try {
        // try to parse as a JSON file. parser throws an exception on errors
        // if successful, temporarily serialize to bytearray for communication
        moldata = json::parse(fp);
        jsondata = json::to_ubjson(moldata);
        jsondata_size = jsondata.size();
        fclose(fp);
      } catch (std::exception &e) {
        fclose(fp);
        error->one(FLERR, fileiarg, "Error parsing JSON file {}: {}", filename, e.what());
      }
    }
    MPI_Bcast(&jsondata_size, 1, MPI_INT, 0, world);

    if (jsondata_size > 0) {

      // broadcast binary JSON data to all processes and deserialize again

      if (comm->me != 0) jsondata.resize(jsondata_size);
      MPI_Bcast(jsondata.data(), jsondata_size, MPI_CHAR, 0, world);

      // convert back to json class on all processors and free temporary storage
      moldata.clear();
      moldata = json::from_ubjson(jsondata);
      jsondata.clear();    // free binary data

      // process JSON data
      Molecule::from_json(id, moldata);
    } else {
      error->all(FLERR, "Molecule file {} does not contain JSON data", filename);
    }

  } else {    // process native molecule file

    if (comm->me == 0) {
      fp = fopen(filename.c_str(), "r");
      if (fp == nullptr)
        error->one(FLERR, fileiarg, "Cannot open molecule file {}: {}", filename,
                   utils::getsyserror());
    }

    // scan file for sizes of all fields and allocate storage for them

    Molecule::read(0);
    Molecule::allocate();

    // read file again to populate all fields

    if (comm->me == 0) rewind(fp);
    Molecule::read(1);
    if (comm->me == 0) fclose(fp);
  }
  Molecule::stats();
}

// ------------------------------------------------------------------------------
//  convert json data structure to molecule data structure
// ------------------------------------------------------------------------------

void Molecule::from_json(const std::string &molid, const json &moldata)
{
  json_format = 1;
  if (!utils::is_id(molid))
    error->all(FLERR, Error::NOLASTLINE,
               "Molecule template ID {} must have only alphanumeric or underscore characters",
               molid);
  delete[] id;
  id = utils::strdup(molid);

  // check required fields if JSON data is compatible

  std::string val;
  if (moldata.contains("application")) {
    if (moldata["application"] != "LAMMPS")
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: JSON data is for incompatible application: {}", id,
                 std::string(moldata["application"]));
  } else {
    error->all(FLERR, Error::NOLASTLINE,
               "Molecule template {}: JSON data does not contain required \"application\" field",
               id);
  }
  if (moldata.contains("format")) {
    if (moldata["format"] != "molecule")
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: JSON data is not for a molecule: {}", id,
                 std::string(moldata["format"]));
  } else {
    error->all(FLERR, Error::NOLASTLINE,
               "Molecule template {}: JSON data does not contain required \"format\" field", id);
  }
  if (moldata.contains("revision")) {
    int rev = moldata["revision"];
    if ((rev < 1) || (rev > 1))
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: JSON molecule data with unsupported revision {}", id, rev);
  } else {
    error->all(FLERR, Error::NOLASTLINE,
               "Molecule template {}: JSON data does not contain required \"revision\" field", id);
  }

  // length of types data list determines the number of atoms in the template and is thus required
  if (!moldata.contains("types"))
    error->all(FLERR, Error::NOLASTLINE,
               "Molecule template {}: JSON data does not contain required \"types\" field", id);

  // optional fields

  // check for compatible units

  if (moldata.contains("units")) {
    bool incompatible_units = true;
    auto jsonunits = std::string(moldata["units"]);
    auto lammpsunits = std::string(update->unit_style);
    if ((jsonunits == "real") || (jsonunits == "metal")) {
      if ((lammpsunits == "real") || (lammpsunits == "metal")) incompatible_units = false;
    } else if (jsonunits == lammpsunits) {
      incompatible_units = false;
    }

    if (incompatible_units)
      error->all(
          FLERR, Error::NOLASTLINE,
          "Molecule template {}: Incompatible units in JSON molecule data: current = {}, JSON = {}",
          id, lammpsunits, jsonunits);
  }
  if (moldata.contains("title")) title = moldata["title"];

  // determine and check sizes

  int dummyvar = 0;

#define JSON_INIT_FIELD(field, sizevar, flagvar, required, sizecheck)                              \
  sizevar = 0;                                                                                     \
  flagvar = 0;                                                                                     \
  if (moldata.contains(#field)) {                                                                  \
    if (!moldata[#field].contains("format"))                                                       \
      error->all(FLERR, Error::NOLASTLINE,                                                         \
                 "Molecule template {}: JSON molecule data does not contain required \"format\" "  \
                 "field for '{}'",                                                                 \
                 id, #field);                                                                      \
    if (moldata[#field].contains("data")) {                                                        \
      flagvar = 1;                                                                                 \
      sizevar = moldata[#field]["data"].size();                                                    \
    } else {                                                                                       \
      error->all(FLERR, Error::NOLASTLINE,                                                         \
                 "Molecule template {}: JSON molecule data does not contain required \"data\" "    \
                 "field for '{}'",                                                                 \
                 id, #field);                                                                      \
    }                                                                                              \
    if (sizevar < 1)                                                                               \
      error->all(FLERR, Error::NOLASTLINE,                                                         \
                 "Molecule template {}: No {} entries in JSON data for molecule", id, #field);     \
  } else {                                                                                         \
    if (required)                                                                                  \
      error->all(                                                                                  \
          FLERR, Error::NOLASTLINE,                                                                \
          "Molecule template {}: JSON data for molecule does not contain required '{}' field", id, \
          #field);                                                                                 \
  }                                                                                                \
  if (flagvar && sizecheck && (sizecheck != sizevar))                                              \
    error->all(FLERR, Error::NOLASTLINE,                                                           \
               "Molecule template {}: Found {} instead of {} data entries for '{}'", id, sizevar,  \
               sizecheck, #field);

  JSON_INIT_FIELD(types, natoms, typeflag, true, 0);
  JSON_INIT_FIELD(coords, dummyvar, xflag, false, natoms);
  JSON_INIT_FIELD(molecules, dummyvar, moleculeflag, false, natoms);
  JSON_INIT_FIELD(fragments, nfragments, fragmentflag, false, 0);
  JSON_INIT_FIELD(charges, dummyvar, qflag, false, natoms);
  JSON_INIT_FIELD(diameters, dummyvar, radiusflag, false, natoms);
  JSON_INIT_FIELD(dipoles, dummyvar, muflag, false, natoms);
  JSON_INIT_FIELD(masses, dummyvar, rmassflag, false, natoms);
  JSON_INIT_FIELD(bonds, nbonds, bondflag, false, 0);
  JSON_INIT_FIELD(angles, nangles, angleflag, false, 0);
  JSON_INIT_FIELD(dihedrals, ndihedrals, dihedralflag, false, 0);
  JSON_INIT_FIELD(impropers, nimpropers, improperflag, false, 0);

#undef JSON_INIT_FIELD
  // special is nested

  if (moldata.contains("special")) {
    if (moldata["special"].contains("counts")) {
      nspecialflag = specialflag_user = 1;
      maxspecial = 0;
      const auto &specialcounts = moldata["special"]["counts"];
      if (!specialcounts.contains("format"))
        error->all(FLERR, Error::NOLASTLINE,
                   "Molecule template {}: JSON molecule data does not contain required \"format\" "
                   "field for 'special:counts'",
                   id);
      if (specialcounts.contains("data")) {
        if ((int) specialcounts["data"].size() != natoms)
          error->all(
              FLERR, Error::NOLASTLINE,
              "Molecule template {}: Found {} instead of {} data entries for 'special:counts'", id,
              specialcounts["data"].size(), natoms);
        for (const auto &item : specialcounts["data"]) {
          if (item.size() != 4)
            error->all(
                FLERR, Error::NOLASTLINE,
                "Molecule template {}: Found {} instead of 4 data entries for 'special:counts'", id,
                item.size());

          const auto &vals = item.get<std::vector<int>>();
          int sumspecial = vals[1] + vals[2] + vals[3];
          maxspecial = MAX(maxspecial, sumspecial);
        }
      } else {
        error->all(FLERR, Error::NOLASTLINE,
                   "Molecule template {}: JSON molecule data does not contain required \"data\" "
                   "field for 'special:counts'",
                   id);
      }
    } else {
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: JSON molecule data does not contain required 'counts' "
                 "field for \"special\"",
                 id);
    }

    if (moldata["special"].contains("bonds")) {
      specialflag = tag_require = 1;
      const auto &specialbonds = moldata["special"]["bonds"];
      if (!specialbonds.contains("format"))
        error->all(FLERR, Error::NOLASTLINE,
                   "Molecule template {}: JSON molecule data does not contain required \"format\" "
                   "field for \"special:bonds\"",
                   id);
      if (specialbonds.contains("data")) {
        if ((int) specialbonds["data"].size() != natoms)
          error->all(
              FLERR, Error::NOLASTLINE,
              "Molecule template {}: Found {} instead of {} data entries for \"special:bonds\"", id,
              specialbonds["data"].size(), natoms);
        if (specialbonds["data"][0].size() != 2)
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: \"special:bonds\" is incorrectly formatted: {}", id,
                     to_string(specialbonds["data"][0]));
        for (int i = 0; i < natoms; ++i) {
          if ((int) specialbonds["data"][i][1].size() > maxspecial)
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: Number of data entries in \"special:bonds\" for atom "
                       "{} exceeds limit: {} vs {}",
                       id, specialbonds["data"][i][1].size(), maxspecial);
        }
      } else {
        error->all(FLERR, Error::NOLASTLINE,
                   "Molecule template {}: JSON molecule data does not contain required \"data\" "
                   "field for \"special:bonds\"",
                   id);
      }
    } else {
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: JSON molecule data does not contain required \"bonds\" "
                 "field for \"special\"",
                 id);
    }
  }

  // shake is nested

  if (moldata.contains("shake")) {
    shakeflag = shakeflagflag = shaketypeflag = shakeatomflag = 0;
    const auto &shakedata = moldata["shake"];

    if (shakedata.contains("flags")) {
      if (!shakedata["flags"].contains("format"))
        error->all(FLERR, Error::NOLASTLINE,
                   "Molecule template {}: JSON molecule data does not contain required \"format\" "
                   "field for \"shake:flags\"",
                   id);
      if (shakedata["flags"].contains("data")) {
        shakeflagflag = 1;
        if ((int) shakedata["flags"]["data"].size() != natoms)
          error->all(
              FLERR, Error::NOLASTLINE,
              "Molecule template {}: Found {} instead of {} data entries for \"shake:flags\"", id,
              shakedata["flags"]["data"].size(), natoms);
      } else {
        error->all(FLERR, Error::NOLASTLINE,
                   "Molecule template {}: JSON molecule data does not contain required \"data\" "
                   "field for \"shake:flags\"",
                   id);
      }
    } else {
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: JSON molecule data does not contain required 'flags' "
                 "field for \"shake\"",
                 id);
    }

    if (shakedata.contains("atoms")) {
      if (!shakedata["atoms"].contains("format"))
        error->all(FLERR, Error::NOLASTLINE,
                   "Molecule template {}: JSON molecule data does not contain required \"format\" "
                   "field for \"shake:atoms\"",
                   id);
      if (shakedata["atoms"].contains("data")) {
        shakeatomflag = 1;
        tag_require = 1;
        if ((int) shakedata["atoms"]["data"].size() != natoms)
          error->all(
              FLERR, Error::NOLASTLINE,
              "Molecule template {}: Found {} instead of {} data entries for \"shake:atoms\"", id,
              shakedata["atoms"]["data"].size(), natoms);
      } else {
        error->all(FLERR, Error::NOLASTLINE,
                   "Molecule template {}: JSON molecule data does not contain required \"data\" "
                   "field for \"shake:atoms\"",
                   id);
      }
    } else {
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: JSON molecule data does not contain required \"atoms\" "
                 "field for \"shake\"",
                 id);
    }

    if (shakedata.contains("types")) {
      if (!shakedata["types"].contains("format"))
        error->all(FLERR, Error::NOLASTLINE,
                   "Molecule template {}: JSON molecule data does not contain required \"format\" "
                   "field for \"shake:types\"",
                   id);
      if (shakedata["types"].contains("data")) {
        shaketypeflag = 1;
        tag_require = 1;
        if ((int) shakedata["types"]["data"].size() != natoms)
          error->all(
              FLERR, Error::NOLASTLINE,
              "Molecule template {}: Found {} instead of {} data entries for \"shake:types\"", id,
              shakedata["types"]["data"].size(), natoms);
      } else {
        error->all(FLERR, Error::NOLASTLINE,
                   "Molecule template {}: JSON molecule data does not contain required \"data\" "
                   "field for \"shake:types\"",
                   id);
      }
    } else {
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: JSON molecule data does not contain required \"types\" "
                 "field for \"shake\"",
                 id);
    }
    if (shakeflagflag && shakeatomflag && shaketypeflag) shakeflag = 1;
  }

  if ((nbonds > 0) || (nangles > 0) || (ndihedrals > 0) || (nimpropers > 0)) tag_require = 1;

  // extract global properties, if present

  if (moldata.contains("masstotal")) {
    massflag = massflag_user = 1;
    masstotal = double(moldata["masstotal"]) * sizescale * sizescale * sizescale;
  }

  if (moldata.contains("com") && (moldata["com"].size() == 3)) {
    comflag = comflag_user = 1;
    com[0] = double(moldata["com"][0]) * sizescale;
    com[1] = double(moldata["com"][1]) * sizescale;
    com[2] = double(moldata["com"][2]) * sizescale;
  }

  if (moldata.contains("inertia") && (moldata["inertia"].size() == 6)) {
    inertiaflag = inertiaflag_user = 1;
    const double scale5 = powint(sizescale, 5);
    itensor[0] = double(moldata["inertia"][0]) * scale5;
    itensor[1] = double(moldata["inertia"][1]) * scale5;
    itensor[2] = double(moldata["inertia"][2]) * scale5;
    itensor[3] = double(moldata["inertia"][3]) * scale5;
    itensor[4] = double(moldata["inertia"][4]) * scale5;
    itensor[5] = double(moldata["inertia"][5]) * scale5;
  }

  if (moldata.contains("body")) {
    avec_body = dynamic_cast<AtomVecBody *>(atom->style_match("body"));
    if (!avec_body)
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: JSON molecule data requires atom style body", id);

    if (moldata["body"].contains("integers") && moldata["body"].contains("doubles")) {
      bodyflag = radiusflag = dbodyflag = ibodyflag = 1;
      nibody = moldata["body"]["integers"].size();
      ndbody = moldata["body"]["doubles"].size();
    } else {
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: JSON molecule \"body\" data requires \"integers\" and "
                 "\"doubles\" sections",
                 id);
    }
  }

  // checks. No checks for < 0 needed since size() is at least 0

  if ((domain->dimension == 2) && (com[2] != 0.0))
    error->all(FLERR, Error::NOLASTLINE,
               "Molecule template {}: Molecule data z center-of-mass must be 0.0 for 2d systems",
               id);

  // allocate required storage

  Molecule::allocate();

  // count = vector for tallying bonds,angles,etc per atom

  memory->create(count, natoms, "molecule:count");

  // process data sections
  std::vector<std::string> secfmt;

  // coords
  if (xflag) {
    for (int i = 0; i < 4; ++i) secfmt.emplace_back(moldata["coords"]["format"][i]);
    if ((secfmt[0] == "atom-id") && (secfmt[1] == "x") && (secfmt[2] == "y") &&
        (secfmt[3] == "z")) {

      memset(count, 0, natoms * sizeof(int));
      for (const auto &c : moldata["coords"]["data"]) {
        if (c.size() < 4)
          error->all(
              FLERR, Error::NOLASTLINE,
              "Molecule template {}: missing data in \"coords\" section of molecule JSON data: {}",
              id, to_string(c));
        if (!c[0].is_number_integer())
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid atom-id in \"coords\" section of molecule JSON "
                     "data: {}",
                     id, to_string(c[0]));

        const int iatom = int(c[0]) - 1;
        if ((iatom < 0) || (iatom >= natoms))
          error->all(
              FLERR, Error::NOLASTLINE,
              "Molecule template {}: invalid atom-id {} in coords section of molecule JSON data",
              id, iatom + 1);
        count[iatom]++;
        x[iatom][0] = c[1];
        x[iatom][1] = c[2];
        x[iatom][2] = c[3];

        x[iatom][0] *= sizescale;
        x[iatom][1] *= sizescale;
        x[iatom][2] *= sizescale;
      }

      // checks
      for (int i = 0; i < natoms; i++) {
        if (count[i] == 0) {
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: atom {} missing in \"coords\" JSON section", id, i + 1);
        }
      }
      if (domain->dimension == 2) {
        for (int i = 0; i < natoms; i++) {
          if (x[i][2] != 0.0) {
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: Z coord for atom {} must be 0.0 for 2d-simulation",
                       id, i + 1);
          }
        }
      }
    } else {
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: Expected \"coords\" format [\"atom-id\",\"x\",\"y\",\"z\"] "
                 "but found [\"{}\",\"{}\",\"{}\",\"{}\"]",
                 id, secfmt[0], secfmt[1], secfmt[2], secfmt[3]);
    }
  }

  // types (is a required section and we tested for it above)

  secfmt.clear();
  for (int i = 0; i < 2; ++i) secfmt.emplace_back(moldata["types"]["format"][i]);
  if ((secfmt[0] == "atom-id") && (secfmt[1] == "type")) {

    memset(count, 0, natoms * sizeof(int));
    for (const auto &c : moldata["types"]["data"]) {
      if (c.size() < 2)
        error->all(
            FLERR, Error::NOLASTLINE,
            "Molecule template {}: missing data in \"types\" section of molecule JSON data: {}", id,
            to_string(c));
      if (!c[0].is_number_integer())
        error->all(
            FLERR, Error::NOLASTLINE,
            "Molecule template {}: invalid atom-id in \"types\" section of molecule JSON data: {}",
            id, to_string(c[0]));
      const int iatom = int(c[0]) - 1;
      if ((iatom < 0) || (iatom >= natoms))
        error->all(
            FLERR, Error::NOLASTLINE,
            "Molecule template {}: invalid atom-id {} in types section of molecule JSON data", id,
            iatom + 1);
      if (c[1].is_number_integer()) {    // numeric type
        type[iatom] = int(c[1]) + toffset;
      } else {
        const auto &typestr = std::string(c[1]);
        if (!atom->labelmapflag)
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid atom type in \"types\" JSON section", id,
                     typestr);
        type[iatom] = atom->lmap->find(typestr, Atom::ATOM);
        if (type[iatom] == -1)
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: Unknown atom type {} in \"types\" JSON section", id,
                     typestr);
      }
      count[iatom]++;
    }
    // checks
    ntypes = 0;
    for (int i = 0; i < natoms; i++) {
      if (count[i] == 0) {
        error->all(FLERR, Error::NOLASTLINE,
                   "Molecule template {}: atom {} missing in \"types\" JSON section", id, i + 1);
      }
      for (int i = 0; i < natoms; i++) {
        if ((type[i] <= 0) || (domain->box_exist && (type[i] > atom->ntypes)))
          error->all(FLERR, fileiarg, "Invalid atom type {} for atom {} in molecule file", type[i],
                     i + 1);
        ntypes = MAX(ntypes, type[i]);
      }
    }
  } else {
    error->all(FLERR, Error::NOLASTLINE,
               "Molecule template {}: Expected \"types\" format [\"atom-id\",\"type\"] but found "
               "[\"{}\",\"{}\"]",
               id, secfmt[0], secfmt[1]);
  }

  // molecules

  if (moleculeflag) {

    secfmt.clear();
    for (int i = 0; i < 2; ++i) secfmt.emplace_back(moldata["molecules"]["format"][i]);
    if ((secfmt[0] == "atom-id") && (secfmt[1] == "molecule-id")) {

      memset(count, 0, natoms * sizeof(int));
      for (const auto &c : moldata["molecules"]["data"]) {
        if (c.size() < 2)
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: missing data in \"molecules\" section of molecule JSON "
                     "data: {}",
                     id, to_string(c));
        if (!c[0].is_number_integer())
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid atom-id in \"molecules\" section of molecule "
                     "JSON data: {}",
                     id, to_string(c[0]));

        const int iatom = int(c[0]) - 1;
        if ((iatom < 0) || (iatom >= natoms))
          error->all(FLERR, Error::NOLASTLINE,
                     "Invalid atom-id {} in \"molecules\" section of molecule JSON data",
                     iatom + 1);
        if (!c[1].is_number_integer())
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid molecule-id in \"molecules\" section of "
                     "molecule JSON data: {}",
                     id, to_string(c[1]));
        molecule[iatom] = int(c[1]);
        if (molecule[iatom] < 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid molecule-id in \"molecules\" section of "
                     "molecule JSON data: {}",
                     id, to_string(c[1]));
        count[iatom]++;
      }
      // checks
      for (int i = 0; i < natoms; i++) {
        if (count[i] == 0) {
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: atom {} missing in \"molecules\" JSON section", id,
                     i + 1);
        }
      }
    } else {
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: Expected \"molecules\" format "
                 "[\"atom-id\",\"molecule-id\"] but found "
                 "[\"{}\",\"{}\"]",
                 id, secfmt[0], secfmt[1]);
    }
  }

  // fragments

  if (fragmentflag) {
    secfmt.clear();
    for (int i = 0; i < 2; ++i) secfmt.emplace_back(moldata["fragments"]["format"][i]);
    if ((secfmt[0] == "fragment-id") && (secfmt[1] == "atom-id-list")) {

      for (int i = 0; i < nfragments; ++i) {
        fragmentnames[i] = to_string(moldata["fragments"]["data"][i][0]);
        for (const auto &c : moldata["fragments"]["data"][i][1]) {
          if (!c.is_number_integer())
            error->all(
                FLERR, Error::NOLASTLINE,
                "Molecule template {}: invalid atom-id in \"fragments\" section  JSON data: {}", id,
                to_string(c));

          const int iatom = int(c) - 1;
          if ((iatom < 0) || (iatom >= natoms))
            error->all(FLERR, Error::NOLASTLINE,
                       "Invalid atom id {} in \"fragments\" section of molecule JSON data",
                       iatom + 1);
          fragmentmask[i][iatom] = 1;
        }
      }
    } else {
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: Expected \"fragments\" format "
                 "[\"fragment-id\",\"atom-id-list\"] but found [\"{}\",\"{}\"]",
                 id, secfmt[0], secfmt[1]);
    }
  }

  // charges

  if (qflag) {

    secfmt.clear();
    for (int i = 0; i < 2; ++i) secfmt.emplace_back(moldata["charges"]["format"][i]);
    if ((secfmt[0] == "atom-id") && (secfmt[1] == "charge")) {

      memset(count, 0, natoms * sizeof(int));
      for (const auto &c : moldata["charges"]["data"]) {
        if (c.size() < 2)
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: missing data in \"charges\" section of molecule JSON "
                     "data: {}",
                     id, to_string(c));
        if (!c[0].is_number_integer())
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid atom-id in \"charges\" section of molecule "
                     "JSON data: {}",
                     id, to_string(c[0]));

        const int iatom = int(c[0]) - 1;
        if ((iatom < 0) || (iatom >= natoms))
          error->all(FLERR, Error::NOLASTLINE,
                     "Invalid atom-id {} in \"charges\" section of molecule JSON data", iatom + 1);
        if (!c[1].is_number())
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid charge in \"charges\" section of "
                     "molecule JSON data: {}",
                     id, to_string(c[1]));
        q[iatom] = double(c[1]);
        count[iatom]++;
      }
      // checks
      for (int i = 0; i < natoms; i++) {
        if (count[i] == 0) {
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: atom {} missing in \"charges\" JSON section", id,
                     i + 1);
        }
      }
    } else {
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: Expected \"charges\" format [\"atom-id\",\"charge\"] but "
                 "found [\"{}\",\"{}\"]",
                 id, secfmt[0], secfmt[1]);
    }
  }

  // diameters

  if (radiusflag && !bodyflag) {
    maxradius = 0.0;
    secfmt.clear();
    for (int i = 0; i < 2; ++i) secfmt.emplace_back(moldata["diameters"]["format"][i]);
    if ((secfmt[0] == "atom-id") && (secfmt[1] == "diameter")) {

      memset(count, 0, natoms * sizeof(int));
      for (const auto &c : moldata["diameters"]["data"]) {
        if (c.size() < 2)
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: missing data in \"diameters\" section of molecule JSON "
                     "data: {}",
                     id, to_string(c));
        if (!c[0].is_number_integer())
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid atom-id in \"diameters\" section of molecule "
                     "JSON data: {}",
                     id, to_string(c[0]));

        const int iatom = int(c[0]) - 1;
        if ((iatom < 0) || (iatom >= natoms))
          error->all(FLERR, Error::NOLASTLINE,
                     "Invalid atom-id {} in \"diameters\" section of molecule JSON data",
                     iatom + 1);
        if (!c[1].is_number())
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid diameter in \"diameters\" section of "
                     "molecule JSON data: {}",
                     id, to_string(c[1]));
        radius[iatom] = double(c[1]) * sizescale * 0.5;
        maxradius = MAX(maxradius, radius[iatom]);
        if (!c[1].is_number())
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid diameter in \"diameters\" section of "
                     "molecule JSON data: {}",
                     id, to_string(c[1]));
        count[iatom]++;
      }
      // checks
      for (int i = 0; i < natoms; i++) {
        if (count[i] == 0) {
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: atom {} missing in \"diameters\" JSON section", id,
                     i + 1);
        }
        if (radius[i] < 0.0)
          error->all(
              FLERR, Error::NOLASTLINE,
              "Molecule template {}: invalid atom diameter {} for atom {} in molecule JSON data",
              id, radius[i] * 2.0 / sizescale, i + 1);
      }
    } else {
      error->all(
          FLERR, Error::NOLASTLINE,
          "Molecule template {}: Expected \"diameters\" format [\"atom-id\",\"diameter\"] but "
          "found [\"{}\",\"{}\"]",
          id, secfmt[0], secfmt[1]);
    }
  }

  // dipoles

  if (muflag) {

    secfmt.clear();
    for (int i = 0; i < 4; ++i) secfmt.emplace_back(moldata["dipoles"]["format"][i]);
    if ((secfmt[0] == "atom-id") && (secfmt[1] == "mux") && (secfmt[2] == "muy") &&
        (secfmt[3] == "muz")) {

      memset(count, 0, natoms * sizeof(int));
      for (const auto &c : moldata["dipoles"]["data"]) {
        if (c.size() < 4)
          error->all(
              FLERR, Error::NOLASTLINE,
              "Molecule template {}: missing data in \"dipoles\" section of molecule JSON data: {}",
              id, to_string(c));
        if (!c[0].is_number_integer())
          error->all(
              FLERR, Error::NOLASTLINE,
              "Molecule template {}: invalid atom-id in \"dipoles\" section of molecule JSON "
              "data: {}",
              id, to_string(c[0]));

        const int iatom = int(c[0]) - 1;
        if ((iatom < 0) || (iatom >= natoms))
          error->all(
              FLERR, Error::NOLASTLINE,
              "Molecule template {}: invalid atom-id {} in dipoles section of molecule JSON data",
              id, iatom + 1);
        count[iatom]++;
        mu[iatom][0] = c[1];
        mu[iatom][1] = c[2];
        mu[iatom][2] = c[3];
        mu[iatom][0] *= sizescale;
        mu[iatom][1] *= sizescale;
        mu[iatom][2] *= sizescale;
      }

      // checks
      for (int i = 0; i < natoms; i++) {
        if (count[i] == 0) {
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: atom {} missing in \"dipoles\" JSON section", id,
                     i + 1);
        }
      }
      if (domain->dimension == 2) {
        for (int i = 0; i < natoms; i++)
          if (mu[i][2] != 0.0)
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: dipole moment z-component in JSON data for atom {} "
                       "must be 0.0 for 2d-simulation",
                       id, i + 1);
      }
    } else {
      error->all(
          FLERR, Error::NOLASTLINE,
          "Molecule template {}: Expected \"dipoles\" format [\"atom-id\",\"mux\",\"muy\",\"muz\"] "
          "but found [\"{}\",\"{}\",\"{}\",\"{}\"]",
          id, secfmt[0], secfmt[1], secfmt[2], secfmt[3]);
    }
  }

  // masses

  if (rmassflag) {
    secfmt.clear();
    for (int i = 0; i < 2; ++i) secfmt.emplace_back(moldata["masses"]["format"][i]);
    if ((secfmt[0] == "atom-id") && (secfmt[1] == "mass")) {

      memset(count, 0, natoms * sizeof(int));
      for (const auto &c : moldata["masses"]["data"]) {
        if (c.size() < 2)
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: missing data in \"masses\" section of molecule JSON "
                     "data: {}",
                     id, to_string(c));
        if (!c[0].is_number_integer())
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid atom-id in \"masses\" section of molecule "
                     "JSON data: {}",
                     id, to_string(c[0]));

        const int iatom = int(c[0]) - 1;
        if ((iatom < 0) || (iatom >= natoms))
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid atom-id {} in \"masses\" section of molecule "
                     "JSON data",
                     iatom + 1);
        if (!c[1].is_number())
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid mass in \"masses\" section of "
                     "molecule JSON data: {}",
                     id, to_string(c[1]));
        rmass[iatom] = double(c[1]) * sizescale * sizescale * sizescale;
        count[iatom]++;
      }
      // checks
      for (int i = 0; i < natoms; i++) {
        if (count[i] == 0) {
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: atom {} missing in \"masses\" JSON section", id, i + 1);
        }
        if (rmass[i] <= 0.0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Invalid atom mass {} for atom {} in molecule JSON data",
                     rmass[i] / sizescale / sizescale / sizescale, i + 1);
      }
    } else {
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: Expected \"masses\" format [\"atom-id\",\"mass\"] but "
                 "found [\"{}\",\"{}\"]",
                 id, secfmt[0], secfmt[1]);
    }
  }

  // bonds

  if (bondflag) {
    int itype;
    tagint m, atom1, atom2;
    const int newton_bond = force->newton_bond;

    // must loop over data twice: first time to count, second time to apply

    for (int flag = 0; flag < 2; ++flag) {
      secfmt.clear();
      for (int i = 0; i < 3; ++i) secfmt.emplace_back(moldata["bonds"]["format"][i]);
      if ((secfmt[0] == "bond-type") && (secfmt[1] == "atom1") && (secfmt[2] == "atom2")) {

        if (flag == 0) {
          memset(count, 0, natoms * sizeof(int));
        } else {
          // must reallocate here in second iteration because bond_per_atom was not set for allocate() .
          memory->destroy(bond_type);
          memory->destroy(bond_atom);
          memory->create(bond_type, natoms, bond_per_atom, "molecule:bond_type");
          memory->create(bond_atom, natoms, bond_per_atom, "molecule:bond_atom");

          memset(num_bond, 0, natoms * sizeof(int));
        }

        for (int i = 0; i < nbonds; ++i) {
          const auto &item = moldata["bonds"]["data"][i];
          if (item.size() < 3)
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: invalid format of JSON data for bond {}: {}", id,
                       i + 1, to_string(item));

          if (item[0].is_number_integer()) {    // numeric type
            itype = int(item[0]) + boffset;
          } else {
            const auto &typestr = std::string(item[0]);
            if (!atom->labelmapflag)
              error->all(FLERR, Error::NOLASTLINE,
                         "Molecule template {}: invalid bond type in \"bonds\" JSON section", id,
                         typestr);
            itype = atom->lmap->find(typestr, Atom::BOND);
            if (itype == -1)
              error->all(FLERR, Error::NOLASTLINE,
                         "Molecule template {}: Unknown bond type {} in \"bonds\" JSON section", id,
                         typestr);
          }

          atom1 = tagint(item[1]);
          atom2 = tagint(item[2]);
          if ((atom1 <= 0) || (atom1 > natoms) || (atom2 <= 0) || (atom2 > natoms) ||
              (atom1 == atom2))
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: invalid atom ID in bond {}: {}", id, i + 1,
                       to_string(item));
          if ((itype <= 0) || (domain->box_exist && (itype > atom->nbondtypes)))
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: invalid bond type in bond {}: {}", id, i + 1,
                       to_string(item));
          if (flag == 0) {
            count[atom1 - 1]++;
            if (newton_bond == 0) count[atom2 - 1]++;
          } else {
            m = atom1 - 1;
            nbondtypes = MAX(nbondtypes, itype);
            bond_type[m][num_bond[m]] = itype;
            bond_atom[m][num_bond[m]] = atom2;
            num_bond[m]++;
            if (newton_bond == 0) {
              m = atom2 - 1;
              bond_type[m][num_bond[m]] = itype;
              bond_atom[m][num_bond[m]] = atom1;
              num_bond[m]++;
            }
          }
        }

        // bond_per_atom = max of count vector

        if (flag == 0) {
          bond_per_atom = 0;
          for (int i = 0; i < natoms; i++) bond_per_atom = MAX(bond_per_atom, count[i]);
        }
      }
    }
  }

  // angles

  if (angleflag) {
    int itype;
    tagint m, atom1, atom2, atom3;
    const int newton_bond = force->newton_bond;

    // must loop over data twice: first time to count, second time to apply

    for (int flag = 0; flag < 2; ++flag) {
      secfmt.clear();
      for (int i = 0; i < 4; ++i) secfmt.emplace_back(moldata["angles"]["format"][i]);
      if ((secfmt[0] == "angle-type") && (secfmt[1] == "atom1") && (secfmt[2] == "atom2") &&
          (secfmt[3] == "atom3")) {

        if (flag == 0) {
          memset(count, 0, natoms * sizeof(int));
        } else {
          // must reallocate here in second iteration because angle_per_atom was not set for allocate() .
          memory->destroy(angle_type);
          memory->destroy(angle_atom1);
          memory->destroy(angle_atom2);
          memory->destroy(angle_atom3);
          memory->create(angle_type, natoms, angle_per_atom, "molecule:angle_type");
          memory->create(angle_atom1, natoms, angle_per_atom, "molecule:angle_atom1");
          memory->create(angle_atom2, natoms, angle_per_atom, "molecule:angle_atom2");
          memory->create(angle_atom3, natoms, angle_per_atom, "molecule:angle_atom3");

          memset(num_angle, 0, natoms * sizeof(int));
        }

        for (int i = 0; i < nangles; ++i) {
          const auto &item = moldata["angles"]["data"][i];
          if (item.size() < 4)
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: invalid format of JSON data for angle {}: {}", id,
                       i + 1, to_string(item));

          if (item[0].is_number_integer()) {    // numeric type
            itype = int(item[0]) + aoffset;
          } else {
            const auto &typestr = std::string(item[0]);
            if (!atom->labelmapflag)
              error->all(FLERR, Error::NOLASTLINE,
                         "Molecule template {}: invalid angle type in \"angles\" JSON section", id,
                         typestr);
            itype = atom->lmap->find(typestr, Atom::ANGLE);
            if (itype == -1)
              error->all(FLERR, Error::NOLASTLINE,
                         "Molecule template {}: Unknown angle type {} in \"angles\" JSON section",
                         id, typestr);
          }

          atom1 = tagint(item[1]);
          atom2 = tagint(item[2]);
          atom3 = tagint(item[3]);

          if ((atom1 <= 0) || (atom1 > natoms) || (atom2 <= 0) || (atom2 > natoms) ||
              (atom3 <= 0) || (atom3 > natoms) || (atom1 == atom2) || (atom1 == atom3) ||
              (atom2 == atom3))
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: invalid atom ID in angle {}: {}", id, i + 1,
                       to_string(item));
          if ((itype <= 0) || (domain->box_exist && (itype > atom->nangletypes)))
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: invalid angle type in angle {}: {}", id, i + 1,
                       to_string(item));
          if (flag == 0) {
            count[atom2 - 1]++;
            if (newton_bond == 0) {
              count[atom1 - 1]++;
              count[atom3 - 1]++;
            }
          } else {
            m = atom2 - 1;
            nangletypes = MAX(nangletypes, itype);
            angle_type[m][num_angle[m]] = itype;
            angle_atom1[m][num_angle[m]] = atom1;
            angle_atom2[m][num_angle[m]] = atom2;
            angle_atom3[m][num_angle[m]] = atom3;
            num_angle[m]++;
            if (newton_bond == 0) {
              m = atom1 - 1;
              angle_type[m][num_angle[m]] = itype;
              angle_atom1[m][num_angle[m]] = atom1;
              angle_atom2[m][num_angle[m]] = atom2;
              angle_atom3[m][num_angle[m]] = atom3;
              num_angle[m]++;
              m = atom3 - 1;
              angle_type[m][num_angle[m]] = itype;
              angle_atom1[m][num_angle[m]] = atom1;
              angle_atom2[m][num_angle[m]] = atom2;
              angle_atom3[m][num_angle[m]] = atom3;
              num_angle[m]++;
            }
          }
        }

        // angle_per_atom = max of count vector

        if (flag == 0) {
          angle_per_atom = 0;
          for (int i = 0; i < natoms; i++) angle_per_atom = MAX(angle_per_atom, count[i]);
        }
      }
    }
  }

  // dihedrals

  if (dihedralflag) {
    int itype;
    tagint m, atom1, atom2, atom3, atom4;
    const int newton_bond = force->newton_bond;

    // must loop over data twice: first time to count, second time to apply

    for (int flag = 0; flag < 2; ++flag) {
      secfmt.clear();
      for (int i = 0; i < 5; ++i) secfmt.emplace_back(moldata["dihedrals"]["format"][i]);
      if ((secfmt[0] == "dihedral-type") && (secfmt[1] == "atom1") && (secfmt[2] == "atom2") &&
          (secfmt[3] == "atom3") && (secfmt[4] == "atom4")) {

        if (flag == 0) {
          memset(count, 0, natoms * sizeof(int));
        } else {
          // must reallocate here in second iteration because dihedral_per_atom was not set for allocate() .
          memory->destroy(dihedral_type);
          memory->destroy(dihedral_atom1);
          memory->destroy(dihedral_atom2);
          memory->destroy(dihedral_atom3);
          memory->destroy(dihedral_atom4);
          memory->create(dihedral_type, natoms, dihedral_per_atom, "molecule:dihedral_type");
          memory->create(dihedral_atom1, natoms, dihedral_per_atom, "molecule:dihedral_atom1");
          memory->create(dihedral_atom2, natoms, dihedral_per_atom, "molecule:dihedral_atom2");
          memory->create(dihedral_atom3, natoms, dihedral_per_atom, "molecule:dihedral_atom3");
          memory->create(dihedral_atom4, natoms, dihedral_per_atom, "molecule:dihedral_atom4");

          memset(num_dihedral, 0, natoms * sizeof(int));
        }

        for (int i = 0; i < ndihedrals; ++i) {
          const auto &item = moldata["dihedrals"]["data"][i];
          if (item.size() < 4)
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: invalid format of JSON data for dihedral {}: {}", id,
                       i + 1, to_string(item));

          if (item[0].is_number_integer()) {    // numeric type
            itype = int(item[0]) + doffset;
          } else {
            const auto &typestr = std::string(item[0]);
            if (!atom->labelmapflag)
              error->all(
                  FLERR, Error::NOLASTLINE,
                  "Molecule template {}: invalid dihedral type in \"dihedrals\" JSON section", id,
                  typestr);
            itype = atom->lmap->find(typestr, Atom::DIHEDRAL);
            if (itype == -1)
              error->all(
                  FLERR, Error::NOLASTLINE,
                  "Molecule template {}: Unknown dihedral type {} in \"dihedrals\" JSON section",
                  id, typestr);
          }

          atom1 = tagint(item[1]);
          atom2 = tagint(item[2]);
          atom3 = tagint(item[3]);
          atom4 = tagint(item[4]);

          if ((atom1 <= 0) || (atom1 > natoms) || (atom2 <= 0) || (atom2 > natoms) ||
              (atom3 <= 0) || (atom3 > natoms) || (atom4 <= 0) || (atom4 > natoms) ||
              (atom1 == atom2) || (atom1 == atom3) || (atom1 == atom4) || (atom2 == atom3) ||
              (atom2 == atom4) || (atom3 == atom4))
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: invalid atom ID in dihedral {}: {}", id, i + 1,
                       to_string(item));
          if ((itype <= 0) || (domain->box_exist && (itype > atom->ndihedraltypes)))
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: invalid dihedral type in dihedral {}: {}", id, i + 1,
                       to_string(item));
          if (flag == 0) {
            count[atom2 - 1]++;
            if (newton_bond == 0) {
              count[atom1 - 1]++;
              count[atom3 - 1]++;
              count[atom4 - 1]++;
            }
          } else {
            m = atom2 - 1;
            ndihedraltypes = MAX(ndihedraltypes, itype);
            dihedral_type[m][num_dihedral[m]] = itype;
            dihedral_atom1[m][num_dihedral[m]] = atom1;
            dihedral_atom2[m][num_dihedral[m]] = atom2;
            dihedral_atom3[m][num_dihedral[m]] = atom3;
            dihedral_atom4[m][num_dihedral[m]] = atom4;
            num_dihedral[m]++;
            if (newton_bond == 0) {
              m = atom1 - 1;
              dihedral_type[m][num_dihedral[m]] = itype;
              dihedral_atom1[m][num_dihedral[m]] = atom1;
              dihedral_atom2[m][num_dihedral[m]] = atom2;
              dihedral_atom3[m][num_dihedral[m]] = atom3;
              dihedral_atom4[m][num_dihedral[m]] = atom4;
              num_dihedral[m]++;
              m = atom3 - 1;
              dihedral_type[m][num_dihedral[m]] = itype;
              dihedral_atom1[m][num_dihedral[m]] = atom1;
              dihedral_atom2[m][num_dihedral[m]] = atom2;
              dihedral_atom3[m][num_dihedral[m]] = atom3;
              dihedral_atom4[m][num_dihedral[m]] = atom4;
              num_dihedral[m]++;
              m = atom4 - 1;
              dihedral_type[m][num_dihedral[m]] = itype;
              dihedral_atom1[m][num_dihedral[m]] = atom1;
              dihedral_atom2[m][num_dihedral[m]] = atom2;
              dihedral_atom3[m][num_dihedral[m]] = atom3;
              dihedral_atom4[m][num_dihedral[m]] = atom4;
              num_dihedral[m]++;
            }
          }
        }

        // dihedral_per_atom = max of count vector

        if (flag == 0) {
          dihedral_per_atom = 0;
          for (int i = 0; i < natoms; i++) dihedral_per_atom = MAX(dihedral_per_atom, count[i]);
        }
      }
    }
  }

  // impropers

  if (improperflag) {
    int itype;
    tagint m, atom1, atom2, atom3, atom4;
    const int newton_bond = force->newton_bond;

    // must loop over data twice: first time to count, second time to apply

    for (int flag = 0; flag < 2; ++flag) {
      secfmt.clear();
      for (int i = 0; i < 5; ++i) secfmt.emplace_back(moldata["impropers"]["format"][i]);
      if ((secfmt[0] == "improper-type") && (secfmt[1] == "atom1") && (secfmt[2] == "atom2") &&
          (secfmt[3] == "atom3") && (secfmt[4] == "atom4")) {

        if (flag == 0) {
          memset(count, 0, natoms * sizeof(int));
        } else {
          // must reallocate here in second iteration because improper_per_atom was not set for allocate() .
          memory->destroy(improper_type);
          memory->destroy(improper_atom1);
          memory->destroy(improper_atom2);
          memory->destroy(improper_atom3);
          memory->destroy(improper_atom4);
          memory->create(improper_type, natoms, improper_per_atom, "molecule:improper_type");
          memory->create(improper_atom1, natoms, improper_per_atom, "molecule:improper_atom1");
          memory->create(improper_atom2, natoms, improper_per_atom, "molecule:improper_atom2");
          memory->create(improper_atom3, natoms, improper_per_atom, "molecule:improper_atom3");
          memory->create(improper_atom4, natoms, improper_per_atom, "molecule:improper_atom4");

          memset(num_improper, 0, natoms * sizeof(int));
        }

        for (int i = 0; i < nimpropers; ++i) {
          const auto &item = moldata["impropers"]["data"][i];
          if (item.size() < 4)
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: invalid format of JSON data for improper {}: {}", id,
                       i + 1, to_string(item));

          if (item[0].is_number_integer()) {    // numeric type
            itype = int(item[0]) + aoffset;
          } else {
            const auto &typestr = std::string(item[0]);
            if (!atom->labelmapflag)
              error->all(
                  FLERR, Error::NOLASTLINE,
                  "Molecule template {}: invalid improper type in \"impropers\" JSON section", id,
                  typestr);
            itype = atom->lmap->find(typestr, Atom::IMPROPER);
            if (itype == -1)
              error->all(
                  FLERR, Error::NOLASTLINE,
                  "Molecule template {}: Unknown improper type {} in \"impropers\" JSON section",
                  id, typestr);
          }

          atom1 = tagint(item[1]);
          atom2 = tagint(item[2]);
          atom3 = tagint(item[3]);
          atom4 = tagint(item[4]);

          if ((atom1 <= 0) || (atom1 > natoms) || (atom2 <= 0) || (atom2 > natoms) ||
              (atom3 <= 0) || (atom3 > natoms) || (atom4 <= 0) || (atom4 > natoms) ||
              (atom1 == atom2) || (atom1 == atom3) || (atom1 == atom4) || (atom2 == atom3) ||
              (atom2 == atom4) || (atom3 == atom4))
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: invalid atom ID in improper {}: {}", id, i + 1,
                       to_string(item));
          if ((itype <= 0) || (domain->box_exist && (itype > atom->nimpropertypes)))
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: invalid improper type in improper {}: {}", id, i + 1,
                       to_string(item));
          if (flag == 0) {
            count[atom2 - 1]++;
            if (newton_bond == 0) {
              count[atom1 - 1]++;
              count[atom3 - 1]++;
              count[atom4 - 1]++;
            }
          } else {
            m = atom2 - 1;
            nimpropertypes = MAX(nimpropertypes, itype);
            improper_type[m][num_improper[m]] = itype;
            improper_atom1[m][num_improper[m]] = atom1;
            improper_atom2[m][num_improper[m]] = atom2;
            improper_atom3[m][num_improper[m]] = atom3;
            improper_atom4[m][num_improper[m]] = atom4;
            num_improper[m]++;
            if (newton_bond == 0) {
              m = atom1 - 1;
              improper_type[m][num_improper[m]] = itype;
              improper_atom1[m][num_improper[m]] = atom1;
              improper_atom2[m][num_improper[m]] = atom2;
              improper_atom3[m][num_improper[m]] = atom3;
              improper_atom4[m][num_improper[m]] = atom4;
              num_improper[m]++;
              m = atom3 - 1;
              improper_type[m][num_improper[m]] = itype;
              improper_atom1[m][num_improper[m]] = atom1;
              improper_atom2[m][num_improper[m]] = atom2;
              improper_atom3[m][num_improper[m]] = atom3;
              improper_atom4[m][num_improper[m]] = atom4;
              num_improper[m]++;
              m = atom4 - 1;
              improper_type[m][num_improper[m]] = itype;
              improper_atom1[m][num_improper[m]] = atom1;
              improper_atom2[m][num_improper[m]] = atom2;
              improper_atom3[m][num_improper[m]] = atom3;
              improper_atom4[m][num_improper[m]] = atom4;
              num_improper[m]++;
            }
          }
        }

        // improper_per_atom = max of count vector

        if (flag == 0) {
          improper_per_atom = 0;
          for (int i = 0; i < natoms; i++) improper_per_atom = MAX(improper_per_atom, count[i]);
        }
      }
    }
  }

  if (specialflag) {

    // process counts

    const auto &specialcounts = moldata["special"]["counts"];
    secfmt.clear();
    for (int i = 0; i < 4; ++i) secfmt.emplace_back(specialcounts["format"][i]);
    if ((secfmt[0] == "atom-id") && (secfmt[1] == "n12") && (secfmt[2] == "n13") &&
        (secfmt[3] == "n14")) {

      memset(count, 0, natoms * sizeof(int));
      for (const auto &item : specialcounts["data"]) {
        if (!item[0].is_number_integer() || !item[1].is_number_integer() ||
            !item[2].is_number_integer() || !item[3].is_number_integer())
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid data in \"special:counts\" section of molecule "
                     "JSON data: {}",
                     id, to_string(item));
        const auto &vals = item.get<std::vector<int>>();
        const int iatom = vals[0] - 1;
        if ((iatom < 0) || (iatom >= natoms))
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid atom-id {} in \"special:counts\" section of "
                     "molecule JSON data",
                     id, iatom + 1);

        nspecial[iatom][0] = vals[1];
        nspecial[iatom][1] = vals[1] + vals[2];
        nspecial[iatom][2] = vals[1] + vals[2] + vals[3];
        count[iatom]++;
      }
      // check
      for (int i = 0; i < natoms; i++) {
        if (count[i] == 0) {
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: atom {} missing in \"special:counts\" JSON section", id,
                     i + 1);
        }
      }
    } else {
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: Expected \"special:counts\" format "
                 "[\"atom-id\",\"n12\",\"n13\",\"n14\"] but found [\"{}\",\"{}\",\"{}\",\"{}\"]",
                 id, secfmt[0], secfmt[1], secfmt[2], secfmt[3]);
    }

    // process bonds

    const auto &specialbonds = moldata["special"]["bonds"];
    secfmt.clear();
    for (int i = 0; i < 2; ++i) secfmt.emplace_back(specialbonds["format"][i]);
    if ((secfmt[0] == "atom-id") && (secfmt[1] == "atom-id-list")) {
      memset(count, 0, natoms * sizeof(int));
      for (int i = 0; i < natoms; ++i) {
        if (!specialbonds["data"][i][0].is_number_integer())
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid atom-id {} for entry {} in \"special:bonds\" "
                     "section of molecule JSON data",
                     id, to_string(specialbonds["data"][i][0]), i + 1);
        const int iatom = int(specialbonds["data"][i][0]) - 1;
        if ((iatom < 0) || (iatom >= natoms))
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid atom-id {} in \"special:bondss\" section of "
                     "molecule JSON data",
                     id, iatom + 1);

        int m = 0;
        for (const auto &item : specialbonds["data"][i][1]) {
          if (!item.is_number_integer())
            error->all(
                FLERR, Error::NOLASTLINE,
                "Molecule template {}: invalid data in \"special:bonds\" section of molecule "
                "JSON data: {}",
                id, to_string(specialbonds["data"][i][1]));

          auto ival = tagint(item);
          if ((ival <= 0) || (ival > natoms) || (ival == iatom + 1))
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: invalid atom index {} in \"special:bonds\" section "
                       "of JSON data",
                       id, ival);
          special[iatom][m++] = tagint(item);
        }
        count[iatom]++;
      }
      // check
      for (int i = 0; i < natoms; i++) {
        if (count[i] == 0) {
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: atom {} missing in \"special:bonds\" JSON section", id,
                     i + 1);
        }
      }
    } else {
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: Expected \"special:bonds\" format "
                 "[\"atom-id\",\"atom-id-list\"] but found [\"{}\",\"{}\"]",
                 id, secfmt[0], secfmt[1]);
    }
  }

  // shake settings

  if (shakeflagflag) {
    const auto &shakedata = moldata["shake"]["flags"];
    secfmt.clear();
    for (int i = 0; i < 2; ++i) secfmt.emplace_back(shakedata["format"][i]);
    if ((secfmt[0] == "atom-id") && (secfmt[1] == "flag")) {

      for (int i = 0; i < natoms; i++) shake_flag[i] = -1;
      memset(count, 0, natoms * sizeof(int));
      for (const auto &c : shakedata["data"]) {
        if (c.size() < 2)
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: missing data in \"shake:flags\" section of molecule "
                     "JSON data: {}",
                     id, to_string(c));
        if (!c[0].is_number_integer())
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid atom-id in \"shake:flags\" section of molecule "
                     "JSON data: {}",
                     id, to_string(c[0]));

        const int iatom = int(c[0]) - 1;
        if ((iatom < 0) || (iatom >= natoms))
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid atom-id {} in \"shake:flags\" section of "
                     "molecule JSON data",
                     id, iatom + 1);
        if (!c[1].is_number_integer())
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid flag in \"shake:flags\" section of "
                     "molecule JSON data: {}",
                     id, to_string(c[1]));
        shake_flag[iatom] = int(c[1]);
        count[iatom]++;
      }
      // checks
      for (int i = 0; i < natoms; i++) {
        if (count[i] == 0) {
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: atom {} missing in \"shake:flags\" JSON section", id,
                     i + 1);
        }
        if ((shake_flag[i] < 0) || (shake_flag[i] > 4))
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid flag value {} in \"shake:flags\" section of "
                     "molecule JSON data",
                     id, shake_flag[i]);
      }
    } else {
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: Expected \"shake:flags\" format [\"atom-id\",\"mass\"] but "
                 "found [\"{}\",\"{}\"]",
                 id, secfmt[0], secfmt[1]);
    }
  }

  // shake_atoms

#define APPLY_SHAKE_ATOMS(ncols)                                                          \
  if (c[1].size() != ncols)                                                               \
    error->all(FLERR, Error::NOLASTLINE,                                                  \
               "Molecule template {}: invalid number of items for atom-id {} in \"shake:" \
               "atoms\" section of molecue JSON data ({} vs {})",                         \
               id, iatom + 1, c[1].size(), ncols);                                        \
  for (int i = 0; i < ncols; ++i) {                                                       \
    if (!c[1][i].is_number_integer())                                                     \
      error->all(FLERR, Error::NOLASTLINE,                                                \
                 "Molecule template {}: invalid atom-id {} in atom-id-list for atom {} "  \
                 "in \"shake:atoms\" section of molecule JSON data",                      \
                 id, to_string(c[1][i]), iatom + 1);                                      \
    shake_atom[iatom][i] = int(c[1][i]);                                                  \
  }

  if (shakeatomflag) {
    const auto &shakedata = moldata["shake"]["atoms"];
    secfmt.clear();
    for (int i = 0; i < 2; ++i) secfmt.emplace_back(shakedata["format"][i]);
    if ((secfmt[0] == "atom-id") && (secfmt[1] == "atom-id-list")) {

      memset(count, 0, natoms * sizeof(int));
      for (const auto &c : shakedata["data"]) {
        if (c.size() < 2)
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: missing data in \"shake:atoms\" section of molecule "
                     "JSON data: {}",
                     id, to_string(c));
        if (!c[0].is_number_integer())
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid atom-id in \"shake:atoms\" section of molecule "
                     "JSON data: {}",
                     id, to_string(c[0]));

        const int iatom = int(c[0]) - 1;
        if ((iatom < 0) || (iatom >= natoms))
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid atom-id {} in \"shake:atoms\" section of "
                     "molecule JSON data",
                     id, iatom + 1);

        switch (shake_flag[iatom]) {
          case 1:
            APPLY_SHAKE_ATOMS(3);
            break;
          case 2:
            APPLY_SHAKE_ATOMS(2);
            break;
          case 3:
            APPLY_SHAKE_ATOMS(3);
            break;
          case 4:
            APPLY_SHAKE_ATOMS(4);
            break;
          case 0:
            APPLY_SHAKE_ATOMS(0);
            break;
          default:
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: Unsupported Shake flag {} for "
                       " atom {} in \"shake:atoms\" section of molecule JSON data",
                       id, shake_flag[iatom], iatom + 1);
        }
        count[iatom]++;
      }
      // checks
      for (int i = 0; i < natoms; i++) {
        if (count[i] == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: atom {} missing in \"shake:atoms\" JSON section", id,
                     i + 1);
      }
    } else {
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: Expected \"shake:atoms\" format "
                 "[\"atom-id\",\"atom-id-list\"] but found [\"{}\",\"{}\"]",
                 id, secfmt[0], secfmt[1]);
    }
  }
#undef APPLY_SHAKE_ATOMS

  // shake_bond_types

#define SET_SHAKE_TYPE(type, idx, ncols, offset)                                          \
  if (c[1].size() < ncols)                                                                \
    error->all(FLERR, Error::NOLASTLINE,                                                  \
               "Molecule template {}: invalid number of items for atom-id {} in \"shake:" \
               "types\" section of molecue JSON data ({} vs {})",                         \
               id, iatom + 1, c[1].size(), ncols);                                        \
  if (c[1][idx].is_number_integer()) {                                                    \
    shake_type[iatom][idx] = int(c[1][idx]) + offset;                                     \
  } else {                                                                                \
    char *subst = utils::expand_type(FLERR, c[1][idx], type, lmp);                        \
    if (subst) {                                                                          \
      shake_type[iatom][idx] = utils::inumeric(FLERR, subst, false, lmp);                 \
      delete[] subst;                                                                     \
    }                                                                                     \
  }

  if (shaketypeflag) {
    const auto &shakedata = moldata["shake"]["types"];
    secfmt.clear();
    for (int i = 0; i < 2; ++i) secfmt.emplace_back(shakedata["format"][i]);
    if ((secfmt[0] == "atom-id") && (secfmt[1] == "type-list")) {

      memset(count, 0, natoms * sizeof(int));
      for (const auto &c : shakedata["data"]) {
        if (c.size() < 2)
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: missing data in \"shake:types\" section of molecule "
                     "JSON data: {}",
                     id, to_string(c));
        if (!c[0].is_number_integer())
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid atom-id in \"shake:types\" section of molecule "
                     "JSON data: {}",
                     id, to_string(c[0]));

        const int iatom = int(c[0]) - 1;
        if ((iatom < 0) || (iatom >= natoms))
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: invalid atom-id {} in \"shake:types\" section of "
                     "molecule JSON data",
                     id, iatom + 1);

        switch (shake_flag[iatom]) {
          case 1:
            SET_SHAKE_TYPE(Atom::BOND, 0, 3, boffset);
            SET_SHAKE_TYPE(Atom::BOND, 1, 3, boffset);
            SET_SHAKE_TYPE(Atom::ANGLE, 2, 3, aoffset);
            break;
          case 2:
            SET_SHAKE_TYPE(Atom::BOND, 0, 1, boffset);
            break;
          case 3:
            SET_SHAKE_TYPE(Atom::BOND, 0, 2, boffset);
            SET_SHAKE_TYPE(Atom::BOND, 1, 2, boffset);
            break;
          case 4:
            SET_SHAKE_TYPE(Atom::BOND, 0, 3, boffset);
            SET_SHAKE_TYPE(Atom::BOND, 1, 3, boffset);
            SET_SHAKE_TYPE(Atom::BOND, 2, 3, boffset);
            break;
          case 0:
            break;
          default:
            error->all(FLERR, Error::NOLASTLINE,
                       "Molecule template {}: Unsupported Shake flag {} for "
                       " atom {} in \"shake:types\" section of molecule JSON data",
                       id, shake_flag[iatom], iatom + 1);
        }
        count[iatom]++;
      }
      // checks
      for (int i = 0; i < natoms; i++) {
        if (count[i] == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Molecule template {}: atom {} missing in \"shake:types\" JSON section", id,
                     i + 1);
      }
    } else {
      error->all(FLERR, Error::NOLASTLINE,
                 "Molecule template {}: Expected \"shake:types\" format "
                 "[\"atom-id\",\"type-list\"] but found [\"{}\",\"{}\"]",
                 id, secfmt[0], secfmt[1]);
    }
  }
#undef SET_SHAKE_TYPE

  // body integers and doubles

  if (bodyflag) {
    for (int i = 0; i < nibody; ++i) ibodyparams[i] = moldata["body"]["integers"][i];
    for (int i = 0; i < ndbody; ++i) dbodyparams[i] = moldata["body"]["doubles"][i];
  }

  // error checks

  if (specialflag && !bondflag)
    error->all(FLERR, fileiarg, "Molecule file has special flags but no bonds");
  if ((shakeflagflag || shakeatomflag || shaketypeflag) && !shakeflag)
    error->all(FLERR, Error::NOLASTLINE,
               "Molecule template {}: \"shake\" info is incomplete in JSON data");
  if (bodyflag && !rmassflag)
    error->all(FLERR, Error::NOLASTLINE,
               "Molecule template {}: \"body\" JSON section requires \"masses\" section", id);

  // auto-generate special bonds if needed and not in file

  if (bondflag && specialflag == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR, fileiarg,
                 "Cannot auto-generate special bonds before simulation box is defined");

    special_generate();
    specialflag_user = 0;
    specialflag = 1;
    nspecialflag = 1;
  }

  // body particle must have natom = 1
  // set radius by having body class compute its own radius

  if (bodyflag) {
    radiusflag = 1;
    if (natoms != 1) error->all(FLERR, fileiarg, "Molecule natoms must be 1 for body particle");
    if (sizescale != 1.0)
      error->all(FLERR, fileiarg, "Molecule sizescale must be 1.0 for body particle");
    radius[0] = avec_body->radius_body(nibody, ndbody, ibodyparams, dbodyparams);
    maxradius = radius[0];
  }

  // clean up

  memory->destroy(count);
}

// clang-format on
// ------------------------------------------------------------------------------
//  convert json data structure to molecule data structure
// ------------------------------------------------------------------------------

json Molecule::to_json() const
{
  json moldata;

  // required global settings
  moldata["application"] = "LAMMPS";
  moldata["format"] = "molecule";
  moldata["revision"] = 1;
  moldata["schema"] = "https://download.lammps.org/json/molecule-schema.json";
  moldata["title"] = title;
  moldata["units"] = update->unit_style;

  // optional global settings
  if (massflag_user) moldata["masstotal"] = masstotal;
  if (comflag_user) {
    moldata["com"][0] = com[0];
    moldata["com"][1] = com[1];
    moldata["com"][2] = com[2];
  }
  if (inertiaflag_user) {
    moldata["inertia"][0] = itensor[0];
    moldata["inertia"][1] = itensor[1];
    moldata["inertia"][2] = itensor[2];
    moldata["inertia"][3] = itensor[3];
    moldata["inertia"][4] = itensor[4];
    moldata["inertia"][5] = itensor[5];
  }

  // fields with format
  if (xflag) {
    moldata["coords"]["format"] = {"atom-id", "x", "y", "z"};
    for (int i = 0; i < natoms; ++i) {
      moldata["coords"]["data"][i] = {i + 1, x[i][0], x[i][1], x[i][2]};
    }
  }

  if (typeflag) {
    moldata["types"]["format"] = {"atom-id", "type"};
    if (atom->labelmapflag && atom->lmap->is_complete(Atom::ATOM)) {
      for (int i = 0; i < natoms; ++i) {
        moldata["types"]["data"][i] = {i + 1, atom->lmap->find(type[i], Atom::ATOM)};
      }
    } else {
      for (int i = 0; i < natoms; ++i) moldata["types"]["data"][i] = {i + 1, type[i]};
    }
  }

  if (moleculeflag) {
    moldata["molecules"]["format"] = {"atom-id", "molecule-id"};
    for (int i = 0; i < natoms; ++i) moldata["molecules"]["data"][i] = {i + 1, molecule[i]};
  }

  if (fragmentflag) {
    moldata["fragments"]["format"] = {"fragment-id", "atom-id-list"};
    for (int i = 0; i < nfragments; ++i) {
      moldata["fragments"]["data"][i][0] = fragmentnames[i];
      int k = 0;
      for (int j = 0; j < natoms; ++j) {
        if (fragmentmask[i][j]) {
          moldata["fragments"]["data"][i][1][k] = j + 1;
          ++k;
        }
      }
    }
  }

  if (qflag) {
    moldata["charges"]["format"] = {"atom-id", "charge"};
    for (int i = 0; i < natoms; ++i) moldata["charges"]["data"][i] = {i + 1, q[i]};
  }

  if (radiusflag) {
    moldata["diameters"]["format"] = {"atom-id", "diameter"};
    for (int i = 0; i < natoms; ++i) moldata["diameters"]["data"][i] = {i + 1, 2.0 * radius[i]};
  }

  if (muflag) {
    moldata["dipoles"]["format"] = {"atom-id", "mux", "muy", "muz"};
    for (int i = 0; i < natoms; ++i)
      moldata["dipoles"]["data"][i] = {i + 1, mu[i][0], mu[i][1], mu[i][2]};
  }

  if (rmassflag) {
    moldata["masses"]["format"] = {"atom-id", "mass"};
    for (int i = 0; i < natoms; ++i) moldata["masses"]["data"][i] = {i + 1, rmass[i]};
  }

  bool has_newton_bond = force->newton_bond > 0;

  if (bondflag) {
    int idx = 0;
    bool has_typelabels = (atom->labelmapflag != 0) && atom->lmap->is_complete(Atom::BOND);
    moldata["bonds"]["format"] = {"bond-type", "atom1", "atom2"};
    for (int i = 0; i < natoms; i++) {
      for (int j = 0; j < num_bond[i]; j++) {
        if (has_newton_bond || (i + 1 < bond_atom[i][j])) {
          if (has_typelabels) {
            moldata["bonds"]["data"][idx] = {atom->lmap->find(bond_type[i][j], Atom::BOND), i + 1,
                                             bond_atom[i][j]};
          } else {
            moldata["bonds"]["data"][idx] = {bond_type[i][j], i + 1, bond_atom[i][j]};
          }
          ++idx;
        }
      }
    }
  }

  if (angleflag) {
    int idx = 0;
    bool has_typelabels = (atom->labelmapflag != 0) && atom->lmap->is_complete(Atom::ANGLE);
    moldata["angles"]["format"] = {"angle-type", "atom1", "atom2", "atom3"};
    for (int i = 0; i < natoms; i++) {
      for (int j = 0; j < num_angle[i]; j++) {
        if (has_newton_bond || (i + 1 == angle_atom2[i][j])) {
          if (has_typelabels) {
            moldata["angles"]["data"][idx] = {atom->lmap->find(angle_type[i][j], Atom::ANGLE),
                                              angle_atom1[i][j], angle_atom2[i][j],
                                              angle_atom3[i][j]};
          } else {
            moldata["angles"]["data"][idx] = {angle_type[i][j], angle_atom1[i][j],
                                              angle_atom2[i][j], angle_atom3[i][j]};
          }
          ++idx;
        }
      }
    }
  }

  if (dihedralflag) {
    int idx = 0;
    bool has_typelabels = (atom->labelmapflag != 0) && atom->lmap->is_complete(Atom::DIHEDRAL);
    moldata["dihedrals"]["format"] = {"dihedral-type", "atom1", "atom2", "atom3", "atom4"};
    for (int i = 0; i < natoms; i++) {
      for (int j = 0; j < num_dihedral[i]; j++) {
        if (has_newton_bond || (i + 1 == dihedral_atom2[i][j])) {
          if (has_typelabels) {
            moldata["dihedrals"]["data"][idx] = {
                atom->lmap->find(dihedral_type[i][j], Atom::DIHEDRAL), dihedral_atom1[i][j],
                dihedral_atom2[i][j], dihedral_atom3[i][j], dihedral_atom4[i][j]};
          } else {
            moldata["dihedrals"]["data"][idx] = {dihedral_type[i][j], dihedral_atom1[i][j],
                                                 dihedral_atom2[i][j], dihedral_atom3[i][j],
                                                 dihedral_atom4[i][j]};
          }
          ++idx;
        }
      }
    }
  }

  if (improperflag) {
    int idx = 0;
    bool has_typelabels = (atom->labelmapflag != 0) && atom->lmap->is_complete(Atom::IMPROPER);
    moldata["impropers"]["format"] = {"improper-type", "atom1", "atom2", "atom3", "atom4"};
    for (int i = 0; i < natoms; i++) {
      for (int j = 0; j < num_improper[i]; j++) {
        if (has_newton_bond || (i + 1 == improper_atom2[i][j])) {
          if (has_typelabels) {
            moldata["impropers"]["data"][idx] = {
                atom->lmap->find(improper_type[i][j], Atom::IMPROPER), improper_atom1[i][j],
                improper_atom2[i][j], improper_atom3[i][j], improper_atom4[i][j]};
          } else {
            moldata["impropers"]["data"][idx] = {improper_type[i][j], improper_atom1[i][j],
                                                 improper_atom2[i][j], improper_atom3[i][j],
                                                 improper_atom4[i][j]};
          }
          ++idx;
        }
      }
    }
  }

  if (specialflag_user) {
    moldata["special"]["counts"]["format"] = {"atom-id", "n12", "n13", "n14"};
    moldata["special"]["bonds"]["format"] = {"atom-id", "atom-id-list"};
    for (int i = 0; i < natoms; i++) {
      moldata["special"]["counts"]["data"][i] = {
          i + 1, nspecial[i][0], nspecial[i][1] - nspecial[i][0], nspecial[i][2] - nspecial[i][1]};
      moldata["special"]["bonds"]["data"][i][0] = i + 1;
      for (int j = 0; j < nspecial[i][2]; ++j)
        moldata["special"]["bonds"]["data"][i][1][j] = special[i][j];
    }
  }

  if (shakeflag) {
    moldata["shake"]["flags"]["format"] = {"atom-id", "flag"};
    moldata["shake"]["atoms"]["format"] = {"atom-id", "atom-id-list"};
    moldata["shake"]["types"]["format"] = {"atom-id", "type-list"};
    for (int i = 0; i < natoms; ++i) {
      bool has_typelabels = (atom->labelmapflag != 0) && atom->lmap->is_complete(Atom::BOND);
      moldata["shake"]["flags"]["data"][i] = {i + 1, shake_flag[i]};
      moldata["shake"]["atoms"]["data"][i][0] = i + 1;
      moldata["shake"]["types"]["data"][i][0] = i + 1;
      switch (shake_flag[i]) {
        case 1:
          has_typelabels = has_typelabels && atom->lmap->is_complete(Atom::ANGLE);
          moldata["shake"]["atoms"]["data"][i][1] = {shake_atom[i][0], shake_atom[i][1],
                                                     shake_atom[i][2]};
          if (has_typelabels) {
            moldata["shake"]["types"]["data"][i][1] = {
                atom->lmap->find(shake_type[i][0], Atom::BOND),
                atom->lmap->find(shake_type[i][1], Atom::BOND),
                atom->lmap->find(shake_type[i][2], Atom::ANGLE)};
          } else {
            moldata["shake"]["types"]["data"][i][1] = {shake_type[i][0], shake_type[i][1],
                                                       shake_type[i][2]};
          }
          break;
        case 2:
          moldata["shake"]["atoms"]["data"][i][1] = {shake_atom[i][0], shake_atom[i][1]};
          if (has_typelabels) {
            moldata["shake"]["types"]["data"][i][1] = {
                atom->lmap->find(shake_type[i][0], Atom::BOND),
                atom->lmap->find(shake_type[i][1], Atom::BOND)};
          } else {
            moldata["shake"]["types"]["data"][i][1] = {shake_type[i][0], shake_type[i][1]};
          }
          break;
        case 3:
          moldata["shake"]["atoms"]["data"][i][1] = {shake_atom[i][0], shake_atom[i][1],
                                                     shake_atom[i][2]};
          if (has_typelabels) {
            moldata["shake"]["types"]["data"][i][1] = {
                atom->lmap->find(shake_type[i][0], Atom::BOND),
                atom->lmap->find(shake_type[i][1], Atom::BOND),
                atom->lmap->find(shake_type[i][2], Atom::BOND)};
          } else {
            moldata["shake"]["types"]["data"][i][1] = {shake_type[i][0], shake_type[i][1],
                                                       shake_type[i][2]};
          }
          break;
        case 4:
          moldata["shake"]["atoms"]["data"][i][1] = {shake_atom[i][0], shake_atom[i][1],
                                                     shake_atom[i][2], shake_atom[i][3]};
          if (has_typelabels) {
            moldata["shake"]["types"]["data"][i][1] = {
                atom->lmap->find(shake_type[i][0], Atom::BOND),
                atom->lmap->find(shake_type[i][1], Atom::BOND),
                atom->lmap->find(shake_type[i][2], Atom::BOND),
                atom->lmap->find(shake_type[i][3], Atom::BOND)};
          } else {
            moldata["shake"]["types"]["data"][i][1] = {shake_type[i][0], shake_type[i][1],
                                                       shake_type[i][2], shake_type[i][3]};
          }
          break;
        case 0:
          moldata["shake"]["atoms"]["data"][i][1] = nullptr;
          moldata["shake"]["types"]["data"][i][1] = nullptr;
          break;
      }
    }
  }

  if (bodyflag) {
    for (int i = 0; i < nibody; ++i) { moldata["body"]["integers"][i] = ibodyparams[i]; }
    for (int i = 0; i < nibody; ++i) { moldata["body"]["doubles"][i] = dbodyparams[i]; }
  }

  return moldata;
}

/* ---------------------------------------------------------------------- */

Molecule::~Molecule()
{
  delete[] id;
  deallocate();
}

/* ----------------------------------------------------------------------
   compute center = geometric center of molecule
   also compute:
     dx = displacement of each atom from center
     molradius = radius of molecule from center
       including finite-size particles or body particles
------------------------------------------------------------------------- */

void Molecule::compute_center()
{
  if (centerflag) return;
  centerflag = 1;

  center[0] = center[1] = center[2] = 0.0;
  for (int i = 0; i < natoms; i++) {
    center[0] += x[i][0];
    center[1] += x[i][1];
    center[2] += x[i][2];
  }
  center[0] /= natoms;
  center[1] /= natoms;
  center[2] /= natoms;

  memory->destroy(dx);
  memory->create(dx, natoms, 3, "molecule:dx");

  for (int i = 0; i < natoms; i++) {
    dx[i][0] = x[i][0] - center[0];
    dx[i][1] = x[i][1] - center[1];
    dx[i][2] = x[i][2] - center[2];
  }

  molradius = 0.0;
  for (int i = 0; i < natoms; i++) {
    double rad = MathExtra::len3(dx[i]);
    if (radiusflag) rad += radius[i];
    molradius = MAX(molradius, rad);
  }
}

/* ----------------------------------------------------------------------
   compute masstotal = total mass of molecule
   could have been set by user, otherwise calculate it
------------------------------------------------------------------------- */

void Molecule::compute_mass()
{
  if (massflag) return;
  massflag = 1;

  atom->check_mass(FLERR);

  masstotal = 0.0;
  for (int i = 0; i < natoms; i++) {
    if (rmassflag)
      masstotal += rmass[i];
    else
      masstotal += atom->mass[type[i]];
  }
}

/* ----------------------------------------------------------------------
   compute com = center of mass of molecule
   could have been set by user, otherwise calculate it
   works for finite size particles assuming no overlap
   also compute:
     dxcom = displacement of each atom from COM
     comatom = which atom (1-Natom) is nearest the COM
     maxextent = furthest any atom in molecule is from comatom (not COM)
------------------------------------------------------------------------- */

void Molecule::compute_com()
{
  if (!comflag) {
    comflag = 1;

    atom->check_mass(FLERR);

    double onemass;
    com[0] = com[1] = com[2] = 0.0;
    for (int i = 0; i < natoms; i++) {
      if (rmassflag)
        onemass = rmass[i];
      else
        onemass = atom->mass[type[i]];
      com[0] += x[i][0] * onemass;
      com[1] += x[i][1] * onemass;
      com[2] += x[i][2] * onemass;
    }
    if (masstotal > 0.0) {
      com[0] /= masstotal;
      com[1] /= masstotal;
      com[2] /= masstotal;
    }
  }

  memory->destroy(dxcom);
  memory->create(dxcom, natoms, 3, "molecule:dxcom");

  for (int i = 0; i < natoms; i++) {
    dxcom[i][0] = x[i][0] - com[0];
    dxcom[i][1] = x[i][1] - com[1];
    dxcom[i][2] = x[i][2] - com[2];
  }

  double rsqmin = BIG;
  for (int i = 0; i < natoms; i++) {
    double rsq = MathExtra::lensq3(dxcom[i]);
    if (rsq < rsqmin) {
      comatom = i;
      rsqmin = rsq;
    }
  }

  double rsqmax = 0.0;
  for (int i = 0; i < natoms; i++) {
    double dx = x[comatom][0] - x[i][0];
    double dy = x[comatom][1] - x[i][1];
    double dz = x[comatom][2] - x[i][2];
    double rsq = dx * dx + dy * dy + dz * dz;
    rsqmax = MAX(rsqmax, rsq);
  }

  comatom++;
  maxextent = sqrt(rsqmax);
}

/* ----------------------------------------------------------------------
   compute itensor = 6 moments of inertia of molecule around xyz axes
   could have been set by user, otherwise calculate it
   accounts for finite size spheres, assuming no overlap
   also compute:
     inertia = 3 principal components of inertia
     ex,ey,ez = principal axes in space coords
     quat = quaternion for orientation of molecule
     dxbody = displacement of each atom from COM in body frame
------------------------------------------------------------------------- */

void Molecule::compute_inertia()
{
  if (!inertiaflag) {
    inertiaflag = 1;

    atom->check_mass(FLERR);

    double onemass, dx, dy, dz;
    for (int i = 0; i < 6; i++) itensor[i] = 0.0;
    for (int i = 0; i < natoms; i++) {
      if (rmassflag)
        onemass = rmass[i];
      else
        onemass = atom->mass[type[i]];
      dx = dxcom[i][0];
      dy = dxcom[i][1];
      dz = dxcom[i][2];
      itensor[0] += onemass * (dy * dy + dz * dz);
      itensor[1] += onemass * (dx * dx + dz * dz);
      itensor[2] += onemass * (dx * dx + dy * dy);
      itensor[3] -= onemass * dy * dz;
      itensor[4] -= onemass * dx * dz;
      itensor[5] -= onemass * dx * dy;
    }

    if (radiusflag && !bodyflag) {
      for (int i = 0; i < natoms; i++) {
        if (rmassflag)
          onemass = rmass[i];
        else
          onemass = atom->mass[type[i]];
        itensor[0] += SINERTIA * onemass * radius[i] * radius[i];
        itensor[1] += SINERTIA * onemass * radius[i] * radius[i];
        itensor[2] += SINERTIA * onemass * radius[i] * radius[i];
      }
    }
  }

  // diagonalize inertia tensor for each body via Jacobi rotations
  // inertia = 3 eigenvalues = principal moments of inertia
  // evectors and exzy = 3 evectors = principal axes of rigid body

  double cross[3];
  double tensor[3][3], evectors[3][3];

  tensor[0][0] = itensor[0];
  tensor[1][1] = itensor[1];
  tensor[2][2] = itensor[2];
  tensor[1][2] = tensor[2][1] = itensor[3];
  tensor[0][2] = tensor[2][0] = itensor[4];
  tensor[0][1] = tensor[1][0] = itensor[5];

  if (MathEigen::jacobi3(tensor, inertia, evectors))
    error->all(FLERR, Error::NOLASTLINE, "Insufficient Jacobi rotations for rigid molecule");

  ex[0] = evectors[0][0];
  ex[1] = evectors[1][0];
  ex[2] = evectors[2][0];
  ey[0] = evectors[0][1];
  ey[1] = evectors[1][1];
  ey[2] = evectors[2][1];
  ez[0] = evectors[0][2];
  ez[1] = evectors[1][2];
  ez[2] = evectors[2][2];

  // if any principal moment < scaled EPSILON, set to 0.0

  double max;
  max = MAX(inertia[0], inertia[1]);
  max = MAX(max, inertia[2]);

  if (inertia[0] < EPSILON * max) inertia[0] = 0.0;
  if (inertia[1] < EPSILON * max) inertia[1] = 0.0;
  if (inertia[2] < EPSILON * max) inertia[2] = 0.0;

  // enforce 3 evectors as a right-handed coordinate system
  // flip 3rd vector if needed

  MathExtra::cross3(ex, ey, cross);
  if (MathExtra::dot3(cross, ez) < 0.0) MathExtra::negate3(ez);

  // create quaternion

  MathExtra::exyz_to_q(ex, ey, ez, quat);

  // compute displacements in body frame defined by quat

  memory->destroy(dxbody);
  memory->create(dxbody, natoms, 3, "molecule:dxbody");

  for (int i = 0; i < natoms; i++) MathExtra::transpose_matvec(ex, ey, ez, dxcom[i], dxbody[i]);
}

// clang-format off

/* ----------------------------------------------------------------------
   read molecule info from file
   flag = 0, just scan for sizes of fields
   flag = 1, read and store fields
------------------------------------------------------------------------- */

void Molecule::read(int flag)
{
  char line[MAXLINE] = {'\0'};
  char *eof;

  // skip 1st line of file

  if (comm->me == 0) {
    eof = fgets(line, MAXLINE, fp);

    // check for units keyword in first line and print warning on mismatch

    auto units = Tokenizer(utils::strfind(line, "units = \\w+")).as_vector();
    if ((flag == 0) && (units.size() > 2)) {
      if (units[2] != update->unit_style)
        error->warning(FLERR, "Inconsistent units in data file: current = {}, data file = {}",
                       update->unit_style, units[2]);
    }

    if (eof == nullptr) error->one(FLERR, fileiarg, "Unexpected end of molecule file");
  }

  if (flag == 0) title = utils::trim(line);

  // read header lines
  // skip blank lines or lines that start with "#"
  // stop when read an unrecognized line
  bool has_atoms = false;

  while (true) {

    readline(line);

    // trim comments. if line is blank or comment, continue

    auto text = utils::trim(utils::trim_comment(line));
    if (text.empty()) continue;
    if (utils::strmatch(text, "^\\s*#")) continue;

    // search line for header keywords and set corresponding variable
    try {
      ValueTokenizer values(text);

      int nmatch = values.count();
      int nwant = 0;
      if (values.matches("^\\s*\\d+\\s+atoms\\s*$")) {
        natoms = values.next_int();
        nwant = 2;
        has_atoms = true;
      } else if (values.matches("^\\s*\\d+\\s+bonds\\s*$")) {
        nbonds = values.next_int();
        nwant = 2;
      } else if (values.matches("^\\s*\\d+\\s+angles\\s*$")) {
        nangles = values.next_int();
        nwant = 2;
      } else if (values.matches("^\\s*\\d+\\s+dihedrals\\s*$")) {
        ndihedrals = values.next_int();
        nwant = 2;
      } else if (values.matches("^\\s*\\d+\\s+impropers\\s*$")) {
        nimpropers = values.next_int();
        nwant = 2;
      } else if (values.matches("^\\s*\\d+\\s+fragments\\s*$")) {
        nfragments = values.next_int();
        nwant = 2;
      } else if (values.matches("^\\s*\\f+\\s+mass\\s*$")) {
        massflag = massflag_user = 1;
        masstotal = values.next_double();
        nwant = 2;
        masstotal *= sizescale * sizescale * sizescale;
      } else if (values.matches("^\\s*\\f+\\s+\\f+\\s+\\f+\\s+com\\s*$")) {
        comflag = comflag_user = 1;
        com[0] = values.next_double();
        com[1] = values.next_double();
        com[2] = values.next_double();
        nwant = 4;
        com[0] *= sizescale;
        com[1] *= sizescale;
        com[2] *= sizescale;
        if ((domain->dimension == 2) && (com[2] != 0.0))
          error->all(FLERR, fileiarg, "Molecule file z center-of-mass must be 0.0 for 2d systems");
      } else if (values.matches("^\\s*\\f+\\s+\\f+\\s+\\f+\\s+\\f+\\s+\\f+\\s+\\f+\\s+inertia\\s*$")) {
        inertiaflag = inertiaflag_user = 1;
        itensor[0] = values.next_double();
        itensor[1] = values.next_double();
        itensor[2] = values.next_double();
        itensor[3] = values.next_double();
        itensor[4] = values.next_double();
        itensor[5] = values.next_double();
        nwant = 7;
        const double scale5 = powint(sizescale, 5);
        itensor[0] *= scale5;
        itensor[1] *= scale5;
        itensor[2] *= scale5;
        itensor[3] *= scale5;
        itensor[4] *= scale5;
        itensor[5] *= scale5;
      } else if (values.matches("^\\s*\\d+\\s+\\d+\\s+body\\s*$")) {
        bodyflag = 1;
        avec_body = dynamic_cast<AtomVecBody *>(atom->style_match("body"));
        if (!avec_body) error->all(FLERR, fileiarg, "Molecule file requires atom style body");
        nibody = values.next_int();
        ndbody = values.next_int();
        nwant = 3;
      } else if (values.matches("^\\s*\\d+\\s+\\S+\\s+types\\s*$")) {
        error->all(FLERR, fileiarg, "Found data file header keyword '{}' in molecule file", text);
      } else if (values.matches("^\\s*\\f+\\s+\\f+\\s+[xyz]lo\\s+[xyz]hi\\s*$")) {
        error->all(FLERR, fileiarg, "Found data file header keyword '{}' in molecule file", text);
      } else {
        // unknown header keyword
        if (values.matches("^\\s*\\f+\\s+\\S+")) {
          error->all(FLERR, fileiarg, "Unknown keyword or incorrectly formatted header line: {}",
                     line);
        } else
          break;
      }
      if (nmatch != nwant)
        error->all(FLERR, fileiarg, "Invalid header line format in molecule file: {}", line);
    } catch (TokenizerException &e) {
      error->all(FLERR, fileiarg, "Invalid header in molecule file: {}", e.what());
    }
  }

  // error checks

  if (!has_atoms)
    error->all(FLERR, fileiarg, "Required \"atoms\" header keyword not found in molecule file");
  if (natoms < 1) error->all(FLERR, fileiarg, "No atoms or invalid atom count in molecule file");
  if (nbonds < 0) error->all(FLERR, fileiarg, "Invalid bond count in molecule file");
  if (nangles < 0) error->all(FLERR, fileiarg, "Invalid angle count in molecule file");
  if (ndihedrals < 0) error->all(FLERR, fileiarg, "Invalid dihedral count in molecule file");
  if (nimpropers < 0) error->all(FLERR, fileiarg, "Invalid improper count in molecule file");

  // count = vector for tallying bonds,angles,etc per atom

  if (flag == 0) memory->create(count, natoms, "molecule:count");

  // grab keyword and skip next line

  std::string keyword = parse_keyword(0, line);
  readline(line);

  // loop over sections of molecule file

  while (!keyword.empty()) {
    if (keyword == "Coords") {
      xflag = 1;
      if (flag)
        coords(line);
      else
        skip_lines(natoms, line, keyword);
    } else if (keyword == "Types") {
      typeflag = 1;
      if (flag)
        types(line);
      else
        skip_lines(natoms, line, keyword);
    } else if (keyword == "Molecules") {
      moleculeflag = 1;
      if (flag)
        molecules(line);
      else
        skip_lines(natoms, line, keyword);
    } else if (keyword == "Fragments") {
      if (nfragments == 0)
        error->all(FLERR, fileiarg, "Found Fragments section but no nfragments setting in header");
      fragmentflag = 1;
      if (flag)
        fragments(line);
      else
        skip_lines(nfragments, line, keyword);
    } else if (keyword == "Charges") {
      qflag = 1;
      if (flag)
        charges(line);
      else
        skip_lines(natoms, line, keyword);
    } else if (keyword == "Diameters") {
      radiusflag = 1;
      if (flag)
        diameters(line);
      else
        skip_lines(natoms, line, keyword);
    } else if (keyword == "Dipoles") {
      muflag = 1;
      if (flag)
        dipoles(line);
      else
        skip_lines(natoms, line, keyword);
    } else if (keyword == "Masses") {
      rmassflag = 1;
      if (flag)
        masses(line);
      else
        skip_lines(natoms, line, keyword);

    } else if (keyword == "Bonds") {
      if (nbonds == 0)
        error->all(FLERR, fileiarg, "Found Bonds section but no nbonds setting in header");
      bondflag = tag_require = 1;
      bonds(flag, line);
    } else if (keyword == "Angles") {
      if (nangles == 0)
        error->all(FLERR, fileiarg, "Found Angles section but no nangles setting in header");
      angleflag = tag_require = 1;
      angles(flag, line);
    } else if (keyword == "Dihedrals") {
      if (ndihedrals == 0)
        error->all(FLERR, fileiarg, "Found Dihedrals section but no ndihedrals setting in header");
      dihedralflag = tag_require = 1;
      dihedrals(flag, line);
    } else if (keyword == "Impropers") {
      if (nimpropers == 0)
        error->all(FLERR, fileiarg, "Found Impropers section but no nimpropers setting in header");
      improperflag = tag_require = 1;
      impropers(flag, line);

    } else if (keyword == "Special Bond Counts") {
      nspecialflag = specialflag_user = 1;
      nspecial_read(flag, line);
    } else if (keyword == "Special Bonds") {
      specialflag = tag_require = 1;
      if (!nspecialflag)
        error->all(FLERR, fileiarg,
                   "Special Bond Counts section must come before Special Bonds section");
      if (flag)
        special_read(line);
      else
        skip_lines(natoms, line, keyword);

    } else if (keyword == "Shake Flags") {
      shakeflagflag = 1;
      if (flag)
        shakeflag_read(line);
      else
        skip_lines(natoms, line, keyword);
    } else if (keyword == "Shake Atoms") {
      shakeatomflag = tag_require = 1;
      if (shaketypeflag) shakeflag = 1;
      if (!shakeflagflag)
        error->all(FLERR, fileiarg, "Shake Flags section must come before Shake Atoms section");
      if (flag)
        shakeatom_read(line);
      else
        skip_lines(natoms, line, keyword);
    } else if (keyword == "Shake Bond Types") {
      shaketypeflag = 1;
      if (shakeatomflag) shakeflag = 1;
      if (!shakeflagflag)
        error->all(FLERR, fileiarg, "Shake Flags section must come before Shake Bonds section");
      if (flag)
        shaketype_read(line);
      else
        skip_lines(natoms, line, keyword);

    } else if (keyword == "Body Integers") {
      if (bodyflag == 0 || nibody == 0)
        error->all(FLERR, fileiarg, "Found Body Integers section but no setting in header");
      ibodyflag = 1;
      body(flag, 0, line);
    } else if (keyword == "Body Doubles") {
      if (bodyflag == 0 || ndbody == 0)
        error->all(FLERR, fileiarg, "Found Body Doubles section but no setting in header");
      dbodyflag = 1;
      body(flag, 1, line);
    } else if ((keyword == "Atoms") || (keyword == "Velocities") || (keyword == "Pair Coeffs") ||
               (keyword == "Bond Coeffs") || (keyword == "Angle Coeffs") ||
               (keyword == "Dihedral Coeffs") || (keyword == "Improper Coeffs")) {
      error->all(FLERR, fileiarg, "Found data file section '{}' in molecule file\n", keyword);
    } else {

      // Error: Either a too long/short section or a typo in the keyword

      if (utils::strmatch(keyword, "^[A-Za-z ]+$"))
        error->all(FLERR, fileiarg, "Unknown section '{}' in molecule file\n", keyword);
      else
        error->all(FLERR, fileiarg,
                   "Unexpected line in molecule file while looking for the next section:\n{}",
                   line);
    }
    keyword = parse_keyword(1, line);
  }

  // error check

  if (flag == 0) {
    if ((nspecialflag && !specialflag) || (!nspecialflag && specialflag))
      error->all(FLERR, fileiarg, "Molecule file needs both Special Bond sections");
    if (specialflag && !bondflag)
      error->all(FLERR, fileiarg, "Molecule file has special flags but no bonds");
    if ((shakeflagflag || shakeatomflag || shaketypeflag) && !shakeflag)
      error->all(FLERR, fileiarg, "Molecule file shake info is incomplete");
    if (bodyflag && !rmassflag)
      error->all(FLERR, fileiarg, "Molecule file must have Masses section for body particle");
    if (bodyflag && nibody && ibodyflag == 0)
      error->all(FLERR, fileiarg, "Molecule file has no Body Integers section");
    if (bodyflag && ndbody && dbodyflag == 0)
      error->all(FLERR, fileiarg, "Molecule file has no Body Doubles section");
    if (nfragments > 0 && !fragmentflag)
      error->all(FLERR, fileiarg, "Molecule file has no Fragments section");
  }

  // auto-generate special bonds if needed and not in file

  if (bondflag && specialflag == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR, fileiarg,
                 "Cannot auto-generate special bonds before simulation box is defined");

    if (flag) {
      special_generate();
      specialflag_user = 0;
      specialflag = 1;
      nspecialflag = 1;
    }
  }

  // body particle must have natom = 1
  // set radius by having body class compute its own radius

  if (bodyflag) {
    radiusflag = 1;
    if (natoms != 1) error->all(FLERR, fileiarg, "Molecule natoms must be 1 for body particle");
    if (sizescale != 1.0)
      error->all(FLERR, fileiarg, "Molecule sizescale must be 1.0 for body particle");
    if (flag) {
      radius[0] = avec_body->radius_body(nibody, ndbody, ibodyparams, dbodyparams);
      maxradius = radius[0];
    }
  }

  // clean up

  if (flag) memory->destroy(count);
}

/* ----------------------------------------------------------------------
   read coords from file
------------------------------------------------------------------------- */

void Molecule::coords(char *line)
{
  for (int i = 0; i < natoms; i++) count[i] = 0;
  try {
    for (int i = 0; i < natoms; i++) {
      readline(line);

      ValueTokenizer values(utils::trim_comment(line));
      if (values.count() != 4)
        error->all(FLERR, fileiarg, "Invalid line in Coords section of molecule file: {}", line);

      int iatom = values.next_int() - 1;
      if (iatom < 0 || iatom >= natoms)
        error->all(FLERR, fileiarg, "Invalid atom index {} in Coords section of molecule file", iatom);
      count[iatom]++;
      x[iatom][0] = values.next_double();
      x[iatom][1] = values.next_double();
      x[iatom][2] = values.next_double();

      x[iatom][0] *= sizescale;
      x[iatom][1] *= sizescale;
      x[iatom][2] *= sizescale;
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, fileiarg, "Invalid line in Coords section of molecule file: {}\n{}", e.what(),
               line);
  }

  for (int i = 0; i < natoms; i++)
    if (count[i] == 0)
      error->all(FLERR, fileiarg, "Atom {} missing in Coords section of molecule file", i + 1);

  if (domain->dimension == 2) {
    for (int i = 0; i < natoms; i++)
      if (x[i][2] != 0.0)
        error->all(FLERR, fileiarg,
                   "Z coord in molecule file for atom {} must be 0.0 for 2d-simulation", i + 1);
  }
}

/* ----------------------------------------------------------------------
   read types from file
   set ntypes = max of any atom type
------------------------------------------------------------------------- */

// clang-format on
void Molecule::types(char *line)
{
  const std::string location = "Types section of molecule file";
  std::string typestr;
  for (int i = 0; i < natoms; i++) count[i] = 0;

  for (int i = 0; i < natoms; i++) {
    readline(line);
    auto values = Tokenizer(utils::trim(line)).as_vector();
    int nwords = values.size();
    for (std::size_t ii = 0; ii < values.size(); ++ii) {
      if (utils::strmatch(values[ii], "^#")) {
        nwords = ii;
        break;
      }
    }
    if (nwords != 2)
      error->all(FLERR, fileiarg, "Invalid format in {}: {}", location, utils::trim(line));

    int iatom = utils::inumeric(FLERR, values[0], false, lmp);
    if (iatom < 1 || iatom > natoms)
      error->all(FLERR, fileiarg, "Invalid atom index {} in {}: {}", iatom, location,
                 utils::trim(line));
    count[--iatom]++;

    typestr = utils::utf8_subst(values[1]);
    switch (utils::is_type(typestr)) {
      case 0: {    // numeric
        type[iatom] = utils::inumeric(FLERR, typestr, false, lmp);
        type[iatom] += toffset;
        break;
      }
      case 1: {    // type label
        if (!atom->labelmapflag)
          error->all(FLERR, fileiarg, "Invalid atom type {} in {}: {}", typestr, location,
                     utils::trim(line));
        type[iatom] = atom->lmap->find(typestr, Atom::ATOM);
        if (type[iatom] == -1)
          error->all(FLERR, fileiarg, "Unknown atom type {} in {}: {}", typestr, location,
                     utils::trim(line));
        break;
      }
      default:    // invalid
        error->one(FLERR, fileiarg, "Invalid format in {}: {}", location, utils::trim(line));
        break;
    }
  }

  for (int i = 0; i < natoms; i++) {
    if (count[i] == 0) error->all(FLERR, fileiarg, "Atom {} missing in {}", i + 1, location);
    if ((type[i] <= 0) || (domain->box_exist && (type[i] > atom->ntypes)))
      error->all(FLERR, fileiarg, "Invalid atom type {} for atom {} in molecule file", type[i],
                 i + 1);
    ntypes = MAX(ntypes, type[i]);
  }
}
// clang-format off
/* ----------------------------------------------------------------------
   read molecules from file
   set nmolecules = max of any molecule type
------------------------------------------------------------------------- */

void Molecule::molecules(char *line)
{
  for (int i = 0; i < natoms; i++) count[i] = 0;
  try {
    for (int i = 0; i < natoms; i++) {
      readline(line);
      ValueTokenizer values(utils::trim_comment(line));
      if (values.count() != 2)
        error->all(FLERR, fileiarg, "Invalid line in Molecules section of molecule file: {}", line);

      int iatom = values.next_int() - 1;
      if (iatom < 0 || iatom >= natoms)
        error->all(FLERR, fileiarg, "Invalid atom index in Molecules section of molecule file");
      count[iatom]++;
      molecule[iatom] = values.next_tagint();
      // molecule[iatom] += moffset; // placeholder for possible molecule offset
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, fileiarg, "Invalid line in Molecules section of molecule file: {}\n{}", e.what(), line);
  }

  for (int i = 0; i < natoms; i++) {
    if (count[i] == 0)
      error->all(FLERR, fileiarg, "Atom {} missing in Molecules section of molecule file", i + 1);
  }
  for (int i = 0; i < natoms; i++) {
    if (molecule[i] < 0)
      error->all(FLERR, fileiarg, "Invalid molecule ID {} for atom {} in molecule file", molecule[i], i + 1);
  }
  for (int i = 0; i < natoms; i++) nmolecules = MAX(nmolecules, molecule[i]);
}

/* ----------------------------------------------------------------------
   read fragments from file
------------------------------------------------------------------------- */

void Molecule::fragments(char *line)
{
  try {
    for (int i = 0; i < nfragments; i++) {
      readline(line);

      ValueTokenizer values(utils::trim_comment(line));

      if ((int) values.count() > natoms + 1)
        error->all(FLERR, fileiarg, "Too many atoms per fragment in Fragments section of molecule file");

      fragmentnames[i] = values.next_string();

      while (values.has_next()) {
        int iatom = values.next_int() - 1;
        if (iatom < 0 || iatom >= natoms)
          error->all(FLERR, fileiarg,
                     "Invalid atom ID {} for fragment {} in Fragments section of molecule file",
                     iatom + 1, fragmentnames[i]);
        fragmentmask[i][iatom] = 1;
      }
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, fileiarg,
               "Invalid atom ID in Fragments section of "
               "molecule file: {}\n{}",
               e.what(), line);
  }
}

/* ----------------------------------------------------------------------
   read charges from file
------------------------------------------------------------------------- */

void Molecule::charges(char *line)
{
  for (int i = 0; i < natoms; i++) count[i] = 0;
  try {
    for (int i = 0; i < natoms; i++) {
      readline(line);

      ValueTokenizer values(utils::trim_comment(line));
      if ((int) values.count() != 2)
        error->all(FLERR, fileiarg, "Invalid line in Charges section of molecule file: {}", line);

      int iatom = values.next_int() - 1;
      if (iatom < 0 || iatom >= natoms)
        error->all(FLERR, fileiarg, "Invalid atom index in Charges section of molecule file");

      count[iatom]++;
      q[iatom] = values.next_double();
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, fileiarg, "Invalid line in Charges section of molecule file: {}\n{}", e.what(), line);
  }

  for (int i = 0; i < natoms; i++) {
    if (count[i] == 0)
      error->all(FLERR, fileiarg, "Atom {} missing in Charges section of molecule file", i + 1);
  }
}

/* ----------------------------------------------------------------------
   read diameters from file and set radii
------------------------------------------------------------------------- */

void Molecule::diameters(char *line)
{
  for (int i = 0; i < natoms; i++) count[i] = 0;
  try {
    maxradius = 0.0;
    for (int i = 0; i < natoms; i++) {
      readline(line);

      ValueTokenizer values(utils::trim_comment(line));
      if (values.count() != 2)
        error->all(FLERR, fileiarg, "Invalid line in Diameters section of molecule file: {}", line);
      int iatom = values.next_int() - 1;
      if (iatom < 0 || iatom >= natoms)
        error->all(FLERR, fileiarg, "Invalid atom index in Diameters section of molecule file");
      count[iatom]++;
      radius[iatom] = values.next_double();
      radius[iatom] *= sizescale;
      radius[iatom] *= 0.5;
      maxradius = MAX(maxradius, radius[iatom]);
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, fileiarg, "Invalid line in Diameters section of molecule file: {}\n{}", e.what(), line);
  }

  for (int i = 0; i < natoms; i++) {
    if (count[i] == 0)
      error->all(FLERR, fileiarg, "Atom {} missing in Diameters section of molecule file", i + 1);
    if (radius[i] < 0.0)
      error->all(FLERR, fileiarg, "Invalid atom diameter {} for atom {} in molecule file", radius[i] * 2.0 / sizescale, i + 1);
  }
}

/* ----------------------------------------------------------------------
   read dipoles from file
------------------------------------------------------------------------- */

void Molecule::dipoles(char *line)
{
  for (int i = 0; i < natoms; i++) count[i] = 0;
  try {
    for (int i = 0; i < natoms; i++) {
      readline(line);

      ValueTokenizer values(utils::trim_comment(line));
      if ((int) values.count() != 4)
        error->all(FLERR, fileiarg, "Invalid line in Dipoles section of molecule file: {}", line);

      int iatom = values.next_int() - 1;
      if (iatom < 0 || iatom >= natoms)
        error->all(FLERR, fileiarg, "Invalid atom index in Dipoles section of molecule file");

      count[iatom]++;
      mu[iatom][0] = values.next_double() * sizescale;
      mu[iatom][1] = values.next_double() * sizescale;
      mu[iatom][2] = values.next_double() * sizescale;
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, fileiarg, "Invalid line in Dipoles section of molecule file: {}\n{}", e.what(), line);
  }

  for (int i = 0; i < natoms; i++) {
    if (count[i] == 0)
      error->all(FLERR, fileiarg, "Atom {} missing in Dipoles section of molecule file", i + 1);
  }
  if (domain->dimension == 2) {
    for (int i = 0; i < natoms; i++)
      if (mu[i][2] != 0.0)
        error->all(FLERR, fileiarg, "Dipole moment z-component in JSON data for atom {} "
                   "must be 0.0 for 2d-simulation", id, i + 1);
  }
}

/* ----------------------------------------------------------------------
   read masses from file
------------------------------------------------------------------------- */

void Molecule::masses(char *line)
{
  for (int i = 0; i < natoms; i++) count[i] = 0;
  try {
    for (int i = 0; i < natoms; i++) {
      readline(line);

      ValueTokenizer values(utils::trim_comment(line));
      if (values.count() != 2)
        error->all(FLERR, fileiarg, "Invalid line in Masses section of molecule file: {}", line);

      int iatom = values.next_int() - 1;
      if (iatom < 0 || iatom >= natoms)
        error->all(FLERR, fileiarg, "Invalid atom index in Masses section of molecule file");
      count[iatom]++;
      rmass[iatom] = values.next_double();
      rmass[iatom] *= sizescale * sizescale * sizescale;
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, fileiarg, "Invalid line in Masses section of molecule file: {}\n{}", e.what(), line);
  }

  for (int i = 0; i < natoms; i++) {
    if (count[i] == 0)
      error->all(FLERR, fileiarg, "Atom {} missing in Masses section of molecule file", i + 1);
    if (rmass[i] <= 0.0)
      error->all(FLERR, fileiarg, "Invalid atom mass {} for atom {} in molecule file", rmass[i] / sizescale / sizescale
    / sizescale, i + 1);
  }
}

/* ----------------------------------------------------------------------
   read bonds from file
   set nbondtypes = max type of any bond
   store each with both atoms if newton_bond = 0
   if flag = 0, just count bonds/atom
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void Molecule::bonds(int flag, char *line)
{
  const std::string location = "Bonds section of molecule file";
  int itype;
  tagint m, atom1, atom2;
  std::string typestr;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_bond[i] = 0;

  for (int i = 0; i < nbonds; i++) {
    readline(line);
    auto values = Tokenizer(utils::trim(line)).as_vector();
    int nwords = values.size();
    for (std::size_t ii = 0; ii < values.size(); ++ii) {
      if (utils::strmatch(values[ii], "^#")) {
        nwords = ii;
        break;
      }
    }
    if (nwords != 4) error->all(FLERR, fileiarg, "Invalid format in {}: {}", location, utils::trim(line));

    typestr = utils::utf8_subst(values[1]);
    switch (utils::is_type(typestr)) {
      case 0: {    // numeric
        itype = utils::inumeric(FLERR, typestr, false, lmp);
        itype += boffset;
        break;
      }
      case 1: {    // type label
        if (!atom->labelmapflag)
          error->all(FLERR, fileiarg, "Invalid bond type {} in {}: {}", typestr, location, utils::trim(line));
        itype = atom->lmap->find(typestr, Atom::BOND);
        if (itype == -1)
          error->all(FLERR, fileiarg, "Unknown bond type {} in {}: {}", typestr, location, utils::trim(line));
        break;
      }
      default:    // invalid
        error->one(FLERR, fileiarg, "Invalid format in {}: {}", location, utils::trim(line));
        break;
    }

    atom1 = utils::tnumeric(FLERR, values[2], false, lmp);
    atom2 = utils::tnumeric(FLERR, values[3], false, lmp);

    if ((atom1 <= 0) || (atom1 > natoms) || (atom2 <= 0) || (atom2 > natoms) || (atom1 == atom2))
      error->all(FLERR, fileiarg, "Invalid atom ID in {}: {}", location, utils::trim(line));
    if ((itype <= 0) || (domain->box_exist && (itype > atom->nbondtypes)))
      error->all(FLERR, fileiarg, "Invalid bond type in {}: {}", location, utils::trim(line));

    if (flag) {
      m = atom1 - 1;
      nbondtypes = MAX(nbondtypes, itype);
      bond_type[m][num_bond[m]] = itype;
      bond_atom[m][num_bond[m]] = atom2;
      num_bond[m]++;
      if (newton_bond == 0) {
        m = atom2 - 1;
        bond_type[m][num_bond[m]] = itype;
        bond_atom[m][num_bond[m]] = atom1;
        num_bond[m]++;
      }
    } else {
      count[atom1 - 1]++;
      if (newton_bond == 0) count[atom2 - 1]++;
    }
  }

  // bond_per_atom = max of count vector

  if (flag == 0) {
    bond_per_atom = 0;
    for (int i = 0; i < natoms; i++) bond_per_atom = MAX(bond_per_atom, count[i]);
  }
}

/* ----------------------------------------------------------------------
   read angles from file
   store each with all 3 atoms if newton_bond = 0
   if flag = 0, just count angles/atom
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void Molecule::angles(int flag, char *line)
{
  const std::string location = "Angles section of molecule file";
  int itype;
  tagint m, atom1, atom2, atom3;
  std::string typestr;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_angle[i] = 0;

  for (int i = 0; i < nangles; i++) {
    readline(line);
    auto values = Tokenizer(utils::trim(line)).as_vector();
    int nwords = values.size();
    for (std::size_t ii = 0; ii < values.size(); ++ii) {
      if (utils::strmatch(values[ii], "^#")) {
        nwords = ii;
        break;
      }
    }
    if (nwords != 5) error->all(FLERR, fileiarg, "Invalid format in {}: {}", location, utils::trim(line));

    typestr = utils::utf8_subst(values[1]);
    switch (utils::is_type(typestr)) {
      case 0: {    // numeric
        itype = utils::inumeric(FLERR, typestr, false, lmp);
        itype += aoffset;
        break;
      }
      case 1: {    // type label
        if (!atom->labelmapflag)
          error->all(FLERR, fileiarg, "Invalid angle type {} in {}: {}", typestr, location, utils::trim(line));
        itype = atom->lmap->find(typestr, Atom::ANGLE);
        if (itype == -1)
          error->all(FLERR, fileiarg, "Unknown angle type {} in {}: {}", typestr, location, utils::trim(line));
        break;
      }
      default:    // invalid
        error->one(FLERR, fileiarg, "Invalid format in {}: {}", location, utils::trim(line));
        break;
    }

    atom1 = utils::tnumeric(FLERR, values[2], false, lmp);
    atom2 = utils::tnumeric(FLERR, values[3], false, lmp);
    atom3 = utils::tnumeric(FLERR, values[4], false, lmp);

    if ((atom1 <= 0) || (atom1 > natoms) || (atom2 <= 0) || (atom2 > natoms) || (atom3 <= 0) ||
        (atom3 > natoms) || (atom1 == atom2) || (atom1 == atom3) || (atom2 == atom3))
      error->all(FLERR, fileiarg, "Invalid atom ID in {}: {}", location, utils::trim(line));
    if ((itype <= 0) || (domain->box_exist && (itype > atom->nangletypes)))
      error->all(FLERR, fileiarg, "Invalid angle type in {}: {}", location, utils::trim(line));

    if (flag) {
      m = atom2 - 1;
      nangletypes = MAX(nangletypes, itype);
      angle_type[m][num_angle[m]] = itype;
      angle_atom1[m][num_angle[m]] = atom1;
      angle_atom2[m][num_angle[m]] = atom2;
      angle_atom3[m][num_angle[m]] = atom3;
      num_angle[m]++;
      if (newton_bond == 0) {
        m = atom1 - 1;
        angle_type[m][num_angle[m]] = itype;
        angle_atom1[m][num_angle[m]] = atom1;
        angle_atom2[m][num_angle[m]] = atom2;
        angle_atom3[m][num_angle[m]] = atom3;
        num_angle[m]++;
        m = atom3 - 1;
        angle_type[m][num_angle[m]] = itype;
        angle_atom1[m][num_angle[m]] = atom1;
        angle_atom2[m][num_angle[m]] = atom2;
        angle_atom3[m][num_angle[m]] = atom3;
        num_angle[m]++;
      }
    } else {
      count[atom2 - 1]++;
      if (newton_bond == 0) {
        count[atom1 - 1]++;
        count[atom3 - 1]++;
      }
    }
  }

  // angle_per_atom = max of count vector

  if (flag == 0) {
    angle_per_atom = 0;
    for (int i = 0; i < natoms; i++) angle_per_atom = MAX(angle_per_atom, count[i]);
  }
}

/* ----------------------------------------------------------------------
   read dihedrals from file
   store each with all 4 atoms if newton_bond = 0
   if flag = 0, just count dihedrals/atom
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void Molecule::dihedrals(int flag, char *line)
{
  const std::string location = "Dihedrals section of molecule file";
  int itype;
  tagint m, atom1, atom2, atom3, atom4;
  std::string typestr;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_dihedral[i] = 0;

  for (int i = 0; i < ndihedrals; i++) {
    readline(line);
    auto values = Tokenizer(utils::trim(line)).as_vector();
    int nwords = values.size();
    for (std::size_t ii = 0; ii < values.size(); ++ii) {
      if (utils::strmatch(values[ii], "^#")) {
        nwords = ii;
        break;
      }
    }
    if (nwords != 6) error->all(FLERR, fileiarg, "Invalid format in {}: {}", location, utils::trim(line));

    typestr = utils::utf8_subst(values[1]);
    switch (utils::is_type(typestr)) {
      case 0: {    // numeric
        itype = utils::inumeric(FLERR, typestr, false, lmp);
        itype += doffset;
        break;
      }
      case 1: {    // type label
        if (!atom->labelmapflag)
          error->all(FLERR, fileiarg, "Invalid dihedral type {} in {}: {}", typestr, location, utils::trim(line));
        itype = atom->lmap->find(typestr, Atom::DIHEDRAL);
        if (itype == -1)
          error->all(FLERR, fileiarg, "Unknown dihedral type {} in {}: {}", typestr, location, utils::trim(line));
        break;
      }
      default:    // invalid
        error->one(FLERR, fileiarg, "Invalid format in {}: {}", location, utils::trim(line));
        break;
    }

    atom1 = utils::tnumeric(FLERR, values[2], false, lmp);
    atom2 = utils::tnumeric(FLERR, values[3], false, lmp);
    atom3 = utils::tnumeric(FLERR, values[4], false, lmp);
    atom4 = utils::tnumeric(FLERR, values[5], false, lmp);

    if ((atom1 <= 0) || (atom1 > natoms) || (atom2 <= 0) || (atom2 > natoms) || (atom3 <= 0) ||
        (atom3 > natoms) || (atom4 <= 0) || (atom4 > natoms) || (atom1 == atom2) ||
        (atom1 == atom3) || (atom1 == atom4) || (atom2 == atom3) || (atom2 == atom4) ||
        (atom3 == atom4))
      error->all(FLERR, fileiarg, "Invalid atom ID in {}: {}", location, utils::trim(line));
    if ((itype <= 0) || (domain->box_exist && (itype > atom->ndihedraltypes)))
      error->all(FLERR, fileiarg, "Invalid dihedral type in {}: {}", location, utils::trim(line));

    if (flag) {
      m = atom2 - 1;
      ndihedraltypes = MAX(ndihedraltypes, itype);
      dihedral_type[m][num_dihedral[m]] = itype;
      dihedral_atom1[m][num_dihedral[m]] = atom1;
      dihedral_atom2[m][num_dihedral[m]] = atom2;
      dihedral_atom3[m][num_dihedral[m]] = atom3;
      dihedral_atom4[m][num_dihedral[m]] = atom4;
      num_dihedral[m]++;
      if (newton_bond == 0) {
        m = atom1 - 1;
        dihedral_type[m][num_dihedral[m]] = itype;
        dihedral_atom1[m][num_dihedral[m]] = atom1;
        dihedral_atom2[m][num_dihedral[m]] = atom2;
        dihedral_atom3[m][num_dihedral[m]] = atom3;
        dihedral_atom4[m][num_dihedral[m]] = atom4;
        num_dihedral[m]++;
        m = atom3 - 1;
        dihedral_type[m][num_dihedral[m]] = itype;
        dihedral_atom1[m][num_dihedral[m]] = atom1;
        dihedral_atom2[m][num_dihedral[m]] = atom2;
        dihedral_atom3[m][num_dihedral[m]] = atom3;
        dihedral_atom4[m][num_dihedral[m]] = atom4;
        num_dihedral[m]++;
        m = atom4 - 1;
        dihedral_type[m][num_dihedral[m]] = itype;
        dihedral_atom1[m][num_dihedral[m]] = atom1;
        dihedral_atom2[m][num_dihedral[m]] = atom2;
        dihedral_atom3[m][num_dihedral[m]] = atom3;
        dihedral_atom4[m][num_dihedral[m]] = atom4;
        num_dihedral[m]++;
      }
    } else {
      count[atom2 - 1]++;
      if (newton_bond == 0) {
        count[atom1 - 1]++;
        count[atom3 - 1]++;
        count[atom4 - 1]++;
      }
    }
  }

  // dihedral_per_atom = max of count vector

  if (flag == 0) {
    dihedral_per_atom = 0;
    for (int i = 0; i < natoms; i++) dihedral_per_atom = MAX(dihedral_per_atom, count[i]);
  }
}

/* ----------------------------------------------------------------------
   read impropers from file
   store each with all 4 atoms if newton_bond = 0
   if flag = 0, just count impropers/atom
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void Molecule::impropers(int flag, char *line)
{
  const std::string location = "Impropers section of molecule file";
  int itype;
  tagint m, atom1, atom2, atom3, atom4;
  std::string typestr;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_improper[i] = 0;

  for (int i = 0; i < nimpropers; i++) {
    readline(line);
    auto values = Tokenizer(utils::trim(line)).as_vector();
    int nwords = values.size();
    for (std::size_t ii = 0; ii < values.size(); ++ii) {
      if (utils::strmatch(values[ii], "^#")) {
        nwords = ii;
        break;
      }
    }
    if (nwords != 6) error->all(FLERR, fileiarg, "Invalid format in {}: {}", location, utils::trim(line));

    typestr = utils::utf8_subst(values[1]);
    switch (utils::is_type(typestr)) {
      case 0: {    // numeric
        itype = utils::inumeric(FLERR, typestr, false, lmp);
        itype += ioffset;
        break;
      }
      case 1: {    // type label
        if (!atom->labelmapflag)
          error->all(FLERR, fileiarg, "Invalid improper type {} in {}: {}", typestr, location, utils::trim(line));
        itype = atom->lmap->find(typestr, Atom::IMPROPER);
        if (itype == -1)
          error->all(FLERR, fileiarg, "Unknown improper type {} in {}: {}", typestr, location, utils::trim(line));
        break;
      }
      default:    // invalid
        error->one(FLERR, fileiarg, "Invalid format in {}: {}", location, utils::trim(line));
        break;
    }

    atom1 = utils::tnumeric(FLERR, values[2], false, lmp);
    atom2 = utils::tnumeric(FLERR, values[3], false, lmp);
    atom3 = utils::tnumeric(FLERR, values[4], false, lmp);
    atom4 = utils::tnumeric(FLERR, values[5], false, lmp);

    if ((atom1 <= 0) || (atom1 > natoms) || (atom2 <= 0) || (atom2 > natoms) || (atom3 <= 0) ||
        (atom3 > natoms) || (atom4 <= 0) || (atom4 > natoms) || (atom1 == atom2) ||
        (atom1 == atom3) || (atom1 == atom4) || (atom2 == atom3) || (atom2 == atom4) ||
        (atom3 == atom4))
      error->all(FLERR, fileiarg, "Invalid atom ID in {}: {}", location, utils::trim(line));
    if ((itype <= 0) || (domain->box_exist && (itype > atom->nimpropertypes)))
      error->all(FLERR, fileiarg, "Invalid improper type in {}: {}", location, utils::trim(line));

    if (flag) {
      m = atom2 - 1;
      nimpropertypes = MAX(nimpropertypes, itype);
      improper_type[m][num_improper[m]] = itype;
      improper_atom1[m][num_improper[m]] = atom1;
      improper_atom2[m][num_improper[m]] = atom2;
      improper_atom3[m][num_improper[m]] = atom3;
      improper_atom4[m][num_improper[m]] = atom4;
      num_improper[m]++;
      if (newton_bond == 0) {
        m = atom1 - 1;
        improper_type[m][num_improper[m]] = itype;
        improper_atom1[m][num_improper[m]] = atom1;
        improper_atom2[m][num_improper[m]] = atom2;
        improper_atom3[m][num_improper[m]] = atom3;
        improper_atom4[m][num_improper[m]] = atom4;
        num_improper[m]++;
        m = atom3 - 1;
        improper_type[m][num_improper[m]] = itype;
        improper_atom1[m][num_improper[m]] = atom1;
        improper_atom2[m][num_improper[m]] = atom2;
        improper_atom3[m][num_improper[m]] = atom3;
        improper_atom4[m][num_improper[m]] = atom4;
        num_improper[m]++;
        m = atom4 - 1;
        improper_type[m][num_improper[m]] = itype;
        improper_atom1[m][num_improper[m]] = atom1;
        improper_atom2[m][num_improper[m]] = atom2;
        improper_atom3[m][num_improper[m]] = atom3;
        improper_atom4[m][num_improper[m]] = atom4;
        num_improper[m]++;
      }
    } else {
      count[atom2 - 1]++;
      if (newton_bond == 0) {
        count[atom1 - 1]++;
        count[atom3 - 1]++;
        count[atom4 - 1]++;
      }
    }
  }

  // improper_per_atom = max of count vector

  if (flag == 0) {
    improper_per_atom = 0;
    for (int i = 0; i < natoms; i++) improper_per_atom = MAX(improper_per_atom, count[i]);
  }
}

/* ----------------------------------------------------------------------
   read 3 special bonds counts from file
   if flag = 0, just tally maxspecial
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void Molecule::nspecial_read(int flag, char *line)
{
  if (flag == 0) maxspecial = 0;

  for (int i = 0; i < natoms; ++i) count[i] = 0;
  for (int i = 0; i < natoms; ++i) {
    readline(line);

    int c0, c1, c2, c3;

    try {
      ValueTokenizer values(utils::trim_comment(line));
      if (values.count() != 4)
        error->all(FLERR, fileiarg, "Invalid line in Special Bond Counts section of molecule file: {}", line);
      c0 = values.next_int();
      c1 = values.next_int();
      c2 = values.next_int();
      c3 = values.next_int();
    } catch (TokenizerException &e) {
      error->all(FLERR, fileiarg, "Invalid line in Special Bond Counts section of molecule file: {}\n{}",
                 e.what(), line);
    }

    if (flag) {
      int iatom = c0 - 1;
      if (iatom < 0 || iatom >= natoms)
        error->all(FLERR, fileiarg, "Invalid atom index in Special Bond Counts section of molecule file");
      count[iatom]++;
      nspecial[iatom][0] = c1;
      nspecial[iatom][1] = c1 + c2;
      nspecial[iatom][2] = c1 + c2 + c3;
    } else {
      maxspecial = MAX(maxspecial, c1 + c2 + c3);
    }
  }

  // check
  if (flag) {
    for (int i = 0; i < natoms; i++) {
      if (count[i] == 0)
        error->all(FLERR, fileiarg, "Atom {} missing in Special Bond Counts section of molecule file", i + 1);
    }
  }
}


/* ----------------------------------------------------------------------
   read special bond indices from file
------------------------------------------------------------------------- */

void Molecule::special_read(char *line)
{
  for (int i = 0; i < natoms; ++i) count[i] = 0;
  try {
    for (int i = 0; i < natoms; i++) {
      readline(line);

      ValueTokenizer values(utils::trim_comment(line));
      int nwords = values.count();

      if (nwords != nspecial[i][2] + 1)
        error->all(FLERR, fileiarg, "Molecule file special list does not match special count");

      int iatom = values.next_int() - 1;
      if (iatom < 0 || iatom >= natoms)
        error->all(FLERR, fileiarg, "Invalid atom index in Special Bonds section of molecule file");

      for (int m = 1; m < nwords; m++) {
        tagint ival = values.next_tagint();
        if ((ival <= 0) || (ival > natoms) || (ival == iatom + 1))
          error->all(FLERR, fileiarg, "Invalid atom index {} in Special Bonds section of "
                     "molecule file", ival);
        special[iatom][m - 1] = ival;
      }
      count[iatom]++;
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, fileiarg, "Invalid line in Special Bonds section of molecule file: {}\n{}",
               e.what(), line);
  }
  for (int i = 0; i < natoms; i++) {
    if (count[i] == 0)
      error->all(FLERR, fileiarg, "Atom {} missing in Special Bonds section of molecule file",
                 i + 1);
  }
}

/* ----------------------------------------------------------------------
   auto generate special bond info
------------------------------------------------------------------------- */

void Molecule::special_generate()
{
  int newton_bond = force->newton_bond;
  tagint atom1, atom2;

  // temporary array for special atoms

  tagint **tmpspecial;
  memory->create(tmpspecial, natoms, atom->maxspecial, "molecule:tmpspecial");
  memset(&tmpspecial[0][0], 0, sizeof(tagint) * natoms * atom->maxspecial);

  for (int i = 0; i < natoms; i++) count[i] = 0;

  // 1-2 neighbors

  if (newton_bond) {
    for (int i = 0; i < natoms; i++) {
      for (int j = 0; j < num_bond[i]; j++) {
        atom2 = bond_atom[i][j] - 1;
        nspecial[i][0]++;
        nspecial[atom2][0]++;
        if (count[i] >= atom->maxspecial || count[atom2] >= atom->maxspecial)
          error->all(FLERR, fileiarg, "Molecule auto special bond generation overflow" + utils::errorurl(23));
        tmpspecial[i][count[i]++] = atom2 + 1;
        tmpspecial[atom2][count[atom2]++] = i + 1;
      }
    }
  } else {
    for (int i = 0; i < natoms; i++) {
      nspecial[i][0] = num_bond[i];
      for (int j = 0; j < num_bond[i]; j++) {
        atom1 = i;
        atom2 = bond_atom[i][j];
        if (count[atom1] >= atom->maxspecial)
          error->all(FLERR, fileiarg, "Molecule auto special bond generation overflow" + utils::errorurl(23));
        tmpspecial[i][count[atom1]++] = atom2;
      }
    }
  }

  // 1-3 neighbors with no duplicates

  for (int i = 0; i < natoms; i++) nspecial[i][1] = nspecial[i][0];

  int dedup;
  for (int i = 0; i < natoms; i++) {
    for (int m = 0; m < nspecial[i][0]; m++) {
      for (int j = 0; j < nspecial[tmpspecial[i][m] - 1][0]; j++) {
        dedup = 0;
        for (int k = 0; k < count[i]; k++) {
          if (tmpspecial[tmpspecial[i][m] - 1][j] == tmpspecial[i][k] ||
              tmpspecial[tmpspecial[i][m] - 1][j] == i + 1) {
            dedup = 1;
          }
        }
        if (!dedup) {
          if (count[i] >= atom->maxspecial)
            error->all(FLERR, fileiarg, "Molecule auto special bond generation overflow" + utils::errorurl(23));
          tmpspecial[i][count[i]++] = tmpspecial[tmpspecial[i][m] - 1][j];
          nspecial[i][1]++;
        }
      }
    }
  }

  // 1-4 neighbors with no duplicates

  for (int i = 0; i < natoms; i++) nspecial[i][2] = nspecial[i][1];

  for (int i = 0; i < natoms; i++) {
    for (int m = nspecial[i][0]; m < nspecial[i][1]; m++) {
      for (int j = 0; j < nspecial[tmpspecial[i][m] - 1][0]; j++) {
        dedup = 0;
        for (int k = 0; k < count[i]; k++) {
          if (tmpspecial[tmpspecial[i][m] - 1][j] == tmpspecial[i][k] ||
              tmpspecial[tmpspecial[i][m] - 1][j] == i + 1) {
            dedup = 1;
          }
        }
        if (!dedup) {
          if (count[i] >= atom->maxspecial)
            error->all(FLERR, fileiarg, "Molecule auto special bond generation overflow" + utils::errorurl(23));
          tmpspecial[i][count[i]++] = tmpspecial[tmpspecial[i][m] - 1][j];
          nspecial[i][2]++;
        }
      }
    }
  }

  maxspecial = 0;
  for (int i = 0; i < natoms; i++) maxspecial = MAX(maxspecial, nspecial[i][2]);

  memory->create(special, natoms, maxspecial, "molecule:special");
  for (int i = 0; i < natoms; i++)
    for (int j = 0; j < nspecial[i][2]; j++) special[i][j] = tmpspecial[i][j];

  memory->destroy(tmpspecial);
}

/* ----------------------------------------------------------------------
   read SHAKE flags from file
------------------------------------------------------------------------- */

void Molecule::shakeflag_read(char *line)
{
  for (int i = 0; i < natoms; i++) count[i] = 0;
  try {
    for (int i = 0; i < natoms; i++) {
      readline(line);

      ValueTokenizer values(utils::trim_comment(line));

      if (values.count() != 2) error->all(FLERR, fileiarg, "Invalid Shake Flags section in molecule file");

      int iatom = values.next_int() - 1;
      if (iatom < 0 || iatom >= natoms)
        error->all(FLERR, fileiarg, "Invalid atom index in Shake Flags section of molecule file");
      count[iatom]++;
      shake_flag[iatom] = values.next_int();
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, fileiarg, "Invalid Shake Flags section in molecule file: {}", e.what());
  }

  for (int i = 0; i < natoms; i++) {
    if (shake_flag[i] < 0 || shake_flag[i] > 4)
      error->all(FLERR, fileiarg, "Invalid shake flag in molecule file");
    if (count[i] == 0)
      error->all(FLERR, fileiarg, "Atom {} missing in Shake Flags section of molecule file", i + 1);
  }
}

/* ----------------------------------------------------------------------
   read SHAKE atom info from file
------------------------------------------------------------------------- */

void Molecule::shakeatom_read(char *line)
{
  int nmatch = 0, nwant = 0;
  for (int i = 0; i < natoms; i++) count[i] = 0;
  try {
    for (int i = 0; i < natoms; i++) {
      readline(line);

      ValueTokenizer values(utils::trim_comment(line));
      nmatch = values.count();
      int iatom = values.next_int() - 1;
      if ((iatom < 0) || (iatom >= natoms))
        throw TokenizerException(fmt::format("Invalid atom-id {} in Shake Atoms section of "
                                             "molecule file", iatom + 1), "");

      switch (shake_flag[iatom]) {
        case 1:
          shake_atom[iatom][0] = values.next_tagint();
          shake_atom[iatom][1] = values.next_tagint();
          shake_atom[iatom][2] = values.next_tagint();
          nwant = 4;
          break;

        case 2:
          shake_atom[iatom][0] = values.next_tagint();
          shake_atom[iatom][1] = values.next_tagint();
          nwant = 3;
          break;

        case 3:
          shake_atom[iatom][0] = values.next_tagint();
          shake_atom[iatom][1] = values.next_tagint();
          shake_atom[iatom][2] = values.next_tagint();
          nwant = 4;
          break;

        case 4:
          shake_atom[iatom][0] = values.next_tagint();
          shake_atom[iatom][1] = values.next_tagint();
          shake_atom[iatom][2] = values.next_tagint();
          shake_atom[iatom][3] = values.next_tagint();
          nwant = 5;
          break;

        case 0:
          nwant = 1;
          break;

        default:
        throw TokenizerException(
          fmt::format("Unexpected Shake flag {} for atom {} in Shake flags "
                      "section of molecule file", shake_flag[iatom], iatom + 1), "");
      }

      if (nmatch != nwant)
        throw TokenizerException(
          fmt::format("Unexpected number of atom-ids ({} vs {}) for atom {} in Shake Atoms "
                      "section of molecule file", nmatch, nwant, iatom + 1), "");
      count[iatom]++;
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, fileiarg, "Invalid Shake Atoms section in molecule file: {}", e.what());
  }

  for (int i = 0; i < natoms; i++) {
    int m = shake_flag[i];
    if (m == 1) m = 3;
    for (int j = 0; j < m; j++)
      if (shake_atom[i][j] <= 0 || shake_atom[i][j] > natoms)
        error->all(FLERR, fileiarg, "Invalid shake atom in molecule file");
    if (count[i] == 0)
      error->all(FLERR, fileiarg, "Atom {} missing in Shake Atoms section of molecule file", i + 1);
  }
}

/* ----------------------------------------------------------------------
   read SHAKE bond type info from file
------------------------------------------------------------------------- */

void Molecule::shaketype_read(char *line)
{
  int nmatch = 0, nwant = 0;
  for (int i = 0; i < natoms; i++) count[i] = 0;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    auto values = Tokenizer(utils::trim(line)).as_vector();
    nmatch = values.size();
    for (std::size_t ii = 0; ii < values.size(); ++ii) {
      if (utils::strmatch(values[ii], "^#")) {
        nmatch = ii;
        break;
      }
    }
    int iatom = utils::inumeric(FLERR, values[0], false, lmp) - 1;
    if ((iatom < 0) || (iatom >= natoms))
      error->all(FLERR, fileiarg, "Invalid atom-id {} in Skake Bond Types section of molecule file",
                 iatom + 1);
    char *subst;
    switch (shake_flag[iatom]) {
      case 1:
        subst = utils::expand_type(FLERR, values[1], Atom::BOND, lmp);
        if (subst) values[1] = subst;
        shake_type[iatom][0] = utils::inumeric(FLERR, values[1], false, lmp) + ((subst) ? 0 : boffset);
        delete[] subst;

        subst = utils::expand_type(FLERR, values[2], Atom::BOND, lmp);
        if (subst) values[2] = subst;
        shake_type[iatom][1] = utils::inumeric(FLERR, values[2], false, lmp) + ((subst) ? 0 : boffset);
        delete[] subst;

        subst = utils::expand_type(FLERR, values[3], Atom::ANGLE, lmp);
        if (subst) values[3] = subst;
        shake_type[iatom][2] = utils::inumeric(FLERR, values[3], false, lmp) + ((subst) ? 0 : aoffset);
        delete[] subst;

        nwant = 4;
        break;

      case 2:
        subst = utils::expand_type(FLERR, values[1], Atom::BOND, lmp);
        if (subst) values[1] = subst;
        shake_type[iatom][0] = utils::inumeric(FLERR, values[1], false, lmp) + ((subst) ? 0 : boffset);
        delete[] subst;

        nwant = 2;
        break;

      case 3:
        subst = utils::expand_type(FLERR, values[1], Atom::BOND, lmp);
        if (subst) values[1] = subst;
        shake_type[iatom][0] = utils::inumeric(FLERR, values[1], false, lmp) + ((subst) ? 0 : boffset);
        delete[] subst;

        subst = utils::expand_type(FLERR, values[2], Atom::BOND, lmp);
        if (subst) values[2] = subst;
        shake_type[iatom][1] = utils::inumeric(FLERR, values[2], false, lmp) + ((subst) ? 0 : boffset);
        delete[] subst;

        nwant = 3;
        break;

      case 4:
        subst = utils::expand_type(FLERR, values[1], Atom::BOND, lmp);
        if (subst) values[1] = subst;
        shake_type[iatom][0] = utils::inumeric(FLERR, values[1], false, lmp) + ((subst) ? 0 : boffset);
        delete[] subst;

        subst = utils::expand_type(FLERR, values[2], Atom::BOND, lmp);
        if (subst) values[2] = subst;
        shake_type[iatom][1] = utils::inumeric(FLERR, values[2], false, lmp) + ((subst) ? 0 : boffset);
        delete[] subst;

        subst = utils::expand_type(FLERR, values[3], Atom::BOND, lmp);
        if (subst) values[3] = subst;
        shake_type[iatom][2] = utils::inumeric(FLERR, values[3], false, lmp) + ((subst) ? 0 : boffset);
        delete[] subst;

        nwant = 4;
        break;

      case 0:
        nwant = 1;
        break;

      default:
        error->all(FLERR, fileiarg, "Invalid shake type values in molecule file");
    }
    if (nmatch != nwant) error->all(FLERR, fileiarg, "Invalid shake type data in molecule file");
    count[iatom]++;
  }

  for (int i = 0; i < natoms; i++) {
    int m = shake_flag[i];
    if (m == 1) m = 3;
    for (int j = 0; j < m - 1; j++)
      if (shake_type[i][j] <= 0) error->all(FLERR, fileiarg, "Invalid shake bond type in molecule file");
    if (shake_flag[i] == 1)
      if (shake_type[i][2] <= 0) error->all(FLERR, fileiarg, "Invalid shake angle type in molecule file");
    if (count[i] == 0)
      error->all(FLERR, fileiarg, "Atom {} missing in Shake Bond Types section of molecule file", i + 1);
  }
}

/* ----------------------------------------------------------------------
   read body params from file
   pflag = 0/1 for integer/double params
------------------------------------------------------------------------- */

void Molecule::body(int flag, int pflag, char *line)
{
  int nparam = nibody;
  if (pflag) nparam = ndbody;

  int nword = 0;
  try {
    while (nword < nparam) {
      readline(line);

      ValueTokenizer values(utils::trim_comment(line));
      int ncount = values.count();

      std::string sec = "Body ";
      sec += pflag ? "Doubles" : "Integers";
      if (ncount == 0)
        error->all(FLERR, fileiarg, "Too few values in {} section of molecule file", sec);
      if (nword + ncount > nparam)
        error->all(FLERR, fileiarg, "Too many values in {} section of molecule file", sec);

      if (flag) {
        if (pflag == 0) {
          while (values.has_next()) ibodyparams[nword++] = values.next_int();
        } else {
          while (values.has_next()) dbodyparams[nword++] = values.next_double();
        }
      } else
        nword += ncount;
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, fileiarg, "Invalid body params in molecule file: {}", e.what());
  }
}

/* ----------------------------------------------------------------------
   return fragment index if name matches existing fragment, -1 if no such fragment
------------------------------------------------------------------------- */

int Molecule::findfragment(const char *name)
{
  for (int i = 0; i < nfragments; i++)
    if (fragmentnames[i] == name) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   error check molecule attributes and topology against system settings
------------------------------------------------------------------------- */

void Molecule::check_attributes()
{
  // check per-atom attributes of molecule
  // warn if not a match

  int mismatch = 0;
  if (qflag && !atom->q_flag) mismatch = 1;
  if (muflag && !atom->mu_flag) mismatch = 1;
  if (radiusflag && !atom->radius_flag) mismatch = 1;
  if (rmassflag && !atom->rmass_flag) mismatch = 1;
  if (bodyflag && !atom->body_flag) mismatch = 1;

  if (mismatch && (comm->me == 0))
    error->warning(FLERR, "Molecule attributes do not match system attributes"
                   + utils::errorurl(26));

  // for all atom styles, check nbondtype,etc

  mismatch = 0;
  if (atom->nbondtypes < nbondtypes) mismatch = 1;
  if (atom->nangletypes < nangletypes) mismatch = 1;
  if (atom->ndihedraltypes < ndihedraltypes) mismatch = 1;
  if (atom->nimpropertypes < nimpropertypes) mismatch = 1;

  if (mismatch)
    error->all(FLERR, fileiarg, "Molecule topology type exceeds system topology type"
               + utils::errorurl(25));

  // for molecular atom styles, check bond_per_atom,etc + maxspecial
  // do not check for atom style template, since nothing stored per atom

  if (atom->molecular == Atom::MOLECULAR) {
    if (atom->avec->bonds_allow && atom->bond_per_atom < bond_per_atom) mismatch = 1;
    if (atom->avec->angles_allow && atom->angle_per_atom < angle_per_atom) mismatch = 1;
    if (atom->avec->dihedrals_allow && atom->dihedral_per_atom < dihedral_per_atom) mismatch = 1;
    if (atom->avec->impropers_allow && atom->improper_per_atom < improper_per_atom) mismatch = 1;
    if (atom->maxspecial < maxspecial) mismatch = 1;

    if (mismatch)
      error->all(FLERR, fileiarg, "Molecule topology/atom exceeds system topology/atom" + utils::errorurl(24));
  }

  // warn if molecule topology defined but no special settings

  if (bondflag && !specialflag)
    if (comm->me == 0) error->warning(FLERR, "Molecule has bond topology but no special bond settings");
}

/* ----------------------------------------------------------------------
   init all data structures to empty
------------------------------------------------------------------------- */

void Molecule::initialize()
{
  title.clear();
  natoms = 0;
  nbonds = nangles = ndihedrals = nimpropers = 0;
  ntypes = 0;
  nmolecules = 1;
  nbondtypes = nangletypes = ndihedraltypes = nimpropertypes = 0;
  nibody = ndbody = 0;
  nfragments = 0;
  masstotal = 0.0;
  maxradius = 0.0;
  molradius = 0.0;
  comatom = 0;
  maxextent = 0.0;

  nset = 0;
  last = 0;
  fileiarg = 0;

  bond_per_atom = angle_per_atom = dihedral_per_atom = improper_per_atom = 0;
  maxspecial = 0;

  xflag = typeflag = moleculeflag = fragmentflag = qflag = radiusflag = muflag = rmassflag = 0;
  bondflag = angleflag = dihedralflag = improperflag = 0;
  nspecialflag = specialflag = 0;
  shakeflag = shakeflagflag = shakeatomflag = shaketypeflag = 0;
  bodyflag = ibodyflag = dbodyflag = 0;

  centerflag = massflag = comflag = inertiaflag = 0;
  massflag_user = comflag_user = inertiaflag_user = specialflag_user = 0;
  tag_require = 0;

  x = nullptr;
  type = nullptr;
  q = nullptr;
  radius = nullptr;
  rmass = nullptr;

  num_bond = nullptr;
  bond_type = nullptr;
  bond_atom = nullptr;

  num_angle = nullptr;
  angle_type = nullptr;
  angle_atom1 = angle_atom2 = angle_atom3 = nullptr;

  num_dihedral = nullptr;
  dihedral_type = nullptr;
  dihedral_atom1 = dihedral_atom2 = dihedral_atom3 = dihedral_atom4 = nullptr;

  num_improper = nullptr;
  improper_type = nullptr;
  improper_atom1 = improper_atom2 = improper_atom3 = improper_atom4 = nullptr;

  nspecial = nullptr;
  special = nullptr;

  shake_flag = nullptr;
  shake_atom = nullptr;
  shake_type = nullptr;

  ibodyparams = nullptr;
  dbodyparams = nullptr;

  dx = nullptr;
  dxcom = nullptr;
  dxbody = nullptr;
}

/* ----------------------------------------------------------------------
   allocate all data structures
   also initialize values for data structures that are always allocated
------------------------------------------------------------------------- */

void Molecule::allocate()
{
  if (xflag) memory->create(x, natoms, 3, "molecule:x");
  if (typeflag) memory->create(type, natoms, "molecule:type");
  if (moleculeflag) memory->create(molecule, natoms, "molecule:molecule");
  if (fragmentflag) {
    fragmentnames.resize(nfragments);
    memory->create(fragmentmask, nfragments, natoms, "molecule:fragmentmask");
    for (int i = 0; i < nfragments; i++)
      for (int j = 0; j < natoms; j++) fragmentmask[i][j] = 0;
  }
  if (qflag) memory->create(q, natoms, "molecule:q");
  if (muflag) memory->create(mu, natoms, 3, "molecule:mu");
  if (radiusflag) memory->create(radius, natoms, "molecule:radius");
  if (rmassflag) memory->create(rmass, natoms, "molecule:rmass");

  // always allocate num_bond,num_angle,etc and nspecial
  // even if not in molecule file, initialize to 0
  // this is so methods that use these arrays don't have to check they exist

  memory->create(num_bond, natoms, "molecule:num_bond");
  for (int i = 0; i < natoms; i++) num_bond[i] = 0;
  memory->create(num_angle, natoms, "molecule:num_angle");
  for (int i = 0; i < natoms; i++) num_angle[i] = 0;
  memory->create(num_dihedral, natoms, "molecule:num_dihedral");
  for (int i = 0; i < natoms; i++) num_dihedral[i] = 0;
  memory->create(num_improper, natoms, "molecule:num_improper");
  for (int i = 0; i < natoms; i++) num_improper[i] = 0;
  memory->create(nspecial, natoms, 3, "molecule:nspecial");
  for (int i = 0; i < natoms; i++) nspecial[i][0] = nspecial[i][1] = nspecial[i][2] = 0;

  if (specialflag) memory->create(special, natoms, maxspecial, "molecule:special");

  if (bondflag) {
    memory->create(bond_type, natoms, bond_per_atom, "molecule:bond_type");
    memory->create(bond_atom, natoms, bond_per_atom, "molecule:bond_atom");
  }

  if (angleflag) {
    memory->create(angle_type, natoms, angle_per_atom, "molecule:angle_type");
    memory->create(angle_atom1, natoms, angle_per_atom, "molecule:angle_atom1");
    memory->create(angle_atom2, natoms, angle_per_atom, "molecule:angle_atom2");
    memory->create(angle_atom3, natoms, angle_per_atom, "molecule:angle_atom3");
  }

  if (dihedralflag) {
    memory->create(dihedral_type, natoms, dihedral_per_atom, "molecule:dihedral_type");
    memory->create(dihedral_atom1, natoms, dihedral_per_atom, "molecule:dihedral_atom1");
    memory->create(dihedral_atom2, natoms, dihedral_per_atom, "molecule:dihedral_atom2");
    memory->create(dihedral_atom3, natoms, dihedral_per_atom, "molecule:dihedral_atom3");
    memory->create(dihedral_atom4, natoms, dihedral_per_atom, "molecule:dihedral_atom4");
  }

  if (improperflag) {
    memory->create(improper_type, natoms, improper_per_atom, "molecule:improper_type");
    memory->create(improper_atom1, natoms, improper_per_atom, "molecule:improper_atom1");
    memory->create(improper_atom2, natoms, improper_per_atom, "molecule:improper_atom2");
    memory->create(improper_atom3, natoms, improper_per_atom, "molecule:improper_atom3");
    memory->create(improper_atom4, natoms, improper_per_atom, "molecule:improper_atom4");
  }

  if (shakeflag) {
    memory->create(shake_flag, natoms, "molecule:shake_flag");
    memory->create(shake_atom, natoms, 4, "molecule:shake_flag");
    memory->create(shake_type, natoms, 3, "molecule:shake_flag");
  }

  if (bodyflag) {
    if (nibody) memory->create(ibodyparams, nibody, "molecule:ibodyparams");
    if (ndbody) memory->create(dbodyparams, ndbody, "molecule:dbodyparams");
  }
}

/* ----------------------------------------------------------------------
   deallocate all data structures
------------------------------------------------------------------------- */

void Molecule::deallocate()
{
  memory->destroy(x);
  memory->destroy(type);
  memory->destroy(molecule);
  memory->destroy(q);
  memory->destroy(mu);
  memory->destroy(radius);
  memory->destroy(rmass);

  memory->destroy(molecule);
  memory->destroy(fragmentmask);

  if (fragmentflag) { fragmentnames.clear(); }

  memory->destroy(num_bond);
  memory->destroy(bond_type);
  memory->destroy(bond_atom);

  memory->destroy(num_angle);
  memory->destroy(angle_type);
  memory->destroy(angle_atom1);
  memory->destroy(angle_atom2);
  memory->destroy(angle_atom3);

  memory->destroy(num_dihedral);
  memory->destroy(dihedral_type);
  memory->destroy(dihedral_atom1);
  memory->destroy(dihedral_atom2);
  memory->destroy(dihedral_atom3);
  memory->destroy(dihedral_atom4);

  memory->destroy(num_improper);
  memory->destroy(improper_type);
  memory->destroy(improper_atom1);
  memory->destroy(improper_atom2);
  memory->destroy(improper_atom3);
  memory->destroy(improper_atom4);

  memory->destroy(nspecial);
  memory->destroy(special);

  memory->destroy(shake_flag);
  memory->destroy(shake_atom);
  memory->destroy(shake_type);

  memory->destroy(dx);
  memory->destroy(dxcom);
  memory->destroy(dxbody);

  memory->destroy(ibodyparams);
  memory->destroy(dbodyparams);
}

/* ----------------------------------------------------------------------
   read and bcast a line
------------------------------------------------------------------------- */

void Molecule::readline(char *line)
{
  int n;
  if (comm->me == 0) {
    if (fgets(line, MAXLINE, fp) == nullptr)
      n = 0;
    else
      n = strlen(line) + 1;
  }
  MPI_Bcast(&n, 1, MPI_INT, 0, world);
  if (n == 0) error->all(FLERR, fileiarg, "Unexpected end of molecule file");
  MPI_Bcast(line, n, MPI_CHAR, 0, world);
}

/* ----------------------------------------------------------------------
   extract keyword from line
   flag = 0, read and bcast line
   flag = 1, line has already been read
------------------------------------------------------------------------- */

std::string Molecule::parse_keyword(int flag, char *line)
{
  char line2[MAXLINE] = {'\0'};
  if (flag) {

    // read upto non-blank line plus 1 following line
    // eof is set to 1 if any read hits end-of-file

    int eof = 0;
    if (comm->me == 0) {
      if (fgets(line, MAXLINE, fp) == nullptr) eof = 1;
      while (eof == 0 && strspn(line, " \t\n\r") == strlen(line)) {
        if (fgets(line, MAXLINE, fp) == nullptr) eof = 1;
      }
      if (fgets(line2, MAXLINE, fp) == nullptr) eof = 1;
    }

    // if eof, set keyword empty and return

    MPI_Bcast(&eof, 1, MPI_INT, 0, world);
    if (eof) return {""};

    // bcast keyword line to all procs

    MPI_Bcast(line, MAXLINE, MPI_CHAR, 0, world);
  }

  // return non-whitespace and non-comment portion of line

  return utils::trim(utils::trim_comment(line));
}

/* ----------------------------------------------------------------------
   skip N lines of file. Check if non-numeric content (e.g. keyword).
------------------------------------------------------------------------- */

void Molecule::skip_lines(int n, char *line, const std::string &section)
{
  for (int i = 0; i < n; i++) {
    readline(line);
    if (utils::strmatch(utils::trim(utils::trim_comment(line)), "^[A-Za-z ]+$"))
      error->one(FLERR, Error::NOLASTLINE,
                 "Unexpected line in molecule file while skipping {} section:\n{}",
                 section, line);
  }
}

/* ------------------------------------------------------------------------------ */

void Molecule::stats()
{
  if (title.empty()) title = "(no title)";
  if (comm->me == 0)
    utils::logmesg(lmp,
                   "Read molecule template {}:\n{}\n"
                   "  {} molecules\n"
                   "  {} fragments\n"
                   "  {} bodies\n"
                   "  {} atoms with max type {}\n"
                   "  {} bonds with max type {}\n"
                   "  {} angles with max type {}\n"
                   "  {} dihedrals with max type {}\n"
                   "  {} impropers with max type {}\n",
                   id, title, nmolecules, nfragments, bodyflag, natoms, ntypes, nbonds, nbondtypes,
                   nangles, nangletypes, ndihedrals, ndihedraltypes, nimpropers, nimpropertypes);
}

// clang-format on
/* ----------------------------------------------------------------------
   print molecule file. may only be called from MPI rank 0
------------------------------------------------------------------------- */

void Molecule::print(FILE *fp)
{
  utils::print(fp, "  {} atoms\n", natoms);
  if (nbonds) utils::print(fp, "  {} bonds\n", nbonds);
  if (nangles) utils::print(fp, "  {} angles\n", nangles);
  if (ndihedrals) utils::print(fp, "  {} dihedrals\n", ndihedrals);
  if (nimpropers) utils::print(fp, "  {} impropers\n", nimpropers);
  if (nfragments) utils::print(fp, "  {} fragments\n", nfragments);
  if (massflag_user) utils::print(fp, "  {} mass\n", masstotal);
  if (bodyflag) utils::print(fp, "  {} {} body\n", nibody, ndbody);
  if (comflag_user) utils::print(fp, "  {} {} {} com\n", com[0], com[1], com[2]);
  if (inertiaflag_user)
    utils::print(fp, "  {} {} {} {} {} {} inertia\n", itensor[0], itensor[1], itensor[2],
                 itensor[3], itensor[4], itensor[5]);

  if (xflag) {
    fputs("\nCoords\n\n", fp);
    for (int i = 0; i < natoms; i++)
      utils::print(fp, " {}  {} {} {}\n", i + 1, x[i][0], x[i][1], x[i][2]);
  }

  if (typeflag) {
    fputs("\nTypes\n\n", fp);
    if (atom->labelmapflag && atom->lmap->is_complete(Atom::ATOM)) {
      for (int i = 0; i < natoms; i++)
        utils::print(fp, " {} {}\n", i + 1, atom->lmap->find(type[i], Atom::ATOM));
    } else {
      for (int i = 0; i < natoms; i++) utils::print(fp, " {}  {}\n", i + 1, type[i]);
    }
  }

  if (moleculeflag) {
    fputs("\nMolecules\n\n", fp);
    for (int i = 0; i < natoms; i++) utils::print(fp, " {}  {}\n", i + 1, molecule[i]);
  }

  if (fragmentflag) {
    fputs("\nFragments\n\n", fp);
    for (int i = 0; i < nfragments; i++) {
      utils::print(fp, " {} ", fragmentnames[i]);
      for (int j = 0; j < natoms; j++) {
        if (fragmentmask[i][j]) utils::print(fp, " {}", j + 1);
      }
      fputs("\n", fp);
    }
  }

  if (qflag) {
    fputs("\nCharges\n\n", fp);
    for (int i = 0; i < natoms; i++) utils::print(fp, " {}  {}\n", i + 1, q[i]);
  }

  if (radiusflag && !bodyflag) {
    fputs("\nDiameters\n\n", fp);
    for (int i = 0; i < natoms; i++) utils::print(fp, " {}  {}\n", i + 1, 2.0 * radius[i]);
  }

  if (muflag) {
    fputs("\nDipoles\n\n", fp);
    for (int i = 0; i < natoms; i++)
      utils::print(fp, " {}  {} {} {}\n", i + 1, mu[i][0], mu[i][1], mu[i][2]);
  }

  if (rmassflag) {
    fputs("\nMasses\n\n", fp);
    for (int i = 0; i < natoms; i++) utils::print(fp, " {}  {}\n", i + 1, rmass[i]);
  }

  bool has_newton_bond = force->newton_bond > 0;

  if (bondflag) {
    fputs("\nBonds\n\n", fp);
    int idx = 0;
    bool has_typelabels = (atom->labelmapflag != 0) && atom->lmap->is_complete(Atom::BOND);

    for (int i = 0; i < natoms; i++) {
      for (int j = 0; j < num_bond[i]; j++) {
        if (has_newton_bond || (i + 1 < bond_atom[i][j])) {
          ++idx;
          if (has_typelabels) {
            utils::print(fp, " {}  {}", idx, atom->lmap->find(bond_type[i][j], Atom::BOND));
          } else {
            utils::print(fp, " {}  {}", idx, bond_type[i][j]);
          }
          utils::print(fp, " {} {}\n", i + 1, bond_atom[i][j]);
        }
      }
    }
  }

  if (angleflag) {
    fputs("\nAngles\n\n", fp);
    int idx = 0;
    bool has_typelabels = (atom->labelmapflag != 0) && atom->lmap->is_complete(Atom::ANGLE);
    for (int i = 0; i < natoms; i++) {
      for (int j = 0; j < num_angle[i]; j++) {
        if (has_newton_bond || (i + 1 == angle_atom2[i][j])) {
          ++idx;
          if (has_typelabels) {
            utils::print(fp, " {}  {}", idx, atom->lmap->find(angle_type[i][j], Atom::ANGLE));
          } else {
            utils::print(fp, " {}  {}", idx, angle_type[i][j]);
          }
          utils::print(fp, " {} {} {}\n", angle_atom1[i][j], angle_atom2[i][j], angle_atom3[i][j]);
        }
      }
    }
  }

  if (dihedralflag) {
    fputs("\nDihedrals\n\n", fp);
    int idx = 0;
    bool has_typelabels = (atom->labelmapflag != 0) && atom->lmap->is_complete(Atom::DIHEDRAL);
    for (int i = 0; i < natoms; i++) {
      for (int j = 0; j < num_dihedral[i]; j++) {
        if (has_newton_bond || (i + 1 == dihedral_atom2[i][j])) {
          ++idx;
          if (has_typelabels) {
            utils::print(fp, " {}  {}", idx, atom->lmap->find(dihedral_type[i][j], Atom::DIHEDRAL));
          } else {
            utils::print(fp, " {}  {}", idx, dihedral_type[i][j]);
          }
          utils::print(fp, " {} {} {} {}\n", dihedral_atom1[i][j], dihedral_atom2[i][j],
                       dihedral_atom3[i][j], dihedral_atom4[i][j]);
        }
      }
    }
  }

  if (improperflag) {
    fputs("\nImpropers\n\n", fp);
    int idx = 0;
    bool has_typelabels = (atom->labelmapflag != 0) && atom->lmap->is_complete(Atom::IMPROPER);
    for (int i = 0; i < natoms; i++) {
      for (int j = 0; j < num_improper[i]; j++) {
        if (has_newton_bond || (i + 1 == improper_atom2[i][j])) {
          ++idx;
          if (has_typelabels) {
            utils::print(fp, " {}  {}", idx, atom->lmap->find(improper_type[i][j], Atom::IMPROPER));
          } else {
            utils::print(fp, " {}  {}", idx, improper_type[i][j]);
          }
          utils::print(fp, " {} {} {} {}\n", improper_atom1[i][j], improper_atom2[i][j],
                       improper_atom3[i][j], improper_atom4[i][j]);
        }
      }
    }
  }

  if (specialflag_user) {
    fputs("\nSpecial Bond Counts\n\n", fp);
    for (int i = 0; i < natoms; i++) {
      utils::print(fp, " {}  {} {} {}\n", i + 1, nspecial[i][0], nspecial[i][1] - nspecial[i][0],
                   nspecial[i][2] - nspecial[i][1]);
    }

    fputs("\nSpecial Bonds\n\n", fp);
    for (int i = 0; i < natoms; i++) {
      utils::print(fp, " {} ", i + 1);
      for (int j = 0; j < nspecial[i][2]; j++) utils::print(fp, " {}", special[i][j]);
      utils::print(fp, "\n");
    }
  }

  if (shakeflag) {
    fputs("\nShake Flags\n\n", fp);
    for (int i = 0; i < natoms; i++) { utils::print(fp, " {}  {}\n", i + 1, shake_flag[i]); }

    fputs("\nShake Atoms\n\n", fp);
    for (int i = 0; i < natoms; i++) {
      utils::print(fp, " {} ", i + 1);
      switch (shake_flag[i]) {
        case 1:
          utils::print(fp, " {} {} {}\n", shake_atom[i][0], shake_atom[i][1], shake_atom[i][2]);
          break;

        case 2:
          utils::print(fp, " {} {}\n", shake_atom[i][0], shake_atom[i][1]);
          break;

        case 3:
          utils::print(fp, " {} {} {}\n", shake_atom[i][0], shake_atom[i][1], shake_atom[i][2]);
          break;

        case 4:
          utils::print(fp, " {} {} {} {}\n", shake_atom[i][0], shake_atom[i][1], shake_atom[i][2],
                       shake_atom[i][3]);
          break;

        case 0:
          fputs("\n", fp);
          break;
      }
    }

    fputs("\nShake Bond Types\n\n", fp);
    for (int i = 0; i < natoms; i++) {
      bool has_typelabels = (atom->labelmapflag != 0) && atom->lmap->is_complete(Atom::BOND);
      utils::print(fp, " {} ", i + 1);
      switch (shake_flag[i]) {
        case 1:
          has_typelabels = has_typelabels && atom->lmap->is_complete(Atom::ANGLE);
          if (has_typelabels) {
            utils::print(fp, " {} {} {}\n", atom->lmap->find(shake_type[i][0], Atom::BOND),
                         atom->lmap->find(shake_type[i][1], Atom::BOND),
                         atom->lmap->find(shake_type[i][2], Atom::ANGLE));
          } else {
            utils::print(fp, " {} {} {}\n", shake_type[i][0], shake_type[i][1], shake_type[i][2]);
          }
          break;

        case 2:
          if (has_typelabels) {
            utils::print(fp, " {} {}\n", atom->lmap->find(shake_type[i][0], Atom::BOND),
                         atom->lmap->find(shake_type[i][1], Atom::BOND));
          } else {
            utils::print(fp, " {} {}\n", shake_type[i][0], shake_type[i][1]);
          }
          break;

        case 3:
          if (has_typelabels) {
            utils::print(fp, " {} {} {}\n", atom->lmap->find(shake_type[i][0], Atom::BOND),
                         atom->lmap->find(shake_type[i][1], Atom::BOND),
                         atom->lmap->find(shake_type[i][2], Atom::BOND));
          } else {
            utils::print(fp, " {} {} {}\n", shake_type[i][0], shake_type[i][1], shake_type[i][2]);
          }
          break;

        case 4:
          if (has_typelabels) {
            utils::print(fp, " {} {} {} {}\n", atom->lmap->find(shake_type[i][0], Atom::BOND),
                         atom->lmap->find(shake_type[i][1], Atom::BOND),
                         atom->lmap->find(shake_type[i][2], Atom::BOND),
                         atom->lmap->find(shake_type[i][3], Atom::BOND));
          } else {
            utils::print(fp, " {} {} {} {}\n", shake_type[i][0], shake_type[i][1], shake_type[i][2],
                         shake_type[i][3]);
          }
          break;

        case 0:
          fputs("\n", fp);
          break;
      }
    }
  }

  if (bodyflag) {
    auto *avec = dynamic_cast<AtomVecBody *>(atom->style_match("body"));
    if (avec) {
      fputs("\nBody Integers\n\n", fp);
      if ((strcmp(avec->bptr->style, "nparticle") == 0) ||
          (strcmp(avec->bptr->style, "rounded/polygon") == 0)) {
        utils::print(fp, " {}\n", ibodyparams[0]);
      }
      if (strcmp(avec->bptr->style, "rounded/polyhedron") == 0) {
        utils::print(fp, " {} {} {}\n", ibodyparams[0], ibodyparams[1], ibodyparams[2]);
      }
      fputs("\nBody Doubles\n\n", fp);
      utils::print(fp, " {} {} {} {} {} {}\n", dbodyparams[0], dbodyparams[1], dbodyparams[2],
                   dbodyparams[3], dbodyparams[4], dbodyparams[5]);
      int idx = 6;
      for (int i = 0; i < ibodyparams[0]; ++i) {
        utils::print(fp, " {} {} {}\n", dbodyparams[idx], dbodyparams[idx + 1],
                     dbodyparams[idx + 2]);
        idx += 3;
      }
      if (strcmp(avec->bptr->style, "rounded/polyhedron") == 0) {
        for (int i = 0; i < ibodyparams[1]; ++i) {
          utils::print(fp, " {} {}\n", dbodyparams[idx], dbodyparams[idx + 1]);
          idx += 2;
        }
        for (int i = 0; i < ibodyparams[2]; ++i) {
          utils::print(fp, " {} {} {} {}\n", dbodyparams[idx], dbodyparams[idx + 1],
                       dbodyparams[idx + 2], dbodyparams[idx + 3]);
          idx += 4;
        }
      }
      utils::print(fp, " {}\n", dbodyparams[idx]);
    }
  }
}
// clang-format off
