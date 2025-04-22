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

#include "molecule_json.h"

#include "json.h"
#include "atom.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "label_map.h"


// PIMPL-ify json class

namespace LAMMPS_NS {
  struct JsonImpl {
    JsonImpl() : jsonf(nullptr) {}
    ~JsonImpl() = default;
    json jsonf;
  };
}

using namespace LAMMPS_NS;

MoleculeJsonCommand::MoleculeJson::MoleculeJson(LAMMPS *lmp, int narg, char **arg, int &index) : Molecule(lmp, narg, arg, index)
{
  json_impl = new JsonImpl();
}

MoleculeJsonCommand::MoleculeJson::~MoleculeJson() 
{
  delete json_impl;
}

/* ----------------------------------------------------------------------
   read molecule info from json file
   flag = 0, just scan for sizes of fields
   flag = 1, read and store fields
------------------------------------------------------------------------- */

void MoleculeJsonCommand::MoleculeJson::read(int flag)
{
  if (flag == 0) {
    json_impl->jsonf = json::parse(fp);

    natoms = json_impl->jsonf["num_atoms"];

    if(json_impl->jsonf.contains("bonds")) {
      bondflag = tag_require = 1;
      nbonds = json_impl->jsonf["bonds"]["num_bonds"];
    }

    if(json_impl->jsonf.contains("angles")) {
      angleflag = tag_require = 1;
      nangles = json_impl->jsonf["angles"]["num_angles"];
    }

    if(json_impl->jsonf.contains("fragments")) {
      nfragments = json_impl->jsonf["fragments"]["num_fragments"];
      fragmentflag = 1;
    }

    if(json_impl->jsonf.contains("dihedrals")) {
      dihedralflag = tag_require = 1;
      ndihedrals = json_impl->jsonf["dihedrals"]["num_dihedrals"];
    }

    if(json_impl->jsonf.contains("impropers")) {
      improperflag = tag_require = 1;
      nimpropers = json_impl->jsonf["impropers"]["num_impropers"];
    }

    if(json_impl->jsonf.contains("coords")) xflag = 1;

    if(json_impl->jsonf.contains("charges")) qflag = 1;

    if(json_impl->jsonf.contains("types")) typeflag = 1;

    if(json_impl->jsonf.contains("molecules")) moleculeflag = 1;

    // error checks
    if (natoms < 1) error->all(FLERR, fileiarg, "No atoms or invalid atom count in molecule file");
    if (nbonds < 0) error->all(FLERR, fileiarg, "Invalid bond count in molecule file");
    if (nangles < 0) error->all(FLERR, fileiarg, "Invalid angle count in molecule file");
    if (ndihedrals < 0) error->all(FLERR, fileiarg, "Invalid dihedral count in molecule file");
    if (nimpropers < 0) error->all(FLERR, fileiarg, "Invalid improper count in molecule file");

    memory->create(count, natoms, "molecule:count");
  }

  if(improperflag) impropers(flag);
  if(bondflag) bonds(flag);
  if(angleflag) angles(flag);
  if(dihedralflag) dihedrals(flag);

  if (flag) {
    // sections of molecule file
    if(xflag) coords();
    if(typeflag) types();
    if(moleculeflag) molecules();
    if(qflag) charges();

    // auto-generate special bonds if needed and not in file
    if (bondflag && specialflag == 0) {
      if (domain->box_exist == 0)
        error->all(FLERR, fileiarg,
                  "Cannot auto-generate special bonds before simulation box is defined");

      special_generate();
      specialflag = 1;
      nspecialflag = 1;
    }
  }
}

/* ----------------------------------------------------------------------
   read coords from json file
------------------------------------------------------------------------- */

void MoleculeJsonCommand::MoleculeJson::coords()
{
  json data = json_impl->jsonf["coords"]["data"];

  for (int i = 0; i < natoms; i++) count[i] = 0;

  for(int iatom = 0; iatom < data.size(); iatom++) {
    x[iatom][0] = data[iatom][1];
    x[iatom][1] = data[iatom][2];
    x[iatom][2] = data[iatom][3];
    count[iatom]++;
    x[iatom][0] *= sizescale;
    x[iatom][1] *= sizescale;
    x[iatom][2] *= sizescale;
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
   read molecules from json file
   set nmolecules = max of any molecule type
------------------------------------------------------------------------- */

void MoleculeJsonCommand::MoleculeJson::molecules()
{
  json data = json_impl->jsonf["molecules"]["data"];

  for (int i = 0; i < natoms; i++) count[i] = 0;

  for (int iatom = 0; iatom < data.size(); iatom++) {
    molecule[iatom] = data[iatom][1];
    count[iatom]++;
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
   read types from json file
------------------------------------------------------------------------- */

void MoleculeJsonCommand::MoleculeJson::types() {
  json data = json_impl->jsonf["types"]["data"];
  std::string typestr;

  for (int i = 0; i < natoms; i++) count[i] = 0;

  for (int iatom = 0; iatom < natoms; iatom++) {
    typestr = data[iatom][1];
    type[iatom] = atom->lmap->find(typestr, Atom::ATOM);
    count[iatom]++;
  }

  for (int i = 0; i < natoms; i++) {
    if (count[i] == 0) error->all(FLERR, fileiarg, "Atom {} missing", i + 1);
    if ((type[i] <= 0) || (domain->box_exist && (type[i] > atom->ntypes)))
      error->all(FLERR, fileiarg, "Invalid atom type {} for atom {} in molecule file", type[i],
                i + 1);
    ntypes = MAX(ntypes, type[i]);
  }
}

/* ----------------------------------------------------------------------
   read charges from json file
------------------------------------------------------------------------- */

void MoleculeJsonCommand::MoleculeJson::charges()
{
  json data = json_impl->jsonf["charges"]["data"];

  for (int i = 0; i < natoms; i++) count[i] = 0;

  for (int iatom = 0; iatom < natoms; iatom++) {
    q[iatom] = data[iatom][1];
    count[iatom]++;
  }

  for (int i = 0; i < natoms; i++) {
    if (count[i] == 0)
      error->all(FLERR, fileiarg, "Atom {} missing in Charges section of molecule file", i + 1);
  }
}

/* ----------------------------------------------------------------------
   read bonds from json file
   set nbondtypes = max type of any bond
   store each with both atoms if newton_bond = 0
   if flag = 0, just count bonds/atom
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void MoleculeJsonCommand::MoleculeJson::bonds(int flag)
{
  json data = json_impl->jsonf["bonds"]["data"];
  std::string typestr;
  tagint atom1, atom2;

  for (int i = 0; i < nbonds; i++) {
    typestr = data[i][1];
    atom1 = data[i][2];
    atom2 = data[i][3];

    store_bond(flag, typestr, atom1, atom2);
  }

  // bond_per_atom = max of count vector

  if (flag == 0) {
    bond_per_atom = 0;
    for (int i = 0; i < natoms; i++) bond_per_atom = MAX(bond_per_atom, count[i]);
  }
}

/* ----------------------------------------------------------------------
   read angles from json file
   set nbondtypes = max type of any bond
   store each with both atoms if newton_bond = 0
   if flag = 0, just count bonds/atom
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void MoleculeJsonCommand::MoleculeJson::angles(int flag)
{
  json data = json_impl->jsonf["angles"]["data"];
  std::string typestr;
  tagint atom1, atom2, atom3;

  for (int i = 0; i < nangles; i++) {
    typestr = data[i][1];
    atom1 = data[i][2];
    atom2 = data[i][3];
    atom3 = data[i][4];

    store_angle(flag, typestr, atom1, atom2, atom3);
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

void MoleculeJsonCommand::MoleculeJson::dihedrals(int flag)
{
  json data = json_impl->jsonf["dihedrals"]["data"];
  std::string typestr;
  tagint atom1, atom2, atom3, atom4;

  for (int i = 0; i < ndihedrals; i++) {
    typestr = data[i][1];
    atom1 = data[i][2];
    atom2 = data[i][3];
    atom3 = data[i][4];
    atom4 = data[i][5];

    store_dihedral(flag, typestr, atom1, atom2, atom3, atom4);
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

void MoleculeJsonCommand::MoleculeJson::impropers(int flag)
{
  json data = json_impl->jsonf["impropers"]["data"];
  std::string typestr;
  tagint atom1, atom2, atom3, atom4;

  for(int i = 0; i < nimpropers; i++) {
    typestr = data[i][1];
    atom1 = data[i][2];
    atom2 = data[i][3];
    atom3 = data[i][4];
    atom4 = data[i][5];

    store_improper(flag, typestr, atom1, atom2, atom3, atom4);
  }

  // improper_per_atom = max of count vector

  if (flag == 0) {
    improper_per_atom = 0;
    for (int i = 0; i < natoms; i++) improper_per_atom = MAX(improper_per_atom, count[i]);
  }
}

void MoleculeJsonCommand::command(int narg, char **arg)
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "molecule", error);

  if (atom->find_molecule(arg[0]) >= 0)
    error->all(FLERR, Error::ARGZERO, "Reuse of molecule template ID {}", arg[0]);

  // 1st molecule in set stores nset = # of mols, others store nset = 0
  // ifile = count of molecules in set
  // index = argument index where next molecule starts, updated by constructor

  int ifile = 1;
  int index = 1;
  while (true) {
    atom->molecules = (Molecule **)
      memory->srealloc(atom->molecules,(atom->nmolecule+1)*sizeof(MoleculeJson *), "atom::molecules");
    atom->molecules[atom->nmolecule] = new MoleculeJson(lmp,narg,arg,index);
    atom->molecules[atom->nmolecule]->scan(arg);
    atom->molecules[atom->nmolecule]->nset = 0;
    atom->molecules[atom->nmolecule-ifile+1]->nset++;
    atom->nmolecule++;
    if (atom->molecules[atom->nmolecule-1]->last) break;
    ifile++;
  }
}
