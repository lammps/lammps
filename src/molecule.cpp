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

#include "molecule.h"

#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_body.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"
#include "math_eigen.h"
#include "memory.h"
#include "tokenizer.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define MAXLINE 256
#define EPSILON 1.0e-7
#define BIG 1.0e20

#define SINERTIA 0.4            // moment of inertia prefactor for sphere

/* ---------------------------------------------------------------------- */

Molecule::Molecule(LAMMPS *lmp, int narg, char **arg, int &index) :
  Pointers(lmp), id(nullptr), x(nullptr), type(nullptr), molecule(nullptr), q(nullptr), radius(nullptr),
  rmass(nullptr), num_bond(nullptr), bond_type(nullptr), bond_atom(nullptr),
  num_angle(nullptr), angle_type(nullptr), angle_atom1(nullptr), angle_atom2(nullptr),
  angle_atom3(nullptr), num_dihedral(nullptr), dihedral_type(nullptr), dihedral_atom1(nullptr),
  dihedral_atom2(nullptr), dihedral_atom3(nullptr), dihedral_atom4(nullptr), num_improper(nullptr),
  improper_type(nullptr), improper_atom1(nullptr), improper_atom2(nullptr),
  improper_atom3(nullptr), improper_atom4(nullptr), nspecial(nullptr), special(nullptr),
  shake_flag(nullptr), shake_atom(nullptr), shake_type(nullptr), avec_body(nullptr), ibodyparams(nullptr),
  dbodyparams(nullptr), fragmentmask(nullptr),
  dx(nullptr), dxcom(nullptr), dxbody(nullptr), quat_external(nullptr), fp(nullptr), count(nullptr)
{
  me = comm->me;

  if (index >= narg) error->all(FLERR,"Illegal molecule command");

  id = utils::strdup(arg[0]);
  if (!utils::is_id(id))
    error->all(FLERR,"Molecule template ID must have only "
               "alphanumeric or underscore characters");

  // parse args until reach unknown arg (next file)

  toffset = 0;
  boffset = aoffset = doffset = ioffset = 0;
  sizescale = 1.0;

  int ifile = index;
  int iarg = ifile+1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"offset") == 0) {
      if (iarg+6 > narg) error->all(FLERR,"Illegal molecule command");
      toffset = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      boffset = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      aoffset = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
      doffset = utils::inumeric(FLERR,arg[iarg+4],false,lmp);
      ioffset = utils::inumeric(FLERR,arg[iarg+5],false,lmp);
      if (toffset < 0 || boffset < 0 || aoffset < 0 ||
          doffset < 0 || ioffset < 0)
        error->all(FLERR,"Illegal molecule command");
      iarg += 6;
    } else if (strcmp(arg[iarg],"toff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal molecule command");
      toffset = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (toffset < 0) error->all(FLERR,"Illegal molecule command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"boff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal molecule command");
      boffset = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (boffset < 0) error->all(FLERR,"Illegal molecule command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"aoff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal molecule command");
      aoffset = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (aoffset < 0) error->all(FLERR,"Illegal molecule command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"doff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal molecule command");
      doffset = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (doffset < 0) error->all(FLERR,"Illegal molecule command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"ioff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal molecule command");
      ioffset = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (ioffset < 0) error->all(FLERR,"Illegal molecule command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal molecule command");
      sizescale = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (sizescale <= 0.0) error->all(FLERR,"Illegal molecule command");
      iarg += 2;
    } else break;
  }

  index = iarg;

  // last molecule if have scanned all args

  if (iarg == narg) last = 1;
  else last = 0;

  // initialize all fields to empty

  initialize();

  // scan file for sizes of all fields and allocate storage for them

  if (me == 0) {
    fp = fopen(arg[ifile],"r");
    if (fp == nullptr)
      error->one(FLERR,"Cannot open molecule file {}: {}",
                                   arg[ifile], utils::getsyserror());
  }
  read(0);
  if (me == 0) fclose(fp);
  allocate();

  // read file again to populate all fields

  if (me == 0) fp = fopen(arg[ifile],"r");
  read(1);
  if (me == 0) fclose(fp);

  // stats

  if (me == 0)
    utils::logmesg(lmp,"Read molecule template {}:\n  {} molecules\n"
                   "  {} fragments\n"
                   "  {} atoms with max type {}\n"
                   "  {} bonds with max type {}\n"
                   "  {} angles with max type {}\n"
                   "  {} dihedrals with max type {}\n"
                   "  {} impropers with max type {}\n", id,nmolecules,
                   nfragments,natoms,ntypes,nbonds,nbondtypes,nangles,
                   nangletypes,ndihedrals,ndihedraltypes,nimpropers,
                   nimpropertypes);
}

/* ---------------------------------------------------------------------- */

Molecule::~Molecule()
{
  delete [] id;
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
  memory->create(dx,natoms,3,"molecule:dx");

  for (int i = 0; i < natoms; i++) {
    dx[i][0] = x[i][0] - center[0];
    dx[i][1] = x[i][1] - center[1];
    dx[i][2] = x[i][2] - center[2];
  }

  molradius = 0.0;
  for (int i = 0; i < natoms; i++) {
    double rad = MathExtra::len3(dx[i]);
    if (radiusflag) rad += radius[i];
    molradius = MAX(molradius,rad);
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
    if (rmassflag) masstotal += rmass[i];
    else masstotal += atom->mass[type[i]];
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
      if (rmassflag) onemass = rmass[i];
      else onemass = atom->mass[type[i]];
      com[0] += x[i][0]*onemass;
      com[1] += x[i][1]*onemass;
      com[2] += x[i][2]*onemass;
    }
    if (masstotal > 0.0) {
      com[0] /= masstotal;
      com[1] /= masstotal;
      com[2] /= masstotal;
    }
  }

  memory->destroy(dxcom);
  memory->create(dxcom,natoms,3,"molecule:dxcom");

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
    double rsq = dx*dx + dy*dy + dz*dz;
    rsqmax = MAX(rsqmax,rsq);
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

    double onemass,dx,dy,dz;
    for (int i = 0; i < 6; i++) itensor[i] = 0.0;
    for (int i = 0; i < natoms; i++) {
      if (rmassflag) onemass = rmass[i];
      else onemass = atom->mass[type[i]];
      dx = dxcom[i][0];
      dy = dxcom[i][1];
      dz = dxcom[i][2];
      itensor[0] += onemass * (dy*dy + dz*dz);
      itensor[1] += onemass * (dx*dx + dz*dz);
      itensor[2] += onemass * (dx*dx + dy*dy);
      itensor[3] -= onemass * dy*dz;
      itensor[4] -= onemass * dx*dz;
      itensor[5] -= onemass * dx*dy;
    }

    if (radiusflag) {
      for (int i = 0; i < natoms; i++) {
        if (rmassflag) onemass = rmass[i];
        else onemass = atom->mass[type[i]];
        itensor[0] += SINERTIA*onemass * radius[i]*radius[i];
        itensor[1] += SINERTIA*onemass * radius[i]*radius[i];
        itensor[2] += SINERTIA*onemass * radius[i]*radius[i];
      }
    }
  }

  // diagonalize inertia tensor for each body via Jacobi rotations
  // inertia = 3 eigenvalues = principal moments of inertia
  // evectors and exzy = 3 evectors = principal axes of rigid body

  double cross[3];
  double tensor[3][3],evectors[3][3];

  tensor[0][0] = itensor[0];
  tensor[1][1] = itensor[1];
  tensor[2][2] = itensor[2];
  tensor[1][2] = tensor[2][1] = itensor[3];
  tensor[0][2] = tensor[2][0] = itensor[4];
  tensor[0][1] = tensor[1][0] = itensor[5];

  if (MathEigen::jacobi3(tensor,inertia,evectors))
    error->all(FLERR,"Insufficient Jacobi rotations for rigid molecule");

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
  max = MAX(inertia[0],inertia[1]);
  max = MAX(max,inertia[2]);

  if (inertia[0] < EPSILON*max) inertia[0] = 0.0;
  if (inertia[1] < EPSILON*max) inertia[1] = 0.0;
  if (inertia[2] < EPSILON*max) inertia[2] = 0.0;

  // enforce 3 evectors as a right-handed coordinate system
  // flip 3rd vector if needed

  MathExtra::cross3(ex,ey,cross);
  if (MathExtra::dot3(cross,ez) < 0.0) MathExtra::negate3(ez);

  // create quaternion

  MathExtra::exyz_to_q(ex,ey,ez,quat);

  // compute displacements in body frame defined by quat

  memory->destroy(dxbody);
  memory->create(dxbody,natoms,3,"molecule:dxbody");

  for (int i = 0; i < natoms; i++)
    MathExtra::transpose_matvec(ex,ey,ez,dxcom[i],dxbody[i]);
}

/* ----------------------------------------------------------------------
   read molecule info from file
   flag = 0, just scan for sizes of fields
   flag = 1, read and store fields
------------------------------------------------------------------------- */

void Molecule::read(int flag)
{
  char line[MAXLINE];
  char *eof;

  // skip 1st line of file

  if (me == 0) {
    eof = fgets(line,MAXLINE,fp);
    if (eof == nullptr) error->one(FLERR,"Unexpected end of molecule file");
  }

  // read header lines
  // skip blank lines or lines that start with "#"
  // stop when read an unrecognized line

  while (1) {

    readline(line);

    // trim comments. if line is blank, continue

    auto text = utils::trim(utils::trim_comment(line));
    if (text.empty()) continue;

    // search line for header keywords and set corresponding variable
    try {
      ValueTokenizer values(text);

      int nmatch = values.count();
      int nwant = 0;
      if (values.contains("atoms")) {
        natoms = values.next_int();
        nwant = 2;
      } else if (values.contains("bonds")) {
        nbonds = values.next_int();
        nwant = 2;
      } else if (values.contains("angles")) {
        nangles = values.next_int();
        nwant = 2;
      } else if (values.contains("dihedrals")) {
        ndihedrals = values.next_int();
        nwant = 2;
      } else if (values.contains("impropers")) {
        nimpropers = values.next_int();
        nwant = 2;
      } else if (values.contains("fragments")) {
        nfragments = values.next_int();
        nwant = 2;
      } else if (values.contains("mass")) {
        massflag = 1;
        masstotal = values.next_double();
        nwant = 2;
        masstotal *= sizescale*sizescale*sizescale;
      } else if (values.contains("com")) {
        comflag = 1;
        com[0] = values.next_double();
        com[1] = values.next_double();
        com[2] = values.next_double();
        nwant = 4;
        com[0] *= sizescale;
        com[1] *= sizescale;
        com[2] *= sizescale;
        if (domain->dimension == 2 && com[2] != 0.0)
          error->all(FLERR,"Molecule file z center-of-mass must be 0.0 for 2d");
      } else if (values.contains("inertia")) {
        inertiaflag = 1;
        itensor[0] = values.next_double();
        itensor[1] = values.next_double();
        itensor[2] = values.next_double();
        itensor[3] = values.next_double();
        itensor[4] = values.next_double();
        itensor[5] = values.next_double();
        nwant = 7;
        const double scale5 = sizescale*sizescale*sizescale*sizescale*sizescale;
        itensor[0] *= scale5;
        itensor[1] *= scale5;
        itensor[2] *= scale5;
        itensor[3] *= scale5;
        itensor[4] *= scale5;
        itensor[5] *= scale5;
      } else if (values.contains("body")) {
        bodyflag = 1;
        avec_body = (AtomVecBody *) atom->style_match("body");
        if (!avec_body)
          error->all(FLERR,"Molecule file requires atom style body");
        nibody = values.next_int();
        ndbody = values.next_int();
        nwant = 3;
      } else {
        // unknown header keyword
        if (utils::strmatch(text,"^\\d+\\s+\\S+")) {
          values.next_int();
          auto keyword = values.next_string();
          error->one(FLERR,"Invalid header keyword: {}",keyword);
        } else break;
      }
      if (nmatch != nwant)
        error->one(FLERR,"Invalid header line format in molecule file");
    } catch (TokenizerException &e) {
      error->one(FLERR, "Invalid header in molecule file\n"
                                    "{}", e.what());
    }
  }

  // error checks

  if (natoms < 1)
    error->all(FLERR,"No atoms or invalid atom count in molecule file");
  if (nbonds < 0) error->all(FLERR,"Invalid bond count in molecule file");
  if (nangles < 0) error->all(FLERR,"Invalid angle count in molecule file");
  if (ndihedrals < 0)
    error->all(FLERR,"Invalid dihedral count in molecule file");
  if (nimpropers < 0)
    error->all(FLERR,"Invalid improper count in molecule file");

  // count = vector for tallying bonds,angles,etc per atom

  if (flag == 0) memory->create(count,natoms,"molecule:count");

  // grab keyword and skip next line

  std::string keyword = parse_keyword(0,line);
  readline(line);

  // loop over sections of molecule file

  while (!keyword.empty()) {
    if (keyword == "Coords") {
      xflag = 1;
      if (flag) coords(line);
      else skip_lines(natoms,line,keyword);
    } else if (keyword == "Types") {
      typeflag = 1;
      if (flag) types(line);
      else skip_lines(natoms,line,keyword);
    } else if (keyword == "Molecules") {
      moleculeflag = 1;
      if (flag) molecules(line);
      else skip_lines(natoms,line,keyword);
    } else if (keyword == "Fragments") {
      if (nfragments == 0)
        error->all(FLERR,"Found Fragments section but no nfragments setting in header");
      fragmentflag = 1;
      if (flag) fragments(line);
      else skip_lines(nfragments,line,keyword);
    } else if (keyword == "Charges") {
      qflag = 1;
      if (flag) charges(line);
      else skip_lines(natoms,line,keyword);
    } else if (keyword == "Diameters") {
      radiusflag = 1;
      if (flag) diameters(line);
      else skip_lines(natoms,line,keyword);
    } else if (keyword == "Masses") {
      rmassflag = 1;
      if (flag) masses(line);
      else skip_lines(natoms,line,keyword);

    } else if (keyword == "Bonds") {
      if (nbonds == 0)
        error->all(FLERR,"Found Bonds section but no nbonds setting in header");
      bondflag = tag_require = 1;
      bonds(flag,line);
    } else if (keyword == "Angles") {
      if (nangles == 0)
        error->all(FLERR,"Found Angles section but no nangles setting in header");
      angleflag = tag_require = 1;
      angles(flag,line);
    } else if (keyword == "Dihedrals") {
      if (ndihedrals == 0) error->all(FLERR,"Found Dihedrals section"
                                      "but no ndihedrals setting in header");
      dihedralflag = tag_require = 1;
      dihedrals(flag,line);
    } else if (keyword == "Impropers") {
      if (nimpropers == 0) error->all(FLERR,"Found Impropers section"
                                      "but no nimpropers setting in header");
      improperflag = tag_require = 1;
      impropers(flag,line);

    } else if (keyword == "Special Bond Counts") {
      nspecialflag = 1;
      nspecial_read(flag,line);
    } else if (keyword == "Special Bonds") {
      specialflag = tag_require = 1;
      if (flag) special_read(line);
      else skip_lines(natoms,line,keyword);

    } else if (keyword == "Shake Flags") {
      shakeflagflag = 1;
      if (flag) shakeflag_read(line);
      else skip_lines(natoms,line,keyword);
    } else if (keyword == "Shake Atoms") {
      shakeatomflag = tag_require = 1;
      if (shaketypeflag) shakeflag = 1;
      if (!shakeflagflag)
        error->all(FLERR,"Shake Flags section must come before Shake Atoms section");
      if (flag) shakeatom_read(line);
      else skip_lines(natoms,line,keyword);
    } else if (keyword == "Shake Bond Types") {
      shaketypeflag = 1;
      if (shakeatomflag) shakeflag = 1;
      if (!shakeflagflag)
        error->all(FLERR,"Shake Flags section must come before Shake Bonds section");
      if (flag) shaketype_read(line);
      else skip_lines(natoms,line,keyword);

    } else if (keyword == "Body Integers") {
      if (bodyflag == 0 || nibody == 0)
        error->all(FLERR,"Found Body Integers section but no setting in header");
      ibodyflag = 1;
      body(flag,0,line);
    } else if (keyword == "Body Doubles") {
      if (bodyflag == 0 || ndbody == 0)
        error->all(FLERR,"Found Body Doubles section but no setting in header");
      dbodyflag = 1;
      body(flag,1,line);
    } else {

      // Error: Either a too long/short section or a typo in the keyword

      if (utils::strmatch(keyword,"^[A-Za-z ]+$"))
        error->one(FLERR,"Unknown section '{}' in molecule "
                                     "file\n",keyword);
      else error->one(FLERR,"Unexpected line in molecule file "
                                        "while looking for the next "
                                        "section:\n{}",line);
    }
    keyword = parse_keyword(1,line);
  }

  // error check

  if (flag == 0) {
    if ((nspecialflag && !specialflag) || (!nspecialflag && specialflag))
      error->all(FLERR,"Molecule file needs both Special Bond sections");
    if (specialflag && !bondflag)
      error->all(FLERR,"Molecule file has special flags but no bonds");
    if ((shakeflagflag || shakeatomflag || shaketypeflag) && !shakeflag)
      error->all(FLERR,"Molecule file shake info is incomplete");
    if (bodyflag && nibody && ibodyflag == 0)
      error->all(FLERR,"Molecule file has no Body Integers section");
    if (bodyflag && ndbody && dbodyflag == 0)
      error->all(FLERR,"Molecule file has no Body Doubles section");
    if (nfragments > 0 && !fragmentflag)
      error->all(FLERR,"Molecule file has no Fragments section");
  }

  // auto-generate special bonds if needed and not in file

  if (bondflag && specialflag == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot auto-generate special bonds before "
                       "simulation box is defined");

    if (flag) {
      special_generate();
      specialflag = 1;
      nspecialflag = 1;
    }
  }

  // body particle must have natom = 1
  // set radius by having body class compute its own radius

  if (bodyflag) {
    radiusflag = 1;
    if (natoms != 1)
      error->all(FLERR,"Molecule natoms must be 1 for body particle");
    if (sizescale != 1.0)
      error->all(FLERR,"Molecule sizescale must be 1.0 for body particle");
    if (flag) {
      radius[0] = avec_body->radius_body(nibody,ndbody,ibodyparams,dbodyparams);
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
        error->all(FLERR,"Invalid line in Coords section of "
                                     "molecule file: {}",line);

      int iatom = values.next_int() - 1;
      if (iatom < 0 || iatom >= natoms)
        error->all(FLERR,"Invalid atom index in Coords section of molecule file");
      count[iatom]++;
      x[iatom][0] = values.next_double();
      x[iatom][1] = values.next_double();
      x[iatom][2] = values.next_double();

      x[iatom][0] *= sizescale;
      x[iatom][1] *= sizescale;
      x[iatom][2] *= sizescale;
    }
  } catch (TokenizerException &e) {
    error->all(FLERR,"Invalid line in Coords section of "
                                 "molecule file: {}\n{}",e.what(),line);
  }

  for (int i = 0; i < natoms; i++)
    if (count[i] == 0) error->all(FLERR,"Atom {} missing in Coords "
                                                    "section of molecule file",i+1);

  if (domain->dimension == 2) {
    for (int i = 0; i < natoms; i++)
      if (x[i][2] != 0.0)
        error->all(FLERR,"Z coord in molecule file for atom {} "
                                     "must be 0.0 for 2d-simulation.",i+1);
  }
}

/* ----------------------------------------------------------------------
   read types from file
   set ntypes = max of any atom type
------------------------------------------------------------------------- */

void Molecule::types(char *line)
{
  for (int i = 0; i < natoms; i++) count[i] = 0;
  try {
    for (int i = 0; i < natoms; i++) {
      readline(line);

      ValueTokenizer values(utils::trim_comment(line));
      if (values.count() != 2)
        error->all(FLERR,"Invalid line in Types section of "
                                     "molecule file: {}",line);

      int iatom = values.next_int() - 1;
      if (iatom < 0 || iatom >= natoms)
        error->all(FLERR,"Invalid atom index in Types section of molecule file");
      count[iatom]++;
      type[iatom] = values.next_int();
      type[iatom] += toffset;
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, "Invalid line in Types section of "
                                  "molecule file: {}\n{}", e.what(),line);
  }

  for (int i = 0; i < natoms; i++) {
    if (count[i] == 0) error->all(FLERR,"Atom {} missing in Types "
                                                    "section of molecule file",i+1);

    if ((type[i] <= 0) || (domain->box_exist && (type[i] > atom->ntypes)))
      error->all(FLERR,"Invalid atom type {} for atom {} "
                                   "in molecule file",type[i],i+1);

    ntypes = MAX(ntypes,type[i]);
  }
}

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
        error->all(FLERR,"Invalid line in Molecules section of "
                                     "molecule file: {}",line);

      int iatom = values.next_int() - 1;
      if (iatom < 0 || iatom >= natoms)
        error->all(FLERR,"Invalid atom index in Molecules section of molecule file");
      count[iatom]++;
      molecule[iatom] = values.next_tagint();
      // molecule[iatom] += moffset; // placeholder for possible molecule offset
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, "Invalid line in Molecules section of "
                                  "molecule file: {}\n{}",e.what(),line);
  }

  for (int i = 0; i < natoms; i++)
    if (count[i] == 0) error->all(FLERR,"Atom {} missing in Molecules "
                                                    "section of molecule file",i+1);

  for (int i = 0; i < natoms; i++)
    if (molecule[i] < 0)
      error->all(FLERR,"Invalid molecule ID {} for atom {} "
                                   "in molecule file",molecule[i],i+1);

  for (int i = 0; i < natoms; i++)
    nmolecules = MAX(nmolecules,molecule[i]);
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

      if ((int)values.count() > natoms+1)
        error->all(FLERR,"Too many atoms per fragment in Fragments "
                   "section of molecule file");

      fragmentnames[i] = values.next_string();

      while (values.has_next()) {
        int iatom = values.next_int()-1;
        if (iatom < 0 || iatom >= natoms)
          error->all(FLERR,"Invalid atom ID {} for fragment {} in "
                                       "Fragments section of molecule file",
                                       iatom+1, fragmentnames[i]);
        fragmentmask[i][iatom] = 1;
      }
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, "Invalid atom ID in Fragments section of "
                                  "molecule file: {}\n{}", e.what(),line);
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
      if ((int)values.count() != 2)
        error->all(FLERR,"Invalid line in Charges section of "
                                     "molecule file: {}",line);

      int iatom = values.next_int() - 1;
      if (iatom < 0 || iatom >= natoms)
        error->all(FLERR,"Invalid atom index in Charges section of molecule file");

      count[iatom]++;
      q[iatom] = values.next_double();
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, "Invalid line in Charges section of "
                                  "molecule file: {}.\n{}",e.what(),line);
  }

  for (int i = 0; i < natoms; i++)
    if (count[i] == 0)  error->all(FLERR,"Atom {} missing in Charges "
                                                    "section of molecule file",i+1);
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
        error->all(FLERR,"Invalid line in Diameters section of "
                                     "molecule file: {}",line);
      int iatom = values.next_int() - 1;
      if (iatom < 0 || iatom >= natoms)
        error->all(FLERR,"Invalid atom index in Diameters section of molecule file");
      count[iatom]++;
      radius[iatom] = values.next_double();
      radius[iatom] *= sizescale;
      radius[iatom] *= 0.5;
      maxradius = MAX(maxradius,radius[iatom]);
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, "Invalid line in Diameters section of "
                                  "molecule file: {}\n{}",e.what(),line);
  }

  for (int i = 0; i < natoms; i++) {
    if (count[i] == 0)  error->all(FLERR,"Atom {} missing in Diameters "
                                                    "section of molecule file",i+1);
    if (radius[i] < 0.0)
      error->all(FLERR,"Invalid atom diameter {} for atom {} "
                                   "in molecule file", radius[i], i+1);
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
        error->all(FLERR,"Invalid line in Masses section of "
                                     "molecule file: {}",line);

      int iatom = values.next_int() - 1;
      if (iatom < 0 || iatom >= natoms)
        error->all(FLERR,"Invalid atom index in Masses section of molecule file");
      count[iatom]++;
      rmass[iatom] = values.next_double();
      rmass[iatom] *= sizescale*sizescale*sizescale;
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, "Invalid line in Masses section of "
                                  "molecule file: {}\n{}",e.what(),line);
  }

  for (int i = 0; i < natoms; i++) {
    if (count[i] == 0)   error->all(FLERR,"Atom {} missing in Masses "
                                                    "section of molecule file",i+1);
    if (rmass[i] <= 0.0)
      error->all(FLERR,"Invalid atom mass {} for atom {} "
                                   "in molecule file", radius[i], i+1);
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
  int itype;
  tagint m,atom1,atom2;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_bond[i] = 0;

  for (int i = 0; i < nbonds; i++) {
    readline(line);

    try {
      ValueTokenizer values(utils::trim_comment(line));
      if (values.count() != 4)
        error->all(FLERR,"Invalid line in Bonds section of "
                                     "molecule file: {}",line);
      values.next_int();
      itype = values.next_int();
      atom1 = values.next_tagint();
      atom2 = values.next_tagint();
    } catch (TokenizerException &e) {
      error->all(FLERR, "Invalid line in Bonds section of "
                                  "molecule file: {}\n{}",e.what(),line);
    }

    itype += boffset;

    if ((atom1 <= 0) || (atom1 > natoms) ||
        (atom2 <= 0) || (atom2 > natoms) || (atom1 == atom2))
      error->all(FLERR,"Invalid atom ID in Bonds section of molecule file");
    if ((itype <= 0) || (domain->box_exist && (itype > atom->nbondtypes)))
      error->all(FLERR,"Invalid bond type in Bonds section of molecule file");

    if (flag) {
      m = atom1-1;
      nbondtypes = MAX(nbondtypes,itype);
      bond_type[m][num_bond[m]] = itype;
      bond_atom[m][num_bond[m]] = atom2;
      num_bond[m]++;
      if (newton_bond == 0) {
        m = atom2-1;
        bond_type[m][num_bond[m]] = itype;
        bond_atom[m][num_bond[m]] = atom1;
        num_bond[m]++;
      }
    } else {
      count[atom1-1]++;
      if (newton_bond == 0) count[atom2-1]++;
    }
  }

  // bond_per_atom = max of count vector

  if (flag == 0) {
    bond_per_atom = 0;
    for (int i = 0; i < natoms; i++)
      bond_per_atom = MAX(bond_per_atom,count[i]);
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
  int itype;
  tagint m,atom1,atom2,atom3;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_angle[i] = 0;

  for (int i = 0; i < nangles; i++) {
    readline(line);

    try {
      ValueTokenizer values(utils::trim_comment(line));
      if (values.count() != 5)
        error->all(FLERR,"Invalid line in Angles section of "
                                     "molecule file: {}",line);
      values.next_int();
      itype = values.next_int();
      atom1 = values.next_tagint();
      atom2 = values.next_tagint();
      atom3 = values.next_tagint();
    } catch (TokenizerException &e) {
      error->all(FLERR, "Invalid line in Angles section of "
                                    "molecule file: {}\n{}",e.what(),line);
    }

    itype += aoffset;

    if ((atom1 <= 0) || (atom1 > natoms) ||
        (atom2 <= 0) || (atom2 > natoms) ||
        (atom3 <= 0) || (atom3 > natoms) ||
        (atom1 == atom2) || (atom1 == atom3) || (atom2 == atom3))
      error->all(FLERR,"Invalid atom ID in Angles section of molecule file");
    if ((itype <= 0) || (domain->box_exist && (itype > atom->nangletypes)))
      error->all(FLERR,"Invalid angle type in Angles section of molecule file");

    if (flag) {
      m = atom2-1;
      nangletypes = MAX(nangletypes,itype);
      angle_type[m][num_angle[m]] = itype;
      angle_atom1[m][num_angle[m]] = atom1;
      angle_atom2[m][num_angle[m]] = atom2;
      angle_atom3[m][num_angle[m]] = atom3;
      num_angle[m]++;
      if (newton_bond == 0) {
        m = atom1-1;
        angle_type[m][num_angle[m]] = itype;
        angle_atom1[m][num_angle[m]] = atom1;
        angle_atom2[m][num_angle[m]] = atom2;
        angle_atom3[m][num_angle[m]] = atom3;
        num_angle[m]++;
        m = atom3-1;
        angle_type[m][num_angle[m]] = itype;
        angle_atom1[m][num_angle[m]] = atom1;
        angle_atom2[m][num_angle[m]] = atom2;
        angle_atom3[m][num_angle[m]] = atom3;
        num_angle[m]++;
      }
    } else {
      count[atom2-1]++;
      if (newton_bond == 0) {
        count[atom1-1]++;
        count[atom3-1]++;
      }
    }
  }

  // angle_per_atom = max of count vector

  if (flag == 0) {
    angle_per_atom = 0;
    for (int i = 0; i < natoms; i++)
      angle_per_atom = MAX(angle_per_atom,count[i]);
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
  int itype;
  tagint m,atom1,atom2,atom3,atom4;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_dihedral[i] = 0;

  for (int i = 0; i < ndihedrals; i++) {
    readline(line);

    try {
      ValueTokenizer values(utils::trim_comment(line));
      if (values.count() != 6)
        error->all(FLERR,"Invalid line in Dihedrals section of "
                                     "molecule file: {}",line);

      values.next_int();
      itype = values.next_int();
      atom1 = values.next_tagint();
      atom2 = values.next_tagint();
      atom3 = values.next_tagint();
      atom4 = values.next_tagint();
    } catch (TokenizerException &e) {
      error->all(FLERR, "Invalid line in Dihedrals section of "
                                    "molecule file: {}\n{}",e.what(),line);
    }

    itype += doffset;

    if ((atom1 <= 0) || (atom1 > natoms) ||
        (atom2 <= 0) || (atom2 > natoms) ||
        (atom3 <= 0) || (atom3 > natoms) ||
        (atom4 <= 0) || (atom4 > natoms) ||
        (atom1 == atom2) || (atom1 == atom3) || (atom1 == atom4) ||
        (atom2 == atom3) || (atom2 == atom4) || (atom3 == atom4))
      error->all(FLERR,
                 "Invalid atom ID in dihedrals section of molecule file");
    if ((itype <= 0) || (domain->box_exist && (itype > atom->ndihedraltypes)))
      error->all(FLERR,"Invalid dihedral type in Dihedrals section of molecule file");

    if (flag) {
      m = atom2-1;
      ndihedraltypes = MAX(ndihedraltypes,itype);
      dihedral_type[m][num_dihedral[m]] = itype;
      dihedral_atom1[m][num_dihedral[m]] = atom1;
      dihedral_atom2[m][num_dihedral[m]] = atom2;
      dihedral_atom3[m][num_dihedral[m]] = atom3;
      dihedral_atom4[m][num_dihedral[m]] = atom4;
      num_dihedral[m]++;
      if (newton_bond == 0) {
        m = atom1-1;
        dihedral_type[m][num_dihedral[m]] = itype;
        dihedral_atom1[m][num_dihedral[m]] = atom1;
        dihedral_atom2[m][num_dihedral[m]] = atom2;
        dihedral_atom3[m][num_dihedral[m]] = atom3;
        dihedral_atom4[m][num_dihedral[m]] = atom4;
        num_dihedral[m]++;
        m = atom3-1;
        dihedral_type[m][num_dihedral[m]] = itype;
        dihedral_atom1[m][num_dihedral[m]] = atom1;
        dihedral_atom2[m][num_dihedral[m]] = atom2;
        dihedral_atom3[m][num_dihedral[m]] = atom3;
        dihedral_atom4[m][num_dihedral[m]] = atom4;
        num_dihedral[m]++;
        m = atom4-1;
        dihedral_type[m][num_dihedral[m]] = itype;
        dihedral_atom1[m][num_dihedral[m]] = atom1;
        dihedral_atom2[m][num_dihedral[m]] = atom2;
        dihedral_atom3[m][num_dihedral[m]] = atom3;
        dihedral_atom4[m][num_dihedral[m]] = atom4;
        num_dihedral[m]++;
      }
    } else {
      count[atom2-1]++;
      if (newton_bond == 0) {
        count[atom1-1]++;
        count[atom3-1]++;
        count[atom4-1]++;
      }
    }
  }

  // dihedral_per_atom = max of count vector

  if (flag == 0) {
    dihedral_per_atom = 0;
    for (int i = 0; i < natoms; i++)
      dihedral_per_atom = MAX(dihedral_per_atom,count[i]);
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
  int itype;
  tagint m,atom1,atom2,atom3,atom4;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_improper[i] = 0;

  for (int i = 0; i < nimpropers; i++) {
    readline(line);

    try {
      ValueTokenizer values(utils::trim_comment(line));
      if (values.count() != 6)
                error->all(FLERR,"Invalid line in Impropers section of "
                                     "molecule file: {}",line);
      values.next_int();
      itype = values.next_int();
      atom1 = values.next_tagint();
      atom2 = values.next_tagint();
      atom3 = values.next_tagint();
      atom4 = values.next_tagint();
    } catch (TokenizerException &e) {
      error->all(FLERR, "Invalid line in Impropers section of "
                                    "molecule file: {}\n{}",e.what(),line);
    }

    itype += ioffset;

    if ((atom1 <= 0) || (atom1 > natoms) ||
        (atom2 <= 0) || (atom2 > natoms) ||
        (atom3 <= 0) || (atom3 > natoms) ||
        (atom4 <= 0) || (atom4 > natoms) ||
        (atom1 == atom2) || (atom1 == atom3) || (atom1 == atom4) ||
        (atom2 == atom3) || (atom2 == atom4) || (atom3 == atom4))
      error->all(FLERR,
                 "Invalid atom ID in impropers section of molecule file");
    if ((itype <= 0) || (domain->box_exist && (itype > atom->nimpropertypes)))
      error->all(FLERR,"Invalid improper type in Impropers section of molecule file");

    if (flag) {
      m = atom2-1;
      nimpropertypes = MAX(nimpropertypes,itype);
      improper_type[m][num_improper[m]] = itype;
      improper_atom1[m][num_improper[m]] = atom1;
      improper_atom2[m][num_improper[m]] = atom2;
      improper_atom3[m][num_improper[m]] = atom3;
      improper_atom4[m][num_improper[m]] = atom4;
      num_improper[m]++;
      if (newton_bond == 0) {
        m = atom1-1;
        improper_type[m][num_improper[m]] = itype;
        improper_atom1[m][num_improper[m]] = atom1;
        improper_atom2[m][num_improper[m]] = atom2;
        improper_atom3[m][num_improper[m]] = atom3;
        improper_atom4[m][num_improper[m]] = atom4;
        num_improper[m]++;
        m = atom3-1;
        improper_type[m][num_improper[m]] = itype;
        improper_atom1[m][num_improper[m]] = atom1;
        improper_atom2[m][num_improper[m]] = atom2;
        improper_atom3[m][num_improper[m]] = atom3;
        improper_atom4[m][num_improper[m]] = atom4;
        num_improper[m]++;
        m = atom4-1;
        improper_type[m][num_improper[m]] = itype;
        improper_atom1[m][num_improper[m]] = atom1;
        improper_atom2[m][num_improper[m]] = atom2;
        improper_atom3[m][num_improper[m]] = atom3;
        improper_atom4[m][num_improper[m]] = atom4;
        num_improper[m]++;
      }
    } else {
      count[atom2-1]++;
      if (newton_bond == 0) {
        count[atom1-1]++;
        count[atom3-1]++;
        count[atom4-1]++;
      }
    }
  }

  // improper_per_atom = max of count vector

  if (flag == 0) {
    improper_per_atom = 0;
    for (int i = 0; i < natoms; i++)
      improper_per_atom = MAX(improper_per_atom,count[i]);
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

  for (int i = 0; i < natoms; i++) {
    readline(line);

    int c1, c2, c3;

    try {
      ValueTokenizer values(utils::trim_comment(line));
      if (values.count() != 4)
        error->all(FLERR,"Invalid line in Special Bond Counts section of "
                                     "molecule file: {}",line);
      values.next_int();
      c1 = values.next_tagint();
      c2 = values.next_tagint();
      c3 = values.next_tagint();
    } catch (TokenizerException &e) {
      error->all(FLERR, "Invalid line in Special Bond Counts section of "
                                    "molecule file: {}\n{}",e.what(),line);
    }

    if (flag) {
      nspecial[i][0] = c1;
      nspecial[i][1] = c1+c2;
      nspecial[i][2] = c1+c2+c3;
    } else maxspecial = MAX(maxspecial,c1+c2+c3);
  }
}

/* ----------------------------------------------------------------------
   read special bond indices from file
------------------------------------------------------------------------- */

void Molecule::special_read(char *line)
{
  try {
    for (int i = 0; i < natoms; i++) {
      readline(line);

      ValueTokenizer values(utils::trim_comment(line));
      int nwords = values.count();

      if (nwords != nspecial[i][2]+1)
        error->all(FLERR,"Molecule file special list "
                        "does not match special count");

      values.next_int(); // ignore

      for (int m = 1; m < nwords; m++) {
        special[i][m-1] = values.next_tagint();
        if (special[i][m-1] <= 0 || special[i][m-1] > natoms ||
            special[i][m-1] == i+1)
          error->all(FLERR,"Invalid atom index in Special Bonds "
                     "section of molecule file");
      }
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, "Invalid line in Special Bonds section of "
                                  "molecule file: {}\n{}",e.what(),line);
  }
}

/* ----------------------------------------------------------------------
   auto generate special bond info
------------------------------------------------------------------------- */

void Molecule::special_generate()
{
  int newton_bond = force->newton_bond;
  tagint atom1,atom2;

  // temporary array for special atoms

  tagint **tmpspecial;
  memory->create(tmpspecial,natoms,atom->maxspecial,"molecule:tmpspecial");
  memset(&tmpspecial[0][0],0,sizeof(tagint)*natoms*atom->maxspecial);

  for (int i = 0; i < natoms; i++) count[i] = 0;

  // 1-2 neighbors

  if (newton_bond) {
    for (int i = 0; i < natoms; i++) {
      for (int j = 0; j < num_bond[i]; j++) {
        atom1 = i;
        atom2 = bond_atom[i][j]-1;
        nspecial[i][0]++;
        nspecial[atom2][0]++;
        if (count[i] >= atom->maxspecial || count[atom2] >= atom->maxspecial)
          error->all(FLERR,"Molecule auto special bond generation overflow");
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
          error->all(FLERR,"Molecule auto special bond generation overflow");
        tmpspecial[i][count[atom1]++] = atom2;
      }
    }
  }

  // 1-3 neighbors with no duplicates

  for (int i = 0; i < natoms; i++) nspecial[i][1] = nspecial[i][0];

  int dedup;
  for (int i = 0; i < natoms; i++) {
    for (int m = 0; m < nspecial[i][0]; m++) {
      for (int j = 0; j < nspecial[tmpspecial[i][m]-1][0]; j++) {
        dedup = 0;
        for (int k =0; k < count[i]; k++) {
          if (tmpspecial[tmpspecial[i][m]-1][j] == tmpspecial[i][k] ||
              tmpspecial[tmpspecial[i][m]-1][j] == i+1) {
            dedup = 1;
          }
        }
        if (!dedup) {
          if (count[i] >= atom->maxspecial)
            error->all(FLERR,"Molecule auto special bond generation overflow");
          tmpspecial[i][count[i]++] = tmpspecial[tmpspecial[i][m]-1][j];
          nspecial[i][1]++;
        }
      }
    }
  }

  // 1-4 neighbors with no duplicates

  for (int i = 0; i < natoms; i++) nspecial[i][2] = nspecial[i][1];

  for (int i = 0; i < natoms; i++) {
    for (int m = nspecial[i][0]; m < nspecial[i][1]; m++) {
      for (int j = 0; j < nspecial[tmpspecial[i][m]-1][0]; j++) {
        dedup = 0;
        for (int k =0; k < count[i]; k++) {
          if (tmpspecial[tmpspecial[i][m]-1][j] == tmpspecial[i][k] ||
              tmpspecial[tmpspecial[i][m]-1][j] == i+1) {
            dedup = 1;
          }
        }
        if (!dedup) {
          if (count[i] >= atom->maxspecial)
            error->all(FLERR,"Molecule auto special bond generation overflow");
          tmpspecial[i][count[i]++] = tmpspecial[tmpspecial[i][m]-1][j];
          nspecial[i][2]++;
        }
      }
    }
  }

  maxspecial = 0;
  for (int i = 0; i < natoms; i++)
    maxspecial = MAX(maxspecial,nspecial[i][2]);

  memory->create(special,natoms,maxspecial,"molecule:special");
  for (int i = 0; i < natoms; i++)
    for (int j = 0; j < nspecial[i][2]; j++)
      special[i][j] = tmpspecial[i][j];

  memory->destroy(tmpspecial);
}

/* ----------------------------------------------------------------------
   read SHAKE flags from file
------------------------------------------------------------------------- */

void Molecule::shakeflag_read(char *line)
{
  try {
    for (int i = 0; i < natoms; i++) {
      readline(line);

      ValueTokenizer values(utils::trim_comment(line));

      if (values.count() != 2)
        error->all(FLERR,"Invalid Shake Flags section in molecule file");

      values.next_int();
      shake_flag[i] = values.next_int();
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, "Invalid Shake Flags section in molecule file\n"
                                  "{}", e.what());
  }

  for (int i = 0; i < natoms; i++)
    if (shake_flag[i] < 0 || shake_flag[i] > 4)
      error->all(FLERR,"Invalid shake flag in molecule file");
}

/* ----------------------------------------------------------------------
   read SHAKE atom info from file
------------------------------------------------------------------------- */

void Molecule::shakeatom_read(char *line)
{
  int nmatch=0, nwant=0;
  try {
    for (int i = 0; i < natoms; i++) {
      readline(line);

      ValueTokenizer values(utils::trim_comment(line));
      nmatch = values.count();

      switch (shake_flag[i]) {
        case 1:
          values.next_int();
          shake_atom[i][0] = values.next_tagint();
          shake_atom[i][1] = values.next_tagint();
          shake_atom[i][2] = values.next_tagint();
          nwant = 4;
          break;

        case 2:
          values.next_int();
          shake_atom[i][0] = values.next_tagint();
          shake_atom[i][1] = values.next_tagint();
          nwant = 3;
          break;

        case 3:
          values.next_int();
          shake_atom[i][0] = values.next_tagint();
          shake_atom[i][1] = values.next_tagint();
          shake_atom[i][2] = values.next_tagint();
          nwant = 4;
          break;

        case 4:
          values.next_int();
          shake_atom[i][0] = values.next_tagint();
          shake_atom[i][1] = values.next_tagint();
          shake_atom[i][2] = values.next_tagint();
          shake_atom[i][3] = values.next_tagint();
          nwant = 5;
          break;

        case 0:
          values.next_int();
          nwant = 1;
          break;

        default:
          error->all(FLERR,"Invalid shake atom in molecule file");
      }

      if (nmatch != nwant)
        error->all(FLERR,"Invalid shake atom in molecule file");
    }

  } catch (TokenizerException &e) {
    error->all(FLERR,"Invalid shake atom in molecule file\n"
                                 "{}", e.what());
  }

  for (int i = 0; i < natoms; i++) {
    int m = shake_flag[i];
    if (m == 1) m = 3;
    for (int j = 0; j < m; j++)
      if (shake_atom[i][j] <= 0 || shake_atom[i][j] > natoms)
        error->all(FLERR,"Invalid shake atom in molecule file");
  }
}

/* ----------------------------------------------------------------------
   read SHAKE bond type info from file
------------------------------------------------------------------------- */

void Molecule::shaketype_read(char *line)
{
  try {
    int nmatch=0, nwant=0;
    for (int i = 0; i < natoms; i++) {
      readline(line);

      ValueTokenizer values(utils::trim_comment(line));
      nmatch = values.count();

      switch (shake_flag[i]) {
        case 1:
          values.next_int();
          shake_type[i][0] = values.next_int();
          shake_type[i][1] = values.next_int();
          shake_type[i][2] = values.next_int();
          nwant = 4;
          break;

        case 2:
          values.next_int();
          shake_type[i][0] = values.next_int();
          nwant = 2;
          break;

        case 3:
          values.next_int();
          shake_type[i][0] = values.next_int();
          shake_type[i][1] = values.next_int();
          nwant = 3;
          break;

        case 4:
          values.next_int();
          shake_type[i][0] = values.next_int();
          shake_type[i][1] = values.next_int();
          shake_type[i][2] = values.next_int();
          nwant = 4;
          break;

        case 0:
          values.next_int();
          nwant = 1;
          break;

        default:
          error->all(FLERR,"Invalid shake type data in molecule file");
      }

      if (nmatch != nwant)
        error->all(FLERR,"Invalid shake type data in molecule file");
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, "Invalid shake type data in molecule file\n",
                                  "{}", e.what());
  }

  for (int i = 0; i < natoms; i++) {
    int m = shake_flag[i];
    if (m == 1) m = 3;
    for (int j = 0; j < m-1; j++)
      if (shake_type[i][j] <= 0)
        error->all(FLERR,"Invalid shake bond type in molecule file");
    if (shake_flag[i] == 1)
      if (shake_type[i][2] <= 0)
        error->all(FLERR,"Invalid shake angle type in molecule file");
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

      if (ncount == 0)
        error->all(FLERR,"Too few values in body section of molecule file");
      if (nword+ncount > nparam)
        error->all(FLERR,"Too many values in body section of molecule file");

      if (flag) {
        if (pflag == 0) {
          while (values.has_next()) {
            ibodyparams[nword++] = values.next_int();
          }
        } else {
          while (values.has_next()) {
            dbodyparams[nword++] = values.next_double();
          }
        }
      } else nword += ncount;
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, "Invalid body params in molecule file\n",
                                  "{}", e.what());
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
   flag = 0, just check this molecule
   flag = 1, check all molecules in set, this is 1st molecule in set
------------------------------------------------------------------------- */

void Molecule::check_attributes(int flag)
{
  int n = 1;
  if (flag) n = nset;
  int imol = atom->find_molecule(id);

  for (int i = imol; i < imol+n; i++) {
    Molecule *onemol = atom->molecules[imol];

    // check per-atom attributes of molecule
    // warn if not a match

    int mismatch = 0;
    if (onemol->qflag && !atom->q_flag) mismatch = 1;
    if (onemol->radiusflag && !atom->radius_flag) mismatch = 1;
    if (onemol->rmassflag && !atom->rmass_flag) mismatch = 1;

    if (mismatch && me == 0)
      error->warning(FLERR,
                     "Molecule attributes do not match system attributes");

    // for all atom styles, check nbondtype,etc

    mismatch = 0;
    if (atom->nbondtypes < onemol->nbondtypes) mismatch = 1;
    if (atom->nangletypes < onemol->nangletypes) mismatch = 1;
    if (atom->ndihedraltypes < onemol->ndihedraltypes) mismatch = 1;
    if (atom->nimpropertypes < onemol->nimpropertypes) mismatch = 1;

    if (mismatch)
      error->all(FLERR,"Molecule topology type exceeds system topology type");

    // for molecular atom styles, check bond_per_atom,etc + maxspecial
    // do not check for atom style template, since nothing stored per atom

    if (atom->molecular == Atom::MOLECULAR) {
      if (atom->avec->bonds_allow &&
          atom->bond_per_atom < onemol->bond_per_atom) mismatch = 1;
      if (atom->avec->angles_allow &&
          atom->angle_per_atom < onemol->angle_per_atom) mismatch = 1;
      if (atom->avec->dihedrals_allow &&
          atom->dihedral_per_atom < onemol->dihedral_per_atom) mismatch = 1;
      if (atom->avec->impropers_allow &&
          atom->improper_per_atom < onemol->improper_per_atom) mismatch = 1;
      if (atom->maxspecial < onemol->maxspecial) mismatch = 1;

      if (mismatch)
        error->all(FLERR,"Molecule topology/atom exceeds system topology/atom");
    }

    // warn if molecule topology defined but no special settings

    if (onemol->bondflag && !onemol->specialflag)
      if (me == 0) error->warning(FLERR,"Molecule has bond topology "
                                  "but no special bond settings");
  }
}

/* ----------------------------------------------------------------------
   init all data structures to empty
------------------------------------------------------------------------- */

void Molecule::initialize()
{
  natoms = 0;
  nbonds = nangles = ndihedrals = nimpropers = 0;
  ntypes = 0;
  nmolecules = 1;
  nbondtypes = nangletypes = ndihedraltypes = nimpropertypes = 0;
  nibody = ndbody = 0;
  nfragments = 0;

  bond_per_atom = angle_per_atom = dihedral_per_atom = improper_per_atom = 0;
  maxspecial = 0;

  xflag = typeflag = moleculeflag = fragmentflag = qflag = radiusflag = rmassflag = 0;
  bondflag = angleflag = dihedralflag = improperflag = 0;
  nspecialflag = specialflag = 0;
  shakeflag = shakeflagflag = shakeatomflag = shaketypeflag = 0;
  bodyflag = ibodyflag = dbodyflag = 0;

  centerflag = massflag = comflag = inertiaflag = 0;
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
  if (xflag) memory->create(x,natoms,3,"molecule:x");
  if (typeflag) memory->create(type,natoms,"molecule:type");
  if (moleculeflag) memory->create(molecule,natoms,"molecule:molecule");
  if (fragmentflag) {
     fragmentnames.resize(nfragments);
     memory->create(fragmentmask,nfragments,natoms,"molecule:fragmentmask");
     for (int i = 0; i < nfragments; i++)
       for (int j = 0; j < natoms; j++) fragmentmask[i][j] = 0;
  }
  if (qflag) memory->create(q,natoms,"molecule:q");
  if (radiusflag) memory->create(radius,natoms,"molecule:radius");
  if (rmassflag) memory->create(rmass,natoms,"molecule:rmass");

  // always allocate num_bond,num_angle,etc and nspecial
  // even if not in molecule file, initialize to 0
  // this is so methods that use these arrays don't have to check they exist

  memory->create(num_bond,natoms,"molecule:num_bond");
  for (int i = 0; i < natoms; i++) num_bond[i] = 0;
  memory->create(num_angle,natoms,"molecule:num_angle");
  for (int i = 0; i < natoms; i++) num_angle[i] = 0;
  memory->create(num_dihedral,natoms,"molecule:num_dihedral");
  for (int i = 0; i < natoms; i++) num_dihedral[i] = 0;
  memory->create(num_improper,natoms,"molecule:num_improper");
  for (int i = 0; i < natoms; i++) num_improper[i] = 0;
  memory->create(nspecial,natoms,3,"molecule:nspecial");
  for (int i = 0; i < natoms; i++)
    nspecial[i][0] = nspecial[i][1] = nspecial[i][2] = 0;

  if (specialflag)
    memory->create(special,natoms,maxspecial,"molecule:special");

  if (bondflag) {
    memory->create(bond_type,natoms,bond_per_atom,
                   "molecule:bond_type");
    memory->create(bond_atom,natoms,bond_per_atom,
                   "molecule:bond_atom");
  }

  if (angleflag) {
    memory->create(angle_type,natoms,angle_per_atom,
                   "molecule:angle_type");
    memory->create(angle_atom1,natoms,angle_per_atom,
                   "molecule:angle_atom1");
    memory->create(angle_atom2,natoms,angle_per_atom,
                   "molecule:angle_atom2");
    memory->create(angle_atom3,natoms,angle_per_atom,
                   "molecule:angle_atom3");
  }

  if (dihedralflag) {
    memory->create(dihedral_type,natoms,dihedral_per_atom,
                   "molecule:dihedral_type");
    memory->create(dihedral_atom1,natoms,dihedral_per_atom,
                   "molecule:dihedral_atom1");
    memory->create(dihedral_atom2,natoms,dihedral_per_atom,
                   "molecule:dihedral_atom2");
    memory->create(dihedral_atom3,natoms,dihedral_per_atom,
                   "molecule:dihedral_atom3");
    memory->create(dihedral_atom4,natoms,dihedral_per_atom,
                   "molecule:dihedral_atom4");
  }

  if (improperflag) {
    memory->create(improper_type,natoms,improper_per_atom,
                   "molecule:improper_type");
    memory->create(improper_atom1,natoms,improper_per_atom,
                   "molecule:improper_atom1");
    memory->create(improper_atom2,natoms,improper_per_atom,
                   "molecule:improper_atom2");
    memory->create(improper_atom3,natoms,improper_per_atom,
                   "molecule:improper_atom3");
    memory->create(improper_atom4,natoms,improper_per_atom,
                   "molecule:improper_atom4");
  }

  if (shakeflag) {
    memory->create(shake_flag,natoms,"molecule:shake_flag");
    memory->create(shake_atom,natoms,4,"molecule:shake_flag");
    memory->create(shake_type,natoms,3,"molecule:shake_flag");
  }

  if (bodyflag) {
    if (nibody) memory->create(ibodyparams,nibody,"molecule:ibodyparams");
    if (ndbody) memory->create(dbodyparams,ndbody,"molecule:dbodyparams");
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
  memory->destroy(radius);
  memory->destroy(rmass);

  memory->destroy(molecule);
  memory->destroy(fragmentmask);

  if (fragmentflag) {
    fragmentnames.clear();
  }

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
  if (me == 0) {
    if (fgets(line,MAXLINE,fp) == nullptr) n = 0;
    else n = strlen(line) + 1;
  }
  MPI_Bcast(&n,1,MPI_INT,0,world);
  if (n == 0) error->all(FLERR,"Unexpected end of molecule file");
  MPI_Bcast(line,n,MPI_CHAR,0,world);
}

/* ----------------------------------------------------------------------
   extract keyword from line
   flag = 0, read and bcast line
   flag = 1, line has already been read
------------------------------------------------------------------------- */

std::string Molecule::parse_keyword(int flag, char *line)
{
  char line2[MAXLINE];
  if (flag) {

    // read upto non-blank line plus 1 following line
    // eof is set to 1 if any read hits end-of-file

    int eof = 0;
    if (me == 0) {
      if (fgets(line,MAXLINE,fp) == nullptr) eof = 1;
      while (eof == 0 && strspn(line," \t\n\r") == strlen(line)) {
        if (fgets(line,MAXLINE,fp) == nullptr) eof = 1;
      }
      if (fgets(line2,MAXLINE,fp) == nullptr) eof = 1;
    }

    // if eof, set keyword empty and return

    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) {
      return std::string("");
    }

    // bcast keyword line to all procs

    MPI_Bcast(line,MAXLINE,MPI_CHAR,0,world);
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
    if (utils::strmatch(utils::trim(utils::trim_comment(line)),"^[A-Za-z ]+$"))
      error->one(FLERR,"Unexpected line in molecule file while "
                                   "skipping {} section:\n{}",section,line);
  }
}

/* ----------------------------------------------------------------------
   proc 0 prints molecule params
------------------------------------------------------------------------- */

/*

void Molecule::print()
{
  printf("MOLECULE %s\n",id);
  printf("  %d natoms\n",natoms);
  if (nbonds) printf("  %d nbonds\n",nbonds);
  if (nangles) printf("  %d nangles\n",nangles);
  if (ndihedrals) printf("  %d ndihedrals\n",ndihedrals);
  if (nimpropers) printf("  %d nimpropers\n",nimpropers);

  if (xflag) {
    printf(  "Coords:\n");
    for (int i = 0; i < natoms; i++)
      printf("    %d %g %g %g\n",i+1,x[i][0],x[i][1],x[i][2]);
  }
  if (typeflag) {
    printf(  "Types:\n");
    for (int i = 0; i < natoms; i++)
      printf("    %d %d\n",i+1,type[i]);
  }
  if (qflag) {
    printf(  "Charges:\n");
    for (int i = 0; i < natoms; i++)
      printf("    %d %g\n",i+1,q[i]);
  }
  if (radiusflag) {
    printf(  "Radii:\n");
    for (int i = 0; i < natoms; i++)
      printf("    %d %g\n",i+1,radius[i]);
  }
  if (rmassflag) {
    printf(  "Masses:\n");
    for (int i = 0; i < natoms; i++)
      printf("    %d %g\n",i+1,rmass[i]);
  }

  if (bondflag) {
    printf(  "Bonds:\n");
    for (int i = 0; i < natoms; i++) {
      printf("    %d %d\n",i+1,num_bond[i]);
      for (int j = 0; j < num_bond[i]; j++)
        printf("      %d %d %d %d\n",j+1,bond_type[i][j],i+1,bond_atom[i][j]);
    }
  }
  if (angleflag) {
    printf(  "Angles:\n");
    for (int i = 0; i < natoms; i++) {
      printf("    %d %d\n",i+1,num_angle[i]);
      for (int j = 0; j < num_angle[i]; j++)
        printf("      %d %d %d %d %d\n",
               j+1,angle_type[i][j],
               angle_atom1[i][j],angle_atom2[i][j],angle_atom3[i][j]);
    }
  }
  if (dihedralflag) {
    printf(  "Dihedrals:\n");
    for (int i = 0; i < natoms; i++) {
      printf("    %d %d\n",i+1,num_dihedral[i]);
      for (int j = 0; j < num_dihedral[i]; j++)
        printf("      %d %d %d %d %d %d\n",
               j+1,dihedral_type[i][j],
               dihedral_atom1[i][j],dihedral_atom2[i][j],
               dihedral_atom3[i][j],dihedral_atom4[i][j]);
    }
  }
  if (improperflag) {
    printf(  "Impropers:\n");
    for (int i = 0; i < natoms; i++) {
      printf("    %d %d\n",i+1,num_improper[i]);
      for (int j = 0; j < num_improper[i]; j++)
        printf("      %d %d %d %d %d %d\n",
               j+1,improper_type[i][j],
               improper_atom1[i][j],improper_atom2[i][j],
               improper_atom3[i][j],improper_atom4[i][j]);
    }
  }

  if (specialflag) {
    printf(  "Special neighs:\n");
    for (int i = 0; i < natoms; i++) {
      printf("    %d %d %d %d\n",i+1,
             nspecial[i][0],nspecial[i][1]-nspecial[i][0],
             nspecial[i][2]-nspecial[i][1]);
      printf("      ");
      for (int j = 0; j < nspecial[i][2]; j++)
        printf(" %d",special[i][j]);
      printf("\n");
    }
  }
}

*/
