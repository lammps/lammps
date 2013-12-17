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

#include "stdlib.h"
#include "string.h"
#include "molecule.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "comm.h"
#include "domain.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 256

/* ---------------------------------------------------------------------- */

Molecule::Molecule(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  me = comm->me;

  if (narg != 2) error->all(FLERR,"Illegal molecule command");

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR,
                 "Molecule ID must be alphanumeric or underscore characters");

  char *file = arg[1];

  // initialize all fields to empty

  initialize();

  // scan file for sizes of all fields and allocate them

  if (me == 0) open(file);
  read(0);
  if (me == 0) fclose(fp);
  allocate();

  // read file again to populate all fields

  if (me == 0) open(file);
  read(1);
  if (me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

Molecule::~Molecule()
{
  delete [] id;
  deallocate();
}

/* ----------------------------------------------------------------------
   compute center = geometric center of molecule
   also compute dx = displacement of each atom from origin
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
}

/* ----------------------------------------------------------------------
   compute xcm = center or mass of molecule
   account for finite size particles
   also compute dxcom = displacement of each atom from COM
------------------------------------------------------------------------- */

void Molecule::compute_xcm()
{
  comflag = 1;

  // NOTE this is not yet right

  xcm[0] = xcm[1] = xcm[2] = 0.0;
  for (int i = 0; i < natoms; i++) {
    xcm[0] += x[i][0];
    xcm[1] += x[i][1];
    xcm[2] += x[i][2];
  }
  xcm[0] /= natoms;
  xcm[1] /= natoms;
  xcm[2] /= natoms;

  memory->destroy(dxcom);
  memory->create(dxcom,natoms,3,"molecule:dxcom");

  for (int i = 0; i < natoms; i++) {
    dxcom[i][0] = x[i][0] - center[0];
    dxcom[i][1] = x[i][1] - center[1];
    dxcom[i][2] = x[i][2] - center[2];
  }
}

/* ----------------------------------------------------------------------
   read molecule info from file
   flag = 0, just scan for sizes of fields
   flag = 1, read and store fields
------------------------------------------------------------------------- */

void Molecule::read(int flag)
{
  char line[MAXLINE],keyword[MAXLINE];
  char *eof,*ptr;

  // skip 1st line of file

  if (me == 0) {
    eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one(FLERR,"Unexpected end of molecule file");
  }

  // read header lines
  // skip blank lines or lines that start with "#"
  // stop when read an unrecognized line

  while (1) {

    readline(line);

    // trim anything from '#' onward
    // if line is blank, continue

    if (ptr = strchr(line,'#')) *ptr = '\0';
    if (strspn(line," \t\n\r") == strlen(line)) continue;

    // search line for header keywords and set corresponding variable

    if (strstr(line,"atoms")) sscanf(line,"%d",&natoms);
    else if (strstr(line,"bonds")) sscanf(line,"%d",&nbonds);
    else if (strstr(line,"angles")) sscanf(line,"%d",&nangles);
    else if (strstr(line,"dihedrals")) sscanf(line,"%d",&ndihedrals);
    else if (strstr(line,"impropers")) sscanf(line,"%d",&nimpropers);
    else break;
  }

  // error check

  if (flag == 0) {
    if (natoms == 0) error->all(FLERR,"No atom count in molecule file");

    if (nbonds && !atom->avec->bonds_allow)
      error->all(FLERR,"Bonds in molecule file not supported by atom style");
    if (nangles && !atom->avec->angles_allow)
      error->all(FLERR,"Angles in molecule file not supported by atom style");
    if (ndihedrals && !atom->avec->dihedrals_allow)
      error->all(FLERR,
                 "Dihedrals in molecule file not supported by atom style");
    if (nimpropers && !atom->avec->impropers_allow)
      error->all(FLERR,
                 "Impropers in molecule file not supported by atom style");
  }

  // count = vector for tallying bonds,angles,etc per atom

  if (flag == 0) memory->create(count,natoms,"molecule:count");
  else count = NULL;
  
  // grab keyword and skip next line

  parse_keyword(0,line,keyword);
  readline(line);

  // loop over sections of molecule file

  while (strlen(keyword)) {
    if (strcmp(keyword,"Coords") == 0) {
      xflag = 1;
      if (flag) coords(line);
      else skip_lines(natoms,line);
    } else if (strcmp(keyword,"Types") == 0) {
      typeflag = 1;
      if (flag) types(line);
      else skip_lines(natoms,line);
    } else if (strcmp(keyword,"Charges") == 0) {
      qflag = 1;
      if (flag) charges(line);
      else skip_lines(natoms,line);
    } else if (strcmp(keyword,"Diameters") == 0) {
      radiusflag = 1;
      if (flag) diameters(line);
      else skip_lines(natoms,line);
    } else if (strcmp(keyword,"Masses") == 0) {
      rmassflag = 1;
      if (flag) masses(line);
      else skip_lines(natoms,line);

    } else if (strcmp(keyword,"Bonds") == 0) {
      if (nbonds == 0)
	error->all(FLERR,"Molecule file has bonds but no nbonds setting");
      bondflag = 1;
      bonds(flag,line);
    } else if (strcmp(keyword,"Angles") == 0) {
      if (nangles == 0) 
	error->all(FLERR,"Molecule file has angles but no nangles setting");
      angleflag = 1;
      angles(flag,line);
    } else if (strcmp(keyword,"Dihedrals") == 0) {
      if (ndihedrals == 0) error->all(FLERR,"Molecule file has dihedrals "
				      "but no ndihedrals setting");
      dihedralflag = 1;
      dihedrals(flag,line);
    } else if (strcmp(keyword,"Impropers") == 0) {
      if (nimpropers == 0) error->all(FLERR,"Molecule file has impropers "
				      "but no nimpropers setting");
      improperflag = 1;
      impropers(flag,line);

    } else if (strcmp(keyword,"Special Bond Counts") == 0) {
      nspecialflag = 1;
      nspecial_read(flag,line);
    } else if (strcmp(keyword,"Special Bonds") == 0) {
      specialflag = 1;
      if (flag) special_read(line);
      else skip_lines(natoms,line);

    } else error->one(FLERR,"Unknown section in molecule file");
	     
    parse_keyword(1,line,keyword);
  }

  // clean up

  memory->destroy(count);

  // error check

  if (flag == 0) {
    if (qflag && !atom->q_flag)
      error->all(FLERR,"Molecule file has undefined atom property");
    if (radiusflag && !atom->radius_flag)
      error->all(FLERR,"Molecule file has undefined atom property");
    if (rmassflag && !atom->rmass_flag)
      error->all(FLERR,"Molecule file has undefined atom property");

    if (bondflag && !atom->avec->bonds_allow)
      error->all(FLERR,"Invalid molecule file section: Bonds");
    if (angleflag && !atom->avec->angles_allow)
      error->all(FLERR,"Invalid molecule file section: Angles");
    if (dihedralflag && !atom->avec->dihedrals_allow)
      error->all(FLERR,"Invalid molecule file section: Dihedrals");
    if (improperflag && !atom->avec->impropers_allow)
      error->all(FLERR,"Invalid molecule file section: Impropers");

    if (bond_per_atom > atom->bond_per_atom ||
	angle_per_atom > atom->angle_per_atom ||
	dihedral_per_atom > atom->dihedral_per_atom ||
	improper_per_atom > atom->improper_per_atom)
      error->all(FLERR,"Molecule file bond/angle/etc counts "
		 "per atom are too large");

    if ((nspecialflag && !specialflag) || (!nspecialflag && specialflag))
      error->all(FLERR,"Molecule file needs both Special Bond sections");
    if (specialflag && !bondflag) 
      error->all(FLERR,"Molecule file has special flags but no bonds");
    if (maxspecial > atom->maxspecial)
      error->all(FLERR,"Molecule file special bond counts are too large");
  }
}

/* ----------------------------------------------------------------------
   read coords from file
------------------------------------------------------------------------- */

void Molecule::coords(char *line)
{
  int tmp;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    sscanf(line,"%d %lg %lg %lg",&tmp,&x[i][0],&x[i][1],&x[i][2]);
  }

  if (domain->dimension == 2) {
    for (int i = 0; i < natoms; i++)
      if (x[i][2] != 0.0) 
        error->all(FLERR,"Molecule file z coord must be 0.0 for 2d");
  }
}

/* ----------------------------------------------------------------------
   read types from file
------------------------------------------------------------------------- */

void Molecule::types(char *line)
{
  int tmp;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    sscanf(line,"%d %d",&tmp,&type[i]);
  }
}

/* ----------------------------------------------------------------------
   read charges from file
------------------------------------------------------------------------- */

void Molecule::charges(char *line)
{
  int tmp;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    sscanf(line,"%d %lg",&tmp,&q[i]);
  }
}

/* ----------------------------------------------------------------------
   read diameters from file and set radii
------------------------------------------------------------------------- */

void Molecule::diameters(char *line)
{
  int tmp;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    sscanf(line,"%d %lg",&tmp,&radius[i]);
    radius[i] *= 0.5;
  }
}

/* ----------------------------------------------------------------------
   read masses from file
------------------------------------------------------------------------- */

void Molecule::masses(char *line)
{
  int tmp;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    sscanf(line,"%d %lg",&tmp,&rmass[i]);
  }
}

/* ----------------------------------------------------------------------
   read bonds from file
   store each with both atoms if newton_bond = 0
   if flag = 0, just count bonds/atom
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void Molecule::bonds(int flag, char *line)
{
  int m,tmp,itype,atom1,atom2;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_bond[i] = 0;

  for (int i = 0; i < nbonds; i++) {
    readline(line);
    sscanf(line,"%d %d %d %d",&tmp,&itype,&atom1,&atom2);

    if (atom1 <= 0 || atom1 > natoms ||
	atom2 <= 0 || atom2 > natoms)
      error->one(FLERR,"Invalid atom ID in Bonds section of molecule file");
    if (itype <= 0 || itype > atom->nbondtypes)
      error->one(FLERR,"Invalid bond type in Bonds section of molecule file");

    if (flag) {
      m = atom1-1;
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
  int m,tmp,itype,atom1,atom2,atom3;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_angle[i] = 0;

  for (int i = 0; i < nangles; i++) {
    readline(line);
    sscanf(line,"%d %d %d %d %d",&tmp,&itype,&atom1,&atom2,&atom3);

    if (atom1 <= 0 || atom1 > natoms ||
        atom2 <= 0 || atom2 > natoms ||
        atom3 <= 0 || atom3 > natoms)
      error->one(FLERR,"Invalid atom ID in Angles section of molecule file");
    if (itype <= 0 || itype > atom->nangletypes)
      error->one(FLERR,"Invalid angle type in Angles section of molecule file");

    if (flag) {
      m = atom2-1;
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
  int m,tmp,itype,atom1,atom2,atom3,atom4;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_dihedral[i] = 0;

  for (int i = 0; i < ndihedrals; i++) {
    readline(line);
    sscanf(line,"%d %d %d %d %d %d",&tmp,&itype,&atom1,&atom2,&atom3,&atom4);

    if (atom1 <= 0 || atom1 > natoms ||
        atom2 <= 0 || atom2 > natoms ||
        atom3 <= 0 || atom3 > natoms ||
        atom4 <= 0 || atom4 > natoms)
      error->one(FLERR,
		 "Invalid atom ID in dihedrals section of molecule file");
    if (itype <= 0 || itype > atom->ndihedraltypes)
      error->one(FLERR,
		 "Invalid dihedral type in dihedrals section of molecule file");

    if (flag) {
      m = atom2-1;
      dihedral_type[m][num_dihedral[m]] = itype;
      dihedral_atom1[m][num_dihedral[m]] = atom1;
      dihedral_atom2[m][num_dihedral[m]] = atom2;
      dihedral_atom3[m][num_dihedral[m]] = atom3;
      num_dihedral[m]++;
      if (newton_bond == 0) {
	m = atom1-1;
	dihedral_type[m][num_dihedral[m]] = itype;
	dihedral_atom1[m][num_dihedral[m]] = atom1;
	dihedral_atom2[m][num_dihedral[m]] = atom2;
	dihedral_atom3[m][num_dihedral[m]] = atom3;
	num_dihedral[m]++;
	m = atom3-1;
	dihedral_type[m][num_dihedral[m]] = itype;
	dihedral_atom1[m][num_dihedral[m]] = atom1;
	dihedral_atom2[m][num_dihedral[m]] = atom2;
	dihedral_atom3[m][num_dihedral[m]] = atom3;
	num_dihedral[m]++;
	m = atom4-1;
	dihedral_type[m][num_dihedral[m]] = itype;
	dihedral_atom1[m][num_dihedral[m]] = atom1;
	dihedral_atom2[m][num_dihedral[m]] = atom2;
	dihedral_atom3[m][num_dihedral[m]] = atom3;
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
  int m,tmp,itype,atom1,atom2,atom3,atom4;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_improper[i] = 0;

  for (int i = 0; i < nimpropers; i++) {
    readline(line);
    sscanf(line,"%d %d %d %d %d %d",&tmp,&itype,&atom1,&atom2,&atom3,&atom4);

    if (atom1 <= 0 || atom1 > natoms ||
        atom2 <= 0 || atom2 > natoms ||
        atom3 <= 0 || atom3 > natoms ||
        atom4 <= 0 || atom4 > natoms)
      error->one(FLERR,
		 "Invalid atom ID in impropers section of molecule file");
    if (itype <= 0 || itype > atom->nimpropertypes)
      error->one(FLERR,
		 "Invalid improper type in impropers section of molecule file");

    if (flag) {
      m = atom2-1;
      improper_type[m][num_improper[m]] = itype;
      improper_atom1[m][num_improper[m]] = atom1;
      improper_atom2[m][num_improper[m]] = atom2;
      improper_atom3[m][num_improper[m]] = atom3;
      num_improper[m]++;
      if (newton_bond == 0) {
	m = atom1-1;
	improper_type[m][num_improper[m]] = itype;
	improper_atom1[m][num_improper[m]] = atom1;
	improper_atom2[m][num_improper[m]] = atom2;
	improper_atom3[m][num_improper[m]] = atom3;
	num_improper[m]++;
	m = atom3-1;
	improper_type[m][num_improper[m]] = itype;
	improper_atom1[m][num_improper[m]] = atom1;
	improper_atom2[m][num_improper[m]] = atom2;
	improper_atom3[m][num_improper[m]] = atom3;
	num_improper[m]++;
	m = atom4-1;
	improper_type[m][num_improper[m]] = itype;
	improper_atom1[m][num_improper[m]] = atom1;
	improper_atom2[m][num_improper[m]] = atom2;
	improper_atom3[m][num_improper[m]] = atom3;
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
  int tmp,c1,c2,c3;

  if (flag == 0) maxspecial = 0;

  for (int i = 0; i < natoms; i++) {
    readline(line);
    sscanf(line,"%d %d %d %d",&tmp,&c1,&c2,&c3);

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
  int m,nwords;
  char **words = new char*[maxspecial+1];

  for (int i = 0; i < natoms; i++) {
    readline(line);
    nwords = parse(line,words,maxspecial+1);
    if (nwords != nspecial[i][2]+1)
      error->all(FLERR,"Molecule file special list "
		 "does not match special count");

    for (m = 1; m < nwords; m++) {
      special[i][m-1] = atoi(words[m]);
      if (special[i][m-1] <= 0 || special[i][m-1] > natoms ||
	  special[i][m-1] == i+1)
	error->all(FLERR,"Invalid special atom index in molecule file");
    }
  }

  delete [] words;
}

/* ----------------------------------------------------------------------
   init all data structures to empty
------------------------------------------------------------------------- */

void Molecule::initialize()
{
  natoms = 0;
  nbonds = nangles = ndihedrals = nimpropers = 0;
  bond_per_atom = angle_per_atom = dihedral_per_atom = improper_per_atom = 0;
  maxspecial = 0;

  xflag = typeflag = qflag = radiusflag = rmassflag = 0;
  bondflag = angleflag = dihedralflag = improperflag = 0;
  nspecialflag = specialflag = 0;

  centerflag = comflag = inertiaflag = 0;

  x = NULL;
  type = NULL;
  q = NULL;
  radius = NULL;
  rmass = NULL;

  num_bond = NULL;
  bond_type = NULL;
  bond_atom = NULL;

  num_angle = NULL;
  angle_type = NULL;
  angle_atom1 = angle_atom2 = angle_atom3 = NULL;

  num_dihedral = NULL;
  dihedral_type = NULL;
  dihedral_atom1 = dihedral_atom2 = dihedral_atom3 = dihedral_atom4 = NULL;

  num_improper = NULL;
  improper_type = NULL;
  improper_atom1 = improper_atom2 = improper_atom3 = improper_atom4 = NULL;

  nspecial = NULL;
  special = NULL;

  dx = NULL;
  dxcom = NULL;
}

/* ----------------------------------------------------------------------
   allocate all data structures
------------------------------------------------------------------------- */

void Molecule::allocate()
{
  if (xflag) memory->create(x,natoms,3,"molecule:x");
  if (typeflag) memory->create(type,natoms,"molecule:type");
  if (qflag) memory->create(q,natoms,"molecule:q");
  if (radiusflag) memory->create(radius,natoms,"molecule:radius");
  if (rmassflag) memory->create(rmass,natoms,"molecule:rmass");
  
  if (bondflag) {
    memory->create(num_bond,natoms,"molecule:num_bond");
    memory->create(bond_type,natoms,bond_per_atom,
		   "molecule:bond_type");
    memory->create(bond_atom,natoms,bond_per_atom,
		   "molecule:bond_atom");
  }

  if (angleflag) {
    memory->create(num_angle,natoms,"molecule:num_angle");
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
    memory->create(num_dihedral,natoms,"molecule:num_dihedral");
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
    memory->create(num_improper,natoms,"molecule:num_improper");
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

  if (nspecialflag)
    memory->create(nspecial,natoms,3,"molecule:nspecial");
  if (specialflag)
    memory->create(special,natoms,maxspecial,"molecule:special");
}

/* ----------------------------------------------------------------------
   deallocate all data structures
------------------------------------------------------------------------- */

void Molecule::deallocate()
{
  memory->destroy(x);
  memory->destroy(type);
  memory->destroy(q);
  memory->destroy(radius);
  memory->destroy(rmass);
  
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

  memory->destroy(dx);
  memory->destroy(dxcom);
}

/* ----------------------------------------------------------------------
   open molecule file
------------------------------------------------------------------------- */

void Molecule::open(char *file)
{
  fp = fopen(file,"r");
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open molecule file %s",file);
    error->one(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   read and bcast a line
------------------------------------------------------------------------- */

void Molecule::readline(char *line)
{
  int n;
  if (me == 0) {
    if (fgets(line,MAXLINE,fp) == NULL) n = 0;
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

void Molecule::parse_keyword(int flag, char *line, char *keyword)
{
  if (flag) {

    // read upto non-blank line plus 1 following line
    // eof is set to 1 if any read hits end-of-file

    int eof = 0;
    if (me == 0) {
      if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
      while (eof == 0 && strspn(line," \t\n\r") == strlen(line)) {
	if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
      }
      if (fgets(keyword,MAXLINE,fp) == NULL) eof = 1;
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
  }

  // copy non-whitespace portion of line into keyword

  int start = strspn(line," \t\n\r");
  int stop = strlen(line) - 1;
  while (line[stop] == ' ' || line[stop] == '\t'
         || line[stop] == '\n' || line[stop] == '\r') stop--;
  line[stop+1] = '\0';
  strcpy(keyword,&line[start]);
}

/* ----------------------------------------------------------------------
   skip N lines of file
------------------------------------------------------------------------- */

void Molecule::skip_lines(int n, char *line)
{
  for (int i = 0; i < n; i++) readline(line);
}

/* ----------------------------------------------------------------------
   parse line into words separated by whitespace
   return # of words
   max = max pointers storable in words
------------------------------------------------------------------------- */

int Molecule::parse(char *line, char **words, int max)
{
  char *ptr;

  int nwords = 0;
  words[nwords++] = strtok(line," \t\n\r\f");

  while (ptr = strtok(NULL," \t\n\r\f")) {
    if (nwords < max) words[nwords] = ptr;
    nwords++;
  }

  return nwords;
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
