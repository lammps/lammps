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

#include "lmptype.h"
#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "ctype.h"
#include "read_data.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
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
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

#define MAXLINE 256
#define LB_FACTOR 1.1
#define CHUNK 1024
#define DELTA 4            // must be 2 or larger
#define MAXBODY 20         // max # of lines in one body, also in Atom class

                           // customize for new sections
#define NSECTIONS 25       // change when add to header::section_keywords

/* ---------------------------------------------------------------------- */

ReadData::ReadData(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  line = new char[MAXLINE];
  keyword = new char[MAXLINE];
  buffer = new char[CHUNK*MAXLINE];
  narg = maxarg = 0;
  arg = NULL;

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
  delete [] keyword;
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

  if (domain->box_exist)
    error->all(FLERR,"Cannot read_data after simulation box is defined");
  if (domain->dimension == 2 && domain->zperiodic == 0)
    error->all(FLERR,"Cannot run 2d simulation with nonperiodic Z dimension");

  // fixes that process data file info

  nfix = 0;
  fix_index = NULL;
  fix_header = NULL;
  fix_section = NULL;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"fix") == 0) {
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
      int n = strlen(arg[iarg+2]) + 1;
      fix_header[nfix] = new char[n];
      strcpy(fix_header[nfix],arg[iarg+2]);
      n = strlen(arg[iarg+3]) + 1;
      fix_section[nfix] = new char[n];
      strcpy(fix_section[nfix],arg[iarg+3]);
      nfix++;
      iarg += 4;
    } else error->all(FLERR,"Illegal read_data command");
  }

  // scan data file to determine max topology needed per atom
  // allocate initial topology arrays

  if (atom->molecular) {
    if (me == 0) {
      if (screen) fprintf(screen,"Scanning data file ...\n");
      open(arg[0]);
      header(0);
      scan(atom->bond_per_atom,atom->angle_per_atom,
           atom->dihedral_per_atom,atom->improper_per_atom);
      if (compressed) pclose(fp);
      else fclose(fp);
      atom->bond_per_atom += atom->extra_bond_per_atom;
    }

    MPI_Bcast(&atom->bond_per_atom,1,MPI_INT,0,world);
    MPI_Bcast(&atom->angle_per_atom,1,MPI_INT,0,world);
    MPI_Bcast(&atom->dihedral_per_atom,1,MPI_INT,0,world);
    MPI_Bcast(&atom->improper_per_atom,1,MPI_INT,0,world);

  } else
    atom->bond_per_atom = atom->angle_per_atom =
      atom->dihedral_per_atom = atom->improper_per_atom = 0;

  // read header info

  if (me == 0) {
    if (screen) fprintf(screen,"Reading data file ...\n");
    open(arg[0]);
  } else fp = NULL;
  header(1);
  domain->box_exist = 1;

  // problem setup using info from header

  update->ntimestep = 0;

  int n;
  if (comm->nprocs == 1) n = static_cast<int> (atom->natoms);
  else n = static_cast<int> (LB_FACTOR * atom->natoms / comm->nprocs);

  atom->allocate_type_arrays();
  atom->avec->grow(n);
  n = atom->nmax;

  domain->print_box("  ");
  domain->set_initial_box();
  domain->set_global_box();
  comm->set_proc_grid();
  domain->set_local_box();

  // customize for new sections
  // read rest of file in free format

  int atomflag = 0;

  while (strlen(keyword)) {

    // allow special fixes first chance to match and process the section
    // if fix matches, continue to next section

    if (nfix) {
      for (n = 0; n < nfix; n++)
        if (strstr(line,fix_section[n])) {
          bigint nlines =
            modify->fix[fix_index[n]]->read_data_skip_lines(keyword);
          fix(n,keyword,nlines);
          parse_keyword(0,1);
          break;
        }
      if (n < nfix) continue;
    }

    if (strcmp(keyword,"Atoms") == 0) {
      atoms();
      atomflag = 1;
    } else if (strcmp(keyword,"Velocities") == 0) {
      if (atomflag == 0) error->all(FLERR,"Must read Atoms before Velocities");
      velocities();

    } else if (strcmp(keyword,"Ellipsoids") == 0) {
      if (!avec_ellipsoid)
        error->all(FLERR,"Invalid data file section: Ellipsoids");
      if (atomflag == 0) error->all(FLERR,"Must read Atoms before Ellipsoids");
      bonus(nellipsoids,(AtomVec *) avec_ellipsoid,"ellipsoids");
    } else if (strcmp(keyword,"Lines") == 0) {
      if (!avec_line)
        error->all(FLERR,"Invalid data file section: Lines");
      if (atomflag == 0) error->all(FLERR,"Must read Atoms before Lines");
      bonus(nlines,(AtomVec *) avec_line,"lines");
    } else if (strcmp(keyword,"Triangles") == 0) {
      if (!avec_tri)
        error->all(FLERR,"Invalid data file section: Triangles");
      if (atomflag == 0) error->all(FLERR,"Must read Atoms before Triangles");
      bonus(ntris,(AtomVec *) avec_tri,"triangles");
    } else if (strcmp(keyword,"Bodies") == 0) {
      if (!avec_body)
        error->all(FLERR,"Invalid data file section: Bodies");
      if (atomflag == 0) error->all(FLERR,"Must read Atoms before Bodies");
      bodies();

    } else if (strcmp(keyword,"Bonds") == 0) {
      if (atom->avec->bonds_allow == 0)
        error->all(FLERR,"Invalid data file section: Bonds");
      if (atomflag == 0) error->all(FLERR,"Must read Atoms before Bonds");
      bonds();
    } else if (strcmp(keyword,"Angles") == 0) {
      if (atom->avec->angles_allow == 0)
        error->all(FLERR,"Invalid data file section: Angles");
      if (atomflag == 0) error->all(FLERR,"Must read Atoms before Angles");
      angles();
    } else if (strcmp(keyword,"Dihedrals") == 0) {
      if (atom->avec->dihedrals_allow == 0)
        error->all(FLERR,"Invalid data file section: Dihedrals");
      if (atomflag == 0) error->all(FLERR,"Must read Atoms before Dihedrals");
      dihedrals();
    } else if (strcmp(keyword,"Impropers") == 0) {
      if (atom->avec->impropers_allow == 0)
        error->all(FLERR,"Invalid data file section: Impropers");
      if (atomflag == 0) error->all(FLERR,"Must read Atoms before Impropers");
      impropers();

    } else if (strcmp(keyword,"Masses") == 0) {
      mass();
    } else if (strcmp(keyword,"Pair Coeffs") == 0) {
      if (force->pair == NULL)
        error->all(FLERR,"Must define pair_style before Pair Coeffs");
      paircoeffs();
    } else if (strcmp(keyword,"PairIJ Coeffs") == 0) {
      if (force->pair == NULL)
        error->all(FLERR,"Must define pair_style before PairIJ Coeffs");
      pairIJcoeffs();
    } else if (strcmp(keyword,"Bond Coeffs") == 0) {
      if (atom->avec->bonds_allow == 0)
        error->all(FLERR,"Invalid data file section: Bond Coeffs");
      if (force->bond == NULL)
        error->all(FLERR,"Must define bond_style before Bond Coeffs");
      bondcoeffs();
    } else if (strcmp(keyword,"Angle Coeffs") == 0) {
      if (atom->avec->angles_allow == 0)
        error->all(FLERR,"Invalid data file section: Angle Coeffs");
      if (force->angle == NULL)
        error->all(FLERR,"Must define angle_style before Angle Coeffs");
      anglecoeffs(0);
    } else if (strcmp(keyword,"Dihedral Coeffs") == 0) {
      if (atom->avec->dihedrals_allow == 0)
        error->all(FLERR,"Invalid data file section: Dihedral Coeffs");
      if (force->dihedral == NULL)
        error->all(FLERR,"Must define dihedral_style before Dihedral Coeffs");
      dihedralcoeffs(0);
    } else if (strcmp(keyword,"Improper Coeffs") == 0) {
      if (atom->avec->impropers_allow == 0)
        error->all(FLERR,"Invalid data file section: Improper Coeffs");
      if (force->improper == NULL)
        error->all(FLERR,"Must define improper_style before Improper Coeffs");
      impropercoeffs(0);

    } else if (strcmp(keyword,"BondBond Coeffs") == 0) {
      if (atom->avec->angles_allow == 0)
        error->all(FLERR,"Invalid data file section: BondBond Coeffs");
      if (force->angle == NULL)
        error->all(FLERR,"Must define angle_style before BondBond Coeffs");
      anglecoeffs(1);
    } else if (strcmp(keyword,"BondAngle Coeffs") == 0) {
      if (atom->avec->angles_allow == 0)
        error->all(FLERR,"Invalid data file section: BondAngle Coeffs");
      if (force->angle == NULL)
        error->all(FLERR,"Must define angle_style before BondAngle Coeffs");
      anglecoeffs(2);

    } else if (strcmp(keyword,"MiddleBondTorsion Coeffs") == 0) {
      if (atom->avec->dihedrals_allow == 0)
        error->all(FLERR,"Invalid data file section: MiddleBondTorsion Coeffs");
      if (force->dihedral == NULL)
        error->all(FLERR,
                   "Must define dihedral_style before "
                   "MiddleBondTorsion Coeffs");
      dihedralcoeffs(1);
    } else if (strcmp(keyword,"EndBondTorsion Coeffs") == 0) {
      if (atom->avec->dihedrals_allow == 0)
        error->all(FLERR,"Invalid data file section: EndBondTorsion Coeffs");
      if (force->dihedral == NULL)
        error->all(FLERR,
                   "Must define dihedral_style before EndBondTorsion Coeffs");
      dihedralcoeffs(2);
    } else if (strcmp(keyword,"AngleTorsion Coeffs") == 0) {
      if (atom->avec->dihedrals_allow == 0)
        error->all(FLERR,"Invalid data file section: AngleTorsion Coeffs");
      if (force->dihedral == NULL)
        error->all(FLERR,
                   "Must define dihedral_style before AngleTorsion Coeffs");
      dihedralcoeffs(3);
    } else if (strcmp(keyword,"AngleAngleTorsion Coeffs") == 0) {
      if (atom->avec->dihedrals_allow == 0)
        error->all(FLERR,"Invalid data file section: AngleAngleTorsion Coeffs");
      if (force->dihedral == NULL)
        error->all(FLERR,
                   "Must define dihedral_style before "
                   "AngleAngleTorsion Coeffs");
      dihedralcoeffs(4);
    } else if (strcmp(keyword,"BondBond13 Coeffs") == 0) {
      if (atom->avec->dihedrals_allow == 0)
        error->all(FLERR,"Invalid data file section: BondBond13 Coeffs");
      if (force->dihedral == NULL)
        error->all(FLERR,"Must define dihedral_style before BondBond13 Coeffs");
      dihedralcoeffs(5);

    } else if (strcmp(keyword,"AngleAngle Coeffs") == 0) {
      if (atom->avec->impropers_allow == 0)
        error->all(FLERR,"Invalid data file section: AngleAngle Coeffs");
      if (force->improper == NULL)
        error->all(FLERR,"Must define improper_style before AngleAngle Coeffs");
      impropercoeffs(1);

    } else {
      char str[128];
      sprintf(str,"Unknown identifier in data file: %s",keyword);
      error->all(FLERR,str);
    }

    parse_keyword(0,1);
  }

  // close file

  if (me == 0) {
    if (compressed) pclose(fp);
    else fclose(fp);
  }

  // error if natoms > 0 yet no atoms were read

  if (atom->natoms > 0 && atomflag == 0)
    error->all(FLERR,"No atoms in data file");

  // create bond topology now that system is defined

  if (atom->molecular) {
    Special special(lmp);
    special.build();
  }
}

/* ----------------------------------------------------------------------
   read free-format header of data file
   if flag = 0, only called by proc 0
   if flag = 1, called by all procs so bcast lines as read them
   1st line and blank lines are skipped
   non-blank lines are checked for header keywords and leading value is read
   header ends with EOF or non-blank line containing no header keyword
     if EOF, line is set to blank line
     else line has first keyword line for rest of file
------------------------------------------------------------------------- */

void ReadData::header(int flag)
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

  // customize for new header lines

  while (1) {

    // read a line and bcast length if flag is set

    if (me == 0) {
      if (fgets(line,MAXLINE,fp) == NULL) n = 0;
      else n = strlen(line) + 1;
    }
    if (flag) MPI_Bcast(&n,1,MPI_INT,0,world);

    // if n = 0 then end-of-file so return with blank line

    if (n == 0) {
      line[0] = '\0';
      return;
    }

    // bcast line if flag is set

    if (flag) MPI_Bcast(line,n,MPI_CHAR,0,world);

    // trim anything from '#' onward
    // if line is blank, continue

    if (ptr = strchr(line,'#')) *ptr = '\0';
    if (strspn(line," \t\n\r") == strlen(line)) continue;

    // allow special fixes first chance to match and process the line
    // if fix matches, continue to next header line

    if (nfix) {
      for (n = 0; n < nfix; n++)
        if (strstr(line,fix_header[n])) {
          modify->fix[fix_index[n]]->read_data_header(line);
          break;
        }
      if (n < nfix) continue;
    }

    // search line for header keyword and set corresponding variable

    if (strstr(line,"atoms")) {
      sscanf(line,BIGINT_FORMAT,&atom->natoms);

    // check for these first
    // otherwise "triangles" will be matched as "angles"

    } else if (strstr(line,"ellipsoids")) {
      if (!avec_ellipsoid && me == 0)
        error->one(FLERR,"No ellipsoids allowed with this atom style");
      sscanf(line,BIGINT_FORMAT,&nellipsoids);
    } else if (strstr(line,"lines")) {
      if (!avec_line && me == 0)
        error->one(FLERR,"No lines allowed with this atom style");
      sscanf(line,BIGINT_FORMAT,&nlines);
    } else if (strstr(line,"triangles")) {
      if (!avec_tri && me == 0)
        error->one(FLERR,"No triangles allowed with this atom style");
      sscanf(line,BIGINT_FORMAT,&ntris);
    } else if (strstr(line,"bodies")) {
      if (!avec_body && me == 0)
        error->one(FLERR,"No bodies allowed with this atom style");
      sscanf(line,BIGINT_FORMAT,&nbodies);
    }

    else if (strstr(line,"bonds")) sscanf(line,BIGINT_FORMAT,&atom->nbonds);
    else if (strstr(line,"angles")) sscanf(line,BIGINT_FORMAT,&atom->nangles);
    else if (strstr(line,"dihedrals")) sscanf(line,BIGINT_FORMAT,
                                              &atom->ndihedrals);
    else if (strstr(line,"impropers")) sscanf(line,BIGINT_FORMAT,
                                              &atom->nimpropers);

    else if (strstr(line,"atom types")) sscanf(line,"%d",&atom->ntypes);
    else if (strstr(line,"bond types")) sscanf(line,"%d",&atom->nbondtypes);
    else if (strstr(line,"angle types")) sscanf(line,"%d",&atom->nangletypes);
    else if (strstr(line,"dihedral types"))
      sscanf(line,"%d",&atom->ndihedraltypes);
    else if (strstr(line,"improper types"))
      sscanf(line,"%d",&atom->nimpropertypes);

    else if (strstr(line,"extra bond per atom"))
      sscanf(line,"%d",&atom->extra_bond_per_atom);

    else if (strstr(line,"xlo xhi"))
      sscanf(line,"%lg %lg",&domain->boxlo[0],&domain->boxhi[0]);
    else if (strstr(line,"ylo yhi"))
      sscanf(line,"%lg %lg",&domain->boxlo[1],&domain->boxhi[1]);
    else if (strstr(line,"zlo zhi"))
      sscanf(line,"%lg %lg",&domain->boxlo[2],&domain->boxhi[2]);
    else if (strstr(line,"xy xz yz")) {
      domain->triclinic = 1;
      sscanf(line,"%lg %lg %lg",&domain->xy,&domain->xz,&domain->yz);
    } else break;
  }

  // error check on total system size

  if (atom->natoms < 0 || atom->natoms > MAXBIGINT ||
      atom->nbonds < 0 || atom->nbonds > MAXBIGINT ||
      atom->nangles < 0 || atom->nangles > MAXBIGINT ||
      atom->ndihedrals < 0 || atom->ndihedrals > MAXBIGINT ||
      atom->nimpropers < 0 || atom->nimpropers > MAXBIGINT) {
    if (me == 0) error->one(FLERR,"System in data file is too big");
  }

  // check that exiting string is a valid section keyword

  parse_keyword(1,flag);
  for (n = 0; n < NSECTIONS; n++)
    if (strcmp(keyword,section_keywords[n]) == 0) break;
  if (n == NSECTIONS && me == 0) {
    char str[128];
    sprintf(str,"Unknown identifier in data file: %s",keyword);
    error->one(FLERR,str);
  }

  // error check on consistency of header values

  if ((atom->nbonds || atom->nbondtypes) &&
      atom->avec->bonds_allow == 0 && me == 0)
    error->one(FLERR,"No bonds allowed with this atom style");
  if ((atom->nangles || atom->nangletypes) &&
      atom->avec->angles_allow == 0 && me == 0)
    error->one(FLERR,"No angles allowed with this atom style");
  if ((atom->ndihedrals || atom->ndihedraltypes) &&
      atom->avec->dihedrals_allow == 0 && me == 0)
    error->one(FLERR,"No dihedrals allowed with this atom style");
  if ((atom->nimpropers || atom->nimpropertypes) &&
      atom->avec->impropers_allow == 0 && me == 0)
    error->one(FLERR,"No impropers allowed with this atom style");

  if (atom->nbonds > 0 && atom->nbondtypes <= 0 && me == 0)
    error->one(FLERR,"Bonds defined but no bond types");
  if (atom->nangles > 0 && atom->nangletypes <= 0 && me == 0)
    error->one(FLERR,"Angles defined but no angle types");
  if (atom->ndihedrals > 0 && atom->ndihedraltypes <= 0 && me == 0)
    error->one(FLERR,"Dihedrals defined but no dihedral types");
  if (atom->nimpropers > 0 && atom->nimpropertypes <= 0 && me == 0)
    error->one(FLERR,"Impropers defined but no improper types");
}

/* ----------------------------------------------------------------------
   read all atoms
------------------------------------------------------------------------- */

void ReadData::atoms()
{
  int i,m,nchunk,eof;

  bigint nread = 0;
  bigint natoms = atom->natoms;

  while (nread < natoms) {
    nchunk = MIN(natoms-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_atoms(nchunk,buffer);
    nread += nchunk;
  }

  // check that all atoms were assigned correctly

  bigint tmp = atom->nlocal;
  MPI_Allreduce(&tmp,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " atoms\n",natoms);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " atoms\n",natoms);
  }

  if (natoms != atom->natoms)
    error->all(FLERR,"Did not assign all atoms correctly");

  // if any atom ID < 0, error
  // if all atom IDs = 0, tag_enable = 0
  // if any atom ID > 0, error if any atom ID == 0
  // not checking if atom IDs > natoms or are unique

  int nlocal = atom->nlocal;
  int *tag = atom->tag;

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (tag[i] < 0) flag = 1;
  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all)
    error->all(FLERR,"Invalid atom ID in Atoms section of data file");

  flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (tag[i] > 0) flag = 1;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_MAX,world);
  if (flag_all == 0) atom->tag_enable = 0;

  if (atom->tag_enable) {
    flag = 0;
    for (int i = 0; i < nlocal; i++)
      if (tag[i] == 0) flag = 1;
    MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
    if (flag_all)
      error->all(FLERR,"Invalid atom ID in Atoms section of data file");
  }

  // create global mapping

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
  int i,m,nchunk,eof;

  int mapflag = 0;
  if (atom->map_style == 0) {
    mapflag = 1;
    atom->map_style = 1;
    atom->map_init();
    atom->map_set();
  }

  bigint nread = 0;
  bigint natoms = atom->natoms;

  while (nread < natoms) {
    nchunk = MIN(natoms-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_vels(nchunk,buffer);
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
   read all bonus data
   to find atoms, must build atom map if not a molecular system
------------------------------------------------------------------------- */

void ReadData::bonus(bigint nbonus, AtomVec *ptr, const char *type)
{
  int i,m,nchunk,eof;

  int mapflag = 0;
  if (atom->map_style == 0) {
    mapflag = 1;
    atom->map_style = 1;
    atom->map_init();
    atom->map_set();
  }

  bigint nread = 0;
  bigint natoms = nbonus;

  while (nread < natoms) {
    nchunk = MIN(natoms-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_bonus(nchunk,buffer,ptr);
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
------------------------------------------------------------------------- */

void ReadData::bodies()
{
  int i,m,nchunk,nmax,ninteger,ndouble,tmp,onebody;
  char *eof;

  int mapflag = 0;
  if (atom->map_style == 0) {
    mapflag = 1;
    atom->map_style = 1;
    atom->map_init();
    atom->map_set();
  }

  // nmax = max # of bodies to read in this chunk
  // nchunk = actual # read

  bigint nread = 0;
  bigint natoms = nbodies;

  while (nread < natoms) {
    nchunk = MIN(natoms-nread,CHUNK);

    if (me == 0) {
      nchunk = 0;
      nlines = 0;
      m = 0;

      while (nchunk < nmax && nlines <= CHUNK-MAXBODY) {
        eof = fgets(&buffer[m],MAXLINE,fp);
        if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
        sscanf(&buffer[m],"%d %d %d",&tmp,&ninteger,&ndouble);
        m += strlen(&buffer[m]);

        onebody = 0;
        if (ninteger) onebody += (ninteger-1)/10 + 1;
        if (ndouble) onebody += (ndouble-1)/10 + 1;
        if (onebody+1 > MAXBODY)
          error->one(FLERR,
                     "Too many lines in one body in data file - boost MAXBODY");

        for (i = 0; i < onebody; i++) {
          eof = fgets(&buffer[m],MAXLINE,fp);
          if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
          m += strlen(&buffer[m]);
        }

        nchunk++;
        nlines += onebody+1;
      }

      if (buffer[m-1] != '\n') strcpy(&buffer[m++],"\n");
      m++;
    }

    MPI_Bcast(&nchunk,1,MPI_INT,0,world);
    MPI_Bcast(&m,1,MPI_INT,0,world);
    MPI_Bcast(buffer,m,MPI_CHAR,0,world);

    atom->data_bodies(nchunk,buffer,avec_body);
    nread += nchunk;
  }

  if (mapflag) {
    atom->map_delete();
    atom->map_style = 0;
  }

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " bodies\n",natoms);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " bodies\n",natoms);
  }
}

/* ---------------------------------------------------------------------- */

void ReadData::bonds()
{
  int i,m,nchunk,eof;

  bigint nread = 0;
  bigint nbonds = atom->nbonds;

  bigint natoms = atom->natoms;

  while (nread < nbonds) {
    nchunk = MIN(nbonds-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_bonds(nchunk,buffer);
    nread += nchunk;
  }

  // check that bonds were assigned correctly

  int nlocal = atom->nlocal;
  bigint sum;
  bigint n = 0;
  for (i = 0; i < nlocal; i++) n += atom->num_bond[i];
  MPI_Allreduce(&n,&sum,1,MPI_LMP_BIGINT,MPI_SUM,world);
  int factor = 1;
  if (!force->newton_bond) factor = 2;

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " bonds\n",sum/factor);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " bonds\n",sum/factor);
  }
  if (sum != factor*atom->nbonds)
    error->all(FLERR,"Bonds assigned incorrectly");
}

/* ---------------------------------------------------------------------- */

void ReadData::angles()
{
  int i,m,nchunk,eof;

  bigint nread = 0;
  bigint nangles = atom->nangles;

  while (nread < nangles) {
    nchunk = MIN(nangles-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_angles(nchunk,buffer);
    nread += nchunk;
  }

  // check that ang

  int nlocal = atom->nlocal;
  bigint sum;
  bigint n = 0;
  for (i = 0; i < nlocal; i++) n += atom->num_angle[i];
  MPI_Allreduce(&n,&sum,1,MPI_LMP_BIGINT,MPI_SUM,world);
  int factor = 1;
  if (!force->newton_bond) factor = 3;

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " angles\n",sum/factor);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " angles\n",sum/factor);
  }
  if (sum != factor*atom->nangles)
    error->all(FLERR,"Angles assigned incorrectly");
}

/* ---------------------------------------------------------------------- */

void ReadData::dihedrals()
{
  int i,m,nchunk,eof;

  bigint nread = 0;
  bigint ndihedrals = atom->ndihedrals;

  while (nread < ndihedrals) {
    nchunk = MIN(ndihedrals-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_dihedrals(nchunk,buffer);
    nread += nchunk;
  }

  // check that dihedrals were assigned correctly

  int nlocal = atom->nlocal;
  bigint sum;
  bigint n = 0;
  for (i = 0; i < nlocal; i++) n += atom->num_dihedral[i];
  MPI_Allreduce(&n,&sum,1,MPI_LMP_BIGINT,MPI_SUM,world);
  int factor = 1;
  if (!force->newton_bond) factor = 4;

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " dihedrals\n",sum/factor);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " dihedrals\n",sum/factor);
  }
  if (sum != factor*atom->ndihedrals)
    error->all(FLERR,"Dihedrals assigned incorrectly");
}

/* ---------------------------------------------------------------------- */

void ReadData::impropers()
{
  int i,m,nchunk,eof;

  bigint nread = 0;
  bigint nimpropers = atom->nimpropers;

  while (nread < nimpropers) {
    nchunk = MIN(nimpropers-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_impropers(nchunk,buffer);
    nread += nchunk;
  }

  // check that impropers were assigned correctly

  int nlocal = atom->nlocal;
  bigint sum;
  bigint n = 0;
  for (i = 0; i < nlocal; i++) n += atom->num_improper[i];
  MPI_Allreduce(&n,&sum,1,MPI_LMP_BIGINT,MPI_SUM,world);
  int factor = 1;
  if (!force->newton_bond) factor = 4;

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " impropers\n",sum/factor);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " impropers\n",sum/factor);
  }
  if (sum != factor*atom->nimpropers)
    error->all(FLERR,"Impropers assigned incorrectly");
}

/* ---------------------------------------------------------------------- */

void ReadData::mass()
{
  int i,m;
  char *next;
  char *buf = new char[atom->ntypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,atom->ntypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (i = 0; i < atom->ntypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    atom->set_mass(buf);
    buf += strlen(buf) + 1;
  }
  delete [] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::paircoeffs()
{
  int i,m;
  char *next;
  char *buf = new char[atom->ntypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,atom->ntypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (i = 0; i < atom->ntypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    m = strlen(buf) + 1;
    parse_coeffs(buf,NULL,1);
    force->pair->coeff(narg,arg);
    buf += m;
  }
  delete [] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::pairIJcoeffs()
{
  int i,j,m;
  char *next;
  
  int nsq = atom->ntypes* (atom->ntypes+1) / 2;
  char *buf = new char[nsq * MAXLINE];

  int eof = comm->read_lines_from_file(fp,nsq,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (i = 0; i < atom->ntypes; i++)
    for (j = i; j < atom->ntypes; j++) {
      next = strchr(buf,'\n');
      *next = '\0';
      m = strlen(buf) + 1;
      parse_coeffs(buf,NULL,0);
      force->pair->coeff(narg,arg);
      buf += m;
    }
  delete [] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::bondcoeffs()
{
  int i,m;
  char *next;
  char *buf = new char[atom->nbondtypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,atom->nbondtypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (i = 0; i < atom->nbondtypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    m = strlen(buf) + 1;
    parse_coeffs(buf,NULL,0);
    force->bond->coeff(narg,arg);
    buf += m;
  }
  delete [] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::anglecoeffs(int which)
{
  int i,m;
  char *next;
  char *buf = new char[atom->nangletypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,atom->nangletypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (i = 0; i < atom->nangletypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    m = strlen(buf) + 1;
    if (which == 0) parse_coeffs(buf,NULL,0);
    else if (which == 1) parse_coeffs(buf,"bb",0);
    else if (which == 2) parse_coeffs(buf,"ba",0);
    force->angle->coeff(narg,arg);
    buf += m;
  }
  delete [] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::dihedralcoeffs(int which)
{
  int i,m;
  char *next;
  char *buf = new char[atom->ndihedraltypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,atom->ndihedraltypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (i = 0; i < atom->ndihedraltypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    m = strlen(buf) + 1;
    if (which == 0) parse_coeffs(buf,NULL,0);
    else if (which == 1) parse_coeffs(buf,"mbt",0);
    else if (which == 2) parse_coeffs(buf,"ebt",0);
    else if (which == 3) parse_coeffs(buf,"at",0);
    else if (which == 4) parse_coeffs(buf,"aat",0);
    else if (which == 5) parse_coeffs(buf,"bb13",0);
    force->dihedral->coeff(narg,arg);
    buf += m;
  }
  delete [] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::impropercoeffs(int which)
{
  int i,m;
  char *next;
  char *buf = new char[atom->nimpropertypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,atom->nimpropertypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (i = 0; i < atom->nimpropertypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    m = strlen(buf) + 1;
    if (which == 0) parse_coeffs(buf,NULL,0);
    else if (which == 1) parse_coeffs(buf,"aa",0);
    force->improper->coeff(narg,arg);
    buf += m;
  }
  delete [] original;
}

/* ----------------------------------------------------------------------
   read fix section, pass lines to fix to process
   n = index of fix
------------------------------------------------------------------------- */

void ReadData::fix(int ifix, char *line, bigint nlines)
{
  int i,m,nchunk,eof;

  bigint nread = 0;
  while (nread < nlines) {
    nchunk = MIN(nlines-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    modify->fix[ifix]->read_data_section(line,nchunk,buffer);
    nread += nchunk;
  }
}

/* ----------------------------------------------------------------------
   proc 0 scans the data file for topology maximums
------------------------------------------------------------------------- */

void ReadData::scan(int &bond_per_atom, int &angle_per_atom,
                    int &dihedral_per_atom, int &improper_per_atom)
{
  int i,tmp1,tmp2,atom1,atom2,atom3,atom4;
  char *eof;

  if (atom->natoms > MAXSMALLINT)
    error->one(FLERR,"Molecular data file has too many atoms");

  // customize for new sections

  int natoms = static_cast<int> (atom->natoms);
  bond_per_atom = angle_per_atom = dihedral_per_atom = improper_per_atom = 0;
  int ellipsoid_flag = 0;
  int line_flag = 0;
  int tri_flag = 0;
  int body_flag = 0;

  // customize for new sections
  // allocate topology counting vector
  // initially, array length = 1 to natoms
  // will grow via reallocate() if atom IDs > natoms

  int cmax = natoms + 1;
  int *count;
  memory->create(count,cmax,"read_data:count");

  while (strlen(keyword)) {

    // allow special fixes first chance to match and process the section
    // if fix matches, continue to next section

    if (nfix) {
      for (i = 0; i < nfix; i++) {
        printf("LINE SECTION %s %s\n",line,fix_section[i]);
        if (strstr(line,fix_section[i])) {
          int n = modify->fix[fix_index[i]]->read_data_skip_lines(keyword);
          printf("NLINES SKIP %d\n",n);
          skip_lines(n);
          parse_keyword(0,0);
          break;
        }
      }
      if (i < nfix) continue;
    }

    if (strcmp(keyword,"Masses") == 0) skip_lines(atom->ntypes);
    else if (strcmp(keyword,"Atoms") == 0) skip_lines(natoms);
    else if (strcmp(keyword,"Velocities") == 0) skip_lines(natoms);

    else if (strcmp(keyword,"Ellipsoids") == 0) {
      if (!avec_ellipsoid)
        error->one(FLERR,"Invalid data file section: Ellipsoids");
      ellipsoid_flag = 1;
      skip_lines(nellipsoids);
    } else if (strcmp(keyword,"Lines") == 0) {
      if (!avec_line) error->one(FLERR,"Invalid data file section: Lines");
      line_flag = 1;
      skip_lines(nlines);
    } else if (strcmp(keyword,"Triangles") == 0) {
      if (!avec_tri) error->one(FLERR,"Invalid data file section: Triangles");
      tri_flag = 1;
      skip_lines(ntris);
    } else if (strcmp(keyword,"Bodies") == 0) {
      if (!avec_body) error->one(FLERR,"Invalid data file section: Bodies");
      body_flag = 1;
      skip_lines(nbodies);

    } else if (strcmp(keyword,"Pair Coeffs") == 0) {
      if (force->pair == NULL)
        error->one(FLERR,"Must define pair_style before Pair Coeffs");
      skip_lines(atom->ntypes);
    } else if (strcmp(keyword,"PairIJ Coeffs") == 0) {
      if (force->pair == NULL)
        error->one(FLERR,"Must define pair_style before Pair Coeffs");
      skip_lines(atom->ntypes*(atom->ntypes+1)/2);
    } else if (strcmp(keyword,"Bond Coeffs") == 0) {
      if (atom->avec->bonds_allow == 0)
        error->one(FLERR,"Invalid data file section: Bond Coeffs");
      if (force->bond == NULL)
        error->one(FLERR,"Must define bond_style before Bond Coeffs");
      skip_lines(atom->nbondtypes);
    } else if (strcmp(keyword,"Angle Coeffs") == 0) {
      if (atom->avec->angles_allow == 0)
        error->one(FLERR,"Invalid data file section: Angle Coeffs");
      if (force->angle == NULL)
        error->one(FLERR,"Must define angle_style before Angle Coeffs");
      skip_lines(atom->nangletypes);
    } else if (strcmp(keyword,"Dihedral Coeffs") == 0) {
      skip_lines(atom->ndihedraltypes);
      if (atom->avec->dihedrals_allow == 0)
        error->one(FLERR,"Invalid data file section: Dihedral Coeffs");
      if (force->dihedral == NULL)
        error->one(FLERR,"Must define dihedral_style before Dihedral Coeffs");
    }  else if (strcmp(keyword,"Improper Coeffs") == 0) {
      if (atom->avec->impropers_allow == 0)
        error->one(FLERR,"Invalid data file section: Improper Coeffs");
      if (force->improper == NULL)
        error->one(FLERR,"Must define improper_style before Improper Coeffs");
      skip_lines(atom->nimpropertypes);

    } else if (strcmp(keyword,"BondBond Coeffs") == 0) {
      if (atom->avec->angles_allow == 0)
        error->one(FLERR,"Invalid data file section: BondBond Coeffs");
      if (force->angle == NULL)
        error->one(FLERR,"Must define angle_style before BondBond Coeffs");
      skip_lines(atom->nangletypes);
    } else if (strcmp(keyword,"BondAngle Coeffs") == 0) {
      if (atom->avec->angles_allow == 0)
        error->one(FLERR,"Invalid data file section: BondAngle Coeffs");
      if (force->angle == NULL)
        error->one(FLERR,"Must define angle_style before BondAngle Coeffs");
      skip_lines(atom->nangletypes);
    } else if (strcmp(keyword,"MiddleBondTorsion Coeffs") == 0) {
      if (atom->avec->dihedrals_allow == 0)
        error->one(FLERR,"Invalid data file section: MiddleBondTorsion Coeffs");
      if (force->dihedral == NULL)
        error->one(FLERR,
                   "Must define dihedral_style before "
                   "MiddleBondTorsion Coeffs");
      skip_lines(atom->ndihedraltypes);
    } else if (strcmp(keyword,"EndBondTorsion Coeffs") == 0) {
      if (atom->avec->dihedrals_allow == 0)
        error->one(FLERR,"Invalid data file section: EndBondTorsion Coeffs");
      if (force->dihedral == NULL)
        error->one(FLERR,
                   "Must define dihedral_style before EndBondTorsion Coeffs");
      skip_lines(atom->ndihedraltypes);
    } else if (strcmp(keyword,"AngleTorsion Coeffs") == 0) {
      if (atom->avec->dihedrals_allow == 0)
        error->one(FLERR,"Invalid data file section: AngleTorsion Coeffs");
      if (force->dihedral == NULL)
        error->one(FLERR,
                   "Must define dihedral_style before AngleTorsion Coeffs");
      skip_lines(atom->ndihedraltypes);
    } else if (strcmp(keyword,"AngleAngleTorsion Coeffs") == 0) {
      if (atom->avec->dihedrals_allow == 0)
        error->one(FLERR,"Invalid data file section: AngleAngleTorsion Coeffs");
      if (force->dihedral == NULL)
        error->one(FLERR,
                   "Must define dihedral_style before "
                   "AngleAngleTorsion Coeffs");
      skip_lines(atom->ndihedraltypes);
    } else if (strcmp(keyword,"BondBond13 Coeffs") == 0) {
      if (atom->avec->dihedrals_allow == 0)
        error->one(FLERR,"Invalid data file section: BondBond13 Coeffs");
      if (force->dihedral == NULL)
        error->one(FLERR,"Must define dihedral_style before BondBond13 Coeffs");
      skip_lines(atom->ndihedraltypes);
    } else if (strcmp(keyword,"AngleAngle Coeffs") == 0) {
      if (atom->avec->impropers_allow == 0)
        error->one(FLERR,"Invalid data file section: AngleAngle Coeffs");
      if (force->improper == NULL)
        error->one(FLERR,"Must define improper_style before AngleAngle Coeffs");
      skip_lines(atom->nimpropertypes);

    } else if (strcmp(keyword,"Bonds") == 0) {
      for (i = 1; i < cmax; i++) count[i] = 0;
      if (force->newton_bond)
        for (i = 0; i < atom->nbonds; i++) {
          eof = fgets(line,MAXLINE,fp);
          if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
          sscanf(line,"%d %d %d %d",&tmp1,&tmp2,&atom1,&atom2);
          if (atom1 >= cmax) cmax = reallocate(&count,cmax,atom1);
          count[atom1]++;
        }
      else
        for (i = 0; i < atom->nbonds; i++) {
          eof = fgets(line,MAXLINE,fp);
          if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
          sscanf(line,"%d %d %d %d",&tmp1,&tmp2,&atom1,&atom2);
          int amax = MAX(atom1,atom2);
          if (amax >= cmax) cmax = reallocate(&count,cmax,amax);
          count[atom1]++;
          count[atom2]++;
        }
      for (i = 1; i < cmax; i++) bond_per_atom = MAX(bond_per_atom,count[i]);
      if (screen) fprintf(screen,"  %d = max bonds/atom\n",bond_per_atom);
      if (logfile) fprintf(logfile,"  %d = max bonds/atom\n",bond_per_atom);

    } else if (strcmp(keyword,"Angles") == 0) {
      for (i = 1; i < cmax; i++) count[i] = 0;
      if (force->newton_bond)
        for (i = 0; i < atom->nangles; i++) {
          eof = fgets(line,MAXLINE,fp);
          if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
          sscanf(line,"%d %d %d %d %d",&tmp1,&tmp2,&atom1,&atom2,&atom3);
          if (atom2 >= cmax) cmax = reallocate(&count,cmax,atom2);
          count[atom2]++;
        }
      else
        for (i = 0; i < atom->nangles; i++) {
          eof = fgets(line,MAXLINE,fp);
          if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
          sscanf(line,"%d %d %d %d %d",&tmp1,&tmp2,&atom1,&atom2,&atom3);
          int amax = MAX(atom1,atom2);
          amax = MAX(amax,atom3);
          if (amax >= cmax) cmax = reallocate(&count,cmax,amax);
          count[atom1]++;
          count[atom2]++;
          count[atom3]++;
        }
      for (i = 1; i < cmax; i++) angle_per_atom = MAX(angle_per_atom,count[i]);
      if (screen) fprintf(screen,"  %d = max angles/atom\n",angle_per_atom);
      if (logfile) fprintf(logfile,"  %d = max angles/atom\n",angle_per_atom);

    } else if (strcmp(keyword,"Dihedrals") == 0) {
      for (i = 1; i < cmax; i++) count[i] = 0;
      if (force->newton_bond)
        for (i = 0; i < atom->ndihedrals; i++) {
          eof = fgets(line,MAXLINE,fp);
          if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
          sscanf(line,"%d %d %d %d %d %d",
                 &tmp1,&tmp2,&atom1,&atom2,&atom3,&atom4);
          if (atom2 >= cmax) cmax = reallocate(&count,cmax,atom2);
          count[atom2]++;
        }
      else
        for (i = 0; i < atom->ndihedrals; i++) {
          eof = fgets(line,MAXLINE,fp);
          if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
          sscanf(line,"%d %d %d %d %d %d",
                 &tmp1,&tmp2,&atom1,&atom2,&atom3,&atom4);
          int amax = MAX(atom1,atom2);
          amax = MAX(amax,atom3);
          amax = MAX(amax,atom4);
          if (amax >= cmax) cmax = reallocate(&count,cmax,amax);
          count[atom1]++;
          count[atom2]++;
          count[atom3]++;
          count[atom4]++;
        }
      for (i = 1; i < cmax; i++)
        dihedral_per_atom = MAX(dihedral_per_atom,count[i]);
      if (screen)
        fprintf(screen,"  %d = max dihedrals/atom\n",dihedral_per_atom);
      if (logfile)
        fprintf(logfile,"  %d = max dihedrals/atom\n",dihedral_per_atom);

    } else if (strcmp(keyword,"Impropers") == 0) {
      for (i = 1; i < cmax; i++) count[i] = 0;
      if (force->newton_bond)
        for (i = 0; i < atom->nimpropers; i++) {
          eof = fgets(line,MAXLINE,fp);
          if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
          sscanf(line,"%d %d %d %d %d %d",
                 &tmp1,&tmp2,&atom1,&atom2,&atom3,&atom4);
          if (atom2 >= cmax) cmax = reallocate(&count,cmax,atom2);
          count[atom2]++;
        }
      else
        for (i = 0; i < atom->nimpropers; i++) {
          eof = fgets(line,MAXLINE,fp);
          if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
          sscanf(line,"%d %d %d %d %d %d",
                 &tmp1,&tmp2,&atom1,&atom2,&atom3,&atom4);
          int amax = MAX(atom1,atom2);
          amax = MAX(amax,atom3);
          amax = MAX(amax,atom4);
          if (amax >= cmax) cmax = reallocate(&count,cmax,amax);
          count[atom1]++;
          count[atom2]++;
          count[atom3]++;
          count[atom4]++;
        }
      for (i = 1; i < cmax; i++)
        improper_per_atom = MAX(improper_per_atom,count[i]);
      if (screen)
        fprintf(screen,"  %d = max impropers/atom\n",improper_per_atom);
      if (logfile)
        fprintf(logfile,"  %d = max impropers/atom\n",improper_per_atom);

    } else {
      char str[128];
      sprintf(str,"Unknown identifier in data file: %s",keyword);
      error->one(FLERR,str);
    }

    parse_keyword(0,0);
  }

  // free topology counting vector

  memory->destroy(count);

  // error check that topology was specified in file

  if ((atom->nbonds && !bond_per_atom) ||
      (atom->nangles && !angle_per_atom) ||
      (atom->ndihedrals && !dihedral_per_atom) ||
      (atom->nimpropers && !improper_per_atom))
    error->one(FLERR,"Needed topology not in data file");

  // customize for new sections
  // error check that Bonus sections were speficied in file

  if (nellipsoids && !ellipsoid_flag)
    error->one(FLERR,"Needed bonus data not in data file");
  if (nlines && !line_flag)
    error->one(FLERR,"Needed bonus data not in data file");
  if (ntris && !tri_flag)
    error->one(FLERR,"Needed bonus data not in data file");
  if (nbodies && !body_flag)
    error->one(FLERR,"Needed bonus data not in data file");
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
    sprintf(gunzip,"gunzip -c %s",file);
    fp = popen(gunzip,"r");
#else
    error->one(FLERR,"Cannot open gzipped file");
#endif
  }

  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",file);
    error->one(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   grab next keyword
   read lines until one is non-blank
   keyword is all text on line w/out leading & trailing white space
   read one additional line (assumed blank)
   if any read hits EOF, set keyword to empty
   if first = 1, line variable holds non-blank line that ended header
   if flag = 0, only proc 0 is calling so no bcast
   else flag = 1, bcast keyword line to all procs
------------------------------------------------------------------------- */

void ReadData::parse_keyword(int first, int flag)
{
  int eof = 0;

  // proc 0 reads upto non-blank line plus 1 following line
  // eof is set to 1 if any read hits end-of-file

  if (me == 0) {
    if (!first) {
      if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
    }
    while (eof == 0 && strspn(line," \t\n\r") == strlen(line)) {
      if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
    }
    if (fgets(buffer,MAXLINE,fp) == NULL) eof = 1;
  }

  // if eof, set keyword empty and return

  if (flag) MPI_Bcast(&eof,1,MPI_INT,0,world);
  if (eof) {
    keyword[0] = '\0';
    return;
  }

  // bcast keyword line to all procs

  if (flag) {
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
   proc 0 reads N lines from file
------------------------------------------------------------------------- */

void ReadData::skip_lines(int n)
{
  char *eof;
  for (int i = 0; i < n; i++) eof = fgets(line,MAXLINE,fp);
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
------------------------------------------------------------------------- */

void ReadData::parse_coeffs(char *line, const char *addstr, int dupflag)
{
  char *ptr;
  if (ptr = strchr(line,'#')) *ptr = '\0';

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
}
