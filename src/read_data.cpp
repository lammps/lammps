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
#include "force.h"
#include "molecule.h"
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
  keyword = new char[MAXLINE];
  style = new char[MAXLINE];
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

  // optional args

  addflag = mergeflag = 0;
  offset[0] = offset[1] = offset[2] = 0.0;
  nfix = 0;
  fix_index = NULL;
  fix_header = NULL;
  fix_section = NULL;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"add") == 0) {
      addflag = 1;
      iarg++;
    } else if (strcmp(arg[iarg],"merge") == 0) {
      mergeflag = 1;
      iarg++;
    } else if (strcmp(arg[iarg],"offset") == 0) {
      if (iarg+4 > narg)
        error->all(FLERR,"Illegal read_data command");
      offset[0] = force->numeric(FLERR,arg[iarg+1]);
      offset[1] = force->numeric(FLERR,arg[iarg+2]);
      offset[2] = force->numeric(FLERR,arg[iarg+3]);
      iarg += 4;
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

  if (domain->box_exist && !addflag && !mergeflag)
    error->all(FLERR,"Cannot read_data after simulation box is defined");
  if (addflag && mergeflag) error->all(FLERR,"Cannot read_data add and merge");
  if (domain->dimension == 2 && offset[2] != 0.0)
    error->all(FLERR,"Cannot use non-zero z offset in read_data "
               "for 2d simulation");

  if (domain->dimension == 2 && domain->zperiodic == 0)
    error->all(FLERR,"Cannot run 2d simulation with nonperiodic Z dimension");

  // perform 1-pass read if no molecular topoogy in file
  // perform 2-pass read if molecular topology,
  //   first pass calculates max topology/atom

  int atomflag,topoflag;
  int bondflag,angleflag,dihedralflag,improperflag;
  int ellipsoidflag,lineflag,triflag,bodyflag;

  atomflag = topoflag = 0;
  bondflag = angleflag = dihedralflag = improperflag = 0;
  ellipsoidflag = lineflag = triflag = bodyflag = 0;

  int firstpass = 1;

  while (1) {

    // open file on proc 0

    if (me == 0) {
      if (firstpass && screen) fprintf(screen,"Reading data file ...\n");
      open(arg[0]);
    } else fp = NULL;
    
    // read header info
    
    header();

    // problem setup using info from header
    // 1st pass only

    if (firstpass) {
      domain->box_exist = 1;
      update->ntimestep = 0;
    
      // insure extra settings are applied before grow(),
      //   even if no topology in file
      // if topology is in file, realloc and another grow() is done below

      atom->bond_per_atom = atom->extra_bond_per_atom;
      atom->angle_per_atom = atom->extra_angle_per_atom;
      atom->dihedral_per_atom = atom->extra_dihedral_per_atom;
      atom->improper_per_atom = atom->extra_improper_per_atom;

      int n;
      if (comm->nprocs == 1) n = static_cast<int> (atom->natoms);
      else n = static_cast<int> (LB_FACTOR * atom->natoms / comm->nprocs);

      atom->allocate_type_arrays();
      atom->avec->grow(n);

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
        } else skip_lines(atom->natoms);
      } else if (strcmp(keyword,"Velocities") == 0) {
        if (atomflag == 0) 
          error->all(FLERR,"Must read Atoms before Velocities");
        if (firstpass) velocities();
        else skip_lines(atom->natoms);

      } else if (strcmp(keyword,"Bonds") == 0) {
        topoflag = bondflag = 1;
        if (atom->nbonds == 0)
          error->all(FLERR,"Invalid data file section: Bonds");
        if (atomflag == 0) error->all(FLERR,"Must read Atoms before Bonds");
        bonds(firstpass);
      } else if (strcmp(keyword,"Angles") == 0) {
        topoflag = angleflag = 1;
        if (atom->nangles == 0)
          error->all(FLERR,"Invalid data file section: Angles");
        if (atomflag == 0) error->all(FLERR,"Must read Atoms before Angles");
        angles(firstpass);
      } else if (strcmp(keyword,"Dihedrals") == 0) {
        topoflag = dihedralflag = 1;
        if (atom->ndihedrals == 0)
          error->all(FLERR,"Invalid data file section: Dihedrals");
        if (atomflag == 0) error->all(FLERR,"Must read Atoms before Dihedrals");
        dihedrals(firstpass);
      } else if (strcmp(keyword,"Impropers") == 0) {
        topoflag = improperflag = 1;
        if (atom->nimpropers == 0)
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
        else skip_lines(atom->ntypes);
      } else if (strcmp(keyword,"Pair Coeffs") == 0) {
        if (force->pair == NULL)
          error->all(FLERR,"Must define pair_style before Pair Coeffs");
        if (firstpass) {
          if (me == 0 && !style_match(style,force->pair_style))
            error->warning(FLERR,"Pair style in data file differs "
                           "from currently defined pair style");
          paircoeffs();
        } else skip_lines(atom->ntypes);
      } else if (strcmp(keyword,"PairIJ Coeffs") == 0) {
        if (force->pair == NULL)
          error->all(FLERR,"Must define pair_style before PairIJ Coeffs");
        if (firstpass) {
          if (me == 0 && !style_match(style,force->pair_style))
            error->warning(FLERR,"Pair style in data file differs "
                           "from currently defined pair style");
          pairIJcoeffs();
        } else skip_lines(atom->ntypes*(atom->ntypes+1)/2);
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
        } else skip_lines(atom->nbondtypes);
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
        } else skip_lines(atom->nangletypes);
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
        } else skip_lines(atom->ndihedraltypes);
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
        } else skip_lines(atom->nimpropertypes);

      } else if (strcmp(keyword,"BondBond Coeffs") == 0) {
        if (atom->avec->angles_allow == 0)
          error->all(FLERR,"Invalid data file section: BondBond Coeffs");
        if (force->angle == NULL)
          error->all(FLERR,"Must define angle_style before BondBond Coeffs");
        if (firstpass) anglecoeffs(1);
        else skip_lines(atom->nangletypes);
      } else if (strcmp(keyword,"BondAngle Coeffs") == 0) {
        if (atom->avec->angles_allow == 0)
          error->all(FLERR,"Invalid data file section: BondAngle Coeffs");
        if (force->angle == NULL)
          error->all(FLERR,"Must define angle_style before BondAngle Coeffs");
        if (firstpass) anglecoeffs(2);
        else skip_lines(atom->nangletypes);

      } else if (strcmp(keyword,"MiddleBondTorsion Coeffs") == 0) {
        if (atom->avec->dihedrals_allow == 0)
          error->all(FLERR,
                     "Invalid data file section: MiddleBondTorsion Coeffs");
        if (force->dihedral == NULL)
          error->all(FLERR,
                     "Must define dihedral_style before "
                     "MiddleBondTorsion Coeffs");
        if (firstpass) dihedralcoeffs(1);
        else skip_lines(atom->ndihedraltypes);
      } else if (strcmp(keyword,"EndBondTorsion Coeffs") == 0) {
        if (atom->avec->dihedrals_allow == 0)
          error->all(FLERR,"Invalid data file section: EndBondTorsion Coeffs");
        if (force->dihedral == NULL)
          error->all(FLERR,
                     "Must define dihedral_style before EndBondTorsion Coeffs");
        if (firstpass) dihedralcoeffs(2);
        else skip_lines(atom->ndihedraltypes);
      } else if (strcmp(keyword,"AngleTorsion Coeffs") == 0) {
        if (atom->avec->dihedrals_allow == 0)
          error->all(FLERR,"Invalid data file section: AngleTorsion Coeffs");
        if (force->dihedral == NULL)
          error->all(FLERR,
                     "Must define dihedral_style before AngleTorsion Coeffs");
        if (firstpass) dihedralcoeffs(3);
        else skip_lines(atom->ndihedraltypes);
      } else if (strcmp(keyword,"AngleAngleTorsion Coeffs") == 0) {
        if (atom->avec->dihedrals_allow == 0)
          error->all(FLERR,
                     "Invalid data file section: AngleAngleTorsion Coeffs");
        if (force->dihedral == NULL)
          error->all(FLERR,
                     "Must define dihedral_style before "
                     "AngleAngleTorsion Coeffs");
        if (firstpass) dihedralcoeffs(4);
        else skip_lines(atom->ndihedraltypes);
      } else if (strcmp(keyword,"BondBond13 Coeffs") == 0) {
        if (atom->avec->dihedrals_allow == 0)
          error->all(FLERR,"Invalid data file section: BondBond13 Coeffs");
        if (force->dihedral == NULL)
          error->all(FLERR,
                     "Must define dihedral_style before BondBond13 Coeffs");
        if (firstpass) dihedralcoeffs(5);
        else skip_lines(atom->ndihedraltypes);
        
      } else if (strcmp(keyword,"AngleAngle Coeffs") == 0) {
        if (atom->avec->impropers_allow == 0)
          error->all(FLERR,"Invalid data file section: AngleAngle Coeffs");
        if (force->improper == NULL)
          error->all(FLERR,
                     "Must define improper_style before AngleAngle Coeffs");
        if (firstpass) impropercoeffs(1);
        else skip_lines(atom->nimpropertypes);

      } else {
        char str[128];
        sprintf(str,"Unknown identifier in data file: %s",keyword);
        error->all(FLERR,str);
      }
      
      parse_keyword(0);
    }

    // error if natoms > 0 yet no atoms were read

    if (atom->natoms > 0 && atomflag == 0)
      error->all(FLERR,"No atoms in data file");
    
    // close file
    
    if (me == 0) {
      if (compressed) pclose(fp);
      else fclose(fp);
    }

    // done if this was 2nd pass

    if (!firstpass) break;

    // at end of 1st pass, error check for required sections
    // customize for new sections

    if ((atom->nbonds && !bondflag) || (atom->nangles && !angleflag) ||
        (atom->ndihedrals && !dihedralflag) ||
        (atom->nimpropers && !improperflag))
      error->one(FLERR,"Needed molecular topology not in data file");

    if ((nellipsoids && !ellipsoidflag) || (nlines && !lineflag) ||
        (ntris && !triflag) || (nbodies && !bodyflag))
      error->one(FLERR,"Needed bonus data not in data file");

    // break out of loop if no molecular topology in file
    // else make 2nd pass
    
    if (!topoflag) break;
    firstpass = 0;

    // reallocate bond,angle,diehdral,improper arrays via grow()
    // use new bond,angle,dihedral,improper per-atom values from 1st pass
    // should leave other atom arrays unchanged, since already nmax in length
    // if bonds/etc not in data file, initialize per-atom size 
    //   with extra settings before grow() of these topology arrays

    if (bondflag) {
      memory->destroy(atom->bond_type);
      memory->destroy(atom->bond_atom);
      atom->bond_type = NULL;
      atom->bond_atom = NULL;
    }

    if (angleflag) {
      memory->destroy(atom->angle_type);
      memory->destroy(atom->angle_atom1);
      memory->destroy(atom->angle_atom2);
      memory->destroy(atom->angle_atom3);
      atom->angle_type = NULL;
      atom->angle_atom1 = atom->angle_atom2 = atom->angle_atom3 = NULL;
    }

    if (dihedralflag) {
      memory->destroy(atom->dihedral_type);
      memory->destroy(atom->dihedral_atom1);
      memory->destroy(atom->dihedral_atom2);
      memory->destroy(atom->dihedral_atom3);
      memory->destroy(atom->dihedral_atom4);
      atom->dihedral_type = NULL;
      atom->dihedral_atom1 = atom->dihedral_atom2 = 
        atom->dihedral_atom3 = atom->dihedral_atom4 = NULL;
    }

    if (improperflag) {
      memory->destroy(atom->improper_type);
      memory->destroy(atom->improper_atom1);
      memory->destroy(atom->improper_atom2);
      memory->destroy(atom->improper_atom3);
      memory->destroy(atom->improper_atom4);
      atom->improper_type = NULL;
      atom->improper_atom1 = atom->improper_atom2 = 
        atom->improper_atom3 = atom->improper_atom4 = NULL;
    }

    atom->avec->grow(atom->nmax);
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
}

/* ----------------------------------------------------------------------
   read free-format header of data file
   1st line and blank lines are skipped
   non-blank lines are checked for header keywords and leading value is read
   header ends with EOF or non-blank line containing no header keyword
     if EOF, line is set to blank line
     else line has first keyword line for rest of file
------------------------------------------------------------------------- */

void ReadData::header()
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

    // read a line and bcast length if flag is set

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

    if (strstr(line,"atoms")) {
      sscanf(line,BIGINT_FORMAT,&atom->natoms);

    // check for these first
    // otherwise "triangles" will be matched as "angles"

    } else if (strstr(line,"ellipsoids")) {
      if (!avec_ellipsoid)
        error->all(FLERR,"No ellipsoids allowed with this atom style");
      sscanf(line,BIGINT_FORMAT,&nellipsoids);
    } else if (strstr(line,"lines")) {
      if (!avec_line)
        error->all(FLERR,"No lines allowed with this atom style");
      sscanf(line,BIGINT_FORMAT,&nlines);
    } else if (strstr(line,"triangles")) {
      if (!avec_tri)
        error->all(FLERR,"No triangles allowed with this atom style");
      sscanf(line,BIGINT_FORMAT,&ntris);
    } else if (strstr(line,"bodies")) {
      if (!avec_body)
        error->all(FLERR,"No bodies allowed with this atom style");
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
    else if (strstr(line,"extra angle per atom"))
      sscanf(line,"%d",&atom->extra_angle_per_atom);
    else if (strstr(line,"extra dihedral per atom"))
      sscanf(line,"%d",&atom->extra_dihedral_per_atom);
    else if (strstr(line,"extra improper per atom"))
      sscanf(line,"%d",&atom->extra_improper_per_atom);
    else if (strstr(line,"extra special per atom"))
      sscanf(line,"%d",&force->special_extra);

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

  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT ||
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
  
  // check that atom IDs are valid

  atom->tag_check();

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
   scan or read all bonds
------------------------------------------------------------------------- */

void ReadData::bonds(int firstpass)
{
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
    for (int i = 0; i < nlocal; i++) count[i] = 0;
  }

  // read and process bonds

  int nchunk,eof;
  bigint nread = 0;
  bigint nbonds = atom->nbonds;

  while (nread < nbonds) {
    nchunk = MIN(nbonds-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_bonds(nchunk,buffer,count);
    nread += nchunk;
  }

  // if firstpass: tally max bond/atom and return

  if (firstpass) {
    int max = 0;
    for (int i = 0; i < nlocal; i++) max = MAX(max,count[i]);
    int maxall;
    MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);

    if (me == 0) {
      if (screen) fprintf(screen,"  %d = max bonds/atom\n",maxall);
      if (logfile) fprintf(logfile,"  %d = max bonds/atom\n",maxall);
    }
    atom->bond_per_atom = maxall + atom->extra_bond_per_atom;
    memory->destroy(count);
    return;
  }

  // if 2nd pass: check that bonds were assigned correctly

  bigint n = 0;
  for (int i = 0; i < nlocal; i++) n += atom->num_bond[i];
  bigint sum;
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

/* ----------------------------------------------------------------------
   scan or read all angles
------------------------------------------------------------------------- */

void ReadData::angles(int firstpass)
{
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
    for (int i = 0; i < nlocal; i++) count[i] = 0;
  }

  // read and process angles

  int nchunk,eof;
  bigint nread = 0;
  bigint nangles = atom->nangles;

  while (nread < nangles) {
    nchunk = MIN(nangles-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_angles(nchunk,buffer,count);
    nread += nchunk;
  }

  // if firstpass: tally max angle/atom and return

  if (firstpass) {
    int max = 0;
    for (int i = 0; i < nlocal; i++) max = MAX(max,count[i]);
    int maxall;
    MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);

    if (me == 0) {
      if (screen) fprintf(screen,"  %d = max angles/atom\n",maxall);
      if (logfile) fprintf(logfile,"  %d = max angles/atom\n",maxall);
    }
    atom->angle_per_atom = maxall + atom->extra_angle_per_atom;
    memory->destroy(count);
    return;
  }

  // if 2nd pass: check that angles were assigned correctly

  bigint n = 0;
  for (int i = 0; i < nlocal; i++) n += atom->num_angle[i];
  bigint sum;
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

/* ----------------------------------------------------------------------
   scan or read all dihedrals
------------------------------------------------------------------------- */

void ReadData::dihedrals(int firstpass)
{
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
    for (int i = 0; i < nlocal; i++) count[i] = 0;
  }

  // read and process dihedrals

  int nchunk,eof;
  bigint nread = 0;
  bigint ndihedrals = atom->ndihedrals;

  while (nread < ndihedrals) {
    nchunk = MIN(ndihedrals-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_dihedrals(nchunk,buffer,count);
    nread += nchunk;
  }

  // if firstpass: tally max dihedral/atom and return

  if (firstpass) {
    int max = 0;
    for (int i = 0; i < nlocal; i++) max = MAX(max,count[i]);
    int maxall;
    MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);

    if (me == 0) {
      if (screen) fprintf(screen,"  %d = max dihedrals/atom\n",maxall);
      if (logfile) fprintf(logfile,"  %d = max dihedrals/atom\n",maxall);
    }
    atom->dihedral_per_atom = maxall + atom->extra_dihedral_per_atom;
    memory->destroy(count);
    return;
  }

  // if 2nd pass: check that dihedrals were assigned correctly

  bigint n = 0;
  for (int i = 0; i < nlocal; i++) n += atom->num_dihedral[i];
  bigint sum;
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

/* ----------------------------------------------------------------------
   scan or read all impropers
------------------------------------------------------------------------- */

void ReadData::impropers(int firstpass)
{
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
    for (int i = 0; i < nlocal; i++) count[i] = 0;
  }

  // read and process impropers

  int nchunk,eof;
  bigint nread = 0;
  bigint nimpropers = atom->nimpropers;

  while (nread < nimpropers) {
    nchunk = MIN(nimpropers-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_impropers(nchunk,buffer,count);
    nread += nchunk;
  }

  // if firstpass: tally max improper/atom and return

  if (firstpass) {
    int max = 0;
    for (int i = 0; i < nlocal; i++) max = MAX(max,count[i]);
    int maxall;
    MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);

    if (me == 0) {
      if (screen) fprintf(screen,"  %d = max impropers/atom\n",maxall);
      if (logfile) fprintf(logfile,"  %d = max impropers/atom\n",maxall);
    }
    atom->improper_per_atom = maxall + atom->extra_improper_per_atom;
    memory->destroy(count);
    return;
  }

  // if 2nd pass: check that impropers were assigned correctly

  bigint n = 0;
  for (int i = 0; i < nlocal; i++) n += atom->num_improper[i];
  bigint sum;
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
   if not firstpass, just read but no processing of data
------------------------------------------------------------------------- */

void ReadData::bodies(int firstpass)
{
  int i,m,nchunk,nline,nmax,ninteger,ndouble,tmp,onebody;
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
        nline += onebody+1;
      }

      if (buffer[m-1] != '\n') strcpy(&buffer[m++],"\n");
      m++;
    }

    MPI_Bcast(&nchunk,1,MPI_INT,0,world);
    MPI_Bcast(&m,1,MPI_INT,0,world);
    MPI_Bcast(buffer,m,MPI_CHAR,0,world);

    if (firstpass) atom->data_bodies(nchunk,buffer,avec_body);
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
  char *buf = new char[atom->ntypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,atom->ntypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < atom->ntypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    atom->set_mass(buf);
    buf = next + 1;
  }
  delete [] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::paircoeffs()
{
  char *next;
  char *buf = new char[atom->ntypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,atom->ntypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < atom->ntypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    parse_coeffs(buf,NULL,1);
    if (narg == 0) error->all(FLERR,"Unexpected end of PairCoeffs section");
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
  
  int nsq = atom->ntypes* (atom->ntypes+1) / 2;
  char *buf = new char[nsq * MAXLINE];

  int eof = comm->read_lines_from_file(fp,nsq,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (i = 0; i < atom->ntypes; i++)
    for (j = i; j < atom->ntypes; j++) {
      next = strchr(buf,'\n');
      *next = '\0';
      parse_coeffs(buf,NULL,0);
      if (narg == 0) error->all(FLERR,"Unexpected end of PairCoeffs section");
      force->pair->coeff(narg,arg);
      buf = next + 1;
    }
  delete [] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::bondcoeffs()
{
  char *next;
  char *buf = new char[atom->nbondtypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,atom->nbondtypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < atom->nbondtypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    parse_coeffs(buf,NULL,0);
    if (narg == 0) error->all(FLERR,"Unexpected end of BondCoeffs section");
    force->bond->coeff(narg,arg);
    buf = next + 1;
  }
  delete [] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::anglecoeffs(int which)
{
  char *next;
  char *buf = new char[atom->nangletypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,atom->nangletypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < atom->nangletypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    if (which == 0) parse_coeffs(buf,NULL,0);
    else if (which == 1) parse_coeffs(buf,"bb",0);
    else if (which == 2) parse_coeffs(buf,"ba",0);
    if (narg == 0) error->all(FLERR,"Unexpected end of AngleCoeffs section");
    force->angle->coeff(narg,arg);
    buf = next + 1;
  }
  delete [] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::dihedralcoeffs(int which)
{
  char *next;
  char *buf = new char[atom->ndihedraltypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,atom->ndihedraltypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < atom->ndihedraltypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    if (which == 0) parse_coeffs(buf,NULL,0);
    else if (which == 1) parse_coeffs(buf,"mbt",0);
    else if (which == 2) parse_coeffs(buf,"ebt",0);
    else if (which == 3) parse_coeffs(buf,"at",0);
    else if (which == 4) parse_coeffs(buf,"aat",0);
    else if (which == 5) parse_coeffs(buf,"bb13",0);
    if (narg == 0) error->all(FLERR,"Unexpected end of DihedralCoeffs section");
    force->dihedral->coeff(narg,arg);
    buf = next + 1;
  }
  delete [] original;
}

/* ---------------------------------------------------------------------- */

void ReadData::impropercoeffs(int which)
{
  char *next;
  char *buf = new char[atom->nimpropertypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,atom->nimpropertypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < atom->nimpropertypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    if (which == 0) parse_coeffs(buf,NULL,0);
    else if (which == 1) parse_coeffs(buf,"aa",0);
    if (narg == 0) error->all(FLERR,"Unexpected end of ImproperCoeffs section");
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
    modify->fix[ifix]->read_data_section(keyword,nchunk,buffer);
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
    sprintf(gunzip,"gzip -c -d %s",file);

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
    sprintf(str,"Cannot open file %s",file);
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
  char *eof;
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
------------------------------------------------------------------------- */

void ReadData::parse_coeffs(char *line, const char *addstr, int dupflag)
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
