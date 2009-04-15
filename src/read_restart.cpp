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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "sys/types.h"
#include "dirent.h"
#include "read_restart.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "group.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "special.h"
#include "universe.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

// same as write_restart.cpp

enum{VERSION,UNITS,NTIMESTEP,DIMENSION,NPROCS,PROCGRID_0,PROCGRID_1,PROCGRID_2,
       NEWTON_PAIR,NEWTON_BOND,XPERIODIC,YPERIODIC,ZPERIODIC,
       BOUNDARY_00,BOUNDARY_01,BOUNDARY_10,BOUNDARY_11,BOUNDARY_20,BOUNDARY_21,
       ATOM_STYLE,NATOMS,NTYPES,
       NBONDS,NBONDTYPES,BOND_PER_ATOM,
       NANGLES,NANGLETYPES,ANGLE_PER_ATOM,
       NDIHEDRALS,NDIHEDRALTYPES,DIHEDRAL_PER_ATOM,
       NIMPROPERS,NIMPROPERTYPES,IMPROPER_PER_ATOM,
       BOXLO_0,BOXHI_0,BOXLO_1,BOXHI_1,BOXLO_2,BOXHI_2,
       SPECIAL_LJ_1,SPECIAL_LJ_2,SPECIAL_LJ_3,
       SPECIAL_COUL_1,SPECIAL_COUL_2,SPECIAL_COUL_3,
       XY,XZ,YZ};
enum{MASS,SHAPE,DIPOLE};
enum{PAIR,BOND,ANGLE,DIHEDRAL,IMPROPER};

#define LB_FACTOR 1.1

/* ---------------------------------------------------------------------- */

ReadRestart::ReadRestart(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void ReadRestart::command(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal read_restart command");

  if (domain->box_exist) 
    error->all("Cannot read_restart after simulation box is defined");

  MPI_Comm_rank(world,&me);

  // if filename contains "*", search dir for latest restart file

  char *file = new char[strlen(arg[0]) + 16];
  if (strchr(arg[0],'*')) file_search(arg[0],file);
  else strcpy(file,arg[0]);

  // check if filename contains "%"

  int multiproc;
  if (strchr(file,'%')) multiproc = 1;
  else multiproc = 0;

  // open single restart file or base file for multiproc case

  if (me == 0) {
    if (screen) fprintf(screen,"Reading restart file ...\n");
    char *hfile;
    if (multiproc) {
      hfile = new char[strlen(file) + 16];
      char *ptr = strchr(file,'%');
      *ptr = '\0';
      sprintf(hfile,"%s%s%s",file,"base",ptr+1);
      *ptr = '%';
    } else hfile = file;
    fp = fopen(hfile,"rb");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open restart file %s",hfile);
      error->one(str);
    }
    if (multiproc) delete [] hfile;
  }

  // read header info and create atom style and simulation box

  header();
  domain->box_exist = 1;

  // problem setup using info from header

  int n;
  if (comm->nprocs == 1) n = static_cast<int> (atom->natoms);
  else n = static_cast<int> (LB_FACTOR * atom->natoms / comm->nprocs);

  atom->allocate_type_arrays();
  atom->avec->grow(n);
  n = atom->nmax;

  domain->print_box("  ");
  domain->set_initial_box();
  domain->set_global_box();
  comm->set_procs();
  domain->set_local_box();

  // read groups, ntype-length arrays, force field, fix info from file
  // nextra = max # of extra quantities stored with each atom

  group->read_restart(fp);
  type_arrays();
  force_fields();

  int nextra = modify->read_restart(fp);
  atom->nextra_store = nextra;
  atom->extra = memory->create_2d_double_array(n,nextra,"atom:extra");

  // if single file:
  //   proc 0 reads atoms from file, one chunk per proc (nprocs_file)
  // else if one file per proc:
  //   proc 0 reads chunks from series of files (nprocs_file)
  // proc 0 bcasts each chunk to other procs
  // each proc unpacks the atoms, saving ones in it's sub-domain
  // check for atom in sub-domain is different for orthogonal vs triclinic box

  AtomVec *avec = atom->avec;
  int triclinic = domain->triclinic;

  int maxbuf = 0;
  double *buf = NULL;

  double *x,lamda[3];
  double *coord,*sublo,*subhi;
  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  int m;
  char *perproc = new char[strlen(file) + 16];
  char *ptr = strchr(file,'%');

  for (int iproc = 0; iproc < nprocs_file; iproc++) {
    if (me == 0) {
      if (multiproc) {
	fclose(fp);
	*ptr = '\0';
	sprintf(perproc,"%s%d%s",file,iproc,ptr+1);
	*ptr = '%';
	fp = fopen(perproc,"rb");
	if (fp == NULL) {
	  char str[128];
	  sprintf(str,"Cannot open restart file %s",perproc);
	  error->one(str);
	}
      }
      fread(&n,sizeof(int),1,fp);
    }

    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n > maxbuf) {
      maxbuf = n;
      delete [] buf;
      buf = new double[maxbuf];
    }

    if (n > 0) {
      if (me == 0) fread(buf,sizeof(double),n,fp);
      MPI_Bcast(buf,n,MPI_DOUBLE,0,world);
    }

    m = 0;
    while (m < n) {
      x = &buf[m+1];
      if (triclinic) {
	domain->x2lamda(x,lamda);
	coord = lamda;
      } else coord = x;

      if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
	  coord[1] >= sublo[1] && coord[1] < subhi[1] &&
	  coord[2] >= sublo[2] && coord[2] < subhi[2])
	m += avec->unpack_restart(&buf[m]);
      else m += static_cast<int> (buf[m]);
    }
  }

  // close restart file and clean-up memory
  
  if (me == 0) fclose(fp);
  delete [] buf;
  delete [] file;
  delete [] perproc;

  // check that all atoms were assigned to procs

  double natoms;
  double rlocal = atom->nlocal;
  MPI_Allreduce(&rlocal,&natoms,1,MPI_DOUBLE,MPI_SUM,world);

  if (me == 0) {
    if (screen) fprintf(screen,"  %.15g atoms\n",natoms);
    if (logfile) fprintf(logfile,"  %.15g atoms\n",natoms);
  }

  if (natoms != atom->natoms) error->all("Did not assign all atoms correctly");

  if (me == 0) {
    if (atom->nbonds) {
      if (screen) fprintf(screen,"  %d bonds\n",atom->nbonds);
      if (logfile) fprintf(logfile,"  %d bonds\n",atom->nbonds);
    }
    if (atom->nangles) {
      if (screen) fprintf(screen,"  %d angles\n",atom->nangles);
      if (logfile) fprintf(logfile,"  %d angles\n",atom->nangles);
    }
    if (atom->ndihedrals) {
      if (screen) fprintf(screen,"  %d dihedrals\n",atom->ndihedrals);
      if (logfile) fprintf(logfile,"  %d dihedrals\n",atom->ndihedrals);
    }
    if (atom->nimpropers) {
      if (screen) fprintf(screen,"  %d impropers\n",atom->nimpropers);
      if (logfile) fprintf(logfile,"  %d impropers\n",atom->nimpropers);
    }
  }

  // check if tags are being used
  // create global mapping and bond topology now that system is defined

  int flag = 0;
  for (int i = 0; i < atom->nlocal; i++)
    if (atom->tag[i] > 0) flag = 1;
  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_MAX,world);
  if (flag_all == 0) atom->tag_enable = 0;

  if (atom->map_style) {
    atom->map_init();
    atom->map_set();
  }
  if (atom->molecular) {
    Special special(lmp);
    special.build();
  }
}

/* ----------------------------------------------------------------------
   search for all files matching infile which contains a "*"
   replace "*" with latest timestep value to create outfile name
   search dir referenced by initial pathname of file
   if infile also contains "%", need to use "base" when search directory
------------------------------------------------------------------------- */

void ReadRestart::file_search(char *infile, char *outfile)
{
  char *ptr;

  // separate infile into dir + filename

  char *dirname = new char[strlen(infile) + 1];
  char *filename = new char[strlen(infile) + 1];

  if (strchr(infile,'/')) {
    ptr = strrchr(infile,'/');
    *ptr = '\0';
    strcpy(dirname,infile);
    strcpy(filename,ptr+1);
    *ptr = '/';
  } else {
    strcpy(dirname,"./");
    strcpy(filename,infile);
  }

  // if filename contains "%" replace "%" with "base"

  char *pattern = new char[strlen(filename) + 16];

  if (ptr = strchr(filename,'%')) {
    *ptr = '\0';
    sprintf(pattern,"%s%s%s",filename,"base",ptr+1);
    *ptr = '%';
  } else strcpy(pattern,filename);

  // scan all files in directory, searching for files that match pattern
  // maxnum = largest int that matches "*"

  int n = strlen(pattern) + 16;
  char *begin = new char[n];
  char *middle = new char[n];
  char *end = new char[n];

  ptr = strchr(pattern,'*');
  *ptr = '\0';
  strcpy(begin,pattern);
  strcpy(end,ptr+1);
  int nbegin = strlen(begin);
  int maxnum = -1;

  if (me == 0) {
    struct dirent *ep;
    DIR *dp = opendir(dirname);
    if (dp == NULL) 
      error->one("Cannot open dir to search for restart file");
    while (ep = readdir(dp)) {
      if (strstr(ep->d_name,begin) != ep->d_name) continue;
      if ((ptr = strstr(&ep->d_name[nbegin],end)) == NULL) continue;
      if (strlen(end) == 0) ptr = ep->d_name + strlen(ep->d_name);
      *ptr = '\0';
      if (strlen(&ep->d_name[nbegin]) < n) {
	strcpy(middle,&ep->d_name[nbegin]);
	if (atoi(middle) > maxnum) maxnum = atoi(middle);
      }
    }
    closedir(dp);
    if (maxnum < 0) error->one("Found no restart file matching pattern");
  }

  // create outfile with maxint substituted for "*"
  // use original infile, not pattern, since need to retain "%" in filename

  ptr = strchr(infile,'*');
  *ptr = '\0';
  sprintf(outfile,"%s%d%s",infile,maxnum,ptr+1);
  *ptr = '*';

  // clean up

  delete [] dirname;
  delete [] filename;
  delete [] pattern;
  delete [] begin;
  delete [] middle;
  delete [] end;
}

/* ----------------------------------------------------------------------
   read header of restart file 
------------------------------------------------------------------------- */

void ReadRestart::header()
{
  int px,py,pz;
  int xperiodic,yperiodic,zperiodic;
  int boundary[3][2];

  // read flags and values until flag = -1

  int flag = read_int();
  while (flag >= 0) {

    // check restart file version, warn if different

    if (flag == VERSION) {
      char *version = read_char();
      if (strcmp(version,universe->version) != 0 && me == 0) {
	error->warning("Restart file version does not match LAMMPS version");
	if (screen) fprintf(screen,"  restart file = %s, LAMMPS = %s\n",
			    version,universe->version);
      }
      delete [] version;

      // reset unit_style only if different
      // so that timestep,neighbor-skin are not changed

    } else if (flag == UNITS) {
      char *style = read_char();
      if (strcmp(style,update->unit_style) != 0) update->set_units(style);
      delete [] style;

    } else if (flag == NTIMESTEP) {
      update->ntimestep = read_int();

      // set dimension from restart file

    } else if (flag == DIMENSION) {
      int dimension = read_int();
      domain->dimension = dimension;
      if (domain->dimension == 2 && domain->zperiodic == 0)
	error->all("Cannot run 2d simulation with nonperiodic Z dimension");

      // read nprocs from restart file, warn if different

    } else if (flag == NPROCS) {
      nprocs_file = read_int();
      if (nprocs_file != comm->nprocs && me == 0)
	error->warning("Restart file used different # of processors");

      // don't set procgrid, warn if different

    } else if (flag == PROCGRID_0) {
      px = read_int();
    } else if (flag == PROCGRID_1) {
      py = read_int();
    } else if (flag == PROCGRID_2) {
      pz = read_int();
      if (comm->user_procgrid[0] != 0 && 
	  (px != comm->user_procgrid[0] || py != comm->user_procgrid[1] || 
	   pz != comm->user_procgrid[2]) && me == 0)
	error->warning("Restart file used different 3d processor grid");

      // don't set newton_pair, leave input script value unchanged
      // set newton_bond from restart file
      // warn if different and input script settings are not default

    } else if (flag == NEWTON_PAIR) {
      int newton_pair_file = read_int();
      if (force->newton_pair != 1) {
	if (newton_pair_file != force->newton_pair && me == 0)
	  error->warning("Restart file used different newton pair setting, "
			 "using input script value");
      }
    } else if (flag == NEWTON_BOND) {
      int newton_bond_file = read_int();
      if (force->newton_bond != 1) {
	if (newton_bond_file != force->newton_bond && me == 0)
	  error->warning("Restart file used different newton bond setting, "
			 "using restart file value");
      }
      force->newton_bond = newton_bond_file;
      if (force->newton_pair || force->newton_bond) force->newton = 1;
      else force->newton = 0;

      // set boundary settings from restart file
      // warn if different and input script settings are not default

    } else if (flag == XPERIODIC) {
      xperiodic = read_int();
    } else if (flag == YPERIODIC) {
      yperiodic = read_int();
    } else if (flag == ZPERIODIC) {
      zperiodic = read_int();
    } else if (flag == BOUNDARY_00) {
      boundary[0][0] = read_int();
    } else if (flag == BOUNDARY_01) {
      boundary[0][1] = read_int();
    } else if (flag == BOUNDARY_10) {
      boundary[1][0] = read_int();
    } else if (flag == BOUNDARY_11) {
      boundary[1][1] = read_int();
    } else if (flag == BOUNDARY_20) {
      boundary[2][0] = read_int();
    } else if (flag == BOUNDARY_21) {
      boundary[2][1] = read_int();

      if (domain->boundary[0][0] || domain->boundary[0][1] || 
	  domain->boundary[1][0] || domain->boundary[1][1] || 
	  domain->boundary[2][0] || domain->boundary[2][1]) {
	if (boundary[0][0] != domain->boundary[0][0] ||
	    boundary[0][1] != domain->boundary[0][1] ||
	    boundary[1][0] != domain->boundary[1][0] ||
	    boundary[1][1] != domain->boundary[1][1] ||
	    boundary[2][0] != domain->boundary[2][0] ||
	    boundary[2][1] != domain->boundary[2][1]) {
	  if (me == 0) 
	    error->warning("Restart file used different boundary settings, "
			   "using restart file values");
	}
      }

      int flag = 0;
      if ((xperiodic != domain->xperiodic || yperiodic != domain->yperiodic ||
	   zperiodic != domain->zperiodic)) flag = 1;

      domain->boundary[0][0] = boundary[0][0];
      domain->boundary[0][1] = boundary[0][1];
      domain->boundary[1][0] = boundary[1][0];
      domain->boundary[1][1] = boundary[1][1];
      domain->boundary[2][0] = boundary[2][0];
      domain->boundary[2][1] = boundary[2][1];

      domain->periodicity[0] = domain->xperiodic = xperiodic;
      domain->periodicity[1] = domain->yperiodic = yperiodic;
      domain->periodicity[2] = domain->zperiodic = zperiodic;
  
      domain->nonperiodic = 0;
      if (xperiodic == 0 || yperiodic == 0 || zperiodic == 0) {
	domain->nonperiodic = 1;
	if (boundary[0][0] >= 2 || boundary[0][1] >= 2 ||
	    boundary[1][0] >= 2 || boundary[1][1] >= 2 ||
	    boundary[2][0] >= 2 || boundary[2][1] >= 2)
	  domain->nonperiodic = 2;
      }

      // create new AtomVec class
      // if style = hybrid, read additional sub-class arguments

    } else if (flag == ATOM_STYLE) {
      char *style = read_char();

      int nwords = 0;
      char **words = NULL;

      if (strcmp(style,"hybrid") == 0) {
	nwords = read_int();
	words = new char*[nwords];
	for (int i = 0; i < nwords; i++) words[i] = read_char();
      }

      atom->create_avec(style,nwords,words);
      for (int i = 0; i < nwords; i++) delete [] words[i];
      delete [] words;
      delete [] style;

    } else if (flag == NATOMS) {
      atom->natoms = read_double();
    } else if (flag == NTYPES) {
      atom->ntypes = read_int();
    } else if (flag == NBONDS) {
      atom->nbonds = read_int();
    } else if (flag == NBONDTYPES) {
      atom->nbondtypes = read_int();
    } else if (flag == BOND_PER_ATOM) {
      atom->bond_per_atom = read_int();
    } else if (flag == NANGLES) {
      atom->nangles = read_int();
    } else if (flag == NANGLETYPES) {
      atom->nangletypes = read_int();
    } else if (flag == ANGLE_PER_ATOM) {
      atom->angle_per_atom = read_int();
    } else if (flag == NDIHEDRALS) {
      atom->ndihedrals = read_int();
    } else if (flag == NDIHEDRALTYPES) {
      atom->ndihedraltypes = read_int();
    } else if (flag == DIHEDRAL_PER_ATOM) {
      atom->dihedral_per_atom = read_int();
    } else if (flag == NIMPROPERS) {
      atom->nimpropers = read_int();
    } else if (flag == NIMPROPERTYPES) {
      atom->nimpropertypes = read_int();
    } else if (flag == IMPROPER_PER_ATOM) {
      atom->improper_per_atom = read_int();

    } else if (flag == BOXLO_0) {
      domain->boxlo[0] = read_double();
    } else if (flag == BOXHI_0) {
      domain->boxhi[0] = read_double();
    } else if (flag == BOXLO_1) {
      domain->boxlo[1] = read_double();
    } else if (flag == BOXHI_1) {
      domain->boxhi[1] = read_double();
    } else if (flag == BOXLO_2) {
      domain->boxlo[2] = read_double();
    } else if (flag == BOXHI_2) {
      domain->boxhi[2] = read_double();

    } else if (flag == SPECIAL_LJ_1) {
      force->special_lj[1] = read_double();
    } else if (flag == SPECIAL_LJ_2) {
      force->special_lj[2] = read_double();
    } else if (flag == SPECIAL_LJ_3) {
      force->special_lj[3] = read_double();
    } else if (flag == SPECIAL_COUL_1) {
      force->special_coul[1] = read_double();
    } else if (flag == SPECIAL_COUL_2) {
      force->special_coul[2] = read_double();
    } else if (flag == SPECIAL_COUL_3) {
      force->special_coul[3] = read_double();

    } else if (flag == XY) {
      domain->triclinic = 1;
      domain->xy = read_double();
    } else if (flag == XZ) {
      domain->triclinic = 1;
      domain->xz = read_double();
    } else if (flag == YZ) {
      domain->triclinic = 1;
      domain->yz = read_double();

    } else error->all("Invalid flag in header section of restart file");

    flag = read_int();
  }
}

/* ---------------------------------------------------------------------- */

void ReadRestart::type_arrays()
{
  int flag = read_int();
  while (flag >= 0) {

    if (flag == MASS) {
      double *mass = new double[atom->ntypes+1];
      if (me == 0) fread(&mass[1],sizeof(double),atom->ntypes,fp);
      MPI_Bcast(&mass[1],atom->ntypes,MPI_DOUBLE,0,world);
      atom->set_mass(mass);
      delete [] mass;

    } else if (flag == SHAPE) {
      double **shape =
	memory->create_2d_double_array(atom->ntypes+1,3,"restart:shape");
      if (me == 0) fread(&shape[1][0],sizeof(double),atom->ntypes*3,fp);
      MPI_Bcast(&shape[1][0],atom->ntypes*3,MPI_DOUBLE,0,world);
      atom->set_shape(shape);
      memory->destroy_2d_double_array(shape);

    } else if (flag == DIPOLE) {
      double *dipole = new double[atom->ntypes+1];
      if (me == 0) fread(&dipole[1],sizeof(double),atom->ntypes,fp);
      MPI_Bcast(&dipole[1],atom->ntypes,MPI_DOUBLE,0,world);
      atom->set_dipole(dipole);
      delete [] dipole;

    } else error->all("Invalid flag in type arrays section of restart file");

    flag = read_int();
  }
}

/* ---------------------------------------------------------------------- */

void ReadRestart::force_fields()
{
  int n;
  char *style;

  int flag = read_int();
  while (flag >= 0) {

    if (flag == PAIR) {
      if (me == 0) fread(&n,sizeof(int),1,fp);
      MPI_Bcast(&n,1,MPI_INT,0,world);
      style = new char[n];
      if (me == 0) fread(style,sizeof(char),n,fp);
      MPI_Bcast(style,n,MPI_CHAR,0,world);

      force->create_pair(style);
      delete [] style;
      force->pair->read_restart(fp);

    } else if (flag == BOND) {
      if (me == 0) fread(&n,sizeof(int),1,fp);
      MPI_Bcast(&n,1,MPI_INT,0,world);
      style = new char[n];
      if (me == 0) fread(style,sizeof(char),n,fp);
      MPI_Bcast(style,n,MPI_CHAR,0,world);
      
      force->create_bond(style);
      delete [] style;
      force->bond->read_restart(fp);

    } else if (flag == ANGLE) {
      if (me == 0) fread(&n,sizeof(int),1,fp);
      MPI_Bcast(&n,1,MPI_INT,0,world);
      style = new char[n];
      if (me == 0) fread(style,sizeof(char),n,fp);
      MPI_Bcast(style,n,MPI_CHAR,0,world);

      force->create_angle(style);
      delete [] style;
      force->angle->read_restart(fp);

    } else if (flag == DIHEDRAL) {
      if (me == 0) fread(&n,sizeof(int),1,fp);
      MPI_Bcast(&n,1,MPI_INT,0,world);
      style = new char[n];
      if (me == 0) fread(style,sizeof(char),n,fp);
      MPI_Bcast(style,n,MPI_CHAR,0,world);

      force->create_dihedral(style);
      delete [] style;
      force->dihedral->read_restart(fp);

    } else if (flag == IMPROPER) {
      if (me == 0) fread(&n,sizeof(int),1,fp);
      MPI_Bcast(&n,1,MPI_INT,0,world);
      style = new char[n];
      if (me == 0) fread(style,sizeof(char),n,fp);
      MPI_Bcast(style,n,MPI_CHAR,0,world);

      force->create_improper(style);
      delete [] style;
      force->improper->read_restart(fp);

    } else error->all("Invalid flag in force field section of restart file");

    flag = read_int();
  }
}

/* ----------------------------------------------------------------------
   read an int from restart file and bcast it
------------------------------------------------------------------------- */

int ReadRestart::read_int()
{
  int value;
  if (me == 0) fread(&value,sizeof(int),1,fp);
  MPI_Bcast(&value,1,MPI_INT,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read a double from restart file and bcast it
------------------------------------------------------------------------- */

double ReadRestart::read_double()
{
  double value;
  if (me == 0) fread(&value,sizeof(double),1,fp);
  MPI_Bcast(&value,1,MPI_DOUBLE,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read a char str from restart file and bcast it
   str is allocated here, ptr is returned, caller must deallocate
------------------------------------------------------------------------- */

char *ReadRestart::read_char()
{
  int n;
  if (me == 0) fread(&n,sizeof(int),1,fp);
  MPI_Bcast(&n,1,MPI_INT,0,world);
  char *value = new char[n];
  if (me == 0) fread(value,sizeof(char),n,fp);
  MPI_Bcast(value,n,MPI_CHAR,0,world);
  return value;
}
