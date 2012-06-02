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

/* ----------------------------------------------------------------------
   Contributing author: Timothy Sirk (U Vermont)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "read_dump.h"
#include "read_dump_native.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "irregular.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

#define CHUNK 1024
#define EPSILON 1.0e-6

enum{ID,TYPE,X,Y,Z,VX,VY,VZ,IX,IY,IZ};
enum{UNSET,UNSCALED,SCALED};
enum{NATIVE};

/* ---------------------------------------------------------------------- */

ReadDump::ReadDump(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
}

/* ---------------------------------------------------------------------- */

void ReadDump::command(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal read_dump command");

  nstep = ATOBIGINT(arg[1]);

  // per-field vectors

  int firstfield = 2;
  fieldtype = new int[narg];
  fieldlabel = new char*[narg];

  // scan ahead to see if "add yes" keyword/value is used
  // requires extra "type" field from from dump file
  // add id and type fields as needed

  int iarg;
  for (iarg = firstfield; iarg < narg; iarg++)
    if (strcmp(arg[iarg],"add") == 0)
      if (iarg < narg-1 && strcmp(arg[iarg+1],"yes") == 0) break;

  nfield = 0;
  fieldtype[nfield++] = ID;
  if (iarg < narg) fieldtype[nfield++] = TYPE;

  // parse fields

  iarg = firstfield;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"x") == 0) fieldtype[nfield++] = X;
    else if (strcmp(arg[iarg],"y") == 0) fieldtype[nfield++] = Y;
    else if (strcmp(arg[iarg],"z") == 0) fieldtype[nfield++] = Z;
    else if (strcmp(arg[iarg],"vx") == 0) fieldtype[nfield++] = VX;
    else if (strcmp(arg[iarg],"vy") == 0) fieldtype[nfield++] = VY;
    else if (strcmp(arg[iarg],"vz") == 0) fieldtype[nfield++] = VZ;
    else if (strcmp(arg[iarg],"ix") == 0) fieldtype[nfield++] = IX;
    else if (strcmp(arg[iarg],"iy") == 0) fieldtype[nfield++] = IY;
    else if (strcmp(arg[iarg],"iz") == 0) fieldtype[nfield++] = IZ;
    else break;
    iarg++;
  }

  dimension = domain->dimension;
  triclinic = domain->triclinic;

  if (fieldtype[nfield-1] == ID || fieldtype[nfield-1] == TYPE)
    error->all(FLERR,"Illegal read_dump command");

  if (dimension == 2) {
    for (int i = 0; i < nfield; i++)
      if (fieldtype[i] == Z || fieldtype[i] == VZ || fieldtype[i] == IZ)
	error->all(FLERR,"Illegal read_dump command");
  }

  for (int i = 0; i < nfield; i++)
    for (int j = i+1; j < nfield; j++)
      if (fieldtype[i] == fieldtype[j])
	error->all(FLERR,"Duplicate fields in read_dump command");

  // parse optional args

  boxflag = 1;
  replaceflag = 1;
  purgeflag = 0;
  trimflag = 0;
  addflag = 0;
  for (int i = 0; i < nfield; i++) fieldlabel[i] = NULL;
  scaledflag = UNSCALED;
  format = NATIVE;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"box") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_dump command");
      if (strcmp(arg[iarg+1],"yes") == 0) boxflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) boxflag = 0;
      else error->all(FLERR,"Illegal read_dump command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"replace") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_dump command");
      if (strcmp(arg[iarg+1],"yes") == 0) replaceflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) replaceflag = 0;
      else error->all(FLERR,"Illegal read_dump command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"purge") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_dump command");
      if (strcmp(arg[iarg+1],"yes") == 0) purgeflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) purgeflag = 0;
      else error->all(FLERR,"Illegal read_dump command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"trim") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_dump command");
      if (strcmp(arg[iarg+1],"yes") == 0) trimflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) trimflag = 0;
      else error->all(FLERR,"Illegal read_dump command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"add") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_dump command");
      if (strcmp(arg[iarg+1],"yes") == 0) addflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) addflag = 0;
      else error->all(FLERR,"Illegal read_dump command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"label") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal read_dump command");
      int i;
      for (i = 0; i < nfield; i++)
	if (strcmp(arg[firstfield+i],arg[iarg+1]) == 0) break;
      if (i == nfield) error->all(FLERR,"Illegal read_dump command");
      fieldlabel[i] = arg[iarg+2];
      iarg += 3;
    } else if (strcmp(arg[iarg],"scaled") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_dump command");
      if (strcmp(arg[iarg+1],"yes") == 0) scaledflag = SCALED;
      else if (strcmp(arg[iarg+1],"no") == 0) scaledflag = UNSCALED;
      else error->all(FLERR,"Illegal read_dump command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_dump command");
      if (strcmp(arg[iarg+1],"native") == 0) format = NATIVE;
      else error->all(FLERR,"Illegal read_dump command");
      iarg += 2;
    } else error->all(FLERR,"Illegal read_dump command");
  }

  if (purgeflag && (replaceflag || trimflag))
    error->all(FLERR,"If read_dump purges it cannot replace or trim");

  // allocate snapshot field buffer

  memory->create(fields,CHUNK,nfield,"read_dump:fields");

  // create reader class
  // could make this a parent class and customize with other readers

  if (format == NATIVE) reader = new ReadDumpNative(lmp);

  // proc 0 opens dump file and scans to correct snapshot
  // after scan these values are set, so Bcast them:
  // nsnapatoms, box[3][3], scaled
  // NOTE: fieldlabel is just ptrs to input args in read_dump command
  //       will not persist if want to use labels in rerun() command

  if (me == 0) {
    if (screen) fprintf(screen,"Scanning dump file ...\n");
    open(arg[0]);
    reader->init(fp);
    reader->scan(nstep,nfield,fieldtype,fieldlabel,scaledflag,
		 nsnapatoms,box,scaled);
  }

  MPI_Bcast(&nsnapatoms,1,MPI_LMP_BIGINT,0,world);
  MPI_Bcast(&box[0][0],9,MPI_DOUBLE,0,world);
  MPI_Bcast(&scaled,1,MPI_INT,0,world);

  // for scaled coords and triclinic box:
  // yindex,zindex = index of Y and Z fields
  // already known to exist because checked in scan()
  // needed for unscaling to absolute coords in xfield(), yfield(), zfield()

  if (scaled == SCALED && triclinic) {
    for (int i = 0; i < nfield; i++) {
      if (fieldtype[i] == Y) yindex = i;
      if (fieldtype[i] == Z) zindex = i;
    }
  }

  // make local copy of snapshot box params

  xlo = box[0][0];
  xhi = box[0][1];
  ylo = box[1][0];
  yhi = box[1][1];
  zlo = box[2][0];
  zhi = box[2][1];
  xprd = xhi - xlo;
  yprd = yhi - ylo;
  zprd = zhi - zlo;
  if (triclinic) {
    xy = box[0][2];
    xz = box[1][2];
    yz = box[2][2];
  }

  // reset timestep to nstep

  char *tstr[1];
  char str[32];
  sprintf(str,BIGINT_FORMAT,nstep);
  tstr[0] = str;
  update->reset_timestep(1,tstr);

  // reset simulation box from snapshot box parameters if requested
  // do it now, so if adding atoms, procs will have correct sub-domains
  // call domain->reset_box() later,
  //   since can't shrink wrap until atom coords change and atoms are added

  if (boxflag) {
    domain->boxlo[0] = xlo;
    domain->boxhi[0] = xhi;
    domain->boxlo[1] = ylo;
    domain->boxhi[1] = yhi;
    if (dimension == 3) {
      domain->boxlo[2] = zlo;
      domain->boxhi[2] = zhi;
    }
    if (triclinic) {
      domain->xy = xy;
      if (dimension == 3) {
	domain->xz = xz;
	domain->yz = yz;
      }
    }

    domain->set_initial_box();
    domain->set_global_box();
    comm->set_proc_grid();
    domain->set_local_box();
  }

  // read in the snapshot

  if (me == 0)
    if (screen) fprintf(screen,"Reading snapshot from dump file ...\n");

  // counters

  bigint natoms_prev = atom->natoms;
  npurge = nreplace = ntrim = nadd = 0;

  // if purgeflag set, delete all current atoms

  if (purgeflag) {
    if (atom->map_style) atom->map_clear();
    npurge = atom->nlocal;
    atom->nlocal = atom->nghost = 0;
    atom->natoms = 0;
  }

  // to match existing atoms to dump atoms:
  // must build map if not a molecular system

  int mapflag = 0;
  if (atom->map_style == 0) {
    mapflag = 1;
    atom->map_style = 1;
    atom->map_init();
    atom->map_set();
  }

  // uflag[i] = 1 for each owned atom appearing in dump
  // ucflag = similar flag for each chunk atom, used in process_atoms()

  int nlocal = atom->nlocal;
  memory->create(uflag,nlocal,"read_dump:uflag");
  for (int i = 0; i < nlocal; i++) uflag[i] = 0;
  memory->create(ucflag,CHUNK,"read_dump:ucflag");
  memory->create(ucflag_all,CHUNK,"read_dump:ucflag");
  
  // read, broadcast, and process atoms from snapshot in chunks

  addproc = -1;

  int nchunk;
  bigint nread = 0;
  while (nread < nsnapatoms) {
    nchunk = MIN(nsnapatoms-nread,CHUNK);
    if (me == 0) reader->read(nchunk,fields);
    MPI_Bcast(&fields[0][0],nchunk*nfield,MPI_DOUBLE,0,world);
    process_atoms(nchunk);
    nread += nchunk;
  }

  // if addflag set, add tags to new atoms if possible

  if (addflag) {
    bigint nblocal = atom->nlocal;
    MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
    if (atom->natoms < 0 || atom->natoms > MAXBIGINT)
      error->all(FLERR,"Too many total atoms");
    if (atom->natoms > MAXTAGINT) atom->tag_enable = 0;
    if (atom->natoms <= MAXTAGINT) atom->tag_extend();
  }

  // if trimflag set, delete atoms not replaced by snapshot atoms

  if (trimflag) {
    delete_atoms(uflag);
    bigint nblocal = atom->nlocal;
    MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }

  // delete atom map if created it above
  // else reinitialize map for current atoms
  // do this before migrating atoms to new procs via Irregular

  if (mapflag) {
    atom->map_delete();
    atom->map_style = 0;
  } else {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // close dump file

  if (me == 0) {
    if (compressed) pclose(fp);
    else fclose(fp);
  }

  // move atoms back inside simulation box and to new processors
  // use remap() instead of pbc() in case atoms moved a long distance
  // adjust image flags of all atoms (old and new) based on current box
  // use irregular() in case atoms moved a long distance

  double **x = atom->x;
  int *image = atom->image;
  nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) domain->remap(x[i],image[i]);

  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->reset_box();
  Irregular *irregular = new Irregular(lmp);
  irregular->migrate_atoms();
  delete irregular;
  if (triclinic) domain->lamda2x(atom->nlocal);

  domain->print_box("  ");

  // clean up

  delete reader;
  delete [] fieldtype;
  delete [] fieldlabel;
  memory->destroy(fields);
  memory->destroy(uflag);
  memory->destroy(ucflag);
  memory->destroy(ucflag_all);

  // print out stats

  bigint npurge_all,nreplace_all,ntrim_all,nadd_all;

  bigint tmp;
  tmp = npurge;
  MPI_Allreduce(&tmp,&npurge_all,1,MPI_LMP_BIGINT,MPI_SUM,world);
  tmp = nreplace;
  MPI_Allreduce(&tmp,&nreplace_all,1,MPI_LMP_BIGINT,MPI_SUM,world);
  tmp = ntrim;
  MPI_Allreduce(&tmp,&ntrim_all,1,MPI_LMP_BIGINT,MPI_SUM,world);
  tmp = nadd;
  MPI_Allreduce(&tmp,&nadd_all,1,MPI_LMP_BIGINT,MPI_SUM,world);

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  " BIGINT_FORMAT " atoms before read\n",natoms_prev);
      fprintf(screen,"  " BIGINT_FORMAT " atoms in snapshot\n",nsnapatoms);
      fprintf(screen,"  " BIGINT_FORMAT " atoms purged\n",npurge_all);
      fprintf(screen,"  " BIGINT_FORMAT " atoms replaced\n",nreplace_all);
      fprintf(screen,"  " BIGINT_FORMAT " atoms trimmed\n",ntrim_all);
      fprintf(screen,"  " BIGINT_FORMAT " atoms added\n",nadd_all);
      fprintf(screen,"  " BIGINT_FORMAT " atoms after read\n",atom->natoms);
    }
    if (logfile) {
      fprintf(logfile,"  " BIGINT_FORMAT " atoms before read\n",natoms_prev);
      fprintf(logfile,"  " BIGINT_FORMAT " atoms in snapshot\n",nsnapatoms);
      fprintf(logfile,"  " BIGINT_FORMAT " atoms purged\n",npurge_all);
      fprintf(logfile,"  " BIGINT_FORMAT " atoms replaced\n",nreplace_all);
      fprintf(logfile,"  " BIGINT_FORMAT " atoms trimmed\n",ntrim_all);
      fprintf(logfile,"  " BIGINT_FORMAT " atoms added\n",nadd_all);
      fprintf(logfile,"  " BIGINT_FORMAT " atoms after read\n",atom->natoms);
    }
  }
}


/* ----------------------------------------------------------------------
   process each of N atoms in chunk read from dump file
   if in replace mode and atom ID matches current atom,
     overwrite atom info with fields from dump file
   if in add mode and atom ID does not match any current atom,
     create new atom with dump file field values,
     and assign to a proc in round-robin manner
   use round-robin method, b/c atom coords may not be inside simulation box
------------------------------------------------------------------------- */

void ReadDump::process_atoms(int n)
{
  int i,m,ifield,itype;
  int xbox,ybox,zbox;

  double **x = atom->x;
  double **v = atom->v;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  for (i = 0; i < n; i++) {
    ucflag[i] = 0;

    // map() call is invalid if purged all atoms
    // setting m = -1 forces new atom not to match

    if (!purgeflag) m = atom->map(static_cast<int> (fields[i][0]));
    else m = -1;
    if (m < 0 || m >= nlocal) continue;

    ucflag[i] = 1;
    uflag[m] = 1;

    if (replaceflag) {
      nreplace++;

      // current image flags
      
      xbox = (image[m] & 1023) - 512;
      ybox = (image[m] >> 10 & 1023) - 512;
      zbox = (image[m] >> 20) - 512;
    
      // overwrite atom attributes with field info
      // start from field 1 since 0 = id, 1 will be skipped if type
      
      for (ifield = 1; ifield < nfield; ifield++) {
	switch (fieldtype[ifield]) {
	case X:
	  x[m][0] = xfield(i,ifield);
	  break;
	case Y:
	  x[m][1] = yfield(i,ifield);
	  break;
	case Z:
	  x[m][2] = zfield(i,ifield);
	  break;
	case VX:
	  v[m][0] = fields[i][ifield];
	  break;
	case VY:
	  v[m][1] = fields[i][ifield];
	  break;
	case VZ:
	  v[m][2] = fields[i][ifield];
	  break;
	case IX:
	  xbox = static_cast<int> (fields[i][ifield]);
	  break;
	case IY:
	  ybox = static_cast<int> (fields[i][ifield]);
	  break;
	case IZ:
	  zbox = static_cast<int> (fields[i][ifield]);
	  break;
	}
      }
      
      // replace image flag in case changed by ix,iy,iz fields
      
      image[m] = (xbox << 20) | (ybox << 10) | zbox;
    }
  }

  // create any atoms in chunk that no processor owned
  // add atoms in round-robin sequence on processors
  // cannot do it geometrically b/c dump coords may not be in simulation box

  if (!addflag) return;

  MPI_Allreduce(ucflag,ucflag_all,n,MPI_INT,MPI_SUM,world);

  double lamda[3],one[3];
  double *coord;

  for (i = 0; i < n; i++) {
    if (ucflag_all[i]) continue;

    // each processor adds every Pth atom

    addproc++;
    if (addproc == nprocs) addproc = 0;
    if (addproc != me) continue;

    // create type and coord fields from dump file
    // coord = 0.0 unless corresponding dump file field was specified

    one[0] = one[1] = one[2] = 0.0;
    for (ifield = 1; ifield < nfield; ifield++) {
      switch (fieldtype[ifield]) {
      case TYPE:
	itype = static_cast<int> (fields[i][ifield]);
	break;
      case X:
	one[0] = xfield(i,ifield);
	break;
      case Y:
	one[1] = yfield(i,ifield);
	break;
      case Z:
	one[2] = zfield(i,ifield);
	break;
      }
    }

    // create the atom on proc that owns it
    // reset v,image ptrs in case they are reallocated

    atom->avec->create_atom(itype,one);
    nadd++;
    
    v = atom->v;
    image = atom->image;
    m = atom->nlocal;

    // set atom attributes from other dump file fields
    // xyzbox = 512 is default value set by create_atom()

    xbox = ybox = zbox = 512;

    for (ifield = 1; ifield < nfield; ifield++) {
      switch (fieldtype[ifield]) {
      case VX:
	v[m][0] = fields[i][ifield];
	break;
      case VY:
	v[m][1] = fields[i][ifield];
	break;
      case VZ:
	v[m][2] = fields[i][ifield];
	break;
      case IX:
	xbox = static_cast<int> (fields[i][ifield]);
	break;
      case IY:
	ybox = static_cast<int> (fields[i][ifield]);
	break;
      case IZ:
	zbox = static_cast<int> (fields[i][ifield]);
	break;
      }

      // replace image flag in case changed by ix,iy,iz fields
      
      image[m] = (xbox << 20) | (ybox << 10) | zbox;
    }
  }
}

/* ----------------------------------------------------------------------
   delete atoms not flagged as replaced by dump atoms
------------------------------------------------------------------------- */

void ReadDump::delete_atoms(int *uflag)
{
  AtomVec *avec = atom->avec;
  int nlocal = atom->nlocal;

  int i = 0;
  while (i < nlocal) {
    if (uflag[i] == 0) {
      avec->copy(nlocal-1,i,1);
      uflag[i] = uflag[nlocal-1];
      nlocal--;
      ntrim++;
    } else i++;
  }

  atom->nlocal = nlocal;
}

/* ----------------------------------------------------------------------
   convert XYZ fields in dump file into absolute, unscaled coordinates
   depends on scaled vs unscaled and triclinic vs orthogonal
   does not depend on wrapped vs unwrapped
------------------------------------------------------------------------- */

double ReadDump::xfield(int i, int j)
{
  if (scaled == UNSCALED) return fields[i][j];
  else if (!triclinic) return fields[i][j]*xprd + xlo;
  else if (dimension == 2) 
    return xprd*fields[i][j] + xy*fields[i][yindex] + xlo;
  return xprd*fields[i][j] + xy*fields[i][yindex] + xz*fields[i][zindex] + xlo;
}

double ReadDump::yfield(int i, int j)
{
  if (scaled == UNSCALED) return fields[i][j];
  else if (!triclinic) return fields[i][j]*yprd + ylo;
  else if (dimension == 2) return yprd*fields[i][j] + ylo;
  return yprd*fields[i][j] + yz*fields[i][zindex] + ylo;
}

double ReadDump::zfield(int i, int j)
{
  if (scaled == UNSCALED) return fields[i][j];
  return fields[i][j]*zprd + zlo;
}

/* ----------------------------------------------------------------------
   proc 0 opens dump file
   test if gzipped
------------------------------------------------------------------------- */

void ReadDump::open(char *file)
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
