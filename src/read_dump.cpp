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
   Contributing author: Timothy Sirk (ARL)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "read_dump.h"
#include "reader.h"
#include "style_reader.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "domain.h"
#include "comm.h"
#include "irregular.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

#define CHUNK 1024
#define EPSILON 1.0e-6

// also in reader_native.cpp

enum{ID,TYPE,X,Y,Z,VX,VY,VZ,Q,IX,IY,IZ};
enum{UNSET,NOSCALE_NOWRAP,NOSCALE_WRAP,SCALE_NOWRAP,SCALE_WRAP};

/* ---------------------------------------------------------------------- */

ReadDump::ReadDump(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  dimension = domain->dimension;
  triclinic = domain->triclinic;

  nfile = 0;
  files = NULL;

  nfield = 0;
  fieldtype = NULL;
  fieldlabel = NULL;
  fields = NULL;

  int n = strlen("native") + 1;
  readerstyle = new char[n];
  strcpy(readerstyle,"native");

  reader = NULL;
  fp = NULL;
}

/* ---------------------------------------------------------------------- */

ReadDump::~ReadDump()
{
  for (int i = 0; i < nfile; i++) delete [] files[i];
  delete [] files;
  for (int i = 0; i < nfield; i++) delete [] fieldlabel[i];
  delete [] fieldlabel;
  delete [] fieldtype;
  delete [] readerstyle;

  memory->destroy(fields);
  delete reader;
}

/* ---------------------------------------------------------------------- */

void ReadDump::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Read_dump command before simulation box is defined");

  if (narg < 2) error->all(FLERR,"Illegal read_dump command");

  store_files(1,&arg[0]);
  bigint nstep = ATOBIGINT(arg[1]);

  int nremain = narg - 2;
  if (nremain) nremain = fields_and_keywords(nremain,&arg[narg-nremain]);
  else nremain = fields_and_keywords(0,NULL);
  if (nremain) setup_reader(nremain,&arg[narg-nremain]);
  else setup_reader(0,NULL);

  // find the snapshot and read/bcast/process header info

  if (me == 0 && screen) fprintf(screen,"Scanning dump file ...\n");

  bigint ntimestep = seek(nstep,1);
  if (ntimestep < 0)
    error->all(FLERR,"Dump file does not contain requested snapshot");
  header(1);

  // reset timestep to nstep

  update->reset_timestep(nstep);

  // counters

  // read in the snapshot and reset system

  if (me == 0 && screen)
    fprintf(screen,"Reading snapshot from dump file ...\n");

  bigint natoms_prev = atom->natoms;
  atoms();

  if (me == 0) reader->close_file();

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

  domain->print_box("  ");

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

/* ---------------------------------------------------------------------- */

void ReadDump::store_files(int nstr, char **str)
{
  nfile = nstr;
  files = new char*[nfile];

  for (int i = 0; i < nfile; i++) {
    int n = strlen(str[i]) + 1;
    files[i] = new char[n];
    strcpy(files[i],str[i]);
  }
}

/* ---------------------------------------------------------------------- */

void ReadDump::setup_reader(int narg, char **arg)
{
  // allocate snapshot field buffer

  memory->create(fields,CHUNK,nfield,"read_dump:fields");

  // create reader class
  // match readerstyle to options in style_reader.h

  if (0) return;        // dummy line to enable else-if macro expansion

#define READER_CLASS
#define ReaderStyle(key,Class) \
  else if (strcmp(readerstyle,#key) == 0) reader = new Class(lmp);
#include "style_reader.h"
#undef READER_CLASS

  // unrecognized style

  else error->all(FLERR,"Invalid dump reader style");

  // pass any arguments to reader

  if (narg > 0) reader->settings(narg,arg);
}

/* ----------------------------------------------------------------------
   seek Nrequest timestep in one or more dump files
   if exact = 1, must find exactly Nrequest
   if exact = 0, find first step >= Nrequest
   return matching ntimestep or -1 if did not find a match
------------------------------------------------------------------------- */

bigint ReadDump::seek(bigint nrequest, int exact)
{
  int ifile,eofflag;
  bigint ntimestep;

  if (me == 0) {

    // exit file loop when dump timestep >= nrequest
    // or files exhausted

    for (ifile = 0; ifile < nfile; ifile++) {
      ntimestep = -1;
      reader->open_file(files[ifile]);
      while (1) {
        eofflag = reader->read_time(ntimestep);
        if (eofflag) break;
        if (ntimestep >= nrequest) break;
        reader->skip();
      }
      if (ntimestep >= nrequest) break;
      reader->close_file();
    }

    currentfile = ifile;
    if (ntimestep < nrequest) ntimestep = -1;
    if (exact && ntimestep != nrequest) ntimestep = -1;
    if (ntimestep < 0) reader->close_file();
  }

  MPI_Bcast(&ntimestep,1,MPI_LMP_BIGINT,0,world);
  return ntimestep;
}

/* ----------------------------------------------------------------------
   find next matching snapshot in one or more dump files
   Ncurrent = current timestep from last snapshot
   Nlast = match no timestep bigger than Nlast
   Nevery = only match timesteps that are a multiple of Nevery
   Nskip = skip every this many timesteps
   return matching ntimestep or -1 if did not find a match
------------------------------------------------------------------------- */

bigint ReadDump::next(bigint ncurrent, bigint nlast, int nevery, int nskip)
{
  int ifile,eofflag;
  bigint ntimestep;

  if (me == 0) {

    // exit file loop when dump timestep matches all criteria
    // or files exhausted

    int iskip = 0;

    for (ifile = currentfile; ifile < nfile; ifile++) {
      ntimestep = -1;
      if (ifile != currentfile) reader->open_file(files[ifile]);
      while (1) {
        eofflag = reader->read_time(ntimestep);
        if (iskip == nskip) iskip = 0;
        iskip++;
        if (eofflag) break;
        if (ntimestep <= ncurrent) break;
        if (ntimestep > nlast) break;
        if (nevery && ntimestep % nevery) reader->skip();
        else if (iskip < nskip) reader->skip();
        else break;
      }
      if (eofflag) reader->close_file();
      else break;
    }

    currentfile = ifile;
    if (eofflag) ntimestep = -1;
    if (ntimestep <= ncurrent) ntimestep = -1;
    if (ntimestep > nlast) ntimestep = -1;
    if (ntimestep < 0) reader->close_file();
  }

  MPI_Bcast(&ntimestep,1,MPI_LMP_BIGINT,0,world);
  return ntimestep;
}

/* ----------------------------------------------------------------------
   read and broadcast and store snapshot header info
   set nsnapatoms = # of atoms in snapshot
------------------------------------------------------------------------- */

void ReadDump::header(int fieldinfo)
{
  int triclinic_snap;
  int fieldflag,xflag,yflag,zflag;

  if (me == 0)
    nsnapatoms = reader->read_header(box,triclinic_snap,
                                     fieldinfo,nfield,fieldtype,fieldlabel,
                                     scaleflag,wrapflag,fieldflag,
                                     xflag,yflag,zflag);

  MPI_Bcast(&nsnapatoms,1,MPI_LMP_BIGINT,0,world);
  MPI_Bcast(&triclinic_snap,1,MPI_INT,0,world);
  MPI_Bcast(&box[0][0],9,MPI_DOUBLE,0,world);

  // local copy of snapshot box parameters
  // used in xfield,yfield,zfield when converting dump atom to absolute coords

  xlo = box[0][0];
  xhi = box[0][1];
  ylo = box[1][0];
  yhi = box[1][1];
  zlo = box[2][0];
  zhi = box[2][1];
  if (triclinic_snap) {
    xy = box[0][2];
    xz = box[1][2];
    yz = box[2][2];
    double xdelta = MIN(0.0,xy);
    xdelta = MIN(xdelta,xz);
    xdelta = MIN(xdelta,xy+xz);
    xlo = xlo - xdelta;
    xdelta = MAX(0.0,xy);
    xdelta = MAX(xdelta,xz);
    xdelta = MAX(xdelta,xy+xz);
    xhi = xhi - xdelta;
    ylo = ylo - MIN(0.0,yz);
    yhi = yhi - MAX(0.0,yz);
  }
  xprd = xhi - xlo;
  yprd = yhi - ylo;
  zprd = zhi - zlo;

  // done if not checking fields

  if (!fieldinfo) return;

  MPI_Bcast(&fieldflag,1,MPI_INT,0,world);
  MPI_Bcast(&xflag,1,MPI_INT,0,world);
  MPI_Bcast(&yflag,1,MPI_INT,0,world);
  MPI_Bcast(&zflag,1,MPI_INT,0,world);

  // error check on current vs new box and fields
  // triclinic_snap < 0 means no box info in file

  if (triclinic_snap < 0 && boxflag > 0)
    error->all(FLERR,"No box information in dump. You have to use 'box no'");
  if (triclinic_snap >= 0) {
    if ((triclinic_snap && !triclinic) ||
        (!triclinic_snap && triclinic))
      error->one(FLERR,"Read_dump triclinic status does not match simulation");
  }

  // error check on requested fields exisiting in dump file

  if (fieldflag < 0)
    error->one(FLERR,"Read_dump field not found in dump file");

  // all explicitly requested x,y,z must have consistent scaling & wrapping

  int value = MAX(xflag,yflag);
  value = MAX(zflag,value);
  if ((xflag != UNSET && xflag != value) ||
      (yflag != UNSET && yflag != value) ||
      (zflag != UNSET && zflag != value))
    error->one(FLERR,
               "Read_dump xyz fields do not have consistent scaling/wrapping");

  // set scaled/wrapped based on xyz flags

  value = UNSET;
  if (xflag != UNSET) value = xflag;
  if (yflag != UNSET) value = yflag;
  if (zflag != UNSET) value = zflag;

  if (value == UNSET) {
    scaled = wrapped = 0;
  } else if (value == NOSCALE_NOWRAP) {
    scaled = wrapped = 0;
  } else if (value == NOSCALE_WRAP) {
    scaled = 0;
    wrapped = 1;
  } else if (value == SCALE_NOWRAP) {
    scaled = 1;
    wrapped = 0;
  } else if (value == SCALE_WRAP) {
    scaled = wrapped = 1;
  }

  // scaled, triclinic coords require all 3 x,y,z fields, to perform unscaling
  // set yindex,zindex = column index of Y and Z fields in fields array
  // needed for unscaling to absolute coords in xfield(), yfield(), zfield()

  if (scaled && triclinic == 1) {
    int flag = 0;
    if (xflag == UNSET) flag = 1;
    if (yflag == UNSET) flag = 1;
    if (dimension == 3 && zflag == UNSET) flag = 1;
    if (flag)
      error->one(FLERR,"All read_dump x,y,z fields must be specified for "
                 "scaled, triclinic coords");

    for (int i = 0; i < nfield; i++) {
      if (fieldtype[i] == Y) yindex = i;
      if (fieldtype[i] == Z) zindex = i;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ReadDump::atoms()
{
  // initialize counters

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
    if (me == 0) reader->read_atoms(nchunk,nfield,fields);
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
    // change these to MAXTAGINT when allow tagint = bigint
    if (atom->natoms > MAXSMALLINT) atom->tag_enable = 0;
    if (atom->natoms <= MAXSMALLINT) atom->tag_extend();
  }

  // if trimflag set, delete atoms not replaced by snapshot atoms

  if (trimflag) {
    delete_atoms();
    bigint nblocal = atom->nlocal;
    MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }

  // can now delete uflag arrays

  memory->destroy(uflag);
  memory->destroy(ucflag);
  memory->destroy(ucflag_all);

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

  // overwrite simulation box with dump snapshot box if requested
  // reallocate processors to box

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
    comm->set_proc_grid(0);
    domain->set_local_box();
  }

  // move atoms back inside simulation box and to new processors
  // use remap() instead of pbc() in case atoms moved a long distance
  // adjust image flags of all atoms (old and new) based on current box
  // use irregular() in case atoms moved a long distance

  double **x = atom->x;
  imageint *image = atom->image;
  nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) domain->remap(x[i],image[i]);

  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->reset_box();
  Irregular *irregular = new Irregular(lmp);
  irregular->migrate_atoms();
  delete irregular;
  if (triclinic) domain->lamda2x(atom->nlocal);
}

/* ----------------------------------------------------------------------
   process arg list for dump file fields and optional keywords
------------------------------------------------------------------------- */

int ReadDump::fields_and_keywords(int narg, char **arg)
{
  // per-field vectors, leave space for ID and TYPE

  fieldtype = new int[narg+2];
  fieldlabel = new char*[narg+2];

  // add id and type fields as needed
  // scan ahead to see if "add yes" keyword/value is used
  // requires extra "type" field from from dump file

  int iarg;
  for (iarg = 0; iarg < narg; iarg++)
    if (strcmp(arg[iarg],"add") == 0)
      if (iarg < narg-1 && strcmp(arg[iarg+1],"yes") == 0) break;

  nfield = 0;
  fieldtype[nfield++] = ID;
  if (iarg < narg) fieldtype[nfield++] = TYPE;

  // parse fields

  iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"x") == 0) fieldtype[nfield++] = X;
    else if (strcmp(arg[iarg],"y") == 0) fieldtype[nfield++] = Y;
    else if (strcmp(arg[iarg],"z") == 0) fieldtype[nfield++] = Z;
    else if (strcmp(arg[iarg],"vx") == 0) fieldtype[nfield++] = VX;
    else if (strcmp(arg[iarg],"vy") == 0) fieldtype[nfield++] = VY;
    else if (strcmp(arg[iarg],"vz") == 0) fieldtype[nfield++] = VZ;
    else if (strcmp(arg[iarg],"q") == 0) {
      if (!atom->q_flag)
        error->all(FLERR,"Read dump of atom property that isn't allocated");
      fieldtype[nfield++] = Q;
    }
    else if (strcmp(arg[iarg],"ix") == 0) fieldtype[nfield++] = IX;
    else if (strcmp(arg[iarg],"iy") == 0) fieldtype[nfield++] = IY;
    else if (strcmp(arg[iarg],"iz") == 0) fieldtype[nfield++] = IZ;
    else break;
    iarg++;
  }

  // check for no fields

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
  scaleflag = 0;
  wrapflag = 1;

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
        if (fieldlabel[i] && strcmp(arg[iarg+1],fieldlabel[i]) == 0) break;
      if (i == nfield) error->all(FLERR,"Illegal read_dump command");
      int n = strlen(arg[iarg+2]) + 1;
      fieldlabel[i] = new char[n];
      strcpy(fieldlabel[i],arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"scaled") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_dump command");
      if (strcmp(arg[iarg+1],"yes") == 0) scaleflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) scaleflag = 0;
      else error->all(FLERR,"Illegal read_dump command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"wrapped") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_dump command");
      if (strcmp(arg[iarg+1],"yes") == 0) wrapflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) wrapflag = 0;
      else error->all(FLERR,"Illegal read_dump command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_dump command");
      delete [] readerstyle;
      int n = strlen(arg[iarg+1]) + 1;
      readerstyle = new char[n];
      strcpy(readerstyle,arg[iarg+1]);
      iarg += 2;
      break;
    } else error->all(FLERR,"Illegal read_dump command");
  }

  if (purgeflag && (replaceflag || trimflag))
    error->all(FLERR,"If read_dump purges it cannot replace or trim");

  return narg-iarg;
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
  int i,m,ifield,itype,itag;;
  int xbox,ybox,zbox;

  double **x = atom->x;
  double **v = atom->v;
  double *q = atom->q;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  int map_tag_max = atom->map_tag_max;

  for (i = 0; i < n; i++) {
    ucflag[i] = 0;

    // check if new atom matches one I own
    // setting m = -1 forces new atom not to match

    itag = static_cast<int> (fields[i][0]);
    if (itag <= map_tag_max) m = atom->map(static_cast<int> (fields[i][0]));
    else m = -1;
    if (m < 0 || m >= nlocal) continue;

    ucflag[i] = 1;
    uflag[m] = 1;

    if (replaceflag) {
      nreplace++;

      // current image flags

      xbox = (image[m] & IMGMASK) - IMGMAX;
      ybox = (image[m] >> IMGBITS & IMGMASK) - IMGMAX;
      zbox = (image[m] >> IMG2BITS) - IMGMAX;

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
        case Q:
          q[m] = fields[i][ifield];
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

      // replace image flag in case changed by ix,iy,iz fields or unwrapping

      if (!wrapped) xbox = ybox = zbox = 0;

      image[m] = ((imageint) (xbox + IMGMAX) & IMGMASK) | 
        (((imageint) (ybox + IMGMAX) & IMGMASK) << IMGBITS) | 
        (((imageint) (zbox + IMGMAX) & IMGMASK) << IMG2BITS);
    }
  }

  // create any atoms in chunk that no processor owned
  // add atoms in round-robin sequence on processors
  // cannot do it geometrically b/c dump coords may not be in simulation box

  if (!addflag) return;

  MPI_Allreduce(ucflag,ucflag_all,n,MPI_INT,MPI_SUM,world);

  int nlocal_previous = atom->nlocal;
  double one[3];

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

    m = atom->nlocal;
    atom->avec->create_atom(itype,one);
    nadd++;

    v = atom->v;
    q = atom->q;
    image = atom->image;

    // set atom attributes from other dump file fields

    xbox = ybox = zbox = 0;

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
      case Q:
        q[m] = fields[i][ifield];
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

      image[m] = ((imageint) (xbox + IMGMAX) & IMGMASK) | 
        (((imageint) (ybox + IMGMAX) & IMGMASK) << IMGBITS) | 
        (((imageint) (zbox + IMGMAX) & IMGMASK) << IMG2BITS);
    }
  }

  // invoke set_arrays() for fixes that need initialization of new atoms
  // same as in CreateAtoms

  nlocal = atom->nlocal;
  for (m = 0; m < modify->nfix; m++) {
    Fix *fix = modify->fix[m];
    if (fix->create_attribute)
      for (i = nlocal_previous; i < nlocal; i++)
        fix->set_arrays(i);
  }
}

/* ----------------------------------------------------------------------
   delete atoms not flagged as replaced by dump atoms
------------------------------------------------------------------------- */

void ReadDump::delete_atoms()
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
  if (!scaled) return fields[i][j];
  else if (!triclinic) return fields[i][j]*xprd + xlo;
  else if (dimension == 2)
    return xprd*fields[i][j] + xy*fields[i][yindex] + xlo;
  return xprd*fields[i][j] + xy*fields[i][yindex] + xz*fields[i][zindex] + xlo;
}

double ReadDump::yfield(int i, int j)
{
  if (!scaled) return fields[i][j];
  else if (!triclinic) return fields[i][j]*yprd + ylo;
  else if (dimension == 2) return yprd*fields[i][j] + ylo;
  return yprd*fields[i][j] + yz*fields[i][zindex] + ylo;
}

double ReadDump::zfield(int i, int j)
{
  if (!scaled) return fields[i][j];
  return fields[i][j]*zprd + zlo;
}
