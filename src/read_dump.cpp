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

/* ----------------------------------------------------------------------
   Contributing author: Timothy Sirk (ARL)
------------------------------------------------------------------------- */

#include "read_dump.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "irregular.h"
#include "memory.h"
#include "reader.h"
#include "style_reader.h"       // IWYU pragma: keep
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

#define CHUNK 16384

// also in reader_native.cpp

enum{ID,TYPE,X,Y,Z,VX,VY,VZ,Q,IX,IY,IZ,FX,FY,FZ};
enum{UNSET,NOSCALE_NOWRAP,NOSCALE_WRAP,SCALE_NOWRAP,SCALE_WRAP};
enum{NOADD,YESADD,KEEPADD};

/* ---------------------------------------------------------------------- */

ReadDump::ReadDump(LAMMPS *lmp) : Command(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  dimension = domain->dimension;
  triclinic = domain->triclinic;

  nfile = 0;
  files = nullptr;

  nnew = maxnew = 0;
  nfield = 0;
  fieldtype = nullptr;
  fieldlabel = nullptr;
  fields = nullptr;
  buf = nullptr;

  readerstyle = utils::strdup("native");

  nreader = 0;
  readers = nullptr;
  nsnapatoms = nullptr;
  clustercomm = MPI_COMM_NULL;
  filereader = 0;
  parallel = 0;
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
  memory->destroy(buf);

  for (int i = 0; i < nreader; i++) delete readers[i];
  delete [] readers;
  delete [] nsnapatoms;

  MPI_Comm_free(&clustercomm);
}

/* ---------------------------------------------------------------------- */

void ReadDump::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Read_dump command before simulation box is defined");

  if (narg < 2) error->all(FLERR,"Illegal read_dump command");

  store_files(1,&arg[0]);
  bigint nstep = utils::bnumeric(FLERR,arg[1],false,lmp);

  int nremain = narg - 2;
  if (nremain) nremain = fields_and_keywords(nremain,&arg[narg-nremain]);
  else nremain = fields_and_keywords(0,nullptr);
  if (nremain) setup_reader(nremain,&arg[narg-nremain]);
  else setup_reader(0,nullptr);

  // find the snapshot and read/bcast/process header info

  if (me == 0) utils::logmesg(lmp,"Scanning dump file ...\n");

  bigint ntimestep = seek(nstep,1);
  if (ntimestep < 0)
    error->all(FLERR,"Dump file does not contain requested snapshot");
  header(1);

  // reset timestep to nstep

  update->reset_timestep(nstep);

  // counters

  // read in the snapshot and reset system

  if (me == 0) utils::logmesg(lmp,"Reading snapshot from dump file ...\n");

  bigint natoms_prev = atom->natoms;
  atoms();

  if (filereader)
    for (int i = 0; i < nreader; i++)
      readers[i]->close_file();

  // print out stats

  bigint nsnap_all,npurge_all,nreplace_all,ntrim_all,nadd_all;

  bigint tmp = 0;
  if (filereader)
    for (int i = 0; i < nreader; i++)
      tmp += nsnapatoms[i];
  MPI_Allreduce(&tmp,&nsnap_all,1,MPI_LMP_BIGINT,MPI_SUM,world);

  tmp = npurge;
  MPI_Allreduce(&tmp,&npurge_all,1,MPI_LMP_BIGINT,MPI_SUM,world);
  tmp = nreplace;
  MPI_Allreduce(&tmp,&nreplace_all,1,MPI_LMP_BIGINT,MPI_SUM,world);
  tmp = ntrim;
  MPI_Allreduce(&tmp,&ntrim_all,1,MPI_LMP_BIGINT,MPI_SUM,world);
  tmp = nadd;
  MPI_Allreduce(&tmp,&nadd_all,1,MPI_LMP_BIGINT,MPI_SUM,world);

  domain->print_box("  ");

  if (me == 0)
    utils::logmesg(lmp,"  {} atoms before read\n"
                   "  {} atoms in snapshot\n"
                   "  {} atoms purged\n"
                   "  {} atoms replaced\n"
                   "  {} atoms trimmed\n"
                   "  {} atoms added\n"
                   "  {} atoms after read\n",natoms_prev,nsnap_all,
                   npurge_all,nreplace_all,ntrim_all,nadd_all,atom->natoms);
}

/* ---------------------------------------------------------------------- */

void ReadDump::store_files(int nstr, char **str)
{
  nfile = nstr;
  files = new char*[nfile];

  // either all or none of files must have '%' wild-card

  for (int i = 0; i < nfile; i++) {
    files[i] = utils::strdup(str[i]);

    if (i == 0) {
      if (strchr(files[i],'%')) multiproc = 1;
      else multiproc = 0;
    } else {
      if (multiproc && !strchr(files[i],'%'))
        error->all(FLERR,"All read_dump files must be serial or parallel");
      if (!multiproc && strchr(files[i],'%'))
        error->all(FLERR,"All read_dump files must be serial or parallel");
    }
  }
}

/* ---------------------------------------------------------------------- */

void ReadDump::setup_reader(int narg, char **arg)
{
  // setup serial or parallel file reading
  // multiproc = 0: only one file to read from, only proc 0 is a reader
  // multiproc_nfile >= nprocs: every proc reads one or more files
  // multiproc_nfile < nprocs: multiproc_nfile readers, create clusters
  // see read_dump.h for explanation of these variables

  if (multiproc == 0) {
    nreader = 1;
    firstfile = -1;
    MPI_Comm_dup(world,&clustercomm);
  } else if (multiproc_nfile >= nprocs) {
    firstfile = static_cast<int> ((bigint) me * multiproc_nfile/nprocs);
    int lastfile = static_cast<int> ((bigint) (me+1) * multiproc_nfile/nprocs);
    nreader = lastfile - firstfile;
    MPI_Comm_split(world,me,0,&clustercomm);
  } else if (multiproc_nfile < nprocs) {
    nreader = 1;
    int icluster = static_cast<int> ((bigint) me * multiproc_nfile/nprocs);
    firstfile = icluster;
    MPI_Comm_split(world,icluster,0,&clustercomm);
  }

  MPI_Comm_rank(clustercomm,&me_cluster);
  MPI_Comm_size(clustercomm,&nprocs_cluster);
  if (me_cluster == 0) filereader = 1;
  else filereader = 0;

  readers = new Reader*[nreader];
  nsnapatoms = new bigint[nreader];
  for (int i=0; i < nreader; ++i) {
    readers[i] = nullptr;
    nsnapatoms[i] = 0;
  }

  // create Nreader reader classes per reader
  // match readerstyle to options in style_reader.h

  if (0) {
    return;        // dummy line to enable else-if macro expansion

#define READER_CLASS
#define ReaderStyle(key,Class) \
  } else if (strcmp(readerstyle,#key) == 0) { \
    for (int i = 0; i < nreader; i++) { \
      readers[i] = new Class(lmp); \
    }
#include "style_reader.h"       // IWYU pragma: keep
#undef READER_CLASS

  // unrecognized style

  } else error->all(FLERR,utils::check_packages_for_style("reader",readerstyle,lmp));

  if (utils::strmatch(readerstyle,"^adios")) {
      // everyone is a reader with adios
      parallel = 1;
      filereader = 1;
  }

  // pass any arguments to readers

  if (narg > 0 && filereader)
    for (int i = 0; i < nreader; i++)
      readers[i]->settings(narg,arg);
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

  // proc 0 finds the timestep in its first reader

  if (me == 0 || parallel) {

    // exit file loop when dump timestep >= nrequest
    // or files exhausted

    for (ifile = 0; ifile < nfile; ifile++) {
      ntimestep = -1;
      if (multiproc) {
        std::string multiname = files[ifile];
        multiname.replace(multiname.find('%'),1,"0");
        readers[0]->open_file(multiname.c_str());
      } else readers[0]->open_file(files[ifile]);

      while (1) {
        eofflag = readers[0]->read_time(ntimestep);
        if (eofflag) break;
        if (ntimestep >= nrequest) break;
        readers[0]->skip();
      }

      if (ntimestep >= nrequest) break;
      readers[0]->close_file();
    }

    currentfile = ifile;
    if (ntimestep < nrequest) ntimestep = -1;
    if (exact && ntimestep != nrequest) ntimestep = -1;
  }

  if (!parallel) {
    // proc 0 broadcasts timestep and currentfile to all procs

    MPI_Bcast(&ntimestep,1,MPI_LMP_BIGINT,0,world);
    MPI_Bcast(&currentfile,1,MPI_INT,0,world);
  }

  // if ntimestep < 0:
  // all filereader procs close all their files and return

  if (ntimestep < 0) {
    if (filereader)
      for (int i = 0; i < nreader; i++)
        readers[i]->close_file();
    return ntimestep;
  }

  // for multiproc mode:
  // all filereader procs search for same ntimestep in currentfile

  if (multiproc && filereader) {
    for (int i = 0; i < nreader; i++) {
      if (me == 0 && i == 0) continue;    // proc 0, reader 0 already found it
      std::string multiname = files[currentfile];
      multiname.replace(multiname.find('%'),1,fmt::format("{}",firstfile+i));
      readers[i]->open_file(multiname.c_str());

      bigint step;
      while (1) {
        eofflag = readers[i]->read_time(step);
        if (eofflag) break;
        if (step == ntimestep) break;
        readers[i]->skip();
      }

      if (eofflag)
        error->one(FLERR,"Read dump parallel files "
                   "do not all have same timestep");
    }
  }

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

  // proc 0 finds the timestep in its first reader

  if (me == 0 || parallel) {

    // exit file loop when dump timestep matches all criteria
    // or files exhausted

    int iskip = 0;

    for (ifile = currentfile; ifile < nfile; ifile++) {
      ntimestep = -1;
      if (ifile != currentfile) {
        if (multiproc) {
          std::string multiname = files[ifile];
          multiname.replace(multiname.find('%'),1,"0");
          readers[0]->open_file(multiname.c_str());
        } else readers[0]->open_file(files[ifile]);
      }

      while (1) {
        eofflag = readers[0]->read_time(ntimestep);
        if (eofflag) break;
        if (ntimestep > nlast) break;
        if (ntimestep <= ncurrent) {
          readers[0]->skip();
          continue;
        }
        if (iskip == nskip) iskip = 0;
        iskip++;
        if (nevery && ntimestep % nevery) readers[0]->skip();
        else if (iskip < nskip) readers[0]->skip();
        else break;
      }

      if (eofflag) readers[0]->close_file();
      else break;
    }

    currentfile = ifile;
    if (eofflag) ntimestep = -1;
    if (ntimestep <= ncurrent) ntimestep = -1;
    if (ntimestep > nlast) ntimestep = -1;
  }

  if (!parallel) {
    // proc 0 broadcasts timestep and currentfile to all procs

    MPI_Bcast(&ntimestep,1,MPI_LMP_BIGINT,0,world);
    MPI_Bcast(&currentfile,1,MPI_INT,0,world);
  }

  // if ntimestep < 0:
  // all filereader procs close all their files and return

  if (ntimestep < 0) {
    if (filereader)
      for (int i = 0; i < nreader; i++)
        readers[i]->close_file();
    return ntimestep;
  }

  // for multiproc mode:
  // all filereader procs search for same ntimestep in currentfile

  if (multiproc && filereader) {
    for (int i = 0; i < nreader; i++) {
      if (me == 0 && i == 0) continue;
      std::string multiname = files[currentfile];
      multiname.replace(multiname.find('%'),1,fmt::format("{}",firstfile+i));
      readers[i]->open_file(multiname.c_str());

      bigint step;
      while (1) {
        eofflag = readers[i]->read_time(step);
        if (eofflag) break;
        if (step == ntimestep) break;
        readers[i]->skip();
      }

      if (eofflag)
        error->one(FLERR,"Read dump parallel files "
                   "do not all have same timestep");
    }
  }

  return ntimestep;
}

/* ----------------------------------------------------------------------
   read and broadcast and store snapshot header info
   set nsnapatoms = # of atoms in snapshot
------------------------------------------------------------------------- */

void ReadDump::header(int fieldinfo)
{
  int boxinfo, triclinic_snap;
  int fieldflag,xflag,yflag,zflag;

  if (filereader) {
    for (int i = 0; i < nreader; i++)
      nsnapatoms[i] = readers[i]->read_header(box,boxinfo,triclinic_snap,fieldinfo,
                                              nfield,fieldtype,fieldlabel,
                                              scaleflag,wrapflag,fieldflag,
                                              xflag,yflag,zflag);
  }

  if (!parallel) {
    MPI_Bcast(nsnapatoms,nreader,MPI_LMP_BIGINT,0,clustercomm);
    MPI_Bcast(&boxinfo,1,MPI_INT,0,clustercomm);
    MPI_Bcast(&triclinic_snap,1,MPI_INT,0,clustercomm);
    MPI_Bcast(&box[0][0],9,MPI_DOUBLE,0,clustercomm);
  }

  // local copy of snapshot box parameters
  // used in xfield,yfield,zfield when converting dump atom to absolute coords

  if (boxinfo) {
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
  }

  // done if not checking fields

  if (!fieldinfo) return;

  MPI_Bcast(&fieldflag,1,MPI_INT,0,clustercomm);
  MPI_Bcast(&xflag,1,MPI_INT,0,clustercomm);
  MPI_Bcast(&yflag,1,MPI_INT,0,clustercomm);
  MPI_Bcast(&zflag,1,MPI_INT,0,clustercomm);

  // error check on current vs new box and fields
  // boxinfo == 0 means no box info in file

  if (boxflag) {
    if (!boxinfo)
      error->all(FLERR,"No box information in dump, must use 'box no'");
    else if ((triclinic_snap && !triclinic) ||
             (!triclinic_snap && triclinic))
      error->one(FLERR,"Read_dump triclinic status does not match simulation");
  }

  // error check on requested fields existing in dump file

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

/* ----------------------------------------------------------------------
   read and process one snapshot of atoms
------------------------------------------------------------------------- */

void ReadDump::atoms()
{
  // initialize counters

  npurge = nreplace = ntrim = nadd = 0;

  // if purgeflag set, delete all current atoms

  if (purgeflag) {
    if (atom->map_style != Atom::MAP_NONE) atom->map_clear();
    npurge = atom->nlocal;
    atom->nlocal = atom->nghost = 0;
    atom->natoms = 0;
  }

  // read all the snapshot atoms into fields
  // each proc will own an arbitrary subset of atoms

  read_atoms();

  // migrate old owned atoms to new procs based on atom IDs
  // not necessary if purged all old atoms or if only 1 proc

  if (!purgeflag && nprocs > 1) migrate_old_atoms();

  // migrate new snapshot atoms to same new procs based on atom IDs
  // not necessary if purged all old atoms or if only 1 proc

  if (!purgeflag && nprocs > 1) migrate_new_atoms();

  // must build map if not a molecular system
  // this will be needed to match new atoms to old atoms

  int mapflag = 0;
  if (atom->map_style == Atom::MAP_NONE) {
    mapflag = 1;
    atom->map_init();
    atom->map_set();
  }

  // each proc now owns both old and new info for same subset of atoms
  // update each local atom with new info

  process_atoms();

  // check that atom IDs are valid

  atom->tag_check();

  // delete atom map if created it above
  // else reinitialize map for current atoms
  // do this before migrating atoms to new procs via Irregular

  if (mapflag) {
    atom->map_delete();
    atom->map_style = Atom::MAP_NONE;
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

  // migrate atoms to their new owing proc, based on atom coords

  migrate_atoms_by_coords();
}

/* ----------------------------------------------------------------------
   read all the snapshot atoms into fields
   done in different ways for multiproc no/yes and # of procs < or >= nprocs
   nnew = # of snapshot atoms this proc stores
------------------------------------------------------------------------- */

void ReadDump::read_atoms()
{
  int count,nread,nsend,nrecv,otherproc;
  bigint nsnap,ntotal,ofirst,olast,rfirst,rlast,lo,hi;
  MPI_Request request;
  MPI_Status status;

  // one reader per cluster of procs
  // each reading proc reads one file and splits data across cluster
  // cluster can be all procs or a subset

  if (!parallel && (!multiproc || multiproc_nfile < nprocs)) {
    nsnap = nsnapatoms[0];

    if (filereader) {
      if (!buf) memory->create(buf,CHUNK,nfield,"read_dump:buf");

      otherproc = 0;
      ofirst = (bigint) otherproc * nsnap/nprocs_cluster;
      olast = (bigint) (otherproc+1) * nsnap/nprocs_cluster;
      if (olast-ofirst > MAXSMALLINT)
        error->one(FLERR,"Read dump snapshot is too large for a proc");
      nnew = static_cast<int> (olast - ofirst);

      if (nnew > maxnew || maxnew == 0) {
        memory->destroy(fields);
        maxnew = MAX(nnew,1);    // avoid null pointer
        memory->create(fields,maxnew,nfield,"read_dump:fields");
      }

      ntotal = 0;
      while (ntotal < nsnap) {
        nread = MIN(CHUNK,nsnap-ntotal);
        readers[0]->read_atoms(nread,nfield,buf);
        rfirst = ntotal;
        rlast = ntotal + nread;

        nsend = 0;
        while (nsend < nread) {
          lo = MAX(ofirst,rfirst);
          hi = MIN(olast,rlast);
          if (otherproc)    // send to otherproc or copy to self
            MPI_Send(&buf[nsend][0],(hi-lo)*nfield,MPI_DOUBLE,
                     otherproc,0,clustercomm);
          else
            memcpy(&fields[rfirst][0],&buf[nsend][0],
                   (hi-lo)*nfield*sizeof(double));
          nsend += hi-lo;
          if (hi == olast) {
            otherproc++;
            ofirst = (bigint) otherproc * nsnap/nprocs_cluster;
            olast = (bigint) (otherproc+1) * nsnap/nprocs_cluster;
          }
        }

        ntotal += nread;
      }

    } else {
      ofirst = (bigint) me_cluster * nsnap/nprocs_cluster;
      olast = (bigint) (me_cluster+1) * nsnap/nprocs_cluster;
      if (olast-ofirst > MAXSMALLINT)
        error->one(FLERR,"Read dump snapshot is too large for a proc");
      nnew = static_cast<int> (olast - ofirst);
      if (nnew > maxnew || maxnew == 0) {
        memory->destroy(fields);
        maxnew = MAX(nnew,1);     // avoid null pointer
        memory->create(fields,maxnew,nfield,"read_dump:fields");
      }

      nrecv = 0;
      while (nrecv < nnew) {
        MPI_Irecv(&fields[nrecv][0],(nnew-nrecv)*nfield,MPI_DOUBLE,0,0,
                  clustercomm,&request);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&count);
        nrecv += count/nfield;
      }
    }

  // every proc is a filereader, reads one or more files
  // each proc keeps all data it reads, no communication required

  } else if (multiproc_nfile >= nprocs || parallel) {
    bigint sum = 0;
    for (int i = 0; i < nreader; i++)
      sum += nsnapatoms[i];
    if (sum > MAXSMALLINT)
      error->one(FLERR,"Read dump snapshot is too large for a proc");
    nnew = static_cast<int> (sum);
    if (nnew > maxnew || maxnew == 0) {
      memory->destroy(fields);
      maxnew = MAX(nnew,1);     // avoid null pointer
      memory->create(fields,maxnew,nfield,"read_dump:fields");
    }

    nnew = 0;
    for (int i = 0; i < nreader; i++) {
      nsnap = nsnapatoms[i];
      ntotal = 0;
      while (ntotal < nsnap) {
        if (parallel) {
          // read the whole thing at once
          nread = nsnap-ntotal;
        } else {
          nread = MIN(CHUNK,nsnap-ntotal);
        }
        readers[i]->read_atoms(nread,nfield,&fields[nnew+ntotal]);
        ntotal += nread;
      }
      nnew += nsnap;
    }
  }
}

/* ----------------------------------------------------------------------
   update info for each old atom I own based on snapshot info
   if in replace mode and atom ID matches current atom,
     overwrite atom info with fields from dump file
   if in add mode and atom ID does not match any old atom,
     create new atom with dump file field values
------------------------------------------------------------------------- */

void ReadDump::process_atoms()
{
  int i,m,ifield,itype;
  int xbox,ybox,zbox;
  tagint mtag;
  int *updateflag,*newflag;

  // updateflag[i] = flag for old atoms, 1 if updated, else 0
  // newflag[i] = flag for new atoms, 0 if used to update old atom, else 1

  int nlocal = atom->nlocal;
  memory->create(updateflag,nlocal,"read_dump:updateflag");
  for (i = 0; i < nlocal; i++) updateflag[i] = 0;
  memory->create(newflag,nnew,"read_dump:newflag");
  for (i = 0; i < nnew; i++) newflag[i] = 1;

  // loop over new atoms

  double **x = atom->x;
  double **v = atom->v;
  double *q = atom->q;
  double **f = atom->f;
  tagint *tag = atom->tag;
  imageint *image = atom->image;
  tagint map_tag_max = atom->map_tag_max;

  for (i = 0; i < nnew; i++) {

    // check if new atom matches one I own
    // setting m = -1 forces new atom not to match
    // NOTE: atom ID in fields is stored as double, not as ubuf
    //       so can only cast it to tagint, thus cannot be full 64-bit ID

    mtag = static_cast<tagint> (fields[i][0]);
    if (mtag <= map_tag_max) m = atom->map(mtag);
    else m = -1;
    if (m < 0 || m >= nlocal) continue;

    updateflag[m] = 1;
    newflag[i] = 0;

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
        case FX:
          f[m][0] = fields[i][ifield];
          break;
        case FY:
          f[m][1] = fields[i][ifield];
          break;
        case FZ:
          f[m][2] = fields[i][ifield];
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

  // if trimflag set, delete atoms not updated by snapshot atoms

  if (trimflag) {
    AtomVec *avec = atom->avec;

    i = 0;
    while (i < nlocal) {
      if (!updateflag[i]) {
        avec->copy(nlocal-1,i,1);
        updateflag[i] = updateflag[nlocal-1];
        nlocal--;
        ntrim++;
      } else i++;
    }

    atom->nlocal = nlocal;
    bigint nblocal = atom->nlocal;
    MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }

  // done if cannot add new atoms

  if (addflag == NOADD) {
    memory->destroy(updateflag);
    memory->destroy(newflag);
    return;
  }

  // ----------------------------------------------------
  // create new atoms for dump file atoms with ID that matches no old atom
  // ----------------------------------------------------

  // first check that dump file snapshot has atom type field

  int tflag = 0;
  for (ifield = 0; ifield < nfield; ifield++)
    if (fieldtype[ifield] == TYPE) tflag = 1;
  if (!tflag)
    error->all(FLERR,"Cannot add atoms if dump file does not store atom type");

  int nlocal_previous = atom->nlocal;
  double one[3];

  for (i = 0; i < nnew; i++) {
    if (!newflag[i]) continue;

    // create type and coord fields from dump file
    // coord = 0.0 unless corresponding dump file field was specified

    itype = 0;
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

    tag = atom->tag;
    v = atom->v;
    q = atom->q;
    image = atom->image;

    // set atom attributes from other dump file fields

    xbox = ybox = zbox = 0;

    for (ifield = 0; ifield < nfield; ifield++) {
      switch (fieldtype[ifield]) {
      case ID:
        if (addflag == KEEPADD)
          tag[m] = static_cast<tagint> (fields[i][ifield]);
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

      // reset image flag in case changed by ix,iy,iz fields

      image[m] = ((imageint) (xbox + IMGMAX) & IMGMASK) |
        (((imageint) (ybox + IMGMAX) & IMGMASK) << IMGBITS) |
        (((imageint) (zbox + IMGMAX) & IMGMASK) << IMG2BITS);
    }
  }

  // if addflag = YESADD or KEEPADD, update total atom count

  if (addflag == YESADD || addflag == KEEPADD) {
    bigint nblocal = atom->nlocal;
    MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }

  // if addflag = YESADD,
  // assign consistent IDs to new snapshot atoms across all procs

  if (addflag == YESADD) {
    if (atom->natoms < 0 || atom->natoms >= MAXBIGINT)
      error->all(FLERR,"Too many total atoms");
    if (atom->tag_enable) atom->tag_extend();
  }

  // init per-atom fix/compute/variable values for created atoms

  atom->data_fix_compute_variable(nlocal_previous,atom->nlocal);

  // free allocated vectors

  memory->destroy(updateflag);
  memory->destroy(newflag);
}

/* ----------------------------------------------------------------------
   migrate old atoms to new procs based on atom IDs
   use migrate_atoms() with explicit processor assignments
------------------------------------------------------------------------- */

void ReadDump::migrate_old_atoms()
{
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int *procassign;
  memory->create(procassign,nlocal,"read_dump:procassign");
  for (int i = 0; i < nlocal; i++)
    procassign[i] = tag[i] % nprocs;

  Irregular *irregular = new Irregular(lmp);
  irregular->migrate_atoms(1,1,procassign);
  delete irregular;

  memory->destroy(procassign);
}

/* ----------------------------------------------------------------------
   migrate new atoms to same new procs based on atom IDs
------------------------------------------------------------------------- */

void ReadDump::migrate_new_atoms()
{
  tagint mtag;
  int *procassign;
  double **newfields;

  memory->create(procassign,nnew,"read_dump:procassign");
  for (int i = 0; i < nnew; i++) {
    mtag = static_cast<tagint> (fields[i][0]);
    procassign[i] = mtag % nprocs;
  }

  Irregular *irregular = new Irregular(lmp);
  int nrecv = irregular->create_data(nnew,procassign,1);
  int newmaxnew = MAX(nrecv,maxnew);
  newmaxnew = MAX(newmaxnew,1);    // avoid null pointer
  memory->create(newfields,newmaxnew,nfield,"read_dump:newfields");
  irregular->exchange_data((char *) &fields[0][0],nfield*sizeof(double),
                           (char *) &newfields[0][0]);
  irregular->destroy_data();
  delete irregular;

  memory->destroy(fields);
  memory->destroy(procassign);

  // point fields at newfields

  fields = newfields;
  maxnew = newmaxnew;
  nnew = nrecv;
}

/* ----------------------------------------------------------------------
   migrate final atoms to new procs based on atom coords
   use migrate_atoms() with implicit processor assignments based on atom coords
   move atoms back inside simulation box and to new processors
   use remap() instead of pbc() in case atoms moved a long distance
   adjust image flags of all atoms (old and new) based on current box
------------------------------------------------------------------------- */

void ReadDump::migrate_atoms_by_coords()
{
  double **x = atom->x;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) domain->remap(x[i],image[i]);

  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->reset_box();
  Irregular *irregular = new Irregular(lmp);
  irregular->migrate_atoms(1);
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
  // scan ahead to see if "add yes/keep" keyword/value is used
  // requires extra "type" field from from dump file

  int iarg;
  for (iarg = 0; iarg < narg; iarg++)
    if (strcmp(arg[iarg],"add") == 0)
      if (iarg < narg-1 && (strcmp(arg[iarg+1],"yes") == 0 ||
                            strcmp(arg[iarg+1],"keep") == 0)) break;

  nfield = 0;
  fieldtype[nfield++] = ID;
  if (iarg < narg) fieldtype[nfield++] = TYPE;

  // parse fields

  iarg = 0;
  while (iarg < narg) {
    int type = whichtype(arg[iarg]);
    if (type < 0) break;
    if (type == Q && !atom->q_flag)
      error->all(FLERR,"Read dump of atom property that isn't allocated");
    fieldtype[nfield++] = type;
    iarg++;
  }

  // check for no fields

  if (fieldtype[nfield-1] == ID || fieldtype[nfield-1] == TYPE)
    error->all(FLERR,"Illegal read_dump command");

  if (dimension == 2) {
    for (int i = 0; i < nfield; i++)
      if (fieldtype[i] == Z || fieldtype[i] == VZ ||
          fieldtype[i] == IZ || fieldtype[i] == FZ)
        error->all(FLERR,"Illegal read_dump command");
  }

  for (int i = 0; i < nfield; i++)
    for (int j = i+1; j < nfield; j++)
      if (fieldtype[i] == fieldtype[j])
        error->all(FLERR,"Duplicate fields in read_dump command");

  // parse optional args

  multiproc_nfile = 0;
  boxflag = 1;
  replaceflag = 1;
  purgeflag = 0;
  trimflag = 0;
  addflag = NOADD;
  for (int i = 0; i < nfield; i++) fieldlabel[i] = nullptr;
  scaleflag = 0;
  wrapflag = 1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"nfile") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_dump command");
      multiproc_nfile = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"box") == 0) {
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
      if (strcmp(arg[iarg+1],"yes") == 0) addflag = YESADD;
      else if (strcmp(arg[iarg+1],"no") == 0) addflag = NOADD;
      else if (strcmp(arg[iarg+1],"keep") == 0) addflag = KEEPADD;
      else error->all(FLERR,"Illegal read_dump command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"label") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal read_dump command");
      int type = whichtype(arg[iarg+1]);
      int i;
      for (i = 0; i < nfield; i++)
        if (type == fieldtype[i]) break;
      if (i == nfield) error->all(FLERR,"Illegal read_dump command");
      fieldlabel[i] = utils::strdup(arg[iarg+2]);
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
      readerstyle = utils::strdup(arg[iarg+1]);
      iarg += 2;
      break;
    } else error->all(FLERR,"Illegal read_dump command");
  }

  if (multiproc == 0 && multiproc_nfile)
    error->all(FLERR,"Dump file is not a multi-proc file");
  if (multiproc && multiproc_nfile == 0)
    error->all(FLERR,"Dump file is a multi-proc file");

  if (purgeflag && (replaceflag || trimflag))
    error->all(FLERR,"If read_dump purges it cannot replace or trim");
  if (addflag == KEEPADD && atom->tag_enable == 0)
    error->all(FLERR,"Read_dump cannot use 'add keep' without atom IDs");

  return narg-iarg;
}

/* ----------------------------------------------------------------------
   check if str is a field argument
   if yes, return index of which
   if not, return -1
------------------------------------------------------------------------- */

int ReadDump::whichtype(char *str)
{
  int type = -1;
  if (strcmp(str,"id") == 0) type = ID;
  else if (strcmp(str,"type") == 0) type = TYPE;
  else if (strcmp(str,"x") == 0) type = X;
  else if (strcmp(str,"y") == 0) type = Y;
  else if (strcmp(str,"z") == 0) type = Z;
  else if (strcmp(str,"vx") == 0) type = VX;
  else if (strcmp(str,"vy") == 0) type = VY;
  else if (strcmp(str,"vz") == 0) type = VZ;
  else if (strcmp(str,"q") == 0) type = Q;
  else if (strcmp(str,"ix") == 0) type = IX;
  else if (strcmp(str,"iy") == 0) type = IY;
  else if (strcmp(str,"iz") == 0) type = IZ;
  else if (strcmp(str,"fx") == 0) type = FX;
  else if (strcmp(str,"fy") == 0) type = FY;
  else if (strcmp(str,"fz") == 0) type = FZ;
  return type;
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
