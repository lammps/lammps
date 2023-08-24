// clang-format off
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

#include "dump.h"

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "group.h"
#include "input.h"
#include "irregular.h"
#include "memory.h"
#include "modify.h"
#include "output.h"
#include "update.h"
#include "variable.h"

#include <cstring>
#include <stdexcept>

using namespace LAMMPS_NS;

#if defined(LMP_QSORT)
// allocate space for static class variable
Dump *Dump::dumpptr;
#endif

#define BIG 1.0e20
#define EPSILON 1.0e-6

enum { ASCEND, DESCEND };

/* ---------------------------------------------------------------------- */

Dump::Dump(LAMMPS *lmp, int /*narg*/, char **arg) :
    Pointers(lmp), multiname(nullptr), refresh(nullptr), skipvar(nullptr), format(nullptr),
    format_default(nullptr), format_line_user(nullptr), format_float_user(nullptr),
    format_int_user(nullptr), format_bigint_user(nullptr), format_column_user(nullptr), fp(nullptr),
    nameslist(nullptr), buf(nullptr), sbuf(nullptr), ids(nullptr), bufsort(nullptr),
    idsort(nullptr), index(nullptr), proclist(nullptr), xpbc(nullptr), vpbc(nullptr),
    imagepbc(nullptr), irregular(nullptr)
{
  MPI_Comm_rank(world, &me);
  MPI_Comm_size(world, &nprocs);

  id = utils::strdup(arg[0]);

  igroup = group->find(arg[1]);
  groupbit = group->bitmask[igroup];

  style = utils::strdup(arg[2]);

  filename = utils::strdup(arg[4]);

  comm_forward = comm_reverse = 0;

  first_flag = 0;
  flush_flag = 1;

  refreshflag = 0;

  clearstep = 0;
  sort_flag = 0;
  balance_flag = 0;
  append_flag = 0;
  buffer_allow = 0;
  buffer_flag = 0;
  padflag = 0;
  pbcflag = 0;
  time_flag = 0;
  unit_flag = 0;
  unit_count = 0;
  delay_flag = 0;
  write_header_flag = 1;
  has_id = 1;

  skipflag = 0;

  maxfiles = -1;
  numfiles = 0;
  fileidx = 0;

  maxbuf = maxids = maxsort = maxproc = 0;
  maxsbuf = 0;

  maxpbc = -1;

  // parse filename for special syntax
  // if contains '%', write one file per proc and replace % with proc-ID
  // if contains '*', write one file per timestep and replace * with timestep
  // check file suffixes
  //   if ends in .bin -> binary file
  //   else if ends in .gz or other known extensions -> compressed text file
  //   else ASCII text file

  singlefile_opened = 0;
  compressed = 0;
  binary = 0;
  multifile = 0;
  size_one = 0;

  multiproc = 0;
  nclusterprocs = nprocs;
  filewriter = 0;
  if (me == 0) filewriter = 1;
  fileproc = 0;

  char *ptr;
  if ((ptr = strchr(filename,'%'))) {
    multiproc = 1;
    nclusterprocs = 1;
    filewriter = 1;
    fileproc = me;
    MPI_Comm_split(world,me,0,&clustercomm);
    *ptr = '\0';
    multiname = utils::strdup(fmt::format("{}{}{}", filename, me, ptr+1));
    *ptr = '%';
  }

  if (strchr(filename,'*')) multifile = 1;

  if (utils::strmatch(filename, "\\.bin$")
      || utils::strmatch(filename, "\\.lammpsbin$")) binary = 1;
  if (platform::has_compress_extension(filename)) compressed = 1;
}

/* ---------------------------------------------------------------------- */

Dump::~Dump()
{
  delete[] id;
  delete[] style;
  delete[] filename;
  delete[] multiname;

  delete[] format;
  delete[] format_default;
  delete[] format_line_user;
  delete[] format_float_user;
  delete[] format_int_user;
  delete[] format_bigint_user;

  delete[] refresh;
  delete[] skipvar;

  // format_column_user is deallocated by child classes that use it

  memory->destroy(buf);
  memory->destroy(bufsort);
  memory->destroy(ids);
  memory->destroy(idsort);
  memory->destroy(index);
  memory->destroy(proclist);
  delete irregular;

  memory->destroy(sbuf);

  if (pbcflag) {
    memory->destroy(xpbc);
    memory->destroy(vpbc);
    memory->destroy(imagepbc);
  }

  if (multiproc) MPI_Comm_free(&clustercomm);

  // delete storage for caching file names

  if (maxfiles > 0) {
    for (int idx=0; idx < numfiles; ++idx)
      delete[] nameslist[idx];
    delete[] nameslist;
  }

  // XTC style sets fp to a null pointer since it closes file in its destructor

  if (multifile == 0 && fp != nullptr) {
    if (compressed) {
      if (filewriter) platform::pclose(fp);
    } else {
      if (filewriter) fclose(fp);
    }
    fp = nullptr;
  }
}

/* ---------------------------------------------------------------------- */

void Dump::init()
{
  init_style();

  if (!sort_flag) {
    memory->destroy(bufsort);
    memory->destroy(ids);
    memory->destroy(idsort);
    memory->destroy(index);
    memory->destroy(proclist);
    delete irregular;

    maxids = maxsort = maxproc = 0;
    bufsort = nullptr;
    ids = idsort = nullptr;
    index = proclist = nullptr;
    irregular = nullptr;
    if ((has_id == 0) && (me == 0))
      error->warning(FLERR,"Dump {} includes no atom IDs and is not sorted by ID. This may complicate "
                     "post-processing tasks or visualization", id);
  }

  if (sort_flag) {
    if (multiproc > 1)
      error->all(FLERR,
                 "Cannot sort dump when 'nfile' or 'fileper' keywords are set to non-default values");
    if (sortcol == 0 && atom->tag_enable == 0)
      error->all(FLERR,"Cannot sort dump on atom IDs with no atom IDs defined");
    if (sortcol && sortcol > size_one)
      error->all(FLERR,"Dump sort column is invalid");
    if ((sortcol != 0) && (has_id == 0) && (me == 0))
      error->warning(FLERR,"Dump {} includes no atom IDs and is not sorted by ID. This may complicate "
                     "post-processing tasks or visualization", id);
    if (nprocs > 1 && irregular == nullptr)
      irregular = new Irregular(lmp);

    bigint size = group->count(igroup);

    // set reorderflag = 1 if can simply reorder local atoms rather than sort
    // criteria: sorting by ID, atom IDs are consecutive from 1 to Natoms
    //           min/max IDs of group match size of group
    // compute ntotal_reorder, nme_reorder, idlo/idhi to test against later

    reorderflag = 0;

    int gcmcflag = 0;
    for (const auto &fix : modify->get_fix_list())
      if (utils::strmatch(fix->style,"^gcmc"))
        gcmcflag = 1;

    if (sortcol == 0 && atom->tag_consecutive() && !gcmcflag) {
      tagint *tag = atom->tag;
      int *mask = atom->mask;
      int nlocal = atom->nlocal;

      tagint min = MAXTAGINT;
      tagint max = 0;
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          min = MIN(min,tag[i]);
          max = MAX(max,tag[i]);
        }
      tagint minall,maxall;
      MPI_Allreduce(&min,&minall,1,MPI_LMP_TAGINT,MPI_MIN,world);
      MPI_Allreduce(&max,&maxall,1,MPI_LMP_TAGINT,MPI_MAX,world);

      if (maxall-minall+1 == size) {
        reorderflag = 1;
        double range = maxall-minall + EPSILON;
        idlo = static_cast<tagint> (range*me/nprocs + minall);
        auto  idhi = static_cast<tagint> (range*(me+1)/nprocs + minall);

        auto  lom1 = static_cast<tagint> ((idlo-1-minall)/range * nprocs);
        auto  lo = static_cast<tagint> ((idlo-minall)/range * nprocs);
        auto  him1 = static_cast<tagint> ((idhi-1-minall)/range * nprocs);
        auto  hi = static_cast<tagint> ((idhi-minall)/range * nprocs);
        if (me && me == lom1) idlo--;
        else if (me && me != lo) idlo++;
        if (me+1 == him1) idhi--;
        else if (me+1 != hi) idhi++;

        nme_reorder = idhi-idlo;
        ntotal_reorder = size;
      }
    }
  }

  // search for refresh compute specified by dump_modify refresh

  if (refreshflag) {
    int icompute;
    for (icompute = 0; icompute < modify->ncompute; icompute++)
      if (strcmp(refresh,modify->compute[icompute]->id) == 0) break;
    if (icompute < modify->ncompute) irefresh = icompute;
    else error->all(FLERR,"Dump could not find refresh compute ID");
  }

  // if skipflag, check skip variable

  if (skipflag) {
    skipindex = input->variable->find(skipvar);
    if (skipindex < 0) error->all(FLERR,"Dump skip variable not found");
    if (!input->variable->equalstyle(skipindex))
      error->all(FLERR,"Variable for dump skip is invalid style");
  }

  // preallocation for PBC copies if requested

  if (pbcflag && atom->nlocal > maxpbc) pbc_allocate();
}

/* ---------------------------------------------------------------------- */

int Dump::count()
{
  if (igroup == 0) return atom->nlocal;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) m++;
  return m;
}

/* ---------------------------------------------------------------------- */

void Dump::write()
{
  imageint *imagehold;
  double **xhold,**vhold;

  // simulation box bounds

  if (domain->triclinic == 0) {
    boxxlo = domain->boxlo[0];
    boxxhi = domain->boxhi[0];
    boxylo = domain->boxlo[1];
    boxyhi = domain->boxhi[1];
    boxzlo = domain->boxlo[2];
    boxzhi = domain->boxhi[2];
  } else {
    boxxlo = domain->boxlo_bound[0];
    boxxhi = domain->boxhi_bound[0];
    boxylo = domain->boxlo_bound[1];
    boxyhi = domain->boxhi_bound[1];
    boxzlo = domain->boxlo_bound[2];
    boxzhi = domain->boxhi_bound[2];
    boxxy = domain->xy;
    boxxz = domain->xz;
    boxyz = domain->yz;
  }

  // nme = # of dump lines this proc contributes to dump

  nme = count();

  // if timestep < delaystep, just return
  // if skip condition is defined and met, just return
  // must do both these tests after count() b/c it invokes computes,
  //   this enables caller to trigger future invocation of needed computes

  if (delay_flag && update->ntimestep < delaystep) return;

  if (skipflag) {
    double value = input->variable->compute_equal(skipindex);
    if (value != 0.0) return;
  }

  // if file per timestep, open new file
  // do this after skip check, so no file is opened if skip occurs

  if (multifile) openfile();
  if (fp) clearerr(fp);

  // ntotal = total # of dump lines in snapshot
  // nmax = max # of dump lines on any proc

  bigint bnme = nme;
  MPI_Allreduce(&bnme,&ntotal,1,MPI_LMP_BIGINT,MPI_SUM,world);

  int nmax;
  if (multiproc != nprocs) MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
  else nmax = nme;

  // ensure buf is sized for packing and communicating
  // use nmax to ensure filewriter proc can receive info from others
  // limit nmax*size_one to int since used as arg in MPI calls

  if (nmax*size_one > maxbuf) {
    if ((bigint) nmax * size_one > MAXSMALLINT)
      error->all(FLERR,"Too much per-proc info for dump");
    maxbuf = nmax * size_one;
    memory->destroy(buf);
    memory->create(buf,maxbuf,"dump:buf");
  }

  // ensure ids buffer is sized for sorting

  if (sort_flag && sortcol == 0 && nmax > maxids) {
    maxids = nmax;
    memory->destroy(ids);
    memory->create(ids,maxids,"dump:ids");
  }

  // apply PBC on copy of x,v,image if requested

  if (pbcflag) {
    int nlocal = atom->nlocal;
    if (nlocal > maxpbc) pbc_allocate();
    if (nlocal) {
      memcpy(&xpbc[0][0],&atom->x[0][0],3*nlocal*sizeof(double));
      memcpy(&vpbc[0][0],&atom->v[0][0],3*nlocal*sizeof(double));
      memcpy(imagepbc,atom->image,nlocal*sizeof(imageint));
    }
    xhold = atom->x;
    vhold = atom->v;
    imagehold = atom->image;
    atom->x = xpbc;
    atom->v = vpbc;
    atom->image = imagepbc;

    // for triclinic, PBC is applied in lamda coordinates

    if (domain->triclinic) domain->x2lamda(nlocal);
    domain->pbc();
    if (domain->triclinic) domain->lamda2x(nlocal);
  }

  // pack my data into buf
  // if sorting on IDs also request ID list from pack()
  // sort buf as needed

  if (sort_flag && (ntotal > 1) && sortcol == 0) pack(ids);
  else pack(nullptr);
  if (sort_flag && (ntotal > 1)) sort();
  if (balance_flag && (ntotal > 1)) balance();

  // write timestep header
  // for multiproc,
  //   nheader = # of lines in this file via Allreduce on clustercomm
  //   must come after sort, which can change nme

  bigint nheader = ntotal;
  if (multiproc) {
    bnme = nme;
    MPI_Allreduce(&bnme,&nheader,1,MPI_LMP_BIGINT,MPI_SUM,clustercomm);
  }

  if (filewriter && write_header_flag) write_header(nheader);

  // if buffering, convert doubles into strings
  // ensure sbuf is sized for communicating
  // cannot buffer if output is to binary file

  if (buffer_flag && !binary) {
    nsme = convert_string(nme,buf);
    int nsmin,nsmax;
    MPI_Allreduce(&nsme,&nsmin,1,MPI_INT,MPI_MIN,world);
    if (nsmin < 0) error->all(FLERR,"Too much buffered per-proc info for dump");
    if (multiproc != nprocs)
      MPI_Allreduce(&nsme,&nsmax,1,MPI_INT,MPI_MAX,world);
    else nsmax = nsme;
    if (nsmax > maxsbuf) {
      maxsbuf = nsmax;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }
  }

  // filewriter = 1 = this proc writes to file
  // ping each proc in my cluster, receive its data, write data to file
  // else wait for ping from fileproc, send my data to fileproc

  int tmp,nlines,nchars;
  MPI_Status status;
  MPI_Request request;

  // comm and output buf of doubles

  if (buffer_flag == 0 || binary) {
    if (filewriter) {
      for (int iproc = 0; iproc < nclusterprocs; iproc++) {
        if (iproc) {
          MPI_Irecv(buf,maxbuf,MPI_DOUBLE,me+iproc,0,world,&request);
          MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
          MPI_Wait(&request,&status);
          MPI_Get_count(&status,MPI_DOUBLE,&nlines);
          nlines /= size_one;
        } else nlines = nme;

        write_data(nlines,buf);
      }
      if (flush_flag && fp) fflush(fp);

    } else {
      MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
      MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,fileproc,0,world);
    }

  // comm and output sbuf = one big string of formatted values per proc

  } else {
    if (filewriter) {
      for (int iproc = 0; iproc < nclusterprocs; iproc++) {
        if (iproc) {
          MPI_Irecv(sbuf,maxsbuf,MPI_CHAR,me+iproc,0,world,&request);
          MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
          MPI_Wait(&request,&status);
          MPI_Get_count(&status,MPI_CHAR,&nchars);
        } else nchars = nsme;

        write_data(nchars,(double *) sbuf);
      }
      if (flush_flag && fp) fflush(fp);

    } else {
      MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
      MPI_Rsend(sbuf,nsme,MPI_CHAR,fileproc,0,world);
    }
  }

  // restore original x,v,image unaltered by PBC

  if (pbcflag) {
    atom->x = xhold;
    atom->v = vhold;
    atom->image = imagehold;
  }

  // trigger post-dump refresh by specified compute
  // currently used for incremental dump files

  if (refreshflag) modify->compute[irefresh]->refresh();

  if (filewriter && fp != nullptr) write_footer();

  if (fp && ferror(fp)) error->one(FLERR,"Error writing dump {}: {}", id, utils::getsyserror());

  // if file per timestep, close file if I am filewriter

  if (multifile) {
    if (compressed) {
      if (filewriter && fp != nullptr) platform::pclose(fp);
    } else {
      if (filewriter && fp != nullptr) fclose(fp);
    }
    fp = nullptr;
  }
}

/* ----------------------------------------------------------------------
   generic opening of a dump file
   ASCII or binary or compressed
   some derived classes override this function
------------------------------------------------------------------------- */

void Dump::openfile()
{
  // single file, already opened, so just return

  if (singlefile_opened) return;
  if (multifile == 0) singlefile_opened = 1;

  unit_count = 0;

  // if one file per timestep, replace '*' with current timestep

  char *filecurrent = filename;
  if (multiproc) filecurrent = multiname;

  if (multifile) {
    filecurrent = utils::strdup(utils::star_subst(filecurrent, update->ntimestep, padflag));
    if (maxfiles > 0) {
      if (numfiles < maxfiles) {
        nameslist[numfiles] = utils::strdup(filecurrent);
        ++numfiles;
      } else {
        remove(nameslist[fileidx]);
        delete[] nameslist[fileidx];
        nameslist[fileidx] = utils::strdup(filecurrent);
        fileidx = (fileidx + 1) % maxfiles;
      }
    }
  }

  // each proc with filewriter = 1 opens a file

  if (filewriter) {
    if (compressed) {
      fp = platform::compressed_write(filecurrent);
    } else if (binary) {
      fp = fopen(filecurrent,"wb");
    } else if (append_flag) {
      fp = fopen(filecurrent,"a");
    } else {
      fp = fopen(filecurrent,"w");
    }

    if (fp == nullptr) error->one(FLERR,"Cannot open dump file");
  } else fp = nullptr;

  // delete string with timestep replaced

  if (multifile) delete[] filecurrent;
}

/* ----------------------------------------------------------------------
   parallel sort of buf across all procs
   changes nme, reorders datums in buf, grows buf if necessary
------------------------------------------------------------------------- */

void Dump::sort()
{
  int i,iproc;
  double value;

  // if single proc, swap ptrs to buf,ids <-> bufsort,idsort

  if (nprocs == 1) {
    if (nme > maxsort) {
      maxsort = nme;
      memory->destroy(bufsort);
      memory->create(bufsort,maxsort*size_one,"dump:bufsort");
      memory->destroy(index);
      memory->create(index,maxsort,"dump:index");
      if (sortcol == 0) {
        memory->destroy(idsort);
        memory->create(idsort,maxsort,"dump:idsort");
      }
    }

    double *dptr = buf;
    buf = bufsort;
    bufsort = dptr;

    if (sortcol == 0) {
      tagint *iptr = ids;
      ids = idsort;
      idsort = iptr;
    }

  // if multiple procs, exchange datums between procs via irregular

  } else {

    // grow proclist if necessary

    if (nme > maxproc) {
      maxproc = nme;
      memory->destroy(proclist);
      memory->create(proclist,maxproc,"dump:proclist");
    }

    // proclist[i] = which proc Ith datum will be sent to

    if (sortcol == 0) {
      tagint min = MAXTAGINT;
      tagint max = 0;
      for (i = 0; i < nme; i++) {
        min = MIN(min,ids[i]);
        max = MAX(max,ids[i]);
      }
      tagint minall,maxall;
      MPI_Allreduce(&min,&minall,1,MPI_LMP_TAGINT,MPI_MIN,world);
      MPI_Allreduce(&max,&maxall,1,MPI_LMP_TAGINT,MPI_MAX,world);

      // use 0.5 instead of EPSILON since atom IDs are integers
      // if use EPSILON, it can be lost if 64-bit maxall-minall is too big
      // then iproc == nprocs for largest ID, causing irregular to crash

      double range = maxall-minall + 0.5;
      for (i = 0; i < nme; i++) {
        iproc = static_cast<int> ((ids[i]-minall)/range * nprocs);
        proclist[i] = iproc;
      }

    } else {
      double min = BIG;
      double max = -BIG;
      for (i = 0; i < nme; i++) {
        value = buf[i*size_one + sortcolm1];
        min = MIN(min,value);
        max = MAX(max,value);
      }
      double minall,maxall;
      MPI_Allreduce(&min,&minall,1,MPI_DOUBLE,MPI_MIN,world);
      MPI_Allreduce(&max,&maxall,1,MPI_DOUBLE,MPI_MAX,world);
      double range = maxall-minall + EPSILON*(maxall-minall);
      if (range == 0.0) range = EPSILON;

      // proc assignment is inverted if sortorder = DESCEND

      for (i = 0; i < nme; i++) {
        value = buf[i*size_one + sortcolm1];
        iproc = static_cast<int> ((value-minall)/range * nprocs);
        if (sortorder == DESCEND) iproc = nprocs-1 - iproc;
        proclist[i] = iproc;
      }
    }

    // create comm plan, grow recv bufs if necessary,
    // exchange datums, destroy plan
    // if sorting on atom IDs, exchange IDs also

    nme = irregular->create_data(nme,proclist);

    if (nme > maxsort) {
      maxsort = nme;
      memory->destroy(bufsort);
      memory->create(bufsort,maxsort*size_one,"dump:bufsort");
      memory->destroy(index);
      memory->create(index,maxsort,"dump:index");
      if (sortcol == 0) {
        memory->destroy(idsort);
        memory->create(idsort,maxsort,"dump:idsort");
      }
    }

    irregular->exchange_data((char *) buf,size_one*sizeof(double),
                             (char *) bufsort);
    if (sortcol == 0)
      irregular->exchange_data((char *) ids,sizeof(tagint),(char *) idsort);
    irregular->destroy_data();
  }

  // if reorder flag is set & total/per-proc counts match pre-computed values,
  // then create index directly from idsort
  // else quicksort of index using IDs or buf column as comparator

  if (reorderflag) {
    if (ntotal != ntotal_reorder) reorderflag = 0;
    int flag = 0;
    if (nme != nme_reorder) flag = 1;
    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
    if (flagall) reorderflag = 0;

    if (reorderflag)
      for (i = 0; i < nme; i++)
        index[idsort[i]-idlo] = i;
  }

#if defined(LMP_QSORT)
  if (!reorderflag) {
    dumpptr = this;
    for (i = 0; i < nme; i++) index[i] = i;
    if (sortcol == 0) qsort(index,nme,sizeof(int),idcompare);
    else if (sortorder == ASCEND) qsort(index,nme,sizeof(int),bufcompare);
    else qsort(index,nme,sizeof(int),bufcompare_reverse);
  }
#else
  if (!reorderflag) {
    for (i = 0; i < nme; i++) index[i] = i;
    if (sortcol == 0) utils::merge_sort(index,nme,(void *)this,idcompare);
    else if (sortorder == ASCEND) utils::merge_sort(index,nme,(void *)this,bufcompare);
    else utils::merge_sort(index,nme,(void *)this,bufcompare_reverse);
  }
#endif

  // reset buf size and maxbuf to largest of any post-sort nme values
  // this ensures proc 0 can receive everyone's info

  int nmax;
  MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);

  if (nmax*size_one > maxbuf) {
    maxbuf = nmax * size_one;
    memory->destroy(buf);
    memory->create(buf,maxbuf,"dump:buf");
  }

  // copy data from bufsort to buf using index

  int nbytes = size_one*sizeof(double);
  for (i = 0; i < nme; i++)
    memcpy(&buf[i*size_one],&bufsort[index[i]*size_one],nbytes);
}

#if defined(LMP_QSORT)

/* ----------------------------------------------------------------------
   compare two atom IDs
   called via qsort() in sort() method
   is a static method so access data via dumpptr
------------------------------------------------------------------------- */

int Dump::idcompare(const void *pi, const void *pj)
{
  tagint *idsort = dumpptr->idsort;

  int i = *((int *) pi);
  int j = *((int *) pj);

  if (idsort[i] < idsort[j]) return -1;
  if (idsort[i] > idsort[j]) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   compare two buffer values with size_one stride
   called via qsort() in sort() method
   is a static method so access data via dumpptr
   sort in ASCENDing order
------------------------------------------------------------------------- */

int Dump::bufcompare(const void *pi, const void *pj)
{
  double *bufsort = dumpptr->bufsort;
  int size_one = dumpptr->size_one;
  int sortcolm1 = dumpptr->sortcolm1;

  int i = *((int *) pi)*size_one + sortcolm1;
  int j = *((int *) pj)*size_one + sortcolm1;

  if (bufsort[i] < bufsort[j]) return -1;
  if (bufsort[i] > bufsort[j]) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   compare two buffer values with size_one stride
   called via qsort() in sort() method
   is a static method so access data via dumpptr
   sort in DESCENDing order
------------------------------------------------------------------------- */

int Dump::bufcompare_reverse(const void *pi, const void *pj)
{
  double *bufsort = dumpptr->bufsort;
  int size_one = dumpptr->size_one;
  int sortcolm1 = dumpptr->sortcolm1;

  int i = *((int *) pi)*size_one + sortcolm1;
  int j = *((int *) pj)*size_one + sortcolm1;

  if (bufsort[i] > bufsort[j]) return -1;
  if (bufsort[i] < bufsort[j]) return 1;
  return 0;
}

#else

/* ----------------------------------------------------------------------
   compare two atom IDs
   called via merge_sort() in sort() method
------------------------------------------------------------------------- */

int Dump::idcompare(const int i, const int j, void *ptr)
{
  tagint *idsort = ((Dump *)ptr)->idsort;
  if (idsort[i] < idsort[j]) return -1;
  else if (idsort[i] > idsort[j]) return 1;
  else return 0;
}

/* ----------------------------------------------------------------------
   compare two buffer values with size_one stride
   called via merge_sort() in sort() method
   sort in ASCENDing order
------------------------------------------------------------------------- */

int Dump::bufcompare(const int i, const int j, void *ptr)
{
  auto dptr = (Dump *) ptr;
  double *bufsort     = dptr->bufsort;
  const int size_one  = dptr->size_one;
  const int sortcolm1 = dptr->sortcolm1;

  const int ii=i*size_one + sortcolm1;
  const int jj=j*size_one + sortcolm1;

  if (bufsort[ii] < bufsort[jj]) return -1;
  else if (bufsort[ii] > bufsort[jj]) return 1;
  else return 0;
}

/* ----------------------------------------------------------------------
   compare two buffer values with size_one stride
   called via merge_sort() in sort() method
   sort in DESCENDing order
------------------------------------------------------------------------- */

int Dump::bufcompare_reverse(const int i, const int j, void *ptr)
{
  auto dptr = (Dump *) ptr;
  double *bufsort     = dptr->bufsort;
  const int size_one  = dptr->size_one;
  const int sortcolm1 = dptr->sortcolm1;

  const int ii=i*size_one + sortcolm1;
  const int jj=j*size_one + sortcolm1;

  if (bufsort[ii] < bufsort[jj]) return 1;
  else if (bufsort[ii] > bufsort[jj]) return -1;
  else return 0;
}

#endif

/* ----------------------------------------------------------------------
   parallel load balance of buf across all procs
   must come after sort
------------------------------------------------------------------------- */

void Dump::balance()
{
  bigint *proc_offsets,*proc_new_offsets;
  memory->create(proc_offsets,nprocs+1,"dump:proc_offsets");
  memory->create(proc_new_offsets,nprocs+1,"dump:proc_new_offsets");

  // compute atom offset for this proc

  bigint offset;
  bigint bnme = nme;
  MPI_Scan(&bnme,&offset,1,MPI_LMP_BIGINT,MPI_SUM,world);

  // gather atom offsets for all procs

  MPI_Allgather(&offset,1,MPI_LMP_BIGINT,&proc_offsets[1],1,MPI_LMP_BIGINT,world);

  proc_offsets[0] = 0;

  // how many atoms should I own after balance

  int nme_balance = static_cast<int>(ntotal/nprocs);

  // include remainder atoms on first procs

  int remainder = ntotal % nprocs;
  if (me < remainder) nme_balance += 1;

  // compute new atom offset for this proc

  bigint offset_balance;
  bigint bnme_balance = nme_balance;
  MPI_Scan(&bnme_balance,&offset_balance,1,MPI_LMP_BIGINT,MPI_SUM,world);

  // gather new atom offsets for all procs

  MPI_Allgather(&offset_balance,1,MPI_LMP_BIGINT,&proc_new_offsets[1],1,MPI_LMP_BIGINT,world);

  proc_new_offsets[0] = 0;

  // reset buf size to largest of any post-balance nme values
  // this ensures proc 0 can receive everyone's info
  // cannot shrink buf to nme_balance, must use previous maxbuf value

  int nmax;
  MPI_Allreduce(&nme_balance,&nmax,1,MPI_INT,MPI_MAX,world);
  if (nmax*size_one > maxbuf) maxbuf = nmax*size_one;

  // allocate a second buffer for balanced data

  double *buf_balance;
  memory->create(buf_balance,maxbuf,"dump:buf_balance");

  // compute from which procs I am receiving atoms
  // post recvs first

  int nswap = 0;
  auto request = new MPI_Request[nprocs];

  // find which proc starting atom belongs to

  int startproc = me;
  while (proc_new_offsets[me] < proc_offsets[startproc]) startproc--;
  while (proc_new_offsets[me] > proc_offsets[startproc+1]-1) startproc++;

  // find which proc ending atom belongs to

  int endproc = me;
  while (proc_new_offsets[me] + nme_balance-1 < proc_offsets[endproc]) endproc--;
  while (proc_new_offsets[me] + nme_balance-1 > proc_offsets[endproc+1]-1) endproc++;

  // loop over procs

  for (int iproc = startproc; iproc <= endproc; iproc++) {
    int istart = MAX(0, proc_offsets[iproc] - proc_new_offsets[me]);
    int iend = MIN(nme_balance-1, proc_offsets[iproc+1]-1 - proc_new_offsets[me]);
    int nrecv = iend - istart + 1;
    if (nrecv == 0) continue;

    // post receive for this proc

    if (iproc != me)
      MPI_Irecv(&buf_balance[istart*size_one],nrecv*size_one,MPI_DOUBLE,
                iproc,0,world,&request[nswap++]);
  }

  // compute which atoms I am sending and to which procs

  // find which proc starting atom belongs to

  startproc = me;
  while (proc_offsets[me] < proc_new_offsets[startproc]) startproc--;
  while (proc_offsets[me] > proc_new_offsets[startproc+1]-1) startproc++;

  // find which proc ending atom belongs to

  endproc = me;
  while (proc_offsets[me] + nme-1 < proc_new_offsets[endproc]) endproc--;
  while (proc_offsets[me] + nme-1 > proc_new_offsets[endproc+1]-1) endproc++;

  // loop over procs

  for (int iproc = startproc; iproc <= endproc; iproc++) {
    int istart = MAX(0,proc_new_offsets[iproc] - proc_offsets[me]);
    int iend = MIN(nme-1,proc_new_offsets[iproc+1]-1 - proc_offsets[me]);
    int nsend = iend - istart + 1;
    if (nsend == 0) continue;

    // send for this proc

    if (iproc != me) {
      MPI_Send(&buf[istart*size_one],nsend*size_one,MPI_DOUBLE,iproc,0,world);
    } else {

      // sending to self, copy buffers

      int offset_me = proc_offsets[me] - proc_new_offsets[me];
      memcpy(&buf_balance[(offset_me + istart)*size_one],&buf[istart*size_one],sizeof(double)*nsend*size_one);
    }
  }

  // wait for all recvs

  for (int n = 0; n < nswap; n++)
    MPI_Wait(&request[n],MPI_STATUS_IGNORE);

  nme = nme_balance;

  // swap buffers

  double *tmp = buf;
  buf = buf_balance;

  // cleanup

  memory->destroy(tmp);
  memory->destroy(proc_offsets);
  memory->destroy(proc_new_offsets);
  delete[] request;
}

/* ----------------------------------------------------------------------
   process params common to all dumps here
   if unknown param, call modify_param specific to the dump
------------------------------------------------------------------------- */

void Dump::modify_params(int narg, char **arg)
{
  if (narg == 0) utils::missing_cmd_args(FLERR, "dump_modify", error);

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"append") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify append", error);
      append_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;

    } else if (strcmp(arg[iarg],"balance") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify balance", error);
      if (nprocs > 1)
        balance_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;

    } else if (strcmp(arg[iarg],"buffer") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify buffer", error);
      buffer_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      if (buffer_flag && buffer_allow == 0)
        error->all(FLERR,"Dump_modify buffer yes not allowed for this style");
      iarg += 2;

    } else if (strcmp(arg[iarg],"colname") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify colname", error);
      if (strcmp(arg[iarg+1],"default") == 0) {
        for (auto &item : keyword_user) item.clear();
        iarg += 2;
      } else {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "dump_modify colname", error);
        int icol = -1;
        if (utils::is_integer(arg[iarg + 1])) {
          icol = utils::inumeric(FLERR,arg[iarg + 1],false,lmp);
          if (icol < 0) icol = keyword_user.size() + icol + 1;
          icol--;
        } else {
          try {
            icol = key2col.at(arg[iarg + 1]);
          } catch (std::out_of_range &) {
            icol = -1;
          }
        }
        if ((icol < 0) || (icol >= (int)keyword_user.size()))
          error->all(FLERR, "Incorrect dump_modify arguments: {} {} {}",
                     arg[iarg], arg[iarg+1], arg[iarg+2]);
        keyword_user[icol] = arg[iarg+2];
        iarg += 3;
      }

    } else if (strcmp(arg[iarg],"delay") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify delay", error);
      delaystep = utils::bnumeric(FLERR,arg[iarg+1],false,lmp);
      if (delaystep >= 0) delay_flag = 1;
      else delay_flag = 0;
      iarg += 2;

    } else if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify every", error);
      int idump;
      for (idump = 0; idump < output->ndump; idump++)
        if (strcmp(id,output->dump[idump]->id) == 0) break;
      int n;
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        delete[] output->var_dump[idump];
        output->var_dump[idump] = utils::strdup(&arg[iarg+1][2]);
        output->last_dump[idump] = -1;
        n = 0;
      } else {
        n = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
        if (n <= 0) error->all(FLERR, "Invalid dump_modify every argument: {}", n);
        output->next_dump[idump] = (update->ntimestep/n)*n+n;
      }
      output->mode_dump[idump] = 0;
      output->every_dump[idump] = n;
      iarg += 2;

    } else if (strcmp(arg[iarg],"every/time") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify every/time", error);
      int idump;
      for (idump = 0; idump < output->ndump; idump++)
        if (strcmp(id,output->dump[idump]->id) == 0) break;
      double delta;
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        delete[] output->var_dump[idump];
        output->var_dump[idump] = utils::strdup(&arg[iarg+1][2]);
        delta = 0.0;
      } else {
        delta = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (delta <= 0.0) error->all(FLERR, "Invalid dump_modify every/time argument: {}", delta);
      }
      output->mode_dump[idump] = 1;
      output->every_time_dump[idump] = delta;
      output->next_dump[idump] = update->ntimestep;
      iarg += 2;

    } else if (strcmp(arg[iarg],"fileper") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify fileper", error);
      if (!multiproc)
        error->all(FLERR,"Cannot use dump_modify fileper without % in dump file name");
      int nper = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nper <= 0) error->all(FLERR, "Invalid dump_modify fileper argument: {}", nper);

      multiproc = nprocs/nper;
      if (nprocs % nper) multiproc++;
      fileproc = me/nper * nper;
      int fileprocnext = MIN(fileproc+nper,nprocs);
      nclusterprocs = fileprocnext - fileproc;
      if (me == fileproc) filewriter = 1;
      else filewriter = 0;
      int icluster = fileproc/nper;

      MPI_Comm_free(&clustercomm);
      MPI_Comm_split(world,icluster,0,&clustercomm);

      delete[] multiname;
      char *ptr = strchr(filename,'%');
      *ptr = '\0';
      multiname = utils::strdup(fmt::format("{}{}{}", filename, icluster, ptr+1));
      *ptr = '%';
      iarg += 2;

    } else if (strcmp(arg[iarg],"first") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify first", error);
      first_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;

    } else if (strcmp(arg[iarg],"flush") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify flush", error);
      flush_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;

    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify format", error);

      if (strcmp(arg[iarg+1],"none") == 0) {
        delete[] format_line_user;
        delete[] format_int_user;
        delete[] format_bigint_user;
        delete[] format_float_user;
        format_line_user = nullptr;
        format_int_user = nullptr;
        format_bigint_user = nullptr;
        format_float_user = nullptr;
        // pass format none to child classes which may use it
        // not an error if they don't
        modify_param(narg-iarg,&arg[iarg]);
        iarg += 2;
        continue;
      }

      if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "dump_modify format", error);

      if (strcmp(arg[iarg+1],"line") == 0) {
        delete[] format_line_user;
        format_line_user = utils::strdup(arg[iarg+2]);
        iarg += 3;
      } else {   // pass other format options to child classes
        int n = modify_param(narg-iarg,&arg[iarg]);
        if (n == 0) error->all(FLERR,"Unknown dump_modify format keyword: {}", arg[iarg+1]);
        iarg += n;
      }

    } else if (strcmp(arg[iarg],"header") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify header", error);
      write_header_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;

    } else if (strcmp(arg[iarg],"maxfiles") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify maxfiles", error);
      if (!multifile)
        error->all(FLERR,"Cannot use dump_modify maxfiles without * in dump file name");
      // wipe out existing storage
      if (maxfiles > 0) {
        for (int idx=0; idx < numfiles; ++idx)
          delete[] nameslist[idx];
        delete[] nameslist;
      }
      maxfiles = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (maxfiles == 0) error->all(FLERR, "Invalid dump_modify maxfiles argument: {}", maxfiles);
      if (maxfiles > 0) {
        nameslist = new char*[maxfiles];
        numfiles = 0;
        for (int idx=0; idx < maxfiles; ++idx)
          nameslist[idx] = nullptr;
        fileidx = 0;
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"nfile") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify nfile", error);
      if (!multiproc)
        error->all(FLERR,"Cannot use dump_modify nfile without % in dump file name");
      int nfile = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nfile <= 0) error->all(FLERR, "Invalid dump_modify nfile argument: {}", nfile);
      nfile = MIN(nfile,nprocs);

      multiproc = nfile;
      int icluster = static_cast<int> ((bigint) me * nfile/nprocs);
      fileproc = static_cast<int> ((bigint) icluster * nprocs/nfile);
      int fcluster = static_cast<int> ((bigint) fileproc * nfile/nprocs);
      if (fcluster < icluster) fileproc++;
      int fileprocnext =
        static_cast<int> ((bigint) (icluster+1) * nprocs/nfile);
      fcluster = static_cast<int> ((bigint) fileprocnext * nfile/nprocs);
      if (fcluster < icluster+1) fileprocnext++;
      nclusterprocs = fileprocnext - fileproc;
      if (me == fileproc) filewriter = 1;
      else filewriter = 0;

      MPI_Comm_free(&clustercomm);
      MPI_Comm_split(world,icluster,0,&clustercomm);

      delete[] multiname;
      char *ptr = strchr(filename,'%');
      *ptr = '\0';
      multiname = utils::strdup(fmt::format("{}{}{}", filename, icluster, ptr+1));
      *ptr = '%';
      iarg += 2;

    } else if (strcmp(arg[iarg],"pad") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify pad", error);
      padflag = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (padflag < 0) error->all(FLERR, "Invalid dump_modify pad argument: {}", padflag);
      iarg += 2;

    } else if (strcmp(arg[iarg],"pbc") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify pbc", error);
      pbcflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;

    } else if (strcmp(arg[iarg],"skip") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      skipflag = 1;
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        delete[] skipvar;
        skipvar = utils::strdup(&arg[iarg+1][2]);
      } else error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"sort") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify sort", error);
      if (strcmp(arg[iarg+1],"off") == 0) sort_flag = 0;
      else if (strcmp(arg[iarg+1],"id") == 0) {
        sort_flag = 1;
        sortcol = 0;
        sortorder = ASCEND;
      } else {
        sort_flag = 1;
        sortcol = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
        sortorder = ASCEND;
        if (sortcol == 0) error->all(FLERR, "Invalid dump_modify sort argument: {}", sortcol);
        if (sortcol < 0) {
          sortorder = DESCEND;
          sortcol = -sortcol;
        }
        sortcolm1 = sortcol - 1;
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"time") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify time", error);
      time_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;

    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "dump_modify units", error);
      unit_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;

    } else {
      int n = modify_param(narg-iarg,&arg[iarg]);
      if (n == 0) error->all(FLERR,"Unknown dump_modify keyword: {}", arg[iarg]);
      iarg += n;
    }
  }
}

/* ---------------------------------------------------------------------- */

double Dump::compute_time()
{
  return update->atime + (update->ntimestep - update->atimestep)*update->dt;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

void Dump::pbc_allocate()
{
  memory->destroy(xpbc);
  memory->destroy(vpbc);
  memory->destroy(imagepbc);
  maxpbc = atom->nmax;
  memory->create(xpbc,maxpbc,3,"dump:xbpc");
  memory->create(vpbc,maxpbc,3,"dump:vbpc");
  memory->create(imagepbc,maxpbc,"dump:imagebpc");
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

double Dump::memory_usage()
{
  double bytes = memory->usage(buf,maxbuf);
  bytes += memory->usage(sbuf,maxsbuf);
  if (sort_flag) {
    if (sortcol == 0) bytes += memory->usage(ids,maxids);
    bytes += memory->usage(bufsort,size_one*maxsort);
    if (sortcol == 0) bytes += memory->usage(idsort,maxsort);
    bytes += memory->usage(index,maxsort);
    bytes += memory->usage(proclist,maxproc);
    if (irregular) bytes += (double)irregular->memory_usage();
  }
  if (pbcflag) {
    bytes += (double)6*maxpbc * sizeof(double);
    bytes += (double)maxpbc * sizeof(imageint);
  }
  return bytes;
}
