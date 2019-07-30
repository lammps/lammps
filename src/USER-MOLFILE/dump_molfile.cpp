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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "dump_molfile.h"
#include <mpi.h>
#include <cstring>
#include <cmath>
#include "domain.h"
#include "atom.h"
#include "comm.h"
#include "update.h"
#include "group.h"
#include "memory.h"
#include "error.h"

#include "molfile_interface.h"

using namespace LAMMPS_NS;
typedef MolfileInterface MFI;

// syntax:
// dump <id> <groupid> molfile <every> <filename> <type> [<path>]
// path defaults to "." -> will look for .so files in CWD.
//
// XXX: potential change: add more options and make them optional
// path <path>
// template <file> <type> (import name and topology information from file)
// bonds <yes|no>         (write out bond information)
// topology <yes|no>      (write out all topology information)

/* ---------------------------------------------------------------------- */

DumpMolfile::DumpMolfile(LAMMPS *lmp, int narg, char **arg)
  : Dump(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal dump molfile command");

  if (binary || compressed || multiproc)
    error->all(FLERR,"Invalid dump molfile filename");

  // required settings

  sort_flag = 1;
  sortcol = 0;

  // storage for collected information

  size_one = 4;
  if (atom->molecule_flag)  ++size_one;
  if (atom->q_flag)         ++size_one;
  if (atom->rmass_flag)     ++size_one;
  if (atom->radius_flag)    ++size_one;

  need_structure = 0;
  unwrap_flag = 0;
  velocity_flag = 0;
  topology_flag = 0;
  ntotal = 0;
  me = comm->me;

  coords = vels = masses = charges = radiuses = NULL;
  types = molids = NULL;
  ntypes = atom->ntypes;
  typenames = NULL;

  // allocate global array for atom coords

  bigint n = group->count(igroup);
  if (n > static_cast<bigint>(MAXSMALLINT/3/sizeof(float)))
    error->all(FLERR,"Too many atoms for dump molfile");
  if (n < 1)
    error->all(FLERR,"Not enough atoms for dump molfile");
  natoms = static_cast<int>(n);

  if (me == 0) {
    memory->create(types,natoms,"dump:types");
    memory->create(coords,3*natoms,"dump:coords");
    if (atom->molecule_flag) memory->create(molids,natoms,"dump:molids");
    if (atom->q_flag) memory->create(charges,natoms,"dump:charges");
    if (atom->rmass_flag) memory->create(masses,natoms,"dump:masses");
    if (atom->radius_flag) memory->create(radiuses,natoms,"dump:radiuses");

    mf = new MolfileInterface(arg[5],MFI::M_WRITE);

    const char *path = (const char *) ".";
    if (narg > 6)
      path=arg[6];

    if (mf->find_plugin(path)!= MFI::E_MATCH)
      error->one(FLERR,"No suitable molfile plugin found");

    if (screen)
      fprintf(screen,"Dump '%s' uses molfile plugin: %s\n",
              id, mf->get_plugin_name());
    if (logfile)
      fprintf(logfile,"Dump '%s' uses molfile plugin: %s\n",
              id,mf->get_plugin_name());
  }
}

/* ---------------------------------------------------------------------- */

DumpMolfile::~DumpMolfile()
{
  if (me == 0) {
    mf->close();
    memory->destroy(types);
    memory->destroy(coords);
    memory->destroy(vels);
    memory->destroy(masses);
    memory->destroy(charges);
    memory->destroy(radiuses);
    delete mf;
  }

  if (typenames) {
    for (int i = 1; i <= ntypes; i++)
      delete [] typenames[i];

    delete [] typenames;
    typenames = NULL;
  }
}

/* ---------------------------------------------------------------------- */

void DumpMolfile::init_style()
{
  if (sort_flag == 0 || sortcol != 0)
    error->all(FLERR,"Dump molfile requires sorting by atom ID");

  if (me == 0) {

    /* initialize typenames array to numeric types by default */
    if (typenames == NULL) {
      typenames = new char*[ntypes+1];
      for (int itype = 1; itype <= ntypes; itype++) {
        /* a 32-bit int can be maximally 10 digits plus sign */
        typenames[itype] = new char[12];
        sprintf(typenames[itype],"%d",itype);
      }
    }

    // open single file, one time only
    if (multifile == 0) openfile();
  }
}

/* ---------------------------------------------------------------------- */

void DumpMolfile::write()
{
  // simulation box dimensions

  if (domain->triclinic == 1) {
    double *h = domain->h;
    double alen = h[0];
    double blen = sqrt(h[5]*h[5] + h[1]*h[1]);
    double clen = sqrt(h[4]*h[4] + h[3]*h[3] + h[2]*h[2]);
    cell[0] = alen;
    cell[1] = blen;
    cell[2] = clen;
    cell[3] = (90.0 - asin((h[5]*h[4] + h[1]*h[3]) / blen/clen)); // alpha
    cell[4] = (90.0 - asin((h[0]*h[4]) / alen/clen));             // beta
    cell[5] = (90.0 - asin((h[0]*h[5]) / alen/blen));             // gamma
  } else {
    cell[0] = domain->xprd;
    cell[1] = domain->yprd;
    cell[2] = domain->zprd;
    cell[3] = cell[4] = cell[5] = 90.0f;
  }

  // nme = # of dump lines this proc will contribute to dump

  nme = count();
  bigint bnme = nme;

  // ntotal = total # of dump lines
  // nmax = max # of dump lines on any proc

  int nmax;
  MPI_Allreduce(&bnme,&ntotal,1,MPI_LMP_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);

  // for single file output, the number of atoms must not change.
  if (natoms != ntotal) {
    if (multifile == 0) {
      error->all(FLERR,"Single file molfile dump needs constant #atoms");
    } else {
      natoms = ntotal;
    }
  }
  ntotal = 0;

  // if file per timestep, open new file

  if (multifile) openfile();

  // insure proc 0 can receive everyone's info
  // limit nmax*size_one to int since used as arg in MPI_Rsend() below
  // pack my data into buf
  // if sorting on IDs also request ID list from pack()
  // sort buf as needed

  if (nmax > maxbuf) {
    if ((bigint) nmax * size_one > MAXSMALLINT)
      error->all(FLERR,"Too much per-proc info for dump");
    maxbuf = nmax;
    memory->destroy(buf);
    memory->create(buf,maxbuf*size_one,"dump:buf");
  }
  if (nmax > maxids) {
    maxids = nmax;
    memory->destroy(ids);
    memory->create(ids,maxids,"dump:ids");
  }

  pack(ids);
  sort();

  int tmp,nlines;

  if (me == 0) {
    MPI_Status status;
    MPI_Request request;
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(buf,maxbuf*size_one,MPI_DOUBLE,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&nlines);
        nlines /= size_one;
      } else nlines = nme;

      write_data(nlines,buf);
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,0,0,world);
  }
}

/* ---------------------------------------------------------------------- */

void DumpMolfile::openfile()
{
  // single file, already opened, so just return

  if (singlefile_opened) return;
  if (multifile == 0) singlefile_opened = 1;
  need_structure = 1;

  if (me == 0) {

    // close open file, if needed.
    if (mf->is_open()) mf->close();

    // if one file per timestep, replace '*' with current timestep

    char *filecurrent = new char[strlen(filename) + 16];
    if (multifile == 0) {
      strcpy(filecurrent,filename);
    } else {
      char *ptr = strchr(filename,'*');
      char *p1 = filename;
      char *p2 = filecurrent;
      while (p1 != ptr)
        *p2++ = *p1++;

      if (padflag == 0) {
        sprintf(p2,BIGINT_FORMAT "%s",update->ntimestep,ptr+1);
      } else {
        char bif[8],pad[16];
        strcpy(bif,BIGINT_FORMAT);
        sprintf(pad,"%%0%d%s%%s",padflag,&bif[1]);
        sprintf(p2,pad,update->ntimestep,ptr+1);
      }
    }

    if (mf->open(filecurrent,&natoms))
      error->one(FLERR,"Cannot open dump file");
    delete[] filecurrent;
  }
}

/* ---------------------------------------------------------------------- */

void DumpMolfile::pack(tagint *ids)
{
  int m,n;

  tagint *tag = atom->tag;
  int *type = atom->type;
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  m = n = 0;
  if (unwrap_flag) {
    double xprd = domain->xprd;
    double yprd = domain->yprd;
    double zprd = domain->zprd;
    double xy = domain->xy;
    double xz = domain->xz;
    double yz = domain->yz;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        int ix = (image[i] & IMGMASK) - IMGMAX;
        int iy = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        int iz = (image[i] >> IMG2BITS) - IMGMAX;

        buf[m++] = type[i];
        if (domain->triclinic) {
          buf[m++] = x[i][0] + ix * xprd + iy * xy + iz * xz;
          buf[m++] = x[i][1] + iy * yprd + iz * yz;
          buf[m++] = x[i][2] + iz * zprd;
        } else {
          buf[m++] = x[i][0] + ix * xprd;
          buf[m++] = x[i][1] + iy * yprd;
          buf[m++] = x[i][2] + iz * zprd;
        }
        if (atom->molecule_flag) buf[m++] = atom->molecule[i];
        if (atom->q_flag)        buf[m++] = atom->q[i];
        if (atom->rmass_flag)    buf[m++] = atom->mass[i];
        if (atom->radius_flag)   buf[m++] = atom->radius[i];
        ids[n++] = tag[i];
      }
    }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        buf[m++] = type[i];
        buf[m++] = x[i][0];
        buf[m++] = x[i][1];
        buf[m++] = x[i][2];
        if (atom->molecule_flag) buf[m++] = atom->molecule[i];
        if (atom->q_flag)        buf[m++] = atom->q[i];
        if (atom->rmass_flag)    buf[m++] = atom->mass[i];
        if (atom->radius_flag)   buf[m++] = atom->radius[i];
        ids[n++] = tag[i];
      }
  }
}

/* ---------------------------------------------------------------------- */

void DumpMolfile::write_data(int n, double *mybuf)
{
  if (me == 0) {
    // copy buf atom coords into global arrays
    int m = 0;
    for (int i = 0; i < n; i++) {
      types[ntotal] = static_cast<int>(mybuf[m++]);
      coords[3*ntotal + 0] = mybuf[m++];
      coords[3*ntotal + 1] = mybuf[m++];
      coords[3*ntotal + 2] = mybuf[m++];
      if (atom->molecule_flag) molids[ntotal]   = static_cast<int>(mybuf[m++]);
      if (atom->q_flag)        charges[ntotal]  = mybuf[m++];
      if (atom->rmass_flag)    masses[ntotal]   = mybuf[m++];
      if (atom->radius_flag)   radiuses[ntotal] = mybuf[m++];
      ++ntotal;
    }

    // if last chunk of atoms in this snapshot, write global arrays to file

    if (ntotal == natoms) {
      ntotal = 0;

      if (need_structure) {
        mf->property(MFI::P_NAME,types,typenames);
        mf->property(MFI::P_TYPE,types,typenames);

        if (atom->molecule_flag)
          mf->property(MFI::P_RESI,molids);

        if (atom->rmass_flag) {
          mf->property(MFI::P_MASS,masses);
        } else {
          mf->property(MFI::P_MASS,types,atom->mass);
        }

        if (atom->q_flag)
          mf->property(MFI::P_CHRG,charges);

        if (atom->radius_flag)
          mf->property(MFI::P_RADS,radiuses);

        // update/write structure information in plugin
        mf->structure();
        need_structure = 0;
      }
      double simtime = update->ntimestep * update->dt;
      mf->timestep(coords,NULL,cell,&simtime);
    }
  }
}

/* ---------------------------------------------------------------------- */

int DumpMolfile::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"unwrap") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"yes") == 0) unwrap_flag = 1;
    else if (strcmp(arg[1],"no") == 0) unwrap_flag = 0;
    else error->all(FLERR,"Illegal dump_modify command");
    return 2;

  } else if (strcmp(arg[0],"element") == 0) {
    if (narg < ntypes+1)
      error->all(FLERR, "Dump modify element names do not match atom types");

    if (typenames) {
      for (int i = 1; i <= ntypes; i++)
        delete [] typenames[i];

      delete [] typenames;
      typenames = NULL;
    }

    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      int n = strlen(arg[itype]) + 1;
      typenames[itype] = new char[n];
      strcpy(typenames[itype],arg[itype]);
    }

    return ntypes+1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory in buf and global coords array
------------------------------------------------------------------------- */

bigint DumpMolfile::memory_usage()
{
  bigint bytes = Dump::memory_usage();
  bytes += memory->usage(coords,natoms*3);
  bytes += sizeof(MFI);
  return bytes;
}
