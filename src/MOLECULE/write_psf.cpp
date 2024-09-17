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

/* ----------------------------------------------------------------------
   Contributing author: Mitch Murphy (alphataubio at gmail.com)
------------------------------------------------------------------------- */

#include "write_psf.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "comm.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "improper.h"
#include "label_map.h"
#include "memory.h"
#include "modify.h"
#include "output.h"
#include "thermo.h"
#include "update.h"
#include "irregular.h"

#include <cstring>

#include <iostream>


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

WritePsf::WritePsf(LAMMPS *lmp) : Command(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  int flag,cols;
  int index_atom_iarray = atom->find_custom("psf",flag,cols);

  if( index_atom_iarray == -1 )
    atom_iarray_psf = nullptr;
  else
    atom_iarray_psf = atom->iarray[index_atom_iarray];
}

/* ----------------------------------------------------------------------
   called as write_psf command in input script
------------------------------------------------------------------------- */

void WritePsf::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"write_psf command before simulation box is defined");

  if (narg < 2) utils::missing_cmd_args(FLERR, "write_psf", error);

  igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Could not find group ID {}", arg[0]);
  groupbit = group->bitmask[igroup];

  // if filename contains a "*", replace with current timestep

  std::string file = arg[1];
  std::size_t found = file.find('*');
  if (found != std::string::npos)
    file.replace(found,1,fmt::format("{}",update->ntimestep));

  lmapflag = 1;
  // store current (default) setting since we may change it.
  int types_style = atom->types_style;
  int noinit = 0;

  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc
  // exception is when called by -r command-line switch
  //   then write_psf immediately follows reading of restart file
  //   assume that read_restart initialized necessary values
  //   if don't make exception:
  //     pair->init() can fail due to various unset values:
  //     e.g. pair hybrid coeffs, dpd ghost-atom velocity setting

  if (noinit == 0) {
    if (comm->me == 0) utils::logmesg(lmp,"System init for write_psf ...\n");
    lmp->init();

    // move atoms to new processors before writing file
    // do setup_pre_exchange to force update of per-atom info if needed
    // enforce PBC in case atoms are outside box
    // call borders() to rebuild atom map since exchange() destroys map

    modify->setup_pre_exchange();
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    domain->reset_box();
    comm->setup();
    comm->exchange();
    comm->borders();
    if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  }

  // set reorderflag = 1 if can simply reorder local atoms rather than sort
  // criteria: sorting by ID, atom IDs are consecutive from 1 to Natoms
  //           min/max IDs of group match size of group
  // compute ntotal_reorder, nme_reorder, idlo/idhi to test against later

  reorderflag = 0;

  if (nprocs > 1 && irregular == nullptr)
    irregular = new Irregular(lmp);


  write(file);
  // restore saved setting
  atom->types_style = types_style;
}

/* ---------------------------------------------------------------------- */

int WritePsf::count()
{
  if (igroup == 0) return atom->nlocal;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) m++;
  return m;
}


/* ----------------------------------------------------------------------
   called from command()
   might later let it be directly called within run/minimize loop
------------------------------------------------------------------------- */

void WritePsf::write(const std::string &file)
{
  // special case where reneighboring is not done in integrator
  //   on timestep data file is written (due to build_once being set)
  // if box is changing, must be reset, else data file will have
  //   wrong box size and atoms will be lost when data file is read
  // other calls to pbc and domain and comm are not made,
  //   b/c they only make sense if reneighboring is actually performed

  //if (neighbor->build_once) domain->reset_box();

  // natoms = sum of nlocal = value to write into data file
  // if unequal and thermo lostflag is "error", don't write data file

  bigint nblocal = atom->nlocal;
  bigint natoms;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms && output->thermo->lostflag == Thermo::ERROR)
    error->all(FLERR,"Atom count is inconsistent, cannot write data file");

  // sum up bond,angle,dihedral,improper counts
  // may be different than atom->nbonds,nangles, etc. if broken/turned-off

  // FIXME: handle Atom::TEMPLATE

  if (atom->molecular == Atom::MOLECULAR && (atom->nbonds || atom->nbondtypes)) {
    nbonds_local = atom->avec->pack_bond(nullptr);
    MPI_Allreduce(&nbonds_local,&nbonds,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }
  if (atom->molecular == Atom::MOLECULAR && (atom->nangles || atom->nangletypes)) {
    nangles_local = atom->avec->pack_angle(nullptr);
    MPI_Allreduce(&nangles_local,&nangles,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }

  if (atom->molecular == Atom::MOLECULAR && (atom->ndihedrals || atom->ndihedraltypes)) {
    ndihedrals_local = atom->avec->pack_dihedral(nullptr);
    MPI_Allreduce(&ndihedrals_local,&ndihedrals,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }

  if (atom->molecular == Atom::MOLECULAR && (atom->nimpropers || atom->nimpropertypes)) {
    nimpropers_local = atom->avec->pack_improper(nullptr);
    MPI_Allreduce(&nimpropers_local,&nimpropers,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }

  // open data file

  if (me == 0) {
    fp = fopen(file.c_str(),"w");
    if (fp == nullptr)
      error->one(FLERR,"Cannot open data file {}: {}",
                                   file, utils::getsyserror());
  }

  // proc 0 writes header, ntype-length arrays, force fields
  // label map must come before coeffs

  if (me == 0) {
    header();
  }

  // per atom info in Atoms and Velocities sections

  if (natoms) atoms();

  // molecular topology info if defined
  // do not write molecular topology for atom_style template

  if ( atom->molecular == Atom::MOLECULAR) {
    if (atom->nbonds && nbonds) bonds();
    if (atom->nangles && nangles) angles();
    if (atom->ndihedrals) dihedrals();
    if (atom->nimpropers) impropers();
  }

  // close data file

  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out data file header
------------------------------------------------------------------------- */

void WritePsf::header()
{

  fmt::print(fp,"PSF EXT XPLOR\n\n         5 !NTITLE\n");

  fmt::print(fp,"* LAMMPS psf file via write_psf, version {}, timestep = {}, units = {}\n",
             lmp->version, update->ntimestep, update->unit_style);

  fmt::print(fp,"* II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I)\n");
  fmt::print(fp,"* expanded format EXT:\n");
  fmt::print(fp,"* (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8) XPLOR\n");
  fmt::print(fp,"* [https://userguide.mdanalysis.org/stable/formats/reference/psf.html]\n");

}




/* ----------------------------------------------------------------------
   write out Bonds section of psf file
------------------------------------------------------------------------- */

void WritePsf::bonds()
{
  // communication buffer for all my Bond info
  // maxrow X ncol = largest buffer needed by any proc

  int ncol = 3;
  int sendrow = static_cast<int> (nbonds_local);
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  tagint **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_psf:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_psf:buf");

  // pack my bond data into buf

  atom->avec->pack_bond(buf);

  // write one chunk of info per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;
  int j = 0;

  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fmt::print(fp,"\n {:9} !NBOND: bonds\n",nbonds);

    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_LMP_TAGINT,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_LMP_TAGINT,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      // void AtomVec::write_bond(FILE *fp, int n, tagint **buf, int index)
      // atom->avec->write_bond(fp,recvrow,buf,index);
      index += recvrow;

      for (int i = 0; i < recvrow; i++) {
        fmt::print(fp, " {:9} {:9}", buf[i][1], buf[i][2]);
        j++;
        if( j % 4 == 0) // newline every 4 bonds
          fmt::print(fp, "\n");
      }

    }

    if( j % 4 != 0) // newline after last bond if row not full
      fmt::print(fp, "\n");


  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_LMP_TAGINT,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Angles section of data file
------------------------------------------------------------------------- */

void WritePsf::angles()
{
  // communication buffer for all my Angle info
  // maxrow X ncol = largest buffer needed by any proc

  int ncol = 4;
  int sendrow = static_cast<int> (nangles_local);
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  tagint **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_psf:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_psf:buf");

  // pack my angle data into buf

  atom->avec->pack_angle(buf);

  // write one chunk of info per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;
  int j = 0;

  int index = 1;
  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fmt::print(fp,"\n {:9} !NTHETA: angles\n",nangles);

    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_LMP_TAGINT,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_LMP_TAGINT,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      index += recvrow;

      for (int i = 0; i < recvrow; i++) {
        fmt::print(fp, " {:9} {:9} {:9}", buf[i][1], buf[i][2], buf[i][3]);
        j++;
        if( j % 3 == 0) // newline every 3 angles
          fmt::print(fp, "\n");
      }

    }

    if( j % 3 != 0) // newline after last angle if row not full
      fmt::print(fp, "\n");

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_LMP_TAGINT,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Dihedrals section of data file
------------------------------------------------------------------------- */

void WritePsf::dihedrals()
{
  // communication buffer for all my Dihedral info
  // maxrow X ncol = largest buffer needed by any proc

  int ncol = 5;
  int sendrow = static_cast<int> (ndihedrals_local);
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  tagint **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_psf:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_psf:buf");

  // pack my dihedral data into buf

  atom->avec->pack_dihedral(buf);

  // write one chunk of info per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;
  int j = 0;

  int index = 1;
  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fmt::print(fp,"\n {:9} !NPHI: dihedrals\n", ndihedrals);

    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_LMP_TAGINT,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_LMP_TAGINT,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      index += recvrow;

      for (int i = 0; i < recvrow; i++) {
        fmt::print(fp, " {:9} {:9} {:9} {:9}", buf[i][1], buf[i][2], buf[i][3], buf[i][4]);
        j++;
        if( j % 2 == 0) // newline every 2 dihedrals
          fmt::print(fp, "\n");
      }

    }

    if( j % 2 != 0) // newline after last dihedral if row not full
      fmt::print(fp, "\n");

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_LMP_TAGINT,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Impropers section of data file
------------------------------------------------------------------------- */

void WritePsf::impropers()
{
  // communication buffer for all my Improper info
  // maxrow X ncol = largest buffer needed by any proc

  int ncol = 5;
  int sendrow = static_cast<int> (nimpropers_local);
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  tagint **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_psf:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_psf:buf");

  // pack my improper data into buf

  atom->avec->pack_improper(buf);

  // write one chunk of info per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;
  int j = 0;

  int index = 1;
  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fmt::print(fp,"\n {:9} !NIMPHI: impropers\n", nimpropers);

    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_LMP_TAGINT,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_LMP_TAGINT,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      index += recvrow;

      for (int i = 0; i < recvrow; i++) {
        fmt::print(fp, " {:9} {:9} {:9} {:9}", buf[i][1], buf[i][2], buf[i][3], buf[i][4]);
        j++;
        if( j % 2 == 0) // newline every 2 dihedrals
          fmt::print(fp, "\n");
      }

    }

    if( j % 2 != 0) // newline after last dihedral if row not full
      fmt::print(fp, "\n");

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_LMP_TAGINT,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Atoms section of psf file
------------------------------------------------------------------------- */

void WritePsf::atoms()
{
  // communication buffer for all my Atom info
  // maxrow X ncol = largest buffer needed by any proc

  int ncol = atom->avec->size_data_atom + 3;
  size_one = atom->avec->size_data_atom;
  int sendrow = atom->nlocal;
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  double **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_psf:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_psf:buf");

  // pack my atom data into buf
  atom->avec->pack_data(buf);

  // write one chunk of atoms per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;

  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fmt::print(fp,"\n {:8} !NATOM\n",count());

    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
          //std::cerr << "ok 1a\n";
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_DOUBLE,iproc,0,world,&request);
           // std::cerr << "ok 1b\n";
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
           // std::cerr << "ok 1c\n";
        MPI_Wait(&request,&status);
           // std::cerr << "ok 1d\n";
        MPI_Get_count(&status,MPI_DOUBLE,&recvrow);
           // std::cerr << "ok 1e\n";
        recvrow /= ncol;
      } else recvrow = sendrow;

          //std::cerr << "ok 1f\n";


      for (int i = 0; i < atom->natoms; i++) {

        // skip if not in group
        if (igroup && !((atom->mask)[i] & groupbit))
            continue;

        // II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I)
        // (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8)

        std::string label;
        auto atom_tag = ubuf(buf[i][0]).i;
        int atom_type = ubuf(buf[i][2]).i;
        // utils::logmesg(lmp, "write_psf... atom_tag={} atom_type={}\n", atom_tag, atom_type);

        fmt::print(fp, "{:10} ", atom_tag );

        if( atom_iarray_psf == nullptr ) {

          // defaults when atom_iarray_psf doesnt exists:
          // - molecule id for Segment ID, molecule ID, Residue ID
          // - numerical type for atom name and atom type
          fmt::print(fp, "{0:<8} {0:<8} {0:<8} {1:<8} {1:<4} ", ubuf(buf[i][1]).i, atom_type );

        } else {

          // atom segment
          int segment_type = atom_iarray_psf[atom_tag-1][0];
          label = atom->lmap->label(segment_type, Atom::SEGMENT);
          fmt::print(fp, "{:<8} ", label );

          // molecule ID
          fmt::print(fp, "{:<8} ", ubuf(buf[i][1]).i );

          // atom residue
          int residue_type = atom_iarray_psf[atom_tag-1][1];
          label = atom->lmap->label(residue_type, Atom::RESIDUE);
          fmt::print(fp, "{:<8} ", label );

          // atom name
          int name_type = atom_iarray_psf[atom_tag-1][2];
          label = atom->lmap->label(name_type, Atom::NAME);
          fmt::print(fp, "{:<8} ", label );

          // atom type
          label = atom->lmap->label(atom_type, Atom::ATOM);
          fmt::print(fp, "{:<4} ", label );

        }

        // charge
        fmt::print(fp, "{:12.6F}      ", buf[i][3] );

        // mass
        fmt::print(fp, "{:8g}           0\n", atom->mass[atom_type] );
      }
    }
  } else {
    std::cerr << "ok 2a\n";
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    std::cerr << "ok 2b\n";
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
    std::cerr << "ok 2c\n";
  }

  memory->destroy(buf);

  std::cerr << "ok 3\n";

}
