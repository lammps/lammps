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

#include "dump_atom_zstd.h"
#include "domain.h"
#include "error.h"
#include "update.h"
#include "force.h"

#include <fmt/format.h>
#include <cstring>

using namespace LAMMPS_NS;

DumpAtomZstd::DumpAtomZstd(LAMMPS *lmp, int narg, char **arg) :
  DumpAtom(lmp, narg, arg)
{
  cctx = nullptr;
  zstdFp = nullptr;
  fp = nullptr;
  out_buffer_size = ZSTD_CStreamOutSize();
  out_buffer = new char[out_buffer_size];

  checksum_flag = 1;
  compression_level = 0; // = default

  if (!compressed)
    error->all(FLERR,"Dump atom/zstd only writes compressed files");
}

/* ---------------------------------------------------------------------- */

DumpAtomZstd::~DumpAtomZstd()
{
  if(cctx && zstdFp) zstd_close();

  delete [] out_buffer;
  out_buffer = nullptr;
  out_buffer_size = 0;
}

/* ----------------------------------------------------------------------
   generic opening of a dump file
   ASCII or binary or zstdipped
   some derived classes override this function
------------------------------------------------------------------------- */

void DumpAtomZstd::openfile()
{
  // single file, already opened, so just return

  if (singlefile_opened) return;
  if (multifile == 0) singlefile_opened = 1;

  // if one file per timestep, replace '*' with current timestep

  char *filecurrent = filename;
  if (multiproc) filecurrent = multiname;

  if (multifile) {
    char *filestar = filecurrent;
    filecurrent = new char[strlen(filestar) + 16];
    char *ptr = strchr(filestar,'*');
    *ptr = '\0';
    if (padflag == 0)
      sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
              filestar,update->ntimestep,ptr+1);
    else {
      char bif[8],pad[16];
      strcpy(bif,BIGINT_FORMAT);
      sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
      sprintf(filecurrent,pad,filestar,update->ntimestep,ptr+1);
    }
    *ptr = '*';
    if (maxfiles > 0) {
      if (numfiles < maxfiles) {
        nameslist[numfiles] = new char[strlen(filecurrent)+1];
        strcpy(nameslist[numfiles],filecurrent);
        ++numfiles;
      } else {
        remove(nameslist[fileidx]);
        delete[] nameslist[fileidx];
        nameslist[fileidx] = new char[strlen(filecurrent)+1];
        strcpy(nameslist[fileidx],filecurrent);
        fileidx = (fileidx + 1) % maxfiles;
      }
    }
  }

  // each proc with filewriter = 1 opens a file

  if (filewriter) {
    if (append_flag) {
      zstdFp = fopen(filecurrent,"ab");
    } else {
      zstdFp = fopen(filecurrent,"wb");
    }

    if (zstdFp == nullptr) error->one(FLERR,"Cannot open dump file");

    cctx = ZSTD_createCCtx();
    ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, compression_level);
    ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, checksum_flag);

    if (cctx == nullptr) error->one(FLERR,"Cannot create Zstd context");
  } else zstdFp = nullptr;

  // delete string with timestep replaced

  if (multifile) delete [] filecurrent;
}

/* ---------------------------------------------------------------------- */

void DumpAtomZstd::write_header(bigint ndump)
{
  std::string header;

  if ((multiproc) || (!multiproc && me == 0)) {
    if (unit_flag && !unit_count) {
      ++unit_count;
      header = fmt::format("ITEM: UNITS\n%s\n",update->unit_style);
    }

    if (time_flag) {
      header += fmt::format("ITEM: TIME\n%.16g\n", compute_time());
    }

    header += fmt::format("ITEM: TIMESTEP\n{}\n", update->ntimestep);
    header += fmt::format("ITEM: NUMBER OF ATOMS\n{}\n", ndump);
    if (domain->triclinic == 0) {
      header += fmt::format("ITEM: BOX BOUNDS {}\n", boundstr);
      header += fmt::format("{0:-1.16e} {1:-1.16e}\n", boxxlo, boxxhi);
      header += fmt::format("{0:-1.16e} {1:-1.16e}\n", boxylo, boxyhi);
      header += fmt::format("{0:-1.16e} {1:-1.16e}\n", boxzlo, boxzhi);
    } else {
      header += fmt::format("ITEM: BOX BOUNDS xy xz yz {}\n", boundstr);
      header += fmt::format("{0:-1.16e} {1:-1.16e} {2:-1.16e}\n", boxxlo, boxxhi, boxxy);
      header += fmt::format("{0:-1.16e} {1:-1.16e} {2:-1.16e}\n", boxylo, boxyhi, boxxz);
      header += fmt::format("{0:-1.16e} {1:-1.16e} {2:-1.16e}\n", boxzlo, boxzhi, boxyz);
    }
    header += fmt::format("ITEM: ATOMS {}\n", columns);

    zstd_write(header.c_str(), header.length());
  }
}

/* ---------------------------------------------------------------------- */

void DumpAtomZstd::write_data(int n, double *mybuf)
{
  ZSTD_inBuffer input = { mybuf, (size_t)n, 0 };
  ZSTD_EndDirective mode = ZSTD_e_continue;

  do {
      ZSTD_outBuffer output = { out_buffer, out_buffer_size, 0 };
      size_t const remaining = ZSTD_compressStream2(cctx, &output, &input, mode);
      fwrite(out_buffer, sizeof(char), output.pos, zstdFp);
  } while(input.pos < input.size);
}

/* ---------------------------------------------------------------------- */

void DumpAtomZstd::write()
{
  DumpAtom::write();
  if (filewriter) {
    if (multifile) {
      zstd_close();
    } else {
      if (flush_flag && zstdFp) {
        zstd_flush();
        fflush(zstdFp);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int DumpAtomZstd::modify_param(int narg, char **arg)
{
  int consumed = DumpAtom::modify_param(narg, arg);
  if(consumed == 0) {
    if (strcmp(arg[0],"checksum") == 0) {
      if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[1],"yes") == 0) checksum_flag = 1;
      else if (strcmp(arg[1],"no") == 0) checksum_flag = 0;
      else error->all(FLERR,"Illegal dump_modify command");
      return 2;
    } else if (strcmp(arg[0],"compression_level") == 0) {
      if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
      compression_level = force->inumeric(FLERR,arg[1]);
      if (compression_level <= 0) error->all(FLERR,"Illegal dump_modify command");
      return 2;
    }
  }
  return consumed;
}

/* ---------------------------------------------------------------------- */

void DumpAtomZstd::zstd_write(const void * buffer, size_t length)
{
  ZSTD_inBuffer input = { buffer, length, 0 };
  ZSTD_EndDirective mode = ZSTD_e_continue;
  
  do {
    ZSTD_outBuffer output = { out_buffer, out_buffer_size, 0 };
    size_t const remaining = ZSTD_compressStream2(cctx, &output, &input, mode);
    fwrite(out_buffer, sizeof(char), output.pos, zstdFp);
  } while(input.pos < input.size);
}

void DumpAtomZstd::zstd_flush() {
  size_t remaining;
  ZSTD_inBuffer input = { nullptr, 0, 0 };
  ZSTD_EndDirective mode = ZSTD_e_flush;

  do {
    ZSTD_outBuffer output = { out_buffer, out_buffer_size, 0 };
    remaining = ZSTD_compressStream2(cctx, &output, &input, mode);
    fwrite(out_buffer, sizeof(char), output.pos, zstdFp);
  } while(remaining);
}

/* ---------------------------------------------------------------------- */

void DumpAtomZstd::zstd_close()
{
  size_t remaining;
  ZSTD_inBuffer input = { nullptr, 0, 0 };
  ZSTD_EndDirective mode = ZSTD_e_end;

  do {
    ZSTD_outBuffer output = { out_buffer, out_buffer_size, 0 };
    remaining = ZSTD_compressStream2(cctx, &output, &input, mode);
    fwrite(out_buffer, sizeof(char), output.pos, zstdFp);
  } while(remaining);

  ZSTD_freeCCtx(cctx);
  cctx = nullptr;
  if (zstdFp) fclose(zstdFp);
  zstdFp = nullptr;
}
