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

#include "restartable.h"
#include "buffer_reader.h"
#include "buffer_reader_root_file.h"
#include "file_writer.h"
#include "file_writer_wrapper.h"
#include "file_writer_sizer.h"
#include "file_writer_buffer.h"

#include "comm.h"
#include "error.h"

namespace LAMMPS_NS {

void Restartable::write_restart_global(FileWriter&& fw) { this->write_restart_global(fw); }
void Restartable::write_restart_global(FileWriter&) const {
  if (comm->me == 0)
    error->warning(FLERR, "write_restart_global on object without support");
}

void Restartable::write_restart_local(FileWriter&& fw) { this->write_restart_local(fw); }
void Restartable::write_restart_local(FileWriter&) const {
  if (comm->me == 0)
    error->warning(FLERR, "write_restart_local on object without support");
}

void Restartable::write_restart_settings(FileWriter&& fw) { this->write_restart_settings(fw); }
void Restartable::write_restart_settings(FileWriter&) const {
  if (comm->me == 0)
    error->warning(FLERR, "write_restart_settings on object without support");
}

void Restartable::read_restart_global(BufferReader&& br) { this->read_restart_global(br); }
void Restartable::read_restart_global(BufferReader&) {
  if (comm->me == 0)
    error->warning(FLERR, "read_restart_global on object without support");
}

void Restartable::read_restart_local(BufferReader&& br) { this->read_restart_local(br); }
void Restartable::read_restart_local(BufferReader&) {
  if (comm->me == 0)
    error->warning(FLERR, "read_restart_local on object without support");
}

void Restartable::read_restart_settings(BufferReader&& br) { this->read_restart_settings(br); }
void Restartable::read_restart_settings(BufferReader&) {
  if (comm->me == 0)
    error->warning(FLERR, "read_restart_settings on object without support");
}

// Simple forward helpers for temporaries -> references
void Restartable::read_restart_local_extra(int n, BufferReader&& br) { this->read_restart_local_extra(n, br); }
void Restartable::read_restart_local_extra(int n, BufferReader& br) { }

/* ---------------------------------------------------------------------- */

static void write_restart_impl(
  const Restartable *obj, FileWriter& fw, int nprocs,
  const std::vector<char>& ldata
) {
  if (obj->restartable_global) {
    obj->write_restart_global(fw);
  }
  if (obj->restartable_local) {
    fw.writev(nprocs);
    fw.writev((bigint)ldata.size());
    fw.writev(ldata.data(), ldata.size());
  }
}

/* ----------------------------------------------------------------------
   Default write_restart implementation for backwards compatibility. Writes all
     data on rank 0.
   For now, assume local_size is consistent on all ranks. If that needs to
     change, add a restartable_props bitfield variable with a variable-size bit.
------------------------------------------------------------------------- */

void Restartable::write_restart(FILE *fp) {
  if (!restartable() && comm->me == 0)
    error->warning(FLERR, "write_restart on object without support");

  std::vector<char> ldata;
  if (restartable_local) {
    FileWriterSizer sizer;
    this->write_restart_local(sizer);
    bigint local_size = sizer.size();

    std::vector<char> m_ldata(local_size);
    FileWriterBuffer fw{m_ldata.data(), local_size};
    this->write_restart_local(fw);

    if(comm->me == 0) ldata.resize(comm->nprocs * local_size);
    MPI_Gather(m_ldata.data(), local_size, MPI_BYTE,
               ldata.data(),   local_size, MPI_BYTE, 0, world);
  }

  if (comm->me == 0) {
    FileWriterWrapper fw{fp};

    if (write_restart_size_prefix) {
      FileWriterSizer sizer;
      write_restart_impl(this, sizer, comm->nprocs, ldata);
      fw.writev<int>(sizer.size());
    }

    write_restart_impl(this, fw, comm->nprocs, ldata);
  }
}

/* ----------------------------------------------------------------------
   Default read_restart implementation for backwards compatibility. Reads all
     data on rank 0 and sends it out.
   Ranks beyond checkpointed nprocs get empty buffers.
------------------------------------------------------------------------- */

void Restartable::read_restart(FILE *fp) {
  if (!restartable()) {
    if (comm->me == 0)
      error->warning(FLERR, "read_restart on object without support");
    return;
  }

  BufferReaderRootFile br(world, fp, error);
  this->read_restart(br);
}

void Restartable::read_restart(BufferReaderRootFile& br) {
  if (!restartable()) {
    // Read with file API, but keep reader consistent
    bigint prepos = 0, postpos = 0;
    if (comm->me == 0) prepos = platform::ftell(br.get_fp());
    this->read_restart(br.get_fp());
    if (comm->me == 0) postpos = platform::ftell(br.get_fp());

    bigint offset = postpos - prepos;
    MPI_Bcast(&offset, 1, MPI_LMP_BIGINT, 0, world);
    br.seek_offset(BRERR, offset);
    return;
  }

  // only important for read_restart(char*)
  if (write_restart_size_prefix) br.read<int>(BRERR);

  if (restartable_global) this->read_restart_global(br);
  if (restartable_local) {
    int nprocs = br.read<int>(BRERR);
    bigint nbytes = br.read<int>(BRERR);
    bigint nper = nbytes / nprocs;

    // Some manual work to avoid broadcasting all of the local data
    std::vector<char> bytes(comm->me == 0 ? nbytes : nper);
    br.read_raw_buf_root_only(BRERR, bytes.data(), nbytes);
    if (comm->nprocs <= nprocs) {
      MPI_Scatter(MPI_IN_PLACE, nper, MPI_CHAR,
                  bytes.data(), nper, MPI_CHAR, 0, world);
    } else {
      std::vector<int> counts(comm->nprocs), displs(comm->nprocs);
      counts[0] = nper;
      for (int i = 1; i < comm->nprocs; i++) {
        counts[i] = i < nprocs ? nper : 0;
        displs[i] = displs[i-1] + counts[i-1];
      }
      MPI_Scatterv(MPI_IN_PLACE, counts.data(), displs.data(), MPI_CHAR,
                   bytes.data(), counts[comm->me], MPI_CHAR, 0, world);
    }

    this->read_restart_local(BufferReader(bytes.data(), nper, error));

    if (comm->nprocs < nprocs) {
      if (comm->me == 0) {
        int n_extra_procs = nprocs - comm->nprocs;
        char* e_bytes = bytes.data() + (comm->nprocs * nper);
        BufferReader br(e_bytes, n_extra_procs*nper, error);
        this->read_restart_local_extra(n_extra_procs, br);
      } else {
        this->read_restart_local_extra(0, BufferReader(nullptr, 0, error));
      }
    }
  }
}

/* ----------------------------------------------------------------------
   As above, but for the Fix API.
------------------------------------------------------------------------- */

void Restartable::restart(char *data) {
  if (!restartable() && comm->me == 0)
    error->warning(FLERR, "read_restart on object without support");

  BufferReader br{data, error};
  if (restartable_global) this->read_restart_global(br);
  if (restartable_local) {
    int nprocs = br.read<int>(BRERR);
    bigint nbytes = br.read<bigint>(BRERR);
    bigint nper = nprocs / nbytes;

    // Reset position to 0, add max length info for safety
    br = br.sub_buf(BRERR, nbytes);

    if (comm->me < nprocs) {
      br.seek(BRERR, nper*comm->me);
      this->read_restart_local(br.sub_buf(BRERR, nper));
    } else {
      this->read_restart_local(BufferReader(nullptr, 0, error));
    }

    if (comm->nprocs < nprocs) {
      if (comm->me == 0) {
        br.seek(BRERR, nper*comm->nprocs);
        this->read_restart_local_extra(nprocs - comm->nprocs, br);
      } else {
        this->read_restart_local_extra(0, BufferReader(nullptr, 0, error));
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Read/write settings APIs, so classes that implement a separate settings
     function can let that be interoperable with the FILE* API.
   Expected to be called within a class's other read/write restart functions,
     not directly by checkpointers.
------------------------------------------------------------------------- */

void Restartable::write_restart_settings(FILE *fp) {
  FileWriterWrapper fw(fp);
  if (comm->me == 0) this->write_restart_settings(fw);
}

void Restartable::read_restart_settings(FILE *fp) {
  this->read_restart_settings(BufferReaderRootFile(world,fp,error));
}

} //namespace LAMMPS_NS
