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
#include "file_writer.h"
#include "file_writer_wrapper.h"
#include "file_writer_sizer.h"
#include "file_writer_buffer.h"

#include "comm.h"
#include "error.h"

namespace LAMMPS_NS {

void Restartable::write_restart_global(FileWriter *) const {
  if (comm->me == 0)
    error->warning(FLERR, "write_restart_global on object without support");
}

void Restartable::write_restart_local(FileWriter *) const {
  if (comm->me == 0)
    error->warning(FLERR, "write_restart_local on object without support");
}

void Restartable::write_restart_settings(FileWriter *) const {
  if (comm->me == 0)
    error->warning(FLERR, "write_restart_settings on object without support");
}

void Restartable::read_restart_global(BufferReader) {
  if (comm->me == 0)
    error->warning(FLERR, "read_restart_global on object without support");
}

void Restartable::read_restart_local(BufferReader) {
  if (comm->me == 0)
    error->warning(FLERR, "read_restart_local on object without support");
}

void Restartable::read_restart_settings(BufferReader *) {
  if (comm->me == 0)
    error->warning(FLERR, "read_restart_settings on object without support");
}

/* ---------------------------------------------------------------------- */

static void write_restart_impl(
  const Restartable *obj, FileWriter *fw, bigint global_size, int nprocs,
  const std::vector<char>& ldata
) {
  if (obj->restartable_global) {
    fw->writev(global_size);
    obj->write_restart_global(fw);
  }
  if (obj->restartable_local) {
    fw->writev(nprocs);
    fw->writev((bigint)ldata.size());
    fw->writev(ldata.data(), ldata.size());
  }
}

/* ----------------------------------------------------------------------
   Default write_restart implementation for backwards compatibility. Writes all
     data on rank 0.
   For now, assume local_size is consistent on all ranks. If that needs to
     change, add a restartable_props bitfield variable with a variable-size bit.
------------------------------------------------------------------------- */

void Restartable::write_restart(FILE *fp) {
  if (!restartable_local && !restartable_global && comm->me == 0)
    error->warning(FLERR, "write_restart on object without support");

  bigint global_size = 0;
  if (restartable_global) {
    FileWriterSizer sizer;
    this->write_restart_global(&sizer);
    global_size = sizer.size();
  }
  std::vector<char> ldata;
  if (restartable_local) {
    FileWriterSizer sizer;
    this->write_restart_local(&sizer);
    bigint local_size = sizer.size();

    std::vector<char> m_ldata(local_size);
    FileWriterBuffer fw{m_ldata.data(), local_size};
    this->write_restart_local(&fw);

    if(comm->me == 0) ldata.resize(comm->nprocs * local_size);
    MPI_Gather(m_ldata.data(), local_size, MPI_BYTE,
               ldata.data(),   local_size, MPI_BYTE, 0, world);
  }
  int total_size = 0;
  if (write_restart_size_prefix) {
    FileWriterSizer sizer;
    write_restart_impl(this, &sizer, global_size, comm->nprocs, ldata);
    int total_size = sizer.size();
  }

  // Only rank 0 does actual write
  if (comm->me == 0) {
    FileWriterWrapper fw{fp};
    if (write_restart_size_prefix) fw.writev(total_size);
    write_restart_impl(this, &fw, global_size, comm->nprocs, ldata);
  }
}

/* ----------------------------------------------------------------------
   Default read_restart implementation for backwards compatibility. Reads all
     data on rank 0 and sends it out.
   For now, assume data from procs beyond the current nprocs should be ignored.
     Otherwise, add bitfield options for how to distribute the rest. Similarly,
     assume ranks beyond checkpointed nprocs should just get empty buffers.
------------------------------------------------------------------------- */

void Restartable::read_restart(FILE *fp) {
  if (!restartable_local && !restartable_global && comm->me == 0)
    error->warning(FLERR, "read_restart on object without support");

  if (write_restart_size_prefix) {
    // write_restart_size_prefix only important for read_restart(char*)
    int size;
    if (comm->me == 0)
      utils::sfread(FLERR, &size, sizeof(size), 1, fp, nullptr, error);
  }

  if (restartable_global) {
    bigint nbytes;
    std::vector<char> bytes;
    if (comm->me == 0) {
      utils::sfread(FLERR, &nbytes, sizeof(nbytes), 1, fp, nullptr, error);
      bytes.resize(nbytes);
      utils::sfread(FLERR, bytes.data(), 1, nbytes, fp, nullptr, error);
    }

    MPI_Bcast(&nbytes, 1, MPI_LMP_BIGINT, 0, world);
    bytes.resize(nbytes);
    MPI_Bcast(bytes.data(), nbytes, MPI_CHAR, 0, world);

    this->read_restart_global(BufferReader(bytes.data(), bytes.size()));
  }
  if (restartable_local) {
    int nprocs = -1;
    std::vector<char> bytes;
    if (comm->me == 0) {
      utils::sfread(FLERR, &nprocs, sizeof(nprocs), 1, fp, nullptr, error);
      bigint nbytes;
      utils::sfread(FLERR, &nbytes, sizeof(nbytes), 1, fp, nullptr, error);
      bytes.resize(nbytes);
      utils::sfread(FLERR, bytes.data(), 1, nbytes, fp, nullptr, error);
    }
    bigint nper = bytes.size() / nprocs;

    MPI_Bcast(&nprocs, 1, MPI_INT, 0, world);
    MPI_Bcast(&nper, 1, MPI_LMP_BIGINT, 0, world);

    if (comm->nprocs <= nprocs) {
      if (comm->me != 0) bytes.resize(nper);
      MPI_Scatter(MPI_IN_PLACE, nper, MPI_CHAR,
                  bytes.data(), nper, MPI_CHAR, 0, world);
    } else {
      std::vector<int> counts(comm->nprocs), displs(comm->nprocs);
      counts[0] = nper;
      for (int i = 1; i < comm->nprocs; i++) {
        counts[i] = i < nprocs ? nper : 0;
        displs[i] = displs[i-1] + counts[i-1];
      }
      if (comm->me != 0) bytes.resize(counts[comm->me]);
      MPI_Scatterv(MPI_IN_PLACE, counts.data(), displs.data(), MPI_CHAR,
                   bytes.data(), counts[comm->me], MPI_CHAR, 0, world);
    }
    if (comm->me == 0) bytes.resize(nper);

    this->read_restart_local(BufferReader(bytes.data(), bytes.size()));
  }
}

/* ----------------------------------------------------------------------
   As above, but for the Fix API.
------------------------------------------------------------------------- */

void Restartable::restart(char *data) {
  if (!restartable_local && !restartable_global && comm->me == 0)
    error->warning(FLERR, "read_restart on object without support");

  BufferReader br{data};
  if (restartable_global) {
    bigint nbytes = br.read<bigint>();
    this->read_restart_global(br.sub_buf(nbytes));
  }
  if (restartable_local) {
    int nprocs = br.read<int>();
    bigint nbytes = br.read<bigint>();
    // This is effectively just adding the max length info
    br = br.sub_buf(nbytes);

    bigint nper = nprocs / nbytes;
    br.skip_bytes(nper * comm->me);
    this->read_restart_global(br.sub_buf(nper));
  }
}

/* ----------------------------------------------------------------------
   Read/write settings APIs, so classes that implement a separate settings
     function can let that be interoperable with the FILE* API.
   Expected to be called within a class's other read/write restart functions,
     not directly by checkpointers.
------------------------------------------------------------------------- */

void Restartable::write_restart_settings(FILE *fp) {
  if (comm->me == 0) {
    FileWriterSizer sizer;
    this->write_restart_settings(&sizer);

    FileWriterWrapper fw{fp};
    fw.writev(sizer.size());
    this->write_restart_settings(&fw);
  }
}

void Restartable::read_restart_settings(FILE *fp) {
  bigint nbytes;
  std::vector<char> bytes;
  if (comm->me == 0) {
      utils::sfread(FLERR, &nbytes, sizeof(nbytes), 1, fp, nullptr, error);
      bytes.resize(nbytes);
      utils::sfread(FLERR, bytes.data(), 1, nbytes, fp, nullptr, error);
  }

  MPI_Bcast(&nbytes, 1, MPI_LMP_BIGINT, 0, world);
  bytes.resize(nbytes);
  MPI_Bcast(bytes.data(), nbytes, MPI_CHAR, 0, world);
  BufferReader br{bytes.data(), nbytes};
  this->read_restart_settings(&br);
}

}
