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
   Contributing author:          Norbert Podhorszki (ORNL)
   ADIOS 2.11.0 (BP5) and C++20: Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "dump_custom_adios.h"

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <algorithm>       // std::shift_left
#include <cmath>
#include <cstddef>         // std::size_t
#include <cstring>         // std::strchr, std::strlen
#include <filesystem>
#include <fstream>
#include <memory>          // std::make_unique, std::unique_ptr
#include <string>
#include <string_view>
#include <vector>

#include "adios2.h"
#include "adios_common.h"

using namespace LAMMPS_NS;

// -------------------------------------------------------------------------
// Pimpl implementation type
// -------------------------------------------------------------------------

namespace LAMMPS_NS {

struct DumpCustomADIOSInternal {

  DumpCustomADIOSInternal() = default;

  // Engine must be closed before the owning ADIOS object is destroyed.
  // See DumpAtomADIOSInternal for the full rationale.
  ~DumpCustomADIOSInternal()
  {
    if (fh) fh.Close();
  }

  DumpCustomADIOSInternal(const DumpCustomADIOSInternal &) = delete;
  DumpCustomADIOSInternal &operator=(const DumpCustomADIOSInternal &) = delete;
  DumpCustomADIOSInternal(DumpCustomADIOSInternal &&) = delete;
  DumpCustomADIOSInternal &operator=(DumpCustomADIOSInternal &&) = delete;

  // Name of the ADIOS2 IO group; must match the XML config and the
  // DeclareIO() call below.
  static constexpr std::string_view ioName{"custom"};

  std::unique_ptr<adios2::ADIOS> ad;    // owning handle – must outlive io/fh
  adios2::IO io;                        // non-owning reference handle
  adios2::Engine fh;                    // non-owning reference handle
  adios2::Variable<double> varAtoms;    // non-owning reference handle

  // Human-readable column labels captured at construction time.
  std::vector<std::string> columnNames;
};

}    // namespace LAMMPS_NS

// -------------------------------------------------------------------------
// Constructor / Destructor
// -------------------------------------------------------------------------

DumpCustomADIOS::DumpCustomADIOS(LAMMPS *lmp, int narg, char **arg)
    : DumpCustom(lmp, narg, arg), internal(std::make_unique<DumpCustomADIOSInternal>())
{
  // Create a default adios2_config.xml if it does not already exist.
  // See DumpAtomADIOS constructor for discussion of the MPI-concurrent race.
  namespace fs = std::filesystem;
  if (!fs::exists("adios2_config.xml")) {
    if (std::ofstream cfg{"adios2_config.xml"}) { cfg << default_config; }
  }

  try {
#if defined(MPI_STUBS)
    internal->ad = std::make_unique<adios2::ADIOS>("adios2_config.xml");
#else
    internal->ad = std::make_unique<adios2::ADIOS>("adios2_config.xml", world);
#endif
  } catch (const std::ios_base::failure &e) {
    error->all(FLERR, "ADIOS initialization failed with error: {}", e.what());
  }

  // Capture the user-requested field names from the dump command arguments.
  // These become the column labels written as ADIOS2 attributes.
  internal->columnNames.reserve(nfield);
  for (int i = 0; i < nfield; ++i) internal->columnNames.emplace_back(earg[i]);
}

// Defined out-of-line so that the compiler sees the complete definition of
// DumpCustomADIOSInternal when it instantiates ~unique_ptr<DumpCustomADIOSInternal>.
DumpCustomADIOS::~DumpCustomADIOS() = default;

// -------------------------------------------------------------------------
// openfile()
// -------------------------------------------------------------------------

void DumpCustomADIOS::openfile()
{
  if (multifile) {
    // One file (directory) per timestep: substitute '*' with the step number.
    const auto filecurrent = utils::star_subst(filename, update->ntimestep, padflag);
#if defined(MPI_STUBS)
    internal->fh = internal->io.Open(filecurrent, adios2::Mode::Write);
#else
    internal->fh = internal->io.Open(filecurrent, adios2::Mode::Write, world);
#endif
    if (!internal->fh) error->one(FLERR, "Cannot open dump file {}", filecurrent);
  } else {
    if (!singlefile_opened) {
#if defined(MPI_STUBS)
      internal->fh = internal->io.Open(filename, adios2::Mode::Write);
#else
      internal->fh = internal->io.Open(filename, adios2::Mode::Write, world);
#endif
      if (!internal->fh) error->one(FLERR, "Cannot open dump file {}", filename);
      singlefile_opened = 1;
    }
  }
}

// -------------------------------------------------------------------------
// write()
// -------------------------------------------------------------------------

void DumpCustomADIOS::write()
{
  // Snapshot box bounds.
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
    boxxy  = domain->xy;
    boxxz  = domain->xz;
    boxyz  = domain->yz;
  }

  // nme = atom count on this rank
  nme = count();

  // ntotal  = global atom count across all ranks
  // atomOffset = exclusive prefix sum of nme (atoms written by ranks 0..me-1)
  const bigint bnme = nme;
  MPI_Allreduce(&bnme, &ntotal, 1, MPI_LMP_BIGINT, MPI_SUM, world);

  bigint atomOffset{0};
  MPI_Scan(&bnme, &atomOffset, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  atomOffset -= bnme;    // MPI_Scan gives inclusive sum; we need exclusive

  // Update the global shape and this rank's hyperslab selection.
  const auto nAtomsGlobal = static_cast<std::size_t>(ntotal);
  const auto startRow     = static_cast<std::size_t>(atomOffset);
  const auto nAtomsLocal  = static_cast<std::size_t>(nme);
  const auto nColumns     = static_cast<std::size_t>(size_one);

  internal->varAtoms.SetShape({nAtomsGlobal, nColumns});
  internal->varAtoms.SetSelection({{startRow, 0}, {nAtomsLocal, nColumns}});

  // Allocate / grow the pack buffer.
  if (nme > maxbuf) {
    if (static_cast<bigint>(nme) * size_one > MAXSMALLINT)
      error->all(FLERR, "Too much per-proc info for dump");
    maxbuf = nme;
    memory->destroy(buf);
    memory->create(buf, maxbuf * size_one, "dump:buf");
  }
  if (sort_flag && sortcol == 0 && nme > maxids) {
    maxids = nme;
    memory->destroy(ids);
    memory->create(ids, maxids, "dump:ids");
  }

  if (sort_flag && sortcol == 0)
    pack(ids);
  else
    pack(nullptr);

  if (sort_flag) sort();

  openfile();
  internal->fh.BeginStep();

  // ------------------------------------------------------------------
  // Scalars – written by rank 0 only (global metadata).
  // ------------------------------------------------------------------
  if (me == 0) {
    internal->fh.Put<uint64_t>("ntimestep", static_cast<uint64_t>(update->ntimestep));
    internal->fh.Put<int>     ("nprocs",    nprocs);

    internal->fh.Put<double>("boxxlo", boxxlo);
    internal->fh.Put<double>("boxxhi", boxxhi);
    internal->fh.Put<double>("boxylo", boxylo);
    internal->fh.Put<double>("boxyhi", boxyhi);
    internal->fh.Put<double>("boxzlo", boxzlo);
    internal->fh.Put<double>("boxzhi", boxzhi);

    if (domain->triclinic) {
      internal->fh.Put<double>("boxxy", boxxy);
      internal->fh.Put<double>("boxxz", boxxz);
      internal->fh.Put<double>("boxyz", boxyz);
    }
  }

  // LocalValueDim variables – every rank writes its own value.
  internal->fh.Put<uint64_t>("natoms",   static_cast<uint64_t>(ntotal));
  internal->fh.Put<int>     ("ncolumns", size_one);
  internal->fh.Put<uint64_t>("nme",      static_cast<uint64_t>(bnme));
  internal->fh.Put<uint64_t>("offset",   static_cast<uint64_t>(atomOffset));

  // BUG FIX: use the Variable object, NOT the string name "atoms".
  //
  // SetShape()/SetSelection() update the internal state of the Variable
  // object 'varAtoms'.  When Put is given a string name, ADIOS2 resolves
  // it to the variable's original placeholder dimensions from
  // DefineVariable({1, nColumns}, ...) — completely ignoring the updated
  // shape.  The engine then calculates the wrong buffer size, reads past
  // the end of 'buf', and corrupts the heap.
  //
  //   WRONG:  internal->fh.Put<double>("atoms", buf);
  //   RIGHT:
  internal->fh.Put<double>(internal->varAtoms, buf);

  internal->fh.EndStep();    // deferred Puts are flushed here

  if (multifile) internal->fh.Close();
}

// -------------------------------------------------------------------------
// init_style()
// -------------------------------------------------------------------------

void DumpCustomADIOS::init_style()
{
  // Assemble the merged column string from defaults and user overrides.
  delete[] columns;
  std::string combined;
  int icol = 0;
  for (const auto &item : utils::split_words(columns_default)) {
    if (!combined.empty()) combined += ' ';
    combined += keyword_user[icol].size() ? keyword_user[icol] : item;
    ++icol;
  }
  columns = utils::strdup(combined);

  domain->boundary_string(boundstr);

  // Strip '%' from the filename (same rationale as dump_atom_adios).
  if (char *pct = std::strchr(filename, '%'); pct != nullptr) {
    const std::size_t tail = std::strlen(pct) + 1;
    std::shift_left(pct, pct + tail, 1);
  }

  // -----------------------------------------------------------------------
  // Validate compute / fix / variable references and region.
  //
  // These loops do two things:
  //   1. Cache the pointers that count() and pack() will use at each step.
  //   2. Emit an early error if an ID is missing or a fix frequency is
  //      incompatible, rather than crashing silently during the first write.
  // They are NOT dead code even though their results are not used further
  // in this function.
  // -----------------------------------------------------------------------
  for (int i = 0; i < ncompute; ++i) {
    compute[i] = modify->get_compute_by_id(id_compute[i]);
    if (!compute[i])
      error->all(FLERR, "Could not find dump custom/adios compute ID {}", id_compute[i]);
  }

  for (int i = 0; i < nfix; ++i) {
    fix[i] = modify->get_fix_by_id(id_fix[i]);
    if (!fix[i])
      error->all(FLERR, "Could not find dump custom/adios fix ID {}", id_fix[i]);
    if (nevery % fix[i]->peratom_freq)
      error->all(FLERR, Error::NOLASTLINE,
                 "dump custom/adios and fix {} with ID {} not "
                 "computed at compatible times{}",
                 fix[i]->style, id_fix[i], utils::errorurl(7));
  }

  for (int i = 0; i < nvariable; ++i) {
    const int ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0)
      error->all(FLERR, "Could not find dump custom/adios variable name {}", id_variable[i]);
    variable[i] = ivariable;
  }

  if (idregion && !domain->get_region_by_id(idregion))
    error->all(FLERR, "Region {} for dump custom/adios does not exist", idregion);

  // ------------------------------------------------------------------
  // Define the ADIOS2 IO group and all variables/attributes once.
  // ------------------------------------------------------------------
  if (!internal->io) {
    internal->io = internal->ad->DeclareIO(std::string{DumpCustomADIOSInternal::ioName});

    if (!internal->io.InConfigFile()) {
      internal->io.SetEngine("BP5");

      const int  num_aggregators = (multiproc > 0) ? multiproc : 1;
      const auto nstreams        = std::to_string(num_aggregators);
      internal->io.SetParameters({{"NumAggregators", nstreams}});

      if (me == 0)
        utils::logmesg(lmp,
                       "ADIOS method for {} is n-to-m "
                       "(BP5 aggregation with {} writers)\n",
                       filename, nstreams);
    }

    // --- Per-step scalar variables -----------------------------------
    internal->io.DefineVariable<uint64_t>("ntimestep");
    internal->io.DefineVariable<uint64_t>("natoms");
    internal->io.DefineVariable<int>     ("nprocs");
    internal->io.DefineVariable<int>     ("ncolumns");

    internal->io.DefineVariable<double>("boxxlo");
    internal->io.DefineVariable<double>("boxxhi");
    internal->io.DefineVariable<double>("boxylo");
    internal->io.DefineVariable<double>("boxyhi");
    internal->io.DefineVariable<double>("boxzlo");
    internal->io.DefineVariable<double>("boxzhi");
    internal->io.DefineVariable<double>("boxxy");
    internal->io.DefineVariable<double>("boxxz");
    internal->io.DefineVariable<double>("boxyz");

    // --- Attributes (written once to the BP5 metadata) ---------------
    internal->io.DefineAttribute<int>("triclinic", domain->triclinic);

    // domain->boundary is int[3][2] – 6 values laid out contiguously.
    const auto *boundaryptr = reinterpret_cast<const int *>(domain->boundary);
    internal->io.DefineAttribute<int>("boundary", boundaryptr, 6);

    const auto nCols = static_cast<std::size_t>(size_one);
    internal->io.DefineAttribute<std::string>("columns",
                                              internal->columnNames.data(), nCols);
    internal->io.DefineAttribute<std::string>("columnstr",   columns);
    internal->io.DefineAttribute<std::string>("boundarystr", boundstr);
    internal->io.DefineAttribute<std::string>("LAMMPS/dump_style", "custom");
    internal->io.DefineAttribute<std::string>("LAMMPS/version",    lmp->version);
    internal->io.DefineAttribute<std::string>("LAMMPS/num_ver",    std::to_string(lmp->num_ver));

    // LocalValueDim variables.
    internal->io.DefineVariable<uint64_t>("nme",    {adios2::LocalValueDim});
    internal->io.DefineVariable<uint64_t>("offset", {adios2::LocalValueDim});

    // Atom table with unknown shape at definition time.
    constexpr std::size_t kUnknown{1};
    internal->varAtoms = internal->io.DefineVariable<double>(
        "atoms",
        {kUnknown, nCols},
        {kUnknown, 0},
        {kUnknown, nCols});
  }
}
