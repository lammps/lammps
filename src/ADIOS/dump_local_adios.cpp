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

#include "dump_local_adios.h"

#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "memory.h"
#include "modify.h"
#include "universe.h"
#include "update.h"

#include <algorithm>       // std::shift_left
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

struct DumpLocalADIOSInternal {

  DumpLocalADIOSInternal() = default;

  ~DumpLocalADIOSInternal()
  {
    if (fh) fh.Close();
  }

  DumpLocalADIOSInternal(const DumpLocalADIOSInternal &) = delete;
  DumpLocalADIOSInternal &operator=(const DumpLocalADIOSInternal &) = delete;
  DumpLocalADIOSInternal(DumpLocalADIOSInternal &&) = delete;
  DumpLocalADIOSInternal &operator=(DumpLocalADIOSInternal &&) = delete;

  // Name of the ADIOS2 IO group; must match the XML config and the
  // DeclareIO() call below.
  static constexpr std::string_view ioName{"local"};

  std::unique_ptr<adios2::ADIOS> ad;    // owning handle – must outlive io/fh
  adios2::IO io;                        // non-owning reference handle
  adios2::Engine fh;                    // non-owning reference handle

  // Local-data table variables – exactly one is valid per dump instance.
  adios2::Variable<double> varLocal;    // used when adios_use_float_ == false
  adios2::Variable<float>  varLocalF;   // used when adios_use_float_ == true

  // Human-readable column labels captured at construction time.
  std::vector<std::string> columnNames;
};

}    // namespace LAMMPS_NS

// -------------------------------------------------------------------------
// Constructor / Destructor
// -------------------------------------------------------------------------

DumpLocalADIOS::DumpLocalADIOS(LAMMPS *lmp, int narg, char **arg)
    : DumpLocal(lmp, narg, arg), internal(std::make_unique<DumpLocalADIOSInternal>())
{
  // Create a default adios2_config.xml if it does not already exist.
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

  // Capture column labels for ADIOS attributes.
  // DumpLocal builds columns_default from the expanded field names in its
  // constructor, then may free earg; split columns_default to get the names.
  for (const auto &col : utils::split_words(columns_default))
    internal->columnNames.emplace_back(col);
}

// Defined out-of-line so that the compiler sees the complete definition of
// DumpLocalADIOSInternal when it instantiates ~unique_ptr<DumpLocalADIOSInternal>.
DumpLocalADIOS::~DumpLocalADIOS() = default;

// -------------------------------------------------------------------------
// modify_param()  – ADIOS-specific dump_modify keywords
// -------------------------------------------------------------------------

int DumpLocalADIOS::modify_param(int narg, char **arg)
{
  // precision fp32|float|fp64|double
  if (strcmp(arg[0], "precision") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "dump_modify precision", error);
    if (strcmp(arg[1], "fp32") == 0 || strcmp(arg[1], "float") == 0)
      adios_use_float_ = true;
    else if (strcmp(arg[1], "fp64") == 0 || strcmp(arg[1], "double") == 0)
      adios_use_float_ = false;
    else
      error->all(FLERR,
                 "Unknown precision value '{}': expected fp32/float or fp64/double",
                 arg[1]);
    return 2;
  }

  // shared yes|no
  if (strcmp(arg[0], "shared") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "dump_modify shared", error);
    if (strcmp(arg[1], "yes") == 0)
      adios_shared_replicas_ = true;
    else if (strcmp(arg[1], "no") == 0)
      adios_shared_replicas_ = false;
    else
      error->all(FLERR,
                 "Unknown shared value '{}': expected yes or no", arg[1]);
    return 2;
  }

  return DumpLocal::modify_param(narg, arg);
}

// -------------------------------------------------------------------------
// openfile()
// -------------------------------------------------------------------------

void DumpLocalADIOS::openfile()
{
  // Guard: DumpLocal::init_style() calls openfile() before our ADIOS IO is
  // set up.  Return early; write() will call openfile() again once IO is ready.
  if (!internal->io) return;

#if !defined(MPI_STUBS)
  const bool sharedMode = adios_shared_replicas_ && (universe->nworlds > 1);
  MPI_Comm writeComm    = sharedMode ? universe->uworld : world;
#endif

  if (multifile) {
    const auto filecurrent = utils::star_subst(filename, update->ntimestep, padflag);
#if defined(MPI_STUBS)
    internal->fh = internal->io.Open(filecurrent, adios2::Mode::Write);
#else
    internal->fh = internal->io.Open(filecurrent, adios2::Mode::Write, writeComm);
#endif
    if (!internal->fh) error->one(FLERR, "Cannot open dump file {}", filecurrent);
  } else {
    if (!singlefile_opened) {
#if defined(MPI_STUBS)
      internal->fh = internal->io.Open(filename, adios2::Mode::Write);
#else
      internal->fh = internal->io.Open(filename, adios2::Mode::Write, writeComm);
#endif
      if (!internal->fh) error->one(FLERR, "Cannot open dump file {}", filename);
      singlefile_opened = 1;
    }
  }
}

// -------------------------------------------------------------------------
// write()
// -------------------------------------------------------------------------

void DumpLocalADIOS::write()
{
  // Snapshot box bounds for metadata.
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

  // nmine = number of local entries contributed by this rank.
  nmine = count();

  // ntotal  = global count across all ranks in this partition.
  // localOffset = exclusive prefix-sum of nmine.
  const bigint bnmine = nmine;
  bigint ntotal_local{0};
  MPI_Allreduce(&bnmine, &ntotal_local, 1, MPI_LMP_BIGINT, MPI_SUM, world);

  bigint localOffset{0};
  MPI_Scan(&bnmine, &localOffset, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  localOffset -= bnmine;    // inclusive → exclusive

  // Per-replica suffix for shared mode.
  const bool sharedMode = adios_shared_replicas_ && (universe->nworlds > 1);
  const std::string sfx =
      sharedMode ? ("_replica" + std::to_string(universe->iworld)) : "";

  // Update the global shape and this rank's hyperslab selection.
  const auto nRowsGlobal = static_cast<std::size_t>(ntotal_local);
  const auto startRow    = static_cast<std::size_t>(localOffset);
  const auto nRowsLocal  = static_cast<std::size_t>(nmine);
  const auto nColumns    = static_cast<std::size_t>(size_one);

  if (adios_use_float_) {
    internal->varLocalF.SetShape({nRowsGlobal, nColumns});
    internal->varLocalF.SetSelection({{startRow, 0}, {nRowsLocal, nColumns}});
  } else {
    internal->varLocal.SetShape({nRowsGlobal, nColumns});
    internal->varLocal.SetSelection({{startRow, 0}, {nRowsLocal, nColumns}});
  }

  // Allocate / grow the pack buffer.
  if (nmine > maxbuf) {
    if (static_cast<bigint>(nmine) * size_one > MAXSMALLINT)
      error->all(FLERR, "Too much per-proc info for dump");
    maxbuf = nmine;
    memory->destroy(buf);
    memory->create(buf, maxbuf * size_one, "dump:buf");
  }

  pack(nullptr);

  openfile();
  internal->fh.BeginStep();

  // ------------------------------------------------------------------
  // Global scalars.
  // ------------------------------------------------------------------
  const bool isGlobalRoot = (me == 0) && (!sharedMode || universe->iworld == 0);
  if (isGlobalRoot) {
    internal->fh.Put<uint64_t>("ntimestep", static_cast<uint64_t>(update->ntimestep));
    internal->fh.Put<int>     ("nprocs",    nprocs);
  }

  // ------------------------------------------------------------------
  // Per-replica box bounds.
  // ------------------------------------------------------------------
  if (me == 0) {
    if (adios_use_float_) {
      internal->fh.Put<float>("boxxlo" + sfx, static_cast<float>(boxxlo));
      internal->fh.Put<float>("boxxhi" + sfx, static_cast<float>(boxxhi));
      internal->fh.Put<float>("boxylo" + sfx, static_cast<float>(boxylo));
      internal->fh.Put<float>("boxyhi" + sfx, static_cast<float>(boxyhi));
      internal->fh.Put<float>("boxzlo" + sfx, static_cast<float>(boxzlo));
      internal->fh.Put<float>("boxzhi" + sfx, static_cast<float>(boxzhi));
      if (domain->triclinic) {
        internal->fh.Put<float>("boxxy" + sfx, static_cast<float>(boxxy));
        internal->fh.Put<float>("boxxz" + sfx, static_cast<float>(boxxz));
        internal->fh.Put<float>("boxyz" + sfx, static_cast<float>(boxyz));
      }
    } else {
      internal->fh.Put<double>("boxxlo" + sfx, boxxlo);
      internal->fh.Put<double>("boxxhi" + sfx, boxxhi);
      internal->fh.Put<double>("boxylo" + sfx, boxylo);
      internal->fh.Put<double>("boxyhi" + sfx, boxyhi);
      internal->fh.Put<double>("boxzlo" + sfx, boxzlo);
      internal->fh.Put<double>("boxzhi" + sfx, boxzhi);
      if (domain->triclinic) {
        internal->fh.Put<double>("boxxy" + sfx, boxxy);
        internal->fh.Put<double>("boxxz" + sfx, boxxz);
        internal->fh.Put<double>("boxyz" + sfx, boxyz);
      }
    }
  }

  // LocalValueDim variables.
  internal->fh.Put<uint64_t>("nlocal"   + sfx, static_cast<uint64_t>(ntotal_local));
  internal->fh.Put<int>     ("ncolumns" + sfx, size_one);
  internal->fh.Put<uint64_t>("nme"      + sfx, static_cast<uint64_t>(bnmine));
  internal->fh.Put<uint64_t>("offset"   + sfx, static_cast<uint64_t>(localOffset));

  // ------------------------------------------------------------------
  // Local-data table.  See dump_custom_adios.cpp for the rationale
  // behind using the Variable object rather than a string name.
  // ------------------------------------------------------------------
  if (adios_use_float_) {
    const std::size_t nElems = nRowsLocal * nColumns;
    std::vector<float> fbuf(nElems);
    for (std::size_t i = 0; i < nElems; ++i) fbuf[i] = static_cast<float>(buf[i]);
    internal->fh.Put<float>(internal->varLocalF, fbuf.data());
  } else {
    internal->fh.Put<double>(internal->varLocal, buf);
  }

  internal->fh.EndStep();

  if (multifile) internal->fh.Close();
}

// -------------------------------------------------------------------------
// init_style()
// -------------------------------------------------------------------------

void DumpLocalADIOS::init_style()
{
  // Run the DumpLocal base initialisation (column assembly, format strings,
  // compute/fix validation, and boundary string).  DumpLocal::init_style()
  // calls openfile() at the end for single-file mode; our openfile() override
  // returns immediately when the ADIOS IO is not yet initialised, so the call
  // is harmless.  The file will be opened properly by the first write().
  DumpLocal::init_style();

  // Strip '%' from the filename: ADIOS2 always produces a single global
  // BP5 directory and does not support per-process file splitting.
  if (char *pct = std::strchr(filename, '%'); pct != nullptr) {
    const std::size_t tail = std::strlen(pct) + 1;
    std::shift_left(pct, pct + tail, 1);
  }

  // ------------------------------------------------------------------
  // Define the ADIOS2 IO group and all variables/attributes once.
  // ------------------------------------------------------------------
  if (!internal->io) {
#if !defined(MPI_STUBS)
    const bool sharedMode = adios_shared_replicas_ && (universe->nworlds > 1);
    if (sharedMode) {
      internal->ad.reset();
      internal->ad =
          std::make_unique<adios2::ADIOS>("adios2_config.xml", universe->uworld);
    }
#else
    const bool sharedMode = false;
#endif

    internal->io = internal->ad->DeclareIO(std::string{DumpLocalADIOSInternal::ioName});

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

    const std::string sfx =
        sharedMode ? ("_replica" + std::to_string(universe->iworld)) : "";

    // --- Global per-step scalars --------------------------------------
    if (!sharedMode || universe->iworld == 0) {
      internal->io.DefineVariable<uint64_t>("ntimestep");
      internal->io.DefineVariable<int>     ("nprocs");
    }

    // --- Per-replica box bounds --------------------------------------
    if (adios_use_float_) {
      internal->io.DefineVariable<float>("boxxlo" + sfx);
      internal->io.DefineVariable<float>("boxxhi" + sfx);
      internal->io.DefineVariable<float>("boxylo" + sfx);
      internal->io.DefineVariable<float>("boxyhi" + sfx);
      internal->io.DefineVariable<float>("boxzlo" + sfx);
      internal->io.DefineVariable<float>("boxzhi" + sfx);
      internal->io.DefineVariable<float>("boxxy"  + sfx);
      internal->io.DefineVariable<float>("boxxz"  + sfx);
      internal->io.DefineVariable<float>("boxyz"  + sfx);
    } else {
      internal->io.DefineVariable<double>("boxxlo" + sfx);
      internal->io.DefineVariable<double>("boxxhi" + sfx);
      internal->io.DefineVariable<double>("boxylo" + sfx);
      internal->io.DefineVariable<double>("boxyhi" + sfx);
      internal->io.DefineVariable<double>("boxzlo" + sfx);
      internal->io.DefineVariable<double>("boxzhi" + sfx);
      internal->io.DefineVariable<double>("boxxy"  + sfx);
      internal->io.DefineVariable<double>("boxxz"  + sfx);
      internal->io.DefineVariable<double>("boxyz"  + sfx);
    }

    // --- Attributes --------------------------------------------------
    internal->io.DefineAttribute<int>("triclinic", domain->triclinic);

    const auto *boundaryptr = reinterpret_cast<const int *>(domain->boundary);
    internal->io.DefineAttribute<int>("boundary", boundaryptr, 6);

    const auto nCols = static_cast<std::size_t>(size_one);
    internal->io.DefineAttribute<std::string>("columns",
                                              internal->columnNames.data(), nCols);
    internal->io.DefineAttribute<std::string>("columnstr",   columns);
    internal->io.DefineAttribute<std::string>("boundarystr", boundstr);
    internal->io.DefineAttribute<std::string>("LAMMPS/dump_style", "local");
    internal->io.DefineAttribute<std::string>("LAMMPS/version",    lmp->version);
    internal->io.DefineAttribute<std::string>("LAMMPS/num_ver",    std::to_string(lmp->num_ver));
    internal->io.DefineAttribute<std::string>("LAMMPS/precision",
                                              adios_use_float_ ? "float32" : "float64");

    if (sharedMode)
      internal->io.DefineAttribute<int>("nreplicas", universe->nworlds);

    // --- LocalValueDim variables -------------------------------------
    internal->io.DefineVariable<uint64_t>("nlocal"   + sfx, {adios2::LocalValueDim});
    internal->io.DefineVariable<int>     ("ncolumns" + sfx, {adios2::LocalValueDim});
    internal->io.DefineVariable<uint64_t>("nme"      + sfx, {adios2::LocalValueDim});
    internal->io.DefineVariable<uint64_t>("offset"   + sfx, {adios2::LocalValueDim});

    // --- Local-data table with unknown shape at definition time ------
    constexpr std::size_t kUnknown{1};
    if (adios_use_float_) {
      internal->varLocalF = internal->io.DefineVariable<float>(
          "local_data" + sfx,
          {kUnknown, nCols},
          {kUnknown, 0},
          {kUnknown, nCols});
    } else {
      internal->varLocal = internal->io.DefineVariable<double>(
          "local_data" + sfx,
          {kUnknown, nCols},
          {kUnknown, 0},
          {kUnknown, nCols});
    }
  }
}
