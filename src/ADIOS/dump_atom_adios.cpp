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

#include "dump_atom_adios.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "memory.h"
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
// Defined in the .cpp so that adios2.h is not part of the public interface.
// All ADIOS2 handle types (IO, Engine, Variable) are reference-counted
// lightweight wrappers; they are NOT copyable in a meaningful way, so the
// struct is non-copyable and non-movable.
// -------------------------------------------------------------------------

namespace LAMMPS_NS {

struct DumpAtomADIOSInternal {

  DumpAtomADIOSInternal() = default;

  // Explicit destructor: Engine must be closed before the owning ADIOS
  // object (held by unique_ptr 'ad') is destroyed.  Member destruction
  // order is reverse-declaration order, so fh destructs before ad;
  // however the Engine wrapper does NOT call Close() on destruction —
  // we must do it here explicitly.
  ~DumpAtomADIOSInternal()
  {
    if (fh) fh.Close();
    // ad unique_ptr destructs automatically after the body, which is
    // correct because the engine is already closed above.
  }

  DumpAtomADIOSInternal(const DumpAtomADIOSInternal &) = delete;
  DumpAtomADIOSInternal &operator=(const DumpAtomADIOSInternal &) = delete;
  DumpAtomADIOSInternal(DumpAtomADIOSInternal &&) = delete;
  DumpAtomADIOSInternal &operator=(DumpAtomADIOSInternal &&) = delete;

  // Name of the ADIOS2 IO group; must match the XML config and the
  // DeclareIO() call below.
  static constexpr std::string_view ioName{"atom"};

  std::unique_ptr<adios2::ADIOS> ad;    // owning handle – must outlive io/fh
  adios2::IO io;                        // non-owning reference handle
  adios2::Engine fh;                    // non-owning reference handle
  adios2::Variable<double> varAtoms;    // non-owning reference handle
};

}    // namespace LAMMPS_NS

// -------------------------------------------------------------------------
// Constructor / Destructor
// -------------------------------------------------------------------------

DumpAtomADIOS::DumpAtomADIOS(LAMMPS *lmp, int narg, char **arg)
    : DumpAtom(lmp, narg, arg), internal(std::make_unique<DumpAtomADIOSInternal>())
{
  // Create a default adios2_config.xml if it does not already exist.
  //
  // In MPI runs every rank may reach this concurrently; since every rank
  // would write identical content to a small file the race is benign in
  // practice (last writer wins, content is the same).  A more rigorous
  // approach would gate the write on rank 0 and then MPI_Barrier, but
  // that would add a coupling dependency here that the original code
  // deliberately avoids.
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
}

// Defined out-of-line so that the compiler sees the complete definition of
// DumpAtomADIOSInternal when it instantiates ~unique_ptr<DumpAtomADIOSInternal>.
DumpAtomADIOS::~DumpAtomADIOS() = default;

// -------------------------------------------------------------------------
// openfile()
// -------------------------------------------------------------------------

void DumpAtomADIOS::openfile()
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

void DumpAtomADIOS::write()
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
  // This must be done every step because ntotal may vary (NPT, grand-canonical).
  const auto nAtomsGlobal = static_cast<std::size_t>(ntotal);
  const auto startRow     = static_cast<std::size_t>(atomOffset);
  const auto nAtomsLocal  = static_cast<std::size_t>(nme);
  const auto nColumns     = static_cast<std::size_t>(size_one);

  internal->varAtoms.SetShape({nAtomsGlobal, nColumns});
  internal->varAtoms.SetSelection({{startRow, 0}, {nAtomsLocal, nColumns}});

  // Allocate / grow the pack buffer.
  // ADIOS2 imposes no per-process data-size ceiling, so nme*size_one is
  // NOT constrained to int (unlike legacy MPI-gather paths).
  if (nme > maxbuf) {
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

  // LocalValueDim variables – every rank writes its own value so that
  // readers can reconstruct per-rank extents without MPI.
  internal->fh.Put<uint64_t>("natoms",   static_cast<uint64_t>(ntotal));
  internal->fh.Put<int>     ("ncolumns", size_one);
  internal->fh.Put<uint64_t>("nme",      static_cast<uint64_t>(bnme));
  internal->fh.Put<uint64_t>("offset",   static_cast<uint64_t>(atomOffset));

  // Write the atom table using the Variable object (not a string name).
  // Using the object is required because SetShape/SetSelection above
  // update the object's internal state; a string-name lookup would
  // resolve to the stale placeholder dimensions from DefineVariable.
  internal->fh.Put<double>(internal->varAtoms, buf);

  internal->fh.EndStep();    // deferred Puts are flushed here

  if (multifile) internal->fh.Close();
}

// -------------------------------------------------------------------------
// init_style()
// -------------------------------------------------------------------------

void DumpAtomADIOS::init_style()
{
  if (image_flag == 0)
    size_one = 5;
  else
    size_one = 8;

  domain->boundary_string(boundstr);

  // ADIOS2 always produces a single global BP5 directory regardless of
  // the '%' multi-processor placeholder.  Strip any '%' from the filename
  // so the path passed to io.Open() is clean.
  // std::shift_left (C++20) shifts the range [pct, end) left by one,
  // effectively deleting the '%' character in-place.
  if (char *pct = std::strchr(filename, '%'); pct != nullptr) {
    const std::size_t tail = std::strlen(pct) + 1;    // include '\0'
    std::shift_left(pct, pct + tail, 1);
  }

  // Build the column name list.
  std::vector<std::string> columnNames;

  // BUG FIX (was heap corruption):
  // The original code assigned string literals directly to 'columns', a
  // char* member that the Dump base class manages with new[]/delete[].
  // Calling delete[] on a string-literal address silently corrupts the
  // allocator's free-list metadata, which manifests later as:
  //   "malloc_consolidate(): invalid chunk size"  /
  //   "free(): corrupted unsorted chunks"
  // Fix: always use utils::strdup() to place a heap-allocated copy.
  delete[] columns;
  columns = nullptr;

  if (scale_flag == 0 && image_flag == 0) {
    columns     = utils::strdup("id type x y z");
    columnNames = {"id", "type", "x", "y", "z"};
  } else if (scale_flag == 0 && image_flag == 1) {
    columns     = utils::strdup("id type x y z ix iy iz");
    columnNames = {"id", "type", "x", "y", "z", "ix", "iy", "iz"};
  } else if (scale_flag == 1 && image_flag == 0) {
    columns     = utils::strdup("id type xs ys zs");
    columnNames = {"id", "type", "xs", "ys", "zs"};
  } else if (scale_flag == 1 && image_flag == 1) {
    columns     = utils::strdup("id type xs ys zs ix iy iz");
    columnNames = {"id", "type", "xs", "ys", "zs", "ix", "iy", "iz"};
  }

  // Apply user-supplied column label overrides.
  for (int icol = 0; icol < static_cast<int>(columnNames.size()); ++icol)
    if (keyword_user[icol].size()) columnNames[icol] = keyword_user[icol];

  // Select the pack function pointer based on scale/image/triclinic flags.
  if      (scale_flag == 1 && image_flag == 0 && domain->triclinic == 0)
    pack_choice = &DumpAtomADIOS::pack_scale_noimage;
  else if (scale_flag == 1 && image_flag == 1 && domain->triclinic == 0)
    pack_choice = &DumpAtomADIOS::pack_scale_image;
  else if (scale_flag == 1 && image_flag == 0 && domain->triclinic == 1)
    pack_choice = &DumpAtomADIOS::pack_scale_noimage_triclinic;
  else if (scale_flag == 1 && image_flag == 1 && domain->triclinic == 1)
    pack_choice = &DumpAtomADIOS::pack_scale_image_triclinic;
  else if (scale_flag == 0 && image_flag == 0)
    pack_choice = &DumpAtomADIOS::pack_noscale_noimage;
  else if (scale_flag == 0 && image_flag == 1)
    pack_choice = &DumpAtomADIOS::pack_noscale_image;

  // ------------------------------------------------------------------
  // Define the ADIOS2 IO group and all variables/attributes once.
  // The guard (operator bool on IO) is idempotent: subsequent calls to
  // init_style() (e.g. after a restart) skip re-definition.
  // ------------------------------------------------------------------
  if (!internal->io) {
    internal->io = internal->ad->DeclareIO(std::string{DumpAtomADIOSInternal::ioName});

    if (!internal->io.InConfigFile()) {
      // Not configured by the user's XML: apply defaults.
      // BP5 is the recommended engine since ADIOS2 2.9 (lower memory,
      // better append performance, node-level shared-memory aggregation).
      internal->io.SetEngine("BP5");

      const int  num_aggregators = (multiproc > 0) ? multiproc : 1;
      const auto nstreams        = std::to_string(num_aggregators);
      // BP5 uses NumAggregators; the old BP4 "substreams" key is ignored.
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
    internal->io.DefineAttribute<int>("scaled",    scale_flag);
    internal->io.DefineAttribute<int>("image",     image_flag);

    // domain->boundary is int[3][2] – 6 values laid out contiguously.
    // reinterpret_cast to int* is well-defined here (same underlying type).
    const auto *boundaryptr = reinterpret_cast<const int *>(domain->boundary);
    internal->io.DefineAttribute<int>("boundary", boundaryptr, 6);

    const auto nCols = static_cast<std::size_t>(size_one);
    internal->io.DefineAttribute<std::string>("columns",     columnNames.data(), nCols);
    internal->io.DefineAttribute<std::string>("columnstr",   columns);
    internal->io.DefineAttribute<std::string>("boundarystr", boundstr);
    internal->io.DefineAttribute<std::string>("LAMMPS/dump_style", "atom");
    internal->io.DefineAttribute<std::string>("LAMMPS/version",    lmp->version);
    internal->io.DefineAttribute<std::string>("LAMMPS/num_ver",    std::to_string(lmp->num_ver));

    // LocalValueDim variables let readers reconstruct the per-rank layout.
    internal->io.DefineVariable<uint64_t>("nme",    {adios2::LocalValueDim});
    internal->io.DefineVariable<uint64_t>("offset", {adios2::LocalValueDim});

    // Atom table: the global shape and local selection are unknown at
    // definition time; they are set correctly in write() each step via
    // SetShape() / SetSelection() on the Variable object.
    constexpr std::size_t kUnknown{1};
    internal->varAtoms = internal->io.DefineVariable<double>(
        "atoms",
        {kUnknown, nCols},    // shape   (updated each step)
        {kUnknown, 0},        // start   (updated each step)
        {kUnknown, nCols});   // count   (updated each step)
  }
}
