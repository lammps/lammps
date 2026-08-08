// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   This file initially came from LIGGGHTS (www.liggghts.com)
   Copyright (2014) DCS Computing GmbH, Linz
   Copyright (2015) Johannes Kepler University Linz

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Daniel Queteschiner (DCS, JKU)
   Christoph Kloss (DCS)
   Richard Berger (JKU)
------------------------------------------------------------------------- */

#include "dump_vtk.h"

#include "arg_info.h"
#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "fix_store_atom.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "platform.h"
#include "region.h"
#include "safe_pointers.h"
#include "update.h"
#include "variable.h"

#include <cstdio>
#include <cstring>
#include <sstream>
#include <utility>
#include <vector>

using namespace LAMMPS_NS;

// customize by
// * adding an enum constant (add vector components in consecutive order)
// * adding a pack_*(int) function for the value
// * adjusting parse_vtk_fields function to add the pack_* function to pack_choice
//   (in case of vectors, adjust identify_vectors as well)
// * adjusting thresh part in modify_param and count functions

enum{X,Y,Z, // required for vtk, must come first
     ID,MOL,PROC,PROCP1,TYPE,ELEMENT,MASS,
     XS,YS,ZS,XSTRI,YSTRI,ZSTRI,XU,YU,ZU,XUTRI,YUTRI,ZUTRI,
     XSU,YSU,ZSU,XSUTRI,YSUTRI,ZSUTRI,
     IX,IY,IZ,
     VX,VY,VZ,FX,FY,FZ,
     Q,MUX,MUY,MUZ,MU,RADIUS,DIAMETER,
     OMEGAX,OMEGAY,OMEGAZ,ANGMOMX,ANGMOMY,ANGMOMZ,
     TQX,TQY,TQZ,
     COMPUTE,FIX,VARIABLE,IVEC,DVEC,IARRAY,DARRAY,
     ATTRIBUTES}; // must come last
enum{LT,LE,GT,GE,EQ,NEQ,XOR};
enum{VTK,VTP,VTU,PVTP,PVTU}; // file formats

// replace the '%' character in a file name with the given text

static std::string percent_subst(const std::string &name, const std::string &subst)
{
  auto pos = name.find('%');
  if (pos == std::string::npos) return name;
  return name.substr(0,pos) + subst + name.substr(pos+1);
}

// insert text in front of the file name extension

static std::string insert_before_extension(const std::string &name, const std::string &text)
{
  auto pos = name.find_last_of('.');
  if (pos == std::string::npos) return name + text;
  return name.substr(0,pos) + text + name.substr(pos);
}

/* ---------------------------------------------------------------------- */

DumpVTK::DumpVTK(LAMMPS *lmp, int narg, char **arg) :
    DumpCustom(lmp, narg, arg), label(nullptr), write_choice(nullptr), n_calls_(0),
    precision_warned(0), writeprec(VTKWriter::SINGLE), filecurrent(nullptr),
    domainfilecurrent(nullptr), parallelfilecurrent(nullptr), multiname_ex(nullptr)
{
  if (narg < 6) utils::missing_cmd_args(FLERR, "dump vtk", error);

  // process attributes
  // ioptional = start of additional optional args
  // only dump image and dump movie styles process optional args

  ioptional = parse_vtk_fields(nargnew,earg);

  if (ioptional < nargnew)
    error->all(FLERR,"Invalid attribute {} in dump vtk command", earg[ioptional]);
  size_one = pack_choice.size();
  current_pack_choice_key = -1;

  if (filewriter) reset_vtk_data_containers();

  vtk_file_format = VTK;

  char *suffix = filename + strlen(filename) - strlen(".vtp");
  if (suffix > filename && strcmp(suffix,".vtp") == 0) {
    if (multiproc) vtk_file_format = PVTP;
    else           vtk_file_format = VTP;
  } else if (suffix > filename && strcmp(suffix,".vtu") == 0) {
    if (multiproc) vtk_file_format = PVTU;
    else           vtk_file_format = VTU;
  }

  if (vtk_file_format == VTK) { // no multiproc support for legacy vtk format

    // the legacy format collects the entire snapshot in a single file.  the
    // '%' character is replaced by 0, as for any other single file writer in
    // LAMMPS, and only proc 0 writes that file

    if (multiproc) {
      MPI_Comm_free(&clustercomm);

      char *ptr = strchr(filename,'%');
      *ptr = '\0';
      char *singlename = utils::strdup(fmt::format("{}0{}", filename, ptr+1));
      delete[] filename;
      filename = singlename;

      if (me == 0)
        error->warning(FLERR,"Cannot write one file per processor with the legacy VTK file "
                       "format: writing all data to file {} instead", filename);
    }

    if (me != 0) filewriter = 0;
    fileproc = 0;
    multiproc = 0;
    nclusterprocs = nprocs;
  }

  // parallel vtp/vtu requires the proc number to be preceded by underscore '_'

  if (multiproc) multiname_ex = utils::strdup(percent_subst(filename, fmt::format("_{}", me)));
}

/* ---------------------------------------------------------------------- */

DumpVTK::~DumpVTK()
{
  delete[] filecurrent;
  delete[] domainfilecurrent;
  delete[] parallelfilecurrent;
  delete[] multiname_ex;
  delete[] label;
}

/* ---------------------------------------------------------------------- */

void DumpVTK::init_style()
{
  // default for element names = C

  if (typenames == nullptr) {
    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      typenames[itype] = new char[2];
      strcpy(typenames[itype],"C");
    }
  }

  // setup boundary string

  domain->boundary_string(boundstr);

  // setup function ptrs

  if (vtk_file_format == VTP || vtk_file_format == PVTP)
    write_choice = &DumpVTK::write_vtp;
  else if (vtk_file_format == VTU || vtk_file_format == PVTU)
    write_choice = &DumpVTK::write_vtu;
  else
    write_choice = &DumpVTK::write_vtk;

  // find current ptr for each compute,fix,variable and custom atom property
  // check that fix frequency is acceptable

  for (int i = 0; i < ncompute; i++) {
    compute[i] = modify->get_compute_by_id(id_compute[i]);
    if (!compute[i]) error->all(FLERR,"Could not find dump vtk compute ID {}", id_compute[i]);
  }

  for (int i = 0; i < nfix; i++) {
    fix[i] = modify->get_fix_by_id(id_fix[i]);
    if (!fix[i]) {
      error->all(FLERR,"Could not find dump vtk fix ID {}", id_fix[i]);
    } else {
      if (nevery % fix[i]->peratom_freq)
        error->all(FLERR,"Dump vtk and fix ID {} not called at compatible times{}", id_fix[i],
                   utils::errorurl(7));
    }
  }

  for (int i = 0; i < nvariable; i++) {
    int ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0)
      error->all(FLERR,"Could not find dump vtk variable name");
    variable[i] = ivariable;
  }

  int icustom,flag,cols;
  for (int i = 0; i < ncustom; i++) {
    icustom = atom->find_custom(id_custom[i],flag,cols);
    if (icustom < 0)
      error->all(FLERR,"Could not find custom per-atom property ID");
    custom[i] = icustom;
    if (!flag && !cols) custom_flag[i] = IVEC;
    else if (flag && !cols) custom_flag[i] = DVEC;
    else if (!flag && cols) custom_flag[i] = IARRAY;
    else if (flag && cols) custom_flag[i] = DARRAY;
  }

  // check validity of region

  if (idregion) {
    if (!domain->get_region_by_id(idregion))
      error->all(FLERR,"Region {} for dump vtk does not exist",idregion);
  }
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_header(bigint)
{
}

/* ---------------------------------------------------------------------- */

int DumpVTK::count()
{
  n_calls_ = 0;

  int i;

  // grow choose and variable vbuf arrays if needed

  const int nlocal = atom->nlocal;
  if (atom->nmax > maxlocal) {
    maxlocal = atom->nmax;

    memory->destroy(choose);
    memory->destroy(dchoose);
    memory->destroy(clist);
    memory->create(choose,maxlocal,"dump:choose");
    memory->create(dchoose,maxlocal,"dump:dchoose");
    memory->create(clist,maxlocal,"dump:clist");

    for (i = 0; i < nvariable; i++) {
      memory->destroy(vbuf[i]);
      memory->create(vbuf[i],maxlocal,"dump:vbuf");
    }
  }

  // invoke Computes for per-atom quantities
  // cannot invoke before first run, otherwise invoke if necessary

  if (ncompute) {
    for (i = 0; i < ncompute; i++) {
      if (!compute[i]->is_initialized())
        error->all(FLERR,"Dump compute ID {} cannot be invoked before initialization by a run",
          compute[i]->id);
      if (!(compute[i]->invoked_flag & Compute::INVOKED_PERATOM)) {
        compute[i]->compute_peratom();
        compute[i]->invoked_flag |= Compute::INVOKED_PERATOM;
      }
    }
  }

  // evaluate atom-style Variables for per-atom quantities

  if (nvariable)
    for (i = 0; i < nvariable; i++)
      input->variable->compute_atom(variable[i],igroup,vbuf[i],1,0);

  // choose all local atoms for output

  for (i = 0; i < nlocal; i++) choose[i] = 1;

  // un-choose if not in group

  if (igroup) {
    int *mask = atom->mask;
    for (i = 0; i < nlocal; i++)
      if (!(mask[i] & groupbit))
        choose[i] = 0;
  }

  // un-choose if not in region

  if (idregion) {
    auto *region = domain->get_region_by_id(idregion);
    if (region) {
      region->prematch();
      double **x = atom->x;
      for (i = 0; i < nlocal; i++)
        if (choose[i] && region->match(x[i][0],x[i][1],x[i][2]) == 0)
          choose[i] = 0;
    }
  }

  // un-choose if any threshold criterion isn't met

  if (nthresh) {
    double *ptr,*ptrhold;
    double *values;
    double value;
    int nstride,lastflag;

    for (int ithresh = 0; ithresh < nthresh; ithresh++) {

      // customize by adding to if statement

      if (thresh_array[ithresh] == ID) {
        tagint *tag = atom->tag;
        for (i = 0; i < nlocal; i++) dchoose[i] = tag[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == MOL) {
        if (!atom->molecule_flag)
          error->all(FLERR, "Threshold for an atom property that isn't allocated");
        tagint *molecule = atom->molecule;
        for (i = 0; i < nlocal; i++) dchoose[i] = molecule[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == PROC) {
        for (i = 0; i < nlocal; i++) dchoose[i] = me;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == PROCP1) {
        for (i = 0; i < nlocal; i++) dchoose[i] = me;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == TYPE) {
        int *type = atom->type;
        for (i = 0; i < nlocal; i++) dchoose[i] = type[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ELEMENT) {
        int *type = atom->type;
        for (i = 0; i < nlocal; i++) dchoose[i] = type[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == MASS) {
        if (atom->rmass) {
          ptr = atom->rmass;
          nstride = 1;
        } else {
          double *mass = atom->mass;
          int *type = atom->type;
          for (i = 0; i < nlocal; i++) dchoose[i] = mass[type[i]];
          ptr = dchoose;
          nstride = 1;
        }

      } else if (thresh_array[ithresh] == X) {
        ptr = &atom->x[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == Y) {
        ptr = &atom->x[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == Z) {
        ptr = &atom->x[0][2];
        nstride = 3;

      } else if (thresh_array[ithresh] == XS) {
        double **x = atom->x;
        double boxxlo = domain->boxlo[0];
        double invxprd = 1.0/domain->xprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (x[i][0] - boxxlo) * invxprd;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == YS) {
        double **x = atom->x;
        double boxylo = domain->boxlo[1];
        double invyprd = 1.0/domain->yprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (x[i][1] - boxylo) * invyprd;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ZS) {
        double **x = atom->x;
        double boxzlo = domain->boxlo[2];
        double invzprd = 1.0/domain->zprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (x[i][2] - boxzlo) * invzprd;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == XSTRI) {
        double **x = atom->x;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[0]*(x[i][0]-boxlo[0]) +
            h_inv[5]*(x[i][1]-boxlo[1]) + h_inv[4]*(x[i][2]-boxlo[2]);
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == YSTRI) {
        double **x = atom->x;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[1]*(x[i][1]-boxlo[1]) +
            h_inv[3]*(x[i][2]-boxlo[2]);
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ZSTRI) {
        double **x = atom->x;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[2]*(x[i][2]-boxlo[2]);
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == XU) {
        double **x = atom->x;
        imageint *image = atom->image;
        double xprd = domain->xprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = x[i][0] + ((image[i] & IMGMASK) - IMGMAX) * xprd;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == YU) {
        double **x = atom->x;
        imageint *image = atom->image;
        double yprd = domain->yprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = x[i][1] +
            ((image[i] >> IMGBITS & IMGMASK) - IMGMAX) * yprd;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ZU) {
        double **x = atom->x;
        imageint *image = atom->image;
        double zprd = domain->zprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = x[i][2] + ((image[i] >> IMG2BITS) - IMGMAX) * zprd;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == XUTRI) {
        double **x = atom->x;
        imageint *image = atom->image;
        double *h = domain->h;
        int xbox,ybox,zbox;
        for (i = 0; i < nlocal; i++) {
          xbox = (image[i] & IMGMASK) - IMGMAX;
          ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
          zbox = (image[i] >> IMG2BITS) - IMGMAX;
          dchoose[i] = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
        }
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == YUTRI) {
        double **x = atom->x;
        imageint *image = atom->image;
        double *h = domain->h;
        int ybox,zbox;
        for (i = 0; i < nlocal; i++) {
          ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
          zbox = (image[i] >> IMG2BITS) - IMGMAX;
          dchoose[i] = x[i][1] + h[1]*ybox + h[3]*zbox;
        }
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ZUTRI) {
        double **x = atom->x;
        imageint *image = atom->image;
        double *h = domain->h;
        int zbox;
        for (i = 0; i < nlocal; i++) {
          zbox = (image[i] >> IMG2BITS) - IMGMAX;
          dchoose[i] = x[i][2] + h[2]*zbox;
        }
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == XSU) {
        double **x = atom->x;
        imageint *image = atom->image;
        double boxxlo = domain->boxlo[0];
        double invxprd = 1.0/domain->xprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (x[i][0] - boxxlo) * invxprd +
            (image[i] & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == YSU) {
        double **x = atom->x;
        imageint *image = atom->image;
        double boxylo = domain->boxlo[1];
        double invyprd = 1.0/domain->yprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] =
            (x[i][1] - boxylo) * invyprd +
            (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == ZSU) {
        double **x = atom->x;
        imageint *image = atom->image;
        double boxzlo = domain->boxlo[2];
        double invzprd = 1.0/domain->zprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (x[i][2] - boxzlo) * invzprd +
            (image[i] >> IMG2BITS) - IMGMAX;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == XSUTRI) {
        double **x = atom->x;
        imageint *image = atom->image;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[0]*(x[i][0]-boxlo[0]) +
            h_inv[5]*(x[i][1]-boxlo[1]) +
            h_inv[4]*(x[i][2]-boxlo[2]) +
            (image[i] & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == YSUTRI) {
        double **x = atom->x;
        imageint *image = atom->image;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[1]*(x[i][1]-boxlo[1]) +
            h_inv[3]*(x[i][2]-boxlo[2]) +
            (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ZSUTRI) {
        double **x = atom->x;
        imageint *image = atom->image;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[2]*(x[i][2]-boxlo[2]) +
            (image[i] >> IMG2BITS) - IMGMAX;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == IX) {
        imageint *image = atom->image;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (image[i] & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == IY) {
        imageint *image = atom->image;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == IZ) {
        imageint *image = atom->image;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (image[i] >> IMG2BITS) - IMGMAX;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == VX) {
        ptr = &atom->v[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == VY) {
        ptr = &atom->v[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == VZ) {
        ptr = &atom->v[0][2];
        nstride = 3;
      } else if (thresh_array[ithresh] == FX) {
        ptr = &atom->f[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == FY) {
        ptr = &atom->f[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == FZ) {
        ptr = &atom->f[0][2];
        nstride = 3;

      } else if (thresh_array[ithresh] == Q) {
        if (!atom->q_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = atom->q;
        nstride = 1;
      } else if (thresh_array[ithresh] == MUX) {
        if (!atom->mu_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][0];
        nstride = 4;
      } else if (thresh_array[ithresh] == MUY) {
        if (!atom->mu_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][1];
        nstride = 4;
      } else if (thresh_array[ithresh] == MUZ) {
        if (!atom->mu_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][2];
        nstride = 4;
      } else if (thresh_array[ithresh] == MU) {
        if (!atom->mu_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][3];
        nstride = 4;

      } else if (thresh_array[ithresh] == RADIUS) {
        if (!atom->radius_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = atom->radius;
        nstride = 1;
      } else if (thresh_array[ithresh] == DIAMETER) {
        if (!atom->radius_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        double *radius = atom->radius;
        for (i = 0; i < nlocal; i++) dchoose[i] = 2.0*radius[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == OMEGAX) {
        if (!atom->omega_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = &atom->omega[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == OMEGAY) {
        if (!atom->omega_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = &atom->omega[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == OMEGAZ) {
        if (!atom->omega_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = &atom->omega[0][2];
        nstride = 3;
      } else if (thresh_array[ithresh] == ANGMOMX) {
        if (!atom->angmom_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = &atom->angmom[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == ANGMOMY) {
        if (!atom->angmom_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = &atom->angmom[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == ANGMOMZ) {
        if (!atom->angmom_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = &atom->angmom[0][2];
        nstride = 3;
      } else if (thresh_array[ithresh] == TQX) {
        if (!atom->torque_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = &atom->torque[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == TQY) {
        if (!atom->torque_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = &atom->torque[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == TQZ) {
        if (!atom->torque_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->torque[0][2];
        nstride = 3;

      } else if (thresh_array[ithresh] == COMPUTE) {
        i = ATTRIBUTES + nfield + ithresh;
        if (argindex[i] == 0) {
          ptr = compute[field2index[i]]->vector_atom;
          nstride = 1;
        } else {
          ptr = &compute[field2index[i]]->array_atom[0][argindex[i]-1];
          nstride = compute[field2index[i]]->size_peratom_cols;
        }

      } else if (thresh_array[ithresh] == FIX) {
        i = ATTRIBUTES + nfield + ithresh;
        if (argindex[i] == 0) {
          ptr = fix[field2index[i]]->vector_atom;
          nstride = 1;
        } else {
          ptr = &fix[field2index[i]]->array_atom[0][argindex[i]-1];
          nstride = fix[field2index[i]]->size_peratom_cols;
        }

      } else if (thresh_array[ithresh] == VARIABLE) {
        i = ATTRIBUTES + nfield + ithresh;
        ptr = vbuf[field2index[i]];
        nstride = 1;

      } else if (thresh_array[ithresh] == IVEC) {
        i = ATTRIBUTES + nfield + ithresh;
        int iwhich = custom[field2index[i]];
        int *ivector = atom->ivector[iwhich];
        for (i = 0; i < nlocal; i++)
          dchoose[i] = ivector[i];
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == DVEC) {
        i = ATTRIBUTES + nfield + ithresh;
        int iwhich = custom[field2index[i]];
        ptr = atom->dvector[iwhich];
        nstride = 1;

      } else if (thresh_array[ithresh] == IARRAY) {
        i = ATTRIBUTES + nfield + ithresh;
        int iwhich = custom[field2index[i]];
        int **iarray = atom->iarray[iwhich];
        int icol = argindex[i] - 1;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = iarray[i][icol];
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == DARRAY) {
        i = ATTRIBUTES + nfield + ithresh;
        int iwhich = custom[field2index[i]];
        double **darray = atom->darray[iwhich];
        ptr = &darray[0][argindex[i]-1];
        nstride = atom->dcols[iwhich];
      }

      // unselect atoms that don't meet threshold criterion
      // compare to single value or values stored in threshfix
      // copy ptr attribute into thresh_fix if this is first comparison

      if (thresh_last[ithresh] < 0) {
        lastflag = 0;
        value = thresh_value[ithresh];
      } else {
        lastflag = 1;
        int ilast = thresh_last[ithresh];
        values = thresh_fix[ilast]->vstore;
        ptrhold = ptr;
        if (thresh_first[ilast]) {
          thresh_first[ilast] = 0;
          for (i = 0; i < nlocal; i++, ptr += nstride) values[i] = *ptr;
          ptr = ptrhold;
        }
      }

      if (thresh_op[ithresh] == LT) {
        if (lastflag) {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr >= values[i]) choose[i] = 0;
        } else {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr >= value) choose[i] = 0;
        }
      } else if (thresh_op[ithresh] == LE) {
        if (lastflag) {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr > values[i]) choose[i] = 0;
        } else {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr > value) choose[i] = 0;
        }
      } else if (thresh_op[ithresh] == GT) {
        if (lastflag) {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr <= values[i]) choose[i] = 0;
        } else {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr <= value) choose[i] = 0;
        }
      } else if (thresh_op[ithresh] == GE) {
        if (lastflag) {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr < values[i]) choose[i] = 0;
        } else {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr < value) choose[i] = 0;
        }
      } else if (thresh_op[ithresh] == EQ) {
        if (lastflag) {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr != values[i]) choose[i] = 0;
        } else {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr != value) choose[i] = 0;
        }
      } else if (thresh_op[ithresh] == NEQ) {
        if (lastflag) {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr == values[i]) choose[i] = 0;
        } else {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if (choose[i] && *ptr == value) choose[i] = 0;
        }
      } else if (thresh_op[ithresh] == XOR) {
        if (lastflag) {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if ((choose[i] && *ptr == 0.0 && values[i] == 0.0) ||
                (*ptr != 0.0 && values[i] != 0.0))
              choose[i] = 0;
        } else {
          for (i = 0; i < nlocal; i++, ptr += nstride)
            if ((choose[i] && *ptr == 0.0 && value == 0.0) ||
                (*ptr != 0.0 && value != 0.0))
              choose[i] = 0;
        }
      }

      // update values stored in threshfix

      if (lastflag) {
        ptr = ptrhold;
        for (i = 0; i < nlocal; i++, ptr += nstride) values[i] = *ptr;
      }
    }
  }

  // compress choose flags into clist
  // nchoose = # of selected atoms
  // clist[i] = local index of each selected atom

  nchoose = 0;
  for (i = 0; i < nlocal; i++)
    if (choose[i]) clist[nchoose++] = i;

  return nchoose;
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write()
{
  // simulation box bounds

  if (domain->triclinic == 0) {
    boxxlo = domain->boxlo[0];
    boxxhi = domain->boxhi[0];
    boxylo = domain->boxlo[1];
    boxyhi = domain->boxhi[1];
    boxzlo = domain->boxlo[2];
    boxzhi = domain->boxhi[2];
  } else {
    domain->box_corners();
  }

  // nme = # of dump lines this proc contributes to dump

  nme = count();

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

  if (nmax > maxbuf) {
    if ((bigint) nmax * size_one > MAXSMALLINT)
      error->all(FLERR,"Too much per processor data info for dump");
    maxbuf = nmax;
    memory->destroy(buf);
    memory->create(buf,maxbuf*size_one,"dump:buf");
  }

  // ensure ids buffer is sized for sorting

  if (sort_flag && sortcol == 0 && nmax > maxids) {
    maxids = nmax;
    memory->destroy(ids);
    memory->create(ids,maxids,"dump:ids");
  }

  // pack my data into buf
  // if sorting on IDs also request ID list from pack()
  // sort buf as needed

  if (sort_flag && sortcol == 0) pack(ids);
  else pack(nullptr);
  if (sort_flag) sort();

  // filewriter = 1 = this proc writes to file
  //   ping each proc in my cluster, receive its data, write data to file
  // else wait for ping from fileproc, send my data to fileproc

  int tmp = 0;
  int nlines;
  MPI_Status status;
  MPI_Request request;

  // comm and output buf of doubles

  if (filewriter) {
    for (int iproc = 0; iproc < nclusterprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(buf,maxbuf*size_one,MPI_DOUBLE,me+iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&nlines);
        nlines /= size_one;
      } else nlines = nme;

      write_data(nlines,buf);
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,&status);
    MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,fileproc,0,world);
  }
}

/* ---------------------------------------------------------------------- */

void DumpVTK::pack(tagint *ids)
{
  int n = 0;
  for (auto &choice : pack_choice) {
      current_pack_choice_key = choice.first; // work-around for pack_vtk_compute, pack_vtk_fix, pack_vtk_variable
      (this->*(choice.second))(n);
      ++n;
  }

  if (ids) {
    tagint *tag = atom->tag;
    for (int i = 0; i < nchoose; i++)
      ids[i] = tag[clist[i]];
  }
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpVTK::setFileCurrent() {

  // per processor piece file.  "%" was replaced with "_<me>" in the
  // constructor; dump_modify fileper or nfile group procs into fewer files,
  // so their file id has to be recomputed for every snapshot

  if (multiproc > 1) {
    delete[] multiname_ex;
    int id;
    if (me + nclusterprocs == nprocs) // last filewriter
      id = multiproc - 1;
    else
      id = me/nclusterprocs;
    multiname_ex = utils::strdup(percent_subst(filename, fmt::format("_{}", id)));
  }

  std::string current = multiproc ? multiname_ex : filename;
  if (multifile) current = utils::star_subst(current, update->ntimestep, padflag);
  delete[] filecurrent;
  filecurrent = utils::strdup(current);

  // domain box data file: common for all procs, so without the "_<id>" part

  std::string domaincurrent = insert_before_extension(percent_subst(filename, ""), "_boundingBox");
  if (multifile) domaincurrent = utils::star_subst(domaincurrent, update->ntimestep, padflag);
  delete[] domainfilecurrent;
  domainfilecurrent = utils::strdup(domaincurrent);

  // parallel summary file: "%" removed and 'p' added to the file extension

  if (multiproc && me == 0) {
    std::string parallelcurrent = percent_subst(filename, "");
    parallelcurrent.insert(parallelcurrent.find_last_of('.') + 1, "p");
    if (multifile) parallelcurrent = utils::star_subst(parallelcurrent, update->ntimestep, padflag);
    delete[] parallelfilecurrent;
    parallelfilecurrent = utils::strdup(parallelcurrent);
  }
}

/* ---------------------------------------------------------------------- */

void DumpVTK::buf2arrays(int n, double *mybuf)
{
  points.reserve(points.size() + (std::size_t) 3*n);
  for (auto &array : myarrays) {
    if (array.type == Dump::STRING) array.strings.reserve(array.strings.size() + n);
    else array.values.reserve(array.values.size() + (std::size_t) n*array.ncomp);
  }

  for (int iatom=0; iatom < n; ++iatom) {
    const double *atom = &mybuf[iatom*size_one];

    points.push_back(atom[0]);
    points.push_back(atom[1]);
    points.push_back(atom[2]);

    int j=3; // 0,1,2 = x,y,z handled just above
    for (auto &array : myarrays) {
      if (array.type == Dump::STRING) {
        array.strings.emplace_back(typenames[static_cast<int>(atom[j])]);
        ++j;
      } else {
        for (int k = 0; k < array.ncomp; ++k) array.values.push_back(atom[j+k]);
        j += array.ncomp;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_domain(VTKWriter::Flavor flavor)
{
  try {
    VTKWriter writer(flavor, binary != 0, writeprec);
    writer.set_title(label ? label : "Generated by LAMMPS");

    if (domain->triclinic == 0) {
      writer.set_rectilinear_grid({boxxlo, boxxhi}, {boxylo, boxyhi}, {boxzlo, boxzhi});
    } else {

      // the corners of a LAMMPS box (from Domain::box_corners() called in
      // write()) and the corners of a VTK hexahedron cell run in a different
      // order around the bottom and top faces

      static const int corner_map[8] = {0, 1, 3, 2, 4, 5, 7, 6};
      double corners[8][3];
      for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 3; ++j) corners[i][j] = domain->corners[corner_map[i]][j];
      writer.set_hexahedron(corners);
    }
    writer.write(domainfilecurrent);
  } catch (VTKWriterException &e) {
    error->one(FLERR,"Cannot write dump vtk box file {}: {}", domainfilecurrent, e.what());
  }
}

/* ----------------------------------------------------------------------
   write the atom data of the current snapshot
------------------------------------------------------------------------- */

void DumpVTK::write_points(VTKWriter::Flavor flavor, bool unstructured)
{
  try {
    VTKWriter writer(flavor, binary != 0, writeprec);
    writer.set_title(label ? label : "Generated by LAMMPS");

    // the collected data is donated to the writer, since it is rebuilt for
    // every snapshot anyway.  only the array names and types remain valid,
    // which is all that write_pvtk() needs afterwards.

    if (unstructured) writer.set_unstructured_grid(std::move(points));
    else writer.set_polydata(std::move(points));

    for (auto &array : myarrays) {
      if (array.type == Dump::STRING) {
        writer.add_point_array(array.name, std::move(array.strings));
      } else if (array.type == Dump::INT) {
        writer.add_point_array(array.name, array.ncomp,
                               std::vector<int>(array.values.begin(), array.values.end()));
      } else {
        writer.add_point_array(array.name, array.ncomp, std::move(array.values));
      }
    }
    writer.write(filecurrent);
    check_coordinate_precision(writer.max_single_precision_value());
  } catch (VTKWriterException &e) {
    error->one(FLERR,"Cannot write dump vtk file {}: {}", filecurrent, e.what());
  }
}

/* ----------------------------------------------------------------------
   atom coordinates are stored in single precision by default, which is
   what visualization programs expect.  warn once when the system has
   grown so large that this no longer resolves the coordinates well
   enough.  with dump_modify double yes nothing is tracked, so the
   warning is skipped automatically.
------------------------------------------------------------------------- */

void DumpVTK::check_coordinate_precision(double maxcoord)
{
  if (precision_warned || (me != 0)) return;
  const double resolution = VTKWriter::single_precision_resolution(maxcoord, force->angstrom);
  if (resolution == 0.0) return;

  precision_warned = 1;
  error->warning(FLERR,"Dump vtk writes atom coordinates in single precision, which resolves "
                 "the largest coordinate of this system, {:.4g}, to only about {:.2g} length "
                 "units. Use dump_modify double yes if your analysis needs more resolution "
                 "than that.", maxcoord, resolution);
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_vtk(int n, double *mybuf)
{
  ++n_calls_;

  buf2arrays(n, mybuf);

  if (n_calls_ < nclusterprocs)
    return; // multiple processors but only proc 0 is a filewriter (-> nclusterprocs procs contribute to the filewriter's output data)

  setFileCurrent();

  write_points(VTKWriter::LEGACY, false);
  write_domain(VTKWriter::LEGACY);

  reset_vtk_data_containers();
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_vtp(int n, double *mybuf)
{
  write_xml_snapshot(n, mybuf, false);
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_vtu(int n, double *mybuf)
{
  write_xml_snapshot(n, mybuf, true);
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_xml_snapshot(int n, double *mybuf, bool unstructured)
{
  ++n_calls_;

  buf2arrays(n, mybuf);

  if (n_calls_ < nclusterprocs)
    return; // multiple processors but not all are filewriters (-> nclusterprocs procs contribute to the filewriter's output data)

  setFileCurrent();

  write_points(VTKWriter::XML, unstructured);

  if (me == 0) {
    if (multiproc) write_pvtk();

    // adjust the filename extension to the dataset type used for the box:
    // rectilinear grid (.vtr) for an orthogonal box, unstructured grid
    // (.vtu) for a triclinic one

    domainfilecurrent[strlen(domainfilecurrent)-1] = (domain->triclinic == 0) ? 'r' : 'u';
    write_domain(VTKWriter::XML);
  }

  reset_vtk_data_containers();
}

/* ----------------------------------------------------------------------
   construct the name of the per-processor piece file for piece "id" as it
   should be referenced from the parallel (.pvtp/.pvtu) summary file.
   mirrors the "%" -> "_<id>" substitution used to build the actual piece
   files (see constructor / setFileCurrent) and returns the name relative to
   the summary file (leading directory removed), since both live in the same
   directory.  filename is guaranteed to contain '%' when multiproc is set.
------------------------------------------------------------------------- */

std::string DumpVTK::pvtk_piece_filename(int id)
{
  std::string fname = percent_subst(filename, fmt::format("_{}", id));
  if (multifile) fname = utils::star_subst(fname, update->ntimestep, padflag);

  // strip the leading directory so the reference is relative to the summary file

  return platform::path_basename(fname);
}

/* ----------------------------------------------------------------------
   write the parallel summary file (.pvtp for PVTP, .pvtu for PVTU) by hand.
   we cannot let the VTK vtkXMLP*Writer classes do this: they derive the
   per-processor piece file names from the summary file name by stripping
   everything after the *first* '.', which yields wrong <Piece> references
   (and an extra bogus piece file) whenever the dump file name contains extra
   '.' characters (e.g. "dump.run.*%.vtp").  LAMMPS already writes the actual
   piece files itself with the correct names, so we emit a matching summary
   file directly.  called on rank 0 only, with multiproc set.
------------------------------------------------------------------------- */

void DumpVTK::write_pvtk()
{
  SafeFilePtr fp(fopen(parallelfilecurrent,"w"));
  if (!fp)
    error->one(FLERR,"Cannot open dump vtk parallel file {}: {}",
               parallelfilecurrent, utils::getsyserror());

  const char *gridtype = (vtk_file_format == PVTP) ? "PPolyData" : "PUnstructuredGrid";

  utils::print(fp, R"(<?xml version="1.0"?>)" "\n");
  utils::print(fp, R"(<VTKFile type="{}" version="1.0" byte_order="{}">)" "\n", gridtype,
               VTKWriter::xml_byte_order());
  utils::print(fp, R"(  <{} GhostLevel="0">)" "\n", gridtype);

  // declare the same arrays that reset_vtk_data_containers() set up for the
  // piece files, so that the two can never get out of step.  the declared
  // types have to match what write_points() makes the writer produce, so
  // floating point arrays follow the selected precision

  const char *real_type = VTKWriter::xml_real_type(writeprec);

  utils::print(fp, R"(    <PPointData>)" "\n");
  for (const auto &array : myarrays) {
    const char *type = real_type;
    if (array.type == Dump::INT) type = "Int32";
    else if (array.type == Dump::STRING) type = "String";
    if (array.ncomp == 3) {
      utils::print(fp, R"(      <PDataArray type="{}" Name="{}" NumberOfComponents="3"/>)" "\n",
                   type, array.name);
    } else {
      utils::print(fp, R"(      <PDataArray type="{}" Name="{}"/>)" "\n", type, array.name);
    }
  }
  utils::print(fp, R"(    </PPointData>)" "\n");

  utils::print(fp, R"(    <PPoints>)" "\n");
  utils::print(fp, R"(      <PDataArray type="{}" Name="Points" NumberOfComponents="3"/>)" "\n",
               real_type);
  utils::print(fp, R"(    </PPoints>)" "\n");

  // one <Piece> entry per per-processor piece file

  int npieces = (multiproc > 1) ? multiproc : nprocs;
  for (int i = 0; i < npieces; ++i)
    utils::print(fp, R"(    <Piece Source="{}"/>)" "\n", pvtk_piece_filename(i));

  utils::print(fp, R"(  </{}>)" "\n", gridtype);
  utils::print(fp, R"(</VTKFile>)" "\n");
}

/* ---------------------------------------------------------------------- */

void DumpVTK::reset_vtk_data_containers()
{
  points.clear();
  myarrays.clear();

  // build the list of data arrays in the order their values appear in buf.
  // the first three fields are the x,y,z coordinates, which become the
  // points of the dataset, and three consecutive fields of a vector
  // attribute become a single array with three components.

  auto it=vtype.begin();
  ++it; ++it; ++it;
  for (; it != vtype.end(); ++it) {
    VTKArray array;
    array.name = name[it->first];
    array.type = it->second;
    array.ncomp = 1;

    if (vector_set.find(it->first) != vector_set.end()) {
      array.ncomp = 3;
      ++it; ++it;
    }
    myarrays.push_back(std::move(array));
  }
}

/* ---------------------------------------------------------------------- */

int DumpVTK::parse_vtk_fields(int narg, char **arg)
{

  pack_choice[X] = &DumpVTK::pack_x;
  vtype[X] = Dump::DOUBLE;
  name[X] = "x";
  pack_choice[Y] = &DumpVTK::pack_y;
  vtype[Y] = Dump::DOUBLE;
  name[Y] = "y";
  pack_choice[Z] = &DumpVTK::pack_z;
  vtype[Z] = Dump::DOUBLE;
  name[Z] = "z";

  // customize by adding to if statement

  for (int iarg = 0; iarg < narg; iarg++) {

    if (strcmp(arg[iarg],"id") == 0) {
      pack_choice[ID] = &DumpVTK::pack_id;
      vtype[ID] = Dump::INT;
      name[ID] = arg[iarg];
    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (!atom->molecule_flag)
        error->all(FLERR,"Dumping atom property 'mol' that isn't allocated");
      pack_choice[MOL] = &DumpVTK::pack_molecule;
      vtype[MOL] = Dump::INT;
      name[MOL] = arg[iarg];
    } else if (strcmp(arg[iarg],"proc") == 0) {
      pack_choice[PROC] = &DumpVTK::pack_proc;
      vtype[PROC] = Dump::INT;
      name[PROC] = arg[iarg];
    } else if (strcmp(arg[iarg],"procp1") == 0) {
      pack_choice[PROCP1] = &DumpVTK::pack_procp1;
      vtype[PROCP1] = Dump::INT;
      name[PROCP1] = arg[iarg];
    } else if (strcmp(arg[iarg],"type") == 0) {
      pack_choice[TYPE] = &DumpVTK::pack_type;
      vtype[TYPE] = Dump::INT;
      name[TYPE] =arg[iarg];
    } else if (strcmp(arg[iarg],"element") == 0) {
      pack_choice[ELEMENT] = &DumpVTK::pack_type;
      vtype[ELEMENT] = Dump::STRING;
      name[ELEMENT] = arg[iarg];
    } else if (strcmp(arg[iarg],"mass") == 0) {
      pack_choice[MASS] = &DumpVTK::pack_mass;
      vtype[MASS] = Dump::DOUBLE;
      name[MASS] = arg[iarg];

    } else if (strcmp(arg[iarg],"x") == 0) {
      // required property
    } else if (strcmp(arg[iarg],"y") == 0) {
      // required property
    } else if (strcmp(arg[iarg],"z") == 0) {
      // required property
    } else if (strcmp(arg[iarg],"xs") == 0) {
      if (domain->triclinic) pack_choice[XS] = &DumpVTK::pack_xs_triclinic;
      else pack_choice[XS] = &DumpVTK::pack_xs;
      vtype[XS] = Dump::DOUBLE;
      name[XS] = arg[iarg];
    } else if (strcmp(arg[iarg],"ys") == 0) {
      if (domain->triclinic) pack_choice[YS] = &DumpVTK::pack_ys_triclinic;
      else pack_choice[YS] = &DumpVTK::pack_ys;
      vtype[YS] = Dump::DOUBLE;
      name[YS] = arg[iarg];
    } else if (strcmp(arg[iarg],"zs") == 0) {
      if (domain->triclinic) pack_choice[ZS] = &DumpVTK::pack_zs_triclinic;
      else pack_choice[ZS] = &DumpVTK::pack_zs;
      vtype[ZS] = Dump::DOUBLE;
      name[ZS] = arg[iarg];
    } else if (strcmp(arg[iarg],"xu") == 0) {
      if (domain->triclinic) pack_choice[XU] = &DumpVTK::pack_xu_triclinic;
      else pack_choice[XU] = &DumpVTK::pack_xu;
      vtype[XU] = Dump::DOUBLE;
      name[XU] = arg[iarg];
    } else if (strcmp(arg[iarg],"yu") == 0) {
      if (domain->triclinic) pack_choice[YU] = &DumpVTK::pack_yu_triclinic;
      else pack_choice[YU] = &DumpVTK::pack_yu;
      vtype[YU] = Dump::DOUBLE;
      name[YU] = arg[iarg];
    } else if (strcmp(arg[iarg],"zu") == 0) {
      if (domain->triclinic) pack_choice[ZU] = &DumpVTK::pack_zu_triclinic;
      else pack_choice[ZU] = &DumpVTK::pack_zu;
      vtype[ZU] = Dump::DOUBLE;
      name[ZU] = arg[iarg];
    } else if (strcmp(arg[iarg],"xsu") == 0) {
      if (domain->triclinic) pack_choice[XSU] = &DumpVTK::pack_xsu_triclinic;
      else pack_choice[XSU] = &DumpVTK::pack_xsu;
      vtype[XSU] = Dump::DOUBLE;
      name[XSU] = arg[iarg];
    } else if (strcmp(arg[iarg],"ysu") == 0) {
      if (domain->triclinic) pack_choice[YSU] = &DumpVTK::pack_ysu_triclinic;
      else pack_choice[YSU] = &DumpVTK::pack_ysu;
      vtype[YSU] = Dump::DOUBLE;
      name[YSU] = arg[iarg];
    } else if (strcmp(arg[iarg],"zsu") == 0) {
      if (domain->triclinic) pack_choice[ZSU] = &DumpVTK::pack_zsu_triclinic;
      else pack_choice[ZSU] = &DumpVTK::pack_zsu;
      vtype[ZSU] = Dump::DOUBLE;
      name[ZSU] = arg[iarg];
    } else if (strcmp(arg[iarg],"ix") == 0) {
      pack_choice[IX] = &DumpVTK::pack_ix;
      vtype[IX] = Dump::INT;
      name[IX] = arg[iarg];
    } else if (strcmp(arg[iarg],"iy") == 0) {
      pack_choice[IY] = &DumpVTK::pack_iy;
      vtype[IY] = Dump::INT;
      name[IY] = arg[iarg];
    } else if (strcmp(arg[iarg],"iz") == 0) {
      pack_choice[IZ] = &DumpVTK::pack_iz;
      vtype[IZ] = Dump::INT;
      name[IZ] = arg[iarg];

    } else if (strcmp(arg[iarg],"vx") == 0) {
      pack_choice[VX] = &DumpVTK::pack_vx;
      vtype[VX] = Dump::DOUBLE;
      name[VX] = arg[iarg];
    } else if (strcmp(arg[iarg],"vy") == 0) {
      pack_choice[VY] = &DumpVTK::pack_vy;
      vtype[VY] = Dump::DOUBLE;
      name[VY] = arg[iarg];
    } else if (strcmp(arg[iarg],"vz") == 0) {
      pack_choice[VZ] = &DumpVTK::pack_vz;
      vtype[VZ] = Dump::DOUBLE;
      name[VZ] = arg[iarg];
    } else if (strcmp(arg[iarg],"fx") == 0) {
      pack_choice[FX] = &DumpVTK::pack_fx;
      vtype[FX] = Dump::DOUBLE;
      name[FX] = arg[iarg];
    } else if (strcmp(arg[iarg],"fy") == 0) {
      pack_choice[FY] = &DumpVTK::pack_fy;
      vtype[FY] = Dump::DOUBLE;
      name[FY] = arg[iarg];
    } else if (strcmp(arg[iarg],"fz") == 0) {
      pack_choice[FZ] = &DumpVTK::pack_fz;
      vtype[FZ] = Dump::DOUBLE;
      name[FZ] = arg[iarg];
    } else if (strcmp(arg[iarg],"q") == 0) {
      if (!atom->q_flag)
        error->all(FLERR,"Dumping atom property '{}' that isn't allocated", arg[iarg]);
      pack_choice[Q] = &DumpVTK::pack_q;
      vtype[Q] = Dump::DOUBLE;
      name[Q] = arg[iarg];
    } else if (strcmp(arg[iarg],"mux") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping atom property '{}' that isn't allocated", arg[iarg]);
      pack_choice[MUX] = &DumpVTK::pack_mux;
      vtype[MUX] = Dump::DOUBLE;
      name[MUX] = arg[iarg];
    } else if (strcmp(arg[iarg],"muy") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping atom property '{}' that isn't allocated", arg[iarg]);
      pack_choice[MUY] = &DumpVTK::pack_muy;
      vtype[MUY] = Dump::DOUBLE;
      name[MUY] = arg[iarg];
    } else if (strcmp(arg[iarg],"muz") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping atom property '{}' that isn't allocated", arg[iarg]);
      pack_choice[MUZ] = &DumpVTK::pack_muz;
      vtype[MUZ] = Dump::DOUBLE;
      name[MUZ] = arg[iarg];
    } else if (strcmp(arg[iarg],"mu") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping atom property '{}' that isn't allocated", arg[iarg]);
      pack_choice[MU] = &DumpVTK::pack_mu;
      vtype[MU] = Dump::DOUBLE;
      name[MU] = arg[iarg];
    } else if (strcmp(arg[iarg],"radius") == 0) {
      if (!atom->radius_flag)
        error->all(FLERR,"Dumping atom property '{}' that isn't allocated", arg[iarg]);
      pack_choice[RADIUS] = &DumpVTK::pack_radius;
      vtype[RADIUS] = Dump::DOUBLE;
      name[RADIUS] = arg[iarg];
    } else if (strcmp(arg[iarg],"diameter") == 0) {
      if (!atom->radius_flag)
        error->all(FLERR,"Dumping atom property '{}' that isn't allocated", arg[iarg]);
      pack_choice[DIAMETER] = &DumpVTK::pack_diameter;
      vtype[DIAMETER] = Dump::DOUBLE;
      name[DIAMETER] = arg[iarg];
    } else if (strcmp(arg[iarg],"omegax") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,"Dumping atom property '{}' that isn't allocated", arg[iarg]);
      pack_choice[OMEGAX] = &DumpVTK::pack_omegax;
      vtype[OMEGAX] = Dump::DOUBLE;
      name[OMEGAX] = arg[iarg];
    } else if (strcmp(arg[iarg],"omegay") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,"Dumping atom property '{}' that isn't allocated", arg[iarg]);
      pack_choice[OMEGAY] = &DumpVTK::pack_omegay;
      vtype[OMEGAY] = Dump::DOUBLE;
      name[OMEGAY] = arg[iarg];
    } else if (strcmp(arg[iarg],"omegaz") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,"Dumping atom property '{}' that isn't allocated", arg[iarg]);
      pack_choice[OMEGAZ] = &DumpVTK::pack_omegaz;
      vtype[OMEGAZ] = Dump::DOUBLE;
      name[OMEGAZ] = arg[iarg];
    } else if (strcmp(arg[iarg],"angmomx") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,"Dumping atom property '{}' that isn't allocated", arg[iarg]);
      pack_choice[ANGMOMX] = &DumpVTK::pack_angmomx;
      vtype[ANGMOMX] = Dump::DOUBLE;
      name[ANGMOMX] = arg[iarg];
    } else if (strcmp(arg[iarg],"angmomy") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,"Dumping atom property '{}' that isn't allocated", arg[iarg]);
      pack_choice[ANGMOMY] = &DumpVTK::pack_angmomy;
      vtype[ANGMOMY] = Dump::DOUBLE;
      name[ANGMOMY] = arg[iarg];
    } else if (strcmp(arg[iarg],"angmomz") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,"Dumping atom property '{}' that isn't allocated", arg[iarg]);
      pack_choice[ANGMOMZ] = &DumpVTK::pack_angmomz;
      vtype[ANGMOMZ] = Dump::DOUBLE;
      name[ANGMOMZ] = arg[iarg];
    } else if (strcmp(arg[iarg],"tqx") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,"Dumping atom property '{}' that isn't allocated", arg[iarg]);
      pack_choice[TQX] = &DumpVTK::pack_tqx;
      vtype[TQX] = Dump::DOUBLE;
      name[TQX] = arg[iarg];
    } else if (strcmp(arg[iarg],"tqy") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,"Dumping atom property '{}' that isn't allocated", arg[iarg]);
      pack_choice[TQY] = &DumpVTK::pack_tqy;
      vtype[TQY] = Dump::DOUBLE;
      name[TQY] = arg[iarg];
    } else if (strcmp(arg[iarg],"tqz") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,"Dumping atom property '{}' that isn't allocated", arg[iarg]);
      pack_choice[TQZ] = &DumpVTK::pack_tqz;
      vtype[TQZ] = Dump::DOUBLE;
      name[TQZ] = arg[iarg];

    } else {

      int n,flag,cols;
      ArgInfo argi(arg[iarg],ArgInfo::COMPUTE|ArgInfo::FIX|ArgInfo::VARIABLE
                   |ArgInfo::DNAME|ArgInfo::INAME);
      argindex[ATTRIBUTES+iarg] = argi.get_index1();
      const auto *aname = argi.get_name();

      switch (argi.get_type()) {

      case ArgInfo::UNKNOWN:
        error->all(FLERR,"Invalid attribute in dump vtk command: {}",arg[iarg]);
        break;

      // compute value = c_ID
      // if no trailing [], then arg is set to 0, else arg is int between []

      case ArgInfo::COMPUTE:
      {
        pack_choice[ATTRIBUTES+iarg] = &DumpVTK::pack_vtk_compute;
        vtype[ATTRIBUTES+iarg] = Dump::DOUBLE;

        auto *icompute = modify->get_compute_by_id(aname);
        if (!icompute) {
          error->all(FLERR,"Could not find dump vtk compute ID: {}",aname);
        } else {
          if (icompute->peratom_flag == 0)
            error->all(FLERR,"Dump vtk compute {} does not compute per-atom info",aname);
          if (argi.get_dim() == 0 && icompute->size_peratom_cols > 0)
            error->all(FLERR,"Dump vtk compute {} does not calculate per-atom vector",aname);
          if (argi.get_dim() > 0 && icompute->size_peratom_cols == 0)
            error->all(FLERR,"Dump vtk compute {} does not calculate per-atom array",aname);
          if (argi.get_dim() > 0 && argi.get_index1() > icompute->size_peratom_cols)
            error->all(FLERR,"Dump vtk compute {} vector is accessed out-of-range{}",
                       aname, utils::errorurl(20));
          field2index[ATTRIBUTES+iarg] = add_compute(aname);
          name[ATTRIBUTES+iarg] = arg[iarg];
        }
        break;
      }

      // fix value = f_ID
      // if no trailing [], then arg is set to 0, else arg is between []

      case ArgInfo::FIX:
      {
        pack_choice[ATTRIBUTES+iarg] = &DumpVTK::pack_vtk_fix;
        vtype[ATTRIBUTES+iarg] = Dump::DOUBLE;

        auto *ifix = modify->get_fix_by_id(aname);
        if (!ifix) {
          error->all(FLERR,"Could not find dump vtk fix ID: {}",aname);
        } else {
          if (ifix->peratom_flag == 0)
            error->all(FLERR,"Dump vtk fix {} does not compute per-atom info",aname);
          if (argi.get_dim() == 0 && ifix->size_peratom_cols > 0)
            error->all(FLERR,"Dump vtk fix {} does not compute per-atom vector",aname);
          if (argi.get_dim() > 0 && ifix->size_peratom_cols == 0)
            error->all(FLERR,"Dump vtk fix {} does not compute per-atom array",aname);
          if (argi.get_dim() > 0 && argi.get_index1() > ifix->size_peratom_cols)
            error->all(FLERR,"Dump vtk fix {} vector is accessed out-of-range{}",
                       aname, utils::errorurl(20));

          field2index[ATTRIBUTES+iarg] = add_fix(aname);
          name[ATTRIBUTES+iarg] = arg[iarg];
        }
        break;
      }

      // variable value = v_name

      case ArgInfo::VARIABLE:
        pack_choice[ATTRIBUTES+iarg] = &DumpVTK::pack_vtk_variable;
        vtype[ATTRIBUTES+iarg] = Dump::DOUBLE;

        n = input->variable->find(aname);
        if (n < 0) error->all(FLERR,"Could not find dump vtk variable name {}",aname);
        if (input->variable->atomstyle(n) == 0)
          error->all(FLERR,"Dump vtk variable {} is not atom-style variable",aname);

        field2index[ATTRIBUTES+iarg] = add_variable(aname);
        name[ATTRIBUTES+iarg] = arg[iarg];
        break;

      // custom per-atom floating point vector or array = d_ID d2_ID

      case ArgInfo::DNAME:
        pack_choice[ATTRIBUTES+iarg] = &DumpVTK::pack_vtk_custom;
        vtype[ATTRIBUTES+iarg] = Dump::DOUBLE;

        n = atom->find_custom(aname,flag,cols);

        if (n < 0)
          error->all(FLERR,"Could not find custom per-atom property ID: {}", aname);
        if (argindex[ATTRIBUTES+iarg] == 0) {
          if (!flag || cols)
            error->all(FLERR,"Property double vector {} for dump vtk does not exist",aname);
        } else {
          if (!flag || !cols)
            error->all(FLERR,"Property double array {} for dump vtk does not exist",aname);
          if (argindex[ATTRIBUTES+iarg] > atom->dcols[n])
            error->all(FLERR,"Dump vtk property array {} is accessed out-of-range{}",aname,
                       utils::errorurl(20));
        }
        field2index[ATTRIBUTES+iarg] = add_custom(aname,1);
        name[ATTRIBUTES+iarg] = arg[iarg];
        break;

      // custom per-atom integer vector or array = i_ID or i2_ID

      case ArgInfo::INAME:
        pack_choice[ATTRIBUTES+iarg] = &DumpVTK::pack_vtk_custom;
        vtype[ATTRIBUTES+iarg] = Dump::INT;

        n = atom->find_custom(aname,flag,cols);

        if (n < 0)
          error->all(FLERR,"Could not find custom per-atom property ID: {}", aname);
        if (argindex[ATTRIBUTES+iarg] == 0) {
          if (flag || cols)
            error->all(FLERR,"Property integer vector {} for dump vtk does not exist",aname);
        } else {
          if (flag || !cols)
            error->all(FLERR,"Property integer array {} for dump vtk does not exist",aname);
          if (argindex[ATTRIBUTES+iarg] > atom->icols[n])
            error->all(FLERR,"Dump vtk property array {} is accessed out-of-range{}",aname,
                       utils::errorurl(20));
        }
        field2index[ATTRIBUTES+iarg] = add_custom(aname,0);
        name[ATTRIBUTES+iarg] = arg[iarg];
        break;

      // no match

      default:
        return iarg;
        break;
      }
    }
  }

  identify_vectors();

  return narg;
}

/* ---------------------------------------------------------------------- */

void DumpVTK::identify_vectors()
{
  // detect vectors
  vector_set.insert(X); // required

  int vector3_starts[] = {XS, XU, XSU, IX, VX, FX, MUX, OMEGAX, ANGMOMX, TQX};
  int num_vector3_starts = sizeof(vector3_starts) / sizeof(int);

  for (int v3s = 0; v3s < num_vector3_starts; v3s++) {
    if (name.count(vector3_starts[v3s]  ) &&
       name.count(vector3_starts[v3s]+1) &&
       name.count(vector3_starts[v3s]+2) )
    {
      std::string vectorName = name[vector3_starts[v3s]];
      std::string::size_type erase_start = vectorName.find_first_of('x');
      if (erase_start == 0) {
        vectorName.erase(0,1);
      } else {
        vectorName.erase(erase_start);
      }
      name[vector3_starts[v3s]] = std::move(vectorName);
      vector_set.insert(vector3_starts[v3s]);
    }
  }

  // compute and fix vectors
  for (auto it=name.begin(); it != name.end(); ++it) {
    if (it->first < ATTRIBUTES) // neither fix nor compute
      continue;

    if (argindex[it->first] == 0) // single value
      continue;

    // assume components are grouped together and in correct order
    if (name.count(it->first + 1) && name.count(it->first + 2)) { // more attributes?
      if (it->second.compare(0,it->second.length()-3,name[it->first + 1],0,it->second.length()-3) == 0  && // same attributes?
         it->second.compare(0,it->second.length()-3,name[it->first + 2],0,it->second.length()-3) == 0)
      {
        it->second.erase(it->second.length()-1);
        std::ostringstream oss;
        oss << "-" << argindex[it->first+2] << "]";
        it->second += oss.str();
        vector_set.insert(it->first);
        ++it; ++it;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int DumpVTK::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"region") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) {
      delete[] idregion;
      idregion = nullptr;
    } else {
      if (!domain->get_region_by_id(arg[1]))
        error->all(FLERR,"Dump_modify region {} does not exist",arg[1]);
      delete[] idregion;
      idregion = utils::strdup(arg[1]);
    }
    return 2;
  }

  if (strcmp(arg[0],"label") == 0) {
     if (narg < 2) utils::missing_cmd_args(FLERR, "dump_modify label", error);
     delete[] label;
     label = utils::strdup(arg[1]);
     return 2;
   }

  if (strcmp(arg[0],"binary") == 0) {
     if (narg < 2) utils::missing_cmd_args(FLERR, "dump_modify binary", error);
     binary = utils::logical(FLERR,arg[1],false,lmp);
     return 2;
  }

  if (strcmp(arg[0],"double") == 0) {
     if (narg < 2) utils::missing_cmd_args(FLERR, "dump_modify double", error);
     writeprec = utils::logical(FLERR,arg[1],false,lmp) ? VTKWriter::DOUBLE : VTKWriter::SINGLE;
     return 2;
  }

  if (strcmp(arg[0],"element") == 0) {
    if (narg < ntypes+1)
      error->all(FLERR,"Number of dump_modify element names does not match number of atom types");

    for (int i = 1; i <= ntypes; i++) delete[] typenames[i];
    delete[] typenames;
    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      typenames[itype] = utils::strdup(arg[itype]);
    }
    return ntypes+1;
  }

  if (strcmp(arg[0],"refresh") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    ArgInfo argi(arg[1],ArgInfo::COMPUTE);
    if ((argi.get_type() != ArgInfo::COMPUTE) || (argi.get_dim() != 0))
      error->all(FLERR,"Illegal dump_modify command");
    if (refreshflag) error->all(FLERR,"Dump_modify can only have one refresh");

    refreshflag = 1;
    idrefresh = argi.copy_name();
    return 2;
  }

  if (strcmp(arg[0],"thresh") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) {
      if (nthresh) {
        memory->destroy(thresh_array);
        memory->destroy(thresh_op);
        memory->destroy(thresh_value);
        thresh_array = nullptr;
        thresh_op = nullptr;
        thresh_value = nullptr;
        thresh_last = nullptr;
        for (int i = 0; i < nthreshlast; i++) {
          modify->delete_fix(thresh_fixID[i]);
          delete[] thresh_fixID[i];
        }
        thresh_fix = nullptr;
        thresh_fixID = nullptr;
        thresh_first = nullptr;
      }
      nthresh = nthreshlast = 0;
      return 2;
    }

    if (narg < 4) error->all(FLERR,"Illegal dump_modify command");

    // grow threshold arrays

    memory->grow(thresh_array,nthresh+1,"dump:thresh_array");
    memory->grow(thresh_op,(nthresh+1),"dump:thresh_op");
    memory->grow(thresh_value,(nthresh+1),"dump:thresh_value");
    memory->grow(thresh_last,(nthresh+1),"dump:thresh_last");

    // set attribute type of threshold
    // customize by adding to if statement

    if (strcmp(arg[1],"id") == 0) thresh_array[nthresh] = ID;
    else if (strcmp(arg[1],"mol") == 0) thresh_array[nthresh] = MOL;
    else if (strcmp(arg[1],"proc") == 0) thresh_array[nthresh] = PROC;
    else if (strcmp(arg[1],"procp1") == 0) thresh_array[nthresh] = PROCP1;
    else if (strcmp(arg[1],"type") == 0) thresh_array[nthresh] = TYPE;
    else if (strcmp(arg[1],"mass") == 0) thresh_array[nthresh] = MASS;

    else if (strcmp(arg[1],"x") == 0) thresh_array[nthresh] = X;
    else if (strcmp(arg[1],"y") == 0) thresh_array[nthresh] = Y;
    else if (strcmp(arg[1],"z") == 0) thresh_array[nthresh] = Z;

    else if (strcmp(arg[1],"xs") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = XS;
    else if (strcmp(arg[1],"xs") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = XSTRI;
    else if (strcmp(arg[1],"ys") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = YS;
    else if (strcmp(arg[1],"ys") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = YSTRI;
    else if (strcmp(arg[1],"zs") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = ZS;
    else if (strcmp(arg[1],"zs") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = ZSTRI;

    else if (strcmp(arg[1],"xu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = XU;
    else if (strcmp(arg[1],"xu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = XUTRI;
    else if (strcmp(arg[1],"yu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = YU;
    else if (strcmp(arg[1],"yu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = YUTRI;
    else if (strcmp(arg[1],"zu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = ZU;
    else if (strcmp(arg[1],"zu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = ZUTRI;

    else if (strcmp(arg[1],"xsu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = XSU;
    else if (strcmp(arg[1],"xsu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = XSUTRI;
    else if (strcmp(arg[1],"ysu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = YSU;
    else if (strcmp(arg[1],"ysu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = YSUTRI;
    else if (strcmp(arg[1],"zsu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = ZSU;
    else if (strcmp(arg[1],"zsu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = ZSUTRI;

    else if (strcmp(arg[1],"ix") == 0) thresh_array[nthresh] = IX;
    else if (strcmp(arg[1],"iy") == 0) thresh_array[nthresh] = IY;
    else if (strcmp(arg[1],"iz") == 0) thresh_array[nthresh] = IZ;
    else if (strcmp(arg[1],"vx") == 0) thresh_array[nthresh] = VX;
    else if (strcmp(arg[1],"vy") == 0) thresh_array[nthresh] = VY;
    else if (strcmp(arg[1],"vz") == 0) thresh_array[nthresh] = VZ;
    else if (strcmp(arg[1],"fx") == 0) thresh_array[nthresh] = FX;
    else if (strcmp(arg[1],"fy") == 0) thresh_array[nthresh] = FY;
    else if (strcmp(arg[1],"fz") == 0) thresh_array[nthresh] = FZ;

    else if (strcmp(arg[1],"q") == 0) thresh_array[nthresh] = Q;
    else if (strcmp(arg[1],"mux") == 0) thresh_array[nthresh] = MUX;
    else if (strcmp(arg[1],"muy") == 0) thresh_array[nthresh] = MUY;
    else if (strcmp(arg[1],"muz") == 0) thresh_array[nthresh] = MUZ;
    else if (strcmp(arg[1],"mu") == 0) thresh_array[nthresh] = MU;

    else if (strcmp(arg[1],"radius") == 0) thresh_array[nthresh] = RADIUS;
    else if (strcmp(arg[1],"diameter") == 0) thresh_array[nthresh] = DIAMETER;
    else if (strcmp(arg[1],"omegax") == 0) thresh_array[nthresh] = OMEGAX;
    else if (strcmp(arg[1],"omegay") == 0) thresh_array[nthresh] = OMEGAY;
    else if (strcmp(arg[1],"omegaz") == 0) thresh_array[nthresh] = OMEGAZ;
    else if (strcmp(arg[1],"angmomx") == 0) thresh_array[nthresh] = ANGMOMX;
    else if (strcmp(arg[1],"angmomy") == 0) thresh_array[nthresh] = ANGMOMY;
    else if (strcmp(arg[1],"angmomz") == 0) thresh_array[nthresh] = ANGMOMZ;
    else if (strcmp(arg[1],"tqx") == 0) thresh_array[nthresh] = TQX;
    else if (strcmp(arg[1],"tqy") == 0) thresh_array[nthresh] = TQY;
    else if (strcmp(arg[1],"tqz") == 0) thresh_array[nthresh] = TQZ;

    // compute or fix or variable or custom vector/array

    else {
      int n,flag,cols;
      ArgInfo argi(arg[1],ArgInfo::COMPUTE|ArgInfo::FIX|ArgInfo::VARIABLE
                   |ArgInfo::DNAME|ArgInfo::INAME);
      argindex[ATTRIBUTES+nfield+nthresh] = argi.get_index1();
      const auto *aname = argi.get_name();

      switch (argi.get_type()) {

      case ArgInfo::UNKNOWN:
        error->all(FLERR,"Invalid attribute in dump modify command");
        break;

      // compute value = c_ID
      // if no trailing [], then arg is set to 0, else arg is between []

      case ArgInfo::COMPUTE:
      {
        thresh_array[nthresh] = COMPUTE;
        auto *icompute = modify->get_compute_by_id(aname);
        if (!icompute) {
          error->all(FLERR,"Could not find dump modify compute ID: {}",aname);
        } else {
          if (icompute->peratom_flag == 0)
            error->all(FLERR,"Dump modify compute ID {} does not compute per-atom info",aname);
          if (argi.get_dim() == 0 && icompute->size_peratom_cols > 0)
            error->all(FLERR,"Dump modify compute ID {} does not compute per-atom vector",aname);
          if (argi.get_index1() > 0 && icompute->size_peratom_cols == 0)
            error->all(FLERR,"Dump modify compute ID {} does not compute per-atom array",aname);
          if (argi.get_index1() > 0 &&
              argi.get_index1() > icompute->size_peratom_cols)
            error->all(FLERR,"Dump modify compute ID {} vector is not large enough",aname);

          field2index[ATTRIBUTES+nfield+nthresh] = add_compute(aname);
        }
        break;
      }

      // fix value = f_ID
      // if no trailing [], then arg is set to 0, else arg is between []

      case ArgInfo::FIX:
      {
        thresh_array[nthresh] = FIX;
        auto *ifix = modify->get_fix_by_id(aname);
        if (!ifix) {
          error->all(FLERR,"Could not find dump modify fix ID: {}",aname);
        } else {
          if (ifix->peratom_flag == 0)
            error->all(FLERR,"Dump modify fix ID {} does not compute per-atom info",aname);
          if (argi.get_dim() == 0 && ifix->size_peratom_cols > 0)
            error->all(FLERR,"Dump modify fix ID {} does not compute per-atom vector",aname);
          if (argi.get_index1() > 0 && ifix->size_peratom_cols == 0)
            error->all(FLERR,"Dump modify fix ID {} does not compute per-atom array",aname);
          if (argi.get_index1() > 0 && argi.get_index1() > ifix->size_peratom_cols)
            error->all(FLERR,"Dump modify fix ID {} vector is not large enough",aname);

          field2index[ATTRIBUTES+nfield+nthresh] = add_fix(aname);
        }
        break;
      }

      // variable value = v_ID

      case ArgInfo::VARIABLE:
        thresh_array[nthresh] = VARIABLE;
        n = input->variable->find(aname);
        if (n < 0) error->all(FLERR,"Could not find dump modify variable name: {}",aname);
        if (input->variable->atomstyle(n) == 0)
          error->all(FLERR,"Dump modify variable {} is not atom-style variable",aname);

        field2index[ATTRIBUTES+nfield+nthresh] = add_variable(aname);
        break;

      // custom per atom floating point vector or array

      case ArgInfo::DNAME:
        n = atom->find_custom(aname,flag,cols);

        if (n < 0)
          error->all(FLERR,"Could not find custom per-atom property ID: {}", aname);
        if (argindex[ATTRIBUTES+nfield+nthresh] == 0) {
          if (!flag || cols)
            error->all(FLERR,"Property double vector for dump custom does not exist");
          thresh_array[nthresh] = DVEC;
        } else {
          if (!flag || !cols)
            error->all(FLERR,"Property double array for dump custom does not exist");
          if (argindex[ATTRIBUTES+nfield+nthresh] > atom->dcols[n])
            error->all(FLERR,"Dump custom property array is accessed out-of-range");
          thresh_array[nthresh] = DARRAY;
        }

        field2index[ATTRIBUTES+nfield+nthresh] = add_custom(aname,thresh_array[nthresh]);
        break;

      // custom per atom integer vector or array

      case ArgInfo::INAME:
        n = atom->find_custom(aname,flag,cols);

        if (n < 0)
          error->all(FLERR,"Could not find custom per-atom property ID: {}", aname);
        if (argindex[ATTRIBUTES+nfield+nthresh] == 0) {
          if (flag || cols)
            error->all(FLERR,"Property integer vector for dump custom does not exist");
          thresh_array[nthresh] = IVEC;
        } else {
          if (flag || !cols)
            error->all(FLERR,"Property integer array for dump custom does not exist");
          if (argindex[ATTRIBUTES+nfield+nthresh] > atom->icols[n])
            error->all(FLERR,"Dump custom property array is accessed out-of-range");
          thresh_array[nthresh] = IARRAY;
        }

        field2index[ATTRIBUTES+nfield+nthresh] = add_custom(aname,thresh_array[nthresh]);
        break;

      // no match

      default:
        error->all(FLERR,"Invalid dump_modify thresh attribute: {}",aname);
        break;
      }
    }

    // set operation type of threshold

    if (strcmp(arg[2],"<") == 0) thresh_op[nthresh] = LT;
    else if (strcmp(arg[2],"<=") == 0) thresh_op[nthresh] = LE;
    else if (strcmp(arg[2],">") == 0) thresh_op[nthresh] = GT;
    else if (strcmp(arg[2],">=") == 0) thresh_op[nthresh] = GE;
    else if (strcmp(arg[2],"==") == 0) thresh_op[nthresh] = EQ;
    else if (strcmp(arg[2],"!=") == 0) thresh_op[nthresh] = NEQ;
    else if (strcmp(arg[2],"|^") == 0) thresh_op[nthresh] = XOR;
    else error->all(FLERR,"Invalid dump_modify thresh operator");

    // set threshold value as number or special LAST keyword
    // create FixStore to hold LAST values, should work with restart
    // id = dump-ID + nthreshlast + DUMP_STORE, fix group = dump group

    if (strcmp(arg[3],"LAST") != 0) {
      thresh_value[nthresh] = utils::numeric(FLERR,arg[3],false,lmp);
      thresh_last[nthresh] = -1;
    } else {
      thresh_fix = (FixStoreAtom **)
        memory->srealloc(thresh_fix,(nthreshlast+1)*sizeof(FixStoreAtom *),"dump:thresh_fix");
      thresh_fixID = (char **)
        memory->srealloc(thresh_fixID,(nthreshlast+1)*sizeof(char *),"dump:thresh_fixID");
      memory->grow(thresh_first,(nthreshlast+1),"dump:thresh_first");

      std::string threshid = fmt::format("{}{}_DUMP_STORE",id,nthreshlast);
      thresh_fixID[nthreshlast] = utils::strdup(threshid);
      threshid += fmt::format(" {} STORE/ATOM 1 0 0 1", group->names[igroup]);
      thresh_fix[nthreshlast] = dynamic_cast<FixStoreAtom *>(modify->add_fix(threshid));

      thresh_last[nthreshlast] = nthreshlast;
      thresh_first[nthreshlast] = 1;
      nthreshlast++;
    }

    nthresh++;
    return 4;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   extraction of Compute, Fix, Variable results
------------------------------------------------------------------------- */

void DumpVTK::pack_vtk_compute(int n)
{
  double *vector = compute[field2index[current_pack_choice_key]]->vector_atom;
  double **array = compute[field2index[current_pack_choice_key]]->array_atom;
  int index = argindex[current_pack_choice_key];

  if (index == 0) {
    for (int i = 0; i < nchoose; i++) {
      buf[n] = vector[clist[i]];
      n += size_one;
    }
  } else {
    index--;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = array[clist[i]][index];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpVTK::pack_vtk_fix(int n)
{
  double *vector = fix[field2index[current_pack_choice_key]]->vector_atom;
  double **array = fix[field2index[current_pack_choice_key]]->array_atom;
  int index = argindex[current_pack_choice_key];

  if (index == 0) {
    for (int i = 0; i < nchoose; i++) {
      buf[n] = vector[clist[i]];
      n += size_one;
    }
  } else {
    index--;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = array[clist[i]][index];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpVTK::pack_vtk_variable(int n)
{
  double *vector = vbuf[field2index[current_pack_choice_key]];

  for (int i = 0; i < nchoose; i++) {
    buf[n] = vector[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpVTK::pack_vtk_custom(int n)
{
  int flag = custom_flag[field2index[current_pack_choice_key]];
  int iwhich = custom[field2index[current_pack_choice_key]];
  int index = argindex[current_pack_choice_key];

  if (flag == IVEC) {
    int *ivector = atom->ivector[iwhich];
    for (int i = 0; i < nchoose; i++) {
      buf[n] = ivector[clist[i]];
      n += size_one;
    }
  } else if (flag == DVEC) {
    double *dvector = atom->dvector[iwhich];
    for (int i = 0; i < nchoose; i++) {
      buf[n] = dvector[clist[i]];
      n += size_one;
    }
  } else if (flag == IARRAY) {
    index--;
    int **iarray = atom->iarray[iwhich];
    for (int i = 0; i < nchoose; i++) {
      buf[n] = iarray[clist[i]][index];
      n += size_one;
    }
  } else if (flag == DARRAY) {
    index--;
    double **darray = atom->darray[iwhich];
    for (int i = 0; i < nchoose; i++) {
      buf[n] = darray[clist[i]][index];
      n += size_one;
    }
  }
}
