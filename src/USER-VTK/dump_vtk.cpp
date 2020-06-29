/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "dump_vtk.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

#include <vector>
#include <sstream>
#include <vtkVersion.h>

#ifndef VTK_MAJOR_VERSION
#include <vtkConfigure.h>
#endif

#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPPolyDataWriter.h>
#include <vtkRectilinearGrid.h>
#include <vtkRectilinearGridWriter.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

using namespace LAMMPS_NS;

// customize by
// * adding an enum constant (add vector components in consecutive order)
// * adding a pack_*(int) function for the value
// * adjusting parse_fields function to add the pack_* function to pack_choice
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
     VARIABLE,COMPUTE,FIX,INAME,DNAME,
     ATTRIBUTES}; // must come last
enum{LT,LE,GT,GE,EQ,NEQ};
enum{VTK,VTP,VTU,PVTP,PVTU}; // file formats

#define INVOKED_PERATOM 8
#define ONEFIELD 32
#define DELTA 1048576

#if (VTK_MAJOR_VERSION < 5) || (VTK_MAJOR_VERSION > 8)
#error This code has only been tested with VTK 5, 6, 7, and 8
#elif VTK_MAJOR_VERSION > 6
#define InsertNextTupleValue InsertNextTypedTuple
#endif

/* ---------------------------------------------------------------------- */

DumpVTK::DumpVTK(LAMMPS *lmp, int narg, char **arg) :
  DumpCustom(lmp, narg, arg)
{
  if (narg == 5) error->all(FLERR,"No dump vtk arguments specified");

  pack_choice.clear();
  vtype.clear();
  name.clear();

  myarrays.clear();
  n_calls_ = 0;

  // process attributes
  // ioptional = start of additional optional args
  // only dump image and dump movie styles process optional args

  ioptional = parse_fields(narg,arg);

  if (ioptional < narg &&
      strcmp(style,"image") != 0 && strcmp(style,"movie") != 0)
    error->all(FLERR,"Invalid attribute in dump vtk command");
  size_one = pack_choice.size();
  current_pack_choice_key = -1;

  if (filewriter) reset_vtk_data_containers();


  label = NULL;

  {
    // parallel vtp/vtu requires proc number to be preceded by underscore '_'
    multiname_ex = NULL;
    char *ptr = strchr(filename,'%');
    if (ptr) {
      multiname_ex = new char[strlen(filename) + 16];
      *ptr = '\0';
      sprintf(multiname_ex,"%s_%d%s",filename,me,ptr+1);
      *ptr = '%';
    }
  }

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
    if (me != 0) filewriter = 0;
    fileproc = 0;
    multiproc = 0;
    nclusterprocs = nprocs;
  }

  filecurrent = NULL;
  domainfilecurrent = NULL;
  parallelfilecurrent = NULL;
  header_choice = NULL;
  write_choice = NULL;
  boxcorners = NULL;
}

/* ---------------------------------------------------------------------- */

DumpVTK::~DumpVTK()
{
  delete [] filecurrent;
  delete [] domainfilecurrent;
  delete [] parallelfilecurrent;
  delete [] multiname_ex;
  delete [] label;
}

/* ---------------------------------------------------------------------- */

void DumpVTK::init_style()
{
  // default for element names = C

  if (typenames == NULL) {
    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      typenames[itype] = new char[2];
      strcpy(typenames[itype],"C");
    }
  }

  // setup boundary string

  domain->boundary_string(boundstr);

  // setup function ptrs

  header_choice = &DumpVTK::header_vtk;

  if (vtk_file_format == VTP || vtk_file_format == PVTP)
    write_choice = &DumpVTK::write_vtp;
  else if (vtk_file_format == VTU || vtk_file_format == PVTU)
    write_choice = &DumpVTK::write_vtu;
  else
    write_choice = &DumpVTK::write_vtk;

  // find current ptr for each compute,fix,variable
  // check that fix frequency is acceptable

  int icompute;
  for (int i = 0; i < ncompute; i++) {
    icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0) error->all(FLERR,"Could not find dump vtk compute ID");
    compute[i] = modify->compute[icompute];
  }

  int ifix;
  for (int i = 0; i < nfix; i++) {
    ifix = modify->find_fix(id_fix[i]);
    if (ifix < 0) error->all(FLERR,"Could not find dump vtk fix ID");
    fix[i] = modify->fix[ifix];
    if (nevery % modify->fix[ifix]->peratom_freq)
      error->all(FLERR,"Dump vtk and fix not computed at compatible times");
  }

  int ivariable;
  for (int i = 0; i < nvariable; i++) {
    ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0)
      error->all(FLERR,"Could not find dump vtk variable name");
    variable[i] = ivariable;
  }

  int icustom;
  for (int i = 0; i < ncustom; i++) {
    icustom = atom->find_custom(id_custom[i],flag_custom[i]);
    if (icustom < 0)
      error->all(FLERR,"Could not find custom per-atom property ID");
  }

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for dump vtk does not exist");
  }
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_header(bigint)
{
}

/* ---------------------------------------------------------------------- */

void DumpVTK::header_vtk(bigint)
{
}

/* ---------------------------------------------------------------------- */

int DumpVTK::count()
{
  n_calls_ = 0;

  int i;

  // grow choose and variable vbuf arrays if needed

  int nlocal = atom->nlocal;
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
  // only if within a run or minimize
  // else require that computes are current
  // this prevents a compute from being invoked by the WriteDump class

  if (ncompute) {
    if (update->whichflag == 0) {
      for (i = 0; i < ncompute; i++)
        if (compute[i]->invoked_peratom != update->ntimestep)
          error->all(FLERR,"Compute used in dump between runs is not current");
    } else {
      for (i = 0; i < ncompute; i++) {
        if (!(compute[i]->invoked_flag & INVOKED_PERATOM)) {
          compute[i]->compute_peratom();
          compute[i]->invoked_flag |= INVOKED_PERATOM;
        }
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

  if (iregion >= 0) {
    Region *region = domain->regions[iregion];
    region->prematch();
    double **x = atom->x;
    for (i = 0; i < nlocal; i++)
      if (choose[i] && region->match(x[i][0],x[i][1],x[i][2]) == 0)
        choose[i] = 0;
  }

  // un-choose if any threshold criterion isn't met

  if (nthresh) {
    double *ptr;
    double value;
    int nstride;
    int nlocal = atom->nlocal;

    for (int ithresh = 0; ithresh < nthresh; ithresh++) {

      // customize by adding to if statement

      if (thresh_array[ithresh] == ID) {
        tagint *tag = atom->tag;
        for (i = 0; i < nlocal; i++) dchoose[i] = tag[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == MOL) {
        if (!atom->molecule_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
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
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][0];
        nstride = 4;
      } else if (thresh_array[ithresh] == MUY) {
        if (!atom->mu_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][1];
        nstride = 4;
      } else if (thresh_array[ithresh] == MUZ) {
        if (!atom->mu_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][2];
        nstride = 4;
      } else if (thresh_array[ithresh] == MU) {
        if (!atom->mu_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][3];
        nstride = 4;

      } else if (thresh_array[ithresh] == RADIUS) {
        if (!atom->radius_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = atom->radius;
        nstride = 1;
      } else if (thresh_array[ithresh] == DIAMETER) {
        if (!atom->radius_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        double *radius = atom->radius;
        for (i = 0; i < nlocal; i++) dchoose[i] = 2.0*radius[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == OMEGAX) {
        if (!atom->omega_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->omega[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == OMEGAY) {
        if (!atom->omega_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->omega[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == OMEGAZ) {
        if (!atom->omega_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->omega[0][2];
        nstride = 3;
      } else if (thresh_array[ithresh] == ANGMOMX) {
        if (!atom->angmom_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->angmom[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == ANGMOMY) {
        if (!atom->angmom_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->angmom[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == ANGMOMZ) {
        if (!atom->angmom_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->angmom[0][2];
        nstride = 3;
      } else if (thresh_array[ithresh] == TQX) {
        if (!atom->torque_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->torque[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == TQY) {
        if (!atom->torque_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
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

      } else if (thresh_array[ithresh] == DNAME) {
        int iwhich,tmp;
        i = ATTRIBUTES + nfield + ithresh;
        iwhich = atom->find_custom(id_custom[field2index[i]],tmp);
        ptr = atom->dvector[iwhich];
        nstride = 1;

      } else if (thresh_array[ithresh] == INAME) {
        int iwhich,tmp;
        i = ATTRIBUTES + nfield + ithresh;
        iwhich = atom->find_custom(id_custom[field2index[i]],tmp);

        int *ivector = atom->ivector[iwhich];
        for (i = 0; i < nlocal; i++)
          dchoose[i] = ivector[i];
        ptr = dchoose;
        nstride = 1;
      }

      // unselect atoms that don't meet threshold criterion

      value = thresh_value[ithresh];

      switch (thresh_op[ithresh]) {
      case LT:
        for (i = 0; i < nlocal; i++, ptr += nstride)
          if (choose[i] && *ptr >= value) choose[i] = 0;
        break;
      case LE:
        for (i = 0; i < nlocal; i++, ptr += nstride)
          if (choose[i] && *ptr > value) choose[i] = 0;
        break;
      case GT:
        for (i = 0; i < nlocal; i++, ptr += nstride)
          if (choose[i] && *ptr <= value) choose[i] = 0;
        break;
      case GE:
        for (i = 0; i < nlocal; i++, ptr += nstride)
          if (choose[i] && *ptr < value) choose[i] = 0;
        break;
      case EQ:
        for (i = 0; i < nlocal; i++, ptr += nstride)
          if (choose[i] && *ptr != value) choose[i] = 0;
        break;
      case NEQ:
        for (i = 0; i < nlocal; i++, ptr += nstride)
          if (choose[i] && *ptr == value) choose[i] = 0;
        break;
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
    boxcorners = domain->corners;
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

  // write timestep header
  // for multiproc,
  //   nheader = # of lines in this file via Allreduce on clustercomm

  bigint nheader = ntotal;
  if (multiproc)
    MPI_Allreduce(&bnme,&nheader,1,MPI_LMP_BIGINT,MPI_SUM,clustercomm);

  if (filewriter) write_header(nheader);

  // insure buf is sized for packing and communicating
  // use nmax to insure filewriter proc can receive info from others
  // limit nmax*size_one to int since used as arg in MPI calls

  if (nmax > maxbuf) {
    if ((bigint) nmax * size_one > MAXSMALLINT)
      error->all(FLERR,"Too much per-proc info for dump");
    maxbuf = nmax;
    memory->destroy(buf);
    memory->create(buf,maxbuf*size_one,"dump:buf");
  }

  // insure ids buffer is sized for sorting

  if (sort_flag && sortcol == 0 && nmax > maxids) {
    maxids = nmax;
    memory->destroy(ids);
    memory->create(ids,maxids,"dump:ids");
  }

  // pack my data into buf
  // if sorting on IDs also request ID list from pack()
  // sort buf as needed

  if (sort_flag && sortcol == 0) pack(ids);
  else pack(NULL);
  if (sort_flag) sort();

  // filewriter = 1 = this proc writes to file
  //   ping each proc in my cluster, receive its data, write data to file
  // else wait for ping from fileproc, send my data to fileproc

  int tmp,nlines;
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
  for (std::map<int,FnPtrPack>::iterator it=pack_choice.begin(); it!=pack_choice.end(); ++it, ++n) {
      current_pack_choice_key = it->first; // work-around for pack_compute, pack_fix, pack_variable
      (this->*(it->second))(n);
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
  delete [] filecurrent;
  filecurrent = NULL;

  char *filestar = filename;
  if (multiproc) {
    if (multiproc > 1) { // if dump_modify fileper or nfile was used
      delete [] multiname_ex;
      multiname_ex = NULL;
      char *ptr = strchr(filename,'%');
      if (ptr) {
        int id;
        if (me + nclusterprocs == nprocs) // last filewriter
          id = multiproc -1;
        else
          id = me/nclusterprocs;
        multiname_ex = new char[strlen(filename) + 16];
        *ptr = '\0';
        sprintf(multiname_ex,"%s_%d%s",filename,id,ptr+1);
        *ptr = '%';
      }
    } // else multiname_ex built in constructor is OK
    filestar = multiname_ex;
  }

  if (multifile == 0) {
    filecurrent = new char[strlen(filestar) + 1];
    strcpy(filecurrent, filestar);
  } else {
    filecurrent = new char[strlen(filestar) + 16];
    char *ptr = strchr(filestar,'*');
    *ptr = '\0';
    if (padflag == 0) {
      sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
              filestar,update->ntimestep,ptr+1);
    } else {
      char bif[8],pad[16];
      strcpy(bif,BIGINT_FORMAT);
      sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
      sprintf(filecurrent,pad,filestar,update->ntimestep,ptr+1);
    }
    *ptr = '*';
  }

  // filename of domain box data file
  delete [] domainfilecurrent;
  domainfilecurrent = NULL;
  if (multiproc) {
    // remove '%' character
    char *ptr = strchr(filename,'%');
    domainfilecurrent = new char[strlen(filename)];
    *ptr = '\0';
    sprintf(domainfilecurrent,"%s%s",filename,ptr+1);
    *ptr = '%';
    // insert "_boundingBox" string
    ptr = strrchr(domainfilecurrent,'.');
    filestar = new char[strlen(domainfilecurrent)+16];
    *ptr = '\0';
    sprintf(filestar,"%s_boundingBox.%s",domainfilecurrent,ptr+1);
    delete [] domainfilecurrent;
    domainfilecurrent = NULL;

    if (multifile == 0) {
      domainfilecurrent = new char[strlen(filestar) + 1];
      strcpy(domainfilecurrent, filestar);
    } else {
      domainfilecurrent = new char[strlen(filestar) + 16];
      char *ptr = strchr(filestar,'*');
      *ptr = '\0';
      if (padflag == 0) {
        sprintf(domainfilecurrent,"%s" BIGINT_FORMAT "%s",
                filestar,update->ntimestep,ptr+1);
      } else {
        char bif[8],pad[16];
        strcpy(bif,BIGINT_FORMAT);
        sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
        sprintf(domainfilecurrent,pad,filestar,update->ntimestep,ptr+1);
      }
      *ptr = '*';
    }
    delete [] filestar;
    filestar = NULL;
  } else {
    domainfilecurrent = new char[strlen(filecurrent) + 16];
    char *ptr = strrchr(filecurrent,'.');
    *ptr = '\0';
    sprintf(domainfilecurrent,"%s_boundingBox.%s",filecurrent,ptr+1);
    *ptr = '.';
  }

  // filename of parallel file
  if (multiproc && me == 0) {
    delete [] parallelfilecurrent;
    parallelfilecurrent = NULL;

    // remove '%' character and add 'p' to file extension
    // -> string length stays the same
    char *ptr = strchr(filename,'%');
    filestar = new char[strlen(filename) + 1];
    *ptr = '\0';
    sprintf(filestar,"%s%s",filename,ptr+1);
    *ptr = '%';
    ptr = strrchr(filestar,'.');
    ptr++;
    *ptr++='p';
    *ptr++='v';
    *ptr++='t';
    *ptr++= (vtk_file_format == PVTP)?'p':'u';
    *ptr++= 0;

    if (multifile == 0) {
      parallelfilecurrent = new char[strlen(filestar) + 1];
      strcpy(parallelfilecurrent, filestar);
    } else {
      parallelfilecurrent = new char[strlen(filestar) + 16];
      char *ptr = strchr(filestar,'*');
      *ptr = '\0';
      if (padflag == 0) {
        sprintf(parallelfilecurrent,"%s" BIGINT_FORMAT "%s",
                filestar,update->ntimestep,ptr+1);
      } else {
        char bif[8],pad[16];
        strcpy(bif,BIGINT_FORMAT);
        sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
        sprintf(parallelfilecurrent,pad,filestar,update->ntimestep,ptr+1);
      }
      *ptr = '*';
    }
    delete [] filestar;
    filestar = NULL;
  }
}

/* ---------------------------------------------------------------------- */

void DumpVTK::buf2arrays(int n, double *mybuf)
{
  for (int iatom=0; iatom < n; ++iatom) {
    vtkIdType pid[1];
    pid[0] = points->InsertNextPoint(mybuf[iatom*size_one],mybuf[iatom*size_one+1],mybuf[iatom*size_one+2]);

    int j=3; // 0,1,2 = x,y,z handled just above
    for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {
      vtkAbstractArray *paa = it->second;
      if (it->second->GetNumberOfComponents() == 3) {
        switch (vtype[it->first]) {
          case Dump::INT:
            {
              int iv3[3] = { static_cast<int>(mybuf[iatom*size_one+j  ]),
                             static_cast<int>(mybuf[iatom*size_one+j+1]),
                             static_cast<int>(mybuf[iatom*size_one+j+2]) };
              vtkIntArray *pia = static_cast<vtkIntArray*>(paa);
              pia->InsertNextTupleValue(iv3);
              break;
            }
          case Dump::DOUBLE:
            {
              vtkDoubleArray *pda = static_cast<vtkDoubleArray*>(paa);
              pda->InsertNextTupleValue(&mybuf[iatom*size_one+j]);
              break;
            }
        }
        j+=3;
      } else {
        switch (vtype[it->first]) {
          case Dump::INT:
            {
              vtkIntArray *pia = static_cast<vtkIntArray*>(paa);
              pia->InsertNextValue(mybuf[iatom*size_one+j]);
              break;
            }
          case Dump::DOUBLE:
            {
              vtkDoubleArray *pda = static_cast<vtkDoubleArray*>(paa);
              pda->InsertNextValue(mybuf[iatom*size_one+j]);
              break;
            }
          case Dump::STRING:
            {
              vtkStringArray *psa = static_cast<vtkStringArray*>(paa);
              psa->InsertNextValue(typenames[static_cast<int>(mybuf[iatom*size_one+j])]);
              break;
            }
        }
        ++j;
      }
    }

    pointsCells->InsertNextCell(1,pid);
  }
}

/* ---------------------------------------------------------------------- */

void DumpVTK::prepare_domain_data(vtkRectilinearGrid *rgrid)
{
  vtkSmartPointer<vtkDoubleArray> xCoords =  vtkSmartPointer<vtkDoubleArray>::New();
  xCoords->InsertNextValue(boxxlo);
  xCoords->InsertNextValue(boxxhi);
  vtkSmartPointer<vtkDoubleArray> yCoords =  vtkSmartPointer<vtkDoubleArray>::New();
  yCoords->InsertNextValue(boxylo);
  yCoords->InsertNextValue(boxyhi);
  vtkSmartPointer<vtkDoubleArray> zCoords =  vtkSmartPointer<vtkDoubleArray>::New();
  zCoords->InsertNextValue(boxzlo);
  zCoords->InsertNextValue(boxzhi);

  rgrid->SetDimensions(2,2,2);
  rgrid->SetXCoordinates(xCoords);
  rgrid->SetYCoordinates(yCoords);
  rgrid->SetZCoordinates(zCoords);
}

/* ---------------------------------------------------------------------- */

void DumpVTK::prepare_domain_data_triclinic(vtkUnstructuredGrid *hexahedronGrid)
{
  vtkSmartPointer<vtkPoints> hexahedronPoints = vtkSmartPointer<vtkPoints>::New();
  hexahedronPoints->SetNumberOfPoints(8);
  hexahedronPoints->InsertPoint(0, boxcorners[0][0], boxcorners[0][1], boxcorners[0][2]);
  hexahedronPoints->InsertPoint(1, boxcorners[1][0], boxcorners[1][1], boxcorners[1][2]);
  hexahedronPoints->InsertPoint(2, boxcorners[3][0], boxcorners[3][1], boxcorners[3][2]);
  hexahedronPoints->InsertPoint(3, boxcorners[2][0], boxcorners[2][1], boxcorners[2][2]);
  hexahedronPoints->InsertPoint(4, boxcorners[4][0], boxcorners[4][1], boxcorners[4][2]);
  hexahedronPoints->InsertPoint(5, boxcorners[5][0], boxcorners[5][1], boxcorners[5][2]);
  hexahedronPoints->InsertPoint(6, boxcorners[7][0], boxcorners[7][1], boxcorners[7][2]);
  hexahedronPoints->InsertPoint(7, boxcorners[6][0], boxcorners[6][1], boxcorners[6][2]);
  vtkSmartPointer<vtkHexahedron> hexahedron = vtkSmartPointer<vtkHexahedron>::New();
  hexahedron->GetPointIds()->SetId(0, 0);
  hexahedron->GetPointIds()->SetId(1, 1);
  hexahedron->GetPointIds()->SetId(2, 2);
  hexahedron->GetPointIds()->SetId(3, 3);
  hexahedron->GetPointIds()->SetId(4, 4);
  hexahedron->GetPointIds()->SetId(5, 5);
  hexahedron->GetPointIds()->SetId(6, 6);
  hexahedron->GetPointIds()->SetId(7, 7);

  hexahedronGrid->Allocate(1, 1);
  hexahedronGrid->InsertNextCell(hexahedron->GetCellType(),
                                  hexahedron->GetPointIds());
  hexahedronGrid->SetPoints(hexahedronPoints);
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_domain_vtk()
{
  vtkSmartPointer<vtkRectilinearGrid> rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
  prepare_domain_data(rgrid.GetPointer());

  vtkSmartPointer<vtkRectilinearGridWriter> gwriter = vtkSmartPointer<vtkRectilinearGridWriter>::New();

  if(label) gwriter->SetHeader(label);
  else      gwriter->SetHeader("Generated by LAMMPS");

  if (binary) gwriter->SetFileTypeToBinary();
  else        gwriter->SetFileTypeToASCII();

#if VTK_MAJOR_VERSION < 6
  gwriter->SetInput(rgrid);
#else
  gwriter->SetInputData(rgrid);
#endif
  gwriter->SetFileName(domainfilecurrent);
  gwriter->Write();
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_domain_vtk_triclinic()
{
  vtkSmartPointer<vtkUnstructuredGrid> hexahedronGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  prepare_domain_data_triclinic(hexahedronGrid.GetPointer());

  vtkSmartPointer<vtkUnstructuredGridWriter> gwriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();

  if(label) gwriter->SetHeader(label);
  else      gwriter->SetHeader("Generated by LAMMPS");

  if (binary) gwriter->SetFileTypeToBinary();
  else        gwriter->SetFileTypeToASCII();

#if VTK_MAJOR_VERSION < 6
  gwriter->SetInput(hexahedronGrid);
#else
  gwriter->SetInputData(hexahedronGrid);
#endif
  gwriter->SetFileName(domainfilecurrent);
  gwriter->Write();
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_domain_vtr()
{
  vtkSmartPointer<vtkRectilinearGrid> rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
  prepare_domain_data(rgrid.GetPointer());

  vtkSmartPointer<vtkXMLRectilinearGridWriter> gwriter = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();

  if (binary) gwriter->SetDataModeToBinary();
  else        gwriter->SetDataModeToAscii();

#if VTK_MAJOR_VERSION < 6
  gwriter->SetInput(rgrid);
#else
  gwriter->SetInputData(rgrid);
#endif
  gwriter->SetFileName(domainfilecurrent);
  gwriter->Write();
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_domain_vtu_triclinic()
{
  vtkSmartPointer<vtkUnstructuredGrid> hexahedronGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  prepare_domain_data_triclinic(hexahedronGrid.GetPointer());

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> gwriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  if (binary) gwriter->SetDataModeToBinary();
  else        gwriter->SetDataModeToAscii();

#if VTK_MAJOR_VERSION < 6
  gwriter->SetInput(hexahedronGrid);
#else
  gwriter->SetInputData(hexahedronGrid);
#endif
  gwriter->SetFileName(domainfilecurrent);
  gwriter->Write();
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_vtk(int n, double *mybuf)
{
  ++n_calls_;

  buf2arrays(n, mybuf);

  if (n_calls_ < nclusterprocs)
    return; // multiple processors but only proc 0 is a filewriter (-> nclusterprocs procs contribute to the filewriter's output data)

  setFileCurrent();

  {
#ifdef UNSTRUCTURED_GRID_VTK
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_VERTEX, pointsCells);

    for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {
      unstructuredGrid->GetPointData()->AddArray(it->second);
    }

    vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
#else
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetVerts(pointsCells);

    for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {
      polyData->GetPointData()->AddArray(it->second);
    }

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
#endif

    if(label) writer->SetHeader(label);
    else      writer->SetHeader("Generated by LAMMPS");

    if (binary) writer->SetFileTypeToBinary();
    else        writer->SetFileTypeToASCII();

#ifdef UNSTRUCTURED_GRID_VTK
  #if VTK_MAJOR_VERSION < 6
    writer->SetInput(unstructuredGrid);
  #else
    writer->SetInputData(unstructuredGrid);
  #endif
#else
  #if VTK_MAJOR_VERSION < 6
    writer->SetInput(polyData);
  #else
    writer->SetInputData(polyData);
  #endif
#endif
    writer->SetFileName(filecurrent);
    writer->Write();

    if (domain->triclinic == 0)
      write_domain_vtk();
    else
      write_domain_vtk_triclinic();
  }

  reset_vtk_data_containers();
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_vtp(int n, double *mybuf)
{
  ++n_calls_;

  buf2arrays(n, mybuf);

  if (n_calls_ < nclusterprocs)
    return; // multiple processors but not all are filewriters (-> nclusterprocs procs contribute to the filewriter's output data)

  setFileCurrent();

  {
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();

    polyData->SetPoints(points);
    polyData->SetVerts(pointsCells);

    for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {
      polyData->GetPointData()->AddArray(it->second);
    }

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    if (binary) writer->SetDataModeToBinary();
    else        writer->SetDataModeToAscii();

#if VTK_MAJOR_VERSION < 6
    writer->SetInput(polyData);
#else
    writer->SetInputData(polyData);
#endif
    writer->SetFileName(filecurrent);
    writer->Write();

    if (me == 0) {
      if (multiproc) {
        vtkSmartPointer<vtkXMLPPolyDataWriter> pwriter = vtkSmartPointer<vtkXMLPPolyDataWriter>::New();
        pwriter->SetFileName(parallelfilecurrent);
        pwriter->SetNumberOfPieces((multiproc > 1)?multiproc:nprocs);
        if (binary) pwriter->SetDataModeToBinary();
        else        pwriter->SetDataModeToAscii();

#if VTK_MAJOR_VERSION < 6
        pwriter->SetInput(polyData);
#else
        pwriter->SetInputData(polyData);
#endif
        pwriter->Write();
      }

      if (domain->triclinic == 0) {
        domainfilecurrent[strlen(domainfilecurrent)-1] = 'r'; // adjust filename extension
        write_domain_vtr();
      } else {
        domainfilecurrent[strlen(domainfilecurrent)-1] = 'u'; // adjust filename extension
        write_domain_vtu_triclinic();
      }
    }
  }

  reset_vtk_data_containers();
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_vtu(int n, double *mybuf)
{
  ++n_calls_;

  buf2arrays(n, mybuf);

  if (n_calls_ < nclusterprocs)
    return; // multiple processors but not all are filewriters (-> nclusterprocs procs contribute to the filewriter's output data)

  setFileCurrent();

  {
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_VERTEX, pointsCells);

    for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {
      unstructuredGrid->GetPointData()->AddArray(it->second);
    }

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    if (binary) writer->SetDataModeToBinary();
    else        writer->SetDataModeToAscii();

#if VTK_MAJOR_VERSION < 6
    writer->SetInput(unstructuredGrid);
#else
    writer->SetInputData(unstructuredGrid);
#endif
    writer->SetFileName(filecurrent);
    writer->Write();

    if (me == 0) {
      if (multiproc) {
        vtkSmartPointer<vtkXMLPUnstructuredGridWriter> pwriter = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
        pwriter->SetFileName(parallelfilecurrent);
        pwriter->SetNumberOfPieces((multiproc > 1)?multiproc:nprocs);
        if (binary) pwriter->SetDataModeToBinary();
        else        pwriter->SetDataModeToAscii();

#if VTK_MAJOR_VERSION < 6
        pwriter->SetInput(unstructuredGrid);
#else
        pwriter->SetInputData(unstructuredGrid);
#endif
        pwriter->Write();
      }

      if (domain->triclinic == 0) {
        domainfilecurrent[strlen(domainfilecurrent)-1] = 'r'; // adjust filename extension
        write_domain_vtr();
      } else {
        write_domain_vtu_triclinic();
      }
    }
  }

  reset_vtk_data_containers();
}

/* ---------------------------------------------------------------------- */

void DumpVTK::reset_vtk_data_containers()
{
  points = vtkSmartPointer<vtkPoints>::New();
  pointsCells = vtkSmartPointer<vtkCellArray>::New();

  std::map<int,int>::iterator it=vtype.begin();
  ++it; ++it; ++it;
  for (; it!=vtype.end(); ++it) {
    switch(vtype[it->first]) {
      case Dump::INT:
        myarrays[it->first] = vtkSmartPointer<vtkIntArray>::New();
        break;
      case Dump::DOUBLE:
        myarrays[it->first] = vtkSmartPointer<vtkDoubleArray>::New();
        break;
      case Dump::STRING:
        myarrays[it->first] = vtkSmartPointer<vtkStringArray>::New();
        break;
    }

    if (vector_set.find(it->first) != vector_set.end()) {
      myarrays[it->first]->SetNumberOfComponents(3);
      myarrays[it->first]->SetName(name[it->first].c_str());
      ++it; ++it;
    } else {
      myarrays[it->first]->SetName(name[it->first].c_str());
    }
  }
}

/* ---------------------------------------------------------------------- */

int DumpVTK::parse_fields(int narg, char **arg)
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
  int i;
  for (int iarg = 5; iarg < narg; iarg++) {
    i = iarg-5;

    if (strcmp(arg[iarg],"id") == 0) {
      pack_choice[ID] = &DumpVTK::pack_id;
      vtype[ID] = Dump::INT;
      name[ID] = arg[iarg];
    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (!atom->molecule_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
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
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[Q] = &DumpVTK::pack_q;
      vtype[Q] = Dump::DOUBLE;
      name[Q] = arg[iarg];
    } else if (strcmp(arg[iarg],"mux") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[MUX] = &DumpVTK::pack_mux;
      vtype[MUX] = Dump::DOUBLE;
      name[MUX] = arg[iarg];
    } else if (strcmp(arg[iarg],"muy") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[MUY] = &DumpVTK::pack_muy;
      vtype[MUY] = Dump::DOUBLE;
      name[MUY] = arg[iarg];
    } else if (strcmp(arg[iarg],"muz") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[MUZ] = &DumpVTK::pack_muz;
      vtype[MUZ] = Dump::DOUBLE;
      name[MUZ] = arg[iarg];
    } else if (strcmp(arg[iarg],"mu") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[MU] = &DumpVTK::pack_mu;
      vtype[MU] = Dump::DOUBLE;
      name[MU] = arg[iarg];

    } else if (strcmp(arg[iarg],"radius") == 0) {
      if (!atom->radius_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[RADIUS] = &DumpVTK::pack_radius;
      vtype[RADIUS] = Dump::DOUBLE;
      name[RADIUS] = arg[iarg];
    } else if (strcmp(arg[iarg],"diameter") == 0) {
      if (!atom->radius_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[DIAMETER] = &DumpVTK::pack_diameter;
      vtype[DIAMETER] = Dump::DOUBLE;
      name[DIAMETER] = arg[iarg];
    } else if (strcmp(arg[iarg],"omegax") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[OMEGAX] = &DumpVTK::pack_omegax;
      vtype[OMEGAX] = Dump::DOUBLE;
      name[OMEGAX] = arg[iarg];
    } else if (strcmp(arg[iarg],"omegay") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[OMEGAY] = &DumpVTK::pack_omegay;
      vtype[OMEGAY] = Dump::DOUBLE;
      name[OMEGAY] = arg[iarg];
    } else if (strcmp(arg[iarg],"omegaz") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[OMEGAZ] = &DumpVTK::pack_omegaz;
      vtype[OMEGAZ] = Dump::DOUBLE;
      name[OMEGAZ] = arg[iarg];
    } else if (strcmp(arg[iarg],"angmomx") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[ANGMOMX] = &DumpVTK::pack_angmomx;
      vtype[ANGMOMX] = Dump::DOUBLE;
      name[ANGMOMX] = arg[iarg];
    } else if (strcmp(arg[iarg],"angmomy") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[ANGMOMY] = &DumpVTK::pack_angmomy;
      vtype[ANGMOMY] = Dump::DOUBLE;
      name[ANGMOMY] = arg[iarg];
    } else if (strcmp(arg[iarg],"angmomz") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[ANGMOMZ] = &DumpVTK::pack_angmomz;
      vtype[ANGMOMZ] = Dump::DOUBLE;
      name[ANGMOMZ] = arg[iarg];
    } else if (strcmp(arg[iarg],"tqx") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[TQX] = &DumpVTK::pack_tqx;
      vtype[TQX] = Dump::DOUBLE;
      name[TQX] = arg[iarg];
    } else if (strcmp(arg[iarg],"tqy") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[TQY] = &DumpVTK::pack_tqy;
      vtype[TQY] = Dump::DOUBLE;
      name[TQY] = arg[iarg];
    } else if (strcmp(arg[iarg],"tqz") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[TQZ] = &DumpVTK::pack_tqz;
      vtype[TQZ] = Dump::DOUBLE;
      name[TQZ] = arg[iarg];

    // compute value = c_ID
    // if no trailing [], then arg is set to 0, else arg is int between []

    } else if (strncmp(arg[iarg],"c_",2) == 0) {
      pack_choice[ATTRIBUTES+i] = &DumpVTK::pack_compute;
      vtype[ATTRIBUTES+i] = Dump::DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump vtk command");
        argindex[ATTRIBUTES+i] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[ATTRIBUTES+i] = 0;

      n = modify->find_compute(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump vtk compute ID");
      if (modify->compute[n]->peratom_flag == 0)
        error->all(FLERR,"Dump vtk compute does not compute per-atom info");
      if (argindex[ATTRIBUTES+i] == 0 && modify->compute[n]->size_peratom_cols > 0)
        error->all(FLERR,
                   "Dump vtk compute does not calculate per-atom vector");
      if (argindex[ATTRIBUTES+i] > 0 && modify->compute[n]->size_peratom_cols == 0)
        error->all(FLERR,\
                   "Dump vtk compute does not calculate per-atom array");
      if (argindex[ATTRIBUTES+i] > 0 &&
          argindex[ATTRIBUTES+i] > modify->compute[n]->size_peratom_cols)
        error->all(FLERR,"Dump vtk compute vector is accessed out-of-range");

      field2index[ATTRIBUTES+i] = add_compute(suffix);
      name[ATTRIBUTES+i] = arg[iarg];
      delete [] suffix;

    // fix value = f_ID
    // if no trailing [], then arg is set to 0, else arg is between []

    } else if (strncmp(arg[iarg],"f_",2) == 0) {
      pack_choice[ATTRIBUTES+i] = &DumpVTK::pack_fix;
      vtype[ATTRIBUTES+i] = Dump::DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump vtk command");
        argindex[ATTRIBUTES+i] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[ATTRIBUTES+i] = 0;

      n = modify->find_fix(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump vtk fix ID");
      if (modify->fix[n]->peratom_flag == 0)
        error->all(FLERR,"Dump vtk fix does not compute per-atom info");
      if (argindex[ATTRIBUTES+i] == 0 && modify->fix[n]->size_peratom_cols > 0)
        error->all(FLERR,"Dump vtk fix does not compute per-atom vector");
      if (argindex[ATTRIBUTES+i] > 0 && modify->fix[n]->size_peratom_cols == 0)
        error->all(FLERR,"Dump vtk fix does not compute per-atom array");
      if (argindex[ATTRIBUTES+i] > 0 &&
          argindex[ATTRIBUTES+i] > modify->fix[n]->size_peratom_cols)
        error->all(FLERR,"Dump vtk fix vector is accessed out-of-range");

      field2index[ATTRIBUTES+i] = add_fix(suffix);
      name[ATTRIBUTES+i] = arg[iarg];
      delete [] suffix;

    // variable value = v_name

    } else if (strncmp(arg[iarg],"v_",2) == 0) {
      pack_choice[ATTRIBUTES+i] = &DumpVTK::pack_variable;
      vtype[ATTRIBUTES+i] = Dump::DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      argindex[ATTRIBUTES+i] = 0;

      n = input->variable->find(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump vtk variable name");
      if (input->variable->atomstyle(n) == 0)
        error->all(FLERR,"Dump vtk variable is not atom-style variable");

      field2index[ATTRIBUTES+i] = add_variable(suffix);
      name[ATTRIBUTES+i] = suffix;
      delete [] suffix;

    // custom per-atom floating point value = d_ID

    } else if (strncmp(arg[iarg],"d_",2) == 0) {
      pack_choice[ATTRIBUTES+i] = &DumpVTK::pack_custom;
      vtype[ATTRIBUTES+i] = Dump::DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);
      argindex[ATTRIBUTES+i] = 0;

      int tmp = -1;
      n = atom->find_custom(suffix,tmp);
      if (n < 0)
        error->all(FLERR,"Could not find custom per-atom property ID");

      if (tmp != 1)
        error->all(FLERR,"Custom per-atom property ID is not floating point");

      field2index[ATTRIBUTES+i] = add_custom(suffix,1);
      name[ATTRIBUTES+i] = suffix;
      delete [] suffix;

    // custom per-atom integer value = i_ID

    } else if (strncmp(arg[iarg],"i_",2) == 0) {
      pack_choice[ATTRIBUTES+i] = &DumpVTK::pack_custom;
      vtype[ATTRIBUTES+i] = Dump::INT;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);
      argindex[ATTRIBUTES+i] = 0;

      int tmp = -1;
      n = atom->find_custom(suffix,tmp);
      if (n < 0)
        error->all(FLERR,"Could not find custom per-atom property ID");

      if (tmp != 0)
        error->all(FLERR,"Custom per-atom property ID is not integer");

      field2index[ATTRIBUTES+i] = add_custom(suffix,0);
      name[ATTRIBUTES+i] = suffix;
      delete [] suffix;

    } else return iarg;
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
    if(name.count(vector3_starts[v3s]  ) &&
       name.count(vector3_starts[v3s]+1) &&
       name.count(vector3_starts[v3s]+2) )
    {
      std::string vectorName = name[vector3_starts[v3s]];
      vectorName.erase(vectorName.find_first_of('x'));
      name[vector3_starts[v3s]] = vectorName;
      vector_set.insert(vector3_starts[v3s]);
    }
  }

  // compute and fix vectors
  for (std::map<int,std::string>::iterator it=name.begin(); it!=name.end(); ++it) {
    if (it->first < ATTRIBUTES) // neither fix nor compute
      continue;

    if(argindex[it->first] == 0) // single value
      continue;

    // assume components are grouped together and in correct order
    if(name.count(it->first + 1) && name.count(it->first + 2) ) { // more attributes?
      if(it->second.compare(0,it->second.length()-3,name[it->first + 1],0,it->second.length()-3) == 0  && // same attributes?
         it->second.compare(0,it->second.length()-3,name[it->first + 2],0,it->second.length()-3) == 0 )
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

/* ----------------------------------------------------------------------
   add Compute to list of Compute objects used by dump
   return index of where this Compute is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpVTK::add_compute(char *id)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(id,id_compute[icompute]) == 0) break;
  if (icompute < ncompute) return icompute;

  id_compute = (char **)
    memory->srealloc(id_compute,(ncompute+1)*sizeof(char *),"dump:id_compute");
  delete [] compute;
  compute = new Compute*[ncompute+1];

  int n = strlen(id) + 1;
  id_compute[ncompute] = new char[n];
  strcpy(id_compute[ncompute],id);
  ncompute++;
  return ncompute-1;
}

/* ----------------------------------------------------------------------
   add Fix to list of Fix objects used by dump
   return index of where this Fix is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpVTK::add_fix(char *id)
{
  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(id,id_fix[ifix]) == 0) break;
  if (ifix < nfix) return ifix;

  id_fix = (char **)
    memory->srealloc(id_fix,(nfix+1)*sizeof(char *),"dump:id_fix");
  delete [] fix;
  fix = new Fix*[nfix+1];

  int n = strlen(id) + 1;
  id_fix[nfix] = new char[n];
  strcpy(id_fix[nfix],id);
  nfix++;
  return nfix-1;
}

/* ----------------------------------------------------------------------
   add Variable to list of Variables used by dump
   return index of where this Variable is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpVTK::add_variable(char *id)
{
  int ivariable;
  for (ivariable = 0; ivariable < nvariable; ivariable++)
    if (strcmp(id,id_variable[ivariable]) == 0) break;
  if (ivariable < nvariable) return ivariable;

  id_variable = (char **)
    memory->srealloc(id_variable,(nvariable+1)*sizeof(char *),
                     "dump:id_variable");
  delete [] variable;
  variable = new int[nvariable+1];
  delete [] vbuf;
  vbuf = new double*[nvariable+1];
  for (int i = 0; i <= nvariable; i++) vbuf[i] = NULL;

  int n = strlen(id) + 1;
  id_variable[nvariable] = new char[n];
  strcpy(id_variable[nvariable],id);
  nvariable++;
  return nvariable-1;
}

/* ----------------------------------------------------------------------
   add custom atom property to list used by dump
   return index of where this property is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpVTK::add_custom(char *id, int flag)
{
  int icustom;
  for (icustom = 0; icustom < ncustom; icustom++)
    if ((strcmp(id,id_custom[icustom]) == 0)
        && (flag == flag_custom[icustom])) break;
  if (icustom < ncustom) return icustom;

  id_custom = (char **)
    memory->srealloc(id_custom,(ncustom+1)*sizeof(char *),"dump:id_custom");
  flag_custom = (int *)
    memory->srealloc(flag_custom,(ncustom+1)*sizeof(int),"dump:flag_custom");

  int n = strlen(id) + 1;
  id_custom[ncustom] = new char[n];
  strcpy(id_custom[ncustom],id);
  flag_custom[ncustom] = flag;

  ncustom++;
  return ncustom-1;
}

/* ---------------------------------------------------------------------- */

int DumpVTK::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"region") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) iregion = -1;
    else {
      iregion = domain->find_region(arg[1]);
      if (iregion == -1)
        error->all(FLERR,"Dump_modify region ID does not exist");
      delete [] idregion;
      int n = strlen(arg[1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[1]);
    }
    return 2;
  }

  if (strcmp(arg[0],"label") == 0) {
     if (narg < 2) error->all(FLERR,"Illegal dump_modify command [label]");
     delete [] label;
     int n = strlen(arg[1]) + 1;
     label = new char[n];
     strcpy(label,arg[1]);
     return 2;
   }

  if (strcmp(arg[0],"binary") == 0) {
     if (narg < 2) error->all(FLERR,"Illegal dump_modify command [binary]");
     if (strcmp(arg[1],"yes") == 0) binary = 1;
     else if (strcmp(arg[1],"no") == 0) binary = 0;
     else error->all(FLERR,"Illegal dump_modify command [binary]");
     return 2;
  }

  if (strcmp(arg[0],"element") == 0) {
    if (narg < ntypes+1)
      error->all(FLERR,"Dump modify: number of element names do not match atom types");

    if (typenames) {
      for (int i = 1; i <= ntypes; i++) delete [] typenames[i];
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

  if (strcmp(arg[0],"thresh") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) {
      if (nthresh) {
        memory->destroy(thresh_array);
        memory->destroy(thresh_op);
        memory->destroy(thresh_value);
        thresh_array = NULL;
        thresh_op = NULL;
        thresh_value = NULL;
      }
      nthresh = 0;
      return 2;
    }

    if (narg < 4) error->all(FLERR,"Illegal dump_modify command");

    // grow threshold arrays

    memory->grow(thresh_array,nthresh+1,"dump:thresh_array");
    memory->grow(thresh_op,(nthresh+1),"dump:thresh_op");
    memory->grow(thresh_value,(nthresh+1),"dump:thresh_value");

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

    // compute value = c_ID
    // if no trailing [], then arg is set to 0, else arg is between []

    else if (strncmp(arg[1],"c_",2) == 0) {
      thresh_array[nthresh] = COMPUTE;
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump modify command");
        argindex[ATTRIBUTES+nfield+nthresh] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[ATTRIBUTES+nfield+nthresh] = 0;

      n = modify->find_compute(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump modify compute ID");

      if (modify->compute[n]->peratom_flag == 0)
        error->all(FLERR,
                   "Dump modify compute ID does not compute per-atom info");
      if (argindex[ATTRIBUTES+nfield+nthresh] == 0 &&
          modify->compute[n]->size_peratom_cols > 0)
        error->all(FLERR,
                   "Dump modify compute ID does not compute per-atom vector");
      if (argindex[ATTRIBUTES+nfield+nthresh] > 0 &&
          modify->compute[n]->size_peratom_cols == 0)
        error->all(FLERR,
                   "Dump modify compute ID does not compute per-atom array");
      if (argindex[ATTRIBUTES+nfield+nthresh] > 0 &&
          argindex[ATTRIBUTES+nfield+nthresh] > modify->compute[n]->size_peratom_cols)
        error->all(FLERR,"Dump modify compute ID vector is not large enough");

      field2index[ATTRIBUTES+nfield+nthresh] = add_compute(suffix);
      delete [] suffix;

    // fix value = f_ID
    // if no trailing [], then arg is set to 0, else arg is between []

    } else if (strncmp(arg[1],"f_",2) == 0) {
      thresh_array[nthresh] = FIX;
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump modify command");
        argindex[ATTRIBUTES+nfield+nthresh] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[ATTRIBUTES+nfield+nthresh] = 0;

      n = modify->find_fix(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump modify fix ID");

      if (modify->fix[n]->peratom_flag == 0)
        error->all(FLERR,"Dump modify fix ID does not compute per-atom info");
      if (argindex[ATTRIBUTES+nfield+nthresh] == 0 &&
          modify->fix[n]->size_peratom_cols > 0)
        error->all(FLERR,"Dump modify fix ID does not compute per-atom vector");
      if (argindex[ATTRIBUTES+nfield+nthresh] > 0 &&
          modify->fix[n]->size_peratom_cols == 0)
        error->all(FLERR,"Dump modify fix ID does not compute per-atom array");
      if (argindex[ATTRIBUTES+nfield+nthresh] > 0 &&
          argindex[ATTRIBUTES+nfield+nthresh] > modify->fix[n]->size_peratom_cols)
        error->all(FLERR,"Dump modify fix ID vector is not large enough");

      field2index[ATTRIBUTES+nfield+nthresh] = add_fix(suffix);
      delete [] suffix;

    // variable value = v_ID

    } else if (strncmp(arg[1],"v_",2) == 0) {
      thresh_array[nthresh] = VARIABLE;
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);

      argindex[ATTRIBUTES+nfield+nthresh] = 0;

      n = input->variable->find(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump modify variable name");
      if (input->variable->atomstyle(n) == 0)
        error->all(FLERR,"Dump modify variable is not atom-style variable");

      field2index[ATTRIBUTES+nfield+nthresh] = add_variable(suffix);
      delete [] suffix;

    } else error->all(FLERR,"Invalid dump_modify threshold operator");

    // set operation type of threshold

    if (strcmp(arg[2],"<") == 0) thresh_op[nthresh] = LT;
    else if (strcmp(arg[2],"<=") == 0) thresh_op[nthresh] = LE;
    else if (strcmp(arg[2],">") == 0) thresh_op[nthresh] = GT;
    else if (strcmp(arg[2],">=") == 0) thresh_op[nthresh] = GE;
    else if (strcmp(arg[2],"==") == 0) thresh_op[nthresh] = EQ;
    else if (strcmp(arg[2],"!=") == 0) thresh_op[nthresh] = NEQ;
    else error->all(FLERR,"Invalid dump_modify threshold operator");

    // set threshold value

    thresh_value[nthresh] = force->numeric(FLERR,arg[3]);

    nthresh++;
    return 4;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory in buf, choose, variable arrays
------------------------------------------------------------------------- */

bigint DumpVTK::memory_usage()
{
  bigint bytes = Dump::memory_usage();
  bytes += memory->usage(choose,maxlocal);
  bytes += memory->usage(dchoose,maxlocal);
  bytes += memory->usage(clist,maxlocal);
  bytes += memory->usage(vbuf,nvariable,maxlocal);
  return bytes;
}

/* ----------------------------------------------------------------------
   extraction of Compute, Fix, Variable results
------------------------------------------------------------------------- */

void DumpVTK::pack_compute(int n)
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

void DumpVTK::pack_fix(int n)
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

void DumpVTK::pack_variable(int n)
{
  double *vector = vbuf[field2index[current_pack_choice_key]];

  for (int i = 0; i < nchoose; i++) {
    buf[n] = vector[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpVTK::pack_custom(int n)
{
  int index = field2index[n];

  if (flag_custom[index] == 0) { // integer
    int iwhich,tmp;
    iwhich = atom->find_custom(id_custom[index],tmp);

    int *ivector = atom->ivector[iwhich];
    for (int i = 0; i < nchoose; i++) {
      buf[n] = ivector[clist[i]];
      n += size_one;
    }
  } else if (flag_custom[index] == 1) { // double
    int iwhich,tmp;
    iwhich = atom->find_custom(id_custom[index],tmp);

    double *dvector = atom->dvector[iwhich];
    for (int i = 0; i < nchoose; i++) {
      buf[n] = dvector[clist[i]];
      n += size_one;
    }
  }
}
