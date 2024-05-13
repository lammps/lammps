// Header file for this class
#include "LammpsInterface.h"
// LAMMPS includes
#include "lammps.h"
#include "atom.h" // x, v, f
#include "atom_vec.h" //for insertion
#include "domain.h" // for basing locations on regions
#include "region.h" // region bounding box and style
#include "force.h" // boltzman constant
#include "group.h" // atom masks
#include "memory.h" // grow atom information
#include "compute.h" // computes
#include "compute_pe_atom.h" // computes potential energy per atom
#include "compute_stress_atom.h" // computes stress per atom
#include "compute_centro_atom.h" // computes centrosymmetry per atom
#include "compute_cna_atom.h" // computes common-neighbor-analysis per atom
#include "compute_coord_atom.h" // computes coordination number per atom
#include "compute_ke_atom.h" // computes kinetic energy per atom
#include "modify.h" //
#include "neighbor.h" // neighbors
#include "neigh_list.h" // neighbor list
#include "update.h" // timestepping information
#include "pair.h" // pair potentials
#include "MANYBODY/pair_eam.h" // pair potentials
#include "lattice.h" // lattice parameters
#include "bond.h" // bond potentials
#include "comm.h" //
#include "fix.h"
#include "utils.h"

// ATC includes
#include "ATC_Error.h"
#include "MatrixLibrary.h"
#include "Utility.h"
using ATC_Utility::to_string;

// Other include files
#include "mpi.h"
#include <cstring>
#include <map>
#include <typeinfo>

using std::max;
using std::stringstream;
using std::copy;
using std::map;
using std::pair;
using std::string;
using std::set;
using LAMMPS_NS::bigint;
using LAMMPS_NS::utils::read_lines_from_file;

namespace ATC {

const static double PI = 3.141592653589793238;

const static int seed_ = 3141592;

const static int MAX_GROUP_BIT = 2147483647; //4294967295; // pow(2,31)-1;

double norm(double * v) {return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); }

LammpsInterface * LammpsInterface::myInstance_ = nullptr;

// -----------------------------------------------------------------
//  instance()
// -----------------------------------------------------------------
LammpsInterface * LammpsInterface::instance()
{
  if (myInstance_ == nullptr) {
    myInstance_ = new LammpsInterface();
  }
  return myInstance_;
}

// -----------------------------------------------------------------
//  Destroy()
// -----------------------------------------------------------------
void LammpsInterface::Destroy()
{
  if (myInstance_) delete myInstance_;
  myInstance_ = nullptr;
}


// -----------------------------------------------------------------
//  constructor
// -----------------------------------------------------------------
LammpsInterface::LammpsInterface()
  : lammps_(nullptr),
    fixPointer_(nullptr),
    commRank_(0),
    atomPE_(nullptr),
    refBoxIsSet_(false),
    random_(nullptr),
    globalrandom_(nullptr)
{
}

// -----------------------------------------------------------------
//  general interface methods
// -----------------------------------------------------------------

MPI_Comm LammpsInterface::world() const { return lammps_->world; }
void LammpsInterface::set_fix_pointer(LAMMPS_NS::Fix * thisFix) { fixPointer_ = thisFix; }
void LammpsInterface::forward_comm_fix() const { lammps_->comm->forward_comm(fixPointer_); }
void LammpsInterface::comm_borders() const { lammps_->comm->borders(); }

#ifndef ISOLATE_FE
void LammpsInterface::sparse_allsum(SparseMatrix<double> &toShare) const
{
  toShare.compress();

  // initialize MPI information
  int nProcs;
  int myRank;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  int error;

  // get numbers of rows, columns, rowsCRS, and
  // sizes (number of nonzero elements in matrix)
  SparseMatInfo *recInfo = new SparseMatInfo[nProcs];
  SparseMatInfo myInfo;
  myInfo.rows    = toShare.nRows();
  myInfo.cols    = toShare.nCols();
  myInfo.rowsCRS = toShare.nRowsCRS();
  myInfo.size    = toShare.size();

  error = MPI_Allgather(&myInfo, 4, MPI_INT,
                        recInfo, 4, MPI_INT, lammps_->world);
  if (error != MPI_SUCCESS) throw ATC_Error("error in sparse_allsum_numrows "+to_string(error));

  // adjust row sendcounts because recRowsCRS is off by one
  int *rowCounts = new int[nProcs];
  int *sizeCounts = new int[nProcs];
  // set up total size of receive buffers for Allgatherv calls
  int totalRowsCRS = 0;
  int totalSize = 0;
  // set up array of displacements for Allgatherv calls
  int *rowOffsets = new int[nProcs];
  rowOffsets[0] = 0;
  int *sizeOffsets = new int[nProcs];
  sizeOffsets[0] = 0;
  for (int i = 0; i < nProcs; i++) {
    // find the total number of entries to share in the mpi calls below
    rowCounts[i] = recInfo[i].rowsCRS + 1;
    sizeCounts[i] = recInfo[i].size;
    totalRowsCRS += rowCounts[i];
    totalSize += recInfo[i].size;
    // these already have their 0th slot filled in
    if (i == 0) continue;
    rowOffsets[i] = rowOffsets[i-1] + rowCounts[i-1];
    sizeOffsets[i] = sizeOffsets[i-1] + sizeCounts[i-1];
  }

  // get actual rows
  INDEX *rec_ia = new INDEX[totalRowsCRS];
  if (toShare.size() == 0) {
    double dummy;
    error = MPI_Allgatherv(&dummy, 0, MPI_INT,
                           rec_ia, rowCounts, rowOffsets, MPI_INT, lammps_->world);
  }
  else
    error = MPI_Allgatherv(toShare.rows(), rowCounts[myRank], MPI_INT,
                           rec_ia, rowCounts, rowOffsets, MPI_INT, lammps_->world);
  if (error != MPI_SUCCESS)
    throw ATC_Error("error in sparse_allsum_rowarray "+to_string(error));

  // get actual cols
  INDEX *rec_ja = new INDEX[totalSize];
  error = MPI_Allgatherv(toShare.cols(), sizeCounts[myRank], MPI_INT,
                         rec_ja, sizeCounts, sizeOffsets, MPI_INT, lammps_->world);
  if (error != MPI_SUCCESS)
    throw ATC_Error("error in sparse_allsum_colarray "+to_string(error));

  // get the array of values
  double *rec_vals = new double[totalSize];
  error = MPI_Allgatherv(toShare.ptr(), sizeCounts[myRank], MPI_DOUBLE,
                         rec_vals, sizeCounts, sizeOffsets, MPI_DOUBLE, lammps_->world);
  if (error != MPI_SUCCESS)
    throw ATC_Error("error in sparse_allsum_valarray "+to_string(error));

  INDEX *rec_ia_proc;
  INDEX *rec_ja_proc;
  double *rec_vals_proc;
  for (int i = 0; i < nProcs; i++) {
    if (myRank != i) {
      // deallocated when tempMat is deleted since it wraps them
      rec_ia_proc = new INDEX[rowCounts[i]];
      rec_ja_proc = new INDEX[sizeCounts[i]];
      rec_vals_proc = new double[sizeCounts[i]];

      // copy the data passed with MPI into the new spots
      copy(rec_ia + rowOffsets[i],
           rec_ia + rowOffsets[i] + rowCounts[i],
           rec_ia_proc);
      copy(rec_ja + sizeOffsets[i],
           rec_ja + sizeOffsets[i] + sizeCounts[i],
           rec_ja_proc);
      copy(rec_vals + sizeOffsets[i],
           rec_vals + sizeOffsets[i] + sizeCounts[i],
           rec_vals_proc);

      // Does anyone know why we have to declare tempMat here (as well as set it equal to
      // something) to avoid segfaults? there are still segfaults, but they happen at a much
      // later stage of the game now (and for less benchmarks overall).
      SparseMatrix<double> tempMat =
        SparseMatrix<double>(rec_ia_proc, rec_ja_proc, rec_vals_proc,
                             recInfo[i].size, recInfo[i].rows,
                             recInfo[i].cols, recInfo[i].rowsCRS);
      toShare += tempMat;
    }
  }
  delete[] rowCounts;
  delete[] sizeCounts;
  delete[] rowOffsets;
  delete[] sizeOffsets;

  delete[] recInfo;
  delete[] rec_ia;
  delete[] rec_ja;
  delete[] rec_vals;
}
#endif

std::string LammpsInterface::read_file(std::string filename) const
{
  FILE *fp = nullptr;
  if (! comm_rank()) {
    fp = fopen(filename.c_str(),"r");
    if (!fp) throw ATC_Error("can't open file: "+filename);
  }
  const int MAXLINE = 256;
  const int CHUNK   = 1024;
  char * buffer = new char[CHUNK*MAXLINE];
  std::stringstream s;
  bool eof = false;
  while ( ! eof) {
    eof = read_lines_from_file(fp,1,MAXLINE,buffer,comm_rank(),lammps_->world);
    s << buffer;
  }
  fclose(fp);
  delete [] buffer;
  return s.str(); // can't copy the stream itself
}


// -----------------------------------------------------------------
//  atom interface methods
// -----------------------------------------------------------------
string LammpsInterface::fix_id() const { return string(fixPointer_->id); }

int LammpsInterface::nlocal() const { return lammps_->atom->nlocal; }

int LammpsInterface::nghost() const { return lammps_->atom->nghost; }

bool LammpsInterface::atoms_sorted() const
{
  int sortfreq = lammps_->atom->sortfreq;
  if (sortfreq > 0) { return true; }
  else              { return false;
}
}


bigint LammpsInterface::natoms() const { return lammps_->atom->natoms; }

int LammpsInterface::nmax() const { return lammps_->atom->nmax; }

int LammpsInterface::ntypes() const { return lammps_->atom->ntypes; }

double ** LammpsInterface::xatom() const { return lammps_->atom->x; }

int LammpsInterface::type_to_charge(int atype) const {
  double *q = lammps_->atom->q;
  if (! q) return 0;
  int nlocal = lammps_->atom->nlocal;
  int *type = lammps_->atom->type;
  double aq = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] == atype) {
      aq = q[i];
      break;
    }
  }
  double pcharge;
  MPI_Allreduce(&aq,&pcharge,1,MPI_DOUBLE,MPI_MAX,world());
  double ncharge;
  MPI_Allreduce(&aq,&ncharge,1,MPI_DOUBLE,MPI_MIN,world());
  double charge = (pcharge == 0.0) ? ncharge : pcharge;
  return charge;
}

//const double ** LammpsInterface::xatom() const { return (const double**)(lammps_->atom->x); }

double ** LammpsInterface::vatom() const { return lammps_->atom->v; }

double ** LammpsInterface::fatom() const { return lammps_->atom->f; }

const int * LammpsInterface::atom_mask() const { return (const int*)lammps_->atom->mask; }

int * LammpsInterface::atom_mask() { return lammps_->atom->mask; }

int * LammpsInterface::atom_type() const { return lammps_->atom->type; }

int * LammpsInterface::atom_tag() const { return lammps_->atom->tag; }

int * LammpsInterface::atom_to_molecule() const { return lammps_->atom->molecule; }


int * LammpsInterface::num_bond() const { return lammps_->atom->num_bond; }

int ** LammpsInterface::bond_atom() const { return lammps_->atom->bond_atom; }

int * LammpsInterface::image() const { return lammps_->atom->image; }

int LammpsInterface::bond_per_atom() const { return lammps_->atom->bond_per_atom; }

int LammpsInterface::newton_bond() const { return lammps_->force->newton_bond; }

int LammpsInterface::local_to_global_map(int global) const { return lammps_->atom->map(global); }

double * LammpsInterface::atom_mass() const { return lammps_->atom->mass; }

double   LammpsInterface::atom_mass(int iType) const { return lammps_->atom->mass[iType]; }

double * LammpsInterface::atom_rmass() const { return lammps_->atom->rmass; }

double * LammpsInterface::atom_charge() const { return lammps_->atom->q; }

double * LammpsInterface::atom_scalar(FundamentalAtomQuantity quantityType) const
{
  if (quantityType==ATOM_MASS) {
    if (atom_mass())
      throw ATC_Error("Atom mass array requested but not defined");
    return atom_rmass();
  }
  else if (quantityType==ATOM_CHARGE) {
    double * atomCharge = atom_charge();
    if (!atomCharge)
      throw ATC_Error("Atom charge array requested but not defined");
    return atomCharge;
  }
  else
    throw ATC_Error("BAD type requested in atom_scalar");
  return nullptr;
}

double ** LammpsInterface::atom_vector(FundamentalAtomQuantity quantityType) const
{
  if (quantityType==ATOM_POSITION)
    return xatom();
  else if (quantityType==ATOM_VELOCITY)
    return vatom();
  else if (quantityType==ATOM_FORCE)
    return fatom();
  else
    throw ATC_Error("BAD type requested in atom_vector");
  return nullptr;
}

int LammpsInterface::atom_quantity_ndof(FundamentalAtomQuantity quantityType) const
{
  if (quantityType==ATOM_MASS || quantityType==ATOM_CHARGE)
    return 1;
  else if (quantityType==ATOM_POSITION || quantityType==ATOM_VELOCITY || quantityType==ATOM_FORCE)
    return dimension();
  else
    throw ATC_Error("BAD type requested in atom_quantity_ndof");
}

double LammpsInterface::atom_quantity_conversion(FundamentalAtomQuantity quantityType) const
{
  if (quantityType==ATOM_MASS || quantityType==ATOM_CHARGE || quantityType==ATOM_POSITION || quantityType==ATOM_VELOCITY)
    return 1;
  else if ( quantityType==ATOM_FORCE)
    return ftm2v();
  else
    throw ATC_Error("BAD type requested in atom_quantity_conversion");
}

// -----------------------------------------------------------------
//  domain interface methods
// -----------------------------------------------------------------

int LammpsInterface::dimension() const { return lammps_->domain->dimension; }

int LammpsInterface::nregion() const { return lammps_->domain->get_region_list().size(); }

void LammpsInterface::box_bounds(double & boxxlo, double & boxxhi,
                                 double & boxylo, double & boxyhi,
                                 double & boxzlo, double &boxzhi) const
{
  if (lammps_->domain->triclinic == 0) {
    boxxlo = lammps_->domain->boxlo[0];
    boxxhi = lammps_->domain->boxhi[0];
    boxylo = lammps_->domain->boxlo[1];
    boxyhi = lammps_->domain->boxhi[1];
    boxzlo = lammps_->domain->boxlo[2];
    boxzhi = lammps_->domain->boxhi[2];
  }
  else {
    boxxlo = lammps_->domain->boxlo_bound[0];
    boxxhi = lammps_->domain->boxhi_bound[0];
    boxylo = lammps_->domain->boxlo_bound[1];
    boxyhi = lammps_->domain->boxhi_bound[1];
    boxzlo = lammps_->domain->boxlo_bound[2];
    boxzhi = lammps_->domain->boxhi_bound[2];
  }
}

bool LammpsInterface::in_box(double * x) const
{
  double xlo,xhi,ylo,yhi,zlo,zhi;
  box_bounds(xlo,xhi,ylo,yhi,zlo,zhi);
  if (x[0] >= xlo && x[0] < xhi &&
      x[1] >= ylo && x[1] < yhi &&
      x[2] >= zlo && x[2] < zhi)
  return true;
  return false;
}

bool LammpsInterface::in_my_processor_box(double * x) const
{
  if (x[0] >= lammps_->domain->sublo[0] && x[0] < lammps_->domain->subhi[0] &&
      x[1] >= lammps_->domain->sublo[1] && x[1] < lammps_->domain->subhi[1] &&
      x[2] >= lammps_->domain->sublo[2] && x[2] < lammps_->domain->subhi[2])
  return true;
  if (! in_box(x))
    throw ATC_Error("point is in no processors box");
  return false;
}


void LammpsInterface::sub_bounds(double & subxlo, double & subxhi,
                                 double & subylo, double & subyhi,
                                 double & subzlo, double & subzhi) const
{
  if (lammps_->domain->triclinic == 0) {
    subxlo = lammps_->domain->sublo[0];
    subxhi = lammps_->domain->subhi[0];
    subylo = lammps_->domain->sublo[1];
    subyhi = lammps_->domain->subhi[1];
    subzlo = lammps_->domain->sublo[2];
    subzhi = lammps_->domain->subhi[2];
  }
  else {
    ATC_Error("Subboxes not accurate when triclinic != 0.");
  }
}

int LammpsInterface::xperiodic() const { return lammps_->domain->xperiodic; }

int LammpsInterface::yperiodic() const { return lammps_->domain->yperiodic; }

int LammpsInterface::zperiodic() const { return lammps_->domain->zperiodic; }

int LammpsInterface::nperiodic() const
{
  int nprd = 0;
  if ( lammps_->domain->xperiodic > 0 ) { nprd++ ; }
  if ( lammps_->domain->yperiodic > 0 ) { nprd++ ; }
  if ( lammps_->domain->zperiodic > 0 ) { nprd++ ; }
  return nprd;
}

// correct posistions for periodic box
void LammpsInterface::periodicity_correction(double * x) const
{
  int* periodicity = lammps_->domain->periodicity;
  if (!refBoxIsSet_) set_reference_box();
  for (int m = 0; m < 3; m++) {
    if ((bool) periodicity[m]) {
      if (x[m] < lower_[m] || x[m] > upper_[m]) {
        x[m] -= length_[m]*floor((x[m]-lower_[m])/length_[m]);
      }
      if (x[m] < lower_[m] || x[m] > upper_[m]) {
        throw ATC_Error("periodicity_correction: still out of box bounds");
      }
    }
  }
}

void LammpsInterface::set_reference_box() const
{
  double * hi = lammps_->domain->boxhi;
  double * lo = lammps_->domain->boxlo;
  double * len = lammps_->domain->prd;
  for (int i = 0; i < 3; i++) {
    upper_[i] = hi[i];
    lower_[i] = lo[i];
    length_[i] = len[i];
  }
  refBoxIsSet_ = true;
}


double LammpsInterface::domain_xprd() const { return lammps_->domain->xprd; }

double LammpsInterface::domain_yprd() const { return lammps_->domain->yprd; }

double LammpsInterface::domain_zprd() const { return lammps_->domain->zprd; }

double LammpsInterface::domain_volume() const
{
  return (lammps_->domain->xprd)*
         (lammps_->domain->yprd)*
         (lammps_->domain->zprd);
}

double LammpsInterface::domain_xy() const { return lammps_->domain->xy; }

double LammpsInterface::domain_xz() const { return lammps_->domain->xz; }

double LammpsInterface::domain_yz() const { return lammps_->domain->yz; }

int LammpsInterface::domain_triclinic() const { return lammps_->domain->triclinic; }

void LammpsInterface::box_periodicity(int & xperiodic,
                                      int & yperiodic,
                                      int & zperiodic) const
{
  xperiodic = lammps_->domain->xperiodic;
  yperiodic = lammps_->domain->yperiodic;
  zperiodic = lammps_->domain->zperiodic;
}

int LammpsInterface::region_id(const char *regionName) const {
  auto regions = lammps_->domain->get_region_list();
  int iregion = 0;
  for (auto reg : regions) {
    if (strcmp(regionName, reg->id) == 0) {
      return iregion;
    }
    ++iregion;
  }
  return -1;
}

bool LammpsInterface::region_bounds(const char * regionName,
  double & xmin, double & xmax,
  double & ymin, double & ymax,
  double & zmin, double & zmax,
  double & xscale, double & yscale, double & zscale) const
{
  int iRegion = region_id(regionName);
  if (iRegion < 0) throw ATC_Error("Unknown region " + to_string(regionName));
  xscale = region_xscale(iRegion);
  yscale = region_yscale(iRegion);
  zscale = region_zscale(iRegion);
  xmin = region_xlo(iRegion);
  xmax = region_xhi(iRegion);
  ymin = region_ylo(iRegion);
  ymax = region_yhi(iRegion);
  zmin = region_zlo(iRegion);
  zmax = region_zhi(iRegion);
  if (strcmp(region_style(iRegion),"block")==0) { return true; }
  else                                          { return false; }
}

void LammpsInterface::minimum_image(double & dx, double & dy, double & dz) const {
  lammps_->domain->minimum_image(dx,dy,dz);
}

void LammpsInterface::closest_image(const double * const xi, const double * const xj, double * const xjImage) const {
  lammps_->domain->closest_image(xi,xj,xjImage);
}


// -----------------------------------------------------------------
//  update interface methods
// -----------------------------------------------------------------
LammpsInterface::UnitsType LammpsInterface::units_style() const
{
    if      (strcmp(lammps_->update->unit_style,"lj") == 0)    return LJ;
    else if (strcmp(lammps_->update->unit_style,"real") == 0)  return REAL;
    else if (strcmp(lammps_->update->unit_style,"metal") == 0) return METAL;
    else return UNKNOWN;
}

double LammpsInterface::convert_units(double value, UnitsType in, UnitsType out, int massExp, int lenExp, int timeExp, int engExp) const
{
    double ps2fs = 1.e3;
    double eV2kcal = 23.069;
    if      (in==REAL) {
      if (out==METAL) {
        return value*pow(ps2fs,-timeExp)*pow(eV2kcal,-engExp);
      }
      else if (out==ATC) {
        if      (units_style()==REAL) {
          return value;
        }
        else if (units_style()==METAL) {
          return convert_units(value, METAL, out, massExp, lenExp, timeExp)*1.0;
        }
      }
      else throw ATC_Error("can't convert");
    }
    else if (in==METAL) {
      if (out==REAL) {
        return value*pow(ps2fs,timeExp)*pow(eV2kcal,engExp);
      }
      else if (out==ATC) {
        if      (units_style()==REAL) {
          return convert_units(value, REAL, out, massExp, lenExp, timeExp)*1.0;
        }
        else if (units_style()==METAL) {
          return value;
        }
      }
      else throw ATC_Error("can't convert");
    }
    else throw ATC_Error("can't convert");
    return value;
}


// -----------------------------------------------------------------
//  lattice interface methods
// -----------------------------------------------------------------

double LammpsInterface::xlattice() const { return lammps_->domain->lattice->xlattice; }

double LammpsInterface::ylattice() const { return lammps_->domain->lattice->ylattice; }

double LammpsInterface::zlattice() const { return lammps_->domain->lattice->zlattice; }

LammpsInterface::LatticeType LammpsInterface::lattice_style() const
{
  if (lammps_->domain->lattice)
    return (LammpsInterface::LatticeType)lammps_->domain->lattice->style;
  else
    throw ATC_Error("Lattice has not been defined");
}

//* returns the number of basis vectors
int LammpsInterface::n_basis() const
{
  return lammps_->domain->lattice->nbasis;
}

//* returns the basis vectors, transformed to the box coords
void LammpsInterface::basis_vectors(double **basis) const
{
  LAMMPS_NS::Lattice *lattice = lammps_->domain->lattice;
  int i,j;
  double origin[3] = {0.0, 0.0, 0.0};
  lattice->lattice2box(origin[0], origin[1], origin[2]);
  for (i=0; i<n_basis(); i++)
  {
    memcpy(basis[i],lattice->basis[i],3*sizeof(double));
    lattice->lattice2box(basis[i][0], basis[i][1], basis[i][2]);
    for (j=0; j<3; j++)  basis[i][j] -= origin[j];
  }
}

//* gets the (max) lattice constant
double LammpsInterface::max_lattice_constant() const
{
  double a1[3], a2[3], a3[3];
  unit_cell(a1,a2,a3);
  double a = norm(a1);
  a = max(a,norm(a2));
  a = max(a,norm(a3));
  return a;
}

//* computes a cutoff distance halfway between 1st and 2nd nearest neighbors
double LammpsInterface::near_neighbor_cutoff() const
{
  double cutoff;
  double alat = LammpsInterface::max_lattice_constant();
  LatticeType type = lattice_style();
  if (type == LammpsInterface::SC) {
    cutoff = 0.5*(1.0+sqrt(2.0))*alat;
  } else if (type == LammpsInterface::BCC) {
    cutoff = 0.5*(0.5*sqrt(3.0)+1.0)*alat;
  } else if (type == LammpsInterface::FCC) {
    cutoff = 0.5*(1.0/sqrt(2.0)+1.0)*alat;
  } else if (type == LammpsInterface::HCP) {
    cutoff = 0.5*(1.0/sqrt(2.0)+1.0)*alat;
  } else if (type == LammpsInterface::DIAMOND) {
    cutoff = 0.5*(0.25*sqrt(3.0)+1.0/sqrt(2.0))*alat;
  } else if (type == LammpsInterface::SQ) {
    cutoff = 0.5*(1.0+sqrt(2.0))*alat;
  } else if (type == LammpsInterface::SQ2) {
    cutoff = 0.5*(1.0/sqrt(2.0)+1.0)*alat;
  } else if (type == LammpsInterface::HEX) {
    cutoff = 0.5*(1.0/sqrt(3.0)+1.0)*alat;
  } else {
    throw ATC_Error("Unknown lattice type");
  }
  return cutoff;
}

//* gets the unit cell vectors
void LammpsInterface::unit_cell(double *a1, double *a2, double *a3) const
{
  int i, j;
  double *a[3] = {a1,a2,a3};
  double origin[3] = {0.0,0.0,0.0};
  LAMMPS_NS::Lattice *lattice = lammps_->domain->lattice;
  // transform origin
  lattice->lattice2box(origin[0], origin[1], origin[2]);

  // copy reference lattice vectors
  memcpy(a[0], lattice->a1, 3*sizeof(double));
  memcpy(a[1], lattice->a2, 3*sizeof(double));
  memcpy(a[2], lattice->a3, 3*sizeof(double));

  for (i=0; i<3; i++)
  {
    lattice->lattice2box(a[i][0], a[i][1], a[i][2]);
    for (j=0; j<3; j++)  a[i][j] -= origin[j];
  }
}

//* gets number of atoms in a unit cell
int LammpsInterface::num_atoms_per_cell() const
{
  int naCell = 0;
  LatticeType type = lattice_style();
  if      (type == LammpsInterface::SC)  naCell = 1;
  else if (type == LammpsInterface::BCC) naCell = 2;
  else if (type == LammpsInterface::FCC) naCell = 4;
  else if (type == LammpsInterface::DIAMOND) naCell = 8;
  else if (comm_rank()==0) {
    //{throw ATC_Error("lattice style not currently supported by ATC");}
    print_msg_once("WARNING:  Cannot get number of atoms per cell from lattice");
    naCell = 1;
  }
  return naCell;
}

//* gets tributary volume for an atom
double LammpsInterface::volume_per_atom() const
{
  double naCell = num_atoms_per_cell();
  double volPerAtom =
    xlattice() * ylattice() * zlattice() / naCell;
  return volPerAtom;
}

//* gets lattice basis
void LammpsInterface::lattice(MATRIX &N, MATRIX &B) const
{
  int nbasis = n_basis();
  double **basis = new double*[nbasis];
  N.reset(3,3);
  B.reset(3,nbasis);
  for (int i=0; i<nbasis; i++)  basis[i] = column(B,i).ptr();
  basis_vectors(basis);
  unit_cell(column(N,0).ptr(),
                                column(N,1).ptr(),
                                column(N,2).ptr());
  delete [] basis;
}

// -----------------------------------------------------------------
//  force interface methods
// -----------------------------------------------------------------

double LammpsInterface::boltz() const{ return lammps_->force->boltz;      }

double LammpsInterface::mvv2e() const{ return lammps_->force->mvv2e;      }

double LammpsInterface::ftm2v()const { return lammps_->force->ftm2v;      }

double LammpsInterface::nktv2p()const{ return lammps_->force->nktv2p;     }

double LammpsInterface::qqr2e() const{ return lammps_->force->qqr2e;      }

double LammpsInterface::qe2f()  const{ return lammps_->force->qe2f;       }
double LammpsInterface::dielectric()const{return lammps_->force->dielectric; }

double LammpsInterface::qqrd2e()const{ return lammps_->force->qqrd2e;     }

double LammpsInterface::qv2e()  const{ return qe2f()*ftm2v();             }

double LammpsInterface::pair_force(int i, int j, double rsq,
  double & fmag_over_rmag) const
{
  int itype = (lammps_->atom->type)[i];
  int jtype = (lammps_->atom->type)[j];
  // return value is the energy
  if (rsq < (lammps_->force->pair->cutsq)[itype][jtype]) {
    return lammps_->force->pair->single(i,j,itype,jtype,
    rsq,1.0,1.0,fmag_over_rmag);
  }
  return 0.0;
}
double LammpsInterface::pair_force(int n, double rsq,
  double & fmag_over_rmag) const
{
  int i = bond_list_i(n);
  int j = bond_list_j(n);
  int type = bond_list_type(n);
  // return value is the energy
  return lammps_->force->bond->single(type,rsq,i,j,fmag_over_rmag);
}
double LammpsInterface::pair_force(
                                   map< std::pair< int,int >,int >::const_iterator itr, double rsq,
  double & fmag_over_rmag, int nbonds) const
{
  int n = itr->second;
  if (n < nbonds) {
    return pair_force(n, rsq,fmag_over_rmag);
  }
  else {
    std::pair <int,int>  ij = itr->first;
    int i = ij.first;
    int j = ij.second;
    return pair_force(i,j, rsq,fmag_over_rmag);
  }
}
double LammpsInterface::pair_force(
  std::pair< std::pair< int,int >,int > apair, double rsq,
  double & fmag_over_rmag, int nbonds) const
{
  int n = apair.second;
  if (n < nbonds) {
    return pair_force(n, rsq,fmag_over_rmag);
  }
  else {
    std::pair <int,int>  ij = apair.first;
    int i = ij.first;
    int j = ij.second;
    return pair_force(i,j, rsq,fmag_over_rmag);
  }
}
double LammpsInterface::bond_stiffness(int i, int j, double rsq0) const
{
  const double perturbation = 1.e-8;
  double rsq1 = sqrt(rsq0)+perturbation;
  rsq1 *= rsq1;
  double f0,f1;
  pair_force(i,j,rsq0,f0);
  pair_force(i,j,rsq1,f1);
  double k =  (f1-f0)/perturbation;
  return k;
}

double LammpsInterface::pair_cutoff() const
{
  return lammps_->force->pair->cutforce;
}

void LammpsInterface::pair_reinit() const
{
  lammps_->force->pair->reinit();
}

int LammpsInterface::single_enable() const
{
  return lammps_->force->pair->single_enable; // all bonds have a single
}

//* insertion/deletion functions : see FixGCMC
//* delete atom
int LammpsInterface::delete_atom(int id) const
{
  LAMMPS_NS::Atom * atom  = lammps_->atom;
  atom->avec->copy(atom->nlocal-1,id,1);
  atom->nlocal--;
  return atom->nlocal;
}

//* insert atom
int LammpsInterface::insert_atom(int atype, int amask,
  double *ax, double *av, double aq) const
{
  LAMMPS_NS::Atom * atom  = lammps_->atom;
  atom->avec->create_atom(atype,ax);
  int m = atom->nlocal - 1;
  atom->mask[m] = amask;
  atom->v[m][0] = av[0];
  atom->v[m][1] = av[1];
  atom->v[m][2] = av[2];
  if (aq != 0) atom->q[m] = aq;

  int nfix = lammps_->modify->nfix;
  LAMMPS_NS::Fix **fix = lammps_->modify->fix;
  for (int j = 0; j < nfix; j++) {
    if (fix[j]->create_attribute) fix[j]->set_arrays(m);
  }
  return m;
}

int LammpsInterface::reset_ghosts(int deln) const
{
  LAMMPS_NS::Atom * atom  = lammps_->atom;
//ATC::LammpsInterface::instance()->print_msg("reset_ghosts del n "+to_string(deln));
  if (atom->tag_enable) {
    atom->natoms += deln;
//ATC::LammpsInterface::instance()->print_msg("reset_ghosts natoms "+to_string(atom->natoms));
    if (deln > 0) { atom->tag_extend(); }
    if (atom->map_style) atom->map_init();
  }
  atom->nghost = 0;
  lammps_->comm->borders();
//ATC::LammpsInterface::instance()->print_msg("reset_ghosts nghosts "+to_string(atom->nghost));
  return atom->nghost;
}

//* energy for interactions within the shortrange cutoff
double LammpsInterface::shortrange_energy(double *coord,
                                          int itype, int id, double /* max */) const
{
  LAMMPS_NS::Atom * atom  = lammps_->atom;
  double **x = atom->x;
  int *type = atom->type;
  int nall = atom->nlocal+ atom->nghost;
  LAMMPS_NS::Pair *pair = lammps_->force->pair;
  double **cutsq =  lammps_->force->pair->cutsq;

  double fpair = 0.0; // an output of single
  double factor_coul = 1.0;
  double factor_lj = 1.0;

  double total_energy = 0.0;
  for (int j = 0; j < nall; j++) {
    if (id == j) continue;
    // factor_lj = special_lj[sbmask(j)];
    // factor_coul = special_coul[sbmask(j)];
    //j &= NEIGHMASK;
    double delx = coord[0] - x[j][0];
    double dely = coord[1] - x[j][1];
    double delz = coord[2] - x[j][2];
    double rsq = delx*delx + dely*dely + delz*delz;
    int jtype = type[j];
    double cut2 = cutsq[itype][jtype];
    if (rsq < cut2) {
      total_energy += pair->single(id,j,itype,jtype,rsq,factor_coul,factor_lj,fpair); //id is for charge lookup
    }
  }
  return total_energy;
}
double LammpsInterface::shortrange_energy(int id, double max) const
{
  double *x = (lammps_->atom->x)[id];
  int type  = (lammps_->atom->type)[id];
  return shortrange_energy(x,type,id,max);
}

POTENTIAL LammpsInterface::potential() const
{
  // find pair style - FRAGILE
  const int nStyles = 4;
  string pairStyles[nStyles] = {"lj/cut",
                          "lj/cut/coul/long",
                          "lj/cut/coul/cut",
                          "lj/charmm/coul/long"};
  LAMMPS_NS::Pair *pair = nullptr;
  for (int i = 0; i < nStyles; i++){
    pair = lammps_->force->pair_match(pairStyles[i].c_str(),1);
    if (pair != nullptr) break;
  }
  return pair;
}

int LammpsInterface::type_to_groupbit(int itype) const
{
  LAMMPS_NS::Atom * atom  = lammps_->atom;
  int groupbit = -1;
  int *type = atom->type;
  int *mask = atom->mask;
  for (int i = 0; i < nlocal(); i++) {
    if (type[i] == itype) {
      groupbit = mask[i];
      break;
    }
  }
  return int_allmax(groupbit);
}

bool LammpsInterface::epsilons(int itype, POTENTIAL pair, double * epsilon0) const
{
  // grab energy parameters
  char * pair_parameter = new char[8];
  strcpy(pair_parameter,"epsilon");
  int dim = 2; // a return value for extract
  double ** epsilons = (double**) ( pair->extract(pair_parameter,dim) );
  delete [] pair_parameter;
  if (epsilons == nullptr) return false;
  //if (epsilons == nullptr) error->all(FLERR,"Fix concentration adapted pair style parameter not supported");
  int i1,i2;
  for (int i=1; i < ntypes()+1; i++) {
    if (i < itype) { i1 = i; i2 = itype; }
    else           { i2 = i; i1 = itype; }
    epsilon0[i-1] = epsilons[i1][i2];
  }
  return true;
}

bool LammpsInterface::set_epsilons(int itype, POTENTIAL pair, double * epsilon) const
{
  // grab energy parameters
  char * pair_parameter = new char[8];
  strcpy(pair_parameter,"epsilon");
  int dim = 2; // a return value for extract
  double ** epsilons = (double**) ( pair->extract(pair_parameter,dim) );
  delete [] pair_parameter;
  if (epsilons == nullptr) return false;
  //if (epsilons == nullptr) error->all(FLERR,"Fix concentration adapted pair style parameter not supported");
  // scale interactions
  int i1,i2;
  for (int i = 1; i < ntypes()+1; i++) {
    if (i < itype) { i1 = i; i2 = itype; }
    else           { i2 = i; i1 = itype; }
    epsilons[i1][i2] = epsilon[i-1];
  }

  return true;
}

int LammpsInterface::set_charge(int itype, double charge) const
{
  int nlocal = lammps_->atom->nlocal;
  int *type = lammps_->atom->type;
  double *q = lammps_->atom->q;
  int count = 0;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] == itype) {
      q[i] = charge;
      count++;
    }
  }
  return count;
}

int LammpsInterface::change_type(int itype, int jtype) const
{
  int nlocal = lammps_->atom->nlocal;
  int *type = lammps_->atom->type;
  int count = 0;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] == itype) {
      type[i] = jtype;
      count++;
    }
  }
  return count;
}

int LammpsInterface::count_type(int itype) const
{
  int nlocal = lammps_->atom->nlocal;
  int *type = lammps_->atom->type;
  int count = 0;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] == itype) { count++; }
  }
  return int_allsum(count);
}

// random number generators
RNG_POINTER LammpsInterface::random_number_generator() const {
  RNG_POINTER p = new LAMMPS_NS::RanPark(lammps_,seed_);
  return p;
}
double LammpsInterface::random_uniform(RNG_POINTER p) const {
  return p->uniform();
}
double LammpsInterface::random_normal (RNG_POINTER p) const {
  return p->gaussian();
}
int LammpsInterface::random_state (RNG_POINTER p) const {
  return p->state();
}
void LammpsInterface::set_random_state (RNG_POINTER p, int seed) const {
  return p->reset(seed);
}
void LammpsInterface::advance_random_generator (RNG_POINTER p, int n) const {
  advance_random_uniform(p,n);
}
void LammpsInterface::advance_random_uniform (RNG_POINTER p, int n) const {
  for (int i = 0; i < n; i++)  p->uniform();
}
void LammpsInterface::advance_random_normal (RNG_POINTER p, int n) const {
  for (int i = 0; i < n; i++)  p->gaussian();
}

//* Boltzmann's constant in M,L,T,t units
double LammpsInterface::kBoltzmann() const  {
  return (lammps_->force->boltz)/(lammps_->force->mvv2e);
}

//* Planck's constant
double LammpsInterface::hbar() const  {
  const int UNITS_STYLE = (int) units_style();
  double hbar = 1.0; // LJ: Dimensionless
  if      (UNITS_STYLE == 2) hbar = 15.1685792814; // Real: KCal/mol-fs
  else if (UNITS_STYLE == 3) hbar = 0.000658212202469; // Metal: eV-ps
  return hbar;
}

//* Dulong-Petit heat capacity
double LammpsInterface::heat_capacity() const  {
  double rhoCp = dimension()*kBoltzmann()/volume_per_atom();
  return rhoCp;
}

//* reference mass density for a *unit cell*
// all that is needed is a unit cell: volume, types, mass per type
double LammpsInterface::mass_density(int* numPerType) const
{
  const double *mass         = lammps_->atom->mass;
  if (!mass) throw ATC_Error("cannot compute a mass density: no mass");
  const int    ntypes        = lammps_->atom->ntypes;
  const int    *mass_setflag = lammps_->atom->mass_setflag;
  const int    *type         = lammps_->atom->type;
  double naCell = num_atoms_per_cell();
  double vol = volume_per_atom();
  if (numPerType) {
    double m = 0.0;
    double n = 0;
    for (int i = 0; i < ntypes; i++) {
       m += numPerType[i]*mass[i+1];
       n += numPerType[i];
    }
    if (n>naCell) throw ATC_Error("cannot compute a mass density: too many basis atoms");
    return m/n/vol;
  }
  // since basis->type map not stored only monatomic lattices are automatic
  // if not given a basis try to guess
  else {
    if (ntypes == 1) {
      if ((this->natoms()==0) && mass_setflag[1]) {
        return mass[1]/vol;
      }
      else {
        if (type)                 return mass[type[0]]/vol;
        else if (mass_setflag[1]) return mass[1]/vol;
      }
    }
  }
  throw ATC_Error("cannot compute a mass density");
  return 0.0;
}

//* permittivity of free space
double LammpsInterface::epsilon0() const
{
  return qe2f()/(4.*PI*qqr2e());
}

//* Coulomb's constant
double LammpsInterface::coulomb_constant() const
{
  return qqr2e()/qe2f();
}

//* special coulombic interactions
double * LammpsInterface::special_coul() const
{
  return lammps_->force->special_coul;
}

//* flag for newton
int LammpsInterface::newton_pair() const
{
  return lammps_->force->newton_pair;
}

// -----------------------------------------------------------------
//  group interface methods
// -----------------------------------------------------------------

int LammpsInterface::ngroup() const { return lammps_->group->ngroup; }

int LammpsInterface::group_bit(string name) const
{
  return group_bit(group_index(name));
}

int LammpsInterface::group_bit(int iGroup) const
{
  int mybit = 0;
  mybit |= lammps_->group->bitmask[iGroup];
  if (mybit < 0 || mybit > MAX_GROUP_BIT) {
    string msg("LammpsInterface::group_bit() lammps group bit "+to_string(mybit)+" is out of range 0:"+to_string(MAX_GROUP_BIT));
    throw ATC_Error(msg);
  }

  return mybit;
}

int LammpsInterface::group_index(string name) const
{
  int igroup = lammps_->group->find(name.c_str());
  if (igroup == -1) {
    string msg("LammpsInterface::group_index() lammps group "+name+" does not exist");
    throw ATC_Error(msg);
  }
  return igroup;
}

int LammpsInterface::group_inverse_mask(int iGroup) const
{
  return lammps_->group->inversemask[iGroup];
}

char * LammpsInterface::group_name(int iGroup) const
{
  return lammps_->group->names[iGroup];
}

void LammpsInterface::group_bounds(int iGroup, double * b) const
{
  lammps_->group->bounds(iGroup, b);
}

// -----------------------------------------------------------------
//  memory interface methods
// -----------------------------------------------------------------


double * LammpsInterface::create_1d_double_array(int length, const char *name) const {
  double * myArray;
  return lammps_->memory->create(myArray, length, name);
}

double *LammpsInterface::grow_1d_double_array(double *array,
                                              int length,
                                              const char *name) const
{
  return lammps_->memory->grow(array, length, name);
}

void LammpsInterface::destroy_1d_double_array(double * d) const {
  lammps_->memory->destroy(d);
}

double ** LammpsInterface::create_2d_double_array(int n1, int n2, const char *name) const {
  double ** myArray;
  return lammps_->memory->create(myArray, n1, n2, name);
}

void LammpsInterface::destroy_2d_double_array(double **d) const {
  lammps_->memory->destroy(d);
}

double **LammpsInterface::grow_2d_double_array(double **array,
                                               int n1,
                                               int n2,
                                               const char *name) const
{
  return lammps_->memory->grow(array, n1, n2, name);
}


int * LammpsInterface::create_1d_int_array(int length, const char *name) const {
  int * myArray;
  return lammps_->memory->create(myArray, length, name);
}

int *LammpsInterface::grow_1d_int_array(int *array,
                                        int length,
                                        const char *name) const
{
  return lammps_->memory->grow(array, length, name);
}

void LammpsInterface::destroy_1d_int_array(int * d) const {
  lammps_->memory->destroy(d);
}

int ** LammpsInterface::create_2d_int_array(int n1, int n2, const char *name) const {
  int ** myArray;
  return lammps_->memory->create(myArray, n1, n2, name);
}

void LammpsInterface::destroy_2d_int_array(int **i) const {
  lammps_->memory->destroy(i);
}

int ** LammpsInterface::grow_2d_int_array(int **array, int n1, int n2, const char *name) const {
  return lammps_->memory->grow(array, n1, n2, name);
}

// -----------------------------------------------------------------
//  update interface methods
// -----------------------------------------------------------------

double LammpsInterface::dt() const        { return lammps_->update->dt; }

bigint LammpsInterface::ntimestep() const { return lammps_->update->ntimestep; }

int LammpsInterface::nsteps() const    { return lammps_->update->nsteps; }


// -----------------------------------------------------------------
//  neighbor list interface methods
// -----------------------------------------------------------------

int LammpsInterface::sbmask(int j) const {return j >> SBBITS & 3; }

void LammpsInterface::set_list(int /* id */, LAMMPS_NS::NeighList *ptr) { list_ = ptr; }

int   LammpsInterface::neighbor_list_inum() const { return list_->inum; }

int * LammpsInterface::neighbor_list_numneigh() const { return list_->numneigh; }

int * LammpsInterface::neighbor_list_ilist() const { return list_->ilist; }

int ** LammpsInterface::neighbor_list_firstneigh() const { return list_->firstneigh; }

int   LammpsInterface::neighbor_ago() const { return lammps_->neighbor->ago; }

int   LammpsInterface::reneighbor_frequency() const {return lammps_->neighbor->every; }

// -----------------------------------------------------------------
//  bond list interface methods
// -----------------------------------------------------------------

int   LammpsInterface::bond_list_length() const { return lammps_->neighbor->nbondlist; }
int**   LammpsInterface::bond_list() const { return lammps_->neighbor->bondlist; }

// -----------------------------------------------------------------
//  region interface methods
// -----------------------------------------------------------------

char * LammpsInterface::region_name(int iRegion)  const
{
  auto regions = lammps_->domain->get_region_list();
  return regions[iRegion]->id;
}

char * LammpsInterface::region_style(int iRegion) const
{
  auto regions = lammps_->domain->get_region_list();
  return regions[iRegion]->style;
}

double LammpsInterface::region_xlo(int iRegion) const
{
  auto regions = lammps_->domain->get_region_list();
  return regions[iRegion]->extent_xlo;
}

double LammpsInterface::region_xhi(int iRegion) const
{
  auto regions = lammps_->domain->get_region_list();
  return regions[iRegion]->extent_xhi;
}

double LammpsInterface::region_ylo(int iRegion) const
{
  auto regions = lammps_->domain->get_region_list();
  return regions[iRegion]->extent_ylo;
}

double LammpsInterface::region_yhi(int iRegion) const
{
  auto regions = lammps_->domain->get_region_list();
  return regions[iRegion]->extent_yhi;
}

double LammpsInterface::region_zlo(int iRegion) const
{
  auto regions = lammps_->domain->get_region_list();
  return regions[iRegion]->extent_zlo;
}

double LammpsInterface::region_zhi(int iRegion) const
{
  auto regions = lammps_->domain->get_region_list();
  return regions[iRegion]->extent_zhi;
}

double LammpsInterface::region_xscale(int iRegion) const
{
  auto regions = lammps_->domain->get_region_list();
  return regions[iRegion]->xscale;
}

double LammpsInterface::region_yscale(int iRegion)  const
{
  auto regions = lammps_->domain->get_region_list();
  return regions[iRegion]->yscale;
}

double LammpsInterface::region_zscale(int iRegion)  const
{
  auto regions = lammps_->domain->get_region_list();
  return regions[iRegion]->zscale;
}

int LammpsInterface::region_match(int iRegion, double x, double y, double z) const {
  auto regions = lammps_->domain->get_region_list();
  return regions[iRegion]->match(x,y,z);
}

// -----------------------------------------------------------------
//  compute methods
// -----------------------------------------------------------------
COMPUTE_POINTER LammpsInterface::compute_pointer(string tag) const
{
  // get the compute id
  char * name = const_cast <char*> (tag.c_str());
  int id = lammps_->modify->find_compute(name);
  if (id < 0) {
    string msg("Could not find compute "+tag);
    msg += tag;
    throw ATC_Error(msg);
  }
  // get the compute
  LAMMPS_NS::Compute* cmpt = lammps_->modify->compute[id];
  // insert it into our set, recall it won't be added if it already exists
  computePointers_.insert(cmpt);
  return cmpt;
}

void LammpsInterface::computes_addstep(int step) const
{
  set<LAMMPS_NS::Compute * >::iterator iter;
  for (iter = computePointers_.begin(); iter != computePointers_.end(); iter++) {
    (*iter)->addstep(step);
  }
}

void LammpsInterface::compute_addstep(COMPUTE_POINTER computePointer, int step) const
{
  LAMMPS_NS::Compute* cmpt = const_to_active(computePointer);
  cmpt->addstep(step);
}

int LammpsInterface::compute_matchstep(COMPUTE_POINTER computePointer, int step) const
{
  LAMMPS_NS::Compute* cmpt = const_to_active(computePointer);
  return cmpt->matchstep(step);
}

void LammpsInterface::reset_invoked_flag(COMPUTE_POINTER computePointer) const
{
  LAMMPS_NS::Compute* cmpt = const_to_active(computePointer);
  cmpt->invoked_flag = 0;
}

int LammpsInterface::compute_ncols_peratom(COMPUTE_POINTER computePointer) const
{
  LAMMPS_NS::Compute* cmpt = const_to_active(computePointer);
  int ndof = cmpt->size_peratom_cols;
  if (ndof == 0 ) ndof = 1;
  return ndof;
}

double*  LammpsInterface::compute_vector_peratom(COMPUTE_POINTER computePointer) const
{
  LAMMPS_NS::Compute* cmpt = const_to_active(computePointer);
  if (!(cmpt->invoked_flag & INVOKED_PERATOM)) {
    cmpt->compute_peratom();
    cmpt->invoked_flag |= INVOKED_PERATOM;
  }
  return cmpt->vector_atom;
}

double**  LammpsInterface::compute_array_peratom(COMPUTE_POINTER computePointer) const
{
  LAMMPS_NS::Compute* cmpt = const_to_active(computePointer);
  if (!(cmpt->invoked_flag & INVOKED_PERATOM)) {
    cmpt->compute_peratom();
    cmpt->invoked_flag |= INVOKED_PERATOM;
  }
  return cmpt->array_atom;
}

LAMMPS_NS::Compute * LammpsInterface::const_to_active(COMPUTE_POINTER computePointer) const
{
  LAMMPS_NS::Compute* cmpt = const_cast<LAMMPS_NS::Compute* >(computePointer);
  set<LAMMPS_NS::Compute * >::iterator cmptPtr;
  cmptPtr = computePointers_.find(cmpt);
  if (cmptPtr != computePointers_.end())
    return *cmptPtr;
  else
    throw ATC_Error("Requested bad computer pointer");
}

// -----------------------------------------------------------------
//  compute pe/atom interface methods
//  - the only compute "owned" by ATC
// -----------------------------------------------------------------
int  LammpsInterface::create_compute_pe_peratom()  const
{
  char **list = new char*[4];
  string atomPeName = compute_pe_name();
  list[0] = (char *) atomPeName.c_str();
  list[1] = (char *) "all";
  list[2] = (char *) "pe/atom";
  list[3] = (char *) "pair";

  int icompute = lammps_->modify->find_compute(list[0]);
  if (icompute < 0) {
    lammps_->modify->add_compute(4,list);
    icompute = lammps_->modify->find_compute(list[0]);
  }
  delete [] list;
  if (! atomPE_ ) {
    atomPE_ = lammps_->modify->compute[icompute];
  }
  computePointers_.insert(atomPE_);
  stringstream ss;
  ss << "peratom PE compute created with ID: " << icompute;
  ATC::LammpsInterface::instance()->print_msg_once(ss.str());
  return icompute;
}

double * LammpsInterface::compute_pe_peratom() const
{
  if (atomPE_) {
    atomPE_->compute_peratom();
    return atomPE_->vector_atom;
  }
  else {
    return nullptr;
  }
}

/* ---------------------------------------------------------------------- */

void LammpsInterface::unwrap_coordinates(int iatom, double* xatom) const
{
  double **x = lammps_->atom->x;
  int *image = lammps_->atom->image;

  double *h   = lammps_->domain->h;
  double xprd = lammps_->domain->xprd;
  double yprd = lammps_->domain->yprd;
  double zprd = lammps_->domain->zprd;
  int xbox,ybox,zbox;

  // for triclinic, need to unwrap current atom coord via h matrix

  if (lammps_->domain->triclinic == 0) {
    xbox = (image[iatom] & 1023) - 512;
    ybox = (image[iatom] >> 10 & 1023) - 512;
    zbox = (image[iatom] >> 20) - 512;
    xatom[0] = x[iatom][0] + xbox*xprd;
    xatom[1] = x[iatom][1] + ybox*yprd;
    xatom[2] = x[iatom][2] + zbox*zprd;
  } else {
    xbox = (image[iatom] & 1023) - 512;
    ybox = (image[iatom] >> 10 & 1023) - 512;
    zbox = (image[iatom] >> 20) - 512;
    xatom[0] = x[iatom][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
    xatom[1] = x[iatom][1] + h[1]*ybox + h[3]*zbox;
    xatom[2] = x[iatom][2] + h[2]*zbox;
  }
}

/* ---------------------------------------------------------------------- */

LAMMPS_NS::PairEAM* LammpsInterface::pair_eam() const
{
  //if (typeid(lammps_->force->pair) == typeid(LAMMPS_NS::PairEAM)) {
  //  return lammps_->force->pair;
  //}
  LAMMPS_NS::PairEAM* pair_eam = dynamic_cast<LAMMPS_NS::PairEAM*> (lammps_->force->pair);
  if (pair_eam != nullptr) {
    return pair_eam;
  }
  else {
    throw ATC_Error("LAMMPS Pair object is not of the derived class type PairEAM");
  }
}

// -----------------------------------------------------------------
//  other methods
// -----------------------------------------------------------------

/** Return lammps pointer -- only as a last resort! */
LAMMPS_NS::LAMMPS * LammpsInterface::lammps_pointer() const { return lammps_; }

}
