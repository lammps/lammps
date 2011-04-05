// Header file for this class
#include "LammpsInterface.h"

// LAMMPS includes
#include "lammps.h"
#include "atom.h" // x, v, f
#include "domain.h" // for basing locations on regions
#include "region.h" // region bounding box and style
#include "force.h" // boltzman constant
#include "group.h" // atom masks
#include "memory.h" // grow atom information
#include "compute.h" // computes
#include "modify.h" // 
#include "neighbor.h" // neighbors
#include "neigh_list.h" // neighbor list
#include "update.h" // timestepping information
#include "pair.h" // pair potentials
#include "lattice.h" // lattice parameters
#include "comm.h" //

// ATC includes
#include "ATC_Error.h"
#include "MatrixLibrary.h"

// Other include files
#include "mpi.h"
#include <cstring>

namespace ATC {


LammpsInterface * LammpsInterface::myInstance_ = NULL;

// -----------------------------------------------------------------
//  instance()
// -----------------------------------------------------------------
LammpsInterface * LammpsInterface::instance()
{
  if (myInstance_ == NULL) {
    myInstance_ = new LammpsInterface();
  }
  return myInstance_;
}

// -----------------------------------------------------------------
//  constructor
// -----------------------------------------------------------------
LammpsInterface::LammpsInterface()
  : lammps_(NULL),
    atomPE_(NULL),
    commRank_(0)
{
}

// -----------------------------------------------------------------
//  general interface methods
// -----------------------------------------------------------------

MPI_Comm LammpsInterface::world() { return lammps_->world; }

// -----------------------------------------------------------------
//  atom interface methods
// -----------------------------------------------------------------

int LammpsInterface::nlocal() { return lammps_->atom->nlocal; }

int LammpsInterface::nghost() { return lammps_->atom->nghost; }

double LammpsInterface::natoms() { return lammps_->atom->natoms; }

int LammpsInterface::nmax() { return lammps_->atom->nmax; }

int LammpsInterface::ntypes() { return lammps_->atom->ntypes; }

double ** LammpsInterface::xatom() { return lammps_->atom->x; }

const double ** LammpsInterface::xatom() const { return (const double**)(lammps_->atom->x); }

double ** LammpsInterface::vatom() { return lammps_->atom->v; }

double ** LammpsInterface::fatom() { return lammps_->atom->f; }

int * LammpsInterface::atom_mask() { return lammps_->atom->mask; }

int * LammpsInterface::atom_type() { return lammps_->atom->type; }

int * LammpsInterface::atom_tag() { return lammps_->atom->tag; }

double * LammpsInterface::atom_mass() { return lammps_->atom->mass; }

double   LammpsInterface::atom_mass(int iType) { return lammps_->atom->mass[iType]; }

double * LammpsInterface::atom_rmass() { return lammps_->atom->rmass; }

double * LammpsInterface::atom_charge() { return lammps_->atom->q; }
  
// -----------------------------------------------------------------
//  domain interface methods
// -----------------------------------------------------------------

int LammpsInterface::dimension() { return lammps_->domain->dimension; }

int LammpsInterface::nregion() { return lammps_->domain->nregion; }

void LammpsInterface::get_box_bounds(double & boxxlo, double & boxxhi,
				     double & boxylo, double & boxyhi,
				     double & boxzlo, double &boxzhi) 
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

int LammpsInterface::xperiodic() { return lammps_->domain->xperiodic; }

int LammpsInterface::yperiodic() { return lammps_->domain->yperiodic; }

int LammpsInterface::zperiodic() { return lammps_->domain->zperiodic; }

int LammpsInterface::nperiodic() 
{ 
  int nprd = 0;
  if ( lammps_->domain->xperiodic > 0 ) { nprd++ ; }
  if ( lammps_->domain->yperiodic > 0 ) { nprd++ ; }
  if ( lammps_->domain->zperiodic > 0 ) { nprd++ ; }
  return nprd;
}

double LammpsInterface::domain_xprd() { return lammps_->domain->xprd; }

double LammpsInterface::domain_yprd() { return lammps_->domain->yprd; }

double LammpsInterface::domain_zprd() { return lammps_->domain->zprd; }

double LammpsInterface::domain_xy() { return lammps_->domain->xy; }

double LammpsInterface::domain_xz() { return lammps_->domain->xz; }

double LammpsInterface::domain_yz() { return lammps_->domain->yz; }

int LammpsInterface::domain_triclinic() { return lammps_->domain->triclinic; }

void LammpsInterface::get_box_periodicity(int & xperiodic,
					  int & yperiodic,
					  int & zperiodic) 
{
  xperiodic = lammps_->domain->xperiodic;
  yperiodic = lammps_->domain->yperiodic;
  zperiodic = lammps_->domain->zperiodic;
}

int LammpsInterface::get_region_id(const char * regionName) {
  int nregion = this->nregion();
  for (int iregion = 0; iregion < nregion; iregion++) {
    if (strcmp(regionName, region_name(iregion)) == 0) {
      return iregion;
    }
  }
  throw ATC_Error(0,"Region has not been defined");
  return -1;
}
// -----------------------------------------------------------------
//  update interface methods
// -----------------------------------------------------------------
LammpsInterface::UnitsType LammpsInterface::units_style(void)
{
    if      (strcmp(lammps_->update->unit_style,"lj") == 0)    return LJ;
    else if (strcmp(lammps_->update->unit_style,"real") == 0)  return REAL;
    else if (strcmp(lammps_->update->unit_style,"metal") == 0) return METAL;
    else return UNKNOWN;
}


// -----------------------------------------------------------------
//  lattice interface methods
// -----------------------------------------------------------------

double LammpsInterface::xlattice() { return lammps_->domain->lattice->xlattice; }

double LammpsInterface::ylattice() { return lammps_->domain->lattice->ylattice; }

double LammpsInterface::zlattice() { return lammps_->domain->lattice->zlattice; }

LammpsInterface::LatticeType LammpsInterface::lattice_style() 
{ 
  if (lammps_->domain->lattice) 
    return (LammpsInterface::LatticeType)lammps_->domain->lattice->style; 
  else
    throw ATC_Error(0,"Lattice has not been defined");
}

//* retuns the number of basis vectors 
int LammpsInterface::get_n_basis()
{
  return lammps_->domain->lattice->nbasis;
}

//* returns the basis vectors, transformed to the box coords
void LammpsInterface::get_basis(double **basis)
{
  LAMMPS_NS::Lattice *lattice = lammps_->domain->lattice;
  int i,j;
  double origin[3] = {0.0, 0.0, 0.0};
  lattice->lattice2box(origin[0], origin[1], origin[2]);
  for (i=0; i<get_n_basis(); i++)
  {
    memcpy(basis[i],lattice->basis[i],3*sizeof(double));
    lattice->lattice2box(basis[i][0], basis[i][1], basis[i][2]);
    for (j=0; j<3; j++)  basis[i][j] -= origin[j];
  }
}

//* gets the unit cell vectors
void LammpsInterface::get_unit_cell(double *a1, double *a2, double *a3)
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
int LammpsInterface::num_atoms_per_cell(void) 
{
  int naCell = 0;
  LatticeType type = lattice_style();
  if      (type == LammpsInterface::SC)  naCell = 1;
  else if (type == LammpsInterface::BCC) naCell = 2;
  else if (type == LammpsInterface::FCC) naCell = 4;
  else if (type == LammpsInterface::DIAMOND) naCell = 8;
  else if (comm_rank()==0) {
    //{throw ATC_Error(0,"lattice style not currently supported by ATC");}
    cout << "ATC WARNING:  Cannot get number of atoms per cell from lattice\n";
    naCell = 1; //HACK to enable us to keep going since this is only used to compute volume per atom
                //     ATC modes with a user specified atomic volume or using only volumetric quantities are fine
  }
  return naCell;
}

//* gets tributary volume for an atom
double LammpsInterface::volume_per_atom(void) 
{
  double naCell = num_atoms_per_cell();
  double volPerAtom =  
    xlattice() * ylattice() * zlattice() / naCell;
  return volPerAtom;
}

//* gets lattice basis
void LammpsInterface::get_lattice(MATRIX &N, MATRIX &B)
{
  int nbasis = get_n_basis();
  double **basis = new double*[nbasis];
  N.reset(3,3);
  B.reset(3,nbasis);
  for (int i=0; i<nbasis; i++)  basis[i] = column(B,i).get_ptr();
  get_basis(basis);
  get_unit_cell(column(N,0).get_ptr(),
                                column(N,1).get_ptr(),
                                column(N,2).get_ptr());
  delete [] basis;
}


// -----------------------------------------------------------------
//  force interface methods
// -----------------------------------------------------------------

double LammpsInterface::boltz()      { return lammps_->force->boltz;      }

double LammpsInterface::mvv2e()      { return lammps_->force->mvv2e;      }     

double LammpsInterface::ftm2v()      { return lammps_->force->ftm2v;      }     

double LammpsInterface::nktv2p()     { return lammps_->force->nktv2p;     }    

double LammpsInterface::qqr2e()      { return lammps_->force->qqr2e;      }     

double LammpsInterface::qe2f()       { return lammps_->force->qe2f;       }      

double LammpsInterface::dielectric() { return lammps_->force->dielectric; }

double LammpsInterface::qqrd2e()     { return lammps_->force->qqrd2e;     }

double LammpsInterface::pair_force(int i, int j, double rsq, 
  double & fmag_over_rmag)     
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

double LammpsInterface::pair_cutoff()
{
  return lammps_->force->pair->cutforce;
}

int LammpsInterface::single_enable()
{
  return lammps_->force->pair->single_enable;
}

//* Boltzmann's constant in M,L,T,t units
double LammpsInterface::kBoltzmann()  { 
  return (lammps_->force->boltz)/(lammps_->force->mvv2e);
}

//* Dulong-Petit heat capacity
double LammpsInterface::heat_capacity()  { 
  double rhoCp = dimension()*kBoltzmann()/volume_per_atom();
  return rhoCp;
}

//* reference mass density
double LammpsInterface::mass_density()  
{ 
  const int    ntypes        = lammps_->atom->ntypes;
  const int    *mass_setflag = lammps_->atom->mass_setflag;
  const int    *type         = lammps_->atom->type;
  const double *mass         = lammps_->atom->mass;
  const double *rmass        = lammps_->atom->rmass;
  // NOTE currently assumes all atoms have same mass and volume
  //      in the future, mass and volume will be different but density should be
  //      an atom indepedent quantity
  if (mass) {
    if (type) return mass[type[0]]/volume_per_atom();
    // no type array - just use first type that has a set mass
    for (int i=1; i<=ntypes; i++) {
      if (mass_setflag[i]) return mass[i]/volume_per_atom();
    }
    // NOTE: no masses specified in input file should we warn the user of this?
    return 0.0;
  }
  // NOTE is this valid - lammps likes to not use 0 index
  if (rmass) return rmass[0]/volume_per_atom();
  return 0.0;
}

// -----------------------------------------------------------------
//  group interface methods
// -----------------------------------------------------------------

int LammpsInterface::ngroup() { return lammps_->group->ngroup; }

int LammpsInterface::group_bit(int iGroup) { return lammps_->group->bitmask[iGroup]; }

int LammpsInterface::find_group(const char * c) { return lammps_->group->find(c); }

int LammpsInterface::group_inverse_mask(int iGroup) 
{
  return lammps_->group->inversemask[iGroup]; 
}

char * LammpsInterface::group_name(int iGroup) 
{ 
  return lammps_->group->names[iGroup]; 
}

void LammpsInterface::group_bounds(int iGroup, double * b) 
{
  lammps_->group->bounds(iGroup, b);
}


// -----------------------------------------------------------------
//  memory interface methods
// -----------------------------------------------------------------

double * LammpsInterface::create_1d_double_array(int nlo, int nhi, const char *name) {
  double *array;
  return lammps_->memory->create1d_offset(array, nlo, nhi, name);
}

void LammpsInterface::destroy_1d_double_array(double * d, int i) {
  lammps_->memory->destroy1d_offset(d, i);
}

double ** LammpsInterface::create_2d_double_array(int n1, int n2, const char *name) {
  double **array;
  return lammps_->memory->create(array, n1, n2, name);
}

void LammpsInterface::destroy_2d_double_array(double **d) {
  lammps_->memory->destroy(d);
}

double **LammpsInterface::grow_2d_double_array(double **array, 
					       int n1, 
					       int n2, 
					       const char *name) 
{
  return lammps_->memory->grow(array, n1, n2, name);
}

int ** LammpsInterface::create_2d_int_array(int n1, int n2, const char *name) {
  int **array;
  return lammps_->memory->create(array, n1, n2, name);
}

void LammpsInterface::destroy_2d_int_array(int **i) {
  lammps_->memory->destroy(i);
}

int ** LammpsInterface::grow_2d_int_array(int **array, int n1, int n2, const char *name) {
  return lammps_->memory->grow(array, n1, n2, name);
}


// -----------------------------------------------------------------
//  update interface methods
// -----------------------------------------------------------------

double LammpsInterface::dt()        { return lammps_->update->dt; }

int LammpsInterface::ntimestep() { return lammps_->update->ntimestep; }

int LammpsInterface::nsteps()    { return lammps_->update->nsteps; }


// -----------------------------------------------------------------
//  neighbor list interface methods
// -----------------------------------------------------------------

void LammpsInterface::init_list(int id, LAMMPS_NS::NeighList *ptr) { list_ = ptr; }

int   LammpsInterface::neighbor_list_inum() { return list_->inum; }

int * LammpsInterface::neighbor_list_numneigh() { return list_->numneigh; }

int * LammpsInterface::neighbor_list_ilist() { return list_->ilist; }

int ** LammpsInterface::neighbor_list_firstneigh() { return list_->firstneigh; }

int   LammpsInterface::neighbor_ago() { return lammps_->neighbor->ago; }

// -----------------------------------------------------------------
//  region interface methods
// -----------------------------------------------------------------

char * LammpsInterface::region_name(int iRegion) 
{
  return lammps_->domain->regions[iRegion]->id; 
}

char * LammpsInterface::region_style(int iRegion)
{ 
  return lammps_->domain->regions[iRegion]->style; 
}

double LammpsInterface::region_xlo(int iRegion) 
{
  return lammps_->domain->regions[iRegion]->extent_xlo; 
}

double LammpsInterface::region_xhi(int iRegion) 
{
  return lammps_->domain->regions[iRegion]->extent_xhi; 
}

double LammpsInterface::region_ylo(int iRegion) 
{
  return lammps_->domain->regions[iRegion]->extent_ylo;
}

double LammpsInterface::region_yhi(int iRegion) 
{
  return lammps_->domain->regions[iRegion]->extent_yhi;
}

double LammpsInterface::region_zlo(int iRegion) 
{
  return lammps_->domain->regions[iRegion]->extent_zlo;
}

double LammpsInterface::region_zhi(int iRegion)
{
  return lammps_->domain->regions[iRegion]->extent_zhi;
}

double LammpsInterface::region_xscale(int iRegion)
{
  return lammps_->domain->regions[iRegion]->xscale; 
}

double LammpsInterface::region_yscale(int iRegion) 
{
  return lammps_->domain->regions[iRegion]->yscale;
}

double LammpsInterface::region_zscale(int iRegion) 
{
  return lammps_->domain->regions[iRegion]->zscale; 
}

int LammpsInterface::region_match(int iRegion, double x, double y, double z) { 
  return lammps_->domain->regions[iRegion]->match(x,y,z);
}

// -----------------------------------------------------------------
//  compute methods
// -----------------------------------------------------------------
int  LammpsInterface::find_compute(const char* tag) 
{
  // a clunky way to safely get rid of the const
  int n = strlen(tag) + 1;
  char* tag_copy = new char[n];
  strcpy(tag_copy,tag);

  int icompute = lammps_->modify->find_compute(tag_copy);
  if (icompute < 0) {
    string msg("Could not find compute ");
    msg += tag;
    throw ATC_Error(0,msg);
  }
  return icompute;
}

LAMMPS_NS::Compute*  LammpsInterface::get_compute(const char* tag) 
{
  int icompute = find_compute(tag);
  LAMMPS_NS::Compute* cmpt = lammps_->modify->compute[icompute];
  return cmpt;
}

double**  LammpsInterface::compute_vector_data(const char* tag) 
{
  LAMMPS_NS::Compute* cmpt = get_compute(tag);
  if (!(cmpt->invoked_flag & INVOKED_PERATOM)) {
    cmpt->compute_peratom();
    cmpt->invoked_flag |= INVOKED_PERATOM;
  }
  return cmpt->array_atom;
}

double*  LammpsInterface::compute_scalar_data(const char* tag) 
{
  LAMMPS_NS::Compute* cmpt = get_compute(tag);
  if (!(cmpt->invoked_flag & INVOKED_PERATOM)) {
    cmpt->compute_peratom();
    cmpt->invoked_flag |= INVOKED_PERATOM;
  }
  return cmpt->vector_atom;
}

int  LammpsInterface::compute_ncols(const char* tag) 
{
  int icompute = find_compute(tag);
  int ncols = lammps_->modify->compute[icompute]->size_peratom_cols;
  if (ncols == 0) ncols = 1; // oddity of lammps, used as flag for scalar
  return ncols;
}

// -----------------------------------------------------------------
//  compute pe/atom interface methods
// -----------------------------------------------------------------
int  LammpsInterface::atomPE_create(void) 
{
  //char * list[3] = {"atcPE","all","pe/atom"};
  char * list[4] = {"atcPE","all","pe/atom","pair"};
  int icompute = lammps_->modify->find_compute(list[0]);
  if (icompute < 0) {
    lammps_->modify->add_compute(3,list);
    icompute = lammps_->modify->find_compute(list[0]);
  }
  if (! atomPE_ ) {
    atomPE_ = lammps_->modify->compute[icompute];
  }
  return icompute;
}

void LammpsInterface::atomPE_init(void)
{
  if (atomPE_) {
    atomPE_->init();
  }
  else {
    throw ATC_Error(0,"no atomPE compute");
  }
}

void LammpsInterface::atomPE_addstep(int step)
{
  atomPE_->addstep(step);
}

double * LammpsInterface::atomPE_compute(void)
{
  if (atomPE_) {
    atomPE_->compute_peratom();
    return atomPE_->vector_atom;
  }
  else {
    return NULL;
  }
}

/* ---------------------------------------------------------------------- */

void LammpsInterface::unwrap_coordinates(int iatom, double* xatom)
{
  double **x = lammps_->atom->x;
  int *mask  = lammps_->atom->mask;
  int *image = lammps_->atom->image;
  int nlocal = lammps_->atom->nlocal;

  double *h   = lammps_->domain->h;
  double xprd = lammps_->domain->xprd;
  double yprd = lammps_->domain->yprd;
  double zprd = lammps_->domain->zprd;
  int xbox,ybox,zbox;

  // for triclinic, need to unwrap current atom coord via h matrix
  // NOTE: Using current box dimensions.
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


// -----------------------------------------------------------------
//  other methods
// -----------------------------------------------------------------

/** Return lammps pointer -- only as a last resort! */
LAMMPS_NS::LAMMPS * LammpsInterface::get_lammps_ptr() { return lammps_; }

}
