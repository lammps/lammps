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

#include "atom.h"
#include <mpi.h>
#include <climits>
#include <cstdlib>
#include <cstring>
#include "style_atom.h"
#include "atom_vec.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "neighbor.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "molecule.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

#ifdef LMP_USER_INTEL
#include "neigh_request.h"
#endif

using namespace LAMMPS_NS;
using namespace MathConst;

#define DELTA 1
#define DELTA_PERATOM 64
#define EPSILON 1.0e-6

enum{DOUBLE,INT,BIGINT};

/* ---------------------------------------------------------------------- */

Atom::Atom(LAMMPS *lmp) : Pointers(lmp)
{
  natoms = 0;
  nlocal = nghost = nmax = 0;
  ntypes = 0;
  nellipsoids = nlines = ntris = nbodies = 0;
  nbondtypes = nangletypes = ndihedraltypes = nimpropertypes = 0;
  nbonds = nangles = ndihedrals = nimpropers = 0;

  firstgroupname = NULL;
  sortfreq = 1000;
  nextsort = 0;
  userbinsize = 0.0;
  maxbin = maxnext = 0;
  binhead = NULL;
  next = permute = NULL;

  // data structure with info on per-atom vectors/arrays

  nperatom = maxperatom = 0;
  peratom = NULL;

  // initialize atom arrays
  // customize by adding new array

  tag = NULL;
  type = mask = NULL;
  image = NULL;
  x = v = f = NULL;

  // charged and dipolar particles

  q = NULL;
  mu = NULL;

  // finite-size particles

  omega = angmom = torque = NULL;
  radius = rmass = NULL;
  ellipsoid = line = tri = body = NULL;

  // molecular systems

  molecule = NULL;
  molindex = molatom = NULL;

  bond_per_atom =  extra_bond_per_atom = 0;
  num_bond = NULL;
  bond_type = NULL;
  bond_atom = NULL;

  angle_per_atom = extra_angle_per_atom = 0;
  num_angle = NULL;
  angle_type = NULL;
  angle_atom1 = angle_atom2 = angle_atom3 = NULL;

  dihedral_per_atom = extra_dihedral_per_atom = 0;
  num_dihedral = NULL;
  dihedral_type = NULL;
  dihedral_atom1 = dihedral_atom2 = dihedral_atom3 = dihedral_atom4 = NULL;

  improper_per_atom = extra_improper_per_atom = 0;
  num_improper = NULL;
  improper_type = NULL;
  improper_atom1 = improper_atom2 = improper_atom3 = improper_atom4 = NULL;

  maxspecial = 1;
  nspecial = NULL;
  special = NULL;

  // PERI package

  vfrac = s0 = NULL;
  x0 = NULL;

  // SPIN package

  sp = fm = fm_long = NULL;

  // USER-EFF and USER-AWPMD packages

  spin = NULL;
  eradius = ervel = erforce = NULL;
  ervelforce = cs = csforce = NULL;
  vforce = NULL;
  etag = NULL;

  // USER-DPD package
  
  uCond = uMech = uChem = uCG = uCGnew = NULL;
  duChem = dpdTheta = NULL;

  // USER-MESO package

  cc = cc_flux = NULL;
  edpd_temp = edpd_flux = edpd_cv = NULL;

  // USER-SMD package

  contact_radius = NULL;
  smd_data_9 = NULL;
  smd_stress = NULL;
  eff_plastic_strain = NULL;
  eff_plastic_strain_rate = NULL;
  damage = NULL;

  // USER-SPH package
  
  rho = drho = e = de = cv = NULL;
  vest = NULL;

  // user-defined molecules

  nmolecule = 0;
  molecules = NULL;

  // custom atom arrays

  nivector = ndvector = 0;
  ivector = NULL;
  dvector = NULL;
  iname = dname = NULL;

  // initialize atom style and array existence flags

  set_atomflag_defaults();

  // initialize peratom data structure

  peratom_create();

  // ntype-length arrays

  mass = NULL;
  mass_setflag = NULL;

  // callback lists & extra restart info

  nextra_grow = nextra_restart = nextra_border = 0;
  extra_grow = extra_restart = extra_border = NULL;
  nextra_grow_max = nextra_restart_max = nextra_border_max = 0;
  nextra_store = 0;
  extra = NULL;

  // default atom ID and mapping values

  tag_enable = 1;
  map_style = map_user = 0;
  map_tag_max = -1;
  map_maxarray = map_nhash = map_nbucket = -1;

  max_same = 0;
  sametag = NULL;
  map_array = NULL;
  map_bucket = NULL;
  map_hash = NULL;

  atom_style = NULL;
  avec = NULL;

  avec_map = new AtomVecCreatorMap();

#define ATOM_CLASS
#define AtomStyle(key,Class) \
  (*avec_map)[#key] = &avec_creator<Class>;
#include "style_atom.h"
#undef AtomStyle
#undef ATOM_CLASS
}

/* ---------------------------------------------------------------------- */

Atom::~Atom()
{
  delete [] atom_style;
  delete avec;
  delete avec_map;

  delete [] firstgroupname;
  memory->destroy(binhead);
  memory->destroy(next);
  memory->destroy(permute);

  // delete peratom data struct

  for (int i = 0; i < nperatom; i++)
    delete [] peratom[i].name;
  memory->sfree(peratom);

  // delete atom arrays
  // customize by adding new array

  memory->destroy(tag);
  memory->destroy(type);
  memory->destroy(mask);
  memory->destroy(image);
  memory->destroy(x);
  memory->destroy(v);
  memory->destroy(f);

  memory->destroy(molecule);
  memory->destroy(molindex);
  memory->destroy(molatom);

  memory->destroy(q);
  memory->destroy(mu);
  memory->destroy(omega);
  memory->destroy(angmom);
  memory->destroy(torque);
  memory->destroy(radius);
  memory->destroy(rmass);
  memory->destroy(ellipsoid);
  memory->destroy(line);
  memory->destroy(tri);
  memory->destroy(body);

  memory->destroy(sp);
  memory->destroy(fm);
  memory->destroy(fm_long);

  memory->destroy(vfrac);
  memory->destroy(s0);
  memory->destroy(x0);

  memory->destroy(spin);
  memory->destroy(eradius);
  memory->destroy(ervel);
  memory->destroy(erforce);
  memory->destroy(ervelforce);
  memory->destroy(cs);
  memory->destroy(csforce);
  memory->destroy(vforce);
  memory->destroy(etag);

  memory->destroy(rho);
  memory->destroy(drho);
  memory->destroy(e);
  memory->destroy(de);
  memory->destroy(cv);
  memory->destroy(vest);

  memory->destroy(contact_radius);
  memory->destroy(smd_data_9);
  memory->destroy(smd_stress);
  memory->destroy(eff_plastic_strain);
  memory->destroy(eff_plastic_strain_rate);
  memory->destroy(damage);

  memory->destroy(dpdTheta);
  memory->destroy(uCond);
  memory->destroy(uMech);
  memory->destroy(uChem);
  memory->destroy(uCG);
  memory->destroy(uCGnew);
  memory->destroy(duChem);

  memory->destroy(cc);
  memory->destroy(cc_flux);
  memory->destroy(edpd_temp);
  memory->destroy(edpd_flux);
  memory->destroy(edpd_cv);

  memory->destroy(nspecial);
  memory->destroy(special);

  memory->destroy(num_bond);
  memory->destroy(bond_type);
  memory->destroy(bond_atom);

  memory->destroy(num_angle);
  memory->destroy(angle_type);
  memory->destroy(angle_atom1);
  memory->destroy(angle_atom2);
  memory->destroy(angle_atom3);

  memory->destroy(num_dihedral);
  memory->destroy(dihedral_type);
  memory->destroy(dihedral_atom1);
  memory->destroy(dihedral_atom2);
  memory->destroy(dihedral_atom3);
  memory->destroy(dihedral_atom4);

  memory->destroy(num_improper);
  memory->destroy(improper_type);
  memory->destroy(improper_atom1);
  memory->destroy(improper_atom2);
  memory->destroy(improper_atom3);
  memory->destroy(improper_atom4);

  // delete custom atom arrays

  for (int i = 0; i < nivector; i++) {
    delete [] iname[i];
    memory->destroy(ivector[i]);
  }
  if (dvector != NULL) {
    for (int i = 0; i < ndvector; i++) {
      delete [] dname[i];
      memory->destroy(dvector[i]);
    }
  }

  memory->sfree(iname);
  memory->sfree(dname);
  memory->sfree(ivector);
  memory->sfree(dvector);

  // delete user-defined molecules

  for (int i = 0; i < nmolecule; i++) delete molecules[i];
  memory->sfree(molecules);

  // delete per-type arrays

  delete [] mass;
  delete [] mass_setflag;

  // delete extra arrays

  memory->destroy(extra_grow);
  memory->destroy(extra_restart);
  memory->destroy(extra_border);
  memory->destroy(extra);

  // delete mapping data structures

  map_delete();
}

/* ----------------------------------------------------------------------
   copy modify settings from old Atom class to current Atom class
------------------------------------------------------------------------- */

void Atom::settings(Atom *old)
{
  tag_enable = old->tag_enable;
  map_user = old->map_user;
  map_style = old->map_style;
  sortfreq = old->sortfreq;
  userbinsize = old->userbinsize;
  if (old->firstgroupname) {
    int n = strlen(old->firstgroupname) + 1;
    firstgroupname = new char[n];
    strcpy(firstgroupname,old->firstgroupname);
  }
}

/* ----------------------------------------------------------------------
   one-time creation of peratom data structure
------------------------------------------------------------------------- */

void Atom::peratom_create()
{
  for (int i = 0; i < nperatom; i++)
    delete [] peratom[i].name;
  memory->sfree(peratom);

  peratom = NULL;
  nperatom = maxperatom = 0;

  // customize: add new peratom variables here, order does not matter
  // register tagint & imageint variables as INT or BIGINT

  int tagintsize = INT;
  if (sizeof(tagint) == 8) tagintsize = BIGINT;
  int imageintsize = INT;
  if (sizeof(imageint) == 8) imageintsize = BIGINT;

  add_peratom("id",&tag,tagintsize,0);
  add_peratom("type",&type,INT,0);
  add_peratom("mask",&mask,INT,0);
  add_peratom("image",&image,imageintsize,0);

  add_peratom("x",&x,DOUBLE,3);
  add_peratom("v",&v,DOUBLE,3);
  add_peratom("f",&f,DOUBLE,3,1);      // set per-thread flag

  add_peratom("rmass",&rmass,DOUBLE,0);
  add_peratom("q",&q,DOUBLE,0);
  add_peratom("mu",&mu,DOUBLE,4);
  add_peratom("mu3",&mu,DOUBLE,3);     // just first 3 values of mu[4]

  // finite size particles

  add_peratom("radius",&radius,DOUBLE,0);
  add_peratom("omega",&omega,DOUBLE,3);
  add_peratom("torque",&torque,DOUBLE,3,1);    // set per-thread flag
  add_peratom("angmom",&angmom,DOUBLE,3);

  add_peratom("ellipsoid",&ellipsoid,INT,0);
  add_peratom("line",&line,INT,0);
  add_peratom("tri",&tri,INT,0);
  add_peratom("body",&body,INT,0);

  // MOLECULE package

  add_peratom("molecule",&molecule,tagintsize,0);
  add_peratom("molindex",&molindex,INT,0);
  add_peratom("molatom",&molatom,INT,0);

  add_peratom("nspecial",&nspecial,INT,3);
  add_peratom_vary("special",&special,tagintsize,&maxspecial,&nspecial,3);

  add_peratom("num_bond",&num_bond,INT,0);
  add_peratom_vary("bond_type",&bond_type,INT,&bond_per_atom,&num_bond);
  add_peratom_vary("bond_atom",&bond_atom,tagintsize,&bond_per_atom,&num_bond);

  add_peratom("num_angle",&num_angle,INT,0);
  add_peratom_vary("angle_type",&angle_type,INT,&angle_per_atom,&num_angle);
  add_peratom_vary("angle_atom1",&angle_atom1,tagintsize,
                   &angle_per_atom,&num_angle);
  add_peratom_vary("angle_atom2",&angle_atom2,tagintsize,
                   &angle_per_atom,&num_angle);
  add_peratom_vary("angle_atom3",&angle_atom3,tagintsize,
                   &angle_per_atom,&num_angle);

  add_peratom("num_dihedral",&num_dihedral,INT,0);
  add_peratom_vary("dihedral_type",&dihedral_type,INT,
                   &dihedral_per_atom,&num_dihedral);
  add_peratom_vary("dihedral_atom1",&dihedral_atom1,tagintsize,
                   &dihedral_per_atom,&num_dihedral);
  add_peratom_vary("dihedral_atom2",&dihedral_atom2,tagintsize,
                   &dihedral_per_atom,&num_dihedral);
  add_peratom_vary("dihedral_atom3",&dihedral_atom3,tagintsize,
                   &dihedral_per_atom,&num_dihedral);
  add_peratom_vary("dihedral_atom4",&dihedral_atom4,tagintsize,
                   &dihedral_per_atom,&num_dihedral);

  add_peratom("num_improper",&num_improper,INT,0);
  add_peratom_vary("improper_type",&improper_type,INT,
                   &improper_per_atom,&num_improper);
  add_peratom_vary("improper_atom1",&improper_atom1,tagintsize,
                   &improper_per_atom,&num_improper);
  add_peratom_vary("improper_atom2",&improper_atom2,tagintsize,
                   &improper_per_atom,&num_improper);
  add_peratom_vary("improper_atom3",&improper_atom3,tagintsize,
                   &improper_per_atom,&num_improper);
  add_peratom_vary("improper_atom4",&improper_atom4,tagintsize,
                   &improper_per_atom,&num_improper);

  // PERI package

  add_peratom("vfrac",&vfrac,DOUBLE,0);
  add_peratom("s0",&s0,DOUBLE,0);
  add_peratom("x0",&x0,DOUBLE,3);

  // SPIN package

  add_peratom("sp",&sp,DOUBLE,4);
  add_peratom("fm",&fm,DOUBLE,3,1);
  add_peratom("fm_long",&fm_long,DOUBLE,3,1);

  // USER-EFF package

  add_peratom("spin",&spin,INT,0);
  add_peratom("eradius",&eradius,DOUBLE,0);
  add_peratom("ervel",&ervel,DOUBLE,0);
  add_peratom("erforce",&erforce,DOUBLE,0,1);     // set per-thread flag

  // USER-AWPMD package

  add_peratom("cs",&cs,DOUBLE,0);
  add_peratom("csforce",&csforce,DOUBLE,0);
  add_peratom("vforce",&vforce,DOUBLE,3);
  add_peratom("ervelforce",&ervelforce,DOUBLE,0);
  add_peratom("etag",&etag,INT,0);

  // USER-DPD package

  add_peratom("dpdTheta",&dpdTheta,DOUBLE,0);
  add_peratom("uCond",&uCond,DOUBLE,0);
  add_peratom("uMech",&uMech,DOUBLE,0);
  add_peratom("uChem",&uChem,DOUBLE,0);
  add_peratom("uCG",&uCG,DOUBLE,0);
  add_peratom("uCGnew",&uCGnew,DOUBLE,0);
  add_peratom("duChem",&duChem,DOUBLE,0);

  // USER-MESO package

  add_peratom("edpd_cv",&edpd_cv,DOUBLE,0);
  add_peratom("edpd_temp",&edpd_temp,DOUBLE,0);
  add_peratom("vest_temp",&vest_temp,DOUBLE,0);
  add_peratom("edpd_flux",&edpd_flux,DOUBLE,0,1);     // set per-thread flag
  add_peratom("cc",&cc,DOUBLE,1);
  add_peratom("cc_flux",&cc_flux,DOUBLE,1,1);         // set per-thread flag

  // USER-SPH package

  add_peratom("rho",&rho,DOUBLE,0);
  add_peratom("drho",&drho,DOUBLE,0,1);               // set per-thread flag
  add_peratom("e",&e,DOUBLE,0);
  add_peratom("de",&de,DOUBLE,0,1);                   // set per-thread flag
  add_peratom("vest",&vest,DOUBLE,3);
  add_peratom("cv",&cv,DOUBLE,0);

  // USER-SMD package

  add_peratom("contact_radius",&contact_radius,DOUBLE,0);
  add_peratom("smd_data_9",&smd_data_9,DOUBLE,1);
  add_peratom("smd_stress",&smd_stress,DOUBLE,1);
  add_peratom("eff_plastic_strain",&eff_plastic_strain,DOUBLE,0);
  add_peratom("eff_plastic_strain_rate",&eff_plastic_strain_rate,DOUBLE,0);
  add_peratom("damage",&damage,DOUBLE,0);
}

/* ----------------------------------------------------------------------
   add info for a single per-atom vector/array to PerAtom data struct
   cols = 0: per-atom vector 
   cols = N: static per-atom array with N columns
   use add_peratom_vary() when column count varies per atom
------------------------------------------------------------------------- */

void Atom::add_peratom(const char *name, void *address, 
                       int datatype, int cols, int threadflag)
{
  if (nperatom == maxperatom) {
    maxperatom += DELTA_PERATOM;
    peratom = (PerAtom *) 
      memory->srealloc(peratom,maxperatom*sizeof(PerAtom),"atom:peratom");
  }

  int n = strlen(name) + 1;
  peratom[nperatom].name = new char[n];
  strcpy(peratom[nperatom].name,name);
  peratom[nperatom].address = address;
  peratom[nperatom].datatype = datatype;
  peratom[nperatom].cols = cols;
  peratom[nperatom].threadflag = threadflag;
  peratom[nperatom].address_length = NULL;

  nperatom++;
}

/* ----------------------------------------------------------------------
   change the column count fof an existing peratom array entry
   allows atom_style to specify column count as an argument
   see atom_style tdpd as an example
------------------------------------------------------------------------- */

void Atom::add_peratom_change_columns(const char *name, int cols)
{
  int i;
  for (int i = 0; i < nperatom; i++)
    if (strcmp(name,peratom[i].name) == 0) peratom[i].cols = cols;
  if (i == nperatom) 
    error->all(FLERR,"Could not find name of peratom array for column change");
}

/* ----------------------------------------------------------------------
   add info for a single per-atom array to PerAtom data struct
   cols = address of int variable with max columns per atom
   for collength = 0:
     length = address of peratom vector with column count per atom
     e.g. num_bond
   for collength = N: 
     length = address of peratom array with column count per atom
     collength = index of column (1 to N) in peratom array with count
     e.g. nspecial
------------------------------------------------------------------------- */

void Atom::add_peratom_vary(const char *name, void *address, 
                            int datatype, int *cols, void *length, int collength)
{
  if (nperatom == maxperatom) {
    maxperatom += DELTA_PERATOM;
    peratom = (PerAtom *) 
      memory->srealloc(peratom,maxperatom*sizeof(PerAtom),"atom:peratom");
  }

  int n = strlen(name) + 1;
  peratom[nperatom].name = new char[n];
  strcpy(peratom[nperatom].name,name);
  peratom[nperatom].address = address;
  peratom[nperatom].datatype = datatype;
  peratom[nperatom].cols = -1;
  peratom[nperatom].threadflag = 0;
  peratom[nperatom].address_maxcols = cols;
  peratom[nperatom].address_length = length;
  peratom[nperatom].collength = collength;

  nperatom++;
}

/* ----------------------------------------------------------------------
   add info for a single per-atom array to PerAtom data struct
   customize by adding new flag, identical list as atom.h 2nd customization
------------------------------------------------------------------------- */

void Atom::set_atomflag_defaults()
{
  sphere_flag = ellipsoid_flag = line_flag = tri_flag = body_flag = 0;
  peri_flag = electron_flag = 0;
  wavepacket_flag = sph_flag = 0;
  molecule_flag = molindex_flag = molatom_flag = 0;
  q_flag = mu_flag = 0;
  rmass_flag = radius_flag = omega_flag = torque_flag = angmom_flag = 0;
  vfrac_flag = spin_flag = eradius_flag = ervel_flag = erforce_flag = 0;
  cs_flag = csforce_flag = vforce_flag = ervelforce_flag = etag_flag = 0;
  rho_flag = e_flag = cv_flag = vest_flag = 0;
  dpd_flag = edpd_flag = tdpd_flag = 0;
  sp_flag = 0;
  x0_flag = 0;
  smd_flag = damage_flag = 0;
  contact_radius_flag = smd_data_9_flag = smd_stress_flag = 0;
  eff_plastic_strain_flag = eff_plastic_strain_rate_flag = 0;

  pdscale = 1.0;
}

/* ----------------------------------------------------------------------
   create an AtomVec style
   called from lammps.cpp, input script, restart file, replicate
------------------------------------------------------------------------- */

void Atom::create_avec(const char *style, int narg, char **arg, int trysuffix)
{
  delete [] atom_style;
  if (avec) delete avec;
  atom_style = NULL;
  avec = NULL;

  // unset atom style and array existence flags
  // may have been set by old avec

  set_atomflag_defaults();

  // create instance of AtomVec
  // use grow() to initialize atom-based arrays to length 1
  //   so that x[0][0] can always be referenced even if proc has no atoms

  int sflag;
  avec = new_avec(style,trysuffix,sflag);
  avec->store_args(narg,arg);
  avec->process_args(narg,arg);
  avec->grow(1);

  if (sflag) {
    char estyle[256];
    if (sflag == 1) snprintf(estyle,256,"%s/%s",style,lmp->suffix);
    else snprintf(estyle,256,"%s/%s",style,lmp->suffix2);
    int n = strlen(estyle) + 1;
    atom_style = new char[n];
    strcpy(atom_style,estyle);
  } else {
    int n = strlen(style) + 1;
    atom_style = new char[n];
    strcpy(atom_style,style);
  }

  // if molecular system:
  // atom IDs must be defined
  // force atom map to be created
  // map style will be reset to array vs hash to by map_init()

  molecular = avec->molecular;
  if (molecular && tag_enable == 0)
    error->all(FLERR,"Atom IDs must be used for molecular systems");
  if (molecular) map_style = 3;
}

/* ----------------------------------------------------------------------
   generate an AtomVec class, first with suffix appended
------------------------------------------------------------------------- */

AtomVec *Atom::new_avec(const char *style, int trysuffix, int &sflag)
{
  if (trysuffix && lmp->suffix_enable) {
    if (lmp->suffix) {
      sflag = 1;
      char estyle[256];
      snprintf(estyle,256,"%s/%s",style,lmp->suffix);
      if (avec_map->find(estyle) != avec_map->end()) {
        AtomVecCreator avec_creator = (*avec_map)[estyle];
        return avec_creator(lmp);
      }
    }

    if (lmp->suffix2) {
      sflag = 2;
      char estyle[256];
      snprintf(estyle,256,"%s/%s",style,lmp->suffix2);
      if (avec_map->find(estyle) != avec_map->end()) {
        AtomVecCreator avec_creator = (*avec_map)[estyle];
        return avec_creator(lmp);
      }
    }
  }

  sflag = 0;
  if (avec_map->find(style) != avec_map->end()) {
    AtomVecCreator avec_creator = (*avec_map)[style];
    return avec_creator(lmp);
  }

  error->all(FLERR,utils::check_packages_for_style("atom",style,lmp).c_str());
  return NULL;
}

/* ----------------------------------------------------------------------
   one instance per AtomVec style in style_atom.h
------------------------------------------------------------------------- */

template <typename T>
AtomVec *Atom::avec_creator(LAMMPS *lmp)
{
  return new T(lmp);
}


/* ---------------------------------------------------------------------- */

void Atom::init()
{
  // delete extra array since it doesn't persist past first run

  if (nextra_store) {
    memory->destroy(extra);
    extra = NULL;
    nextra_store = 0;
  }

  // check arrays that are atom type in length

  check_mass(FLERR);

  // setup of firstgroup

  if (firstgroupname) {
    firstgroup = group->find(firstgroupname);
    if (firstgroup < 0)
      error->all(FLERR,"Could not find atom_modify first group ID");
  } else firstgroup = -1;

  // init AtomVec

  avec->init();
}

/* ---------------------------------------------------------------------- */

void Atom::setup()
{
  // setup bins for sorting
  // cannot do this in init() because uses neighbor cutoff

  if (sortfreq > 0) setup_sort_bins();
}

/* ----------------------------------------------------------------------
   return ptr to AtomVec class if matches style or to matching hybrid sub-class
   return NULL if no match
------------------------------------------------------------------------- */

AtomVec *Atom::style_match(const char *style)
{
  if (strcmp(atom_style,style) == 0) return avec;
  else if (strcmp(atom_style,"hybrid") == 0) {
    AtomVecHybrid *avec_hybrid = (AtomVecHybrid *) avec;
    for (int i = 0; i < avec_hybrid->nstyles; i++)
      if (strcmp(avec_hybrid->keywords[i],style) == 0)
        return avec_hybrid->styles[i];
  }
  return NULL;
}

/* ----------------------------------------------------------------------
   modify parameters of the atom style
   some options can only be invoked before simulation box is defined
   first and sort options cannot be used together
------------------------------------------------------------------------- */

void Atom::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal atom_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"id") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal atom_modify command");
      if (domain->box_exist)
        error->all(FLERR,
                   "Atom_modify id command after simulation box is defined");
      if (strcmp(arg[iarg+1],"yes") == 0) tag_enable = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) tag_enable = 0;
      else error->all(FLERR,"Illegal atom_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"map") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal atom_modify command");
      if (domain->box_exist)
        error->all(FLERR,
                   "Atom_modify map command after simulation box is defined");
      if (strcmp(arg[iarg+1],"array") == 0) map_user = 1;
      else if (strcmp(arg[iarg+1],"hash") == 0) map_user = 2;
      else if (strcmp(arg[iarg+1],"yes") == 0) map_user = 3;
      else error->all(FLERR,"Illegal atom_modify command");
      map_style = map_user;
      iarg += 2;
    } else if (strcmp(arg[iarg],"first") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal atom_modify command");
      if (strcmp(arg[iarg+1],"all") == 0) {
        delete [] firstgroupname;
        firstgroupname = NULL;
      } else {
        int n = strlen(arg[iarg+1]) + 1;
        firstgroupname = new char[n];
        strcpy(firstgroupname,arg[iarg+1]);
        sortfreq = 0;
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"sort") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal atom_modify command");
      sortfreq = force->inumeric(FLERR,arg[iarg+1]);
      userbinsize = force->numeric(FLERR,arg[iarg+2]);
      if (sortfreq < 0 || userbinsize < 0.0)
        error->all(FLERR,"Illegal atom_modify command");
      if (sortfreq >= 0 && firstgroupname)
        error->all(FLERR,"Atom_modify sort and first options "
                   "cannot be used together");
      iarg += 3;
    } else error->all(FLERR,"Illegal atom_modify command");
  }
}

/* ----------------------------------------------------------------------
   check that atom IDs are valid
   error if any atom ID < 0 or atom ID = MAXTAGINT
   if any atom ID > 0, error if any atom ID == 0
   if any atom ID > 0, error if tag_enable = 0
   if all atom IDs = 0, tag_enable must be 0
   if max atom IDs < natoms, must be duplicates
   OK if max atom IDs > natoms
   NOTE: not fully checking that atom IDs are unique
------------------------------------------------------------------------- */

void Atom::tag_check()
{
  tagint min = MAXTAGINT;
  tagint max = 0;

  for (int i = 0; i < nlocal; i++) {
    min = MIN(min,tag[i]);
    max = MAX(max,tag[i]);
  }

  tagint minall,maxall;
  MPI_Allreduce(&min,&minall,1,MPI_LMP_TAGINT,MPI_MIN,world);
  MPI_Allreduce(&max,&maxall,1,MPI_LMP_TAGINT,MPI_MAX,world);

  if (minall < 0) error->all(FLERR,"One or more Atom IDs is negative");
  if (maxall >= MAXTAGINT) error->all(FLERR,"One or more atom IDs is too big");
  if (maxall > 0 && minall == 0)
    error->all(FLERR,"One or more atom IDs is zero");
  if (maxall > 0 && tag_enable == 0)
    error->all(FLERR,"Non-zero atom IDs with atom_modify id = no");
  if (maxall == 0 && natoms && tag_enable)
    error->all(FLERR,"All atom IDs = 0 but atom_modify id = yes");
  if (tag_enable && maxall < natoms)
    error->all(FLERR,"Duplicate atom IDs exist");
}

/* ----------------------------------------------------------------------
   add unique tags to any atoms with tag = 0
   new tags are grouped by proc and start after max current tag
   called after creating new atoms
   error if new tags will exceed MAXTAGINT
------------------------------------------------------------------------- */

void Atom::tag_extend()
{
  // maxtag_all = max tag for all atoms

  tagint maxtag = 0;
  for (int i = 0; i < nlocal; i++) maxtag = MAX(maxtag,tag[i]);
  tagint maxtag_all;
  MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);

  // DEBUG: useful for generating 64-bit IDs even for small systems
  // use only when LAMMPS is compiled with BIGBIG

  //maxtag_all += 1000000000000;

  // notag = # of atoms I own with no tag (tag = 0)
  // notag_sum = # of total atoms on procs <= me with no tag

  bigint notag = 0;
  for (int i = 0; i < nlocal; i++) if (tag[i] == 0) notag++;

  bigint notag_total;
  MPI_Allreduce(&notag,&notag_total,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (notag_total >= MAXTAGINT)
    error->all(FLERR,"New atom IDs exceed maximum allowed ID");

  bigint notag_sum;
  MPI_Scan(&notag,&notag_sum,1,MPI_LMP_BIGINT,MPI_SUM,world);

  // itag = 1st new tag that my untagged atoms should use

  tagint itag = maxtag_all + notag_sum - notag + 1;
  for (int i = 0; i < nlocal; i++) if (tag[i] == 0) tag[i] = itag++;
}

/* ----------------------------------------------------------------------
   check that atom IDs span range from 1 to Natoms inclusive
   return 0 if mintag != 1 or maxtag != Natoms
   return 1 if OK
   doesn't actually check if all tag values are used
------------------------------------------------------------------------- */

int Atom::tag_consecutive()
{
  tagint idmin = MAXTAGINT;
  tagint idmax = 0;

  for (int i = 0; i < nlocal; i++) {
    idmin = MIN(idmin,tag[i]);
    idmax = MAX(idmax,tag[i]);
  }
  tagint idminall,idmaxall;
  MPI_Allreduce(&idmin,&idminall,1,MPI_LMP_TAGINT,MPI_MIN,world);
  MPI_Allreduce(&idmax,&idmaxall,1,MPI_LMP_TAGINT,MPI_MAX,world);

  if (idminall != 1 || idmaxall != natoms) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   check that bonus data settings are valid
   error if number of atoms with ellipsoid/line/tri/body flags
   are consistent with global setting.
------------------------------------------------------------------------- */

void Atom::bonus_check()
{
  bigint local_ellipsoids = 0, local_lines = 0, local_tris = 0;
  bigint local_bodies = 0, num_global;

  for (int i = 0; i < nlocal; ++i) {
    if (ellipsoid && (ellipsoid[i] >=0)) ++local_ellipsoids;
    if (line && (line[i] >=0)) ++local_lines;
    if (tri && (tri[i] >=0)) ++local_tris;
    if (body && (body[i] >=0)) ++local_bodies;
  }

  MPI_Allreduce(&local_ellipsoids,&num_global,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (nellipsoids != num_global)
    error->all(FLERR,"Inconsistent 'ellipsoids' header value and number of "
               "atoms with enabled ellipsoid flags");

  MPI_Allreduce(&local_lines,&num_global,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (nlines != num_global)
    error->all(FLERR,"Inconsistent 'lines' header value and number of "
               "atoms with enabled line flags");

  MPI_Allreduce(&local_tris,&num_global,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (ntris != num_global)
    error->all(FLERR,"Inconsistent 'tris' header value and number of "
               "atoms with enabled tri flags");

  MPI_Allreduce(&local_bodies,&num_global,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (nbodies != num_global)
    error->all(FLERR,"Inconsistent 'bodies' header value and number of "
               "atoms with enabled body flags");
}

/* ----------------------------------------------------------------------
   count and return words in a single line
   make copy of line before using strtok so as not to change line
   trim anything from '#' onward
------------------------------------------------------------------------- */

int Atom::count_words(const char *line)
{
  int n = strlen(line) + 1;
  char *copy;
  memory->create(copy,n,"atom:copy");
  strcpy(copy,line);

  char *ptr;
  if ((ptr = strchr(copy,'#'))) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->destroy(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  memory->destroy(copy);
  return n;
}

/* ----------------------------------------------------------------------
   count and return words in a single line using provided copy buf
   make copy of line before using strtok so as not to change line
   trim anything from '#' onward
------------------------------------------------------------------------- */

int Atom::count_words(const char *line, char *copy)
{
  strcpy(copy,line);

  char *ptr;
  if ((ptr = strchr(copy,'#'))) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->destroy(copy);
    return 0;
  }
  int n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  return n;
}

/* ----------------------------------------------------------------------
   deallocate molecular topology arrays
   done before realloc with (possibly) new 2nd dimension set to
     correctly initialized per-atom values, e.g. bond_per_atom
   needs to be called whenever 2nd dimensions are changed
     and these arrays are already pre-allocated,
     e.g. due to grow(1) in create_avec()
------------------------------------------------------------------------- */

void Atom::deallocate_topology()
{
  memory->destroy(atom->bond_type);
  memory->destroy(atom->bond_atom);
  atom->bond_type = NULL;
  atom->bond_atom = NULL;

  memory->destroy(atom->angle_type);
  memory->destroy(atom->angle_atom1);
  memory->destroy(atom->angle_atom2);
  memory->destroy(atom->angle_atom3);
  atom->angle_type = NULL;
  atom->angle_atom1 = atom->angle_atom2 = atom->angle_atom3 = NULL;

  memory->destroy(atom->dihedral_type);
  memory->destroy(atom->dihedral_atom1);
  memory->destroy(atom->dihedral_atom2);
  memory->destroy(atom->dihedral_atom3);
  memory->destroy(atom->dihedral_atom4);
  atom->dihedral_type = NULL;
  atom->dihedral_atom1 = atom->dihedral_atom2 =
    atom->dihedral_atom3 = atom->dihedral_atom4 = NULL;

  memory->destroy(atom->improper_type);
  memory->destroy(atom->improper_atom1);
  memory->destroy(atom->improper_atom2);
  memory->destroy(atom->improper_atom3);
  memory->destroy(atom->improper_atom4);
  atom->improper_type = NULL;
  atom->improper_atom1 = atom->improper_atom2 =
    atom->improper_atom3 = atom->improper_atom4 = NULL;
}

/* ----------------------------------------------------------------------
   unpack N lines from Atom section of data file
   call style-specific routine to parse line
------------------------------------------------------------------------- */

void Atom::data_atoms(int n, char *buf, tagint id_offset, tagint mol_offset,
                      int type_offset, int shiftflag, double *shift)
{
  int m,xptr,iptr;
  imageint imagedata;
  double xdata[3],lamda[3];
  double *coord;
  char *next;

  next = strchr(buf,'\n');
  *next = '\0';
  int nwords = count_words(buf);
  *next = '\n';

  if (nwords != avec->size_data_atom && nwords != avec->size_data_atom + 3)
    error->all(FLERR,"Incorrect atom format in data file");

  char **values = new char*[nwords];

  // set bounds for my proc
  // if periodic and I am lo/hi proc, adjust bounds by EPSILON
  // insures all data atoms will be owned even with round-off

  int triclinic = domain->triclinic;

  double epsilon[3];
  if (triclinic) epsilon[0] = epsilon[1] = epsilon[2] = EPSILON;
  else {
    epsilon[0] = domain->prd[0] * EPSILON;
    epsilon[1] = domain->prd[1] * EPSILON;
    epsilon[2] = domain->prd[2] * EPSILON;
  }

  double sublo[3],subhi[3];
  if (triclinic == 0) {
    sublo[0] = domain->sublo[0]; subhi[0] = domain->subhi[0];
    sublo[1] = domain->sublo[1]; subhi[1] = domain->subhi[1];
    sublo[2] = domain->sublo[2]; subhi[2] = domain->subhi[2];
  } else {
    sublo[0] = domain->sublo_lamda[0]; subhi[0] = domain->subhi_lamda[0];
    sublo[1] = domain->sublo_lamda[1]; subhi[1] = domain->subhi_lamda[1];
    sublo[2] = domain->sublo_lamda[2]; subhi[2] = domain->subhi_lamda[2];
  }

  if (comm->layout != Comm::LAYOUT_TILED) {
    if (domain->xperiodic) {
      if (comm->myloc[0] == 0) sublo[0] -= epsilon[0];
      if (comm->myloc[0] == comm->procgrid[0]-1) subhi[0] += epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->myloc[1] == 0) sublo[1] -= epsilon[1];
      if (comm->myloc[1] == comm->procgrid[1]-1) subhi[1] += epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->myloc[2] == 0) sublo[2] -= epsilon[2];
      if (comm->myloc[2] == comm->procgrid[2]-1) subhi[2] += epsilon[2];
    }

  } else {
    if (domain->xperiodic) {
      if (comm->mysplit[0][0] == 0.0) sublo[0] -= epsilon[0];
      if (comm->mysplit[0][1] == 1.0) subhi[0] += epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->mysplit[1][0] == 0.0) sublo[1] -= epsilon[1];
      if (comm->mysplit[1][1] == 1.0) subhi[1] += epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->mysplit[2][0] == 0.0) sublo[2] -= epsilon[2];
      if (comm->mysplit[2][1] == 1.0) subhi[2] += epsilon[2];
    }
  }

  // xptr = which word in line starts xyz coords
  // iptr = which word in line starts ix,iy,iz image flags

  xptr = avec->xcol_data - 1;
  int imageflag = 0;
  if (nwords > avec->size_data_atom) imageflag = 1;
  if (imageflag) iptr = nwords - 3;

  // loop over lines of atom data
  // tokenize the line into values
  // extract xyz coords and image flags
  // remap atom into simulation box
  // if atom is in my sub-domain, unpack its values

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');

    values[0] = strtok(buf," \t\n\r\f");
    if (values[0] == NULL)
      error->all(FLERR,"Incorrect atom format in data file");
    for (m = 1; m < nwords; m++) {
      values[m] = strtok(NULL," \t\n\r\f");
      if (values[m] == NULL)
        error->all(FLERR,"Incorrect atom format in data file");
    }

    if (imageflag)
      imagedata = ((imageint) (atoi(values[iptr]) + IMGMAX) & IMGMASK) |
        (((imageint) (atoi(values[iptr+1]) + IMGMAX) & IMGMASK) << IMGBITS) |
        (((imageint) (atoi(values[iptr+2]) + IMGMAX) & IMGMASK) << IMG2BITS);
    else imagedata = ((imageint) IMGMAX << IMG2BITS) |
           ((imageint) IMGMAX << IMGBITS) | IMGMAX;

    xdata[0] = atof(values[xptr]);
    xdata[1] = atof(values[xptr+1]);
    xdata[2] = atof(values[xptr+2]);
    if (shiftflag) {
      xdata[0] += shift[0];
      xdata[1] += shift[1];
      xdata[2] += shift[2];
    }

    domain->remap(xdata,imagedata);
    if (triclinic) {
      domain->x2lamda(xdata,lamda);
      coord = lamda;
    } else coord = xdata;

    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
        coord[1] >= sublo[1] && coord[1] < subhi[1] &&
        coord[2] >= sublo[2] && coord[2] < subhi[2]) {
      avec->data_atom(xdata,imagedata,values);
      if (id_offset) tag[nlocal-1] += id_offset;
      if (mol_offset) molecule[nlocal-1] += mol_offset;
      if (type_offset) {
        type[nlocal-1] += type_offset;
        if (type[nlocal-1] > ntypes)
          error->one(FLERR,"Invalid atom type in Atoms section of data file");
      }
    }

    buf = next + 1;
  }

  delete [] values;
}

/* ----------------------------------------------------------------------
   unpack N lines from Velocity section of data file
   check that atom IDs are > 0 and <= map_tag_max
   call style-specific routine to parse line
------------------------------------------------------------------------- */

void Atom::data_vels(int n, char *buf, tagint id_offset)
{
  int j,m;
  tagint tagdata;
  char *next;

  next = strchr(buf,'\n');
  *next = '\0';
  int nwords = count_words(buf);
  *next = '\n';

  if (nwords != avec->size_data_vel)
    error->all(FLERR,"Incorrect velocity format in data file");

  char **values = new char*[nwords];

  // loop over lines of atom velocities
  // tokenize the line into values
  // if I own atom tag, unpack its values

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');

    values[0] = strtok(buf," \t\n\r\f");
    for (j = 1; j < nwords; j++)
      values[j] = strtok(NULL," \t\n\r\f");

    tagdata = ATOTAGINT(values[0]) + id_offset;
    if (tagdata <= 0 || tagdata > map_tag_max)
      error->one(FLERR,"Invalid atom ID in Velocities section of data file");
    if ((m = map(tagdata)) >= 0) avec->data_vel(m,&values[1]);

    buf = next + 1;
  }

  delete [] values;
}

/* ----------------------------------------------------------------------
   process N bonds read into buf from data files
   if count is non-NULL, just count bonds per atom
   else store them with atoms
   check that atom IDs are > 0 and <= map_tag_max
------------------------------------------------------------------------- */

void Atom::data_bonds(int n, char *buf, int *count, tagint id_offset,
                      int type_offset)
{
  int m,tmp,itype,rv;
  tagint atom1,atom2;
  char *next;
  int newton_bond = force->newton_bond;

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    rv = sscanf(buf,"%d %d " TAGINT_FORMAT " " TAGINT_FORMAT,
                &tmp,&itype,&atom1,&atom2);
    if (rv != 4)
      error->one(FLERR,"Incorrect format of Bonds section in data file");
    if (id_offset) {
      atom1 += id_offset;
      atom2 += id_offset;
    }
    itype += type_offset;

    if ((atom1 <= 0) || (atom1 > map_tag_max) ||
        (atom2 <= 0) || (atom2 > map_tag_max) || (atom1 == atom2))
      error->one(FLERR,"Invalid atom ID in Bonds section of data file");
    if (itype <= 0 || itype > nbondtypes)
      error->one(FLERR,"Invalid bond type in Bonds section of data file");
    if ((m = map(atom1)) >= 0) {
      if (count) count[m]++;
      else {
        bond_type[m][num_bond[m]] = itype;
        bond_atom[m][num_bond[m]] = atom2;
        num_bond[m]++;
      }
    }
    if (newton_bond == 0) {
      if ((m = map(atom2)) >= 0) {
        if (count) count[m]++;
        else {
          bond_type[m][num_bond[m]] = itype;
          bond_atom[m][num_bond[m]] = atom1;
          num_bond[m]++;
        }
      }
    }
    buf = next + 1;
  }
}

/* ----------------------------------------------------------------------
   process N angles read into buf from data files
   if count is non-NULL, just count angles per atom
   else store them with atoms
   check that atom IDs are > 0 and <= map_tag_max
------------------------------------------------------------------------- */

void Atom::data_angles(int n, char *buf, int *count, tagint id_offset,
                       int type_offset)
{
  int m,tmp,itype,rv;
  tagint atom1,atom2,atom3;
  char *next;
  int newton_bond = force->newton_bond;

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    rv = sscanf(buf,"%d %d " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT,
                &tmp,&itype,&atom1,&atom2,&atom3);
    if (rv != 5)
      error->one(FLERR,"Incorrect format of Angles section in data file");
    if (id_offset) {
      atom1 += id_offset;
      atom2 += id_offset;
      atom3 += id_offset;
    }
    itype += type_offset;

    if ((atom1 <= 0) || (atom1 > map_tag_max) ||
        (atom2 <= 0) || (atom2 > map_tag_max) ||
        (atom3 <= 0) || (atom3 > map_tag_max) ||
        (atom1 == atom2) || (atom1 == atom3) || (atom2 == atom3))
      error->one(FLERR,"Invalid atom ID in Angles section of data file");
    if (itype <= 0 || itype > nangletypes)
      error->one(FLERR,"Invalid angle type in Angles section of data file");
    if ((m = map(atom2)) >= 0) {
      if (count) count[m]++;
      else {
        angle_type[m][num_angle[m]] = itype;
        angle_atom1[m][num_angle[m]] = atom1;
        angle_atom2[m][num_angle[m]] = atom2;
        angle_atom3[m][num_angle[m]] = atom3;
        num_angle[m]++;
      }
    }
    if (newton_bond == 0) {
      if ((m = map(atom1)) >= 0) {
        if (count) count[m]++;
        else {
          angle_type[m][num_angle[m]] = itype;
          angle_atom1[m][num_angle[m]] = atom1;
          angle_atom2[m][num_angle[m]] = atom2;
          angle_atom3[m][num_angle[m]] = atom3;
          num_angle[m]++;
        }
      }
      if ((m = map(atom3)) >= 0) {
        if (count) count[m]++;
        else {
          angle_type[m][num_angle[m]] = itype;
          angle_atom1[m][num_angle[m]] = atom1;
          angle_atom2[m][num_angle[m]] = atom2;
          angle_atom3[m][num_angle[m]] = atom3;
          num_angle[m]++;
        }
      }
    }
    buf = next + 1;
  }
}

/* ----------------------------------------------------------------------
   process N dihedrals read into buf from data files
   if count is non-NULL, just count diihedrals per atom
   else store them with atoms
   check that atom IDs are > 0 and <= map_tag_max
------------------------------------------------------------------------- */

void Atom::data_dihedrals(int n, char *buf, int *count, tagint id_offset,
                          int type_offset)
{
  int m,tmp,itype,rv;
  tagint atom1,atom2,atom3,atom4;
  char *next;
  int newton_bond = force->newton_bond;

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    rv = sscanf(buf,"%d %d " TAGINT_FORMAT " " TAGINT_FORMAT
                " " TAGINT_FORMAT " " TAGINT_FORMAT,
                &tmp,&itype,&atom1,&atom2,&atom3,&atom4);
    if (rv != 6)
      error->one(FLERR,"Incorrect format of Dihedrals section in data file");
    if (id_offset) {
      atom1 += id_offset;
      atom2 += id_offset;
      atom3 += id_offset;
      atom4 += id_offset;
    }
    itype += type_offset;

    if ((atom1 <= 0) || (atom1 > map_tag_max) ||
        (atom2 <= 0) || (atom2 > map_tag_max) ||
        (atom3 <= 0) || (atom3 > map_tag_max) ||
        (atom4 <= 0) || (atom4 > map_tag_max) ||
        (atom1 == atom2) || (atom1 == atom3) || (atom1 == atom4) ||
        (atom2 == atom3) || (atom2 == atom4) || (atom3 == atom4))
      error->one(FLERR,"Invalid atom ID in Dihedrals section of data file");
    if (itype <= 0 || itype > ndihedraltypes)
      error->one(FLERR,
                 "Invalid dihedral type in Dihedrals section of data file");
    if ((m = map(atom2)) >= 0) {
      if (count) count[m]++;
      else {
        dihedral_type[m][num_dihedral[m]] = itype;
        dihedral_atom1[m][num_dihedral[m]] = atom1;
        dihedral_atom2[m][num_dihedral[m]] = atom2;
        dihedral_atom3[m][num_dihedral[m]] = atom3;
        dihedral_atom4[m][num_dihedral[m]] = atom4;
        num_dihedral[m]++;
      }
    }
    if (newton_bond == 0) {
      if ((m = map(atom1)) >= 0) {
        if (count) count[m]++;
        else {
          dihedral_type[m][num_dihedral[m]] = itype;
          dihedral_atom1[m][num_dihedral[m]] = atom1;
          dihedral_atom2[m][num_dihedral[m]] = atom2;
          dihedral_atom3[m][num_dihedral[m]] = atom3;
          dihedral_atom4[m][num_dihedral[m]] = atom4;
          num_dihedral[m]++;
        }
      }
      if ((m = map(atom3)) >= 0) {
        if (count) count[m]++;
        else {
          dihedral_type[m][num_dihedral[m]] = itype;
          dihedral_atom1[m][num_dihedral[m]] = atom1;
          dihedral_atom2[m][num_dihedral[m]] = atom2;
          dihedral_atom3[m][num_dihedral[m]] = atom3;
          dihedral_atom4[m][num_dihedral[m]] = atom4;
          num_dihedral[m]++;
        }
      }
      if ((m = map(atom4)) >= 0) {
        if (count) count[m]++;
        else {
          dihedral_type[m][num_dihedral[m]] = itype;
          dihedral_atom1[m][num_dihedral[m]] = atom1;
          dihedral_atom2[m][num_dihedral[m]] = atom2;
          dihedral_atom3[m][num_dihedral[m]] = atom3;
          dihedral_atom4[m][num_dihedral[m]] = atom4;
          num_dihedral[m]++;
        }
      }
    }
    buf = next + 1;
  }
}

/* ----------------------------------------------------------------------
   process N impropers read into buf from data files
   if count is non-NULL, just count impropers per atom
   else store them with atoms
   check that atom IDs are > 0 and <= map_tag_max
------------------------------------------------------------------------- */

void Atom::data_impropers(int n, char *buf, int *count, tagint id_offset,
                          int type_offset)
{
  int m,tmp,itype,rv;
  tagint atom1,atom2,atom3,atom4;
  char *next;
  int newton_bond = force->newton_bond;

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    rv = sscanf(buf,"%d %d "
                TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT,
                &tmp,&itype,&atom1,&atom2,&atom3,&atom4);
    if (rv != 6)
      error->one(FLERR,"Incorrect format of Impropers section in data file");
    if (id_offset) {
      atom1 += id_offset;
      atom2 += id_offset;
      atom3 += id_offset;
      atom4 += id_offset;
    }
    itype += type_offset;

    if ((atom1 <= 0) || (atom1 > map_tag_max) ||
        (atom2 <= 0) || (atom2 > map_tag_max) ||
        (atom3 <= 0) || (atom3 > map_tag_max) ||
        (atom4 <= 0) || (atom4 > map_tag_max) ||
        (atom1 == atom2) || (atom1 == atom3) || (atom1 == atom4) ||
        (atom2 == atom3) || (atom2 == atom4) || (atom3 == atom4))
      error->one(FLERR,"Invalid atom ID in Impropers section of data file");
    if (itype <= 0 || itype > nimpropertypes)
      error->one(FLERR,
                 "Invalid improper type in Impropers section of data file");
    if ((m = map(atom2)) >= 0) {
      if (count) count[m]++;
      else {
        improper_type[m][num_improper[m]] = itype;
        improper_atom1[m][num_improper[m]] = atom1;
        improper_atom2[m][num_improper[m]] = atom2;
        improper_atom3[m][num_improper[m]] = atom3;
        improper_atom4[m][num_improper[m]] = atom4;
        num_improper[m]++;
      }
    }
    if (newton_bond == 0) {
      if ((m = map(atom1)) >= 0) {
        if (count) count[m]++;
        else {
          improper_type[m][num_improper[m]] = itype;
          improper_atom1[m][num_improper[m]] = atom1;
          improper_atom2[m][num_improper[m]] = atom2;
          improper_atom3[m][num_improper[m]] = atom3;
          improper_atom4[m][num_improper[m]] = atom4;
          num_improper[m]++;
        }
      }
      if ((m = map(atom3)) >= 0) {
        if (count) count[m]++;
        else {
          improper_type[m][num_improper[m]] = itype;
          improper_atom1[m][num_improper[m]] = atom1;
          improper_atom2[m][num_improper[m]] = atom2;
          improper_atom3[m][num_improper[m]] = atom3;
          improper_atom4[m][num_improper[m]] = atom4;
          num_improper[m]++;
        }
      }
      if ((m = map(atom4)) >= 0) {
        if (count) count[m]++;
        else {
          improper_type[m][num_improper[m]] = itype;
          improper_atom1[m][num_improper[m]] = atom1;
          improper_atom2[m][num_improper[m]] = atom2;
          improper_atom3[m][num_improper[m]] = atom3;
          improper_atom4[m][num_improper[m]] = atom4;
          num_improper[m]++;
        }
      }
    }
    buf = next + 1;
  }
}

/* ----------------------------------------------------------------------
   unpack N lines from atom-style specific bonus section of data file
   check that atom IDs are > 0 and <= map_tag_max
   call style-specific routine to parse line
------------------------------------------------------------------------- */

void Atom::data_bonus(int n, char *buf, AtomVec *avec_bonus, tagint id_offset)
{
  int j,m,tagdata;
  char *next;

  next = strchr(buf,'\n');
  *next = '\0';
  int nwords = count_words(buf);
  *next = '\n';

  if (nwords != avec_bonus->size_data_bonus)
    error->all(FLERR,"Incorrect bonus data format in data file");

  char **values = new char*[nwords];

  // loop over lines of bonus atom data
  // tokenize the line into values
  // if I own atom tag, unpack its values

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');

    values[0] = strtok(buf," \t\n\r\f");
    for (j = 1; j < nwords; j++)
      values[j] = strtok(NULL," \t\n\r\f");

    tagdata = ATOTAGINT(values[0]) + id_offset;
    if (tagdata <= 0 || tagdata > map_tag_max)
      error->one(FLERR,"Invalid atom ID in Bonus section of data file");

    // ok to call child's data_atom_bonus() method thru parent avec_bonus,
    // since data_bonus() was called with child ptr, and method is virtual

    if ((m = map(tagdata)) >= 0) avec_bonus->data_atom_bonus(m,&values[1]);

    buf = next + 1;
  }

  delete [] values;
}

/* ----------------------------------------------------------------------
   unpack N bodies from Bodies section of data file
   each body spans multiple lines
   check that atom IDs are > 0 and <= map_tag_max
   call style-specific routine to parse line
------------------------------------------------------------------------- */

void Atom::data_bodies(int n, char *buf, AtomVec *avec_body, tagint id_offset)
{
  int j,m,nvalues,tagdata,ninteger,ndouble;

  int maxint = 0;
  int maxdouble = 0;
  int *ivalues = NULL;
  double *dvalues = NULL;

  // loop over lines of body data
  // if I own atom tag, tokenize lines into ivalues/dvalues, call data_body()
  // else skip values

  for (int i = 0; i < n; i++) {
    if (i == 0) tagdata = ATOTAGINT(strtok(buf," \t\n\r\f")) + id_offset;
    else tagdata = ATOTAGINT(strtok(NULL," \t\n\r\f")) + id_offset;

    if (tagdata <= 0 || tagdata > map_tag_max)
      error->one(FLERR,"Invalid atom ID in Bodies section of data file");

    ninteger = force->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
    ndouble = force->inumeric(FLERR,strtok(NULL," \t\n\r\f"));

    if ((m = map(tagdata)) >= 0) {
      if (ninteger > maxint) {
        delete [] ivalues;
        maxint = ninteger;
        ivalues = new int[maxint];
      }
      if (ndouble > maxdouble) {
        delete [] dvalues;
        maxdouble = ndouble;
        dvalues = new double[maxdouble];
      }

      for (j = 0; j < ninteger; j++)
        ivalues[j] = force->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
      for (j = 0; j < ndouble; j++)
        dvalues[j] = force->numeric(FLERR,strtok(NULL," \t\n\r\f"));

      avec_body->data_body(m,ninteger,ndouble,ivalues,dvalues);

    } else {
      nvalues = ninteger + ndouble;    // number of values to skip
      for (j = 0; j < nvalues; j++)
        strtok(NULL," \t\n\r\f");
    }
  }

  delete [] ivalues;
  delete [] dvalues;
}

/* ----------------------------------------------------------------------
   init per-atom fix/compute/variable values for newly created atoms
   called from create_atoms, read_data, read_dump,
     lib::lammps_create_atoms()
   fixes, computes, variables may or may not exist when called
------------------------------------------------------------------------- */

void Atom::data_fix_compute_variable(int nprev, int nnew)
{
  for (int m = 0; m < modify->nfix; m++) {
    Fix *fix = modify->fix[m];
    if (fix->create_attribute)
      for (int i = nprev; i < nnew; i++)
        fix->set_arrays(i);
  }

  for (int m = 0; m < modify->ncompute; m++) {
    Compute *compute = modify->compute[m];
    if (compute->create_attribute)
      for (int i = nprev; i < nnew; i++)
        compute->set_arrays(i);
  }

  for (int i = nprev; i < nnew; i++)
    input->variable->set_arrays(i);
}

/* ----------------------------------------------------------------------
   allocate arrays of length ntypes
   only done after ntypes is set
------------------------------------------------------------------------- */

void Atom::allocate_type_arrays()
{
  if (avec->mass_type) {
    mass = new double[ntypes+1];
    mass_setflag = new int[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) mass_setflag[itype] = 0;
  }
}

/* ----------------------------------------------------------------------
   set a mass and flag it as set
   called from reading of data file
   type_offset may be used when reading multiple data files
------------------------------------------------------------------------- */

void Atom::set_mass(const char *file, int line, const char *str, int type_offset)
{
  if (mass == NULL) error->all(file,line,"Cannot set mass for this atom style");

  int itype;
  double mass_one;
  int n = sscanf(str,"%d %lg",&itype,&mass_one);
  if (n != 2) error->all(file,line,"Invalid mass line in data file");
  itype += type_offset;

  if (itype < 1 || itype > ntypes)
    error->all(file,line,"Invalid type for mass set");

  mass[itype] = mass_one;
  mass_setflag[itype] = 1;

  if (mass[itype] <= 0.0) error->all(file,line,"Invalid mass value");
}

/* ----------------------------------------------------------------------
   set a mass and flag it as set
   called from EAM pair routine
------------------------------------------------------------------------- */

void Atom::set_mass(const char *file, int line, int itype, double value)
{
  if (mass == NULL) error->all(file,line,"Cannot set mass for this atom style");
  if (itype < 1 || itype > ntypes)
    error->all(file,line,"Invalid type for mass set");

  mass[itype] = value;
  mass_setflag[itype] = 1;

  if (mass[itype] <= 0.0) error->all(file,line,"Invalid mass value");
}

/* ----------------------------------------------------------------------
   set one or more masses and flag them as set
   called from reading of input script
------------------------------------------------------------------------- */

void Atom::set_mass(const char *file, int line, int /*narg*/, char **arg)
{
  if (mass == NULL) error->all(file,line,"Cannot set mass for this atom style");

  int lo,hi;
  force->bounds(file,line,arg[0],ntypes,lo,hi);
  if (lo < 1 || hi > ntypes) error->all(file,line,"Invalid type for mass set");

  for (int itype = lo; itype <= hi; itype++) {
    mass[itype] = atof(arg[1]);
    mass_setflag[itype] = 1;

    if (mass[itype] <= 0.0) error->all(file,line,"Invalid mass value");
  }
}

/* ----------------------------------------------------------------------
   set all masses
   called from reading of restart file, also from ServerMD
------------------------------------------------------------------------- */

void Atom::set_mass(double *values)
{
  for (int itype = 1; itype <= ntypes; itype++) {
    mass[itype] = values[itype];
    mass_setflag[itype] = 1;
  }
}

/* ----------------------------------------------------------------------
   check that all per-atom-type masses have been set
------------------------------------------------------------------------- */

void Atom::check_mass(const char *file, int line)
{
  if (mass == NULL) return;
  if (rmass_flag) return;
  for (int itype = 1; itype <= ntypes; itype++)
    if (mass_setflag[itype] == 0)
      error->all(file,line,"Not all per-type masses are set");
}

/* ----------------------------------------------------------------------
   check that radii of all particles of itype are the same
   return 1 if true, else return 0
   also return the radius value for that type
------------------------------------------------------------------------- */

int Atom::radius_consistency(int itype, double &rad)
{
  double value = -1.0;
  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] != itype) continue;
    if (value < 0.0) value = radius[i];
    else if (value != radius[i]) flag = 1;
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) return 0;

  MPI_Allreduce(&value,&rad,1,MPI_DOUBLE,MPI_MAX,world);
  return 1;
}

/* ----------------------------------------------------------------------
   check that shape of all particles of itype are the same
   return 1 if true, else return 0
   also return the 3 shape params for itype
------------------------------------------------------------------------- */

int Atom::shape_consistency(int itype,
                            double &shapex, double &shapey, double &shapez)
{
  double zero[3] = {0.0, 0.0, 0.0};
  double one[3] = {-1.0, -1.0, -1.0};
  double *shape;

  AtomVecEllipsoid *avec_ellipsoid =
    (AtomVecEllipsoid *) style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec_ellipsoid->bonus;

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] != itype) continue;
    if (ellipsoid[i] < 0) shape = zero;
    else shape = bonus[ellipsoid[i]].shape;

    if (one[0] < 0.0) {
      one[0] = shape[0];
      one[1] = shape[1];
      one[2] = shape[2];
    } else if (one[0] != shape[0] || one[1] != shape[1] || one[2] != shape[2])
      flag = 1;
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) return 0;

  double oneall[3];
  MPI_Allreduce(one,oneall,3,MPI_DOUBLE,MPI_MAX,world);
  shapex = oneall[0];
  shapey = oneall[1];
  shapez = oneall[2];
  return 1;
}

/* ----------------------------------------------------------------------
   add a new molecule template = set of molecules
------------------------------------------------------------------------- */

void Atom::add_molecule(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal molecule command");

  if (find_molecule(arg[0]) >= 0)
    error->all(FLERR,"Reuse of molecule template ID");

  // 1st molecule in set stores nset = # of mols, others store nset = 0
  // ifile = count of molecules in set
  // index = argument index where next molecule starts, updated by constructor

  int ifile = 1;
  int index = 1;
  while (1) {
    molecules = (Molecule **)
      memory->srealloc(molecules,(nmolecule+1)*sizeof(Molecule *),
                       "atom::molecules");
    molecules[nmolecule] = new Molecule(lmp,narg,arg,index);
    molecules[nmolecule]->nset = 0;
    molecules[nmolecule-ifile+1]->nset++;
    nmolecule++;
    if (molecules[nmolecule-1]->last) break;
    ifile++;
  }
}

/* ----------------------------------------------------------------------
   find first molecule in set with template ID
   return -1 if does not exist
------------------------------------------------------------------------- */

int Atom::find_molecule(char *id)
{
  if(id == NULL) return -1;
  int imol;
  for (imol = 0; imol < nmolecule; imol++)
    if (strcmp(id,molecules[imol]->id) == 0) return imol;
  return -1;
}

/* ----------------------------------------------------------------------
   add info to current atom ilocal from molecule template onemol and its iatom
   offset = atom ID preceding IDs of atoms in this molecule
   called by fixes and commands that add molecules
------------------------------------------------------------------------- */

void Atom::add_molecule_atom(Molecule *onemol, int iatom,
                             int ilocal, tagint offset)
{
  if (onemol->qflag && q_flag) q[ilocal] = onemol->q[iatom];
  if (onemol->radiusflag && radius_flag) radius[ilocal] = onemol->radius[iatom];
  if (onemol->rmassflag && rmass_flag) rmass[ilocal] = onemol->rmass[iatom];
  else if (rmass_flag)
    rmass[ilocal] = 4.0*MY_PI/3.0 *
      radius[ilocal]*radius[ilocal]*radius[ilocal];
  if (onemol->bodyflag) {
    body[ilocal] = 0;     // as if a body read from data file
    onemol->avec_body->data_body(ilocal,onemol->nibody,onemol->ndbody,
                                 onemol->ibodyparams,onemol->dbodyparams);
    onemol->avec_body->set_quat(ilocal,onemol->quat_external);
  }

  if (molecular != 1) return;

  // add bond topology info
  // for molecular atom styles, but not atom style template

  if (avec->bonds_allow) {
    num_bond[ilocal] = onemol->num_bond[iatom];
    for (int i = 0; i < num_bond[ilocal]; i++) {
      bond_type[ilocal][i] = onemol->bond_type[iatom][i];
      bond_atom[ilocal][i] = onemol->bond_atom[iatom][i] + offset;
    }
  }

  if (avec->angles_allow) {
    num_angle[ilocal] = onemol->num_angle[iatom];
    for (int i = 0; i < num_angle[ilocal]; i++) {
      angle_type[ilocal][i] = onemol->angle_type[iatom][i];
      angle_atom1[ilocal][i] = onemol->angle_atom1[iatom][i] + offset;
      angle_atom2[ilocal][i] = onemol->angle_atom2[iatom][i] + offset;
      angle_atom3[ilocal][i] = onemol->angle_atom3[iatom][i] + offset;
    }
  }

  if (avec->dihedrals_allow) {
    num_dihedral[ilocal] = onemol->num_dihedral[iatom];
    for (int i = 0; i < num_dihedral[ilocal]; i++) {
      dihedral_type[ilocal][i] = onemol->dihedral_type[iatom][i];
      dihedral_atom1[ilocal][i] = onemol->dihedral_atom1[iatom][i] + offset;
      dihedral_atom2[ilocal][i] = onemol->dihedral_atom2[iatom][i] + offset;
      dihedral_atom3[ilocal][i] = onemol->dihedral_atom3[iatom][i] + offset;
      dihedral_atom4[ilocal][i] = onemol->dihedral_atom4[iatom][i] + offset;
    }
  }

  if (avec->impropers_allow) {
    num_improper[ilocal] = onemol->num_improper[iatom];
    for (int i = 0; i < num_improper[ilocal]; i++) {
      improper_type[ilocal][i] = onemol->improper_type[iatom][i];
      improper_atom1[ilocal][i] = onemol->improper_atom1[iatom][i] + offset;
      improper_atom2[ilocal][i] = onemol->improper_atom2[iatom][i] + offset;
      improper_atom3[ilocal][i] = onemol->improper_atom3[iatom][i] + offset;
      improper_atom4[ilocal][i] = onemol->improper_atom4[iatom][i] + offset;
    }
  }

  if (onemol->specialflag) {
    nspecial[ilocal][0] = onemol->nspecial[iatom][0];
    nspecial[ilocal][1] = onemol->nspecial[iatom][1];
    int n = nspecial[ilocal][2] = onemol->nspecial[iatom][2];
    for (int i = 0; i < n; i++)
      special[ilocal][i] = onemol->special[iatom][i] + offset;
  }
}

/* ----------------------------------------------------------------------
   reorder owned atoms so those in firstgroup appear first
   called by comm->exchange() if atom_modify first group is set
   only owned atoms exist at this point, no ghost atoms
------------------------------------------------------------------------- */

void Atom::first_reorder()
{
  // insure there is one extra atom location at end of arrays for swaps

  if (nlocal == nmax) avec->grow(0);

  // loop over owned atoms
  // nfirst = index of first atom not in firstgroup
  // when find firstgroup atom out of place, swap it with atom nfirst

  int bitmask = group->bitmask[firstgroup];
  nfirst = 0;
  while (nfirst < nlocal && mask[nfirst] & bitmask) nfirst++;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & bitmask && i > nfirst) {
      avec->copy(i,nlocal,0);
      avec->copy(nfirst,i,0);
      avec->copy(nlocal,nfirst,0);
      while (nfirst < nlocal && mask[nfirst] & bitmask) nfirst++;
    }
  }
}

/* ----------------------------------------------------------------------
   perform spatial sort of atoms within my sub-domain
   always called between comm->exchange() and comm->borders()
   don't have to worry about clearing/setting atom->map since done in comm
------------------------------------------------------------------------- */

void Atom::sort()
{
  int i,m,n,ix,iy,iz,ibin,empty;

  // set next timestep for sorting to take place

  nextsort = (update->ntimestep/sortfreq)*sortfreq + sortfreq;

  // re-setup sort bins if needed

  if (domain->box_change) setup_sort_bins();
  if (nbins == 1) return;

  // reallocate per-atom vectors if needed

  if (nlocal > maxnext) {
    memory->destroy(next);
    memory->destroy(permute);
    maxnext = atom->nmax;
    memory->create(next,maxnext,"atom:next");
    memory->create(permute,maxnext,"atom:permute");
  }

  // insure there is one extra atom location at end of arrays for swaps

  if (nlocal == nmax) avec->grow(0);

  // bin atoms in reverse order so linked list will be in forward order

  for (i = 0; i < nbins; i++) binhead[i] = -1;

  for (i = nlocal-1; i >= 0; i--) {
    ix = static_cast<int> ((x[i][0]-bboxlo[0])*bininvx);
    iy = static_cast<int> ((x[i][1]-bboxlo[1])*bininvy);
    iz = static_cast<int> ((x[i][2]-bboxlo[2])*bininvz);
    ix = MAX(ix,0);
    iy = MAX(iy,0);
    iz = MAX(iz,0);
    ix = MIN(ix,nbinx-1);
    iy = MIN(iy,nbiny-1);
    iz = MIN(iz,nbinz-1);
    ibin = iz*nbiny*nbinx + iy*nbinx + ix;
    next[i] = binhead[ibin];
    binhead[ibin] = i;
  }

  // permute = desired permutation of atoms
  // permute[I] = J means Ith new atom will be Jth old atom

  n = 0;
  for (m = 0; m < nbins; m++) {
    i = binhead[m];
    while (i >= 0) {
      permute[n++] = i;
      i = next[i];
    }
  }

  // current = current permutation, just reuse next vector
  // current[I] = J means Ith current atom is Jth old atom

  int *current = next;
  for (i = 0; i < nlocal; i++) current[i] = i;

  // reorder local atom list, when done, current = permute
  // perform "in place" using copy() to extra atom location at end of list
  // inner while loop processes one cycle of the permutation
  // copy before inner-loop moves an atom to end of atom list
  // copy after inner-loop moves atom at end of list back into list
  // empty = location in atom list that is currently empty

  for (i = 0; i < nlocal; i++) {
    if (current[i] == permute[i]) continue;
    avec->copy(i,nlocal,0);
    empty = i;
    while (permute[empty] != i) {
      avec->copy(permute[empty],empty,0);
      empty = current[empty] = permute[empty];
    }
    avec->copy(nlocal,empty,0);
    current[empty] = permute[empty];
  }

  // sanity check that current = permute

  //int flag = 0;
  //for (i = 0; i < nlocal; i++)
  //  if (current[i] != permute[i]) flag = 1;
  //int flagall;
  //MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  //if (flagall) error->all(FLERR,"Atom sort did not operate correctly");
}

/* ----------------------------------------------------------------------
   setup bins for spatial sorting of atoms
------------------------------------------------------------------------- */

void Atom::setup_sort_bins()
{
  // binsize:
  // user setting if explicitly set
  // default = 1/2 of neighbor cutoff
  // check if neighbor cutoff = 0.0
  // and in that case, disable sorting

  double binsize = 0.0;
  if (userbinsize > 0.0) binsize = userbinsize;
  else if (neighbor->cutneighmax > 0.0) binsize = 0.5 * neighbor->cutneighmax;

  if ((binsize == 0.0) && (sortfreq > 0)) {
    sortfreq = 0;
    if (comm->me == 0)
          error->warning(FLERR,"No pairwise cutoff or binsize set. "
                         "Atom sorting therefore disabled.");
    return;
  }

  double bininv = 1.0/binsize;

  // nbin xyz = local bins
  // bbox lo/hi = bounding box of my sub-domain

  if (domain->triclinic)
    domain->bbox(domain->sublo_lamda,domain->subhi_lamda,bboxlo,bboxhi);
  else {
    bboxlo[0] = domain->sublo[0];
    bboxlo[1] = domain->sublo[1];
    bboxlo[2] = domain->sublo[2];
    bboxhi[0] = domain->subhi[0];
    bboxhi[1] = domain->subhi[1];
    bboxhi[2] = domain->subhi[2];
  }

  nbinx = static_cast<int> ((bboxhi[0]-bboxlo[0]) * bininv);
  nbiny = static_cast<int> ((bboxhi[1]-bboxlo[1]) * bininv);
  nbinz = static_cast<int> ((bboxhi[2]-bboxlo[2]) * bininv);
  if (domain->dimension == 2) nbinz = 1;

  if (nbinx == 0) nbinx = 1;
  if (nbiny == 0) nbiny = 1;
  if (nbinz == 0) nbinz = 1;

  bininvx = nbinx / (bboxhi[0]-bboxlo[0]);
  bininvy = nbiny / (bboxhi[1]-bboxlo[1]);
  bininvz = nbinz / (bboxhi[2]-bboxlo[2]);

  #ifdef LMP_USER_INTEL
  int intel_neigh = 0;
  if (neighbor->nrequest) {
    if (neighbor->requests[0]->intel) intel_neigh = 1;
  } else if (neighbor->old_nrequest)
    if (neighbor->old_requests[0]->intel) intel_neigh = 1;
  if (intel_neigh && userbinsize == 0.0) {
    if (neighbor->binsizeflag) bininv = 1.0/neighbor->binsize_user;

    double nx_low = neighbor->bboxlo[0];
    double ny_low = neighbor->bboxlo[1];
    double nz_low = neighbor->bboxlo[2];
    double nxbbox = neighbor->bboxhi[0] - nx_low;
    double nybbox = neighbor->bboxhi[1] - ny_low;
    double nzbbox = neighbor->bboxhi[2] - nz_low;
    int nnbinx = static_cast<int> (nxbbox * bininv);
    int nnbiny = static_cast<int> (nybbox * bininv);
    int nnbinz = static_cast<int> (nzbbox * bininv);
    if (domain->dimension == 2) nnbinz = 1;

    if (nnbinx == 0) nnbinx = 1;
    if (nnbiny == 0) nnbiny = 1;
    if (nnbinz == 0) nnbinz = 1;

    double binsizex = nxbbox/nnbinx;
    double binsizey = nybbox/nnbiny;
    double binsizez = nzbbox/nnbinz;

    bininvx = 1.0 / binsizex;
    bininvy = 1.0 / binsizey;
    bininvz = 1.0 / binsizez;

    int lxo = (bboxlo[0] - nx_low) * bininvx;
    int lyo = (bboxlo[1] - ny_low) * bininvy;
    int lzo = (bboxlo[2] - nz_low) * bininvz;
    bboxlo[0] = nx_low + static_cast<double>(lxo) / bininvx;
    bboxlo[1] = ny_low + static_cast<double>(lyo) / bininvy;
    bboxlo[2] = nz_low + static_cast<double>(lzo) / bininvz;
    nbinx = static_cast<int>((bboxhi[0] - bboxlo[0]) * bininvx) + 1;
    nbiny = static_cast<int>((bboxhi[1] - bboxlo[1]) * bininvy) + 1;
    nbinz = static_cast<int>((bboxhi[2] - bboxlo[2]) * bininvz) + 1;
    bboxhi[0] = bboxlo[0] + static_cast<double>(nbinx) / bininvx;
    bboxhi[1] = bboxlo[1] + static_cast<double>(nbiny) / bininvy;
    bboxhi[2] = bboxlo[2] + static_cast<double>(nbinz) / bininvz;
  }
  #endif

  if (1.0*nbinx*nbiny*nbinz > INT_MAX)
    error->one(FLERR,"Too many atom sorting bins");

  nbins = nbinx*nbiny*nbinz;

  // reallocate per-bin memory if needed

  if (nbins > maxbin) {
    memory->destroy(binhead);
    maxbin = nbins;
    memory->create(binhead,maxbin,"atom:binhead");
  }
}

/* ----------------------------------------------------------------------
   register a callback to a fix so it can manage atom-based arrays
   happens when fix is created
   flag = 0 for grow, 1 for restart, 2 for border comm
------------------------------------------------------------------------- */

void Atom::add_callback(int flag)
{
  int ifix;

  // find the fix
  // if find NULL ptr:
  //   it's this one, since it is being replaced and has just been deleted
  //   at this point in re-creation
  // if don't find NULL ptr:
  //   i is set to nfix = new one currently being added at end of list

  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (modify->fix[ifix] == NULL) break;

  // add callback to lists, reallocating if necessary

  if (flag == 0) {
    if (nextra_grow == nextra_grow_max) {
      nextra_grow_max += DELTA;
      memory->grow(extra_grow,nextra_grow_max,"atom:extra_grow");
    }
    extra_grow[nextra_grow] = ifix;
    nextra_grow++;
  } else if (flag == 1) {
    if (nextra_restart == nextra_restart_max) {
      nextra_restart_max += DELTA;
      memory->grow(extra_restart,nextra_restart_max,"atom:extra_restart");
    }
    extra_restart[nextra_restart] = ifix;
    nextra_restart++;
  } else if (flag == 2) {
    if (nextra_border == nextra_border_max) {
      nextra_border_max += DELTA;
      memory->grow(extra_border,nextra_border_max,"atom:extra_border");
    }
    extra_border[nextra_border] = ifix;
    nextra_border++;
  }
}

/* ----------------------------------------------------------------------
   unregister a callback to a fix
   happens when fix is deleted, called by its destructor
   flag = 0 for grow, 1 for restart
------------------------------------------------------------------------- */

void Atom::delete_callback(const char *id, int flag)
{
  if (id == NULL) return;

  int ifix = modify->find_fix(id);

  // compact the list of callbacks

  if (flag == 0) {
    int match;
    for (match = 0; match < nextra_grow; match++)
      if (extra_grow[match] == ifix) break;
    if ((nextra_grow == 0) || (match == nextra_grow))
      error->all(FLERR,"Trying to delete non-existent Atom::grow() callback");
    for (int i = match; i < nextra_grow-1; i++)
      extra_grow[i] = extra_grow[i+1];
    nextra_grow--;

  } else if (flag == 1) {
    int match;
    for (match = 0; match < nextra_restart; match++)
      if (extra_restart[match] == ifix) break;
    if ((nextra_restart == 0) || (match == nextra_restart))
      error->all(FLERR,"Trying to delete non-existent Atom::restart() callback");
    for (int i = match; i < nextra_restart-1; i++)
      extra_restart[i] = extra_restart[i+1];
    nextra_restart--;

  } else if (flag == 2) {
    int match;
    for (match = 0; match < nextra_border; match++)
      if (extra_border[match] == ifix) break;
    if ((nextra_border == 0) || (match == nextra_border))
      error->all(FLERR,"Trying to delete non-existent Atom::border() callback");
    for (int i = match; i < nextra_border-1; i++)
      extra_border[i] = extra_border[i+1];
    nextra_border--;
  }
}

/* ----------------------------------------------------------------------
   decrement ptrs in callback lists to fixes beyond the deleted ifix
   happens after fix is deleted
------------------------------------------------------------------------- */

void Atom::update_callback(int ifix)
{
  for (int i = 0; i < nextra_grow; i++)
    if (extra_grow[i] > ifix) extra_grow[i]--;
  for (int i = 0; i < nextra_restart; i++)
    if (extra_restart[i] > ifix) extra_restart[i]--;
  for (int i = 0; i < nextra_border; i++)
    if (extra_border[i] > ifix) extra_border[i]--;
}

/* ----------------------------------------------------------------------
   find custom per-atom vector with name
   return index if found, and flag = 0/1 for int/double
   return -1 if not found
------------------------------------------------------------------------- */

int Atom::find_custom(const char *name, int &flag)
{
  if(name == NULL) return -1;

  for (int i = 0; i < nivector; i++)
    if (iname[i] && strcmp(iname[i],name) == 0) {
      flag = 0;
      return i;
    }

  for (int i = 0; i < ndvector; i++)
    if (dname[i] && strcmp(dname[i],name) == 0) {
      flag = 1;
      return i;
    }

  return -1;
}

/* ----------------------------------------------------------------------
   add a custom variable with name of type flag = 0/1 for int/double
   assumes name does not already exist
   return index in ivector or dvector of its location
------------------------------------------------------------------------- */

int Atom::add_custom(const char *name, int flag)
{
  int index;

  if (flag == 0) {
    index = nivector;
    nivector++;
    iname = (char **) memory->srealloc(iname,nivector*sizeof(char *),
                                       "atom:iname");
    int n = strlen(name) + 1;
    iname[index] = new char[n];
    strcpy(iname[index],name);
    ivector = (int **) memory->srealloc(ivector,nivector*sizeof(int *),
                                        "atom:ivector");
    memory->create(ivector[index],nmax,"atom:ivector");
  } else {
    index = ndvector;
    ndvector++;
    dname = (char **) memory->srealloc(dname,ndvector*sizeof(char *),
                                       "atom:dname");
    int n = strlen(name) + 1;
    dname[index] = new char[n];
    strcpy(dname[index],name);
    dvector = (double **) memory->srealloc(dvector,ndvector*sizeof(double *),
                                           "atom:dvector");
    memory->create(dvector[index],nmax,"atom:dvector");
  }

  return index;
}

/* ----------------------------------------------------------------------
   remove a custom variable of type flag = 0/1 for int/double at index
   free memory for vector and name and set ptrs to NULL
   ivector/dvector and iname/dname lists never shrink
------------------------------------------------------------------------- */

void Atom::remove_custom(int flag, int index)
{
  if (flag == 0) {
    memory->destroy(ivector[index]);
    ivector[index] = NULL;
    delete [] iname[index];
    iname[index] = NULL;
  } else {
    memory->destroy(dvector[index]);
    dvector[index] = NULL;
    delete [] dname[index];
    dname[index] = NULL;
  }
}

/* ----------------------------------------------------------------------
   return a pointer to a named internal variable
   if don't recognize name, return NULL
   customize by adding names
------------------------------------------------------------------------- */

void *Atom::extract(char *name)
{
  if (strcmp(name,"mass") == 0) return (void *) mass;

  if (strcmp(name,"id") == 0) return (void *) tag;
  if (strcmp(name,"type") == 0) return (void *) type;
  if (strcmp(name,"mask") == 0) return (void *) mask;
  if (strcmp(name,"image") == 0) return (void *) image;
  if (strcmp(name,"x") == 0) return (void *) x;
  if (strcmp(name,"v") == 0) return (void *) v;
  if (strcmp(name,"f") == 0) return (void *) f;
  if (strcmp(name,"molecule") == 0) return (void *) molecule;
  if (strcmp(name,"q") == 0) return (void *) q;
  if (strcmp(name,"mu") == 0) return (void *) mu;
  if (strcmp(name,"omega") == 0) return (void *) omega;
  if (strcmp(name,"angmom") == 0) return (void *) angmom;
  if (strcmp(name,"torque") == 0) return (void *) torque;
  if (strcmp(name,"radius") == 0) return (void *) radius;
  if (strcmp(name,"rmass") == 0) return (void *) rmass;
  if (strcmp(name,"ellipsoid") == 0) return (void *) ellipsoid;
  if (strcmp(name,"line") == 0) return (void *) line;
  if (strcmp(name,"tri") == 0) return (void *) tri;

  if (strcmp(name,"vfrac") == 0) return (void *) vfrac;
  if (strcmp(name,"s0") == 0) return (void *) s0;
  if (strcmp(name,"x0") == 0) return (void *) x0;

  if (strcmp(name,"spin") == 0) return (void *) spin;
  if (strcmp(name,"eradius") == 0) return (void *) eradius;
  if (strcmp(name,"ervel") == 0) return (void *) ervel;
  if (strcmp(name,"erforce") == 0) return (void *) erforce;
  if (strcmp(name,"ervelforce") == 0) return (void *) ervelforce;
  if (strcmp(name,"cs") == 0) return (void *) cs;
  if (strcmp(name,"csforce") == 0) return (void *) csforce;
  if (strcmp(name,"vforce") == 0) return (void *) vforce;
  if (strcmp(name,"etag") == 0) return (void *) etag;

  if (strcmp(name,"rho") == 0) return (void *) rho;
  if (strcmp(name,"drho") == 0) return (void *) drho;
  if (strcmp(name,"e") == 0) return (void *) e;
  if (strcmp(name,"de") == 0) return (void *) de;
  if (strcmp(name,"cv") == 0) return (void *) cv;
  if (strcmp(name,"vest") == 0) return (void *) vest;

  if (strcmp(name, "contact_radius") == 0) return (void *) contact_radius;
  if (strcmp(name, "smd_data_9") == 0) return (void *) smd_data_9;
  if (strcmp(name, "smd_stress") == 0) return (void *) smd_stress;
  if (strcmp(name, "eff_plastic_strain") == 0)
    return (void *) eff_plastic_strain;
  if (strcmp(name, "eff_plastic_strain_rate") == 0)
    return (void *) eff_plastic_strain_rate;
  if (strcmp(name, "damage") == 0) return (void *) damage;

  if (strcmp(name,"dpdTheta") == 0) return (void *) dpdTheta;
  if (strcmp(name,"edpd_temp") == 0) return (void *) edpd_temp;

  return NULL;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   call to avec tallies per-atom vectors
   add in global to local mapping storage
------------------------------------------------------------------------- */

bigint Atom::memory_usage()
{
  bigint bytes = avec->memory_usage();

  bytes += max_same*sizeof(int);
  if (map_style == 1)
    bytes += memory->usage(map_array,map_maxarray);
  else if (map_style == 2) {
    bytes += map_nbucket*sizeof(int);
    bytes += map_nhash*sizeof(HashElem);
  }
  if (maxnext) {
    bytes += memory->usage(next,maxnext);
    bytes += memory->usage(permute,maxnext);
  }

  return bytes;
}

