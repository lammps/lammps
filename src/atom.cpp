// clang-format off
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

#include "atom.h"
#include "atom_vec.h"
#include "style_atom.h"  // IWYU pragma: keep

#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "label_map.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "neighbor.h"
#include "tokenizer.h"
#include "update.h"
#include "variable.h"

#include "library.h"

#include <algorithm>
#include <cstring>

#ifdef LMP_GPU
#include "fix_gpu.h"
#include <cmath>
#endif

using namespace LAMMPS_NS;
using namespace MathConst;

#define DELTA 1
#define EPSILON 1.0e-6
#define MAXLINE 256

/* ----------------------------------------------------------------------
   one instance per AtomVec style in style_atom.h
------------------------------------------------------------------------- */

template <typename T> static AtomVec *avec_creator(LAMMPS *_lmp)
{
  return new T(_lmp);
}

/* ---------------------------------------------------------------------- */

/** \class LAMMPS_NS::Atom
 *  \brief Class to provide access to atom data

\verbatim embed:rst
The Atom class provides access to atom style related global settings and
per-atom data that is stored with atoms and migrates with them from
sub-domain to sub-domain as atoms move around.  This includes topology
data, which is stored with either one specific atom or all atoms involved
depending on the settings of the :doc:`newton command <newton>`.

The actual per-atom data is allocated and managed by one of the various
classes derived from the AtomVec class as determined by
the :doc:`atom_style command <atom_style>`.  The pointers in the Atom class
are updated by the AtomVec class as needed.
\endverbatim
 */

/** Atom class constructor
 *
 * This resets and initialized all kinds of settings,
 * parameters, and pointer variables for per-atom arrays.
 * This also initializes the factory for creating
 * instances of classes derived from the AtomVec base
 * class, which correspond to the selected atom style.
 *
 * \param  _lmp  pointer to the base LAMMPS class */

Atom::Atom(LAMMPS *_lmp) : Pointers(_lmp), atom_style(nullptr), avec(nullptr), avec_map(nullptr)
{
  natoms = 0;
  nlocal = nghost = nmax = 0;
  ntypes = 0;
  nellipsoids = nlines = ntris = nbodies = 0;
  nbondtypes = nangletypes = ndihedraltypes = nimpropertypes = 0;
  nbonds = nangles = ndihedrals = nimpropers = 0;

  firstgroupname = nullptr;
  sortfreq = 1000;
  nextsort = 0;
  userbinsize = 0.0;
  maxbin = maxnext = 0;
  binhead = nullptr;
  next = permute = nullptr;

  // --------------------------------------------------------------------
  // 1st customization section: customize by adding new per-atom variables

  tag = nullptr;
  type = mask = nullptr;
  image = nullptr;
  x = v = f = nullptr;

  // charged and dipolar particles

  q = nullptr;
  mu = nullptr;

  // finite-size particles

  omega = angmom = torque = nullptr;
  radius = rmass = nullptr;
  ellipsoid = line = tri = body = nullptr;
  quat = nullptr;
  temperature = nullptr;
  heatflow = nullptr;

  // molecular systems

  molecule = nullptr;
  molindex = molatom = nullptr;

  bond_per_atom =  extra_bond_per_atom = 0;
  num_bond = nullptr;
  bond_type = nullptr;
  bond_atom = nullptr;

  angle_per_atom = extra_angle_per_atom = 0;
  num_angle = nullptr;
  angle_type = nullptr;
  angle_atom1 = angle_atom2 = angle_atom3 = nullptr;

  dihedral_per_atom = extra_dihedral_per_atom = 0;
  num_dihedral = nullptr;
  dihedral_type = nullptr;
  dihedral_atom1 = dihedral_atom2 = dihedral_atom3 = dihedral_atom4 = nullptr;

  improper_per_atom = extra_improper_per_atom = 0;
  num_improper = nullptr;
  improper_type = nullptr;
  improper_atom1 = improper_atom2 = improper_atom3 = improper_atom4 = nullptr;

  maxspecial = 1;
  nspecial = nullptr;
  special = nullptr;

  // PERI package

  vfrac = s0 = nullptr;
  x0 = nullptr;

  // SPIN package

  sp = fm = fm_long = nullptr;

  // EFF and AWPMD packages

  spin = nullptr;
  eradius = ervel = erforce = nullptr;
  ervelforce = nullptr;
  cs = csforce = vforce = nullptr;
  etag = nullptr;

  // CG-DNA package

  id5p = nullptr;

  // DPD-REACT package

  uCond = uMech = uChem = uCG = uCGnew = nullptr;
  duChem = dpdTheta = nullptr;

  // MESO package

  cc = cc_flux = nullptr;
  edpd_temp = edpd_flux = edpd_cv = vest_temp = nullptr;

  // MACHDYN package

  contact_radius = nullptr;
  smd_data_9 = nullptr;
  smd_stress = nullptr;
  eff_plastic_strain = nullptr;
  eff_plastic_strain_rate = nullptr;
  damage = nullptr;

  // SPH package

  rho = drho = esph = desph = cv = nullptr;
  vest = nullptr;

  // AMOEBA package

  maxspecial15 = 1;
  nspecial15 = nullptr;
  special15 = nullptr;

  // DIELECTRIC package

  area = ed = em = epsilon = curvature = q_scaled = nullptr;

  // end of customization section
  // --------------------------------------------------------------------

  // user-defined molecules

  nmolecule = 0;
  molecules = nullptr;

  // type labels

  lmap = nullptr;
  types_style = NUMERIC;

  // custom atom arrays

  nivector = ndvector = niarray = ndarray = 0;
  ivector = nullptr;
  dvector = nullptr;
  iarray = nullptr;
  darray = nullptr;
  icols = dcols = nullptr;
  ivname = dvname = ianame = daname = nullptr;

  // initialize atom style and array existence flags

  set_atomflag_defaults();

  // initialize peratom data structure

  peratom_create();

  // ntype-length arrays

  mass = nullptr;
  mass_setflag = nullptr;

  // callback lists & extra restart info

  nextra_grow = nextra_restart = nextra_border = 0;
  extra_grow = extra_restart = extra_border = nullptr;
  nextra_grow_max = nextra_restart_max = nextra_border_max = 0;
  nextra_store = 0;
  extra = nullptr;

  // default atom ID and mapping values

  tag_enable = 1;
  map_style = map_user = MAP_NONE;
  map_tag_max = -1;
  map_maxarray = map_nhash = map_nbucket = -1;

  max_same = 0;
  sametag = nullptr;
  map_array = nullptr;
  map_bucket = nullptr;
  map_hash = nullptr;

  unique_tags = nullptr;
  reset_image_flag[0] = reset_image_flag[1] = reset_image_flag[2] = false;

  avec_map = new AtomVecCreatorMap();

#define ATOM_CLASS
#define AtomStyle(key,Class) \
  (*avec_map)[#key] = &avec_creator<Class>;
#include "style_atom.h"  // IWYU pragma: keep
#undef AtomStyle
#undef ATOM_CLASS
}

/* ---------------------------------------------------------------------- */

Atom::~Atom()
{
  delete[] atom_style;
  delete avec;
  delete avec_map;

  delete[] firstgroupname;
  memory->destroy(binhead);
  memory->destroy(next);
  memory->destroy(permute);

  memory->destroy(tag);
  memory->destroy(type);
  memory->destroy(mask);
  memory->destroy(image);
  memory->destroy(x);
  memory->destroy(v);
  memory->destroy(f);

  // delete custom atom arrays

  for (int i = 0; i < nivector; i++) {
    delete[] ivname[i];
    memory->destroy(ivector[i]);
  }
  for (int i = 0; i < ndvector; i++) {
    delete[] dvname[i];
    if (dvector) // (needed for Kokkos)
      memory->destroy(dvector[i]);
  }
  for (int i = 0; i < niarray; i++) {
    delete[] ianame[i];
    memory->destroy(iarray[i]);
  }
  for (int i = 0; i < ndarray; i++) {
    delete[] daname[i];
    memory->destroy(darray[i]);
  }

  memory->sfree(ivname);
  memory->sfree(dvname);
  memory->sfree(ianame);
  memory->sfree(daname);
  memory->sfree(ivector);
  memory->sfree(dvector);
  memory->sfree(iarray);
  memory->sfree(darray);
  memory->sfree(icols);
  memory->sfree(dcols);

  // delete user-defined molecules

  for (int i = 0; i < nmolecule; i++) delete molecules[i];
  memory->sfree(molecules);

  // delete label map

  delete lmap;

  // delete per-type arrays

  delete[] mass;
  delete[] mass_setflag;

  // delete extra arrays

  memory->destroy(extra_grow);
  memory->destroy(extra_restart);
  memory->destroy(extra_border);
  memory->destroy(extra);

  // delete mapping data structures

  Atom::map_delete();

  delete unique_tags;
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
  if (old->firstgroupname)
    firstgroupname = utils::strdup(old->firstgroupname);
}

/* ----------------------------------------------------------------------
   one-time creation of peratom data structure
------------------------------------------------------------------------- */

void Atom::peratom_create()
{
  peratom.clear();

  // --------------------------------------------------------------------
  // 2nd customization section: add peratom variables here, order does not matter
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

  add_peratom("temperature",&temperature,DOUBLE,0);
  add_peratom("heatflow",&heatflow,DOUBLE,0);

  // BPM package

  add_peratom("quat",&quat,DOUBLE,4);

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

  // EFF package

  add_peratom("espin",&spin,INT,0);
  add_peratom("eradius",&eradius,DOUBLE,0);
  add_peratom("ervel",&ervel,DOUBLE,0);
  add_peratom("erforce",&erforce,DOUBLE,0,1);     // set per-thread flag

  // AWPMD package

  add_peratom("cs",&cs,DOUBLE,2);
  add_peratom("csforce",&csforce,DOUBLE,2);
  add_peratom("vforce",&vforce,DOUBLE,3);
  add_peratom("ervelforce",&ervelforce,DOUBLE,0);
  add_peratom("etag",&etag,INT,0);

  // CG-DNA package

  add_peratom("id5p",&id5p,tagintsize,0);

  // DPD-REACT package

  add_peratom("dpdTheta",&dpdTheta,DOUBLE,0);
  add_peratom("uCond",&uCond,DOUBLE,0);
  add_peratom("uMech",&uMech,DOUBLE,0);
  add_peratom("uChem",&uChem,DOUBLE,0);
  add_peratom("uCG",&uCG,DOUBLE,0);
  add_peratom("uCGnew",&uCGnew,DOUBLE,0);
  add_peratom("duChem",&duChem,DOUBLE,0);

  // MESO package

  add_peratom("edpd_cv",&edpd_cv,DOUBLE,0);
  add_peratom("edpd_temp",&edpd_temp,DOUBLE,0);
  add_peratom("vest_temp",&vest_temp,DOUBLE,0);
  add_peratom("edpd_flux",&edpd_flux,DOUBLE,0,1);     // set per-thread flag
  add_peratom("cc",&cc,DOUBLE,1);
  add_peratom("cc_flux",&cc_flux,DOUBLE,1,1);         // set per-thread flag

  // SPH package

  add_peratom("rho",&rho,DOUBLE,0);
  add_peratom("drho",&drho,DOUBLE,0,1);               // set per-thread flag
  add_peratom("esph",&esph,DOUBLE,0);
  add_peratom("desph",&desph,DOUBLE,0,1);             // set per-thread flag
  add_peratom("vest",&vest,DOUBLE,3);
  add_peratom("cv",&cv,DOUBLE,0);

  // MACHDYN package

  add_peratom("contact_radius",&contact_radius,DOUBLE,0);
  add_peratom("smd_data_9",&smd_data_9,DOUBLE,1);
  add_peratom("smd_stress",&smd_stress,DOUBLE,1);
  add_peratom("eff_plastic_strain",&eff_plastic_strain,DOUBLE,0);
  add_peratom("eff_plastic_strain_rate",&eff_plastic_strain_rate,DOUBLE,0);
  add_peratom("damage",&damage,DOUBLE,0);

  // AMOEBA package

  add_peratom("nspecial15",&nspecial15,INT,0);
  add_peratom_vary("special15",&special15,tagintsize,&maxspecial15,&nspecial15,0);

  // DIELECTRIC package

  add_peratom("area",&area,DOUBLE,0);
  add_peratom("ed",&ed,DOUBLE,0);
  add_peratom("em",&em,DOUBLE,0);
  add_peratom("epsilon",&epsilon,DOUBLE,0);
  add_peratom("curvature",&curvature,DOUBLE,0);
  add_peratom("q_scaled",&q_scaled,DOUBLE,0);

  // end of customization section
  // --------------------------------------------------------------------
}

/* ----------------------------------------------------------------------
   add info for a single per-atom vector/array to PerAtom data struct
   cols = 0: per-atom vector
   cols = N: static per-atom array with N columns
   use add_peratom_vary() when column count varies per atom
------------------------------------------------------------------------- */

void Atom::add_peratom(const std::string &name, void *address,
                       int datatype, int cols, int threadflag)
{
  PerAtom item = {name, address, nullptr, nullptr, datatype, cols, 0, threadflag};
  peratom.push_back(item);
}

/* ----------------------------------------------------------------------
   change the column count of an existing peratom array entry
   allows atom_style to specify column count as an argument
   see atom_style tdpd as an example
------------------------------------------------------------------------- */

void Atom::add_peratom_change_columns(const std::string &name, int cols)
{
  auto match = std::find_if(peratom.begin(), peratom.end(),
                            [&name] (const PerAtom &p) { return p.name == name; });

  if (match != peratom.end()) (*match).cols = cols;
  else error->all(FLERR,"Could not find per-atom array name {} for column change", name);
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

void Atom::add_peratom_vary(const std::string &name, void *address,
                            int datatype, int *cols, void *length, int collength)
{
  PerAtom item = {name, address, length, cols, datatype, -1, collength, 0};
  peratom.push_back(item);
}

/* ----------------------------------------------------------------------
   add info for a single per-atom array to PerAtom data struct
------------------------------------------------------------------------- */

void Atom::set_atomflag_defaults()
{
  // --------------------------------------------------------------------
  // 3rd customization section: customize by adding new flag
  // identical list as 2nd customization in atom.h

  labelmapflag = 0;
  sphere_flag = ellipsoid_flag = line_flag = tri_flag = body_flag = 0;
  quat_flag = 0;
  peri_flag = electron_flag = 0;
  wavepacket_flag = sph_flag = 0;
  molecule_flag = molindex_flag = molatom_flag = 0;
  q_flag = mu_flag = 0;
  rmass_flag = radius_flag = omega_flag = torque_flag = angmom_flag = 0;
  temperature_flag = heatflow_flag = 0;
  vfrac_flag = spin_flag = eradius_flag = ervel_flag = erforce_flag = 0;
  cs_flag = csforce_flag = vforce_flag = ervelforce_flag = etag_flag = 0;
  rho_flag = esph_flag = cv_flag = vest_flag = 0;
  dpd_flag = edpd_flag = tdpd_flag = 0;
  sp_flag = 0;
  x0_flag = 0;
  smd_flag = damage_flag = 0;
  mesont_flag = 0;
  contact_radius_flag = smd_data_9_flag = smd_stress_flag = 0;
  eff_plastic_strain_flag = eff_plastic_strain_rate_flag = 0;
  nspecial15_flag = 0;

  pdscale = 1.0;
}

/* ----------------------------------------------------------------------
   create an AtomVec style
   called from lammps.cpp, input script, restart file, replicate
------------------------------------------------------------------------- */

void Atom::create_avec(const std::string &style, int narg, char **arg, int trysuffix)
{
  delete[] atom_style;
  if (avec) delete avec;
  atom_style = nullptr;
  avec = nullptr;

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
    std::string estyle = style + "/";
    if (sflag == 1) estyle += lmp->suffix;
    else if (sflag == 2) estyle += lmp->suffix2;
    else if (sflag == 3) estyle += lmp->non_pair_suffix();
    atom_style = utils::strdup(estyle);
  } else {
    atom_style = utils::strdup(style);
  }

  // if molecular system:
  // atom IDs must be defined
  // force atom map to be created
  // map style will be reset to array vs hash to by map_init()

  molecular = avec->molecular;
  if ((molecular != Atom::ATOMIC) && (tag_enable == 0))
    error->all(FLERR,"Atom IDs must be used for molecular systems");
  if (molecular != Atom::ATOMIC) map_style = MAP_YES;
}

/* ----------------------------------------------------------------------
   generate an AtomVec class, first with suffix appended
------------------------------------------------------------------------- */

AtomVec *Atom::new_avec(const std::string &style, int trysuffix, int &sflag)
{
  if (trysuffix && lmp->suffix_enable) {
    if (lmp->non_pair_suffix()) {
      sflag = 1 + 2*lmp->pair_only_flag;
      std::string estyle = style + "/" + lmp->non_pair_suffix();
      if (avec_map->find(estyle) != avec_map->end()) {
        AtomVecCreator &avec_creator = (*avec_map)[estyle];
        return avec_creator(lmp);
      }
    }

    if (lmp->suffix2) {
      sflag = 2;
      std::string estyle = style + "/" + lmp->suffix2;
      if (avec_map->find(estyle) != avec_map->end()) {
        AtomVecCreator &avec_creator = (*avec_map)[estyle];
        return avec_creator(lmp);
      }
    }
  }

  sflag = 0;
  if (avec_map->find(style) != avec_map->end()) {
    AtomVecCreator &avec_creator = (*avec_map)[style];
    return avec_creator(lmp);
  }

  error->all(FLERR,utils::check_packages_for_style("atom",style,lmp));
  return nullptr;
}

/* ---------------------------------------------------------------------- */

void Atom::init()
{
  // delete extra array since it doesn't persist past first run

  if (nextra_store) {
    memory->destroy(extra);
    extra = nullptr;
    nextra_store = 0;
  }

  // check arrays that are atom type in length

  check_mass(FLERR);

  // setup of firstgroup

  if (firstgroupname) {
    firstgroup = group->find(firstgroupname);
    if (firstgroup < 0)
      error->all(FLERR,"Could not find atom_modify first group ID {}", firstgroupname);
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

/* ---------------------------------------------------------------------- */

std::string Atom::get_style()
{
  std::string retval = atom_style;
  if (retval == "hybrid") {
    auto avec_hybrid = dynamic_cast<AtomVecHybrid *>(avec);
    if (avec_hybrid) {
      for (int i = 0; i < avec_hybrid->nstyles; i++) {
        retval += ' ';
        retval += avec_hybrid->keywords[i];
      }
    }
  }
  return retval;
}

/* ----------------------------------------------------------------------
   return ptr to AtomVec class if matches style or to matching hybrid sub-class
   return nullptr if no match
------------------------------------------------------------------------- */

AtomVec *Atom::style_match(const char *style)
{
  if (strcmp(atom_style,style) == 0) return avec;
  else if (strcmp(atom_style,"hybrid") == 0) {
    auto avec_hybrid = dynamic_cast<AtomVecHybrid *>(avec);
    for (int i = 0; i < avec_hybrid->nstyles; i++)
      if (strcmp(avec_hybrid->keywords[i],style) == 0)
        return avec_hybrid->styles[i];
  }
  return nullptr;
}

/* ----------------------------------------------------------------------
   modify parameters of the atom style
   some options can only be invoked before simulation box is defined
   first and sort options cannot be used together
------------------------------------------------------------------------- */

void Atom::modify_params(int narg, char **arg)
{
  if (narg == 0) utils::missing_cmd_args(FLERR, "atom_modify", error);

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"id") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "atom_modify id", error);
      if (domain->box_exist)
        error->all(FLERR,"Atom_modify id command after simulation box is defined");
      tag_enable = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"map") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "atom_modify map", error);
      if (domain->box_exist)
        error->all(FLERR,"Atom_modify map command after simulation box is defined");
      if (strcmp(arg[iarg+1],"array") == 0) map_user = MAP_ARRAY;
      else if (strcmp(arg[iarg+1],"hash") == 0) map_user = MAP_HASH;
      else if (strcmp(arg[iarg+1],"yes") == 0) map_user = MAP_YES;
      else error->all(FLERR,"Illegal atom_modify map command argument {}", arg[iarg+1]);
      map_style = map_user;
      iarg += 2;
    } else if (strcmp(arg[iarg],"first") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "atom_modify first", error);
      if (strcmp(arg[iarg+1],"all") == 0) {
        delete[] firstgroupname;
        firstgroupname = nullptr;
      } else {
        firstgroupname = utils::strdup(arg[iarg+1]);
        sortfreq = 0;
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"sort") == 0) {
      if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "atom_modify sort", error);
      sortfreq = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      userbinsize = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      if (sortfreq < 0) error->all(FLERR,"Illegal atom_modify sort frequency {}", sortfreq);
      if (userbinsize < 0.0) error->all(FLERR,"Illegal atom_modify sort bin size {}", userbinsize);
      if ((sortfreq >= 0) && firstgroupname)
        error->all(FLERR,"Atom_modify sort and first options cannot be used together");
      iarg += 3;
    } else error->all(FLERR,"Illegal atom_modify command argument: {}", arg[iarg]);
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

  if (minall < 0) error->all(FLERR,"One or more Atom IDs are negative");
  if (maxall >= MAXTAGINT) error->all(FLERR,"One or more atom IDs are too big");
  if (maxall > 0 && minall == 0)
    error->all(FLERR,"One or more atom IDs are zero");
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
    error->all(FLERR,"New atom IDs exceed maximum allowed ID {}", MAXTAGINT);

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
    if (ellipsoid && (ellipsoid[i] >= 0)) ++local_ellipsoids;
    if (line && (line[i] >= 0)) ++local_lines;
    if (tri && (tri[i] >= 0)) ++local_tris;
    if (body && (body[i] >= 0)) ++local_bodies;
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
  atom->bond_type = nullptr;
  atom->bond_atom = nullptr;

  memory->destroy(atom->angle_type);
  memory->destroy(atom->angle_atom1);
  memory->destroy(atom->angle_atom2);
  memory->destroy(atom->angle_atom3);
  atom->angle_type = nullptr;
  atom->angle_atom1 = atom->angle_atom2 = atom->angle_atom3 = nullptr;

  memory->destroy(atom->dihedral_type);
  memory->destroy(atom->dihedral_atom1);
  memory->destroy(atom->dihedral_atom2);
  memory->destroy(atom->dihedral_atom3);
  memory->destroy(atom->dihedral_atom4);
  atom->dihedral_type = nullptr;
  atom->dihedral_atom1 = atom->dihedral_atom2 =
    atom->dihedral_atom3 = atom->dihedral_atom4 = nullptr;

  memory->destroy(atom->improper_type);
  memory->destroy(atom->improper_atom1);
  memory->destroy(atom->improper_atom2);
  memory->destroy(atom->improper_atom3);
  memory->destroy(atom->improper_atom4);
  atom->improper_type = nullptr;
  atom->improper_atom1 = atom->improper_atom2 =
    atom->improper_atom3 = atom->improper_atom4 = nullptr;
}

/* ----------------------------------------------------------------------
   unpack N lines from Atom section of data file
   call style-specific routine to parse line
------------------------------------------------------------------------- */

void Atom::data_atoms(int n, char *buf, tagint id_offset, tagint mol_offset,
                      int type_offset, int shiftflag, double *shift,
                      int labelflag, int *ilabel)
{
  int xptr,iptr;
  imageint imagedata;
  double xdata[3],lamda[3];
  double *coord;
  char *next;
  std::string typestr;
  auto location = "Atoms section of data file";

  // use the first line to detect and validate the number of words/tokens per line
  next = strchr(buf,'\n');
  if (!next) error->all(FLERR, "Missing data in {}", location);
  *next = '\0';
  auto values = Tokenizer(buf).as_vector();
  int nwords = values.size();
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (utils::strmatch(values[i], "^#")) {
      nwords = i;
      break;
    }
  }

  if ((nwords != avec->size_data_atom) && (nwords != avec->size_data_atom + 3))
    error->all(FLERR,"Incorrect format in {}: {}", location, utils::trim(buf));

  *next = '\n';
  // set bounds for my proc
  // if periodic and I am lo/hi proc, adjust bounds by EPSILON
  // ensures all data atoms will be owned even with round-off

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
    if (!next) error->all(FLERR, "Missing data in {}", location);
    *next = '\0';
    auto values = Tokenizer(buf).as_vector();
    int nvalues = values.size();
    if ((nvalues == 0) || (utils::strmatch(values[0],"^#.*")))  {
      // skip over empty or comment lines
    } else if ((nvalues < nwords) ||
               ((nvalues > nwords) && (!utils::strmatch(values[nwords],"^#")))) {
      error->all(FLERR, "Incorrect format in {}: {}", location, utils::trim(buf));
    } else {
      int imx = 0, imy = 0, imz = 0;
      if (imageflag) {
        imx = utils::inumeric(FLERR,values[iptr],false,lmp);
        imy = utils::inumeric(FLERR,values[iptr+1],false,lmp);
        imz = utils::inumeric(FLERR,values[iptr+2],false,lmp);
        if ((domain->dimension == 2) && (imz != 0))
          error->all(FLERR,"Z-direction image flag must be 0 for 2d-systems");
        if ((!domain->xperiodic) && (imx != 0)) { reset_image_flag[0] = true; imx = 0; }
        if ((!domain->yperiodic) && (imy != 0)) { reset_image_flag[1] = true; imy = 0; }
        if ((!domain->zperiodic) && (imz != 0)) { reset_image_flag[2] = true; imz = 0; }
      }
      imagedata = ((imageint) (imx + IMGMAX) & IMGMASK) |
        (((imageint) (imy + IMGMAX) & IMGMASK) << IMGBITS) |
        (((imageint) (imz + IMGMAX) & IMGMASK) << IMG2BITS);

      xdata[0] = utils::numeric(FLERR,values[xptr],false,lmp);
      xdata[1] = utils::numeric(FLERR,values[xptr+1],false,lmp);
      xdata[2] = utils::numeric(FLERR,values[xptr+2],false,lmp);
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
        avec->data_atom(xdata,imagedata,values,typestr);
        typestr = utils::utf8_subst(typestr);
        if (id_offset) tag[nlocal-1] += id_offset;
        if (mol_offset) molecule[nlocal-1] += mol_offset;

        switch (utils::is_type(typestr)) {
          case 0: {    // numeric
            int itype = utils::inumeric(FLERR, typestr, true, lmp) + type_offset;
            if ((itype < 1) || (itype > ntypes))
              error->one(FLERR, "Invalid atom type {} in {}: {}", itype, location,
                         utils::trim(buf));
            type[nlocal - 1] = itype;
            if (labelflag) type[nlocal - 1] = ilabel[itype - 1];
            break;
          }
          case 1: {    // type label
            if (!labelmapflag)
              error->one(FLERR, "Invalid line in {}: {}", location, utils::trim(buf));
            type[nlocal - 1] = lmap->find(typestr, Atom::ATOM);
            if (type[nlocal - 1] == -1)
              error->one(FLERR, "Invalid line in {}: {}", location, utils::trim(buf));
            break;
          }
          default:    // invalid
            error->one(FLERR, "Invalid line in {}: {}", location, utils::trim(buf));
            break;
        }

        if (type[nlocal-1] <= 0 || type[nlocal-1] > ntypes)
          error->one(FLERR,"Invalid atom type {} in {}", location, typestr);
      }
    }
    buf = next + 1;
  }
}

/* ----------------------------------------------------------------------
   unpack N lines from Velocity section of data file
   check that atom IDs are > 0 and <= map_tag_max
   call style-specific routine to parse line
------------------------------------------------------------------------- */

void Atom::data_vels(int n, char *buf, tagint id_offset)
{
  int m;
  char *next;

  // loop over lines of atom velocities
  // tokenize the line into values
  // if I own atom tag, unpack its values

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    if (!next) error->all(FLERR, "Missing data in Velocities section of data file");
    *next = '\0';
    auto values = Tokenizer(utils::trim_comment(buf)).as_vector();
    if (values.size() == 0) {
      // skip over empty or comment lines
    } else if ((int)values.size() != avec->size_data_vel) {
      error->all(FLERR, "Incorrect velocity format in data file: {}", utils::trim(buf));
    } else {
      tagint tagdata = utils::tnumeric(FLERR,values[0],false,lmp) + id_offset;
      if (tagdata <= 0 || tagdata > map_tag_max)
        error->one(FLERR,"Invalid atom ID {} in Velocities section of data file: {}", tagdata, buf);
      if ((m = map(tagdata)) >= 0) avec->data_vel(m,values);
    }
    buf = next + 1;
  }
}

/* ----------------------------------------------------------------------
   process N bonds read into buf from data files
   if count is non-nullptr, just count bonds per atom
   else store them with atoms
   check that atom IDs are > 0 and <= map_tag_max
------------------------------------------------------------------------- */

void Atom::data_bonds(int n, char *buf, int *count, tagint id_offset,
                      int type_offset, int labelflag, int *ilabel)
{
  int m,itype;
  tagint atom1,atom2;
  char *next;
  std::string typestr;
  int newton_bond = force->newton_bond;
  auto location = "Bonds section of data file";

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    if (!next) error->all(FLERR, "Missing data in {}", location);
    *next = '\0';
    auto values = Tokenizer(buf).as_vector();
    int nwords = values.size();
    for (std::size_t ii = 0; ii < values.size(); ++ii) {
      if (utils::strmatch(values[ii], "^#")) {
        nwords = ii;
        break;
      }
    }

    // skip over empty or comment lines
    // Bonds line is: number(ignored), bond type, atomID 1, atomID 2
    if (nwords > 0) {
      if (nwords != 4) error->all(FLERR, "Incorrect format in {}: {}", location, utils::trim(buf));
      typestr = utils::utf8_subst(values[1]);
      atom1 = utils::tnumeric(FLERR, values[2], false, lmp);
      atom2 = utils::tnumeric(FLERR, values[3], false, lmp);
      if (id_offset) {
        atom1 += id_offset;
        atom2 += id_offset;
      }

      switch (utils::is_type(typestr)) {
        case 0: {    // numeric
          itype = utils::inumeric(FLERR, typestr, false, lmp) + type_offset;
          if ((itype < 1) || (itype > nbondtypes))
            error->all(FLERR, "Invalid bond type {} in {}: {}", itype, location, utils::trim(buf));
          if (labelflag) itype = ilabel[itype - 1];
          break;
        }
        case 1: {    // type label
          if (!atom->labelmapflag) error->all(FLERR, "Invalid {}: {}", location, utils::trim(buf));
          itype = lmap->find(typestr, Atom::BOND);
          if (itype == -1) error->all(FLERR, "Invalid {}: {}", location, utils::trim(buf));
          break;
        }
        default:    // invalid
          error->one(FLERR, "Invalid {}: {}", location, utils::trim(buf));
          break;
      }

      if ((atom1 <= 0) || (atom1 > map_tag_max) ||
          (atom2 <= 0) || (atom2 > map_tag_max) || (atom1 == atom2))
        error->all(FLERR,"Invalid atom ID in {}: {}", location, utils::trim(buf));
      if ((itype <= 0) || (itype > nbondtypes))
        error->all(FLERR,"Invalid bond type {} in {}: {}", itype, location, utils::trim(buf));
      if ((m = map(atom1)) >= 0) {
        if (count) count[m]++;
        else {
          bond_type[m][num_bond[m]] = itype;
          bond_atom[m][num_bond[m]] = atom2;
          num_bond[m]++;
          avec->data_bonds_post(m, num_bond[m], atom1, atom2, id_offset);
        }
      }
      if (newton_bond == 0) {
        if ((m = map(atom2)) >= 0) {
          if (count) count[m]++;
          else {
            bond_type[m][num_bond[m]] = itype;
            bond_atom[m][num_bond[m]] = atom1;
            num_bond[m]++;
            avec->data_bonds_post(m, num_bond[m], atom1, atom2, id_offset);
          }
        }
      }
    }
    buf = next + 1;
  }
}

/* ----------------------------------------------------------------------
   process N angles read into buf from data files
   if count is non-nullptr, just count angles per atom
   else store them with atoms
   check that atom IDs are > 0 and <= map_tag_max
------------------------------------------------------------------------- */

void Atom::data_angles(int n, char *buf, int *count, tagint id_offset,
                       int type_offset, int labelflag, int *ilabel)
{
  int m,itype;
  tagint atom1,atom2,atom3;
  char *next;
  std::string typestr;
  int newton_bond = force->newton_bond;
  auto location = "Angles section of data file";

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    if (!next) error->all(FLERR, "Missing data in {}", location);
    *next = '\0';
    auto values = Tokenizer(buf).as_vector();
    int nwords = values.size();
    for (std::size_t ii = 0; ii < values.size(); ++ii) {
      if (utils::strmatch(values[ii], "^#")) {
        nwords = ii;
        break;
      }
    }

    // skip over empty or comment lines
    // Angles line is: number(ignored), angle type, atomID 1, atomID 2, atomID 3
    if (nwords > 0) {
      if (nwords != 5) error->all(FLERR, "Incorrect format in {}: {}", location, utils::trim(buf));
      typestr = utils::utf8_subst(values[1]);
      atom1 = utils::tnumeric(FLERR, values[2], false, lmp);
      atom2 = utils::tnumeric(FLERR, values[3], false, lmp);
      atom3 = utils::tnumeric(FLERR, values[4], false, lmp);
      if (id_offset) {
        atom1 += id_offset;
        atom2 += id_offset;
        atom3 += id_offset;
      }

      switch (utils::is_type(typestr)) {
        case 0: {    // numeric
          itype = utils::inumeric(FLERR, typestr, false, lmp) + type_offset;
          if ((itype < 1) || (itype > nangletypes))
            error->all(FLERR, "Invalid angle type {} in {}: {}", itype, location, utils::trim(buf));
          if (labelflag) itype = ilabel[itype - 1];
          break;
        }
        case 1: {    // type label
          if (!atom->labelmapflag) error->all(FLERR, "Invalid {}: {}", location, utils::trim(buf));
          itype = lmap->find(typestr, Atom::ANGLE);
          if (itype == -1) error->all(FLERR, "Invalid {}: {}", location, utils::trim(buf));
          break;
        }
        default:    // invalid
          error->one(FLERR, "Invalid {}: {}", location, utils::trim(buf));
          break;
      }

      if ((atom1 <= 0) || (atom1 > map_tag_max) ||
          (atom2 <= 0) || (atom2 > map_tag_max) ||
          (atom3 <= 0) || (atom3 > map_tag_max) ||
          (atom1 == atom2) || (atom1 == atom3) || (atom2 == atom3))
        error->one(FLERR,"Invalid atom ID in {}: {}", location, utils::trim(buf));
      if (itype <= 0 || itype > nangletypes)
        error->one(FLERR,"Invalid angle type {} in {}: {}", itype, location, utils::trim(buf));
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
    }
    buf = next + 1;
  }
}

/* ----------------------------------------------------------------------
   process N dihedrals read into buf from data files
   if count is non-nullptr, just count diihedrals per atom
   else store them with atoms
   check that atom IDs are > 0 and <= map_tag_max
------------------------------------------------------------------------- */

void Atom::data_dihedrals(int n, char *buf, int *count, tagint id_offset,
                          int type_offset, int labelflag, int *ilabel)
{
  int m,itype;
  tagint atom1,atom2,atom3,atom4;
  char *next;
  std::string typestr;
  int newton_bond = force->newton_bond;
  auto location = "Dihedrals section of data file";

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    if (!next) error->all(FLERR, "Missing data in {}", location);
    *next = '\0';
    auto values = Tokenizer(buf).as_vector();
    int nwords = values.size();
    for (std::size_t ii = 0; ii < values.size(); ++ii) {
      if (utils::strmatch(values[ii], "^#")) {
        nwords = ii;
        break;
      }
    }

    // skip over empty or comment lines
    // Dihedrals line is: number(ignored), bond type, atomID 1, atomID 2, atomID 3, atomID 4
    if (nwords > 0) {
      if (nwords != 6) error->all(FLERR, "Incorrect format in {}: {}", location, utils::trim(buf));
      typestr = utils::utf8_subst(values[1]);
      atom1 = utils::tnumeric(FLERR, values[2], false, lmp);
      atom2 = utils::tnumeric(FLERR, values[3], false, lmp);
      atom3 = utils::tnumeric(FLERR, values[4], false, lmp);
      atom4 = utils::tnumeric(FLERR, values[5], false, lmp);
      if (id_offset) {
        atom1 += id_offset;
        atom2 += id_offset;
        atom3 += id_offset;
        atom4 += id_offset;
      }

      switch (utils::is_type(typestr)) {
        case 0: {    // numeric
          itype = utils::inumeric(FLERR, typestr, false, lmp) + type_offset;
          if ((itype < 1) || (itype > ndihedraltypes))
            error->all(FLERR, "Invalid dihedral type {} in {}: {}", itype, location,
                       utils::trim(buf));
          if (labelflag) itype = ilabel[itype - 1];
          break;
        }
        case 1: {    // type label
          if (!atom->labelmapflag) error->all(FLERR, "Invalid {}: {}", location, utils::trim(buf));
          itype = lmap->find(typestr, Atom::DIHEDRAL);
          if (itype == -1) error->all(FLERR, "Invalid {}: {}", location, utils::trim(buf));
          break;
        }
        default:    // invalid
          error->one(FLERR, "Invalid {}: {}", location, utils::trim(buf));
          break;
      }

      if ((atom1 <= 0) || (atom1 > map_tag_max) ||
          (atom2 <= 0) || (atom2 > map_tag_max) ||
          (atom3 <= 0) || (atom3 > map_tag_max) ||
          (atom4 <= 0) || (atom4 > map_tag_max) ||
          (atom1 == atom2) || (atom1 == atom3) || (atom1 == atom4) ||
          (atom2 == atom3) || (atom2 == atom4) || (atom3 == atom4))
        error->one(FLERR, "Invalid atom ID in {}: {}", location, utils::trim(buf));
      if (itype <= 0 || itype > ndihedraltypes)
        error->one(FLERR, "Invalid dihedral type {} in {}: {}", itype, location, utils::trim(buf));
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
    }
    buf = next + 1;
  }
}

/* ----------------------------------------------------------------------
   process N impropers read into buf from data files
   if count is non-nullptr, just count impropers per atom
   else store them with atoms
   check that atom IDs are > 0 and <= map_tag_max
------------------------------------------------------------------------- */

void Atom::data_impropers(int n, char *buf, int *count, tagint id_offset,
                          int type_offset, int labelflag, int *ilabel)
{
  int m,itype;
  tagint atom1,atom2,atom3,atom4;
  char *next;
  std::string typestr;
  int newton_bond = force->newton_bond;
  auto location = "Impropers section of data file";

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    if (!next) error->all(FLERR, "Missing data in {}", location);
    *next = '\0';
    auto values = Tokenizer(buf).as_vector();
    int nwords = values.size();
    for (std::size_t ii = 0; ii < values.size(); ++ii) {
      if (utils::strmatch(values[ii], "^#")) {
        nwords = ii;
        break;
      }
    }

    // skip over empty or comment lines
    // Impropers line is: number(ignored), bond type, atomID 1, atomID 2, atomID 3, atomID 4
    if (nwords > 0) {
      if (nwords != 6) error->all(FLERR, "Incorrect format in {}: {}", location, utils::trim(buf));
      typestr = utils::utf8_subst(values[1]);
      atom1 = utils::tnumeric(FLERR, values[2], false, lmp);
      atom2 = utils::tnumeric(FLERR, values[3], false, lmp);
      atom3 = utils::tnumeric(FLERR, values[4], false, lmp);
      atom4 = utils::tnumeric(FLERR, values[5], false, lmp);
      if (id_offset) {
        atom1 += id_offset;
        atom2 += id_offset;
        atom3 += id_offset;
        atom4 += id_offset;
      }

      switch (utils::is_type(typestr)) {
        case 0: {    // numeric
          itype = utils::inumeric(FLERR, typestr, false, lmp) + type_offset;
          if ((itype < 1) || (itype > nimpropertypes))
            error->all(FLERR, "Invalid improper type {} in {}: {}", itype, location,
                       utils::trim(buf));
          if (labelflag) itype = ilabel[itype - 1];
          break;
        }
        case 1: {    // type label
          if (!atom->labelmapflag) error->all(FLERR, "Invalid {}: {}", location, utils::trim(buf));
          itype = lmap->find(typestr, Atom::IMPROPER);
          if (itype == -1) error->all(FLERR, "Invalid {}: {}", location, utils::trim(buf));
          break;
        }
        default:    // invalid
          error->one(FLERR, "Invalid {}: {}", location, utils::trim(buf));
          break;
      }

      if ((atom1 <= 0) || (atom1 > map_tag_max) ||
          (atom2 <= 0) || (atom2 > map_tag_max) ||
          (atom3 <= 0) || (atom3 > map_tag_max) ||
          (atom4 <= 0) || (atom4 > map_tag_max) ||
          (atom1 == atom2) || (atom1 == atom3) || (atom1 == atom4) ||
          (atom2 == atom3) || (atom2 == atom4) || (atom3 == atom4))
        error->one(FLERR, "Invalid atom ID in {}: {}", location, utils::trim(buf));
      if (itype <= 0 || itype > nimpropertypes)
        error->one(FLERR, "Invalid improper type {} in {}: {}", itype, location, utils::trim(buf));
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
  int m;
  char *next;

  // loop over lines of bonus atom data
  // tokenize the line into values
  // if I own atom tag, unpack its values

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    if (!next) error->all(FLERR, "Missing data in Bonus section of data file");
    *next = '\0';
    auto values = Tokenizer(utils::trim_comment(buf)).as_vector();
    if (values.size() == 0) {
      // skip over empty or comment lines
    } else if ((int)values.size() != avec_bonus->size_data_bonus) {
      error->all(FLERR, "Incorrect bonus data format in data file: {}", utils::trim(buf));
    } else {
      tagint tagdata = utils::tnumeric(FLERR,values[0],false,lmp) + id_offset;
      if (tagdata <= 0 || tagdata > map_tag_max)
        error->one(FLERR,"Invalid atom ID in Bonus section of data file");

      // ok to call child's data_atom_bonus() method thru parent avec_bonus,
      // since data_bonus() was called with child ptr, and method is virtual

      if ((m = map(tagdata)) >= 0) avec_bonus->data_atom_bonus(m,values);
    }
    buf = next + 1;
  }
}

/* ----------------------------------------------------------------------
   unpack N bodies from Bodies section of data file
   each body spans multiple lines
   check that atom IDs are > 0 and <= map_tag_max
   call style-specific routine to parse line
------------------------------------------------------------------------- */

void Atom::data_bodies(int n, char *buf, AtomVec *avec_body, tagint id_offset)
{
  std::vector<int> ivalues;
  std::vector<double> dvalues;

  if (!unique_tags) unique_tags = new std::set<tagint>;

  // loop over lines of body data
  // if I own atom tag, tokenize lines into ivalues/dvalues, call data_body()
  // else skip values

  for (int i = 0; i < n; i++) {
    char *next = strchr(buf,'\n');
    if (!next) error->all(FLERR, "Missing data in Bodies section of data file");
    *next = '\0';

    auto values = Tokenizer(utils::trim_comment(buf)).as_vector();
    if (values.size()) {
      tagint tagdata = utils::tnumeric(FLERR,values[0],false,lmp) + id_offset;
      int ninteger = utils::inumeric(FLERR,values[1],false,lmp);
      int ndouble = utils::inumeric(FLERR,values[2],false,lmp);

      if (unique_tags->find(tagdata) == unique_tags->end())
        unique_tags->insert(tagdata);
      else
        error->one(FLERR,"Duplicate atom ID {} in Bodies section of data file", tagdata);

      buf = next + 1;
      int m = map(tagdata);
      if (m >= 0) {
        ivalues.resize(ninteger);
        dvalues.resize(ndouble);

        for (int j = 0; j < ninteger; j++) {
          buf += strspn(buf," \t\n\r\f");
          buf[strcspn(buf," \t\n\r\f")] = '\0';
          ivalues[j] = utils::inumeric(FLERR,buf,false,lmp);
          buf += strlen(buf)+1;
        }

        for (int j = 0; j < ndouble; j++) {
          buf += strspn(buf," \t\n\r\f");
          buf[strcspn(buf," \t\n\r\f")] = '\0';
          dvalues[j] = utils::numeric(FLERR,buf,false,lmp);
          buf += strlen(buf)+1;
        }

        avec_body->data_body(m,ninteger,ndouble,ivalues.data(),dvalues.data());

      } else {
        int nvalues = ninteger + ndouble;    // number of values to skip
        for (int j = 0; j < nvalues; j++) {
          buf += strspn(buf," \t\n\r\f");
          buf[strcspn(buf," \t\n\r\f")] = '\0';
          buf += strlen(buf)+1;
        }
      }
    }
    buf += strspn(buf," \t\n\r\f");
  }
}

/* ----------------------------------------------------------------------
   init per-atom fix/compute/variable values for newly created atoms
   called from create_atoms, read_data, read_dump,
     lib::lammps_create_atoms()
   fixes, computes, variables may or may not exist when called
------------------------------------------------------------------------- */

void Atom::data_fix_compute_variable(int nprev, int nnew)
{
  for (const auto &ifix : modify->get_fix_list()) {
    if (ifix->create_attribute)
      for (int i = nprev; i < nnew; i++)
        ifix->set_arrays(i);
  }

  for (const auto &icompute : modify->get_compute_list()) {
    if (icompute->create_attribute)
      for (int i = nprev; i < nnew; i++)
        icompute->set_arrays(i);
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
  if (avec->mass_type == AtomVec::PER_TYPE) {
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
// clang-format on

void Atom::set_mass(const char *file, int line, const char *str, int type_offset, int labelflag,
                    int *ilabel)
{
  if (mass == nullptr) error->all(file, line, "Cannot set mass for atom style {}", atom_style);

  int itype;
  double mass_one;
  auto location = "Masses section of data file";
  auto values = Tokenizer(str).as_vector();
  int nwords = values.size();
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (utils::strmatch(values[i], "^#")) {
      nwords = i;
      break;
    }
  }
  if (nwords != 2) error->all(file, line, "Invalid format in {}: {}", location, str);
  auto typestr = utils::utf8_subst(values[0]);

  switch (utils::is_type(typestr)) {

    case 0: {    // numeric
      itype = utils::inumeric(file, line, typestr, false, lmp);
      if ((itype < 1) || (itype > ntypes))
        error->all(file, line, "Invalid atom type {} in {}: {}", itype, location, utils::trim(str));
      if (labelflag) itype = ilabel[itype - 1];
      break;
    }

    case 1: {    // type label
      if (!atom->labelmapflag)
        error->all(file, line, "Invalid atom type in {}: {}", location, utils::trim(str));
      itype = lmap->find(typestr, Atom::ATOM);
      if (itype == -1)
        error->all(file, line, "Unknown atom type {} in {}: {}", typestr, location,
                   utils::trim(str));
      break;
    }

    default:    // invalid
      itype = -1000000000;
      error->one(file, line, "Invalid {}: {}", location, utils::trim(str));
      break;
  }
  itype += type_offset;
  mass_one = utils::numeric(file, line, values[1], false, lmp);

  if (itype < 1 || itype > ntypes)
    error->all(file, line, "Invalid atom type {} in {}: {}", itype, location, utils::trim(str));

  if (mass_one <= 0.0)
    error->all(file, line, "Invalid mass value {} in {}: {}", mass_one, location, utils::trim(str));
  mass[itype] = mass_one;
  mass_setflag[itype] = 1;
}
// clang-format off

/* ----------------------------------------------------------------------
   set a mass and flag it as set
   called from EAM pair routine
------------------------------------------------------------------------- */

void Atom::set_mass(const char *file, int line, int itype, double value)
{
  if (mass == nullptr)
    error->all(file,line, "Cannot set per-type mass for atom style {}", atom_style);
  if (itype < 1 || itype > ntypes)
    error->all(file,line,"Invalid type {} for atom mass {}", itype, value);
  if (value <= 0.0) {
    if (comm->me == 0)
      error->warning(file,line,"Ignoring invalid mass value {} for atom type {}", value, itype);
    return;
  }
  mass[itype] = value;
  mass_setflag[itype] = 1;
}

/* ----------------------------------------------------------------------
   set one or more masses and flag them as set
   called from reading of input script
------------------------------------------------------------------------- */

void Atom::set_mass(const char *file, int line, int /*narg*/, char **arg)
{
  if (mass == nullptr)
    error->all(file,line, "Cannot set per-type atom mass for atom style {}", atom_style);

  char *typestr = utils::expand_type(file, line, arg[0], Atom::ATOM, lmp);
  const std::string str = typestr ? typestr : arg[0];
  delete[] typestr;

  int lo, hi;
  utils::bounds(file, line, str, 1, ntypes, lo, hi, error);
  if ((lo < 1) || (hi > ntypes))
    error->all(file, line, "Invalid atom type {} for atom mass", str);

  const double value = utils::numeric(FLERR, arg[1], false, lmp);
  if (value <= 0.0)
    error->all(file, line, "Invalid atom mass value {} for type {}", value, str);

  for (int itype = lo; itype <= hi; itype++) {
    mass[itype] = value;
    mass_setflag[itype] = 1;
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
  if (mass == nullptr) return;
  if (rmass_flag) return;
  for (int itype = 1; itype <= ntypes; itype++)
    if (mass_setflag[itype] == 0)
      error->all(file,line,"Not all per-type masses are set. Type {} is missing.", itype);
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

int Atom::shape_consistency(int itype, double &shapex, double &shapey, double &shapez)
{
  double zero[3] = {0.0, 0.0, 0.0};
  double one[3] = {-1.0, -1.0, -1.0};
  double *shape;

  auto avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(style_match("ellipsoid"));
  auto bonus = avec_ellipsoid->bonus;

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] != itype) continue;
    if (ellipsoid[i] < 0) shape = zero;
    else shape = bonus[ellipsoid[i]].shape;

    if (one[0] < 0.0) {
      one[0] = shape[0];
      one[1] = shape[1];
      one[2] = shape[2];
    } else if ((one[0] != shape[0]) || (one[1] != shape[1]) || (one[2] != shape[2]))
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
  if (narg < 1) utils::missing_cmd_args(FLERR, "molecule", error);

  if (find_molecule(arg[0]) >= 0)
    error->all(FLERR,"Reuse of molecule template ID {}", arg[0]);

  // 1st molecule in set stores nset = # of mols, others store nset = 0
  // ifile = count of molecules in set
  // index = argument index where next molecule starts, updated by constructor

  int ifile = 1;
  int index = 1;
  while (true) {
    molecules = (Molecule **)
      memory->srealloc(molecules,(nmolecule+1)*sizeof(Molecule *), "atom::molecules");
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

int Atom::find_molecule(const char *id)
{
  if (id == nullptr) return -1;
  for (int imol = 0; imol < nmolecule; imol++)
    if (strcmp(id,molecules[imol]->id) == 0) return imol;
  return -1;
}

/* ----------------------------------------------------------------------
   return vector of molecules which match template ID
------------------------------------------------------------------------- */

std::vector<Molecule *>Atom::get_molecule_by_id(const std::string &id)
{
  std::vector<Molecule *> result;
  for (int imol = 0; imol < nmolecule; ++imol)
    if (id == molecules[imol]->id) result.push_back(molecules[imol]);
  return result;
}

/* ----------------------------------------------------------------------
   add info to current atom ilocal from molecule template onemol and its iatom
   offset = atom ID preceding IDs of atoms in this molecule
   called by fixes and commands that add molecules
------------------------------------------------------------------------- */

void Atom::add_molecule_atom(Molecule *onemol, int iatom, int ilocal, tagint offset)
{
  if (onemol->qflag && q_flag) q[ilocal] = onemol->q[iatom];
  if (onemol->radiusflag && radius_flag) radius[ilocal] = onemol->radius[iatom];
  if (onemol->rmassflag && rmass_flag) rmass[ilocal] = onemol->rmass[iatom];
  else if (rmass_flag)
    rmass[ilocal] = 4.0*MY_PI/3.0 * radius[ilocal]*radius[ilocal]*radius[ilocal];
  if (onemol->bodyflag) {
    body[ilocal] = 0;     // as if a body read from data file
    onemol->avec_body->data_body(ilocal,onemol->nibody,onemol->ndbody,
                                 onemol->ibodyparams,onemol->dbodyparams);
    onemol->avec_body->set_quat(ilocal,onemol->quat_external);
  }

  // initialize custom per-atom properties to zero if present

  for (int i = 0; i < nivector; ++i)
    if (ivname[i]) ivector[i][ilocal] = 0;
  for (int i = 0; i < ndvector; ++i)
    if (dvname[i]) dvector[i][ilocal] = 0.0;
  for (int i = 0; i < niarray; ++i)
    if (ianame[i])
      for (int j = 0; j < icols[i]; ++j)
        iarray[i][ilocal][j] = 0;
  for (int i = 0; i < ndarray; ++i)
    if (daname[i])
      for (int j = 0; j < dcols[i]; ++j)
        darray[i][ilocal][j] = 0.0;

  if (molecular != Atom::MOLECULAR) return;

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
   allocate space for type label map
------------------------------------------------------------------------- */

void Atom::add_label_map()
{
  if (lmp->kokkos)
    error->all(FLERR, "Label maps are currently not supported with Kokkos");
  labelmapflag = 1;
  lmap = new LabelMap(lmp,ntypes,nbondtypes,nangletypes,ndihedraltypes,nimpropertypes);
}

/* ----------------------------------------------------------------------
   reorder owned atoms so those in firstgroup appear first
   called by comm->exchange() if atom_modify first group is set
   only owned atoms exist at this point, no ghost atoms
------------------------------------------------------------------------- */

void Atom::first_reorder()
{
  // ensure there is one extra atom location at end of arrays for swaps

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

  // ensure there is one extra atom location at end of arrays for swaps

  if (nlocal == nmax) avec->grow(0);

  // bin atoms in reverse order so linked list will be in forward order

  for (i = 0; i < nbins; i++) binhead[i] = -1;

  // for triclinic, atoms must be in box coords (not lamda) to match bbox

  if (domain->triclinic) domain->lamda2x(nlocal);

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

  // convert back to lamda coords

  if (domain->triclinic) domain->x2lamda(nlocal);

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
      error->warning(FLERR,"No pairwise cutoff or binsize set. Atom sorting therefore disabled.");
    return;
  }

#ifdef LMP_GPU
  if (userbinsize == 0.0) {
    auto ifix = dynamic_cast<FixGPU *>(modify->get_fix_by_id("package_gpu"));
    if (ifix) {
      const double subx = domain->subhi[0] - domain->sublo[0];
      const double suby = domain->subhi[1] - domain->sublo[1];
      const double subz = domain->subhi[2] - domain->sublo[2];
      binsize = ifix->binsize(subx, suby, subz, atom->nlocal, 0.5 * neighbor->cutneighmax);
    }
  }
#endif

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

#ifdef LMP_INTEL
  if (neighbor->has_intel_request() && userbinsize == 0.0) {
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

#ifdef LMP_GPU
  if (userbinsize == 0.0) {
    FixGPU *fix = dynamic_cast<FixGPU *>(modify->get_fix_by_id("package_gpu"));
    if (fix) {
      const double subx = domain->subhi[0] - domain->sublo[0];
      const double suby = domain->subhi[1] - domain->sublo[1];
      const double subz = domain->subhi[2] - domain->sublo[2];

      binsize = fix->binsize(subx, suby, subz, atom->nlocal,neighbor->cutneighmax);
      bininv = 1.0 / binsize;

      nbinx = static_cast<int> (ceil(subx * bininv));
      nbiny = static_cast<int> (ceil(suby * bininv));
      nbinz = static_cast<int> (ceil(subz * bininv));
      if (domain->dimension == 2) nbinz = 1;

      if (nbinx == 0) nbinx = 1;
      if (nbiny == 0) nbiny = 1;
      if (nbinz == 0) nbinz = 1;

      bininvx = bininv;
      bininvy = bininv;
      bininvz = bininv;
    }
  }
#endif

  if (1.0*nbinx*nbiny*nbinz > INT_MAX) error->one(FLERR,"Too many atom sorting bins");

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
   flag = Atom::GROW for grow, Atom::RESTART for restart, Atom::BORDER for border comm
------------------------------------------------------------------------- */

void Atom::add_callback(int flag)
{
  int ifix;

  // find the fix
  // if find null pointer:
  //   it's this one, since it is being replaced and has just been deleted
  //   at this point in re-creation
  // if don't find null pointer:
  //   i is set to nfix = new one currently being added at end of list

  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (modify->fix[ifix] == nullptr) break;

  // add callback to lists and sort, reallocating if necessary
  // sorting is required in cases where fixes were replaced as it ensures atom
  // data is read/written/transfered in the same order that fixes are called

  if (flag == GROW) {
    if (nextra_grow == nextra_grow_max) {
      nextra_grow_max += DELTA;
      memory->grow(extra_grow,nextra_grow_max,"atom:extra_grow");
    }
    extra_grow[nextra_grow] = ifix;
    nextra_grow++;
    std::sort(extra_grow, extra_grow + nextra_grow);
  } else if (flag == RESTART) {
    if (nextra_restart == nextra_restart_max) {
      nextra_restart_max += DELTA;
      memory->grow(extra_restart,nextra_restart_max,"atom:extra_restart");
    }
    extra_restart[nextra_restart] = ifix;
    nextra_restart++;
    std::sort(extra_restart, extra_restart + nextra_restart);
  } else if (flag == BORDER) {
    if (nextra_border == nextra_border_max) {
      nextra_border_max += DELTA;
      memory->grow(extra_border,nextra_border_max,"atom:extra_border");
    }
    extra_border[nextra_border] = ifix;
    nextra_border++;
    std::sort(extra_border, extra_border + nextra_border);
  }
}

/* ----------------------------------------------------------------------
   unregister a callback to a fix
   happens when fix is deleted, called by its destructor
   flag = 0 for grow, 1 for restart
------------------------------------------------------------------------- */

void Atom::delete_callback(const char *id, int flag)
{
  if (id == nullptr) return;

  int ifix = modify->find_fix(id);

  // compact the list of callbacks

  if (flag == GROW) {
    int match;
    for (match = 0; match < nextra_grow; match++)
      if (extra_grow[match] == ifix) break;
    if ((nextra_grow == 0) || (match == nextra_grow))
      error->all(FLERR,"Trying to delete non-existent Atom::grow() callback");
    for (int i = match; i < nextra_grow-1; i++)
      extra_grow[i] = extra_grow[i+1];
    nextra_grow--;

  } else if (flag == RESTART) {
    int match;
    for (match = 0; match < nextra_restart; match++)
      if (extra_restart[match] == ifix) break;
    if ((nextra_restart == 0) || (match == nextra_restart))
      error->all(FLERR,"Trying to delete non-existent Atom::restart() callback");
    for (int i = match; i < nextra_restart-1; i++)
      extra_restart[i] = extra_restart[i+1];
    nextra_restart--;

  } else if (flag == BORDER) {
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
   return index if found, -1 if not found
     lists of names can have NULL entries if previously removed
   return flag = 0/1 for int/double
   return cols = 0/N for vector/array where N = # of columns
------------------------------------------------------------------------- */

int Atom::find_custom(const char *name, int &flag, int &cols)
{
  if (name == nullptr) return -1;

  for (int i = 0; i < nivector; i++)
    if (ivname[i] && strcmp(ivname[i],name) == 0) {
      flag = 0;
      cols = 0;
      return i;
    }

  for (int i = 0; i < ndvector; i++)
    if (dvname[i] && strcmp(dvname[i],name) == 0) {
      flag = 1;
      cols = 0;
      return i;
    }

  for (int i = 0; i < niarray; i++)
    if (ianame[i] && strcmp(ianame[i],name) == 0) {
      flag = 0;
      cols = icols[i];
      return i;
    }

  for (int i = 0; i < ndarray; i++)
    if (daname[i] && strcmp(daname[i],name) == 0) {
      flag = 1;
      cols = dcols[i];
      return i;
    }

  return -1;
}

/** \brief Add a custom per-atom property with the given name and type and size
\verbatim embed:rst

This function will add a custom per-atom property with one or more values
with the name "name" to the list of custom properties.
This function is called, e.g. from :doc:`fix property/atom <fix_property_atom>`.
\endverbatim
 * \param name Name of the property (w/o a "i_" or "d_" or "i2_" or "d2_" prefix)
 * \param flag Data type of property: 0 for int, 1 for double
 * \param cols Number of values: 0 for a single value, 1 or more for a vector of values
 * \return index of property in the respective list of properties
 */
int Atom::add_custom(const char *name, int flag, int cols)
{
  int index = -1;

  if ((flag == 0) && (cols == 0)) {
    index = nivector;
    nivector++;
    ivname = (char **) memory->srealloc(ivname,nivector*sizeof(char *),"atom:ivname");
    ivname[index] = utils::strdup(name);
    ivector = (int **) memory->srealloc(ivector,nivector*sizeof(int *),"atom:ivector");
    memory->create(ivector[index],nmax,"atom:ivector");

  } else if ((flag == 1) && (cols == 0)) {
    index = ndvector;
    ndvector++;
    dvname = (char **) memory->srealloc(dvname,ndvector*sizeof(char *),"atom:dvname");
    dvname[index] = utils::strdup(name);
    dvector = (double **) memory->srealloc(dvector,ndvector*sizeof(double *),"atom:dvector");
    memory->create(dvector[index],nmax,"atom:dvector");

  } else if ((flag == 0) && (cols > 0)) {
    index = niarray;
    niarray++;
    ianame = (char **) memory->srealloc(ianame,niarray*sizeof(char *),"atom:ianame");
    ianame[index] = utils::strdup(name);
    iarray = (int ***) memory->srealloc(iarray,niarray*sizeof(int **),"atom:iarray");
    memory->create(iarray[index],nmax,cols,"atom:iarray");
    icols = (int *) memory->srealloc(icols,niarray*sizeof(int),"atom:icols");
    icols[index] = cols;

  } else if ((flag == 1) && (cols > 0)) {
    index = ndarray;
    ndarray++;
    daname = (char **) memory->srealloc(daname,ndarray*sizeof(char *),"atom:daname");
    daname[index] = utils::strdup(name);
    darray = (double ***) memory->srealloc(darray,ndarray*sizeof(double **),"atom:darray");
    memory->create(darray[index],nmax,cols,"atom:darray");
    dcols = (int *) memory->srealloc(dcols,ndarray*sizeof(int),"atom:dcols");
    dcols[index] = cols;
  }

  if (index < 0)
    error->all(FLERR,"Invalid call to Atom::add_custom()");
  return index;
}

/*! \brief Remove a custom per-atom property of a given type and size
 *
\verbatim embed:rst
This will remove a property that was requested, e.g. by the
:doc:`fix property/atom <fix_property_atom>` command.  It frees the
allocated memory and sets the pointer to ``nullptr`` for the entry in
the list so it can be reused. The lists of these pointers are never
compacted or shrink, so that indices to name mappings remain valid.
\endverbatim
 * \param index Index of property in the respective list of properties
 * \param flag Data type of property: 0 for int, 1 for double
 * \param cols Number of values: 0 for a single value, 1 or more for a vector of values
 */
void Atom::remove_custom(int index, int flag, int cols)
{
  if (flag == 0 && cols == 0) {
    memory->destroy(ivector[index]);
    ivector[index] = nullptr;
    delete[] ivname[index];
    ivname[index] = nullptr;

  } else if (flag == 1 && cols == 0) {
    memory->destroy(dvector[index]);
    dvector[index] = nullptr;
    delete[] dvname[index];
    dvname[index] = nullptr;

  } else if (flag == 0 && cols) {
    memory->destroy(iarray[index]);
    iarray[index] = nullptr;
    delete[] ianame[index];
    ianame[index] = nullptr;

  } else if (flag == 1 && cols) {
    memory->destroy(darray[index]);
    darray[index] = nullptr;
    delete[] daname[index];
    daname[index] = nullptr;
  }
}

/** Provide access to internal data of the Atom class by keyword
 *
\verbatim embed:rst

This function is a way to access internal per-atom data.  This data is
distributed across MPI ranks and thus only the data for "local" atoms
are expected to be available.  Whether also data for "ghost" atoms is
stored and up-to-date depends on various simulation settings.

This table lists a large part of the supported names, their data types,
length of the data area, and a short description.

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Name
     - Type
     - Items per atom
     - Description
   * - mass
     - double
     - 1
     - per-type mass. This array is **NOT** a per-atom array
       but of length ``ntypes+1``, element 0 is ignored.
   * - id
     - tagint
     - 1
     - atom ID of the particles
   * - type
     - int
     - 1
     - atom type of the particles
   * - mask
     - int
     - 1
     - bitmask for mapping to groups. Individual bits are set
       to 0 or 1 for each group.
   * - image
     - imageint
     - 1
     - 3 image flags encoded into a single integer.
       See :cpp:func:`lammps_encode_image_flags`.
   * - x
     - double
     - 3
     - x-, y-, and z-coordinate of the particles
   * - v
     - double
     - 3
     - x-, y-, and z-component of the velocity of the particles
   * - f
     - double
     - 3
     - x-, y-, and z-component of the force on the particles
   * - molecule
     - int
     - 1
     - molecule ID of the particles
   * - q
     - double
     - 1
     - charge of the particles
   * - mu
     - double
     - 3
     - dipole moment of the particles
   * - omega
     - double
     - 3
     - x-, y-, and z-component of rotational velocity of the particles
   * - angmom
     - double
     - 3
     - x-, y-, and z-component of angular momentum of the particles
   * - torque
     - double
     - 3
     - x-, y-, and z-component of the torque on the particles
   * - radius
     - double
     - 1
     - radius of the (extended) particles
   * - rmass
     - double
     - 1
     - per-atom mass of the particles. ``nullptr`` if per-type masses are
       used. See the :cpp:func:`rmass_flag<lammps_extract_setting>` setting.
   * - ellipsoid
     - int
     - 1
     - 1 if the particle is an ellipsoidal particle, 0 if not
   * - line
     - int
     - 1
     - 1 if the particle is a line particle, 0 if not
   * - tri
     - int
     - 1
     - 1 if the particle is a triangulated particle, 0 if not
   * - body
     - int
     - 1
     - 1 if the particle is a body particle, 0 if not
   * - quat
     - double
     - 4
     - four quaternion components of the particles
   * - temperature
     - double
     - 1
     - temperature of the particles
   * - heatflow
     - double
     - 1
     - heatflow of the particles
   * - i_name
     - int
     - 1
     - single integer value defined by fix property/atom vector name
   * - d_name
     - double
     - 1
     - single double value defined by fix property/atom vector name
   * - i2_name
     - int
     - n
     - N integer values defined by fix property/atom array name
   * - d2_name
     - double
     - n
     - N double values defined by fix property/atom array name

*See also*
   :cpp:func:`lammps_extract_atom`

\endverbatim
 *
 * \sa extract_datatype
 *
 * \param  name  string with the keyword of the desired property.
                 Typically the name of the pointer variable returned
 * \return       pointer to the requested data cast to ``void *`` or ``nullptr`` */

void *Atom::extract(const char *name)
{
  // --------------------------------------------------------------------
  // 4th customization section: customize by adding new variable name
  // if new variable is from a package, add package comment

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
  if (strcmp(name,"body") == 0) return (void *) body;
  if (strcmp(name,"quat") == 0) return (void *) quat;
  if (strcmp(name,"temperature") == 0) return (void *) temperature;
  if (strcmp(name,"heatflow") == 0) return (void *) heatflow;

  // PERI PACKAGE

  if (strcmp(name,"vfrac") == 0) return (void *) vfrac;
  if (strcmp(name,"s0") == 0) return (void *) s0;
  if (strcmp(name,"x0") == 0) return (void *) x0;

  // SPIN PACKAGE

  if (strcmp(name,"sp") == 0) return (void *) sp;

  // EFF and AWPMD packages

  if (strcmp(name,"espin") == 0) return (void *) spin;
  if (strcmp(name,"spin") == 0) return (void *) spin;  // backward compatibility
  if (strcmp(name,"eradius") == 0) return (void *) eradius;
  if (strcmp(name,"ervel") == 0) return (void *) ervel;
  if (strcmp(name,"erforce") == 0) return (void *) erforce;
  if (strcmp(name,"ervelforce") == 0) return (void *) ervelforce;
  if (strcmp(name,"cs") == 0) return (void *) cs;
  if (strcmp(name,"csforce") == 0) return (void *) csforce;
  if (strcmp(name,"vforce") == 0) return (void *) vforce;
  if (strcmp(name,"etag") == 0) return (void *) etag;

  // SPH package
  if (strcmp(name,"rho") == 0) return (void *) rho;
  if (strcmp(name,"drho") == 0) return (void *) drho;
  if (strcmp(name,"esph") == 0) return (void *) esph;
  if (strcmp(name,"desph") == 0) return (void *) desph;
  if (strcmp(name,"cv") == 0) return (void *) cv;
  if (strcmp(name,"vest") == 0) return (void *) vest;

  // MACHDYN package

  if (strcmp(name, "contact_radius") == 0) return (void *) contact_radius;
  if (strcmp(name, "smd_data_9") == 0) return (void *) smd_data_9;
  if (strcmp(name, "smd_stress") == 0) return (void *) smd_stress;
  if (strcmp(name, "eff_plastic_strain") == 0)
    return (void *) eff_plastic_strain;
  if (strcmp(name, "eff_plastic_strain_rate") == 0)
    return (void *) eff_plastic_strain_rate;
  if (strcmp(name, "damage") == 0) return (void *) damage;

  // DPD-REACT pakage

  if (strcmp(name,"dpdTheta") == 0) return (void *) dpdTheta;

  // DPD-MESO package

  if (strcmp(name,"edpd_temp") == 0) return (void *) edpd_temp;

  // DIELECTRIC package

  if (strcmp(name,"area") == 0) return (void *) area;
  if (strcmp(name,"ed") == 0) return (void *) ed;
  if (strcmp(name,"em") == 0) return (void *) em;
  if (strcmp(name,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(name,"curvature") == 0) return (void *) curvature;
  if (strcmp(name,"q_scaled") == 0) return (void *) q_scaled;

  // end of customization section
  // --------------------------------------------------------------------

  // custom vectors and arrays

  if (utils::strmatch(name,"^[id]2?_")) {
    int which = 0, array = 0;
    if (name[0] == 'd') which = 1;
    if (name[1] == '2') array = 1;

    int index,flag,cols;
    if (!array) index = find_custom(&name[2],flag,cols);
    else index = find_custom(&name[3],flag,cols);

    if (index < 0) return nullptr;
    if (which != flag) return nullptr;
    if ((!array && cols) || (array && !cols)) return nullptr;

    if (!which && !array) return (void *) ivector[index];
    if (which && !array) return (void *) dvector[index];
    if (!which && array) return (void *) iarray[index];
    if (which && array) return (void *) darray[index];
  }

  return nullptr;
}

/** Provide data type info about internal data of the Atom class
 *
\verbatim embed:rst

.. versionadded:: 18Sep2020

\endverbatim
 *
 * \sa extract
 *
 * \param  name  string with the keyword of the desired property.
 * \return       data type constant for desired property or -1 */

int Atom::extract_datatype(const char *name)
{
  // --------------------------------------------------------------------
  // 5th customization section: customize by adding new variable name

  if (strcmp(name,"mass") == 0) return LAMMPS_DOUBLE;

  if (strcmp(name,"id") == 0) return LAMMPS_TAGINT;
  if (strcmp(name,"type") == 0) return LAMMPS_INT;
  if (strcmp(name,"mask") == 0) return LAMMPS_INT;
  if (strcmp(name,"image") == 0) return LAMMPS_TAGINT;
  if (strcmp(name,"x") == 0) return LAMMPS_DOUBLE_2D;
  if (strcmp(name,"v") == 0) return LAMMPS_DOUBLE_2D;
  if (strcmp(name,"f") == 0) return LAMMPS_DOUBLE_2D;
  if (strcmp(name,"molecule") == 0) return LAMMPS_TAGINT;
  if (strcmp(name,"q") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"mu") == 0) return LAMMPS_DOUBLE_2D;
  if (strcmp(name,"omega") == 0) return LAMMPS_DOUBLE_2D;
  if (strcmp(name,"angmom") == 0) return LAMMPS_DOUBLE_2D;
  if (strcmp(name,"torque") == 0) return LAMMPS_DOUBLE_2D;
  if (strcmp(name,"radius") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"rmass") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"ellipsoid") == 0) return LAMMPS_INT;
  if (strcmp(name,"line") == 0) return LAMMPS_INT;
  if (strcmp(name,"tri") == 0) return LAMMPS_INT;
  if (strcmp(name,"body") == 0) return LAMMPS_INT;
  if (strcmp(name,"quat") == 0) return LAMMPS_DOUBLE_2D;
  if (strcmp(name,"temperature") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"heatflow") == 0) return LAMMPS_DOUBLE;

  if (strcmp(name,"vfrac") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"s0") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"x0") == 0) return LAMMPS_DOUBLE_2D;

  if (strcmp(name,"espin") == 0) return LAMMPS_INT;
  if (strcmp(name,"spin") == 0) return LAMMPS_INT;   // backwards compatibility
  if (strcmp(name,"eradius") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"ervel") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"erforce") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"ervelforce") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"cs") == 0) return LAMMPS_DOUBLE_2D;
  if (strcmp(name,"csforce") == 0) return LAMMPS_DOUBLE_2D;
  if (strcmp(name,"vforce") == 0) return LAMMPS_DOUBLE_2D;
  if (strcmp(name,"etag") == 0) return LAMMPS_INT;

  if (strcmp(name,"rho") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"drho") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"esph") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"desph") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"cv") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"vest") == 0) return LAMMPS_DOUBLE_2D;

  // MACHDYN package

  if (strcmp(name, "contact_radius") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name, "smd_data_9") == 0) return LAMMPS_DOUBLE_2D;
  if (strcmp(name, "smd_stress") == 0) return LAMMPS_DOUBLE_2D;
  if (strcmp(name, "eff_plastic_strain") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name, "eff_plastic_strain_rate") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name, "damage") == 0) return LAMMPS_DOUBLE;

  // DPD-REACT package

  if (strcmp(name,"dpdTheta") == 0) return LAMMPS_DOUBLE;

  // DPD-MESO package

  if (strcmp(name,"edpd_temp") == 0) return LAMMPS_DOUBLE;

  // DIELECTRIC package

  if (strcmp(name,"area") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"ed") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"em") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"epsilon") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"curvature") == 0) return LAMMPS_DOUBLE;
  if (strcmp(name,"q_unscaled") == 0) return LAMMPS_DOUBLE;

  // end of customization section
  // --------------------------------------------------------------------

  // custom vectors and arrays

  if (utils::strmatch(name,"^[id]2?_")) {
    int which = 0, array = 0;
    if (name[0] == 'd') which = 1;
    if (name[1] == '2') array = 1;

    int index,flag,cols;
    if (!array) index = find_custom(&name[2],flag,cols);
    else index = find_custom(&name[3],flag,cols);

    if (index < 0) return -1;
    if (which != flag) return -1;
    if ((!array && cols) || (array && !cols)) return -1;

    if (which == 0) return LAMMPS_INT;
    else return LAMMPS_DOUBLE;
  }

  return -1;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   call to avec tallies per-atom vectors
   add in global to local mapping storage
------------------------------------------------------------------------- */

double Atom::memory_usage()
{
  double bytes = avec->memory_usage();

  bytes += (double)max_same*sizeof(int);
  if (map_style == MAP_ARRAY)
    bytes += memory->usage(map_array,map_maxarray);
  else if (map_style == MAP_HASH) {
    bytes += (double)map_nbucket*sizeof(int);
    bytes += (double)map_nhash*sizeof(HashElem);
  }
  if (maxnext) {
    bytes += memory->usage(next,maxnext);
    bytes += memory->usage(permute,maxnext);
  }

  return bytes;
}
