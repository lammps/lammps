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
Contributing Author: Jacob Gissinger (jgissing@stevens.edu)
------------------------------------------------------------------------- */

#include "reaction_constraints.h"

#include "reaction_parser.h"
#include "superpose3d.h"

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "random_mars.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
evaluate constraints: return 0 if any aren't satisfied
------------------------------------------------------------------------- */

ReactionConstraints::ReactionConstraints(LAMMPS *lmp) : Pointers (lmp)
{
  // reaction functions used by 'custom' constraint
  nrxnfunction = 3;
  rxnfunclist.resize(nrxnfunction);
  peratomflag.resize(nrxnfunction);
  rxnfunclist[0] = "rxnsum";
  peratomflag[0] = 1;
  rxnfunclist[1] = "rxnave";
  peratomflag[1] = 1;
  rxnfunclist[2] = "rxnbond";
  peratomflag[2] = 0;
  nvvec = 0;
  ncustomvars = 0;
  vvec = nullptr;
}

ReactionConstraints::~ReactionConstraints()
{
 if (vvec != nullptr) memory->destroy(vvec);
}

int ReactionConstraints::check(Reaction &rxn, std::vector<tagint> &glove)
{
  double x1[3],x2[3],x3[3],x4[3];
  double delx,dely,delz,rsq;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double rsq1,rsq2,r1,r2,c,t,prrhob;
  // for computation of dihedrals
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
  double ax,ay,az,bx,by,bz,rasq,rbsq,rgsq,rg,ra2inv,rb2inv,rabinv;
  double s,phi;
  int ANDgate;

  tagint atom1,atom2;
  double **x = atom->x;

  for (auto &constraint : rxn.constraints) constraint.satisfied = true;

  for (auto &constraint : rxn.constraints) {
    if (constraint.type == Reaction::Constraint::Type::DISTANCE) {
      get_IDcoords(constraint.idtypes[0], constraint.ids[0], x1, rxn.reactant, glove);
      get_IDcoords(constraint.idtypes[1], constraint.ids[1], x2, rxn.reactant, glove);
      delx = x1[0] - x2[0];
      dely = x1[1] - x2[1];
      delz = x1[2] - x2[2];
      domain->minimum_image(FLERR, delx,dely,delz); // ghost location fix
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < constraint.distance.rminsq || rsq > constraint.distance.rmaxsq) constraint.satisfied = false;
    } else if (constraint.type == Reaction::Constraint::Type::ANGLE) {
      get_IDcoords(constraint.idtypes[0], constraint.ids[0], x1, rxn.reactant, glove);
      get_IDcoords(constraint.idtypes[1], constraint.ids[1], x2, rxn.reactant, glove);
      get_IDcoords(constraint.idtypes[2], constraint.ids[2], x3, rxn.reactant, glove);

      // 1st bond
      delx1 = x1[0] - x2[0];
      dely1 = x1[1] - x2[1];
      delz1 = x1[2] - x2[2];
      domain->minimum_image(FLERR, delx1,dely1,delz1); // ghost location fix
      rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
      r1 = sqrt(rsq1);

      // 2nd bond
      delx2 = x3[0] - x2[0];
      dely2 = x3[1] - x2[1];
      delz2 = x3[2] - x2[2];
      domain->minimum_image(FLERR, delx2,dely2,delz2); // ghost location fix
      rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
      r2 = sqrt(rsq2);

      // angle (cos and sin)
      c = delx1*delx2 + dely1*dely2 + delz1*delz2;
      c /= r1*r2;
      if (c > 1.0) c = 1.0;
      if (c < -1.0) c = -1.0;
      if (acos(c) < constraint.angle.amin || acos(c) > constraint.angle.amax) constraint.satisfied = false;
    } else if (constraint.type == Reaction::Constraint::Type::DIHEDRAL) {
      // phi calculation from dihedral style harmonic
      get_IDcoords(constraint.idtypes[0], constraint.ids[0], x1, rxn.reactant, glove);
      get_IDcoords(constraint.idtypes[1], constraint.ids[1], x2, rxn.reactant, glove);
      get_IDcoords(constraint.idtypes[2], constraint.ids[2], x3, rxn.reactant, glove);
      get_IDcoords(constraint.idtypes[3], constraint.ids[3], x4, rxn.reactant, glove);

      vb1x = x1[0] - x2[0];
      vb1y = x1[1] - x2[1];
      vb1z = x1[2] - x2[2];
      domain->minimum_image(FLERR, vb1x,vb1y,vb1z);

      vb2x = x3[0] - x2[0];
      vb2y = x3[1] - x2[1];
      vb2z = x3[2] - x2[2];
      domain->minimum_image(FLERR, vb2x,vb2y,vb2z);

      vb2xm = -vb2x;
      vb2ym = -vb2y;
      vb2zm = -vb2z;
      domain->minimum_image(FLERR, vb2xm,vb2ym,vb2zm);

      vb3x = x4[0] - x3[0];
      vb3y = x4[1] - x3[1];
      vb3z = x4[2] - x3[2];
      domain->minimum_image(FLERR, vb3x,vb3y,vb3z);

      ax = vb1y*vb2zm - vb1z*vb2ym;
      ay = vb1z*vb2xm - vb1x*vb2zm;
      az = vb1x*vb2ym - vb1y*vb2xm;
      bx = vb3y*vb2zm - vb3z*vb2ym;
      by = vb3z*vb2xm - vb3x*vb2zm;
      bz = vb3x*vb2ym - vb3y*vb2xm;

      rasq = ax*ax + ay*ay + az*az;
      rbsq = bx*bx + by*by + bz*bz;
      rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
      rg = sqrt(rgsq);

      ra2inv = rb2inv = 0.0;
      if (rasq > 0) ra2inv = 1.0/rasq;
      if (rbsq > 0) rb2inv = 1.0/rbsq;
      rabinv = sqrt(ra2inv*rb2inv);

      c = (ax*bx + ay*by + az*bz)*rabinv;
      s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

      if (c > 1.0) c = 1.0;
      if (c < -1.0) c = -1.0;
      phi = atan2(s,c);

      ANDgate = 0;
      if (constraint.dihedral.amin < constraint.dihedral.amax) {
        if (phi > constraint.dihedral.amin && phi < constraint.dihedral.amax) ANDgate = 1;
      } else {
        if (phi > constraint.dihedral.amin || phi < constraint.dihedral.amax) ANDgate = 1;
      }
      if (constraint.dihedral.amin2 < constraint.dihedral.amax2) {
        if (phi > constraint.dihedral.amin2 && phi < constraint.dihedral.amax2) ANDgate = 1;
      } else {
        if (phi > constraint.dihedral.amin2 || phi < constraint.dihedral.amax2) ANDgate = 1;
      }
      if (ANDgate != 1) constraint.satisfied = false;
    } else if (constraint.type == Reaction::Constraint::Type::ARRHENIUS) {
      std::vector<tagint> myglove(glove.begin(), glove.begin() + rxn.reactant->natoms);
      t = get_temperature(myglove);
      prrhob = constraint.arrhenius.A*pow(t,constraint.arrhenius.n)*
        exp(-constraint.arrhenius.E_a/(force->boltz*t));
      if (prrhob < constraint.arrhenius.rrhandom->uniform()) constraint.satisfied = false;
    } else if (constraint.type == Reaction::Constraint::Type::RMSD) {
      // call superpose
      int iatom;
      int iref = -1; // choose first atom as reference
      int n2superpose = 0;
      double **xfrozen; // coordinates for the "frozen" target molecule
      double **xmobile; // coordinates for the "mobile" molecule
      int ifragment = constraint.ids[0];
      if (ifragment >= 0) {
        for (int j = 0; j < rxn.reactant->natoms; j++)
          if (rxn.reactant->fragmentmask[ifragment][j]) n2superpose++;
        memory->create(xfrozen,n2superpose,3,"bond/react:xfrozen");
        memory->create(xmobile,n2superpose,3,"bond/react:xmobile");
        int myincr = 0;
        for (int j = 0; j < rxn.reactant->natoms; j++) {
          if (rxn.reactant->fragmentmask[ifragment][j]) {
            iatom = atom->map(glove[j]);
            if (iref == -1) iref = iatom;
            iatom = domain->closest_image(iref,iatom);
            for (int k = 0; k < 3; k++) {
              xfrozen[myincr][k] = x[iatom][k];
              xmobile[myincr][k] = rxn.reactant->x[j][k];
            }
            myincr++;
          }
        }
      } else {
        int iatom;
        int iref = -1; // choose first atom as reference
        n2superpose = rxn.reactant->natoms;
        memory->create(xfrozen,n2superpose,3,"bond/react:xfrozen");
        memory->create(xmobile,n2superpose,3,"bond/react:xmobile");
        for (int j = 0; j < n2superpose; j++) {
          iatom = atom->map(glove[j]);
          if (iref == -1) iref = iatom;
          iatom = domain->closest_image(iref,iatom);
          for (int k = 0; k < 3; k++) {
            xfrozen[j][k] = x[iatom][k];
            xmobile[j][k] = rxn.reactant->x[j][k];
          }
        }
      }
      Superpose3D<double, double **> superposer(n2superpose);
      double rmsd = superposer.Superpose(xfrozen, xmobile);
      memory->destroy(xfrozen);
      memory->destroy(xmobile);
      if (rmsd > constraint.rmsd.rmsdmax) constraint.satisfied = false;
    } else if (constraint.type == Reaction::Constraint::Type::CUSTOM) {
      constraint.satisfied = custom_constraint(constraint.custom.str, rxn, glove); // NOLINT
    }
  }

  if (!rxn.constraints.empty()) {
    std::string evalstr = rxn.constraintstr;
    for (auto &constraint : rxn.constraints) {
      evalstr.replace(evalstr.find('C'), 1, constraint.satisfied ? "1" : "0");
    }
    std::vector<char> buffer(evalstr.begin(), evalstr.end());
    buffer.push_back('\0');
    double verdict = input->variable->evaluate_boolean(buffer.data());
    if (verdict == 0.0) return 0;
  }

  // let's also check chirality within 'check' routine

  // full special lists - may need correction for unusual special bond settings
  int **nxspecial = atom->nspecial;
  tagint **xspecial = atom->special;

  for (int i = 0; i < rxn.reactant->natoms; i++) {
    if (rxn.atoms[i].chiral[0] == 1) {
      double my4coords[12] = {0.0};
      // already ensured, by transitive property, that chiral simulation atom has four neighs
      for (int j = 0; j < 4; j++) {
        atom1 = atom->map(glove[i]);
        // loop over known types involved in chiral center
        for (int jj = 0; jj < 4; jj++) {
          if (atom->type[atom->map(xspecial[atom1][j])] == rxn.atoms[i].chiral[jj+2]) {
            atom2 = atom->map(xspecial[atom1][j]);
            atom2 = domain->closest_image(atom1,atom2);
            for (int k = 0; k < 3; k++) {
              my4coords[3*jj+k] = x[atom2][k];
            }
            break;
          }
        }
      }
      if (ReactionParser(lmp).get_chirality(my4coords) != rxn.atoms[i].chiral[1]) return 0;
    }
  }

  return 1;
}

/* ----------------------------------------------------------------------
return pre-reaction atom or fragment location
fragment: given pre-reacted molID (reactant) and fragID,
          return geometric center (of mapped simulation atoms)
------------------------------------------------------------------------- */

void ReactionConstraints::get_IDcoords(Reaction::Constraint::IDType idtype, int myID,
                                double *center, Molecule *mol, std::vector<tagint> &glove)
{
  double **x = atom->x;
  if (idtype == Reaction::Constraint::IDType::ATOM) {
    int iatom = atom->map(glove[myID-1]);
    for (int i = 0; i < 3; i++)
      center[i] = x[iatom][i];
  } else {
    int iref = -1; // choose first atom as reference
    int iatom;
    int nfragatoms = 0;
    for (int i = 0; i < 3; i++)
      center[i] = 0;

    for (int i = 0; i < mol->natoms; i++) {
      if (mol->fragmentmask[myID][i]) {
        if (iref == -1) iref = atom->map(glove[i]);
        iatom = atom->map(glove[i]);
        iatom = domain->closest_image(iref,iatom);
        for (int j = 0; j < 3; j++)
          center[j] += x[iatom][j];
        nfragatoms++;
      }
    }
    if (nfragatoms > 0)
      for (int i = 0; i < 3; i++) center[i] /= nfragatoms;
  }
}

/* ----------------------------------------------------------------------
compute local temperature: average over all atoms in reaction template
------------------------------------------------------------------------- */

double ReactionConstraints::get_temperature(std::vector<tagint> &glove)
{
  double adof = domain->dimension;

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;

  double t = 0.0;

  if (rmass) {
    for (const auto &g : glove) {
      auto ilocal = atom->map(g);
      t += (v[ilocal][0]*v[ilocal][0] + v[ilocal][1]*v[ilocal][1] +
            v[ilocal][2]*v[ilocal][2]) * rmass[ilocal];
    }
  } else {
    for (const auto &g : glove) {
      auto ilocal = atom->map(g);
      t += (v[ilocal][0]*v[ilocal][0] + v[ilocal][1]*v[ilocal][1] +
            v[ilocal][2]*v[ilocal][2]) * mass[type[ilocal]];
    }
  }

  // final temperature
  double dof = adof*glove.size();
  double tfactor = force->mvv2e / (dof * force->boltz);
  t *= tfactor;
  return t;
}

/* ----------------------------------------------------------------------
evaulate expression for variable constraint
------------------------------------------------------------------------- */

bool ReactionConstraints::custom_constraint(const std::string &varstr, Reaction &rxn, std::vector<tagint> &glove)
{
  std::size_t pos,pos1,pos2,pos3;
  int irxnfunc;
  int prev3 = -1;
  std::string argstr,varid,fragid,evlcat;
  std::vector<std::string> evlstr;

  // search varstr for special 'rxn' functions
  while (true) {
    // find next reaction special function occurrence
    pos1 = std::string::npos;
    for (int i = 0; i < nrxnfunction; i++) {
      pos = varstr.find(rxnfunclist[i],prev3+1);
      if (pos == std::string::npos) continue;
      if (pos < pos1) {
        pos1 = pos;
        irxnfunc = i;
      }
    }
    if (pos1 == std::string::npos) break;

    fragid = "all"; // operate over entire reaction site by default
    pos2 = varstr.find('(',pos1);
    pos3 = varstr.find(')',pos2);
    if (pos2 == std::string::npos || pos3 == std::string::npos)
      error->one(FLERR,"Fix bond/react: Illegal rxn function syntax\n");
    evlstr.push_back(varstr.substr(prev3+1,pos1-(prev3+1)));
    prev3 = pos3;
    argstr = varstr.substr(pos2+1,pos3-pos2-1);
    argstr.erase(remove_if(argstr.begin(), argstr.end(), isspace), argstr.end()); // remove whitespace
    pos2 = argstr.find(',');
    if (pos2 != std::string::npos) {
      varid = argstr.substr(0,pos2);
      fragid = argstr.substr(pos2+1);
    } else varid = argstr;
    evlstr.push_back(std::to_string(rxnfunction(rxnfunclist[irxnfunc], varid, fragid, rxn.reactant, glove)));
  }
  evlstr.push_back(varstr.substr(prev3+1));

  for (auto & evl : evlstr) evlcat += evl;
  return static_cast<bool>(input->variable->compute_equal(evlcat));
}

/* ----------------------------------------------------------------------
currently three 'rxn' functions: rxnsum, rxnave, and rxnbond
------------------------------------------------------------------------- */

double ReactionConstraints::rxnfunction(const std::string& rxnfunc, const std::string& varid,
                                 const std::string& fragid, Molecule *mol, std::vector<tagint> &glove)
{
  int ifrag = -1;
  if (fragid != "all") {
    ifrag = mol->findfragment(fragid.c_str());
    if (ifrag < 0) error->one(FLERR,"Bond/react: Molecule fragment "
                              "in reaction special function does not exist");
  }

  // start with 'rxnbond' per-bond function
  // for 'rxnbond', varid corresponds to 'compute bond/local' name,
  //                and fragid is a pre-reaction fragment containing the two atoms in the bond
  if (rxnfunc == "rxnbond") {
    int ibond;
    double perbondval;
    std::set<tagint> aset;
    std::string computeid = varid;

    if (computeid.substr(0,2) != "c_") error->one(FLERR,"Bond/react: Reaction special function compute "
                                         "name should begin with 'c_'");
    computeid = computeid.substr(2);
    cperbond = modify->get_compute_by_id(computeid);
    if (!cperbond) error->one(FLERR,"Bond/react: Reaction special function compute name does not exist");
    std::string compute_style = cperbond->style;
    if (compute_style != "bond/local") error->one(FLERR,"Bond/react: Compute used by reaction "
                                         "special function 'rxnbond' must be of style 'bond/local'");
    if (cperbond->size_local_cols > 0) error->one(FLERR,"Bond/react: 'Compute bond/local' used by reaction "
                                         "special function 'rxnbond' must compute one value");

    if (atoms2bondflag == 0) {
      atoms2bondflag = 1;
      get_atoms2bond(cperbond->groupbit);
    }

    for (int i = 0; i < mol->natoms; i++) {
      if (mol->fragmentmask[ifrag][i]) {
        aset.insert(glove[i]);
      }
    }
    if (aset.size() != 2) error->one(FLERR,"Bond/react: Molecule fragment of reaction special function 'rxnbond' "
                     "must contain exactly two atoms");

    if (cperbond->invoked_local != lmp->update->ntimestep)
      cperbond->compute_local();

    auto it = atoms2bond.find(aset);
    if (it == atoms2bond.end()) error->one(FLERR,"Bond/react: Unable to locate bond referenced by "
                                            "reaction special function 'rxnbond'");
    ibond = it->second;
    perbondval = cperbond->vector_local[ibond];
    return perbondval;
  }

  int ivar = -1;
  for (int i = 0; i < ncustomvars; i++) {
    if (varid == customvarstrs[i]) {
      ivar = i;
      break;
    }
  }
  // variable name should always be found, at this point
  // however, let's double check for completeness
  if (ivar < 0)
    error->one(FLERR,"Fix bond/react: Reaction special function variable "
                                 "name does not exist");

  int iatom;
  int nsum = 0;
  double sumvvec = 0;
  if (rxnfunc == "rxnsum" || rxnfunc == "rxnave") {
    if (fragid == "all") {
      for (int i = 0; i < mol->natoms; i++) {
        iatom = atom->map(glove[i]);
        sumvvec += vvec[iatom][ivar];
      }
      nsum = mol->natoms;
    } else {
      for (int i = 0; i < mol->natoms; i++) {
        if (mol->fragmentmask[ifrag][i]) {
          iatom = atom->map(glove[i]);
          sumvvec += vvec[iatom][ivar];
          nsum++;
        }
      }
    }
  }

  if (rxnfunc == "rxnsum") return sumvvec;
  if (rxnfunc == "rxnave") return sumvvec/nsum;
  return 0.0;
}

/* ----------------------------------------------------------------------
populate map to get bond index from atom IDs
------------------------------------------------------------------------- */

void ReactionConstraints::get_atoms2bond(int cgroupbit)
{
  int i,m,atom1,atom2,btype,nb;
  std::set<tagint> aset;

  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *mask = atom->mask;

  m = 0;
  atoms2bond.clear();
  for (atom1 = 0; atom1 < nlocal; atom1++) {
    if (!(mask[atom1] & cgroupbit)) continue;
    nb = num_bond[atom1];
    for (i = 0; i < nb; i++) {
      btype = bond_type[atom1][i];
      atom2 = atom->map(bond_atom[atom1][i]);
      if (atom2 < 0 || !(mask[atom2] & cgroupbit)) continue;
      if (force->newton_bond == 0 && tag[atom1] > tag[atom2]) continue;
      if (btype == 0) continue;
      aset = {tag[atom1], tag[atom2]};
      atoms2bond.insert(std::make_pair(aset,m++));
    }
  }
}

/* ----------------------------------------------------------------------
get per-atom variable names used by custom constraint
initialize RNG for Arrhenius constraint
------------------------------------------------------------------------- */

void ReactionConstraints::customvarnames(std::vector<Reaction> & rxns)
{
  std::size_t pos,pos1,pos2,pos3;
  int prev3;
  std::string varstr,argstr,varid;

  // search all constraints' varstr for special 'rxn' functions
  //   add variable names to customvarstrs
  //   add values to customvars

  for (auto &rxn : rxns) {
    for (auto &constraint : rxn.constraints) {
      if (constraint.type == Reaction::Constraint::Type::CUSTOM) {
        varstr = constraint.custom.str;
        prev3 = -1;
        while (true) {
          // find next reaction special function occurrence
          pos1 = std::string::npos;
          for (int i = 0; i < nrxnfunction; i++) {
            if (peratomflag[i] == 0) continue;
            pos = varstr.find(rxnfunclist[i],prev3+1);
            if (pos == std::string::npos) continue;
            if (pos < pos1) pos1 = pos;
          }
          if (pos1 == std::string::npos) break;

          pos2 = varstr.find('(',pos1);
          pos3 = varstr.find(')',pos2);
          if (pos2 == std::string::npos || pos3 == std::string::npos)
            error->all(FLERR,"Fix bond/react: Illegal rxn function syntax\n");
          prev3 = (int)pos3;
          argstr = varstr.substr(pos2+1,pos3-pos2-1);
          argstr.erase(remove_if(argstr.begin(), argstr.end(), isspace), argstr.end()); // remove whitespace
          pos2 = argstr.find(',');
          if (pos2 != std::string::npos) varid = argstr.substr(0,pos2);
          else varid = argstr;
          // check if we already know about this variable
          int varidflag = 0;
          for (int j = 0; j < ncustomvars; j++) {
            if (customvarstrs[j] == varid) {
              varidflag = 1;
              break;
            }
          }
          if (!varidflag) {
            customvarstrs.resize(ncustomvars+1);
            customvarstrs[ncustomvars++] = varid;
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
evaluate per-atom variables needed for custom constraint
------------------------------------------------------------------------- */

void ReactionConstraints::get_customvars(int igroup)
{
  double *tempvvec;
  std::string varid;
  int nall = atom->nlocal + atom->nghost;

  memory->create(tempvvec,nall,"bond/react:tempvvec");
  if (vvec == nullptr) {
    memory->create(vvec,nall,ncustomvars,"bond/react:vvec");
    nvvec = nall;
  }
  if (nvvec < nall) {
    memory->grow(vvec,nall,ncustomvars,"bond/react:vvec");
    nvvec = nall;
  }
  for (int i = 0; i < ncustomvars; i++) {
    varid = customvarstrs[i];
    if (varid.substr(0,2) != "v_") error->all(FLERR,"Fix bond/react: Reaction special function variable "
                                     "name should begin with 'v_'");
    varid = varid.substr(2);
    int ivar = input->variable->find(varid.c_str());
    if (ivar < 0)
      error->all(FLERR,"Fix bond/react: Reaction special function variable "
                                   "name does not exist");
    if (!input->variable->atomstyle(ivar))
      error->all(FLERR,"Fix bond/react: Reaction special function must "
                                   "reference an atom-style variable");

    input->variable->compute_atom(ivar,igroup,tempvvec,1,0);
    for (int j = 0; j < nall; j++) vvec[j][i] = tempvvec[j];
  }
  memory->destroy(tempvvec);
}
