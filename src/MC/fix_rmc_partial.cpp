/* -*- c++ -*- ----------------------------------------------------------
 * LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 * https://www.lammps.org/, Sandia National Laboratories
 * LAMMPS development team: developers@lammps.org
 *          
 * Copyright (2003) Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software.  This software is distributed under
 * the GNU General Public License.
 * 
 * See the README file in the top-level LAMMPS directory.
 * ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Vishnu Raghuraman (UIUC)
   Email: vishnura@illinois.edu
------------------------------------------------------------------------- */

#include "fix_rmc_partial.h"
#include "compute.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "comm.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "domain.h"
#include "error.h"
#include "modify.h"
#include "neighbor.h"
#include "json.h"
#include "text_file_reader.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixRMCPartial::FixRMCPartial(LAMMPS *lmp, int narg, char **arg) : rng_acc(device_acc()), rng_atom(device_atom()), 
rng_type_source(device_type_source()), rng_type_destination(device_type_destination()), Fix(lmp, narg, arg)
{
  perform_step = 0;
  nmcsteps=0;

  // JSON input filename
  std::string jsonfilename = arg[3];

  // Read the data from rank 0
  json rmcdata;
  std::vector <std::uint8_t> rmcdata_byte;
  int rmcdatasize = 0; 
  FILE *fp;

  if (comm->me == 0){
    if (!utils::strmatch(jsonfilename, "\\.json$")){
      error->all(FLERR, "Filename {} does not have .json extension", jsonfilename);    
    }
    fp = fopen(jsonfilename.c_str(), "r");
    rmcdata = json::parse(fp);
    rmcdata_byte = json::to_ubjson(rmcdata);
    rmcdatasize = rmcdata_byte.size();
    fclose(fp);
  }

  // The problem here is only rank 0 has the data
  // We have to send it to other ranks. First, we send the size
  MPI_Bcast(&rmcdatasize, 1, MPI_INT, 0, world);

  // Ensure all ranks have the correct size data structure for receiving data
  if (comm->me != 0) rmcdata_byte.resize(rmcdatasize);

  // Now send the actual data
  MPI_Bcast(rmcdata_byte.data(), rmcdatasize, MPI_CHAR, 0, world);

  // clear out the rmcdata variable and repopulate on all ranks from rmcdata_byte
  rmcdata.clear();
  rmcdata = json::from_ubjson(rmcdata_byte);
  rmcdata_byte.clear();

  // All ranks have the json data
  // Now we get to business - examining the json data

  // Get the system name, important for restart files
  if (rmcdata.contains("sysname")){
     sysname = std::string(rmcdata["sysname"]);
  } else {
   error->all(FLERR, "Sysname not defined in JSON input file\n");
  }

  // Let's get the temperature
  if (rmcdata.contains("temperature")){
    temperature = rmcdata["temperature"];
  }  else {
    error->all(FLERR, "Temperature not defined in JSON input file\n");
  }
  beta = 1.0/(force->boltz * temperature);

  // Determing if new run or restart
  if (rmcdata.contains("restart")) {
    if (std::string(rmcdata["restart"]) == "y") {
      restart = 1;
    } else {
      restart = 0;
    }
 } else {
   restart = 0;
 }

  // Getting dopant/semiconductor size, number of molecules, type threshold and other parameters
  if (rmcdata.contains("system")){
     if (rmcdata["system"].contains("dopant_num_atoms")){
        dopant_size = rmcdata["system"]["dopant_num_atoms"];
     } else {
       error->all(FLERR, "Dopant_num_atoms missing in JSON input system field\n");
     }
     if (rmcdata["system"].contains("semiconductor_num_atoms")){
        semiconductor_size = rmcdata["system"]["semiconductor_num_atoms"];
     } else {
       error->all(FLERR, "Semiconductor_num_atoms missing in JSON input system field\n");
     }
     size_limit = std::max(dopant_size, semiconductor_size);
     if (rmcdata["system"].contains("num_dopants")){
       n_dopant = rmcdata["system"]["num_dopants"];
     } else {
       error->all(FLERR, "num_dopants missing in JSON input system field");
     }
     if (rmcdata["system"].contains("num_semiconductors")){
       n_semiconductor = rmcdata["system"]["num_semiconductors"];
     } else {
       error->all(FLERR, "num_semiconductors missing in JSON input system field\n");
     }
     n_molecules = n_dopant+n_semiconductor;
  } else {
    error->all(FLERR, "System field missing in JSON input file\n");
  }

  if (rmcdata.contains("cycle")){
     if (rmcdata["cycle"].contains("mdsteps")){
        periodicity = rmcdata["cycle"]["mdsteps"];
     } else {
        error->all(FLERR, "Periodicity missing in JSON input cycle field\n");
     }
     if (rmcdata["cycle"].contains("mcsteps")){
        nmoves = rmcdata["cycle"]["mcsteps"];
     } else {
        error->all(FLERR, "Mcsteps missing in JSON input cycle field\n");
     }
  } else {
   error->all(FLERR, "Cycle field not defined in JSON input file\n");
  }

  if (rmcdata.contains("type_threshold")){
     type_threshold = rmcdata["type_threshold"];
  } else {
    error->all(FLERR, "Type threshold missing in JSON input file");
  }
  
  if (rmcdata.contains("num_charge_states")){
    num_charge_states = rmcdata["num_charge_states"];
  } else {
    error->all(FLERR, "Num charge states missing in JSON input file");
  }

  delta_g_list = new double [num_charge_states];
  if (rmcdata.contains("barrier")) {
    // Create and populate the delta_g list for various charge states
   if (rmcdata["barrier"].size() != num_charge_states){
      error->all(FLERR, "Number of barriers provided inconsistent with number of charge states.");
   }
    for (int i=0;i<num_charge_states;i++){
       delta_g_list[i] = rmcdata["barrier"][i];
    }
  } else {
    error->all(FLERR, "Barrier missing in JSON input file");
  }

  
  // Data structure to store semiconductor and dopant charges
  semiconductor_charges = new double*[num_charge_states];
  dopant_charges = new double*[num_charge_states];
  for (int i=0;i<num_charge_states;i++) {
     semiconductor_charges[i] = new double [semiconductor_size];
     dopant_charges[i] = new double [dopant_size];
  }
  if (rmcdata.contains("semiconductor_charges")){
    // Populate the charge 0 and charge +-1 state
    if (rmcdata["semiconductor_charges"].contains("neutral") && rmcdata["semiconductor_charges"].contains("charged")) {
      for (int i=0;i<semiconductor_size;i++){
         semiconductor_charges[0][i] = rmcdata["semiconductor_charges"]["neutral"][i];
         semiconductor_charges[num_charge_states-1][i] = rmcdata["semiconductor_charges"]["charged"][i];
      }
    } else {
      error->all(FLERR, "The neutral and/or charged values missing from JSON input semiconductor_charges field");
    }
  } else {
    error->all(FLERR, "Semiconductor charges field missing in JSON input file.");
  }

  if (rmcdata.contains("dopant_charges")){
    if (rmcdata["dopant_charges"].contains("neutral") && rmcdata["dopant_charges"].contains("charged")){
      for (int i=0;i<dopant_size;i++){
         dopant_charges[0][i] = rmcdata["dopant_charges"]["neutral"][i];
         dopant_charges[num_charge_states-1][i] = rmcdata["dopant_charges"]["charged"][i];
      }
    } else {
      error->all(FLERR, "The neutral and/or charged values missing from JSON input dopant_charges field");
    } 
 } else {
   error->all(FLERR, "Dopant charges field missing in JSON input file");
 }

  // Populate the intermediate charge values for each semiconductor atom using linear interpolation
  double temp_delta_semiconductor, temp_delta_dopant;
  
  for (int i=0;i<semiconductor_size;i++) {
     temp_delta_semiconductor = (semiconductor_charges[num_charge_states-1][i] - semiconductor_charges[0][i])/(num_charge_states-1);
     for (int j=0;j<num_charge_states-1;j++)
     {
        semiconductor_charges[j+1][i] = semiconductor_charges[j][i] + temp_delta_semiconductor;
     }
  }
  
  for (int i=0;i<dopant_size;i++) {
     temp_delta_dopant = (dopant_charges[num_charge_states-1][i] - dopant_charges[0][i])/(num_charge_states-1);
     for (int j=0;j<num_charge_states-1;j++)
     {
        dopant_charges[j+1][i] = dopant_charges[j][i] + temp_delta_dopant;
     }
  }

  // Separation between charge states
  deltaQ = 1.0/(num_charge_states-1.0);
  
  // List of charges
  charges = new double [num_charge_states];
  charges[0] = 0;
  for (int i=0;i<num_charge_states-1;i++) {
     charges[i+1] = charges[i] + deltaQ;
  }

  // Initializing random number generators
  atom_dist = std::uniform_int_distribution<> (1, n_molecules);
  type_dist = std::uniform_int_distribution<> (0,num_charge_states-1);
  acc_dist = std::uniform_real_distribution<> (0,1);
  

  // Check if we are doing dihedral/angle/bond modifications
  if (rmcdata.contains("dihedral_modification")) {
    if (std::string(rmcdata["dihedral_modification"]) == "y"){
      do_dihedral = 1;
      if (!rmcdata.contains("dihedral_types") || !rmcdata.contains("dihedral_list")){
         error->all(FLERR, "The dihedral_types and/or dihedral_list field are missing in the JSON input file");
      }
    } else {
      do_dihedral = 0;
    } 
  } else {
    do_dihedral = 0;
  }

  if (rmcdata.contains("angle_modification")) {
    if (std::string(rmcdata["angle_modification"]) == "y"){
      do_angle = 1;
      if (!rmcdata.contains("angle_types") || !rmcdata.contains("angle_list")){
         error->all(FLERR, "The angle_types and/or angle_list field are missing in the JSON input file");
      }
    } else {
      do_angle = 0;
    } 
  } else {
    do_angle = 0;
  }

  if (rmcdata.contains("bond_modification")) {
    if (std::string(rmcdata["bond_modification"]) == "y"){
      do_bond = 1;
      if (!rmcdata.contains("bond_types") || !rmcdata.contains("bond_list")){
         error->all(FLERR, "The bond_types and/or bond_list field are missing in the JSON input file");
      }
    } else {
      do_bond = 0;
    } 
  } else {
    do_bond = 0;
  }

  
  // Read in the dihedrals that need to be altered
  if (do_dihedral == 1){
     num_dihedrals = rmcdata["dihedral_list"].size();
     dihedral_types = new int [num_charge_states];
     if (rmcdata["dihedral_types"].size() != num_charge_states){
        error->all(FLERR, "Number of dihedral types doesn't match the number of charge states\n");
     }
     for (int i=0;i<num_charge_states;i++){
        dihedral_types[i] = rmcdata["dihedral_types"][i];
     }
     dihedral_list = new int*[num_dihedrals];
     for (int i=0;i<num_dihedrals;i++){
        dihedral_list[i] = new int[5];
        for (int j=0;j<4;j++){
           dihedral_list[i][j] = rmcdata["dihedral_list"][i][j];
        }
        dihedral_list[i][4] = determine_molecule(dihedral_list[i][0]);
     }
  }

  // Read in the angles that need to be altered
  if (do_angle == 1){
     num_angles = rmcdata["angle_list"].size();
     angle_types = new int [num_charge_states];
     if (rmcdata["angle_types"].size() != num_charge_states){
        error->all(FLERR, "Number of angle types doesn't match the number of charge states\n");
     }
     for (int i=0;i<num_charge_states;i++){
        angle_types[i] = rmcdata["angle_types"][i];
     }
     angle_list = new int*[num_angles];
     for (int i=0;i<num_angles;i++){
        angle_list[i] = new int[4];
        for (int j=0;j<3;j++){
           angle_list[i][j] = rmcdata["angle_list"][i][j];
        }
        angle_list[i][4] = determine_molecule(angle_list[i][0]);
     }
  }

  // Read in the bonds that need to be altered
  if (do_bond == 1){
     num_bonds = rmcdata["bond_list"].size();
     bond_types = new int [num_charge_states];
     if (rmcdata["bond_types"].size() != num_charge_states){
        error->all(FLERR, "Number of bond types doesn't match the number of charge states\n");
     }
     for (int i=0;i<num_charge_states;i++){
        bond_types[i] = rmcdata["bond_types"][i];
     }
     bond_list = new int*[num_bonds];
     for (int i=0;i<num_bonds;i++){
        bond_list[i] = new int[3];
        for (int j=0;j<2;j++){
           bond_list[i][j] = rmcdata["bond_list"][i][j];
        }
        bond_list[i][2] = determine_molecule(bond_list[i][0]);
     }
  }

   // Initialize the dynamic doping efficiency array
   dde = new double [num_charge_states];
   doping_efficiency = new double [num_charge_states];
   for (int i=0;i<num_charge_states;i++)
   {
      dde[i] = 0;
      doping_efficiency[i] = 0;
   }
   
   // Initialize acceptances/rejections
   acceptances = 0;
   rejections = 0;
   
   // get pointer for compute class, which will allow us to 
   // retrieve the potential energy
   c_pe = modify->get_compute_by_id("thermo_pe");
   

   // store the charge state for each molecule
   // also store the type of molecule (i.e is it dopant or semiconductor)

  molecule_charge_states = new double[n_molecules];
  molecule_type = new double[n_molecules];
  num_dopant_charge = new int [num_charge_states];
  num_semiconductor_charge = new int[num_charge_states];
  for (int i=0;i<num_charge_states;i++) {
      num_dopant_charge[i] = 0;
      num_semiconductor_charge[i] = 0;
  }

  
  if (restart == 0) {
     for (int i=0;i<n_molecules;i++) {
        molecule_charge_states[i] = 0.0;
        molecule_type[i] = determine_dopant_or_semiconductor(i+1);
     }
     num_dopant_charge[0] = n_dopant;
     num_semiconductor_charge[0] = n_semiconductor;
  } else {
     std::string moltype = sysname + "_type.dat";
     std::string molcharge = sysname + "_charge.dat";
     TextFileReader molecule_type_handle(moltype.c_str(), "MType");
     TextFileReader molecule_charge_handle(molcharge.c_str(), "MCharge");
     molecule_charge_handle.next_dvector(molecule_charge_states, n_molecules);
     molecule_type_handle.next_dvector(molecule_type, n_molecules);
     
     for (int i=0;i<n_molecules;i++) {
         Mol molecule = get_molecule(i+1, size_limit);
         int cstate = determine_charge_state(&molecule, molecule_type[i]);

         if (molecule_type[i] == 0)
         {
            num_semiconductor_charge[cstate] = num_semiconductor_charge[cstate] + 1;
         }
         else if (molecule_type[i] == 1)
         {
            num_dopant_charge[cstate] = num_dopant_charge[cstate] + 1;
         }
         delete_molecule(&molecule);
     }
  }
  // Combined data structures for self-energy calculations
    // Initialize them
    osc_mol_c = new double* [semiconductor_size];
    dopant_mol_c = new double* [dopant_size];
    for (int i=0;i<semiconductor_size;i++) {
      osc_mol_c[i] = new double [6];
      if (i < dopant_size) {
         dopant_mol_c[i] = new double [6];
      }
    }

  if (comm->me == 0) {
    std::string mesg =
      "###############################################\n"
      "        RMC PARTIAL CHARGE INITIALIZATION      \n"
      "###############################################\n";

    mesg += (restart == 1) ? "Continuing Run\n" : "New Run\n";

    mesg += fmt::format("RMC frequency: {}\n", periodicity);
    mesg += fmt::format("Number of RMC moves per turn: {}\n",nmoves);
    mesg += fmt::format("Size of dopant molecule: {}\n",dopant_size);
    mesg += fmt::format("Size of semiconductor molecule: {}\n",semiconductor_size);
    mesg += fmt::format("Number of dopant molecules: {}\n",n_dopant);
    mesg += fmt::format("Number of semiconductor molecules: {}\n",n_semiconductor);
    mesg += fmt::format("Temperature: {}\n",temperature);
    mesg += fmt::format("Number of charge states: {}\n", num_charge_states);
    mesg += "Charges: ";
    for (int i=0;i<num_charge_states;i++) {
      mesg += fmt::format("{} ",charges[i]);
    }
    mesg += "\nSemiconductor numbers for each charge type: ";
    for (int i=0;i<num_charge_states;i++) {
      mesg += fmt::format("{} ", num_semiconductor_charge[i]);
    }
    mesg += "\nDopant numbers for each charge type: ";
    for (int i=0;i<num_charge_states;i++) {
      mesg += fmt::format("{} ", num_dopant_charge[i]);
    }
    mesg += "\nReaction Energies: ";
    for (int i=0; i<num_charge_states;i++) {
      mesg += fmt::format("{} ",delta_g_list[i]);
    }
    mesg += "\nbetaDelta_G: ";
    for (int i=0;i<num_charge_states;i++) {
      mesg += fmt::format("{} ", beta*delta_g_list[i]);
    }
    mesg += fmt::format("\nDihedral modification?: {}", do_dihedral ? "yes" : "no");
    if (do_dihedral == 1){
      mesg += "Dihedral types: ";
      for (int i=0; i < num_charge_states;i++) {
        mesg += fmt::format("{} ", dihedral_types[i]);
      }
      mesg += fmt::format("\nNumber of dihedrals: {}\n", num_dihedrals);
    }
    mesg += fmt::format("Angle modification?: {}", do_angle ? "yes" : "no");
    if (do_angle == 1){
      mesg += "Angle types: ";
      for (int i=0;i<num_charge_states;i++){
        mesg += fmt::format("{} ", angle_types[i]);
      }
      mesg += fmt::format("\nNumber of angles: {}\n", num_angles);
    }
    mesg += fmt::format("Bond modification?: {}\n", do_bond ? "yes" : "no");
    if (do_bond == 1){
      mesg += "Bond types: ";
      for (int i=0;i<num_charge_states;i++){
        mesg += fmt::format("{} ", bond_types[i]);
      }
      mesg += fmt::format("\nNumber of bonds: {}\n", num_bonds);
    }
    mesg += fmt::format("Type threshold: {}\n",type_threshold);
    mesg += "###############################################\n";
    utils::logmesg(lmp, mesg);
  }

}

int FixRMCPartial::determine_molecule(int global_id)
{ 
  int mol_id = 0;
  int global_mol_id = 0;
  for (int i=0;i<atom->nlocal;i++) {  
     if (atom->tag[i] == global_id) { 
       mol_id = atom->molecule[i];
     }
  }
  MPI_Allreduce(&mol_id, &global_mol_id, 1, MPI_INT, MPI_SUM, world);
  return global_mol_id;
}

void FixRMCPartial::calculateMoleculeCOM(double* com, struct Mol* molecule)
{
   double local_com[3];
   double unwrap[3];
   double global_mass=0;
   double local_mass=0;
   local_com[0] = 0.0, com[0] = 0.0;
   local_com[1] = 0.0, com[1] = 0.0;
   local_com[2] = 0.0, com[2] = 0.0;
   for (int i=0;i < molecule->local_atoms;i++) {
      // unwrap so that the COM isn't complete garbage
      domain->unmap(molecule->pos[i], molecule->image[i], unwrap);
      for (int j=0;j<3;j++) {
         local_com[j] = local_com[j] + molecule->mass[i]*unwrap[j];
      }
      local_mass = local_mass + molecule->mass[i];
   }
   MPI_Allreduce(local_com, com, 3, MPI_DOUBLE, MPI_SUM, world);
   MPI_Allreduce(&local_mass, &global_mass, 1, MPI_DOUBLE, MPI_SUM, world);
   for (int j=0;j<3;j++) {
      com[j] = com[j] * (1.0/global_mass);
   }
   // map it back into the box
   domain->remap(com);
}

double FixRMCPartial::calDistance(double *pos1, double *pos2)
{
   double dr[3], period[3], dmin[3];
   double dx=0.0, dy=0.0, dz=0.0;
   double dist=0.0;
   for (int j=0;j<3;j++) {
      dr[j] = pos1[j] - pos2[j];
      period[j] = domain->boxhi[j] - domain->boxlo[j];
      dmin[j] = std::min(fabs(dr[j]), std::min(fabs(dr[j] - period[j]), fabs(dr[j] + period[j])));
      dist = dist + dmin[j]*dmin[j];
   }
   return pow(dist, 0.5);
}

double FixRMCPartial::calculateColoumbSelf(double **combined_struct, int molsize)
{
   double cpair_old=0.0;
   double cpair_new =0.0;
   double dist=0.0;
   double mol_self_energy=0.0;
   double sbcoeff;
   for (int i=0;i<molsize;i++) {
      for (int j=i+1;j<molsize;j++) {
         dist = calDistance(combined_struct[i], combined_struct[j]);
         sbcoeff = getSpecialBondCoefficient(combined_struct[i][5], combined_struct[j][5]);
         cpair_old = force->qqrd2e * (sbcoeff*combined_struct[i][3]*combined_struct[j][3])/dist;
         cpair_new = force->qqrd2e * (sbcoeff*combined_struct[i][4]*combined_struct[j][4])/dist;
         mol_self_energy = mol_self_energy + cpair_new - cpair_old;
      }
   }
   return mol_self_energy;
}

void FixRMCPartial::calculateColoumbCross(double **osc_struct, int osc_size, double **dopant_struct, int dopant_size, double *result)
{
   double cpair_old=0.0;
   double cpair_new=0.0;
   result[0] = 0.0;
   result[1] = 0.0;
   double dist=0.0;
   for (int i=0;i<osc_size;i++) {
      for (int j=0;j<dopant_size;j++) {
         dist = calDistance(osc_struct[i], dopant_struct[j]);
         cpair_old = force->qqrd2e * osc_struct[i][3] * dopant_struct[j][3]/dist;
         cpair_new = force->qqrd2e * osc_struct[i][4] * dopant_struct[j][4]/dist;
         result[0] = result[0] + cpair_old;
         result[1] = result[1] + cpair_new;
      }
   }
}

double FixRMCPartial::getSpecialBondCoefficient(int iglobal, int jglobal)
{
   double *factor, *global_factor;
   double final_factor;
   factor = new double [comm->nprocs];
   global_factor = new double[comm->nprocs];
   for (int i=0;i<comm->nprocs;i++) {
      factor[i] = 0.0;
      global_factor[i] = 0.0;
   }

   for (int i=0;i<atom->nlocal;i++) {
      if (atom->tag[i] == iglobal) {
         // Found our "first" atom, let's look at its list of specials
         for (int j=0;j<atom->nspecial[i][0];j++) {// check if 1-2
            if (atom->special[i][j] == jglobal) {
               factor[comm->me]=force->special_coul[1]-1;
            }
         }
         for (int j=atom->nspecial[i][0];j<atom->nspecial[i][1];j++) { // check if 1-3
            if (atom->special[i][j] == jglobal) {
               factor[comm->me]=force->special_coul[2]-1;
            }
         }
         for (int j=atom->nspecial[i][1];j<atom->nspecial[i][2];j++) {// check if 1-4
            if (atom->special[i][j] == jglobal) {
               factor[comm->me]=force->special_coul[3]-1;
            }
         }
         if (factor[comm->me] != force->special_coul[1]-1 && factor[comm->me] != force->special_coul[2]-1 && factor[comm->me] != force->special_coul[3]-1) {
            factor[comm->me]=-10; // These two atoms are confirmed to not be special bonds
         }
      }
   }

   MPI_Allreduce(factor, global_factor, comm->nprocs, MPI_DOUBLE, MPI_SUM, world);
   
   int count=0;
   for (int i=0;i<comm->nprocs;i++) {
      if (global_factor[i] != 0) {
         final_factor = global_factor[i];
         count = count + 1;
      }
   }
   if (count > 1) {
      error->all(FLERR, "The same atom seems to be on multiple processors, how is that possible?\n");
   }

   if (final_factor == -10) {
      final_factor = 1;
   } else {
      final_factor = final_factor+1;
   }

   delete[] factor; 
   delete[] global_factor;

   return final_factor;
}

void FixRMCPartial::bringMoleculeTogether(struct Mol* molecule, double **combined_struct, int molsize)
{
   // Set combined_struct to zero
   double **local_struct;
   double unwrap[3];
   local_struct = new double*[molsize];
   for (int i=0;i<molsize;i++) {
      local_struct[i] = new double [6];
      for (int j=0;j<6;j++) {
         local_struct[i][j] = 0.0;
      }
   }

   // First identify how much of the molecule is on each rank, and bring it together
   int *molsize_across_ranks, *molsize_combined;
   molsize_across_ranks = new int [comm->nprocs];
   molsize_combined = new int [comm->nprocs];

   memset(molsize_across_ranks, 0, comm->nprocs*sizeof(int));
   molsize_across_ranks[comm->me] = molecule->local_atoms;
   MPI_Allreduce(molsize_across_ranks, molsize_combined, comm->nprocs, MPI_INT, MPI_SUM, world);


   // Create a cumulative list of atoms per rank
   int *cumulative_molsize;
   cumulative_molsize = new int[comm->nprocs];
   cumulative_molsize[0] = 0;
   for (int i=1;i<comm->nprocs;i++) {
      cumulative_molsize[i] = cumulative_molsize[i-1] + molsize_combined[i-1];
   }

   int base=cumulative_molsize[comm->me];
   int index;

   // Now copy the position and charge data 
   for (int i=0;i<molecule->local_atoms;i++) {
      index=base+i;
      //utils::logmesg(lmp, "{}\n", index);
      for (int j=0;j<3;j++) {
         local_struct[index][j] = molecule->pos[i][j];
      }
      local_struct[index][3] = molecule->charge[i];
      local_struct[index][4] = molecule->new_charge[i];
      local_struct[index][5] = molecule->global_tag[i];
   }
   for (int i=0;i<molsize;i++) {
      MPI_Allreduce(local_struct[i], combined_struct[i], 6, MPI_DOUBLE, MPI_SUM, world);
   }
   for (int i=0;i<molsize;i++) {
      delete[] local_struct[i];
   }
   delete[] local_struct;
   delete[] molsize_across_ranks;
   delete[] molsize_combined;
   delete[] cumulative_molsize;
}

FixRMCPartial::Mol FixRMCPartial::initialize_molecule(int num_atoms_max)
{  
   Mol molecule;
   molecule.pos = new double*[num_atoms_max];
   for (int i=0;i<num_atoms_max;i++) { 
     molecule.pos[i] = new double[3];
   }
   molecule.charge = new double[num_atoms_max];
   molecule.new_charge = new double [num_atoms_max];
   molecule.type = new int[num_atoms_max];
   molecule.mass = new double[num_atoms_max];
   molecule.image = new imageint[num_atoms_max];
   molecule.global_tag = new int[num_atoms_max];
   molecule.local_tag = new int[num_atoms_max];
   molecule.local_index = new int[num_atoms_max];
   return molecule;
}

 
int FixRMCPartial::setmask()
{ 
  int mask = 0;
  mask |= FixConst::INITIAL_INTEGRATE;
  return mask;
}

void FixRMCPartial::initial_integrate(int /*vflag*/)
{
  if (perform_step == 0 || perform_step == update->ntimestep) {
     for (int move=1;move<=nmoves;move++) {
        if (comm->me == 0) utils::logmesg(lmp, "Move {} out of {}\n", move, nmoves);
        make_move();
     }
     nmcsteps = nmcsteps + nmoves;
     // Update when next to perform ReactiveMC
     perform_step = update->ntimestep + periodicity;
  }
}

double FixRMCPartial::energy_full()
{
  
  int eflag = 1;
  int vflag = 0;

  if (modify->n_pre_force) modify->pre_force(vflag);

  if (force->pair) force->pair->compute(eflag, vflag);

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) force->bond->compute(eflag, vflag);
    if (force->angle) force->angle->compute(eflag, vflag);
    if (force->dihedral) force->dihedral->compute(eflag, vflag);
    if (force->improper) force->improper->compute(eflag, vflag);
  }

  if (force->kspace) force->kspace->compute(eflag, vflag);

  if (modify->n_post_force_any) modify->post_force(vflag);
  
  double total_energy = c_pe->compute_scalar();
  update->eflag_global = update->ntimestep;

  return total_energy;
}

int FixRMCPartial::determine_dopant_or_semiconductor(int mol_id)
{
   int *indicator, *global_indicator;
   indicator = new int[comm->nprocs];
   global_indicator = new int[comm->nprocs];
   int final_indicator;
   memset(indicator, 0.0, comm->nprocs*sizeof(int));
   memset(global_indicator, 0.0, comm->nprocs*sizeof(int));
   for (int i=0;i<atom->nlocal;i++) {
      indicator[comm->me] = -1;
      if (atom->molecule[i] == mol_id) {
        if (atom->type[i] > type_threshold) {
           indicator[comm->me] = 1; // dopant, above type_threshold
           break;
        } else {
           indicator[comm->me] = 0; // semiconductor, below type_threshold
           break;
        }
      }
   }
   MPI_Allreduce(indicator, global_indicator, comm->nprocs, MPI_INT, MPI_SUM, world);
   for (int procs=0;procs<comm->nprocs;procs++) {
      if (global_indicator[procs] != -1) {
         final_indicator = global_indicator[procs];
      }
   }
   delete[] indicator;
   delete[] global_indicator;
   return final_indicator;
}

int FixRMCPartial::determine_charge_state(struct Mol* molecule, double d_or_s)
{
    int final_charge_indicator=-1;
    int *charge_indicator, *global_charge_indicator;
    charge_indicator = new int [comm->nprocs];
    global_charge_indicator = new int [comm->nprocs];
    memset(charge_indicator, 0, comm->nprocs*sizeof(int));
    memset(global_charge_indicator, 0, comm->nprocs*sizeof(int));
    if (molecule->local_atoms != 0) {
       double test_charge = molecule->charge[0];
       int local_tag = molecule->local_tag[0];
       
       if (d_or_s == 0) {
          charge_indicator[comm->me] = -1;
          for (int j=0;j<num_charge_states;j++) {
             
             if (fabs(test_charge - semiconductor_charges[j][local_tag-1]) < 1e-7) {
                charge_indicator[comm->me] = j;
             }
          }
          if (charge_indicator[comm->me] == -1) {
             error->all(FLERR, "No charge match for this semiconductor! Check the input file");
          }
       } else if (d_or_s == 1) {
          charge_indicator[comm->me] = -1;
          for (int j=0;j<num_charge_states;j++) {
            if (fabs(test_charge - dopant_charges[j][local_tag-1]) < 1e-7) {
               charge_indicator[comm->me] = j;
            }
          }
          if (charge_indicator[comm->me] == -1) {
             error->all(FLERR, "No charge match for this dopant! Check the input file");
          }
       }
    } else {
       charge_indicator[comm->me] = -1;
    }

    MPI_Allreduce(charge_indicator, global_charge_indicator, comm->nprocs, MPI_INT, MPI_SUM, world);
    for (int procs=0;procs<comm->nprocs;procs++) {
       if (global_charge_indicator[procs] != -1) {
          final_charge_indicator = global_charge_indicator[procs];
       }
    }
    delete[] charge_indicator;
    delete[] global_charge_indicator;
    return final_charge_indicator;
}

double FixRMCPartial::change_dihedral_parameters(int molecule_id, int ending_state)
{
   // go to ending_state
   
   double pre_energy = energy_full();
   
   // Circle through the relevant dihedrals and see if any need to be modified.
   for (int d=0;d<num_dihedrals;d++) {
      if (dihedral_list[d][4] == molecule_id) { // this dihedral type needs to be modified
        // First we need to find the dihedral in the main data structure
        for (int i=0;i<atom->nlocal;i++) {
           if (atom->tag[i] == dihedral_list[d][1]){
             for (int j=0;j<atom->num_dihedral[i];j++){
                if (atom->dihedral_atom1[i][j] == dihedral_list[d][0] &&
                    atom->dihedral_atom2[i][j] == dihedral_list[d][1] &&
                    atom->dihedral_atom3[i][j] == dihedral_list[d][2] &&
                    atom->dihedral_atom4[i][j] == dihedral_list[d][3]){
                  // switch to new type
                  atom->dihedral_type[i][j] = dihedral_types[ending_state];
                }
             } 
           }
        }
      }
   }

   // Update all the dihedral information in the neighbor lists
   // since that is what is used in the energy calculation
   neighbor->build_topology();

   // Now calculate the new energy
   double new_energy = energy_full();
   double energy_diff = new_energy - pre_energy;
   return energy_diff;
}

double FixRMCPartial::change_angle_parameters(int molecule_id, int ending_state)
{
   // go to ending_state
   
   double pre_energy = energy_full();
   
   // Circle through the relevant dihedrals and see if any need to be modified.
   for (int a=0;a<num_angles;a++) {
      if (angle_list[a][3] == molecule_id) { // this dihedral type needs to be modified
        // First we need to find the dihedral in the main data structure
        for (int i=0;i<atom->nlocal;i++) {
           if (atom->tag[i] == angle_list[a][1]){
             for (int j=0;j<atom->num_angle[i];j++){
                if (atom->angle_atom1[i][j] == angle_list[a][0] &&
                    atom->angle_atom2[i][j] == angle_list[a][1] &&
                    atom->angle_atom3[i][j] == angle_list[a][2]){
                  // switch to new type
                  atom->angle_type[i][j] = angle_types[ending_state];
                }
             } 
           }
        }
      }
   }


   // Update all the dihedral information in the neighbor lists
   // since that is what is used in the energy calculation
   neighbor->build_topology();

   // Now calculate the new energy
   double new_energy = energy_full();
   double energy_diff = new_energy - pre_energy;
   return energy_diff;
}

double FixRMCPartial::change_bond_parameters(int molecule_id, int ending_state)
{
   // go to ending_state
   
   double pre_energy = energy_full();
   
   // Circle through the relevant dihedrals and see if any need to be modified.
   for (int b=0;b<num_bonds;b++) {
      if (bond_list[b][2] == molecule_id) { // this dihedral type needs to be modified
        // First we need to find the dihedral in the main data structure
        for (int i=0;i<atom->nlocal;i++) {
           if (atom->tag[i] == bond_list[b][0]){
             for (int j=0;j<atom->num_bond[i];j++){
                if (atom->bond_atom[i][j] == bond_list[b][1]){
                  // switch to new type 
                  atom->bond_type[i][j] = bond_types[ending_state];
                }
             } 
           }
        }
      }
   }

   // Update all the dihedral information in the neighbor lists
   // since that is what is used in the energy calculation
   neighbor->build_topology();

   // Now calculate the new energy
   double new_energy = energy_full();
   double energy_diff = new_energy - pre_energy;
   return energy_diff;
}

FixRMCPartial::Mol FixRMCPartial::get_molecule(int mol_id, int num_atoms_max) 
{
  // Initialize memory for the molecule
  Mol molecule = initialize_molecule(num_atoms_max);

  // Get all the atoms that belong to that molecule
  // This structure will look different across MPI ranks,
  // if atoms of a single molecule are split across ranks.
  int atom_counter=0; 
  int atom_counter_global=0;
  int *lcount_across_ranks,*lmin_tag, *gmin_tag, *gcount_across_ranks;
  lcount_across_ranks = new int[comm->nprocs];
  lmin_tag = new int[comm->nprocs];
  gmin_tag = new int[comm->nprocs];
  gcount_across_ranks = new int[comm->nprocs];
  memset(lmin_tag, 0.0, comm->nprocs * sizeof(int));
  memset(lcount_across_ranks, 0.0, comm->nprocs * sizeof(int));
  int mintag = INT_MAX;
  for (int i=0;i< atom->nlocal;i++){
    if (atom->molecule[i] == (double) mol_id){ 
      molecule.charge[atom_counter] = atom->q[i];
      molecule.new_charge[atom_counter] = atom->q[i];
      molecule.type[atom_counter] = atom->type[i];
      molecule.mass[atom_counter] = atom->mass[atom->type[i]];
      molecule.image[atom_counter] = atom->image[i];
      molecule.global_tag[atom_counter] = atom->tag[i];
      molecule.local_index[atom_counter] = i;
      molecule.local_tag[atom_counter] = 0.0;
    
      if (molecule.global_tag[atom_counter] < mintag){
         mintag = molecule.global_tag[atom_counter];
      }
      for (int j=0;j<3;j++) {
        molecule.pos[atom_counter][j] = atom->x[i][j];
      }
      atom_counter = atom_counter+1;
   }
  }
  // Add up all the atoms across the different MPI ranks
  molecule.local_atoms = atom_counter;
  lmin_tag[comm->me] = mintag;
  lcount_across_ranks[comm->me] = atom_counter;
  MPI_Allreduce(lcount_across_ranks, gcount_across_ranks, comm->nprocs, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(lmin_tag, gmin_tag, comm->nprocs, MPI_INT, MPI_SUM, world);
  MPI_Barrier(world);

  // Determine overall minimum global index
  int gmintagval = INT_MAX;
  for (int i=0;i<comm->nprocs;i++){
    atom_counter_global = atom_counter_global + gcount_across_ranks[i];
    if (gmin_tag[i] != INT_MAX){

      if (gmin_tag[i] < gmintagval){
        gmintagval = gmin_tag[i];
      }
    }
  }
  // Subtract off the minimum global index so we get a 
  // new "local" index. This is not local to processor, but local to a molecule
  // This will be very helpful for charge manipulation
  for (int i=0;i<atom_counter;i++){
    molecule.local_tag[i] = molecule.global_tag[i] - gmintagval + 1;
  }

  delete[] lcount_across_ranks;
  delete[] lmin_tag;
  delete[] gmin_tag;
  delete[] gcount_across_ranks;
  return molecule;
}

void FixRMCPartial::modify_charge(struct Mol *molecule, double *charge_list)
{
  
  for (int i=0;i<molecule->local_atoms;i++){
     for (int j=0;j<atom->nlocal;j++){
        if (atom->tag[j] == molecule->global_tag[i]){
           atom->q[j] = charge_list[molecule->local_tag[i]-1];
           molecule->new_charge[i] = charge_list[molecule->local_tag[i]-1];
        }
     }
  }
}

void FixRMCPartial::restore_charge(struct Mol *molecule)
{
  for (int i=0;i<molecule->local_atoms;i++){
     for (int j=0;j<atom->nlocal;j++){
        if (atom->tag[j] == molecule->global_tag[i]){
           atom->q[j] = molecule->charge[i];
           molecule->new_charge[i] = molecule->charge[i];
        }
     }
  }

}


void FixRMCPartial::make_move()
{
    double reaction_energy = 0.0;
    double prefactor_num = 0.0;
    double prefactor_den = 0.0;
    double prefactor = 0.0;
    double transition_probability = 0.0;
    double edihedral = 0.0;
    double ebond = 0.0;
    double eangle = 0.0;

    //Clear the osc dopant combined structures
   // this is probably unnecessary and adds more time, but I'm doing it for peace of mind

    for (int i=0;i<semiconductor_size;i++){
      for (int j=0;j<6;j++){
         osc_mol_c[i][j] = 0.0;
      }
    }
    for (int i=0;i<dopant_size;i++){
      for (int j=0;j<6;j++){
         dopant_mol_c[i][j] = 0.0;
      }
    }

    // Centre of mass pointers
    double osc_com[3], dopant_com[3];
    double com_diff;

    // Self energy terms
    double e_osc_diff, e_dopant_diff;


    // Calculate the energy before we do any mischief
    double starting_energy = energy_full();

    if (comm->me == 0) utils::logmesg(lmp, "The starting energy is {}\n", starting_energy);

    int rand_semi, rand_dope, rand_charge;
    int indicator=-1;
    int c_indicator=-1;


    // Find a semiconductor - this is not good of course - you're not sampling properly
    // What you need to do is pick a random CHARGE, and then choose a molecule with that charge
    if (comm->me == 0) {
      while (c_indicator != 0){
         rand_charge = type_dist(rng_type_source);
         if (num_semiconductor_charge[rand_charge] != 0 && num_dopant_charge[rand_charge] != 0){
            c_indicator = 0;
         }
      }
      utils::logmesg(lmp, "We have chosen charge state {}\n", charges[rand_charge]);
      int tries=0;
      while (indicator != 0) {
         rand_semi = atom_dist(rng_atom);
         if (molecule_type[rand_semi-1] == 0 && fabs(molecule_charge_states[rand_semi-1] - charges[rand_charge]) < 0.001){
            indicator = 0;
         }
         tries=tries+1;
         if (tries > 10000){
            error->all(FLERR, "Too many attemps at finding a semiconductor, something is going wrong!");
         }
      }
   }
    MPI_Bcast(&rand_semi, 1, MPI_INT, 0, world);
    
    // Retrieve information on the chosen molecule and populate it in this struct
    Mol semiconductor = get_molecule(rand_semi, size_limit);

    int charge_state = determine_charge_state(&semiconductor, 0);
    semiconductor.charge_state = charge_state;

    
    indicator=-1;
    int tries=0;
    
    if (comm->me == 0) {
      while (indicator != 1){
          rand_dope = atom_dist(rng_atom);
          if (molecule_type[rand_dope-1] == 1) {
             if (fabs(molecule_charge_states[rand_dope-1] + molecule_charge_states[rand_semi-1]) < 0.001){
               utils::logmesg(lmp, "FOUND EQUALITY\n"
                              "Semiconductor charge, {} = Dopant charge {}\n",
                              molecule_charge_states[rand_semi-1],
                              molecule_charge_states[rand_dope-1]);
               indicator = 1;
             }
          }
          tries=tries+1;
          if (tries > 10000){
            error->all(FLERR, "Max tries reached to find dopant of same charge state");
          }
       }
    }
    MPI_Bcast(&rand_dope, 1, MPI_INT, 0, world);

    
    //Retrieve information on the chosen dopant and populate it in this struct
    Mol dopant = get_molecule(rand_dope, size_limit);
    charge_state = determine_charge_state(&dopant, 1);
    dopant.charge_state = charge_state;
  
    // Verify charge states are the same
    if (dopant.charge_state != semiconductor.charge_state){
       utils::logmesg(lmp, "Dopant charge state is {} and semiconductor charge state is {}\n",
                      dopant.charge_state, semiconductor.charge_state);
       error->all(FLERR, Error::NOLASTLINE, "Charge states don't match!");
    }

    // Identify a destination charge state, picked randomly
    // of course it has to be different from the starting state
    int destination_charge_state = dopant.charge_state;

    if (comm->me == 0){
      while (destination_charge_state == dopant.charge_state){
         destination_charge_state = type_dist(rng_type_destination);
      }
    }

    MPI_Bcast(&destination_charge_state, 1, MPI_INT, 0, world);
    
    MPI_Bcast(&destination_charge_state, 1, MPI_INT, 0, world);

    if (comm->me == 0)
      utils::logmesg(lmp, "Going from charge state {} to {}\n",
                     charges[semiconductor.charge_state], charges[destination_charge_state]);

    // Calculate centre of mass
    calculateMoleculeCOM(osc_com, &semiconductor);
    calculateMoleculeCOM(dopant_com, &dopant);
    
    com_diff = calDistance(osc_com, dopant_com);

    // Change the dihedral/angle/bond parameters to destination type 
    // and capture the dihedral/angle/bond energy
    if (do_dihedral == 1){
      edihedral = change_dihedral_parameters(rand_semi, destination_charge_state);
    }
    if (do_angle == 1){
      eangle = change_angle_parameters(rand_semi, destination_charge_state);
    }
    if (do_bond == 1){
      ebond = change_bond_parameters(rand_semi, destination_charge_state);
    }

    // Modify charge to new type
    modify_charge(&semiconductor, semiconductor_charges[destination_charge_state]);
    modify_charge(&dopant, dopant_charges[destination_charge_state]);



   // Capture the molecule self-energies difference
    bringMoleculeTogether(&semiconductor, osc_mol_c, semiconductor_size);
    bringMoleculeTogether(&dopant, dopant_mol_c, dopant_size);
    
    e_osc_diff = calculateColoumbSelf(osc_mol_c, semiconductor_size);
    e_dopant_diff = calculateColoumbSelf(dopant_mol_c, dopant_size);

    // Capture the pair coloumb interaction
    double interaction_energy[2];
    calculateColoumbCross(osc_mol_c, semiconductor_size, dopant_mol_c, dopant_size, interaction_energy);


    if (comm->me == 0){
      std::string mesg;
      if (do_dihedral == 1){
        mesg += fmt::format("The dihedral energy change is {}\n", edihedral);
      }
      if (do_angle == 1){
        mesg+= fmt::format("The angle energy change is {}\n", eangle);
      }
      if (do_bond == 1){
        mesg += fmt::format("The bond energy change is {}\n", ebond);
      }
      mesg += fmt::format("The semiconductor self energy difference is {}\n", e_osc_diff);
      mesg += fmt::format("The dopant self-energy difference is {}\n", e_dopant_diff);
      mesg += fmt::format("The cross-energy before and after doping are {} {}\n", interaction_energy[0], interaction_energy[1]);
      mesg += fmt::format("The COM distance is {}\n", com_diff);
      utils::logmesg(lmp, mesg);
    }


    // Reinitialize Ewald
    force->kspace->init();
    reaction_energy = delta_g_list[destination_charge_state] - delta_g_list[semiconductor.charge_state];
    prefactor_num = (double)(num_semiconductor_charge[semiconductor.charge_state]*num_dopant_charge[dopant.charge_state]);
    prefactor_den = (double)(num_semiconductor_charge[destination_charge_state]+1.0)*(num_dopant_charge[destination_charge_state]+1.0);
    prefactor = prefactor_num/prefactor_den;

    // Recalculate energy after these changes
    // subtract the dihedral energy, that has to be put back in with more thought later
    double new_energy = energy_full();
    double energy_diff = new_energy - starting_energy + reaction_energy - e_osc_diff - e_dopant_diff;
    
    if (do_dihedral==1){
      energy_diff = energy_diff - edihedral;
    }
    if (do_angle==1){
      energy_diff = energy_diff - eangle;
    }
    if (do_bond==1){
      energy_diff = energy_diff - ebond;
    }
    double total_electrostatic_diff = energy_diff - reaction_energy;

    if (comm->me == 0)
       utils::logmesg(lmp, "Total Electrostatic Energy Difference: {}\n", total_electrostatic_diff);
    // Calculate acceptance probability
    transition_probability = prefactor*exp(-beta*energy_diff);

    if (comm->me == 0)
       utils::logmesg(lmp, "The new energy is {} and the difference is {}\n",  new_energy, energy_diff);
    

    // Generate random number in one rank
    double rand_number=0;
    if (comm->me == 0){
      rand_number=acc_dist(rng_acc);
    }
    MPI_Bcast(&rand_number, 1, MPI_DOUBLE, 0, world);

    // Determine whether to accept or reject
    if (transition_probability > rand_number){
       // Move accepted
       acceptances = acceptances+1;
       if (comm->me == 0) utils::logmesg(lmp, "MOVE ACCEPTED\n");

       num_semiconductor_charge[destination_charge_state] = num_semiconductor_charge[destination_charge_state]+1;
       num_dopant_charge[destination_charge_state] = num_dopant_charge[destination_charge_state]+1;
       num_semiconductor_charge[semiconductor.charge_state] = num_semiconductor_charge[semiconductor.charge_state]-1;
       num_dopant_charge[dopant.charge_state] = num_dopant_charge[dopant.charge_state]-1;
       molecule_charge_states[rand_semi-1] = charges[destination_charge_state];
       molecule_charge_states[rand_dope-1] = -charges[destination_charge_state];
    } else {
       // Move rejected
       rejections = rejections+1;

       // Revert dihedral coefficients
       if (do_dihedral == 1){
          edihedral = change_dihedral_parameters(rand_semi, semiconductor.charge_state);
       }
       if (do_angle == 1){
         eangle = change_angle_parameters(rand_semi, semiconductor.charge_state);
       }
       if (do_bond == 1){
         ebond = change_bond_parameters(rand_semi, semiconductor.charge_state);
       }

       if (comm->me == 0) utils::logmesg(lmp, "MOVE REJECTED\n");

       // Restore the charges to what they were
       restore_charge(&semiconductor);
       restore_charge(&dopant);
       force->kspace->init();
    }

    // The step is done, so we can delete the molecules
    delete_molecule(&semiconductor);
    delete_molecule(&dopant);

    // Calculate dynamic doping efficiency
    for (int i=0;i<num_charge_states;i++){
         dde[i] = (double) num_dopant_charge[i]/(double)n_dopant;
    }

     if (comm->me == 0){
       std::string mesg = "Dynamic Doping Efficiency: ";
        for (int i=0;i<num_charge_states;i++){
          mesg += fmt::format("{} ", dde[i]);
        }
        mesg += "\n";
        utils::logmesg(lmp, mesg);
     }

}

void FixRMCPartial::delete_molecule(struct Mol *molecule)
{
    delete[] molecule->type;
    delete[] molecule->charge;
    delete[] molecule->mass;
    delete[] molecule->global_tag;
    delete[] molecule->local_tag;
    delete[] molecule->image;
    for (int i=0;i<molecule->local_atoms;i++){
       delete[] molecule->pos[i];
    } 
    delete[] molecule->pos;
    delete[] molecule->local_index;
    delete[] molecule->new_charge;
}

void FixRMCPartial::post_mortem()
{
   acceptance_rate = (double)acceptances/(double)nmcsteps;
   double total_charged_dopants=0;
   double total_charged_semiconductors=0;
   for (int i=0;i<num_charge_states;i++){
      doping_efficiency[i] = (double)num_dopant_charge[i]/(double)n_dopant;
      total_charged_dopants = total_charged_dopants+num_dopant_charge[i];
      total_charged_semiconductors = total_charged_semiconductors+num_semiconductor_charge[i];
   }
   double total_doping_efficiency = (double)total_charged_dopants/(double)n_dopant;

   // Write out COM for all molecules, for RDF/morphology or other CG uses
   std::string comname = sysname + "_com.xyz";
   if (comm->me == 0){
      FILE *fcom = fopen(comname.c_str(), "w");
      fprintf(fcom, "%d\n",n_molecules);
      for (int j=0;j<3;j++){
         fprintf(fcom, "%f %f\n", domain->boxlo[j], domain->boxhi[j]);
      }
      fprintf(fcom, "LAMMPS_RMC_output\n");
      fclose(fcom);
   }

   double *com;
   com = new double[3];
   for (int i=0;i<n_molecules;i++){
      Mol molecule = get_molecule(i+1,size_limit);
      calculateMoleculeCOM(com, &molecule);
      if (comm->me == 0){
         FILE *fcom = fopen(comname.c_str(), "a");
         fprintf(fcom, "%d %f %f %f\n", i+1,com[0],com[1],com[2]);
         fclose(fcom);
      }
      delete_molecule(&molecule);  
   }
 
   if (comm->me == 0){
      // Write the charge of each molecule to a file for morphology analysis
      std::string molcharge = sysname + "_charge.dat";
      std::string moltype = sysname + "_type.dat";
      FILE *fp = fopen(molcharge.c_str(), "w");
      for (int i=0;i<n_molecules;i++){
         fprintf(fp, "%f\n", molecule_charge_states[i]);
      }
      fclose(fp);

      // Write the type of each molecule to a file for morphology analysis
      FILE *ftype = fopen(moltype.c_str(), "w");
      for (int i=0;i<n_molecules;i++){
         fprintf(ftype, "%f\n", molecule_type[i]);
      }
      fclose(ftype);

      std::string mesg =
        "###############################################\n"
        "              RMC OUTPUT SUMMARY               \n"
        "###############################################\n";
      mesg += fmt::format("Number of RMC moves: {}\n", nmcsteps);
      mesg += "Final charged dopant number: ";
      for (int i=0;i<num_charge_states;i++){
        mesg += fmt::format("{} ", num_dopant_charge[i]);
      }
      mesg += "\nFinal charged semiconductor number: ";
      for (int i=0;i<num_charge_states;i++){
        mesg += fmt::format("{} ", num_semiconductor_charge[i]);
      }
      mesg += fmt::format("\nFinal charged dopant number: {}\n", total_charged_dopants);
      mesg += fmt::format("Final charged semiconductor number: {}\n", total_charged_semiconductors);
      mesg += "Doping efficiency: ";
      for (int i=0;i<num_charge_states;i++){
        mesg += fmt::format("{} ", doping_efficiency[i]);
      }
      mesg += fmt::format("\nFinal acceptances: {}\n", acceptances);
      mesg += fmt::format("Final rejections: {}\n", rejections);
      mesg += fmt::format("Acceptance rate: {}\n",acceptance_rate);
      mesg += "###############################################\n";
      utils::logmesg(lmp, mesg);
   }
}


FixRMCPartial::~FixRMCPartial()
{
   post_mortem();
   delete[] dde;
   delete[] delta_g_list;
   delete[] charges;
   delete[] doping_efficiency; 
   delete[] num_dopant_charge;
   delete[] num_semiconductor_charge;
   delete[] molecule_charge_states;
   delete[] molecule_type;

    // Delete the combined structs
    for (int i=0;i<semiconductor_size;i++){
      delete[] osc_mol_c[i];
    }
    for (int i=0;i<dopant_size;i++){
      delete[] dopant_mol_c[i];
    }

    // delete the charge arrays
    for (int i=0;i<num_charge_states;i++){
      delete[] semiconductor_charges[i];
      delete[] dopant_charges[i];
    }
    if (do_dihedral == 1){
      delete[] dihedral_types;
      for (int i=0;i<num_dihedrals;i++){
        delete[] dihedral_list[i];
      }
      delete[] dihedral_list;
    }
    if (do_angle == 1){
      delete[] angle_types;
      for (int i=0;i<num_angles;i++){
         delete[] angle_list[i];
      }
      delete[] angle_list;
    }
    if (do_bond == 1){
      delete[] bond_types;
      for (int i=0;i<num_bonds;i++){
         delete[] bond_list[i];
      }
      delete[] bond_list;
    }
 
    delete[] semiconductor_charges;
    delete[] dopant_charges;
    delete[] osc_mol_c;
    delete[] dopant_mol_c;
}
