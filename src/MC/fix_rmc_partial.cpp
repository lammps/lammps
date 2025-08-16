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
#include "text_file_reader.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixRMCPartial::FixRMCPartial(LAMMPS *lmp, int narg, char **arg) : rng_acc(device_acc()), rng_atom(device_atom()), 
rng_type_source(device_type_source()), rng_type_destination(device_type_destination()), Fix(lmp, narg, arg)
{
  perform_step = 0;

  // After how many MD steps should RMC be performed 
  periodicity = utils::inumeric(FLERR, arg[3], false, lmp);
  
  // Second argument :- number of MC moves to perform per "turn"
  nmoves = utils::inumeric(FLERR, arg[4], false, lmp);

  nmcsteps=0;
 
  // Third argument :- number of atoms in dopant molecule
  dopant_size= utils::inumeric(FLERR, arg[5], false, lmp);
  
  // Fourth argument :- number of atoms in semiconductor molecule
  semiconductor_size = utils::inumeric(FLERR, arg[6], false, lmp);

  // Determine the max of the sizes, to use for memory allocation
  size_limit = std::max(dopant_size, semiconductor_size);
  
  // Fifth argument :- Number of dopant molecules
  n_dopant = utils::inumeric(FLERR, arg[7], false, lmp);
  
  // Sixth argument :- Number of semiconductor molecules
  n_semiconductor = utils::inumeric(FLERR, arg[8], false, lmp);
  
  // Defining the sum of them as a separate variable for convienience
  n_molecules = n_dopant+n_semiconductor;

  // Seventh argument :- Temperature
  temperature = utils::numeric(FLERR, arg[9], false, lmp);
  
  // Calculate beta
  beta = 1.0/(force->boltz * temperature);
  
  // Ninth argument :- Type threshold (any atom type <= this is a semiconductor , > is a dopant)
  // Critical variable to determine if a molecule is an semiconductor or a dopant
  type_threshold = utils::inumeric(FLERR, arg[10], false, lmp);

  // Tenth argument (only for partial) :- How many charge states to study?
  num_charge_states = utils::inumeric(FLERR, arg[11], false, lmp);

  // Eleventh argument (is this fresh or are we restarting)
  restart = utils::inumeric(FLERR, arg[12], false, lmp);

  // Initializing random number generators
  atom_dist = std::uniform_int_distribution<> (1, n_molecules);
  type_dist = std::uniform_int_distribution<> (0,num_charge_states-1);
  acc_dist = std::uniform_real_distribution<> (0,1);

  // Separation between charge states
  deltaQ = 1.0/(num_charge_states-1.0);
  
  // List of charges
  charges = new double [num_charge_states];
  charges[0] = 0;
  for (int i=0;i<num_charge_states-1;i++)
  {
     charges[i+1] = charges[i] + deltaQ;
  }
  
  // Data structure to store semiconductor and dopant charges
  semiconductor_charges = new double*[num_charge_states];
  dopant_charges = new double*[num_charge_states];
  for (int i=0;i<num_charge_states;i++)
  {
     semiconductor_charges[i] = new double [semiconductor_size];
     dopant_charges[i] = new double [dopant_size];
  }

  
  // Populate the neutral and integer charge state from files
  TextFileReader semi_integer("semiconductor_charge.dat", "SCharges");
  TextFileReader semi_neutral("semiconductor_neutral.dat", "SNeutral");
  TextFileReader dope_integer("dopant_charge.dat", "DCharges");
  TextFileReader dope_neutral("dopant_neutral.dat", "DNeutral");
  semi_neutral.next_dvector(semiconductor_charges[0], semiconductor_size);
  semi_integer.next_dvector(semiconductor_charges[num_charge_states-1], semiconductor_size);
  dope_neutral.next_dvector(dopant_charges[0], dopant_size);
  dope_integer.next_dvector(dopant_charges[num_charge_states-1], dopant_size);
  
  // Populate the intermediate charge values for each semiconductor atom using linear interpolation
  double temp_delta_semiconductor, temp_delta_dopant;
  
  for (int i=0;i<semiconductor_size;i++)
  {
     temp_delta_semiconductor = (semiconductor_charges[num_charge_states-1][i] - semiconductor_charges[0][i])/(num_charge_states-1);
     for (int j=0;j<num_charge_states-1;j++)
     {
        semiconductor_charges[j+1][i] = semiconductor_charges[j][i] + temp_delta_semiconductor;
     }
  }
  
  for (int i=0;i<dopant_size;i++)
  {
     temp_delta_dopant = (dopant_charges[num_charge_states-1][i] - dopant_charges[0][i])/(num_charge_states-1);
     for (int j=0;j<num_charge_states-1;j++)
     {
        dopant_charges[j+1][i] = dopant_charges[j][i] + temp_delta_dopant;
     }
  }
 
 
  // Create and populate the delta_g list for various charge states 
  delta_g_list = new double [num_charge_states];
  TextFileReader barrier_data("barrier_list.dat", "barriers");
  char *barriers = barrier_data.next_line(num_charge_states);
  ValueTokenizer vt_barrier(barriers);
  for (int i=0;i<num_charge_states;i++)
  {
     delta_g_list[i] =  vt_barrier.next_double();
  }


  // Read in the dihedrals that need to be altered
  TextFileReader dihedral_data("dihedral_list.dat", "dihedrals");

  // get number of dihedrals
  char *ndihedrals = dihedral_data.next_line(1);
  ValueTokenizer vt(ndihedrals, "\n");
  num_dihedrals = vt.next_int();

  // get dihedrals types
  dihedral_types = new int [num_charge_states];
  char *ndtypes = dihedral_data.next_line(num_charge_states);
  ValueTokenizer vt1(ndtypes);
  for (int i=0;i<num_charge_states;i++)
  {
     dihedral_types[i] = vt1.next_int();  
  }

  // Get the dihedral atom indices 
  dihedral_list = new int*[num_dihedrals];
  
  for (int i=0;i<num_dihedrals;i++)
  {
     dihedral_list[i] = new int[5];
     char *dihedral_line = dihedral_data.next_line(4);
     ValueTokenizer vt(dihedral_line);
     for (int j=0;j<4;j++)
     {
        dihedral_list[i][j] = vt.next_int();
     }
     dihedral_list[i][4] = determine_molecule(dihedral_list[i][0]);
  }

  
  // Read in the angles that need to be altered
  TextFileReader angle_data("angle_list.dat", "angles");

  // get number of angles
  char *nangles = angle_data.next_line(1);
  ValueTokenizer vt_angle(nangles, "\n");
  num_angles = vt_angle.next_int();

  // get angle types
  angle_types = new int [num_charge_states];
  char *natypes = angle_data.next_line(num_charge_states);
  ValueTokenizer vt1_angle(natypes);
  for (int i=0;i<num_charge_states;i++)
  {
     angle_types[i] = vt1_angle.next_int();  
  }

  // Get the angles atom indices
  //fmt::print(screen, "{} {}\n", "The number of dihedrals to process is ", num_dihedrals);  
  angle_list = new int*[num_angles];
  
  for (int i=0;i<num_angles;i++)
  {
     angle_list[i] = new int[4];
     char *angle_line = angle_data.next_line(3);
     ValueTokenizer vt_angle(angle_line);
     for (int j=0;j<3;j++)
     {
        angle_list[i][j] = vt_angle.next_int();
     }
     angle_list[i][3] = determine_molecule(angle_list[i][0]);
  }

  
  // Read in the bonds that need to be altered
  TextFileReader bond_data("bond_list.dat", "bonds");

  // get number of bonds
  char *nbonds = bond_data.next_line(1);
  ValueTokenizer vt_bond(nbonds, "\n");
  num_bonds = vt_bond.next_int();

  // get bonds types
  bond_types = new int [num_charge_states];
  char *nbtypes = bond_data.next_line(num_charge_states);
  ValueTokenizer vt1_bond(nbtypes);
  for (int i=0;i<num_charge_states;i++)
  {
     bond_types[i] = vt1_bond.next_int();  
  }

  // Get the dihedral atom indices
  //fmt::print(screen, "{} {}\n", "The number of dihedrals to process is ", num_dihedrals);  
  bond_list = new int*[num_bonds];
  
  for (int i=0;i<num_bonds;i++)
  {
     bond_list[i] = new int[3];
     char *bond_line = bond_data.next_line(2);
     ValueTokenizer vt_bond(bond_line);
     for (int j=0;j<2;j++)
     {
        bond_list[i][j] = vt_bond.next_int();
     }
     bond_list[i][2] = determine_molecule(bond_list[i][0]);
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
  for (int i=0;i<num_charge_states;i++)
  {
      num_dopant_charge[i] = 0;
      num_semiconductor_charge[i] = 0;
  }

  
  if (restart == 0)
  {
     for (int i=0;i<n_molecules;i++)
     {
        molecule_charge_states[i] = 0.0;
        molecule_type[i] = determine_dopant_or_semiconductor(i+1);
     }
     num_dopant_charge[0] = n_dopant;
     num_semiconductor_charge[0] = n_semiconductor;
  }
  else
  {
     TextFileReader molecule_type_handle("molecule_type.dat", "MType");
     TextFileReader molecule_charge_handle("molecule_charge.dat", "MCharge");
     molecule_charge_handle.next_dvector(molecule_charge_states, n_molecules);
     molecule_type_handle.next_dvector(molecule_type, n_molecules);
     
     for (int i=0;i<n_molecules;i++)
     {
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
    for (int i=0;i<semiconductor_size;i++)
    {
      osc_mol_c[i] = new double [6];
      if (i < dopant_size)
      {
         dopant_mol_c[i] = new double [6];
      }
    }


  if (comm->me == 0)
  {   
      fmt::print(screen, "{}\n", "###############################################");
      fmt::print(screen, "{}\n", "        RMC PARTIAL CHARGE INITIALIZATION      ");
      fmt::print(screen, "{}\n", "###############################################");
      if (restart == 1)
      {
         fmt::print(screen, "{}\n", "Continuing Run");
      }
      else
      {
         fmt::print(screen, "{}\n", "New Run");
      }
      fmt::print(screen,"{} {}\n","RMC frequency: ", periodicity);
      fmt::print(screen,"{} {}\n","Number of RMC moves per turn: ",nmoves);
      fmt::print(screen,"{} {}\n","Size of dopant molecule: ",dopant_size);
      fmt::print(screen,"{} {}\n","Size of semiconductor molecule: ",semiconductor_size);
      fmt::print(screen,"{} {}\n","Number of dopant molecules: ",n_dopant);
      fmt::print(screen,"{} {}\n","Number of semiconductor molecules: ",n_semiconductor);
      fmt::print(screen,"{} {}\n","Temperature: ",temperature);
      fmt::print(screen,"{} {}\n","Number of charge states: ", num_charge_states);
      fmt::print(screen,"{} ", "Charges: ");
      for (int i=0;i<num_charge_states;i++)
      {
         fmt::print(screen,"{} ",charges[i]);
      }
      fmt::print(screen, "\n{} ", "Semiconductor numbers for each charge type: ");
      for (int i=0;i<num_charge_states;i++)
      {
         fmt::print(screen, "{} ", num_semiconductor_charge[i]);
      }
      fmt::print(screen, "\n{} ", "Dopant numbers for each charge type: ");
      for (int i=0;i<num_charge_states;i++)
      {
         fmt::print(screen, "{} ", num_dopant_charge[i]);
      }
      fmt::print(screen,"\n{} ", "Reaction Energies: ");
      for (int i=0; i<num_charge_states;i++)
      {
         fmt::print(screen,"{} ",delta_g_list[i]);
      }
      fmt::print(screen,"\n{} ", "betaDelta_G: ");
      for (int i=0;i<num_charge_states;i++)
      {
         fmt::print(screen,"{} ", beta*delta_g_list[i]);
      }
      fmt::print(screen,"\n{} ", "Dihedral types: ");
      for (int i=0;i<num_charge_states;i++)
      {
         fmt::print(screen,"{} ", dihedral_types[i]);
      }
      fmt::print(screen,"\n{} {}\n","Type threshold: ",type_threshold);
      fmt::print(screen, "{}\n", "###############################################");
  }
     
}

int FixRMCPartial::determine_molecule(int global_id)
{ 
  int mol_id = 0;
  int global_mol_id = 0;
  for (int i=0;i<atom->nlocal;i++)
  {  
     if (atom->tag[i] == global_id)
     { 
       mol_id = atom->molecule[i];
     }
  }
  MPI_Allreduce(&mol_id, &global_mol_id, 1, MPI_INT, MPI_SUM, world);
  return global_mol_id;
}

void FixRMCPartial::calculateMoleculeCOM(double* com, struct Mol* molecule)
{
   double *local_com;
   double unwrap[3];
   double global_mass=0;
   double local_mass=0;
   local_com = new double [3];
   local_com[0] = 0.0, com[0] = 0.0;
   local_com[1] = 0.0, com[1] = 0.0;
   local_com[2] = 0.0, com[2] = 0.0;
   for (int i=0;i < molecule->local_atoms;i++)
   {
      // unwrap so that the COM isn't complete garbage
      domain->unmap(molecule->pos[i], molecule->image[i], unwrap);
      for (int j=0;j<3;j++)
      {
         local_com[j] = local_com[j] + molecule->mass[i]*unwrap[j];
      }
      local_mass = local_mass + molecule->mass[i];
   }
   MPI_Allreduce(local_com, com, 3, MPI_DOUBLE, MPI_SUM, world);
   MPI_Allreduce(&local_mass, &global_mass, 1, MPI_DOUBLE, MPI_SUM, world);
   for (int j=0;j<3;j++)
   {
      com[j] = com[j] * (1.0/global_mass);
   }
   // map it back into the box
   domain->remap(com);
   //fmt::print(screen, "{} {} {}\n", com[0], com[1], com[2]);
}

double FixRMCPartial::calDistance(double *pos1, double *pos2)
{
   double dr[3], period[3], dmin[3];
   double dx=0.0, dy=0.0, dz=0.0;
   double dist=0.0;
   for (int j=0;j<3;j++)
   {
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
   for (int i=0;i<molsize;i++)
   {
      for (int j=i+1;j<molsize;j++)
      {
         
         dist = calDistance(combined_struct[i], combined_struct[j]);
         sbcoeff = getSpecialBondCoefficient(combined_struct[i][5], combined_struct[j][5]);
         cpair_old = force->qqrd2e * (sbcoeff*combined_struct[i][3]*combined_struct[j][3])/dist;
         cpair_new = force->qqrd2e * (sbcoeff*combined_struct[i][4]*combined_struct[j][4])/dist;
         if (comm->me == 0)
         { 
         }
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
   for (int i=0;i<osc_size;i++)
   {
      for (int j=0;j<dopant_size;j++)
      {
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
   for (int i=0;i<comm->nprocs;i++)
   {
      factor[i] = 0.0;
      global_factor[i] = 0.0;
   }

   for (int i=0;i<atom->nlocal;i++)
   {
      if (atom->tag[i] == iglobal)
      {
         // Found our "first" atom, let's look at its list of specials
         for (int j=0;j<atom->nspecial[i][0];j++) // check if 1-2
         {
            if (atom->special[i][j] == jglobal)
            {
               factor[comm->me]=force->special_coul[1]-1;
            }
         }
         for (int j=atom->nspecial[i][0];j<atom->nspecial[i][1];j++) // check if 1-3
         {
            if (atom->special[i][j] == jglobal)
            {
               factor[comm->me]=force->special_coul[2]-1;
            }
         }
         for (int j=atom->nspecial[i][1];j<atom->nspecial[i][2];j++) // check if 1-4
         {
            if (atom->special[i][j] == jglobal)
            {
               factor[comm->me]=force->special_coul[3]-1;
            }
         }
         if (factor[comm->me] != force->special_coul[1]-1 && factor[comm->me] != force->special_coul[2]-1 && factor[comm->me] != force->special_coul[3]-1)
         {
            factor[comm->me]=-10; // These two atoms are confirmed to not be special bonds
         }
      }
   }

   MPI_Allreduce(factor, global_factor, comm->nprocs, MPI_DOUBLE, MPI_SUM, world);
   
   int count=0;
   for (int i=0;i<comm->nprocs;i++)
   {
      if (global_factor[i] != 0)
      {
         final_factor = global_factor[i];
         count = count + 1;
      }
   }
   if (count > 1)
   {
      error->all(FLERR, "The same atom seems to be on multiple processors, how is that possible?\n");
   }

   if (final_factor == -10)
   {
      final_factor = 1;
   }
   else
   {
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
   for (int i=0;i<molsize;i++)
   {
      local_struct[i] = new double [6];
      for (int j=0;j<6;j++)
      {
         local_struct[i][j] = 0.0;
      }
   }

   // First identify how much of the molecule is on each rank, and bring it together
   int molsize_across_ranks[comm->nprocs];
   int molsize_combined[comm->nprocs];

   memset(molsize_across_ranks, 0, comm->nprocs*sizeof(int));
   molsize_across_ranks[comm->me] = molecule->local_atoms;
   MPI_Allreduce(molsize_across_ranks, molsize_combined, comm->nprocs, MPI_INT, MPI_SUM, world);


   // Create a cumulative list of atoms per rank
   int cumulative_molsize[comm->nprocs];
   cumulative_molsize[0] = 0;
   for (int i=1;i<comm->nprocs;i++)
   {
      cumulative_molsize[i] = cumulative_molsize[i-1] + molsize_combined[i-1];
   }

   int base=cumulative_molsize[comm->me];
   int index;

   // Now copy the position and charge data 
   for (int i=0;i<molecule->local_atoms;i++)
   {
      index=base+i;
      //fmt::print(screen, "{}\n", index);
      for (int j=0;j<3;j++)
      {
         local_struct[index][j] = molecule->pos[i][j];
      }
      local_struct[index][3] = molecule->charge[i];
      local_struct[index][4] = molecule->new_charge[i];
      local_struct[index][5] = molecule->global_tag[i];
   }
   for (int i=0;i<molsize;i++)
   {
      MPI_Allreduce(local_struct[i], combined_struct[i], 6, MPI_DOUBLE, MPI_SUM, world);
   }
   for (int i=0;i<molsize;i++)
   {
      delete[] local_struct[i];
   }
   delete[] local_struct;
}

FixRMCPartial::Mol FixRMCPartial::initialize_molecule(int num_atoms_max)
{  
   Mol molecule;
   molecule.pos = new double*[num_atoms_max];
   for (int i=0;i<num_atoms_max;i++)
   { 
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
  if (perform_step == 0 || perform_step == update->ntimestep)
  {  
     for (int move=1;move<=nmoves;move++)
     {  
        if (comm->me == 0)
        {  
           fmt::print(screen, "{} {} {} {}\n", "Move ", move, " out of ", nmoves);
        }
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
   int indicator[comm->nprocs];
   int global_indicator[comm->nprocs];
   int final_indicator;
   memset(indicator, 0.0, comm->nprocs*sizeof(int));
   memset(global_indicator, 0.0, comm->nprocs*sizeof(int));
   for (int i=0;i<atom->nlocal;i++)
   {
      indicator[comm->me] = -1;
      if (atom->molecule[i] == mol_id)
      {
        if (atom->type[i] > type_threshold)
        {
           indicator[comm->me] = 1; // dopant, above type_threshold
           break;
        }
        else
        {
           indicator[comm->me] = 0; // semiconductor, below type_threshold
           break;
        }
      }
   }
   MPI_Allreduce(indicator, global_indicator, comm->nprocs, MPI_INT, MPI_SUM, world);
   for (int procs=0;procs<comm->nprocs;procs++)
   {
      if (global_indicator[procs] != -1)
      {
         final_indicator = global_indicator[procs];
      }
   }
   return final_indicator;
}

int FixRMCPartial::determine_charge_state(struct Mol* molecule, double d_or_s)
{
    int final_charge_indicator=-1;
    int charge_indicator[comm->nprocs];
    int global_charge_indicator[comm->nprocs];
    memset(charge_indicator, 0, comm->nprocs*sizeof(int));
    memset(global_charge_indicator, 0, comm->nprocs*sizeof(int));
    if (molecule->local_atoms != 0)
    {
       double test_charge = molecule->charge[0];
       int local_tag = molecule->local_tag[0];
       
       if (d_or_s == 0)
       {
          charge_indicator[comm->me] = -1;
          for (int j=0;j<num_charge_states;j++)
          {
             
             if (fabs(test_charge - semiconductor_charges[j][local_tag-1]) < 1e-7)
             {
                charge_indicator[comm->me] = j;
             }
          }
          if (charge_indicator[comm->me] == -1)
          {
             error->all(FLERR, "No charge match for this semiconductor! Check the input file");
          }
       }
       else if (d_or_s == 1)
       {
          charge_indicator[comm->me] = -1;
          for (int j=0;j<num_charge_states;j++)
          {
            if (fabs(test_charge - dopant_charges[j][local_tag-1]) < 1e-7)
            {
               charge_indicator[comm->me] = j;
            }
          }
          if (charge_indicator[comm->me] == -1)
          {
             error->all(FLERR, "No charge match for this dopant! Check the input file");
          }
       }
    }
    else
    {
       charge_indicator[comm->me] = -1;
    }

    MPI_Allreduce(charge_indicator, global_charge_indicator, comm->nprocs, MPI_INT, MPI_SUM, world);
    for (int procs=0;procs<comm->nprocs;procs++)
    {
       if (global_charge_indicator[procs] != -1)
       {
          final_charge_indicator = global_charge_indicator[procs];
       }
    }
    return final_charge_indicator;
}

double FixRMCPartial::change_dihedral_parameters(int molecule_id, int ending_state)
{
   // go to ending_state
   
   double pre_energy = energy_full();
   
   // Circle through the relevant dihedrals and see if any need to be modified.
   for (int d=0;d<num_dihedrals;d++)
   {
      if (dihedral_list[d][4] == molecule_id)
      { // this dihedral type needs to be modified
        // First we need to find the dihedral in the main data structure
        for (int i=0;i<atom->nlocal;i++) 
        {
           if (atom->tag[i] == dihedral_list[d][1])
           {
             for (int j=0;j<atom->num_dihedral[i];j++)
             {
                if (atom->dihedral_atom1[i][j] == dihedral_list[d][0] &&
                    atom->dihedral_atom2[i][j] == dihedral_list[d][1] &&
                    atom->dihedral_atom3[i][j] == dihedral_list[d][2] &&
                    atom->dihedral_atom4[i][j] == dihedral_list[d][3])
                {
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
   for (int a=0;a<num_angles;a++)
   {
      if (angle_list[a][3] == molecule_id)
      { // this dihedral type needs to be modified
        // First we need to find the dihedral in the main data structure
        for (int i=0;i<atom->nlocal;i++) 
        {
           if (atom->tag[i] == angle_list[a][1])
           {
             for (int j=0;j<atom->num_angle[i];j++)
             {
                if (atom->angle_atom1[i][j] == angle_list[a][0] &&
                    atom->angle_atom2[i][j] == angle_list[a][1] &&
                    atom->angle_atom3[i][j] == angle_list[a][2])
                {
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
   for (int b=0;b<num_bonds;b++)
   {
      if (bond_list[b][2] == molecule_id)
      { // this dihedral type needs to be modified
        // First we need to find the dihedral in the main data structure
        for (int i=0;i<atom->nlocal;i++) 
        {
           if (atom->tag[i] == bond_list[b][0])
           {
             for (int j=0;j<atom->num_bond[i];j++)
             {
                if (atom->bond_atom[i][j] == bond_list[b][1])
                {
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
  int lcount_across_ranks[comm->nprocs];
  int lmin_tag[comm->nprocs];
  int gmin_tag[comm->nprocs];
  int gcount_across_ranks[comm->nprocs];
  memset(lmin_tag, 0.0, comm->nprocs * sizeof(int));
  memset(lcount_across_ranks, 0.0, comm->nprocs * sizeof(int));
  int mintag = INT_MAX;
  for (int i=0;i< atom->nlocal;i++)
  {
    if (atom->molecule[i] == (double) mol_id)
    { 
      molecule.charge[atom_counter] = atom->q[i];
      molecule.new_charge[atom_counter] = atom->q[i];
      molecule.type[atom_counter] = atom->type[i];
      molecule.mass[atom_counter] = atom->mass[atom->type[i]];
      molecule.image[atom_counter] = atom->image[i];
      molecule.global_tag[atom_counter] = atom->tag[i];
      molecule.local_index[atom_counter] = i;
      molecule.local_tag[atom_counter] = 0.0;
    
      if (molecule.global_tag[atom_counter] < mintag)
      {
         mintag = molecule.global_tag[atom_counter];
      }
      for (int j=0;j<3;j++) 
      {
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
  for (int i=0;i<comm->nprocs;i++)
  {
    atom_counter_global = atom_counter_global + gcount_across_ranks[i];
    if (gmin_tag[i] != INT_MAX)
    {

      if (gmin_tag[i] < gmintagval)
      {
        gmintagval = gmin_tag[i];
      }
    }
  }
  // Subtract off the minimum global index so we get a 
  // new "local" index. This is not local to processor, but local to a molecule
  // This will be very helpful for charge manipulation
  for (int i=0;i<atom_counter;i++)
  {
    molecule.local_tag[i] = molecule.global_tag[i] - gmintagval + 1;
  }

  // Now, arguably the most important step
  // Determine the charge state of the molecule

  return molecule;
}

void FixRMCPartial::modify_charge(struct Mol *molecule, double *charge_list)
{
  
  for (int i=0;i<molecule->local_atoms;i++)
  {
     for (int j=0;j<atom->nlocal;j++)
     {
        if (atom->tag[j] == molecule->global_tag[i])
        {
           atom->q[j] = charge_list[molecule->local_tag[i]-1];
           molecule->new_charge[i] = charge_list[molecule->local_tag[i]-1];
        }
     }
  }
}

void FixRMCPartial::restore_charge(struct Mol *molecule)
{
  for (int i=0;i<molecule->local_atoms;i++)
  {
     for (int j=0;j<atom->nlocal;j++)
     {
        if (atom->tag[j] == molecule->global_tag[i])
        {
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

    for (int i=0;i<semiconductor_size;i++)
    {
      for (int j=0;j<6;j++)
      {
         osc_mol_c[i][j] = 0.0;
      }
    }
    for (int i=0;i<dopant_size;i++)
    {
      for (int j=0;j<6;j++)
      {
         dopant_mol_c[i][j] = 0.0;
      }
    }

    // Centre of mass pointers
    double *osc_com, *dopant_com;
    osc_com = new double [3];
    dopant_com = new double [3];
    double com_diff;

    // Self energy terms
    double e_osc_diff, e_dopant_diff;


    // Calculate the energy before we do any mischief
    double starting_energy = energy_full(); 
    
    if (comm->me == 0)
    {
       fmt::print(screen, "{}, {}\n", "The starting energy is ", starting_energy);
    }
    
    int rand_semi, rand_dope, rand_charge;
    int indicator=-1;
    int c_indicator=-1;
    
    
    // Find a semiconductor - this is not good of course - you're not sampling properly
    // What you need to do is pick a random CHARGE, and then choose a molecule with that charge
    if (comm->me == 0) {
      while (c_indicator != 0)
      {
         rand_charge = type_dist(rng_type_source);
         if (num_semiconductor_charge[rand_charge] != 0 && num_dopant_charge[rand_charge] != 0)
         {
            c_indicator = 0;
         }
      }
      fmt::print(screen, "{}, {}\n", "We have chosen charge state", charges[rand_charge]);
      int tries=0;
      while (indicator != 0) {
         rand_semi = atom_dist(rng_atom);
         if (molecule_type[rand_semi-1] == 0 && fabs(molecule_charge_states[rand_semi-1] - charges[rand_charge]) < 0.001)
         {
            indicator = 0;
         }
         tries=tries+1;
         if (tries > 10000)
         {
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
          if (molecule_type[rand_dope-1] == 1) 
          {
             if (fabs(molecule_charge_states[rand_dope-1] + molecule_charge_states[rand_semi-1]) < 0.001)
             {
               fmt::print(screen, "{}\n", "FOUND EQUALITY");
               fmt::print(screen, "{} {} {} {}\n", "Semiconductor charge, ", molecule_charge_states[rand_semi-1], " = Dopant charge ", molecule_charge_states[rand_dope-1]);
               indicator = 1;
             }
          }
          tries=tries+1;
          if (tries > 10000)
          {
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
    if (dopant.charge_state != semiconductor.charge_state)
    {
       fmt::print(screen, "{} {} {} {}\n", "Dopant charge state is ", dopant.charge_state, " and semiconductor charge state is ", semiconductor.charge_state); 
       error->all(FLERR, "Charge states don't match!");
    }

    // Identify a destination charge state, picked randomly
    // of course it has to be different from the starting state
    int destination_charge_state = dopant.charge_state;

    if (comm->me == 0)
    {
      while (destination_charge_state == dopant.charge_state)
      {
         destination_charge_state = type_dist(rng_type_destination);
      }
    }

    MPI_Bcast(&destination_charge_state, 1, MPI_INT, 0, world);
    
    if (comm->me == 0)
    {
      fmt::print(screen, "{} {} {} {}\n", "Going from charge state", charges[semiconductor.charge_state], "to ", charges[destination_charge_state]);
    }

    // Calculate centre of mass
    calculateMoleculeCOM(osc_com, &semiconductor);
    calculateMoleculeCOM(dopant_com, &dopant);
    
    com_diff = calDistance(osc_com, dopant_com);

    // Change the dihedral parameters to destination type 
    // and capture the dihedral energy
    edihedral = change_dihedral_parameters(rand_semi, destination_charge_state);
    eangle = change_angle_parameters(rand_semi, destination_charge_state);
    ebond = change_bond_parameters(rand_semi, destination_charge_state);

    // Modify charge to new type
    modify_charge(&semiconductor, semiconductor_charges[destination_charge_state]);
    modify_charge(&dopant, dopant_charges[destination_charge_state]);



   // Capture the molecule self-energies difference
    bringMoleculeTogether(&semiconductor, osc_mol_c, semiconductor_size);
    bringMoleculeTogether(&dopant, dopant_mol_c, dopant_size);


    // calculate special bonds matrix
    
    e_osc_diff = calculateColoumbSelf(osc_mol_c, semiconductor_size);
    e_dopant_diff = calculateColoumbSelf(dopant_mol_c, dopant_size);

    // Capture the pair coloumb interaction
    double interaction_energy[2];
    calculateColoumbCross(osc_mol_c, semiconductor_size, dopant_mol_c, dopant_size, interaction_energy);


    if (comm->me == 0)
    {
       fmt::print(screen, "{} {}\n", "The dihedral energy change is ", edihedral);
       fmt::print(screen, "{} {}\n", "The angle energy change is", eangle);
       fmt::print(screen, "{} {}\n", "The bond energy change is", ebond);
       fmt::print(screen, "{} {}\n", "The semiconductor self energy difference is ", e_osc_diff);
       fmt::print(screen, "{} {}\n", "The dopant self-energy difference is ", e_dopant_diff);
       fmt::print(screen, "{} {} {}\n", "The cross-energy before and after doping are ", interaction_energy[0], interaction_energy[1]);
       fmt::print(screen, "{} {}\n", "The COM distance is ", com_diff);
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
    double energy_diff = new_energy - starting_energy + reaction_energy - e_osc_diff - e_dopant_diff - edihedral;
    double total_electrostatic_diff = new_energy - starting_energy - e_osc_diff - e_dopant_diff - edihedral;

    if (comm->me == 0)
    {
       fmt::print(screen, "{} {}\n",  "Total Electrostatic Energy Difference: ", total_electrostatic_diff);
    }
    // Calculate acceptance probability
    transition_probability = prefactor*exp(-beta*energy_diff);

    if (comm->me == 0)
    {
       fmt::print(screen, "{} {} {} {}\n", "The new energy is ", new_energy, " and the difference is ", energy_diff);

    }

    // Generate random number in one rank
    double rand_number=0;
    if (comm->me == 0)
    {
      rand_number=acc_dist(rng_acc);
    }
    MPI_Bcast(&rand_number, 1, MPI_DOUBLE, 0, world);

    // Determine whether to accept or reject
    if (transition_probability > rand_number)
    {
       // Move accepted
       acceptances = acceptances+1;
       if (comm->me == 0)
       {
          fmt::print(screen, "{}\n", "MOVE ACCEPTED");
       }
     
       num_semiconductor_charge[destination_charge_state] = num_semiconductor_charge[destination_charge_state]+1;
       num_dopant_charge[destination_charge_state] = num_dopant_charge[destination_charge_state]+1;
       num_semiconductor_charge[semiconductor.charge_state] = num_semiconductor_charge[semiconductor.charge_state]-1;
       num_dopant_charge[dopant.charge_state] = num_dopant_charge[dopant.charge_state]-1;
       molecule_charge_states[rand_semi-1] = charges[destination_charge_state];
       molecule_charge_states[rand_dope-1] = -charges[destination_charge_state];
    }
    else
    {
       // Move rejected
       rejections = rejections+1;

       // Revert dihedral coefficients
       edihedral = change_dihedral_parameters(rand_semi, semiconductor.charge_state);
       eangle = change_angle_parameters(rand_semi, semiconductor.charge_state);
       ebond = change_bond_parameters(rand_semi, semiconductor.charge_state);

       if (comm->me == 0)
       { 
          fmt::print(screen, "{}\n", "MOVE REJECTED");
       }

       // Restore the charges to what they were
       restore_charge(&semiconductor);
       restore_charge(&dopant);
       force->kspace->init();
    }
    
    // The step is done, so we can delete the molecules
    delete_molecule(&semiconductor);
    delete_molecule(&dopant);

    // Delete the centre of mass structs
    delete[] osc_com; 
    delete[] dopant_com;

    // Calculate dynamic doping efficiency
    for (int i=0;i<num_charge_states;i++)
    {
         dde[i] = (double) num_dopant_charge[i]/(double)n_dopant;
    }

     if (comm->me == 0)
     {
        fmt::print(screen, "{}", "Dynamic Doping Efficiency: ");
        for (int i=0;i<num_charge_states;i++)
        {
            fmt::print(screen, "{} ", dde[i]);
        }
        fmt::print(screen, "{}\n", " ");
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
    for (int i=0;i<molecule->local_atoms;i++)
    {
       delete[] molecule->pos[i];
    } 
    delete[] molecule->pos;
    delete[] molecule->local_index;
    delete[] molecule->new_charge;
    //delete molecule;
}

void FixRMCPartial::post_mortem()
{
   acceptance_rate = (double)acceptances/(double)nmcsteps;
   double total_charged_dopants=0;
   double total_charged_semiconductors=0;
   for (int i=0;i<num_charge_states;i++)
   {
      doping_efficiency[i] = (double)num_dopant_charge[i]/(double)n_dopant;
      total_charged_dopants = total_charged_dopants+num_dopant_charge[i];
      total_charged_semiconductors = total_charged_semiconductors+num_semiconductor_charge[i];
   }
   double total_doping_efficiency = (double)total_charged_dopants/(double)n_dopant;

   // Write out COM for all molecules, for RDF/morphology or other CG uses

   if (comm->me == 0)
   {
      FILE *fcom = fopen("com.xyz", "w");
      fprintf(fcom, "%d\n",n_molecules);
      for (int j=0;j<3;j++)
      {
         fprintf(fcom, "%f %f\n", domain->boxlo[j], domain->boxhi[j]);
      }
      fprintf(fcom, "LAMMPS_RMC_output\n");
      fclose(fcom);
   }

   double *com;
   com = new double[3];
   for (int i=0;i<n_molecules;i++)
   {
      Mol molecule = get_molecule(i+1,size_limit);
      calculateMoleculeCOM(com, &molecule);
      if (comm->me == 0)
      {
         FILE *fcom = fopen("com.xyz", "a");
         fprintf(fcom, "%d %f %f %f\n", i+1,com[0],com[1],com[2]);
         fclose(fcom);
      }
      delete_molecule(&molecule);  
   }
 
   if (comm->me == 0)
   {
      // Write the charge of each molecule to a file for morphology analysis
      FILE *fp = fopen("molecule_charge.dat", "w");
      for (int i=0;i<n_molecules;i++)
      {
         fprintf(fp, "%f\n", molecule_charge_states[i]);
      }
      fclose(fp);

      // Write the type of each molecule to a file for morphology analysis
      FILE *ftype = fopen("molecule_type.dat", "w");
      for (int i=0;i<n_molecules;i++)
      {
         fprintf(ftype, "%f\n", molecule_type[i]);
      }
      fclose(ftype);

      fmt::print(screen, "{}\n", "###############################################");
      fmt::print(screen, "{}\n", "              RMC OUTPUT SUMMARY               ");
      fmt::print(screen, "{}\n", "###############################################");
      fmt::print(screen,"{} {}\n", "Number of RMC moves: ", nmcsteps);
      fmt::print(screen, "{}", "Final charged dopant number: ");
      for (int i=0;i<num_charge_states;i++)
      {
         fmt::print(screen, "{} ", num_dopant_charge[i]);
      }
      fmt::print(screen, "\n{}", "Final charged semiconductor number: ");
      for (int i=0;i<num_charge_states;i++)
      {
         fmt::print(screen, "{} ", num_semiconductor_charge[i]);
      }
      fmt::print(screen,"\n{} {}\n", "Final charged dopant number: ", total_charged_dopants);
      fmt::print(screen,"{} {}\n", "Final charged semiconductor number: ", total_charged_semiconductors);
      fmt::print(screen,"{}", "Doping efficiency: ");
      for (int i=0;i<num_charge_states;i++)
      {
         fmt::print(screen,"{} ", doping_efficiency[i]);
      }
      fmt::print(screen,"\n{} {}\n", "Final acceptances: ", acceptances);
      fmt::print(screen,"{} {}\n", "Final rejections: ", rejections);
      fmt::print(screen,"{} {}\n", "Acceptance rate: ",acceptance_rate);
      fmt::print(screen, "{}\n", "###############################################");
   }
}


FixRMCPartial::~FixRMCPartial()
{
   post_mortem();
   delete[] dde;

    // Delete the combined structs
    for (int i=0;i<semiconductor_size;i++)
    {
      delete[] osc_mol_c[i];
    }
    for (int i=0;i<dopant_size;i++)
    {
      delete[] dopant_mol_c[i];
    }
    delete[] osc_mol_c;
    delete[] dopant_mol_c;
}
