#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarbias.h"


colvarbias::colvarbias (std::string const &conf, char const *key)
  : colvarparse(), has_data (false)
{
  cvm::log ("Initializing a new \""+std::string (key)+"\" instance.\n");

  size_t rank = 1;
  std::string const key_str (key);

  if (to_lower_cppstr (key_str) == std::string ("abf")) {
    rank = cvm::n_abf_biases+1;
  }
  if (to_lower_cppstr (key_str) == std::string ("harmonic")) {
    rank = cvm::n_harm_biases+1;
  }
  if (to_lower_cppstr (key_str) == std::string ("histogram")) {
    rank = cvm::n_histo_biases+1;
  }
  if (to_lower_cppstr (key_str) == std::string ("metadynamics")) {
    rank = cvm::n_meta_biases+1;
  }

  get_keyval (conf, "name", name, key_str+cvm::to_str (rank));

  for (std::vector<colvarbias *>::iterator bi = cvm::biases.begin();
       bi != cvm::biases.end();
       bi++) {
    if ((*bi)->name == this->name)
      cvm::fatal_error ("Error: this bias cannot have the same name, \""+this->name+
                        "\", as another bias.\n");
  }

  // lookup the associated colvars
  std::vector<std::string> colvars_str;
  if (get_keyval (conf, "colvars", colvars_str)) {
    for (size_t i = 0; i < colvars_str.size(); i++) {
      add_colvar (colvars_str[i]);
    }
  }
  if (!colvars.size()) {
    cvm::fatal_error ("Error: no collective variables specified.\n");
  }

  get_keyval (conf, "outputEnergy", b_output_energy, false);
}


colvarbias::colvarbias()
  : colvarparse(), has_data (false)
{}


void colvarbias::add_colvar (std::string const &cv_name)
{
  if (colvar *cvp = cvm::colvar_p (cv_name)) {
    cvp->enable (colvar::task_gradients);
    if (cvm::debug())
      cvm::log ("Applying this bias to collective variable \""+
                cvp->name+"\".\n");
    colvars.push_back (cvp);
    colvar_forces.push_back (colvarvalue (cvp->type()));
  } else {
    cvm::fatal_error ("Error: cannot find a colvar named \""+
                      cv_name+"\".\n");
  }
}


void colvarbias::communicate_forces()
{
  for (size_t i = 0; i < colvars.size(); i++) {
    if (cvm::debug()) {
      cvm::log ("Communicating a force to colvar \""+
                colvars[i]->name+"\", of type \""+
                colvarvalue::type_desc[colvars[i]->type()]+"\".\n");
    }
    colvars[i]->add_bias_force (colvar_forces[i]);
  }
}


void colvarbias::change_configuration(std::string const &conf)
{
  cvm::fatal_error ("Error: change_configuration() not implemented.\n");
}


cvm::real colvarbias::energy_difference(std::string const &conf)
{
  cvm::fatal_error ("Error: energy_difference() not implemented.\n");
  return 0.;
}


std::ostream & colvarbias::write_traj_label (std::ostream &os)
{
  os << " ";
  if (b_output_energy)
    os << " E_"
       << cvm::wrap_string (this->name, cvm::en_width-2);
  return os;
}


std::ostream & colvarbias::write_traj (std::ostream &os)
{
  os << " ";
  if (b_output_energy)
    os << " "
       << bias_energy;
  return os;
}




colvarbias_harmonic::colvarbias_harmonic (std::string const &conf,
                                          char const *key)
  : colvarbias (conf, key),
    target_nstages (0),
    target_nsteps (0)
{
  get_keyval (conf, "forceConstant", force_k, 1.0);
  for (size_t i = 0; i < colvars.size(); i++) {
    if (colvars[i]->width != 1.0)
      cvm::log ("The force constant for colvar \""+colvars[i]->name+
                "\" will be rescaled to "+
                cvm::to_str (force_k/(colvars[i]->width*colvars[i]->width))+
                " according to the specified width.\n");
  }

  // get the initial restraint centers
  colvar_centers.resize (colvars.size());
  colvar_centers_raw.resize (colvars.size());
  for (size_t i = 0; i < colvars.size(); i++) {
    colvar_centers[i].type (colvars[i]->type());
    colvar_centers_raw[i].type (colvars[i]->type());
  }
  if (get_keyval (conf, "centers", colvar_centers, colvar_centers)) {
    for (size_t i = 0; i < colvars.size(); i++) {
      colvar_centers[i].apply_constraints();
      colvar_centers_raw[i] = colvar_centers[i];
    }
  } else {
    colvar_centers.clear();
    cvm::fatal_error ("Error: must define the initial centers of the restraints.\n");
  }

  if (colvar_centers.size() != colvars.size())
    cvm::fatal_error ("Error: number of harmonic centers does not match "
                      "that of collective variables.\n");

  if (get_keyval (conf, "targetCenters", target_centers, colvar_centers)) {
    b_chg_centers = true;
    for (size_t i = 0; i < target_centers.size(); i++) {
      target_centers[i].apply_constraints();
    }
  } else {
    b_chg_centers = false;
    target_centers.clear();
  }

  if (get_keyval (conf, "targetForceConstant", target_force_k, 0.0)) {
    if (b_chg_centers)
      cvm::fatal_error ("Error: cannot specify both targetCenters and targetForceConstant.\n");

    starting_force_k = force_k;
    b_chg_force_k = true;

    get_keyval (conf, "targetEquilSteps", target_equil_steps, 0);

    get_keyval (conf, "lambdaSchedule", lambda_schedule, lambda_schedule);
    if (lambda_schedule.size()) {
      // There is one more lambda-point than stages
      target_nstages = lambda_schedule.size() - 1;
    }
  } else {
    b_chg_force_k = false;
  }

  if (b_chg_centers || b_chg_force_k) {
    get_keyval (conf, "targetNumSteps", target_nsteps, 0);
    if (!target_nsteps)
      cvm::fatal_error ("Error: targetNumSteps must be non-zero.\n");

    if (get_keyval (conf, "targetNumStages", target_nstages, target_nstages) &&
        lambda_schedule.size()) {
      cvm::fatal_error ("Error: targetNumStages and lambdaSchedule are incompatible.\n");
    }

    if (target_nstages) {
      // This means that either numStages of lambdaSchedule has been provided
      stage = 0;
      restraint_FE = 0.0;
    }

    if (get_keyval (conf, "targetForceExponent", force_k_exp, 1.0)) {
      if (! b_chg_force_k)
        cvm::log ("Warning: not changing force constant: targetForceExponent will be ignored\n");
      if (force_k_exp < 1.0)
        cvm::log ("Warning: for all practical purposes, targetForceExponent should be 1.0 or greater.\n");
    }
  }

  get_keyval (conf, "outputCenters", b_output_centers, false);
  get_keyval (conf, "outputAccumulatedWork", b_output_acc_work, false);
  acc_work = 0.0;

  if (cvm::debug())
    cvm::log ("Done initializing a new harmonic restraint bias.\n");
}

colvarbias_harmonic::~colvarbias_harmonic ()
{
  if (cvm::n_harm_biases > 0)
    cvm::n_harm_biases -= 1;
}


void colvarbias_harmonic::change_configuration (std::string const &conf)
{
  get_keyval (conf, "forceConstant", force_k, force_k);
  if (get_keyval (conf, "centers", colvar_centers, colvar_centers)) {
    for (size_t i = 0; i < colvars.size(); i++) {
      colvar_centers[i].apply_constraints();
      colvar_centers_raw[i] = colvar_centers[i];
    }
  }
}


cvm::real colvarbias_harmonic::energy_difference (std::string const &conf)
{
  std::vector<colvarvalue> alt_colvar_centers;
  cvm::real alt_force_k;
  cvm::real alt_bias_energy = 0.0;

  get_keyval (conf, "forceConstant", alt_force_k, force_k);

  alt_colvar_centers.resize (colvars.size());
  for (size_t i = 0; i < colvars.size(); i++) {
    alt_colvar_centers[i].type (colvars[i]->type());
  }
  if (get_keyval (conf, "centers", alt_colvar_centers, colvar_centers)) {
    for (size_t i = 0; i < colvars.size(); i++) {
      colvar_centers[i].apply_constraints();
    }
  }

  for (size_t i = 0; i < colvars.size(); i++) {
    alt_bias_energy += 0.5 * alt_force_k / (colvars[i]->width * colvars[i]->width) *
      colvars[i]->dist2 (colvars[i]->value(), alt_colvar_centers[i]);
  }

  return alt_bias_energy - bias_energy;
}


cvm::real colvarbias_harmonic::update()
{
  bias_energy = 0.0;

  if (cvm::debug())
    cvm::log ("Updating the harmonic bias \""+this->name+"\".\n");

  // Setup first stage of staged variable force constant calculation
  if (b_chg_force_k && target_nstages && cvm::step_absolute() == 0) {
    cvm::real lambda;
    if (lambda_schedule.size()) {
      lambda = lambda_schedule[0];
    } else {
      lambda = 0.0;
    }
    force_k = starting_force_k + (target_force_k - starting_force_k)
              * std::pow (lambda, force_k_exp);
    cvm::log ("Harmonic restraint " + this->name + ", stage " +
        cvm::to_str(stage) + " : lambda = " + cvm::to_str(lambda));
    cvm::log ("Setting force constant to " + cvm::to_str (force_k));
  }

  if (b_chg_centers) {

    if (!centers_incr.size()) {
      // if this is the first calculation, calculate the advancement
      // at each simulation step (or stage, if applicable)
      // (take current stage into account: it can be non-zero
      //  if we are restarting a staged calculation)
      centers_incr.resize (colvars.size());
      for (size_t i = 0; i < colvars.size(); i++) {
        centers_incr[i].type (colvars[i]->type());
        centers_incr[i] = (target_centers[i] - colvar_centers_raw[i]) /
          cvm::real ( target_nstages ? (target_nstages - stage) :
                                      (target_nsteps - cvm::step_absolute()));
      }
      if (cvm::debug())
        cvm::log ("Center increment for the harmonic bias \""+
                  this->name+"\": "+cvm::to_str (centers_incr)+" at stage "+cvm::to_str (stage)+ ".\n");

    }

    if (target_nstages) {
      if ((cvm::step_relative() > 0)
            && (cvm::step_absolute() % target_nsteps) == 0
            && stage < target_nstages) {

          for (size_t i = 0; i < colvars.size(); i++) {
            colvar_centers_raw[i] += centers_incr[i];
            colvar_centers[i] = colvar_centers_raw[i];
            colvars[i]->wrap(colvar_centers[i]);
            colvar_centers[i].apply_constraints();
          }
          stage++;
          cvm::log ("Moving restraint stage " + cvm::to_str(stage) +
              " : setting centers to " + cvm::to_str (colvar_centers) +
              " at step " +  cvm::to_str (cvm::step_absolute()));
      }
    } else if ((cvm::step_relative() > 0) && (cvm::step_absolute() <= target_nsteps)) {
      // move the restraint centers in the direction of the targets
      // (slow growth)
      for (size_t i = 0; i < colvars.size(); i++) {
        colvar_centers_raw[i] += centers_incr[i];
        colvar_centers[i] = colvar_centers_raw[i];
        colvars[i]->wrap(colvar_centers[i]);
        colvar_centers[i].apply_constraints();
      }
    }

    if (cvm::debug())
      cvm::log ("Current centers for the harmonic bias \""+
                this->name+"\": "+cvm::to_str (colvar_centers)+".\n");
  }

  if (b_chg_force_k) {
    // Coupling parameter, between 0 and 1
    cvm::real lambda;

    if (target_nstages) {
      // TI calculation: estimate free energy derivative
      // need current lambda
      if (lambda_schedule.size()) {
        lambda = lambda_schedule[stage];
      } else {
        lambda = cvm::real(stage) / cvm::real(target_nstages);
      }

      if (target_equil_steps == 0 || cvm::step_absolute() % target_nsteps >= target_equil_steps) {
        // Start averaging after equilibration period, if requested

        // Square distance normalized by square colvar width
        cvm::real dist_sq = 0.0;
        for (size_t i = 0; i < colvars.size(); i++) {
          dist_sq += colvars[i]->dist2 (colvars[i]->value(), colvar_centers[i])
            / (colvars[i]->width * colvars[i]->width);
        }

        restraint_FE += 0.5 * force_k_exp * std::pow(lambda, force_k_exp - 1.0)
          * (target_force_k - starting_force_k) * dist_sq;
      }

      // Finish current stage...
      if (cvm::step_absolute() % target_nsteps == 0 &&
          cvm::step_absolute() > 0) {

          cvm::log ("Lambda= " + cvm::to_str (lambda) + " dA/dLambda= "
              + cvm::to_str (restraint_FE / cvm::real(target_nsteps - target_equil_steps)));

        //  ...and move on to the next one
        if (stage < target_nstages) {

          restraint_FE = 0.0;
          stage++;
          if (lambda_schedule.size()) {
            lambda = lambda_schedule[stage];
          } else {
            lambda = cvm::real(stage) / cvm::real(target_nstages);
          }
          force_k = starting_force_k + (target_force_k - starting_force_k)
                    * std::pow (lambda, force_k_exp);
          cvm::log ("Harmonic restraint " + this->name + ", stage " +
              cvm::to_str(stage) + " : lambda = " + cvm::to_str(lambda));
          cvm::log ("Setting force constant to " + cvm::to_str (force_k));
        }
      }
    } else if (cvm::step_absolute() <= target_nsteps) {
      // update force constant (slow growth)
      lambda = cvm::real(cvm::step_absolute()) / cvm::real(target_nsteps);
      force_k = starting_force_k + (target_force_k - starting_force_k)
          * std::pow (lambda, force_k_exp);
    }
  }

  if (cvm::debug())
    cvm::log ("Done updating the harmonic bias \""+this->name+"\".\n");

  // Force and energy calculation
  for (size_t i = 0; i < colvars.size(); i++) {
    colvar_forces[i] =
      (-0.5) * force_k /
      (colvars[i]->width * colvars[i]->width) *
      colvars[i]->dist2_lgrad (colvars[i]->value(),
                               colvar_centers[i]);
    bias_energy += 0.5 * force_k / (colvars[i]->width * colvars[i]->width) *
              colvars[i]->dist2(colvars[i]->value(), colvar_centers[i]);
    if (cvm::debug())
      cvm::log ("dist_grad["+cvm::to_str (i)+
                "] = "+cvm::to_str (colvars[i]->dist2_lgrad (colvars[i]->value(),
                               colvar_centers[i]))+"\n");
  }

  if (b_output_acc_work) {
    if ((cvm::step_relative() > 0) || (cvm::step_absolute() == 0)) {
      for (size_t i = 0; i < colvars.size(); i++) {
        // project forces on the calculated increments at this step
        acc_work += colvar_forces[i] * centers_incr[i];
      }
    }
  }

  if (cvm::debug())
    cvm::log ("Current forces for the harmonic bias \""+
              this->name+"\": "+cvm::to_str (colvar_forces)+".\n");

  return bias_energy;
}


std::istream & colvarbias_harmonic::read_restart (std::istream &is)
{
  size_t const start_pos = is.tellg();

  cvm::log ("Restarting harmonic bias \""+
            this->name+"\".\n");

  std::string key, brace, conf;
  if ( !(is >> key)   || !(key == "harmonic") ||
       !(is >> brace) || !(brace == "{") ||
       !(is >> colvarparse::read_block ("configuration", conf)) ) {

    cvm::log ("Error: in reading restart configuration for harmonic bias \""+
              this->name+"\" at position "+
              cvm::to_str (is.tellg())+" in stream.\n");
    is.clear();
    is.seekg (start_pos, std::ios::beg);
    is.setstate (std::ios::failbit);
    return is;
  }

//   int id = -1;
  std::string name = "";
//   if ( ( (colvarparse::get_keyval (conf, "id", id, -1, colvarparse::parse_silent)) &&
//          (id != this->id) ) ||
  if ( (colvarparse::get_keyval (conf, "name", name, std::string (""), colvarparse::parse_silent)) &&
       (name != this->name) )
    cvm::fatal_error ("Error: in the restart file, the "
                      "\"harmonic\" block has a wrong name\n");
//   if ( (id == -1) && (name == "") ) {
  if (name.size() == 0) {
    cvm::fatal_error ("Error: \"harmonic\" block in the restart file "
                      "has no identifiers.\n");
  }

  if (b_chg_centers) {
//    cvm::log ("Reading the updated restraint centers from the restart.\n");
    if (!get_keyval (conf, "centers", colvar_centers))
      cvm::fatal_error ("Error: restraint centers are missing from the restart.\n");
    if (!get_keyval (conf, "centers_raw", colvar_centers_raw))
      cvm::fatal_error ("Error: \"raw\" restraint centers are missing from the restart.\n");
  }

  if (b_chg_force_k) {
//    cvm::log ("Reading the updated force constant from the restart.\n");
    if (!get_keyval (conf, "forceConstant", force_k))
      cvm::fatal_error ("Error: force constant is missing from the restart.\n");
  }

  if (target_nstages) {
//    cvm::log ("Reading current stage from the restart.\n");
    if (!get_keyval (conf, "stage", stage))
      cvm::fatal_error ("Error: current stage is missing from the restart.\n");
  }

  if (b_output_acc_work) {
    if (!get_keyval (conf, "accumulatedWork", acc_work))
      cvm::fatal_error ("Error: accumulatedWork is missing from the restart.\n");
  }

  is >> brace;
  if (brace != "}") {
    cvm::fatal_error ("Error: corrupt restart information for harmonic bias \""+
                      this->name+"\": no matching brace at position "+
                      cvm::to_str (is.tellg())+" in the restart file.\n");
    is.setstate (std::ios::failbit);
  }
  return is;
}


std::ostream & colvarbias_harmonic::write_restart (std::ostream &os)
{
  os << "harmonic {\n"
     << "  configuration {\n"
    //      << "    id " << this->id << "\n"
     << "    name " << this->name << "\n";

  if (b_chg_centers) {
    os << "    centers ";
    for (size_t i = 0; i < colvars.size(); i++) {
      os << " " << colvar_centers[i];
    }
    os << "\n";
    os << "    centers_raw ";
    for (size_t i = 0; i < colvars.size(); i++) {
      os << " " << colvar_centers_raw[i];
    }
    os << "\n";
  }

  if (b_chg_force_k) {
    os << "    forceConstant "
       << std::setprecision (cvm::en_prec)
       << std::setw (cvm::en_width) << force_k << "\n";
  }

  if (target_nstages) {
    os << "    stage " << std::setw (cvm::it_width)
       << stage << "\n";
  }

  if (b_output_acc_work) {
    os << "    accumulatedWork " << acc_work << "\n";
  }

  os << "  }\n"
     << "}\n\n";

  return os;
}


std::ostream & colvarbias_harmonic::write_traj_label (std::ostream &os)
{
  os << " ";

  if (b_output_energy)
    os << " E_"
       << cvm::wrap_string (this->name, cvm::en_width-2);

  if (b_output_centers)
    for (size_t i = 0; i < colvars.size(); i++) {
      size_t const this_cv_width = (colvars[i]->value()).output_width (cvm::cv_width);
      os << " x0_"
         << cvm::wrap_string (colvars[i]->name, this_cv_width-3);
    }

  if (b_output_acc_work)
    os << " W_"
       << cvm::wrap_string (this->name, cvm::en_width-2);

  return os;
}


std::ostream & colvarbias_harmonic::write_traj (std::ostream &os)
{
  os << " ";

  if (b_output_energy)
    os << " "
       << std::setprecision (cvm::en_prec) << std::setw (cvm::en_width)
       << bias_energy;

  if (b_output_centers)
    for (size_t i = 0; i < colvars.size(); i++) {
      os << " "
         << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
         << colvar_centers[i];
    }

  if (b_output_acc_work)
    os << " "
       << std::setprecision (cvm::en_prec) << std::setw (cvm::en_width)
       << acc_work;

  return os;
}

