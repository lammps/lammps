/// -*- c++ -*-

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarbias_restraint.h"


colvarbias_restraint::colvarbias_restraint(std::string const &conf,
                                           char const *key)
  : colvarbias(conf, key),
    target_nstages(0),
    target_nsteps(0)
{
  get_keyval(conf, "forceConstant", force_k, 1.0);

  {
    // get the initial restraint centers
    colvar_centers.resize(colvars.size());
    colvar_centers_raw.resize(colvars.size());
    size_t i;
    for (i = 0; i < colvars.size(); i++) {
      colvar_centers[i].type(colvars[i]->value());
      colvar_centers_raw[i].type(colvars[i]->value());
      if (cvm::debug()) {
        cvm::log("colvarbias_restraint: center size = "+
                 cvm::to_str(colvar_centers[i].vector1d_value.size())+"\n");
      }
    }
    if (get_keyval(conf, "centers", colvar_centers, colvar_centers)) {
      for (i = 0; i < colvars.size(); i++) {
        if (cvm::debug()) {
          cvm::log("colvarbias_restraint: parsing initial centers, i = "+cvm::to_str(i)+".\n");
        }

        colvar_centers[i].apply_constraints();
        colvar_centers_raw[i] = colvar_centers[i];
      }
    } else {
      colvar_centers.clear();
      cvm::error("Error: must define the initial centers of the restraints.\n");
    }

    if (colvar_centers.size() != colvars.size()) {
      cvm::error("Error: number of centers does not match "
                 "that of collective variables.\n");
    }
  }

  {
    if (cvm::debug()) {
      cvm::log("colvarbias_restraint: parsing target centers.\n");
    }

    size_t i;
    if (get_keyval(conf, "targetCenters", target_centers, colvar_centers)) {
      b_chg_centers = true;
      for (i = 0; i < target_centers.size(); i++) {
        target_centers[i].apply_constraints();
      }
    } else {
      b_chg_centers = false;
      target_centers.clear();
    }
  }

  if (get_keyval(conf, "targetForceConstant", target_force_k, 0.0)) {
    if (b_chg_centers)
      cvm::error("Error: cannot specify both targetCenters and targetForceConstant.\n");

    starting_force_k = force_k;
    b_chg_force_k = true;

    get_keyval(conf, "targetEquilSteps", target_equil_steps, 0);

    get_keyval(conf, "lambdaSchedule", lambda_schedule, lambda_schedule);
    if (lambda_schedule.size()) {
      // There is one more lambda-point than stages
      target_nstages = lambda_schedule.size() - 1;
    }
  } else {
    b_chg_force_k = false;
  }

  if (b_chg_centers || b_chg_force_k) {
    get_keyval(conf, "targetNumSteps", target_nsteps, 0);
    if (!target_nsteps)
      cvm::error("Error: targetNumSteps must be non-zero.\n");

    if (get_keyval(conf, "targetNumStages", target_nstages, target_nstages) &&
        lambda_schedule.size()) {
      cvm::error("Error: targetNumStages and lambdaSchedule are incompatible.\n");
    }

    if (target_nstages) {
      // This means that either numStages of lambdaSchedule has been provided
      stage = 0;
      restraint_FE = 0.0;
    }

    if (get_keyval(conf, "targetForceExponent", force_k_exp, 1.0)) {
      if (! b_chg_force_k)
        cvm::log("Warning: not changing force constant: targetForceExponent will be ignored\n");
      if (force_k_exp < 1.0)
        cvm::log("Warning: for all practical purposes, targetForceExponent should be 1.0 or greater.\n");
    }
  }

  get_keyval(conf, "outputCenters", b_output_centers, false);
  if (b_chg_centers) {
    get_keyval(conf, "outputAccumulatedWork", b_output_acc_work, false);
  } else {
    b_output_acc_work = false;
  }
  acc_work = 0.0;

  if (cvm::debug())
    cvm::log("Done initializing a new restraint bias.\n");
}

colvarbias_restraint::~colvarbias_restraint()
{
  if (cvm::n_rest_biases > 0)
    cvm::n_rest_biases -= 1;
}


void colvarbias_restraint::change_configuration(std::string const &conf)
{
  get_keyval(conf, "forceConstant", force_k, force_k);
  if (get_keyval(conf, "centers", colvar_centers, colvar_centers)) {
    for (size_t i = 0; i < colvars.size(); i++) {
      colvar_centers[i].type(colvars[i]->value());
      colvar_centers[i].apply_constraints();
      colvar_centers_raw[i].type(colvars[i]->value());
      colvar_centers_raw[i] = colvar_centers[i];
    }
  }
}


cvm::real colvarbias_restraint::energy_difference(std::string const &conf)
{
  std::vector<colvarvalue> alt_colvar_centers;
  cvm::real alt_force_k;
  cvm::real alt_bias_energy = 0.0;

  get_keyval(conf, "forceConstant", alt_force_k, force_k);

  alt_colvar_centers.resize(colvars.size());
  size_t i;
  for (i = 0; i < colvars.size(); i++) {
    alt_colvar_centers[i].type(colvars[i]->value());
  }
  if (get_keyval(conf, "centers", alt_colvar_centers, colvar_centers)) {
    for (i = 0; i < colvars.size(); i++) {
      colvar_centers[i].type(colvars[i]->value());
      colvar_centers[i].apply_constraints();
    }
  }

  for (i = 0; i < colvars.size(); i++) {
    alt_bias_energy += restraint_potential(restraint_convert_k(alt_force_k, colvars[i]->width),
					   colvars[i],
					   alt_colvar_centers[i]);
  }

  return alt_bias_energy - bias_energy;
}


cvm::real colvarbias_restraint::update()
{
  bias_energy = 0.0;

  if (cvm::debug())
    cvm::log("Updating the restraint bias \""+this->name+"\".\n");

  // Setup first stage of staged variable force constant calculation
  if (b_chg_force_k && target_nstages && cvm::step_absolute() == 0) {
    cvm::real lambda;
    if (lambda_schedule.size()) {
      lambda = lambda_schedule[0];
    } else {
      lambda = 0.0;
    }
    force_k = starting_force_k + (target_force_k - starting_force_k)
              * std::pow(lambda, force_k_exp);
    cvm::log("Restraint " + this->name + ", stage " +
        cvm::to_str(stage) + " : lambda = " + cvm::to_str(lambda));
    cvm::log("Setting force constant to " + cvm::to_str(force_k));
  }

  if (b_chg_centers) {

    if (!centers_incr.size()) {
      // if this is the first calculation, calculate the advancement
      // at each simulation step (or stage, if applicable)
      // (take current stage into account: it can be non-zero
      //  if we are restarting a staged calculation)
      centers_incr.resize(colvars.size());
      for (size_t i = 0; i < colvars.size(); i++) {
        centers_incr[i].type(colvars[i]->value());
        centers_incr[i] = (target_centers[i] - colvar_centers_raw[i]) /
          cvm::real( target_nstages ? (target_nstages - stage) :
                                      (target_nsteps - cvm::step_absolute()));
      }
      if (cvm::debug()) {
        cvm::log("Center increment for the restraint bias \""+
                  this->name+"\": "+cvm::to_str(centers_incr)+" at stage "+cvm::to_str(stage)+ ".\n");
      }
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
          cvm::log("Moving restraint \"" + this->name +
              "\" stage " + cvm::to_str(stage) +
              " : setting centers to " + cvm::to_str(colvar_centers) +
              " at step " +  cvm::to_str(cvm::step_absolute()));
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
      cvm::log("Current centers for the restraint bias \""+
                this->name+"\": "+cvm::to_str(colvar_centers)+".\n");
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
          dist_sq += colvars[i]->dist2(colvars[i]->value(), colvar_centers[i])
            / (colvars[i]->width * colvars[i]->width);
        }

        restraint_FE += 0.5 * force_k_exp * std::pow(lambda, force_k_exp - 1.0)
          * (target_force_k - starting_force_k) * dist_sq;
      }

      // Finish current stage...
      if (cvm::step_absolute() % target_nsteps == 0 &&
          cvm::step_absolute() > 0) {

          cvm::log("Lambda= " + cvm::to_str(lambda) + " dA/dLambda= "
              + cvm::to_str(restraint_FE / cvm::real(target_nsteps - target_equil_steps)));

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
                    * std::pow(lambda, force_k_exp);
          cvm::log("Restraint " + this->name + ", stage " +
              cvm::to_str(stage) + " : lambda = " + cvm::to_str(lambda));
          cvm::log("Setting force constant to " + cvm::to_str(force_k));
        }
      }
    } else if (cvm::step_absolute() <= target_nsteps) {
      // update force constant (slow growth)
      lambda = cvm::real(cvm::step_absolute()) / cvm::real(target_nsteps);
      force_k = starting_force_k + (target_force_k - starting_force_k)
          * std::pow(lambda, force_k_exp);
    }
  }

  if (cvm::debug())
    cvm::log("Done updating the restraint bias \""+this->name+"\".\n");

  // Force and energy calculation
  for (size_t i = 0; i < colvars.size(); i++) {
    colvar_forces[i].type(colvars[i]->value());
    colvar_forces[i] = -1.0 * restraint_force(restraint_convert_k(force_k, colvars[i]->width),
                                              colvars[i],
                                              colvar_centers[i]);
    bias_energy += restraint_potential(restraint_convert_k(force_k, colvars[i]->width),
				       colvars[i],
				       colvar_centers[i]);
    if (cvm::debug()) {
      cvm::log("dist_grad["+cvm::to_str(i)+
                "] = "+cvm::to_str(colvars[i]->dist2_lgrad(colvars[i]->value(),
                                                             colvar_centers[i]))+"\n");
    }
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
    cvm::log("Current forces for the restraint bias \""+
              this->name+"\": "+cvm::to_str(colvar_forces)+".\n");

  return bias_energy;
}


std::istream & colvarbias_restraint::read_restart(std::istream &is)
{
  size_t const start_pos = is.tellg();

  cvm::log("Restarting restraint bias \""+
            this->name+"\".\n");

  std::string key, brace, conf;
  if ( !(is >> key)   || !(key == "restraint" || key == "harmonic") ||
       !(is >> brace) || !(brace == "{") ||
       !(is >> colvarparse::read_block("configuration", conf)) ) {

    cvm::log("Error: in reading restart configuration for restraint bias \""+
              this->name+"\" at position "+
              cvm::to_str(is.tellg())+" in stream.\n");
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
    return is;
  }

//   int id = -1;
  std::string name = "";
//   if ( ( (colvarparse::get_keyval (conf, "id", id, -1, colvarparse::parse_silent)) &&
//          (id != this->id) ) ||
  if ( (colvarparse::get_keyval(conf, "name", name, std::string(""), colvarparse::parse_silent)) &&
       (name != this->name) )
    cvm::error("Error: in the restart file, the "
                      "\"restraint\" block has a wrong name\n");
//   if ( (id == -1) && (name == "") ) {
  if (name.size() == 0) {
    cvm::error("Error: \"restraint\" block in the restart file "
                      "has no identifiers.\n");
  }

  if (b_chg_centers) {
//    cvm::log ("Reading the updated restraint centers from the restart.\n");
    if (!get_keyval(conf, "centers", colvar_centers))
      cvm::error("Error: restraint centers are missing from the restart.\n");
    if (!get_keyval(conf, "centers_raw", colvar_centers_raw))
      cvm::error("Error: \"raw\" restraint centers are missing from the restart.\n");
  }

  if (b_chg_force_k) {
//    cvm::log ("Reading the updated force constant from the restart.\n");
    if (!get_keyval(conf, "forceConstant", force_k))
      cvm::error("Error: force constant is missing from the restart.\n");
  }

  if (target_nstages) {
//    cvm::log ("Reading current stage from the restart.\n");
    if (!get_keyval(conf, "stage", stage))
      cvm::error("Error: current stage is missing from the restart.\n");
  }

  if (b_output_acc_work) {
    if (!get_keyval(conf, "accumulatedWork", acc_work))
      cvm::error("Error: accumulatedWork is missing from the restart.\n");
  }

  is >> brace;
  if (brace != "}") {
    cvm::error("Error: corrupt restart information for restraint bias \""+
                      this->name+"\": no matching brace at position "+
                      cvm::to_str(is.tellg())+" in the restart file.\n");
    is.setstate(std::ios::failbit);
  }
  return is;
}


std::ostream & colvarbias_restraint::write_restart(std::ostream &os)
{
  os << "restraint {\n"
     << "  configuration {\n"
    //      << "    id " << this->id << "\n"
     << "    name " << this->name << "\n";

  if (b_chg_centers) {
    size_t i;
    os << "    centers ";
    for (i = 0; i < colvars.size(); i++) {
      os << " " << colvar_centers[i];
    }
    os << "\n";
    os << "    centers_raw ";
    for (i = 0; i < colvars.size(); i++) {
      os << " " << colvar_centers_raw[i];
    }
    os << "\n";
  }

  if (b_chg_force_k) {
    os << "    forceConstant "
       << std::setprecision(cvm::en_prec)
       << std::setw(cvm::en_width) << force_k << "\n";
  }

  if (target_nstages) {
    os << "    stage " << std::setw(cvm::it_width)
       << stage << "\n";
  }

  if (b_output_acc_work) {
    os << "    accumulatedWork " << acc_work << "\n";
  }

  os << "  }\n"
     << "}\n\n";

  return os;
}


std::ostream & colvarbias_restraint::write_traj_label(std::ostream &os)
{
  os << " ";

  if (b_output_energy)
    os << " E_"
       << cvm::wrap_string(this->name, cvm::en_width-2);

  if (b_output_centers)
    for (size_t i = 0; i < colvars.size(); i++) {
      size_t const this_cv_width = (colvars[i]->value()).output_width(cvm::cv_width);
      os << " x0_"
         << cvm::wrap_string(colvars[i]->name, this_cv_width-3);
    }

  if (b_output_acc_work)
    os << " W_"
       << cvm::wrap_string(this->name, cvm::en_width-2);

  return os;
}


std::ostream & colvarbias_restraint::write_traj(std::ostream &os)
{
  os << " ";

  if (b_output_energy)
    os << " "
       << std::setprecision(cvm::en_prec) << std::setw(cvm::en_width)
       << bias_energy;

  if (b_output_centers)
    for (size_t i = 0; i < colvars.size(); i++) {
      os << " "
         << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
         << colvar_centers[i];
    }

  if (b_output_acc_work)
    os << " "
       << std::setprecision(cvm::en_prec) << std::setw(cvm::en_width)
       << acc_work;

  return os;
}

colvarbias_restraint_harmonic::colvarbias_restraint_harmonic(std::string const &conf, char const *key) :
  colvarbias_restraint(conf, key) {
  for (size_t i = 0; i < colvars.size(); i++) {
    if (colvars[i]->width != 1.0)
      cvm::log("The force constant for colvar \""+colvars[i]->name+
                "\" will be rescaled to "+
                cvm::to_str(restraint_convert_k(force_k, colvars[i]->width))+
                " according to the specified width.\n");
  }
}

cvm::real colvarbias_restraint_harmonic::restraint_potential(cvm::real k,  colvar* x,  const colvarvalue &xcenter) const
{
  return 0.5 * k * x->dist2(x->value(), xcenter);
}

colvarvalue colvarbias_restraint_harmonic::restraint_force(cvm::real k,  colvar* x,  const colvarvalue &xcenter) const
{
  return 0.5 * k * x->dist2_lgrad(x->value(), xcenter);
}

cvm::real colvarbias_restraint_harmonic::restraint_convert_k(cvm::real k, cvm::real dist_measure) const
{
  return k / (dist_measure * dist_measure);
}


colvarbias_restraint_linear::colvarbias_restraint_linear(std::string const &conf, char const *key) :
  colvarbias_restraint(conf, key) {
  for (size_t i = 0; i < colvars.size(); i++) {
    if (colvars[i]->width != 1.0)
      cvm::log("The force constant for colvar \""+colvars[i]->name+
                "\" will be rescaled to "+
                cvm::to_str(restraint_convert_k(force_k, colvars[i]->width))+
                " according to the specified width.\n");
  }
}

cvm::real colvarbias_restraint_linear::restraint_potential(cvm::real k,  colvar* x,  const colvarvalue &xcenter) const
{
  return k * (x->value() - xcenter);
}

colvarvalue colvarbias_restraint_linear::restraint_force(cvm::real k,  colvar* x,  const colvarvalue &xcenter) const
{
  return k;
}

cvm::real colvarbias_restraint_linear::restraint_convert_k(cvm::real k, cvm::real dist_measure) const
{
  return k / dist_measure;
}
