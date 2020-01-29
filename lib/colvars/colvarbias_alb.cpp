// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <cstdlib>

#include "colvarmodule.h"
#include "colvarbias.h"
#include "colvarbias_alb.h"

#ifdef _MSC_VER
#if _MSC_VER <= 1700
#define copysign(A,B) _copysign(A,B)
double fmax(double A, double B) { return ( A > B ? A : B ); }
double fmin(double A, double B) { return ( A < B ? A : B ); }
#endif
#endif

/* Note about nomenclature. Force constant is called a coupling
 * constant here to emphasize its changing in the code. Outwards,
 * everything is called a force constant to keep it consistent with
 * the rest of colvars.
 *
 */

colvarbias_alb::colvarbias_alb(char const *key)
  : colvarbias(key), update_calls(0), b_equilibration(true)
{
}


int colvarbias_alb::init(std::string const &conf)
{
  colvarbias::init(conf);

  enable(f_cvb_scalar_variables);

  size_t i;

  // get the initial restraint centers
  colvar_centers.resize(num_variables());

  means.resize(num_variables());
  ssd.resize(num_variables()); //sum of squares of differences from mean

  //setup force vectors
  max_coupling_range.resize(num_variables());
  max_coupling_rate.resize(num_variables());
  coupling_accum.resize(num_variables());
  set_coupling.resize(num_variables());
  current_coupling.resize(num_variables());
  coupling_rate.resize(num_variables());

  enable(f_cvb_apply_force);

  for (i = 0; i < num_variables(); i++) {
    colvar_centers[i].type(colvars[i]->value());
    //zero moments
    means[i] = ssd[i] = 0;

    //zero force some of the force vectors that aren't initialized
    coupling_accum[i] = current_coupling[i] = 0;

  }
  if (get_keyval(conf, "centers", colvar_centers, colvar_centers)) {
    for (i = 0; i < num_variables(); i++) {
      colvar_centers[i].apply_constraints();
    }
  } else {
    colvar_centers.clear();
    cvm::fatal_error("Error: must define the initial centers of adaptive linear bias .\n");
  }

  if (colvar_centers.size() != num_variables())
    cvm::fatal_error("Error: number of centers does not match "
                      "that of collective variables.\n");

  if (!get_keyval(conf, "UpdateFrequency", update_freq, 0))
    cvm::fatal_error("Error: must set updateFrequency for adaptive linear bias.\n");

  //we split the time between updating and equilibrating
  update_freq /= 2;

  if (update_freq <= 1)
    cvm::fatal_error("Error: must set updateFrequency to greater than 2.\n");

  enable(f_cvb_history_dependent);

  get_keyval(conf, "outputCenters", b_output_centers, false);
  get_keyval(conf, "outputGradient", b_output_grad, false);
  get_keyval(conf, "outputCoupling", b_output_coupling, true);
  get_keyval(conf, "hardForceRange", b_hard_coupling_range, true);

  //initial guess
  if (!get_keyval(conf, "forceConstant", set_coupling, set_coupling))
    for (i =0 ; i < num_variables(); i++)
      set_coupling[i] = 0.;

  //how we're going to increase to that point
  for (i = 0; i < num_variables(); i++)
    coupling_rate[i] = (set_coupling[i] - current_coupling[i]) / update_freq;


  if (!get_keyval(conf, "forceRange", max_coupling_range, max_coupling_range)) {
    //set to default
    for (i = 0; i < num_variables(); i++) {
      if (cvm::temperature() > 0)
        max_coupling_range[i] =   3 * cvm::temperature() * cvm::boltzmann();
      else
        max_coupling_range[i] =   3 * cvm::boltzmann();
    }
  }

  if (!get_keyval(conf, "rateMax", max_coupling_rate, max_coupling_rate)) {
    //set to default
    for (i = 0; i < num_variables(); i++) {
      max_coupling_rate[i] =   max_coupling_range[i] / (10 * update_freq);
    }
  }


  if (cvm::debug())
    cvm::log(" bias.\n");

  return COLVARS_OK;
}


colvarbias_alb::~colvarbias_alb()
{
}


int colvarbias_alb::update()
{

  bias_energy = 0.0;
  update_calls++;

  if (cvm::debug())
    cvm::log("Updating the adaptive linear bias \""+this->name+"\".\n");

  //log the moments of the CVs
  // Force and energy calculation
  bool finished_equil_flag = 1;
  cvm::real delta;
  for (size_t i = 0; i < num_variables(); i++) {
    colvar_forces[i] = -1.0 * restraint_force(restraint_convert_k(current_coupling[i], colvars[i]->width),
                                              colvars[i],
                                              colvar_centers[i]);
    bias_energy += restraint_potential(restraint_convert_k(current_coupling[i], colvars[i]->width),
                                       colvars[i],
                                       colvar_centers[i]);

    if (!b_equilibration) {
      //Welford, West, and Hanso online variance method

      delta = static_cast<cvm::real>(colvars[i]->value())  - means[i];
      means[i] += delta / update_calls;
      ssd[i] += delta * (static_cast<cvm::real>(colvars[i]->value())  - means[i]);

    } else {
      //check if we've reached the setpoint
      cvm::real const coupling_diff = current_coupling[i] - set_coupling[i];
      if ((coupling_rate[i] == 0) ||
          ((coupling_diff*coupling_diff) < (coupling_rate[i]*coupling_rate[i]))) {
        finished_equil_flag &= 1; //we continue equilibrating as long as we haven't reached all the set points
      }
      else {
        current_coupling[i] += coupling_rate[i];
        finished_equil_flag = 0;
      }


      //update max_coupling_range
      if (!b_hard_coupling_range && fabs(current_coupling[i]) > max_coupling_range[i]) {
        std::ostringstream logStream;
        logStream << "Coupling constant for "
                  << colvars[i]->name
                  << " has exceeded coupling range of "
                  << max_coupling_range[i]
                  << ".\n";

        max_coupling_range[i] *= 1.25;
        logStream << "Expanding coupling range to "  << max_coupling_range[i] << ".\n";
        cvm::log(logStream.str());
      }


    }
  }

  if (b_equilibration && finished_equil_flag) {
    b_equilibration = false;
    update_calls = 0;
  }


  //now we update coupling constant, if necessary
  if (!b_equilibration && update_calls == update_freq) {

    //use estimated variance to take a step
    cvm::real step_size = 0;
    cvm::real temp;

    //reset means and sum of squares of differences
    for (size_t i = 0; i < num_variables(); i++) {

      temp = 2. * (means[i] / (static_cast<cvm::real> (colvar_centers[i])) - 1) * ssd[i] / (update_calls - 1);

      if (cvm::temperature() > 0)
        step_size = temp / (cvm::temperature()  * cvm::boltzmann());
      else
        step_size = temp / cvm::boltzmann();

      means[i] = 0;
      ssd[i] = 0;

      //stochastic if we do that update or not
      if (num_variables() == 1 || rand() < RAND_MAX / ((int) num_variables())) {
        coupling_accum[i] += step_size * step_size;
        current_coupling[i] = set_coupling[i];
        set_coupling[i] += max_coupling_range[i] / sqrt(coupling_accum[i]) * step_size;
        coupling_rate[i] = (set_coupling[i] - current_coupling[i]) / update_freq;
        //set to the minimum rate and then put the sign back on it
        coupling_rate[i] = copysign(fmin(fabs(coupling_rate[i]), max_coupling_rate[i]), coupling_rate[i]);
      } else {
        coupling_rate[i] = 0;
      }

    }

    update_calls = 0;
    b_equilibration = true;

  }

  return COLVARS_OK;
}


int colvarbias_alb::set_state_params(std::string const &conf)
{
  if (!get_keyval(conf, "setCoupling", set_coupling))
    cvm::fatal_error("Error: current setCoupling  is missing from the restart.\n");

  if (!get_keyval(conf, "currentCoupling", current_coupling))
    cvm::fatal_error("Error: current setCoupling  is missing from the restart.\n");

  if (!get_keyval(conf, "maxCouplingRange", max_coupling_range))
    cvm::fatal_error("Error: maxCouplingRange  is missing from the restart.\n");

  if (!get_keyval(conf, "couplingRate", coupling_rate))
    cvm::fatal_error("Error: current setCoupling  is missing from the restart.\n");

  if (!get_keyval(conf, "couplingAccum", coupling_accum))
    cvm::fatal_error("Error: couplingAccum is missing from the restart.\n");

  if (!get_keyval(conf, "mean", means))
    cvm::fatal_error("Error: current mean is missing from the restart.\n");

  if (!get_keyval(conf, "ssd", ssd))
    cvm::fatal_error("Error: current ssd is missing from the restart.\n");

  if (!get_keyval(conf, "updateCalls", update_calls))
    cvm::fatal_error("Error: current updateCalls is missing from the restart.\n");

  if (!get_keyval(conf, "b_equilibration", b_equilibration))
    cvm::fatal_error("Error: current updateCalls is missing from the restart.\n");

  return COLVARS_OK;
}


std::string const colvarbias_alb::get_state_params() const
{
  std::ostringstream os;
  os << "    setCoupling ";
  size_t i;
  for (i = 0; i < num_variables(); i++) {
    os << std::setprecision(cvm::en_prec)
       << std::setw(cvm::en_width) << set_coupling[i] << "\n";
  }
  os << "    currentCoupling ";
  for (i = 0; i < num_variables(); i++) {
    os << std::setprecision(cvm::en_prec)
       << std::setw(cvm::en_width) << current_coupling[i] << "\n";
  }
  os << "    maxCouplingRange ";
  for (i = 0; i < num_variables(); i++) {
    os << std::setprecision(cvm::en_prec)
       << std::setw(cvm::en_width) << max_coupling_range[i] << "\n";
  }
  os << "    couplingRate ";
  for (i = 0; i < num_variables(); i++) {
    os << std::setprecision(cvm::en_prec)
       << std::setw(cvm::en_width) << coupling_rate[i] << "\n";
  }
  os << "    couplingAccum ";
  for (i = 0; i < num_variables(); i++) {
    os << std::setprecision(cvm::en_prec)
       << std::setw(cvm::en_width) << coupling_accum[i] << "\n";
  }
  os << "    mean ";
  for (i = 0; i < num_variables(); i++) {
    os << std::setprecision(cvm::en_prec)
       << std::setw(cvm::en_width) << means[i] << "\n";
  }
  os << "    ssd ";
  for (i = 0; i < num_variables(); i++) {
    os << std::setprecision(cvm::en_prec)
       << std::setw(cvm::en_width) << ssd[i] << "\n";
  }
  os << "    updateCalls " << update_calls << "\n";
  if (b_equilibration)
    os << "    b_equilibration yes\n";
  else
    os << "    b_equilibration no\n";

  return os.str();
}


std::ostream & colvarbias_alb::write_traj_label(std::ostream &os)
{
  os << " ";

  if (b_output_energy)
    os << " E_"
       << cvm::wrap_string(this->name, cvm::en_width-2);

  if (b_output_coupling)
    for (size_t i = 0; i < current_coupling.size(); i++) {
      os << " ForceConst_" << i
         <<std::setw(cvm::en_width - 6 - (i / 10 + 1))
         << "";
    }

  if (b_output_grad)
    for (size_t i = 0; i < means.size(); i++) {
      os << "Grad_"
         << cvm::wrap_string(colvars[i]->name, cvm::cv_width - 4);
    }

  if (b_output_centers)
    for (size_t i = 0; i < num_variables(); i++) {
      size_t const this_cv_width = (colvars[i]->value()).output_width(cvm::cv_width);
      os << " x0_"
         << cvm::wrap_string(colvars[i]->name, this_cv_width-3);
    }

  return os;
}


std::ostream & colvarbias_alb::write_traj(std::ostream &os)
{
  os << " ";

  if (b_output_energy)
    os << " "
       << std::setprecision(cvm::en_prec) << std::setw(cvm::en_width)
       << bias_energy;

  if (b_output_coupling)
    for (size_t i = 0; i < current_coupling.size(); i++) {
      os << " "
         << std::setprecision(cvm::en_prec) << std::setw(cvm::en_width)
         << current_coupling[i];
    }


  if (b_output_centers)
    for (size_t i = 0; i < num_variables(); i++) {
      os << " "
         << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
         << colvar_centers[i];
    }

  if (b_output_grad)
    for (size_t i = 0; i < means.size(); i++) {
      os << " "
         << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
         << -2.0 * (means[i] / (static_cast<cvm::real>(colvar_centers[i])) - 1) * ssd[i] / (fmax(update_calls, 2.0) - 1);

    }

  return os;
}


cvm::real colvarbias_alb::restraint_potential(cvm::real k,
                                              colvar const *x,
                                              colvarvalue const &xcenter) const
{
  return k * (x->value() - xcenter);
}


colvarvalue colvarbias_alb::restraint_force(cvm::real k,
                                            colvar const * /* x */,
                                            colvarvalue const & /* xcenter */) const
{
  return k;
}


cvm::real colvarbias_alb::restraint_convert_k(cvm::real k,
                                              cvm::real dist_measure) const
{
  return k / dist_measure;
}

