#include "colvarmodule.h"
#include "colvarbias_alb.h"
#include "colvarbias.h"
#include <stdio.h>
#include <stdlib.h>

/* Note about nomenclature. Force constant is called a coupling
 * constant here to emphasize its changing in the code. Outwards,
 * everything is called a force constant to keep it consistent with
 * the rest of colvars.
 *
 */

colvarbias_alb::colvarbias_alb(std::string const &conf, char const *key) :
  colvarbias(conf, key), update_calls(0), b_equilibration(true) {
  size_t i;

  // get the initial restraint centers
  colvar_centers.resize (colvars.size());

  means.resize(colvars.size());
  ssd.resize(colvars.size()); //sum of squares of differences from mean

  //setup force vectors
  max_coupling_range.resize(colvars.size());
  max_coupling_rate.resize(colvars.size());
  coupling_accum.resize(colvars.size());
  set_coupling.resize(colvars.size());
  current_coupling.resize(colvars.size());
  coupling_rate.resize(colvars.size());


  for (i = 0; i < colvars.size(); i++) {
    colvar_centers[i].type (colvars[i]->type());
    //zero moments
    means[i] = ssd[i] = 0;

    //zero force some of the force vectors that aren't initialized
    coupling_accum[i] = current_coupling[i] = 0;

  }
  if (get_keyval (conf, "centers", colvar_centers, colvar_centers)) {
    for (i = 0; i < colvars.size(); i++) {
      colvar_centers[i].apply_constraints();
    }
  } else {
    colvar_centers.clear();
    cvm::fatal_error ("Error: must define the initial centers of adaptive linear bias .\n");
  }

  if (colvar_centers.size() != colvars.size())
    cvm::fatal_error ("Error: number of centers does not match "
                      "that of collective variables.\n");

  if(!get_keyval (conf, "UpdateFrequency", update_freq, 0))
    cvm::fatal_error("Error: must set updateFrequency for adaptive linear bias.\n");

  //we split the time between updating and equilibrating
  update_freq /= 2;

  if(update_freq == 0)
    cvm::fatal_error("Error: must set updateFrequency to greater than 2.\n");

  get_keyval (conf, "outputCenters", b_output_centers, false);
  get_keyval (conf, "outputGradient", b_output_grad, false);
  get_keyval (conf, "outputCoupling", b_output_coupling, true);
  get_keyval (conf, "hardForceRange", b_hard_coupling_range, true);

  //initial guess
  if(!get_keyval (conf, "forceConstant", set_coupling, set_coupling))
    for(i =0 ; i < colvars.size(); i++)
      set_coupling[i] = 0.;

  //how we're going to increase to that point
  for(i = 0; i < colvars.size(); i++)
    coupling_rate[i] = (set_coupling[i] - current_coupling[i]) / update_freq;


  if(!get_keyval (conf, "forceRange", max_coupling_range, max_coupling_range)) {
    //set to default
    for(i = 0; i < colvars.size(); i++) {
      if(cvm::temperature() > 0)
        max_coupling_range[i] =   3 * cvm::temperature() * cvm::boltzmann();
      else
        max_coupling_range[i] =   3 * cvm::boltzmann();
    }
  }



  if(!get_keyval (conf, "rateMax", max_coupling_rate, max_coupling_rate)) {
    //set to default
    for(i = 0; i < colvars.size(); i++) {
      max_coupling_rate[i] =   max_coupling_range[i] / (10 * update_freq);
    }
  }


  if (cvm::debug())
    cvm::log (" bias.\n");

}

colvarbias_alb::~colvarbias_alb() {

  if (cvm::n_rest_biases > 0)
    cvm::n_rest_biases -= 1;

}

cvm::real colvarbias_alb::update() {

  bias_energy = 0.0;
  update_calls++;

  if (cvm::debug())
    cvm::log ("Updating the adaptive linear bias \""+this->name+"\".\n");




  //log the moments of the CVs
  // Force and energy calculation
  bool finished_equil_flag = 1;
  cvm::real delta;
  for (size_t i = 0; i < colvars.size(); i++) {
    colvar_forces[i] = -1.0 * restraint_force(restraint_convert_k(current_coupling[i], colvars[i]->width),
                                              colvars[i],
                                              colvar_centers[i]);
    bias_energy += restraint_potential(restraint_convert_k(current_coupling[i], colvars[i]->width),
				       colvars[i],
				       colvar_centers[i]);

    if(!b_equilibration) {
      //Welford, West, and Hanso online variance method

      delta = static_cast<cvm::real>(colvars[i]->value())  - means[i];
      means[i] += delta / update_calls;
      ssd[i] += delta * (static_cast<cvm::real>(colvars[i]->value())  - means[i]);

    } else {
      //check if we've reached the setpoint
      if(coupling_rate[i] == 0 || pow(current_coupling[i] - set_coupling[i],2)  < pow(coupling_rate[i],2)) {
	finished_equil_flag &= 1; //we continue equilibrating as long as we haven't reached all the set points
      }
      else {
	current_coupling[i] += coupling_rate[i];
	finished_equil_flag = 0;
      }


      //update max_coupling_range
      if(!b_hard_coupling_range && fabs(current_coupling[i]) > max_coupling_range[i]) {
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

  if(b_equilibration && finished_equil_flag) {
    b_equilibration = false;
    update_calls = 0;
  }


  //now we update coupling constant, if necessary
  if(!b_equilibration && update_calls == update_freq) {

    //use estimated variance to take a step
    cvm::real step_size = 0;
    cvm::real temp;

    //reset means and sum of squares of differences
    for(size_t i = 0; i < colvars.size(); i++) {

      temp = 2. * (means[i] / (static_cast<cvm::real> (colvar_centers[i])) - 1) * ssd[i] / (update_calls - 1);

      if(cvm::temperature() > 0)
	step_size = temp / (cvm::temperature()  * cvm::boltzmann());
      else
	step_size = temp / cvm::boltzmann();

      means[i] = 0;
      ssd[i] = 0;

      //stochastic if we do that update or not
      if(colvars.size() == 1 || rand() < RAND_MAX / ((int) colvars.size())) {
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

  return bias_energy;

}


std::istream & colvarbias_alb::read_restart (std::istream &is)
{
  size_t const start_pos = is.tellg();

  cvm::log ("Restarting adaptive linear bias \""+
            this->name+"\".\n");

  std::string key, brace, conf;
  if ( !(is >> key)   || !(key == "ALB") ||
       !(is >> brace) || !(brace == "{") ||
       !(is >> colvarparse::read_block ("configuration", conf)) ) {

    cvm::log ("Error: in reading restart configuration for restraint bias \""+
              this->name+"\" at position "+
              cvm::to_str (is.tellg())+" in stream.\n");
    is.clear();
    is.seekg (start_pos, std::ios::beg);
    is.setstate (std::ios::failbit);
    return is;
  }

  std::string name = "";
  if ( (colvarparse::get_keyval (conf, "name", name, std::string (""), colvarparse::parse_silent)) &&
       (name != this->name) )
    cvm::fatal_error ("Error: in the restart file, the "
                      "\"ALB\" block has a wrong name\n");
  if (name.size() == 0) {
    cvm::fatal_error ("Error: \"ALB\" block in the restart file "
                      "has no identifiers.\n");
  }

  if (!get_keyval (conf, "setCoupling", set_coupling))
    cvm::fatal_error ("Error: current setCoupling  is missing from the restart.\n");

  if (!get_keyval (conf, "currentCoupling", current_coupling))
    cvm::fatal_error ("Error: current setCoupling  is missing from the restart.\n");

  if (!get_keyval (conf, "maxCouplingRange", max_coupling_range))
    cvm::fatal_error ("Error: maxCouplingRange  is missing from the restart.\n");


  if (!get_keyval (conf, "couplingRate", coupling_rate))
    cvm::fatal_error ("Error: current setCoupling  is missing from the restart.\n");

  if (!get_keyval (conf, "couplingAccum", coupling_accum))
    cvm::fatal_error ("Error: couplingAccum is missing from the restart.\n");


  if (!get_keyval (conf, "mean", means))
    cvm::fatal_error ("Error: current mean is missing from the restart.\n");

  if (!get_keyval (conf, "ssd", ssd))
    cvm::fatal_error ("Error: current ssd is missing from the restart.\n");

  if (!get_keyval (conf, "updateCalls", update_calls))
    cvm::fatal_error ("Error: current updateCalls is missing from the restart.\n");

  if (!get_keyval (conf, "b_equilibration", b_equilibration))
    cvm::fatal_error ("Error: current updateCalls is missing from the restart.\n");

  is >> brace;
  if (brace != "}") {
    cvm::fatal_error ("Error: corrupt restart information for adaptive linear bias \""+
                      this->name+"\": no matching brace at position "+
                      cvm::to_str (is.tellg())+" in the restart file.\n");
    is.setstate (std::ios::failbit);
  }

  return is;
}


std::ostream & colvarbias_alb::write_restart (std::ostream &os)
{
  os << "ALB {\n"
     << "  configuration {\n"
     << "    name " << this->name << "\n";
  os << "    setCoupling ";
  size_t i;
  for(i = 0; i < colvars.size(); i++) {
    os << std::setprecision (cvm::en_prec)
       << std::setw (cvm::en_width) << set_coupling[i] << "\n";
  }
  os << "    currentCoupling ";
  for(i = 0; i < colvars.size(); i++) {
    os << std::setprecision (cvm::en_prec)
       << std::setw (cvm::en_width) << current_coupling[i] << "\n";
  }
  os << "    maxCouplingRange ";
  for(i = 0; i < colvars.size(); i++) {
    os << std::setprecision (cvm::en_prec)
       << std::setw (cvm::en_width) << max_coupling_range[i] << "\n";
  }
  os << "    couplingRate ";
  for(i = 0; i < colvars.size(); i++) {
    os << std::setprecision (cvm::en_prec)
       << std::setw (cvm::en_width) << coupling_rate[i] << "\n";
  }
  os << "    couplingAccum ";
  for(i = 0; i < colvars.size(); i++) {
    os << std::setprecision (cvm::en_prec)
       << std::setw (cvm::en_width) << coupling_accum[i] << "\n";
  }
  os << "    mean ";
  for(i = 0; i < colvars.size(); i++) {
    os << std::setprecision (cvm::en_prec)
       << std::setw (cvm::en_width) << means[i] << "\n";
  }
  os << "    ssd ";
  for(i = 0; i < colvars.size(); i++) {
    os << std::setprecision (cvm::en_prec)
       << std::setw (cvm::en_width) << ssd[i] << "\n";
  }
  os << "    updateCalls " << update_calls << "\n";
  if(b_equilibration)
    os << "    b_equilibration yes\n";
  else
    os << "    b_equilibration no\n";

  os << "  }\n"
     << "}\n\n";

  return os;
}


std::ostream & colvarbias_alb::write_traj_label (std::ostream &os)
{
  os << " ";

  if (b_output_energy)
    os << " E_"
       << cvm::wrap_string (this->name, cvm::en_width-2);

  if (b_output_coupling)
    for(size_t i = 0; i < current_coupling.size(); i++) {
      os << " ForceConst_" << i
	 <<std::setw(cvm::en_width - 6 - (i / 10 + 1))
	 << "";
    }

  if(b_output_grad)
    for(size_t i = 0; i < means.size(); i++) {
      os << "Grad_"
	 << cvm::wrap_string(colvars[i]->name, cvm::cv_width - 4);
    }

  if (b_output_centers)
    for (size_t i = 0; i < colvars.size(); i++) {
      size_t const this_cv_width = (colvars[i]->value()).output_width (cvm::cv_width);
      os << " x0_"
         << cvm::wrap_string (colvars[i]->name, this_cv_width-3);
    }

  return os;
}


std::ostream & colvarbias_alb::write_traj (std::ostream &os)
{
  os << " ";

  if (b_output_energy)
    os << " "
       << std::setprecision (cvm::en_prec) << std::setw (cvm::en_width)
       << bias_energy;

  if(b_output_coupling)
    for(size_t i = 0; i < current_coupling.size(); i++) {
      os << " "
	 << std::setprecision (cvm::en_prec) << std::setw (cvm::en_width)
	 << current_coupling[i];
    }


  if (b_output_centers)
    for (size_t i = 0; i < colvars.size(); i++) {
      os << " "
         << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
         << colvar_centers[i];
    }

  if(b_output_grad)
    for(size_t i = 0; i < means.size(); i++) {
      os << " "
	 << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
	 << -2. * (means[i] / (static_cast<cvm::real> (colvar_centers[i])) - 1) * ssd[i] / (fmax(update_calls,2) - 1);

    }

  return os;
}


cvm::real colvarbias_alb::restraint_potential(cvm::real k,  const colvar* x,  const colvarvalue &xcenter) const
{
  return k * (x->value() - xcenter);
}

colvarvalue colvarbias_alb::restraint_force(cvm::real k,  const colvar* x,  const colvarvalue &xcenter) const
{
  return k;
}

cvm::real colvarbias_alb::restraint_convert_k(cvm::real k, cvm::real dist_measure) const
{
  return k / dist_measure;
}

