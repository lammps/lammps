// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvar.h"
#include "colvarcomp.h"
#include "colvarmodule.h"
#include "colvarparse.h"
#include "colvarvalue.h"

#include "colvarcomp_torchann.h"


#ifdef COLVARS_TORCH

colvar::torchANN::torchANN()
{
  set_function_type("torchANN");
  provide(f_cvc_periodic);
}

colvar::torchANN::~torchANN() {}


int colvar::torchANN::init(std::string const &conf) {

  int error_code = linearCombination::init(conf);

  std::string model_file ;
  get_keyval(conf, "modelFile", model_file, std::string(""));
  try {
    nn = torch::jit::load(model_file);
    nn.to(torch::kCPU);
    cvm::log("torch model loaded.") ;
  } catch (const std::exception & e) {
    return cvm::error("Error: couldn't load libtorch model (see below).\n" + cvm::to_str(e.what()),
                      COLVARS_INPUT_ERROR);
  }

  auto const legacy_keyword = get_keyval(conf, "m_output_index", m_output_index, m_output_index);
  if (legacy_keyword) {
    cvm::log("Warning: m_output_index is a deprecated keyword, please use output_component instead.\n");
  }
  get_keyval(conf, "output_component", m_output_index, m_output_index);

  get_keyval(conf, "doubleInputTensor", use_double_input, use_double_input);
  //get_keyval(conf, "useGPU", use_gpu, false);

  cvc_indices.resize(cv.size(),0);

  size_t num_inputs = 0;
  // compute total number of inputs of neural network
  for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv)
  {
    num_inputs += cv[i_cv]->value().size() ;
    if (i_cv < cv.size() - 1)
      cvc_indices[i_cv+1] = num_inputs;
  }
  cvm::log("Input dimension of model: " + cvm::to_str(num_inputs));

  // initialize the input tensor
  auto options = torch::TensorOptions().dtype(torch::kFloat32).requires_grad(true);

  /*
  if (use_gpu) {
    if (torch::cuda::is_available()) {
      try {
        nn.to(torch::kCUDA);
      } catch(const std::exception & e) {
        cvm::error("Failed to move model to GPU.");
        use_gpu = false;
      }
    } else {
      use_gpu = false;
      cvm::log("GPU not available.");
    }
  }

  if (use_gpu) {
    options = options.device(torch::kCUDA);
    if (use_double_input) {
      cvm::log("Data type reset to Float for GPU computation!");
      use_double_input = false;
    }
  }
  */

  if (use_double_input) {  // set type to double
    options = options.dtype(torch::kFloat64);
    nn.to(torch::kFloat64);
    cvm::log("Model's dtype: kFloat64.");
  } else {
    cvm::log("Model's dtype: kFloat32.");
  }

  input_tensor = torch::zeros({1,(long int) num_inputs}, options);

  try { // test the model
    std::vector<torch::jit::IValue> inputs={input_tensor};
    nn_outputs = nn.forward(inputs).toTensor()[0][m_output_index];
    cvm::log("Evaluating model with zero tensor succeeded.");
  } catch (const std::exception & e) {
    error_code |= cvm::error("Error: evaluating model with zero tensor failed (see below).\n" +
                                 cvm::to_str(e.what()),
                             COLVARS_INPUT_ERROR);
  }

  x.type(colvarvalue::type_scalar);

  return error_code;
}


void colvar::torchANN::calc_value() {

  for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv)
    cv[i_cv]->calc_value();

  /*
  if (use_gpu)
    input_tensor = input_tensor.to(torch::kCPU);
  */

  // set input tensor with no_grad
  {
    torch::NoGradGuard no_grad;
    size_t l = 0;
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
      const colvarvalue& current_cv_value = cv[i_cv]->value();
      if (current_cv_value.type() == colvarvalue::type_scalar) {
        input_tensor[0][l++] = cv[i_cv]->sup_coeff * (cvm::pow(current_cv_value.real_value, cv[i_cv]->sup_np));
      } else {
        for (size_t j_elem = 0; j_elem < current_cv_value.size(); ++j_elem)
          input_tensor[0][l++] = cv[i_cv]->sup_coeff * current_cv_value[j_elem];
      }
    }
  }

  /*
  if (use_gpu)
    input_tensor = input_tensor.to(torch::kCUDA);
  */

  std::vector<torch::jit::IValue> inputs={input_tensor};

  // evaluate the value of function
  nn_outputs = nn.forward(inputs).toTensor()[0][m_output_index];

  input_grad = torch::autograd::grad({nn_outputs}, {input_tensor})[0][0];

  /*
  if (use_gpu)
    input_grad = input_grad.to(torch::kCPU);
  */

  x = nn_outputs.item<double>() ;

  this->wrap(x);

}

void colvar::torchANN::calc_gradients() {
  for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
    cv[i_cv]->calc_gradients();
    if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
      const cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
      // get the initial index of this cvc
      size_t l = cvc_indices[i_cv];
      for (size_t j_elem = 0; j_elem < cv[i_cv]->value().size(); ++j_elem) {
        // get derivative of neural network wrt its input
        const cvm::real factor = input_grad[l+j_elem].item<double>();
        for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
          for (size_t l_atom = 0; l_atom < (cv[i_cv]->atom_groups)[k_ag]->size(); ++l_atom) {
            (cv[i_cv]->atom_groups)[k_ag]->grad_x(l_atom) *= factor_polynomial * factor;
            (cv[i_cv]->atom_groups)[k_ag]->grad_y(l_atom) *= factor_polynomial * factor;
            (cv[i_cv]->atom_groups)[k_ag]->grad_z(l_atom) *= factor_polynomial * factor;
          }
        }
      }
    }
  }
}

void colvar::torchANN::apply_force(colvarvalue const &force) {

  for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
    // If this CV uses explicit gradients, then atomic gradients is already calculated
    // We can apply the force to atom groups directly
    if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
      for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
        (cv[i_cv]->atom_groups)[k_ag]->apply_colvar_force(force.real_value);
      }
    } else {
      const colvarvalue& current_cv_value = cv[i_cv]->value();
      colvarvalue cv_force(current_cv_value);
      cv_force.reset();
      const cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
      // get the initial index of this cvc
      size_t l = cvc_indices[i_cv];
      for (size_t j_elem = 0; j_elem < current_cv_value.size(); ++j_elem) {
        cv_force[j_elem] = factor_polynomial * input_grad[l+j_elem].item<double>() * force.real_value;
      }
      cv[i_cv]->apply_force(cv_force);
    }
  }
}


#else

colvar::torchANN::torchANN()
{
  set_function_type("torchANN");
}

colvar::torchANN::~torchANN() {}

int colvar::torchANN::init(std::string const &conf) {

  return cvm::error(
          "torchANN requires the libtorch library, but it is not enabled during compilation.\n"
          "Please refer to the Compilation Notes section of the Colvars manual for more "
          "information.\n",
          COLVARS_NOT_IMPLEMENTED);

}

void colvar::torchANN::calc_value()
{
}

#endif
