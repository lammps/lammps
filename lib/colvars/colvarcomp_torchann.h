// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.
//
#ifndef COLVARCOMP_TORCH_H
#define COLVARCOMP_TORCH_H

// Declaration of torchann

#include <memory>

#include "colvar.h"
#include "colvarcomp.h"
#include "colvarmodule.h"

#ifdef COLVARS_TORCH

#include <torch/torch.h>
#include <torch/script.h>

class colvar::torchANN
  : public colvar::linearCombination
{
protected:
    torch::jit::script::Module nn;
    /// the index of nn output component
    size_t m_output_index = 0;
    bool use_double_input = false;
    //bool use_gpu;
    // 1d tensor, concatenation of values of sub-cvcs
    torch::Tensor input_tensor;
    torch::Tensor nn_outputs;
    torch::Tensor input_grad;
    // record the initial index of of sub-cvcs in input_tensor
    std::vector<int> cvc_indices;
public:
    torchANN();
    virtual ~torchANN();
    virtual int init(std::string const &conf);
    virtual void calc_value();
    virtual void calc_gradients();
    virtual void apply_force(colvarvalue const &force);
};

#else

class colvar::torchANN
  : public colvar::cvc
{
public:
    torchANN();
    virtual ~torchANN();
    virtual int init(std::string const &conf);
    virtual void calc_value();
};
#endif // COLVARS_TORCH checking

#endif
