// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#if (__cplusplus >= 201103L)
#ifndef NEURALNETWORKCOMPUTE_H
#define NEURALNETWORKCOMPUTE_H

#include <vector>
#include <functional>
#include <string>
#include <cmath>
#include <memory>
#include <map>

#ifdef LEPTON
#include "Lepton.h"
#endif

namespace neuralnetworkCV {
/// mapping from a string to the activation function and its derivative
extern std::map<std::string, std::pair<std::function<double(double)>, std::function<double(double)>>> activation_function_map;

#ifdef LEPTON
// allow to define a custom activation function
class customActivationFunction {
public:
    /// empty constructor
    customActivationFunction();
    /// construct by an mathematical expression
    customActivationFunction(const std::string& expression_string);
    /// copy constructor
    customActivationFunction(const customActivationFunction& source);
    /// overload assignment operator
    customActivationFunction& operator=(const customActivationFunction& source);
    /// setter for the custom expression
    void setExpression(const std::string& expression_string);
    /// getter for the custom expression
    std::string getExpression() const;
    /// evaluate the value of an expression
    double evaluate(double x) const;
    /// evaluate the gradient of an expression
    double derivative(double x) const;
private:
    std::string expression;
    std::unique_ptr<Lepton::CompiledExpression> value_evaluator;
    std::unique_ptr<Lepton::CompiledExpression> gradient_evaluator;
    double* input_reference;
    double* derivative_reference;
};
#endif

class denseLayer {
private:
    size_t m_input_size;
    size_t m_output_size;
    std::function<double(double)> m_activation_function;
    std::function<double(double)> m_activation_function_derivative;
#ifdef LEPTON
    bool m_use_custom_activation;
    customActivationFunction m_custom_activation_function;
#else
    static const bool m_use_custom_activation = false;
#endif
    /// weights[i][j] is the weight of the i-th output and the j-th input
    std::vector<std::vector<double>> m_weights;
    /// bias of each node
    std::vector<double> m_biases;
public:
    /// empty constructor
    denseLayer() {}
    /*! @param[in]  weights_file    filename of the weights file
     *  @param[in]  biases_file     filename of the biases file
     *  @param[in]  f               activation function
     *  @param[in]  df              derivative of the activation function
     */
    denseLayer(const std::string& weights_file, const std::string& biases_file, const std::function<double(double)>& f, const std::function<double(double)>& df);
#ifdef LEPTON
    /*! @param[in]  weights_file                 filename of the weights file
     *  @param[in]  biases_file                  filename of the biases file
     *  @param[in]  custom_activation_expression the expression of the custom activation function
     */
    denseLayer(const std::string& weights_file, const std::string& biases_file, const std::string& custom_activation_expression);
#endif
    /// read data from file
    void readFromFile(const std::string& weights_file, const std::string& biases_file);
    /// setup activation function
    void setActivationFunction(const std::function<double(double)>& f, const std::function<double(double)>& df);
    /// compute the value of this layer
    void compute(const std::vector<double>& input, std::vector<double>& output) const;
    /// compute the gradient of i-th output wrt j-th input
    double computeGradientElement(const std::vector<double>& input, const size_t i, const size_t j) const;
    /// output[i][j] is the gradient of i-th output wrt j-th input
    void computeGradient(const std::vector<double>& input, std::vector<std::vector<double>>& output_grad) const;
    /// get the input size
    size_t getInputSize() const {
        return m_input_size;
    }
    /// get the output size
    size_t getOutputSize() const {
        return m_output_size;
    }
    /// getter for weights and biases
    double getWeight(size_t i, size_t j) const {
        return m_weights[i][j];
    }
    double getBias(size_t i) const {
        return m_biases[i];
    }
    ~denseLayer() {}
};

class neuralNetworkCompute {
private:
    std::vector<denseLayer> m_dense_layers;
    std::vector<double> m_input;
    /// temporary output for each layer, useful to speedup the gradients' calculation
    std::vector<std::vector<double>> m_layers_output;
    std::vector<std::vector<std::vector<double>>> m_grads_tmp;
    std::vector<std::vector<double>> m_chained_grad;
private:
    /// helper function: multiply two matrix constructed from 2D vector
    static std::vector<std::vector<double>> multiply_matrix(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B);
public:
    neuralNetworkCompute(): m_dense_layers(0), m_layers_output(0) {}
    neuralNetworkCompute(const std::vector<denseLayer>& dense_layers);
    bool addDenseLayer(const denseLayer& layer);
    // for faster computation
    const std::vector<double>& input() const {return m_input;}
    std::vector<double>& input() {return m_input;}
    /// compute the values and the gradients of all output nodes
    void compute();
    double getOutput(const size_t i) const {return m_layers_output.back()[i];}
    double getGradient(const size_t i, const size_t j) const {return m_chained_grad[i][j];}
    /// get a specified layer
    const denseLayer& getLayer(const size_t i) const {return m_dense_layers[i];}
    /// get the number of layers
    size_t getNumberOfLayers() const {return m_dense_layers.size();}
};

}
#endif
#endif
