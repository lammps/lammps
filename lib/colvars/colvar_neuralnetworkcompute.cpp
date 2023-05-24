// -*- Mode:c++; c-basic-offset: 4; -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <iostream>
#include <fstream>

#if (__cplusplus >= 201103L)
#include "colvar_neuralnetworkcompute.h"
#include "colvarparse.h"
#include "colvarproxy.h"

namespace neuralnetworkCV {
std::map<std::string, std::pair<std::function<double(double)>, std::function<double(double)>>> activation_function_map
{
    {"tanh",     {[](double x){return std::tanh(x);},
                  [](double x){return 1.0 - std::tanh(x) * std::tanh(x);}}},
    {"sigmoid",  {[](double x){return 1.0 / (1.0 + std::exp(-x));},
                  [](double x){return std::exp(-x) / ((1.0 + std::exp(-x)) * (1.0 + std::exp(-x)));}}},
    {"linear",   {[](double x){return x;},
                  [](double /*x*/){return 1.0;}}},
    {"relu",     {[](double x){return x < 0. ? 0. : x;},
                  [](double x){return x < 0. ? 0. : 1.;}}},
    {"lrelu100", {[](double x){return x < 0. ? 0.01 * x : x;},
                  [](double x){return x < 0. ? 0.01     : 1.;}}},
    {"elu",      {[](double x){return x < 0. ? std::exp(x)-1. : x;},
                  [](double x){return x < 0. ? std::exp(x)    : 1.;}}}
};

#ifdef LEPTON
customActivationFunction::customActivationFunction():
expression(), value_evaluator(nullptr), gradient_evaluator(nullptr),
input_reference(nullptr), derivative_reference(nullptr) {}

customActivationFunction::customActivationFunction(const std::string& expression_string):
expression(), value_evaluator(nullptr), gradient_evaluator(nullptr),
input_reference(nullptr), derivative_reference(nullptr) {
    setExpression(expression_string);
}

customActivationFunction::customActivationFunction(const customActivationFunction& source):
expression(), value_evaluator(nullptr), gradient_evaluator(nullptr),
input_reference(nullptr), derivative_reference(nullptr) {
    // check if the source object is initialized
    if (source.value_evaluator != nullptr) {
        this->setExpression(source.expression);
    }
}

customActivationFunction& customActivationFunction::operator=(const customActivationFunction& source) {
    if (source.value_evaluator != nullptr) {
        this->setExpression(source.expression);
    } else {
        expression = std::string();
        value_evaluator = nullptr;
        gradient_evaluator = nullptr;
        input_reference = nullptr;
        derivative_reference = nullptr;
    }
    return *this;
}

void customActivationFunction::setExpression(const std::string& expression_string) {
    expression = expression_string;
    Lepton::ParsedExpression parsed_expression;
    // the variable must be "x" for the input of an activation function
    const std::string activation_input_variable{"x"};
    // parse the expression
    try {
        parsed_expression = Lepton::Parser::parse(expression);
    } catch (...) {
        cvm::error("Error parsing or compiling expression \"" + expression + "\".\n", COLVARS_INPUT_ERROR);
    }
    // compile the expression
    try {
        value_evaluator = std::unique_ptr<Lepton::CompiledExpression>(new Lepton::CompiledExpression(parsed_expression.createCompiledExpression()));
    } catch (...) {
        cvm::error("Error compiling expression \"" + expression + "\".\n", COLVARS_INPUT_ERROR);
    }
    // create a compiled expression for the derivative
    try {
        gradient_evaluator = std::unique_ptr<Lepton::CompiledExpression>(new Lepton::CompiledExpression(parsed_expression.differentiate(activation_input_variable).createCompiledExpression()));
    } catch (...) {
        cvm::error("Error creating compiled expression for variable \"" + activation_input_variable + "\".\n", COLVARS_INPUT_ERROR);
    }
    // get the reference to the input variable in the compiled expression
    try {
        input_reference = &(value_evaluator->getVariableReference(activation_input_variable));
    } catch (...) {
        cvm::error("Error on getting the reference to variable \"" + activation_input_variable + "\" in the compiled expression.\n", COLVARS_INPUT_ERROR);
    }
    // get the reference to the input variable in the compiled derivative expression
    try {
        derivative_reference = &(gradient_evaluator->getVariableReference(activation_input_variable));
    } catch (...) {
        cvm::error("Error on getting the reference to variable \"" + activation_input_variable + "\" in the compiled derivative exprssion.\n", COLVARS_INPUT_ERROR);
    }
}

std::string customActivationFunction::getExpression() const {
    return expression;
}

double customActivationFunction::evaluate(double x) const {
    *input_reference = x;
    return value_evaluator->evaluate();
}

double customActivationFunction::derivative(double x) const {
    *derivative_reference = x;
    return gradient_evaluator->evaluate();
}
#endif

denseLayer::denseLayer(const std::string& weights_file, const std::string& biases_file, const std::function<double(double)>& f, const std::function<double(double)>& df): m_activation_function(f), m_activation_function_derivative(df) {
#ifdef LEPTON
    m_use_custom_activation = false;
#endif
    readFromFile(weights_file, biases_file);
}

#ifdef LEPTON
denseLayer::denseLayer(const std::string& weights_file, const std::string& biases_file, const std::string& custom_activation_expression) {
    m_use_custom_activation = true;
    m_custom_activation_function = customActivationFunction(custom_activation_expression);
    readFromFile(weights_file, biases_file);
}
#endif

void denseLayer::readFromFile(const std::string& weights_file, const std::string& biases_file) {
    // parse weights file
    m_weights.clear();
    m_biases.clear();
    std::string line;
    colvarproxy *proxy = cvm::main()->proxy;
    auto &ifs_weights = proxy->input_stream(weights_file, "weights file");
    while (std::getline(ifs_weights, line)) {
        if (!ifs_weights) {
            throw std::runtime_error("I/O error while reading " + weights_file);
        }
        std::vector<std::string> splitted_data;
        colvarparse::split_string(line, std::string{" "}, splitted_data);
        if (splitted_data.size() > 0) {
            std::vector<double> weights_tmp(splitted_data.size());
            for (size_t i = 0; i < splitted_data.size(); ++i) {
                try {
                    weights_tmp[i] = std::stod(splitted_data[i]);
                } catch (...) {
                    throw std::runtime_error("Cannot convert " + splitted_data[i] + " to a number while reading file " + weights_file);
                }
            }
            m_weights.push_back(weights_tmp);
        }
    }
    proxy->close_input_stream(weights_file);

    // parse biases file
    auto &ifs_biases = proxy->input_stream(biases_file, "biases file");
    while (std::getline(ifs_biases, line)) {
        if (!ifs_biases) {
            throw std::runtime_error("I/O error while reading " + biases_file);
        }
        std::vector<std::string> splitted_data;
        colvarparse::split_string(line, std::string{" "}, splitted_data);
        if (splitted_data.size() > 0) {
            double bias = 0;
            try {
                bias = std::stod(splitted_data[0]);
            } catch (...) {
                throw std::runtime_error("Cannot convert " + splitted_data[0] + " to a number while reading file " + biases_file);
            }
            m_biases.push_back(bias);
        }
    }
    proxy->close_input_stream(biases_file);

    m_input_size = m_weights[0].size();
    m_output_size = m_weights.size();
}

void denseLayer::setActivationFunction(const std::function<double(double)>& f, const std::function<double(double)>& df) {
    m_activation_function = f;
    m_activation_function_derivative = df;
}

void denseLayer::compute(const std::vector<double>& input, std::vector<double>& output) const {
    for (size_t i = 0; i < m_output_size; ++i) {
        output[i] = 0;
        for (size_t j = 0; j < m_input_size; ++j) {
            output[i] += input[j] * m_weights[i][j];
        }
        output[i] += m_biases[i];
#ifdef LEPTON
        if (m_use_custom_activation) {
            output[i] = m_custom_activation_function.evaluate(output[i]);
        } else {
#endif
            output[i] = m_activation_function(output[i]);
#ifdef LEPTON
        }
#endif
    }
}

double denseLayer::computeGradientElement(const std::vector<double>& input, const size_t i, const size_t j) const {
    double sum_with_bias = 0;
    for (size_t j_in = 0; j_in < m_input_size; ++j_in) {
        sum_with_bias += input[j_in] * m_weights[i][j_in];
    }
    sum_with_bias += m_biases[i];
#ifdef LEPTON
    if (m_use_custom_activation) {
        const double grad_ij = m_custom_activation_function.derivative(sum_with_bias) * m_weights[i][j];
        return grad_ij;
    } else {
#endif
        const double grad_ij = m_activation_function_derivative(sum_with_bias) * m_weights[i][j];
        return grad_ij;
#ifdef LEPTON
    }
#endif
}

void denseLayer::computeGradient(const std::vector<double>& input, std::vector<std::vector<double>>& output_grad) const {
    for (size_t j = 0; j < m_input_size; ++j) {
        for (size_t i = 0; i < m_output_size; ++i) {
            output_grad[i][j] = computeGradientElement(input, i, j);
        }
    }
}

neuralNetworkCompute::neuralNetworkCompute(const std::vector<denseLayer>& dense_layers): m_dense_layers(dense_layers) {
    m_layers_output.resize(m_dense_layers.size());
    m_grads_tmp.resize(m_dense_layers.size());
    for (size_t i_layer = 0; i_layer < m_layers_output.size(); ++i_layer) {
        m_layers_output[i_layer].assign(m_dense_layers[i_layer].getOutputSize(), 0);
        m_grads_tmp[i_layer].assign(m_dense_layers[i_layer].getOutputSize(), std::vector<double>(m_dense_layers[i_layer].getInputSize(), 0));
    }
}

bool neuralNetworkCompute::addDenseLayer(const denseLayer& layer) {
    if (m_dense_layers.empty()) {
        // add layer to this ann directly if m_dense_layers is empty
        m_dense_layers.push_back(layer);
        m_layers_output.push_back(std::vector<double>(layer.getOutputSize()));
        m_grads_tmp.push_back(std::vector<std::vector<double>>(layer.getOutputSize(), std::vector<double>(layer.getInputSize(), 0)));
        return true;
    } else {
        // otherwise, we need to check if the output of last layer in m_dense_layers matches the input of layer to be added
        if (m_dense_layers.back().getOutputSize() == layer.getInputSize()) {
            m_dense_layers.push_back(layer);
            m_layers_output.push_back(std::vector<double>(layer.getOutputSize()));
            m_grads_tmp.push_back(std::vector<std::vector<double>>(layer.getOutputSize(), std::vector<double>(layer.getInputSize(), 0)));
            return true;
        } else {
            return false;
        }
    }
}

std::vector<std::vector<double>> neuralNetworkCompute::multiply_matrix(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
    const size_t m = A.size();
    const size_t n = B.size();
    if (A[0].size() != n) {
        std::cerr << "Error on multiplying matrices!\n";
    }
    const size_t t = B[0].size();
    std::vector<std::vector<double>> C(m, std::vector<double>(t, 0.0));
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < t; ++j) {
            for (size_t k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

void neuralNetworkCompute::compute() {
    if (m_dense_layers.empty()) {
        return;
    }
    size_t i_layer;
    m_dense_layers[0].compute(m_input, m_layers_output[0]);
    for (i_layer = 1; i_layer < m_dense_layers.size(); ++i_layer) {
        m_dense_layers[i_layer].compute(m_layers_output[i_layer - 1], m_layers_output[i_layer]);
    }
    // gradients of each layer
    m_dense_layers[0].computeGradient(m_input, m_grads_tmp[0]);
    for (i_layer = 1; i_layer < m_dense_layers.size(); ++i_layer) {
        m_dense_layers[i_layer].computeGradient(m_layers_output[i_layer - 1], m_grads_tmp[i_layer]);
    }
    // chain rule
    if (m_dense_layers.size() > 1) {
        m_chained_grad = multiply_matrix(m_grads_tmp[1], m_grads_tmp[0]);
        for (i_layer = 2; i_layer < m_dense_layers.size(); ++i_layer) {
            m_chained_grad = multiply_matrix(m_grads_tmp[i_layer], m_chained_grad);
        }
    } else {
        m_chained_grad = m_grads_tmp[0];
    }
}
}

#endif
