#if (__cplusplus >= 201103L)

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvarcomp.h"
#include "colvar_neuralnetworkcompute.h"

using namespace neuralnetworkCV;

colvar::neuralNetwork::neuralNetwork(std::string const &conf): linearCombination(conf) {
    set_function_type("neuralNetwork");
    // the output of neural network consists of multiple values
    // read "output_component" key to determine it
    get_keyval(conf, "output_component", m_output_index);
    // read weight files
    bool has_weight_files = true;
    size_t num_layers_weight = 0;
    std::vector<std::string> weight_files;
    while (has_weight_files) {
        std::string lookup_key = std::string{"layer"} + cvm::to_str(num_layers_weight + 1) + std::string{"_WeightsFile"};
        if (key_lookup(conf, lookup_key.c_str())) {
            std::string weight_filename;
            get_keyval(conf, lookup_key.c_str(), weight_filename, std::string(""));
            weight_files.push_back(weight_filename);
            cvm::log(std::string{"Will read layer["} + cvm::to_str(num_layers_weight + 1) + std::string{"] weights from "} + weight_filename + '\n');
            ++num_layers_weight;
        } else {
            has_weight_files = false;
        }
    }
    // read bias files
    bool has_bias_files = true;
    size_t num_layers_bias = 0;
    std::vector<std::string> bias_files;
    while (has_bias_files) {
        std::string lookup_key = std::string{"layer"} + cvm::to_str(num_layers_bias + 1) + std::string{"_BiasesFile"};
        if (key_lookup(conf, lookup_key.c_str())) {
            std::string bias_filename;
            get_keyval(conf, lookup_key.c_str(), bias_filename, std::string(""));
            bias_files.push_back(bias_filename);
            cvm::log(std::string{"Will read layer["} + cvm::to_str(num_layers_bias + 1) + std::string{"] biases from "} + bias_filename + '\n');
            ++num_layers_bias;
        } else {
            has_bias_files = false;
        }
    }
    // read activation function strings
    bool has_activation_functions = true;
    size_t num_activation_functions = 0;
    // pair(is_custom_function, function_string)
    std::vector<std::pair<bool, std::string>> activation_functions;
    while (has_activation_functions) {
        std::string lookup_key = std::string{"layer"} + cvm::to_str(num_activation_functions + 1) + std::string{"_activation"};
        std::string lookup_key_custom = std::string{"layer"} + cvm::to_str(num_activation_functions + 1) + std::string{"_custom_activation"};
        if (key_lookup(conf, lookup_key.c_str())) {
            // Ok, this is not a custom function
            std::string function_name;
            get_keyval(conf, lookup_key.c_str(), function_name, std::string(""));
            if (activation_function_map.find(function_name) == activation_function_map.end()) {
                cvm::error("Unknown activation function name: \"" + function_name + "\".\n");
                return;
            }
            activation_functions.push_back(std::make_pair(false, function_name));
            cvm::log(std::string{"The activation function for layer["} + cvm::to_str(num_activation_functions + 1) + std::string{"] is "} + function_name + '\n');
            ++num_activation_functions;
#ifdef LEPTON
        } else if (key_lookup(conf, lookup_key_custom.c_str())) {
            std::string function_expression;
            get_keyval(conf, lookup_key_custom.c_str(), function_expression, std::string(""));
            activation_functions.push_back(std::make_pair(true, function_expression));
            cvm::log(std::string{"The custom activation function for layer["} + cvm::to_str(num_activation_functions + 1) + std::string{"] is "} + function_expression + '\n');
            ++num_activation_functions;
#endif
        } else {
            has_activation_functions = false;
        }
    }
    // expect the three numbers are equal
    if ((num_layers_weight != num_layers_bias) || (num_layers_bias != num_activation_functions)) {
        cvm::error("Error: the numbers of weights, biases and activation functions do not match.\n");
        return;
    }
//     nn = std::make_unique<neuralnetworkCV::neuralNetworkCompute>();
    // std::make_unique is only available in C++14
    nn = std::unique_ptr<neuralnetworkCV::neuralNetworkCompute>(new neuralnetworkCV::neuralNetworkCompute());
    for (size_t i_layer = 0; i_layer < num_layers_weight; ++i_layer) {
        denseLayer d;
#ifdef LEPTON
        if (activation_functions[i_layer].first) {
            // use custom function as activation function
            try {
                d = denseLayer(weight_files[i_layer], bias_files[i_layer], activation_functions[i_layer].second);
            } catch (std::exception &ex) {
                cvm::error("Error on initializing layer " + cvm::to_str(i_layer) + " (" + ex.what() + ")\n", COLVARS_INPUT_ERROR);
                return;
            }
        } else {
#endif
            // query the map of supported activation functions
            const auto& f = activation_function_map[activation_functions[i_layer].second].first;
            const auto& df = activation_function_map[activation_functions[i_layer].second].second;
            try {
                d = denseLayer(weight_files[i_layer], bias_files[i_layer], f, df);
            } catch (std::exception &ex) {
                cvm::error("Error on initializing layer " + cvm::to_str(i_layer) + " (" + ex.what() + ")\n", COLVARS_INPUT_ERROR);
                return;
            }
#ifdef LEPTON
        }
#endif
        // add a new dense layer to network
        if (nn->addDenseLayer(d)) {
            if (cvm::debug()) {
                // show information about the neural network
                cvm::log("Layer " + cvm::to_str(i_layer) + " : has " + cvm::to_str(d.getInputSize()) + " input nodes and " + cvm::to_str(d.getOutputSize()) + " output nodes.\n");
                for (size_t i_output = 0; i_output < d.getOutputSize(); ++i_output) {
                    for (size_t j_input = 0; j_input < d.getInputSize(); ++j_input) {
                        cvm::log("    weights[" + cvm::to_str(i_output) + "][" + cvm::to_str(j_input) + "] = " + cvm::to_str(d.getWeight(i_output, j_input)));
                    }
                    cvm::log("    biases[" + cvm::to_str(i_output) + "] = " + cvm::to_str(d.getBias(i_output)) + "\n");
                }
            }
        } else {
            cvm::error("Error: error on adding a new dense layer.\n");
            return;
        }
    }
    nn->input().resize(cv.size());
}

colvar::neuralNetwork::~neuralNetwork() {
}

void colvar::neuralNetwork::calc_value() {
    x.reset();
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_value();
        const colvarvalue& current_cv_value = cv[i_cv]->value();
        // for current nn implementation we have to assume taht types are always scaler
        if (current_cv_value.type() == colvarvalue::type_scalar) {
            nn->input()[i_cv] = cv[i_cv]->sup_coeff * (cvm::pow(current_cv_value.real_value, cv[i_cv]->sup_np));
        } else {
            cvm::error("Error: using of non-scaler component.\n");
            return;
        }
    }
    nn->compute();
    x = nn->getOutput(m_output_index);
}

void colvar::neuralNetwork::calc_gradients() {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_gradients();
        if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            const cvm::real factor = nn->getGradient(m_output_index, i_cv);
            const cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            for (size_t j_elem = 0; j_elem < cv[i_cv]->value().size(); ++j_elem) {
                for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                    for (size_t l_atom = 0; l_atom < (cv[i_cv]->atom_groups)[k_ag]->size(); ++l_atom) {
                        (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad = factor_polynomial * factor * (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad;
                    }
                }
            }
        }
    }
}

void colvar::neuralNetwork::apply_force(colvarvalue const &force) {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        // If this CV us explicit gradients, then atomic gradients is already calculated
        // We can apply the force to atom groups directly
        if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                (cv[i_cv]->atom_groups)[k_ag]->apply_colvar_force(force.real_value);
            }
        } else {
            // Compute factors for polynomial combinations
            const cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            const cvm::real factor = nn->getGradient(m_output_index, i_cv);;
            colvarvalue cv_force = force.real_value * factor * factor_polynomial;
            cv[i_cv]->apply_force(cv_force);
        }
    }
}

#endif
