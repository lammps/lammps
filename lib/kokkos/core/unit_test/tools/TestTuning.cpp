//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

// This file tests the primitives of the Tuning system

#include <iostream>
#include <Kokkos_Core.hpp>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

static size_t expectedNumberOfContextVariables;
static int64_t expectedContextVariableValue;

int main() {
  Kokkos::initialize();
  {
    auto context = Kokkos::Tools::Experimental::get_new_context_id();

    Kokkos::Tools::Experimental::VariableInfo contextVariableInfo;

    contextVariableInfo.category = Kokkos::Tools::Experimental::
        StatisticalCategory::kokkos_value_categorical;
    contextVariableInfo.type =
        Kokkos::Tools::Experimental::ValueType::kokkos_value_int64;
    contextVariableInfo.valueQuantity =
        Kokkos::Tools::Experimental::CandidateValueType::kokkos_value_unbounded;

    Kokkos::Tools::Experimental::VariableInfo tuningVariableInfo;

    tuningVariableInfo.category = Kokkos::Tools::Experimental::
        StatisticalCategory::kokkos_value_categorical;
    tuningVariableInfo.type =
        Kokkos::Tools::Experimental::ValueType::kokkos_value_int64;
    tuningVariableInfo.valueQuantity =
        Kokkos::Tools::Experimental::CandidateValueType::kokkos_value_set;

    std::vector<int64_t> candidate_value_vector = {0, 1, 2, 3, 4,
                                                   5, 6, 7, 8, 9};

    Kokkos::Tools::Experimental::SetOrRange allowed_values =
        Kokkos::Tools::Experimental::make_candidate_set(
            candidate_value_vector.size(), candidate_value_vector.data());
    // test that ID's are transmitted to the tool
    Kokkos::Tools::Experimental::set_declare_output_type_callback(
        [](const char*, const size_t,
           Kokkos::Tools::Experimental::VariableInfo* info) {
          if (info->type !=
              Kokkos::Tools::Experimental::ValueType::kokkos_value_int64) {
            throw(std::runtime_error("Tuning Variable has wrong type"));
          }
        });
    Kokkos::Tools::Experimental::set_declare_input_type_callback(
        [](const char*, const size_t,
           Kokkos::Tools::Experimental::VariableInfo* info) {
          if (info->type !=
              Kokkos::Tools::Experimental::ValueType::kokkos_value_int64) {
            throw(std::runtime_error("Context Variable has wrong type"));
          }
        });
    tuningVariableInfo.candidates = allowed_values;
    auto contextVariableId = Kokkos::Tools::Experimental::declare_input_type(
        "kokkos.testing.context_variable", contextVariableInfo);
    auto tuningVariableId = Kokkos::Tools::Experimental::declare_output_type(
        "kokkos.testing.tuning_variable", tuningVariableInfo);

    // test that we correctly pass context values, and receive tuning variables
    // back in return
    Kokkos::Tools::Experimental::VariableValue contextValues[] = {
        Kokkos::Tools::Experimental::make_variable_value(contextVariableId,
                                                         int64_t(0))};
    Kokkos::Tools::Experimental::set_input_values(context, 1, contextValues);

    Kokkos::Tools::Experimental::set_request_output_values_callback(
        [](const size_t, const size_t,
           const Kokkos::Tools::Experimental::VariableValue* context_values,
           const size_t,
           Kokkos::Tools::Experimental::VariableValue* tuning_values) {
          auto candidate_values = tuning_values[0].metadata->candidates;
          if (context_values[0].value.int_value !=
              expectedContextVariableValue) {
            throw std::runtime_error(
                "Context variables not correctly passed to tuning callbacks");
          }
          int tuningVariableSetSize = candidate_values.set.size;
          std::cout << "Set of size " << tuningVariableSetSize << std::endl;
          // tuning methodology via https://xkcd.com/221/
          tuning_values[0].value.int_value =
              candidate_values.set.values.int_value[4 % tuningVariableSetSize];
        });

    Kokkos::Tools::Experimental::VariableValue tuningValues[] = {
        Kokkos::Tools::Experimental::make_variable_value(tuningVariableId,
                                                         int64_t(0))};

    Kokkos::Tools::Experimental::request_output_values(context, 1,
                                                       tuningValues);
    std::cout << tuningValues[0].value.int_value << ","
              << candidate_value_vector[4] << std::endl;
    if (tuningValues[0].value.int_value != candidate_value_vector[4]) {
      throw std::runtime_error("Tuning value return is incorrect");
    }

    Kokkos::Tools::Experimental::end_context(context);

    // test nested contexts
    auto outerContext = Kokkos::Tools::Experimental::get_new_context_id();
    auto innerContext = Kokkos::Tools::Experimental::get_new_context_id();

    Kokkos::Tools::Experimental::VariableInfo secondContextVariableInfo;

    secondContextVariableInfo.category = Kokkos::Tools::Experimental::
        StatisticalCategory::kokkos_value_categorical;
    secondContextVariableInfo.type =
        Kokkos::Tools::Experimental::ValueType::kokkos_value_int64;
    secondContextVariableInfo.valueQuantity =
        Kokkos::Tools::Experimental::CandidateValueType::kokkos_value_unbounded;
    auto secondContextVariableId =
        Kokkos::Tools::Experimental::declare_output_type(
            "kokkos.testing.second_context_variable",
            secondContextVariableInfo);

    Kokkos::Tools::Experimental::VariableValue contextValueTwo[] = {
        Kokkos::Tools::Experimental::make_variable_value(
            secondContextVariableId, int64_t(1))};

    Kokkos::Tools::Experimental::set_request_output_values_callback(
        [](const size_t, const size_t num_context_variables,
           const Kokkos::Tools::Experimental::VariableValue*, const size_t,
           Kokkos::Tools::Experimental::VariableValue*) {
          std::cout << "Expect " << expectedNumberOfContextVariables
                    << ", have " << num_context_variables << std::endl;
          if (num_context_variables != expectedNumberOfContextVariables) {
            throw(
                std::runtime_error("Incorrect number of context variables in "
                                   "nested tuning contexts"));
          }
        });
    Kokkos::Tools::Experimental::set_input_values(outerContext, 1,
                                                  contextValues);
    expectedNumberOfContextVariables = 1;
    Kokkos::Tools::Experimental::request_output_values(outerContext, 1,
                                                       tuningValues);
    Kokkos::Tools::Experimental::set_input_values(innerContext, 1,
                                                  contextValueTwo);
    expectedNumberOfContextVariables = 2;
    Kokkos::Tools::Experimental::request_output_values(innerContext, 1,
                                                       tuningValues);
  }  // end Kokkos block

  Kokkos::finalize();
}
