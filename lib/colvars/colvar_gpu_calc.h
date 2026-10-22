#ifndef COLVAR_GPU_CALC_H
#define COLVAR_GPU_CALC_H

#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)

#include <vector>
#include "colvar_gpu_support.h"

class colvardeps;
class colvar;
class colvarmodule;

/**
 * @file colvar_gpu_calc.h
 * @brief Declaration of the class for GPU calculation of CVCs
 */

namespace colvars_gpu {
/**
 * @brief Class for managing the GPU calculation of CVCs
 *
 * @note This class is basically a GPU version of colvarmodule::calc_colvars().
 * In the future, the atom groups are supposed to be shared among CVCs,
 * so instead of iterating over all colvars and computing them one by one,
 * this class will flatten the colvardeps* tree, then compute all
 * atom groups in one pass with CUDA graphs, then all CVCs, and finally
 * collect the colvar values.
 */
class colvarmodule_gpu_calc {
public:
  /**
   * @brief A struct for holding a compute node in the CUDA graph
   * and its corresponding atom group or CVC (child of colvardeps)
   *
   * @note There should be an independent data structure for AST
   * For the time being I have to hack the colvardeps which serves
   * partially as an AST (without any type information). There are two
   * kinds of dependencies, namely (i) the dependencies between computational
   * operations (cv-cvc, cvc-cvc and cvc-atom groups), and (ii) the
   * dependencies of different features. (i) should be similar to an AST
   * in a compiler, while (ii) should be a checker walking over the AST.
   */
  struct compute_node_t {
    /// \brief Pointer to the corresponding colvardeps object
    colvardeps* colvar_node;
    /// \brief CUDA graph node
    cudaGraphNode_t child_graph_node;
    /// \brief Whether the node requires CPU buffers
    bool require_cpu_buffers;
  };
  /**
   * @brief A struct for holding a CUDA graph and its execution object
   */
  class compute_gpu_graph_t {
  public:
    /// \brief Constructor
    compute_gpu_graph_t();
    /// \brief (Re)initialize the CUDA graph
    int init();
    /// \brief Destructor
    ~compute_gpu_graph_t();
    /// \brief Dump the CUDA graph to a dot file for debugging
    void dump_graph(const char* filename);
    /// \brief Flag to describe whether the graph execution instance has been initialized
    bool graph_exec_initialized;
    /// \brief CUDA graph object
    cudaGraph_t graph;
    /// \brief CUDA graph execution instance object
    cudaGraphExec_t graph_exec;
    /// \brief List of compute nodes
    std::vector<compute_node_t> nodes;
  };
  /**
   * @brief Constructor
   */
  colvarmodule_gpu_calc();
  /**
   * @brief Destructor
   */
  ~colvarmodule_gpu_calc() {}
  /**
   * @brief Initialize all the GPU compute graphs
   *
   * @return COLVARS_OK if succeeded
   */
  int init();
  /**
   * @brief Calculate the `colvar {...}` blocks
   *
   * This function performs the following operations:
   * 1. Call `cv_update_flags` to update the flags of all colvars and cvcs
   * 2. If the current step is greater than 0, call `cvc_calc_total_force`
   *    to calculate the total force on each CVC
   * 3. Call `atom_group_read_data_gpu` to read the atom positions to GPU
   *    buffers and calculate the required properties
   * 4. Call `cvc_calc_value` to calculate the values of all CVCs
   * 5. Call `cvc_calc_gradients` to calculate the gradients of all CVCs
   * 6. Call `atom_group_calc_fit_gradients` to calculate the fit gradients
   * 7. Call `cvc_calc_Jacobian_derivative` to calculate the Jacobian
   *    derivatives of all CVCs
   * 8. Call `cv_collect_cvc_data` to sum up the colvar values
   * On the first call, the CUDA graphs will be constructed and
   * instantiated. On subsequent calls, the existing CUDA graphs will be
   * launched for "atom_group_read_data_gpu" and "atom_group_calc_fit_gradients".
   *
   * @param[in] colvars A vector of all colvar objects
   * @param[in] cvmodule The main colvarmodule object
   * @return COLVARS_OK if succeeded
   */
  int calc_cvs(const std::vector<colvar*>& colvars, colvarmodule* cvmodule);
  /**
   * @brief Apply the forces to the atom groups from the CVCs
   *
   * On the first call, the CUDA graph for apply atom group forces will be constructed and
   * instantiated. On subsequent calls, the existing CUDA graph will be launched.
   *
   * @param[in] colvars A vector of all colvar objects
   * @param[in] cvmodule The main colvarmodule object
   * @return COLVARS_OK if succeeded
   */
  int apply_forces(const std::vector<colvar*>& colvars, colvarmodule* cvmodule);
private:
  /// \brief CUDA graph for reading data to GPU and calculating required properties
  compute_gpu_graph_t read_data_compute;
  /// \brief CUDA graph for calculating CVCs
  compute_gpu_graph_t calc_value_compute;
  /// \brief CUDA graph for calculating CVC gradients
  compute_gpu_graph_t calc_gradients_compute;
  /// \brief CUDA graph for calculating CVC Jacobians
  compute_gpu_graph_t calc_Jacobian_derivative_compute;
  /// \brief CUDA graph for calculating CVC total forces
  compute_gpu_graph_t calc_total_force_compute;
  /// \brief CUDA graph for calculating fit gradients
  compute_gpu_graph_t calc_fit_gradients_compute;
  /// \brief CUDA graph for applying forces to atom groups
  compute_gpu_graph_t apply_forces_compute;
  /**
   * @brief A list of atom groups that require forces to be applied
   *
   * @note I have no way to forward declare cvm::atom_group here so I
   * have to use colvardeps as a workaround.
   */
  std::vector<colvardeps*> forced_atom_groups;

  #if defined (COLVARS_NVTX_PROFILING)
  /// @name Timers for profiling
  /// @{
  colvars_gpu::colvar_nvtx_prof ag_read_data_prof;
  colvars_gpu::colvar_nvtx_prof cvc_calc_value_prof;
  colvars_gpu::colvar_nvtx_prof cvc_calc_gradients_prof;
  colvars_gpu::colvar_nvtx_prof ag_calc_fit_gradients_prof;
  colvars_gpu::colvar_nvtx_prof cvc_calc_Jacobian_derivative_prof;
  colvars_gpu::colvar_nvtx_prof cvc_calc_total_force_prof;
  colvars_gpu::colvar_nvtx_prof cv_collect_cvc_data_prof;
  colvars_gpu::colvar_nvtx_prof apply_forces_prof;
  /// @}
  #endif

  /// @name Internal functions used in calc_cvs()
  /// @{
  int cv_update_flags(const std::vector<colvar*>& colvars);
  int cvc_calc_total_force(
    const std::vector<colvar*>& colvars,
    colvarmodule_gpu_calc::compute_gpu_graph_t& g,
    colvarmodule* cvmodule,
    bool use_current_step = false);
  int atom_group_read_data_gpu(
    const std::vector<colvar*>& colvars,
    colvarmodule_gpu_calc::compute_gpu_graph_t& g,
    colvarmodule* cvmodule);
  int cvc_calc_value(
    const std::vector<colvar*>& colvars,
    colvarmodule_gpu_calc::compute_gpu_graph_t& g,
    colvarmodule* cvmodule);
  int cvc_calc_gradients(
    const std::vector<colvar*>& colvars,
    colvarmodule_gpu_calc::compute_gpu_graph_t& g,
    colvarmodule* cvmodule);
  int atom_group_calc_fit_gradients(
    const std::vector<colvar*>& colvars,
    colvarmodule_gpu_calc::compute_gpu_graph_t& g,
    colvarmodule* cvmodule);
  int cvc_debug_gradients(
    const std::vector<colvar*>& colvars,
    colvarmodule* cvmodule);
  int cvc_calc_Jacobian_derivative(
    const std::vector<colvar*>& colvars,
    colvarmodule_gpu_calc::compute_gpu_graph_t& g,
    colvarmodule* cvmodule);
  int cv_collect_cvc_data(
    const std::vector<colvar*>& colvars,
    colvarmodule* cvmodule);
  /// @}
};
}

#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)
#endif // COLVAR_GPU_CALC_H
