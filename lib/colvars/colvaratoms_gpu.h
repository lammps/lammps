#ifndef COLVARATOMS_GPU_H
#define COLVARATOMS_GPU_H

#include "colvarmodule.h"
#include "colvartypes.h"
#include "colvar_rotation_derivative.h"

/**
 * @file colvaratoms_gpu.h
 * @brief Declaration of the class for calculating atom group properties on GPU
 */

namespace colvars_gpu {

#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
/**
 * @brief  A struct for holding GPU atom group buffers
 */
struct colvaratoms_gpu_buffer_t {
  /// \brief GPU atom proxy indices (size: num_atoms)
  int* d_atoms_index;
  /// \brief GPU atom positions (size: 3 * num_atoms)
  cvm::real* d_atoms_pos;
  /// \brief GPU atom charges (size: num_atoms)
  cvm::real* d_atoms_charge;
  /// \brief GPU atom velocities (size: 3 * num_atoms)
  cvm::real* d_atoms_vel;
  /// \brief GPU atom mass (size: num_atoms)
  cvm::real* d_atoms_mass;
  /// \brief GPU atom gradients (size: 3 * num_atoms)
  cvm::real* d_atoms_grad;
  /// \brief GPU atom total forces (size: 3 * num_atoms)
  cvm::real* d_atoms_total_force;
  /// \brief Atom masses divided by total mass (size: num_atoms)
  cvm::real* d_atoms_weight;
  /// \brief GPU atom applied force
  cvm::real* d_atoms_applied_force;
  /// \brief GPU fit gradients
  cvm::real* d_fit_gradients;
  /// \brief GPU reference coordinates for f_ag_center or f_ag_rotate
  cvm::real* d_ref_pos;
  /// \brief GPU atom positions (size: 3 * num_atoms)
  cvm::real* d_atoms_pos_unrotated;
  /// \brief GPU center-of-mass
  cvm::rvector* d_com;
  /// \brief GPU temporary buffer for COM, used for avoiding memset
  cvm::rvector* d_com_tmp;
  /// \brief GPU center-of-geometry
  cvm::rvector* d_cog;
  /// \brief GPU temporary buffer for COG, used for avoiding memset
  cvm::rvector* d_cog_tmp;
  /// \brief GPU center of geometry before any fitting
  cvm::rvector* d_cog_orig;
  /// \brief GPU atomic counter for block reduction
  unsigned int* d_com_cog_tbcount;
  /// \brief Center-of-mass on the host-pinned memory for CPU compatibility
  cvm::rvector* h_com;
  /// \brief Center-of-geometry on the host-pinned memory for CPU compatibility
  cvm::rvector* h_cog;
  /// \brief Center-of-geometry before any fitting on the host-pinned memory for CPU compatibility
  cvm::rvector* h_cog_orig;
  /// \brief GPU center of geometry of the reference coordinates
  cvm::rvector* d_ref_pos_cog;
};

/**
 * @brief  A struct for temporary variables for calculating the fit gradients
 */
struct colvaratoms_gpu_calc_fit_info_t {
  /// \brief Fit gradients due to centering
  double3* d_atom_grad;
  /// \brief Gradients of the CV with respect to the quaternion
  double* d_sum_dxdq;
  /// \brief Gradients of the CV with respect to the correlation matrix
  cvm::rmatrix* d_dxdC;
  /// \brief GPU atomic counter for block reduction
  unsigned int* d_tbcount;
};

/**
 * @brief  A struct for the graph for debug gradients
 */
struct colvaratoms_gpu_debug_graph_t {
  /// \brief Flag to describe whether the graph has been initialized
  bool initialized;
  /// \brief CUDA Graph of calc_required_properties
  cudaGraph_t graph_calc_required_properties;
  /// \brief CUDA Graph execution instance of calc_required_properties
  cudaGraphExec_t graph_exec_calc_required_properties;
};

/**
 * @brief The main class for calculating the atom group properties on GPU
 */
class colvaratoms_gpu {
public:
  /**
   * @brief Constructor
   *
   * All variables are expected to initialize to 0.
   */
  colvaratoms_gpu(colvarmodule *cvmodule_in);
  /**
   * @brief Destructor
   */
  ~colvaratoms_gpu();
  /**
   * @brief Initialize the object
   *
   * This function re-allocates the GPU buffers for COM, COG, atomic counters,
   * reference COG, and the temporary variables used for fit gradients. However,
   * it does not (re)allocate the atom-wise data fields like positions. The
   * atom-wise data fields are expected to be updated when sync_to_gpu_buffers
   * is called.
   *
   * @return COLVARS_OK if succeeded
   */
  int init_gpu();
  /**
   * @brief Destructor
   *
   * @return COLVARS_OK if succeeded
   */
  int destroy_gpu();
  /**
   * @brief Synchronize atom-wise data fields from the CPU buffers
   *
   * This function (re)allocates the GPU buffers for proxy indices, charges,
   * masses, weights, positions, velocities, total forces, gradients, and
   * forces that will be applied. Then the proxy indices, charges, masses,
   * positions, velocities, gradients and total forces will be updated from
   * the corresponding CPU buffers.
   *
   * @param[in] cpu_atoms CPU atom group class
   * @return COLVARS_OK if succeeded
   */
  int sync_to_gpu_buffers(const cvm::atom_group* cpu_atoms);
  /**
   * @brief Clear GPU atom-wise data fields
   *
   * This function clears the GPU buffers for proxy indices, charges, masses,
   * weights, positions, gradients, total forces, velocities, and applied forces.
   *
   * @param[in] cpu_atoms CPU atom group class
   * @return COLVARS_OK if succeeded
   */
  int clear_gpu_buffers(const cvm::atom_group* cpu_atoms);
  /**
   * @brief Add a node for resetting some atom-wise data fields to the CUDA graph
   *
   * This function adds cudaMemset nodes to clear gradients, total forces and
   * applied forces to the CUDA graph. The positions and velocities are not
   * cleared since they are expected to be overwritten.
   *
   * @param[in] cpu_atoms CPU atom group class
   * @param[in] graph CUDA graph object
   * @param[in,out] nodes_map A map to store the added nodes with operation name
   * @return COLVARS_OK if succeeded
   */
  int add_reset_atoms_data_nodes(
    const cvm::atom_group* cpu_atoms, cudaGraph_t& graph,
    std::unordered_map<std::string, cudaGraphNode_t>& nodes_map);
  /**
   * @brief Add a node for reading positions and velocities to the CUDA graph
   *
   * This function adds cudaMemcpy nodes to read positions and velocities from
   * the CPU buffers to the GPU buffers to the CUDA graph.
   *
   * @param[in] cpu_atoms CPU atom group class
   * @param[in] graph CUDA graph object
   * @param[in,out] nodes_map A map to store the added nodes with operation name
   * @return COLVARS_OK if succeeded
   */
  int add_read_positions_nodes(
    const cvm::atom_group* cpu_atoms, cudaGraph_t& graph,
    std::unordered_map<std::string, cudaGraphNode_t>& nodes_map);
  /**
   * @brief Add nodes for calculating the required properties to the CUDA graph
   *
   * This function adds nodes for calculating COM, COG, rotation matrix,
   * rotated positions, and reference COG to the CUDA graph.
   *
   * @param[in] cpu_atoms CPU atom group class
   * @param[in] graph CUDA graph object
   * @param[in,out] nodes_map A map to store the added nodes with operation name
   * @param[in] extra_initial_dependencies Additional dependencies to be added
   *        to the initial nodes
   * @return COLVARS_OK if succeeded
   */
  int add_calc_required_properties_nodes(
    const cvm::atom_group* cpu_atoms, cudaGraph_t& graph,
    std::unordered_map<std::string, cudaGraphNode_t>& nodes_map,
    const std::vector<cudaGraphNode_t>& extra_initial_dependencies = {});
  /**
   * @brief Add nodes for updating the CPU buffers after the GPU calculation
   *
   * This function adds cudaMemcpy nodes to update COM, COG, original COG,
   * rotation matrix, and fit gradients to the CPU buffers.
   *
   * @param[in] cpu_atoms CPU atom group class
   * @param[in] graph CUDA graph object
   * @param[in,out] nodes_map A map to store the added nodes with operation name
   * @return COLVARS_OK if succeeded
   */
  int add_update_cpu_buffers_nodes(
    cvm::atom_group* cpu_atoms, cudaGraph_t& graph,
    std::unordered_map<std::string, cudaGraphNode_t>& nodes_map);
  /**
   * @brief Update the CPU COM, COG and rotation object after GPU synchronization
   *
   * @param[in] cpu_atoms CPU atom group class
   * @param[in] copy_to_cpu If true, copy the COM, COG and rotation object to CPU buffers
   * @param[in] stream CUDA stream to be synchronized
   * @return COLVARS_OK if succeeded
   */
  int after_read_data_sync(
    cvm::atom_group* cpu_atoms, bool copy_to_cpu, cudaStream_t stream);
  /**
   * @brief Clear the CPU force buffer for scalar components before applying forces on GPU
   */
  int begin_apply_force_gpu();
  /**
   * @brief Add nodes for applying forces to the CUDA graph
   *
   * @param[in] cpu_atoms CPU atom group class
   * @param[in] graph CUDA graph object
   * @param[in,out] nodes_map A map to store the added nodes with operation name
   * @param[in] extra_initial_dependencies Additional dependencies to be added
   *        to the initial nodes
   * @return COLVARS_OK if succeeded
   */
  int add_apply_force_nodes(
    const cvm::atom_group* cpu_atoms, cudaGraph_t& graph,
    std::unordered_map<std::string, cudaGraphNode_t>& nodes_map,
    const std::vector<cudaGraphNode_t>& extra_initial_dependencies = {});
  /**
   * @brief Add nodes for calculating the fit gradients to the CUDA graph
   *
   * @param[in] cpu_atoms CPU atom group class
   * @param[in] graph CUDA graph object
   * @param[in,out] nodes_map A map to store the added nodes with operation name
   * @param[in] use_cpu_buffers If true, add a node for copying fit gradients from GPU to CPU
   * @return COLVARS_OK if succeeded
   */
  int add_calc_fit_gradients_nodes(
    cvm::atom_group* cpu_atoms, cudaGraph_t& graph,
    std::unordered_map<std::string, cudaGraphNode_t>& nodes_map,
    bool use_cpu_buffers = false);
  /**
   * @brief Function for reading atom positions used for debug gradients
   *
   * This function reads the atom positions from the proxy buffer to the atom group GPU buffer,
   * and optionally modifies one of the coordinates for testing the numerical gradients.
   *
   * @param[in] cpu_atoms CPU atom group class
   * @param[in] change_fitting_group If true, change the fitting group atom position
   * @param[in] change_atom_i Index of the atom to be changed
   * @param[in] xyz Coordinate index (0, 1, or 2) to be changed
   * @param[in] to_cpu If true, copy the changed positions to CPU buffers
   * @param[in] sign Sign of the change (+1 or -1)
   * @param[in] stream CUDA stream to be used
   * @return COLVARS_OK if succeeded
   */
  int read_positions_gpu_debug(
    cvm::atom_group* cpu_atoms, bool change_fitting_group, size_t change_atom_i,
    int xyz, bool to_cpu, double sign, cudaStream_t stream);
  /**
   * @brief Function for calculating the required properties used for debug gradients
   *
   * @param[in] cpu_atoms CPU atom group class
   * @param[in] to_cpu If true, copy the calculated properties to CPU buffers
   * @param[in] stream CUDA stream to be used
   * @return COLVARS_OK if succeeded
   */
  int calc_required_properties_gpu_debug(
    cvm::atom_group* cpu_atoms, bool to_cpu, cudaStream_t stream);
  /**
   * @brief Function to be called when a colvardeps feature is enabled
   *
   * This function performs any side effects needed when a feature is enabled.
   * Specifically, it allocates the GPU buffers for unrotated positions if
   * f_ag_fit_gradients is enabled, and sets up the rotation object if
   * f_ag_rotate is enabled.
   *
   * @param[in] cpu_atoms CPU atom group class
   * @param[in] id Feature ID that has been enabled
   */
  void do_feature_side_effects_gpu(
    cvm::atom_group* cpu_atoms, int id);
  /**
   * @brief Setup the rotation object and copy the reference positions to GPU
   *
   * @param[in] cpu_atoms CPU atom group class
   * @return COLVARS_OK if succeeded
   */
  int setup_rotation(const cvm::atom_group* cpu_atoms);
  /**
   * @brief Setup the rotation derivative object on GPU
   *
   * @param[in] cpu_atoms CPU atom group class
   * @return COLVARS_OK if succeeded
   */
  int setup_rotation_derivative(const cvm::atom_group* cpu_atoms);
  /**
   * @brief Get the GPU buffers
   *
   * @return Reference to the GPU buffers struct
   */
  colvaratoms_gpu_buffer_t& get_gpu_buffers() { return gpu_buffers; }
  /**
   * @brief Read the total forces from the proxy buffer to the GPU buffer
   *
   * @param[in] cpu_atoms CPU atom group class
   * @return COLVARS_OK if succeeded
   */
  int read_total_forces(cvm::atom_group* cpu_atoms);
  /**
   * @brief Function to intercept the forces applied from the CPU interface
   *
   * This function is called when apply_colvar_force() is called on CPU.
   * It stores the applied force in a temporary variable, which will be
   * added to the GPU buffer.
   *
   * @param[in] cpu_force Force applied from the CPU interface
   */
  void apply_colvar_force_from_cpu(cvm::real const& cpu_force);
  /**
   * @brief Set whether to use the CPU atom group force
   *
   * This function is called when group_force_object is used on CPU.
   * It sets a flag to indicate that the group force should be
   * added to the GPU buffer.
   *
   * @param[in] yesno True to use the CPU group force, false otherwise
   */
  void set_use_cpu_group_force(bool yesno) { use_group_force = yesno; }
  /**
   * @brief Getter of the internal GPU rotation object
   */
  colvars_gpu::rotation_gpu& get_rot_gpu() { return rot_gpu; }
  const colvars_gpu::rotation_gpu& get_rot_gpu() const { return rot_gpu; }
  /**
   * @brief Getter of the internal GPU rotation derivative object
   */
  colvars_gpu::rotation_derivative_gpu* get_rot_deriv_gpu() { return rot_deriv_gpu; }
  const colvars_gpu::rotation_derivative_gpu* get_rot_deriv_gpu() const { return rot_deriv_gpu; }
private:
  colvaratoms_gpu_buffer_t gpu_buffers;
  /// \brief Temporary variables for calc_fit_gradients GPU kernel
  colvaratoms_gpu_calc_fit_info_t calc_fit_gradients_gpu_info;
  /// \brief Temporary variables for calc_fit_forces (or "calc_fit_gradients" for vector CVCs) GPU kernel
  colvaratoms_gpu_calc_fit_info_t calc_fit_forces_gpu_info;
  /// \brief Separate CUDA graphs for supporting debug gradients
  colvaratoms_gpu_debug_graph_t debug_graphs;
  /// \brief For intercepting the forces applied from the CPU interface
  cvm::real* h_sum_applied_colvar_force;
  /// \brief If the CPU code path use apply_colvar_force(),
  /// this will be set to true, and then reset to false in begin_apply_force_gpu()
  bool use_apply_colvar_force;
  /// \brief If the CPU code path use group_force_object,
  /// this will be set to true, and then reset to false in begin_apply_force_gpu()
  bool use_group_force;
  /// \brief GPU rotation object
  colvars_gpu::rotation_gpu rot_gpu;
  /// \brief GPU Rotation derivative;
  colvars_gpu::rotation_derivative_gpu* rot_deriv_gpu;
  /// \brief Pointer to the parent colvarmodule
  colvarmodule *cvmodule;
};

#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)
}

#endif // COLVARATOMS_GPU_H
