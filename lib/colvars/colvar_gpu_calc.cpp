#include "colvar_gpu_calc.h"
#include "colvar_gpu_support.h"
#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvar.h"
#include "colvarcomp.h"

#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)

using namespace colvars_gpu;

int colvarmodule_gpu_calc::cv_update_flags(const std::vector<colvar*>& colvars) {
  int error_code = COLVARS_OK;
  for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
    error_code |= (*cvi)->update_cvc_flags();
    if (error_code != COLVARS_OK) return error_code;
  }
  return error_code;
}

int colvarmodule_gpu_calc::cvc_calc_total_force(
  const std::vector<colvar*>& colvars,
  colvarmodule_gpu_calc::compute_gpu_graph_t& g,
  colvarmodule* cvmodule,
  bool use_current_step) {
  int error_code = COLVARS_OK;
  const bool total_force_valid = cvmodule->proxy ? cvmodule->proxy->total_forces_valid() : false;
  if (!total_force_valid) {
    return error_code;
  }
  colvarproxy* p = cvmodule->proxy;
  cudaStream_t stream = p->get_default_stream();
#if defined (COLVARS_NVTX_PROFILING)
  cvc_calc_total_force_prof.start();
#endif // defined (COLVARS_NVTX_PROFILING)
  if (!g.graph_exec_initialized) {
    using node_map_t = std::unordered_map<std::string, cudaGraphNode_t>;
    for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
      // Calculate CVC total force
      if (!(*cvi)->is_enabled(colvardeps::f_cv_total_force_calc)) continue;
      const bool do_total_force =
      use_current_step ?
        (*cvi)->is_enabled(colvardeps::f_cv_total_force_current_step) :
        !(*cvi)->is_enabled(colvardeps::f_cv_total_force_current_step);
      if (do_total_force) {
        const auto all_cvcs = (*cvi)->get_cvcs();
        for (auto cvc = all_cvcs.begin(); cvc != all_cvcs.end(); ++cvc) {
          if (!(*cvc)->is_enabled(colvardeps::f_cvc_active)) continue;
          node_map_t cvc_node_map;
          cudaGraph_t cvc_graph;
          if ((*cvc)->has_gpu_implementation()) {
            error_code |= checkGPUError(cudaGraphCreate(&cvc_graph, 0));
            error_code |= (*cvc)->add_calc_force_invgrads_node(cvc_graph, cvc_node_map);
            cudaGraphNode_t child_graph_node;
            error_code |= checkGPUError(cudaGraphAddChildGraphNode(
              &child_graph_node, g.graph, NULL, 0, cvc_graph));
            error_code |= checkGPUError(cudaGraphDestroy(cvc_graph));
            if (error_code != COLVARS_OK) return error_code;
            // Let's assume that the outermost scope (aka colvars {...} blocks)
            // always require CPU buffers
            g.nodes.push_back({cvc->get(), child_graph_node, true});
          }
        }
      }
    }
    if (!g.nodes.empty()) {
      cudaGraphInstantiateParams params{0};
      params.flags = cudaGraphInstantiateFlagUpload;
      params.uploadStream = stream;
      error_code |= checkGPUError(cudaGraphInstantiateWithParams(&g.graph_exec, g.graph, &params));
      if (params.result_out != cudaGraphInstantiateSuccess) {
        error_code |= cvmodule->error("Failed to instantiate CUDA graph!", COLVARS_ERROR);
      }
      g.graph_exec_initialized = true;
    }
  }
  if (cvmodule->debug()) {
    if (g.graph_exec_initialized) {
      // Debug graph
      const std::string filename = cvmodule->output_prefix() + "_calc_total_force.dot";
      cvmodule->log("Writing calc_total_force graph to " + filename);
      g.dump_graph(filename.c_str());
    }
  }
  if (g.graph_exec_initialized) {
    // Launch the graph
    error_code |= checkGPUError(cudaGraphLaunch(g.graph_exec, stream));
    if (error_code != COLVARS_OK) return error_code;
  }
  for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
    // Calculate CVC total force
    if (!(*cvi)->is_enabled(colvardeps::f_cv_total_force_calc)) continue;
    const bool do_total_force =
      use_current_step ?
       (*cvi)->is_enabled(colvardeps::f_cv_total_force_current_step) :
      !(*cvi)->is_enabled(colvardeps::f_cv_total_force_current_step);
    if (do_total_force) {
      if (cvmodule->debug()) {
        cvmodule->log("Calculating total force of colvar \""+(*cvi)->name+"\".\n");
      }
      const auto all_cvcs = (*cvi)->get_cvcs();
      for (auto cvc = all_cvcs.begin(); cvc != all_cvcs.end(); ++cvc) {
        if (!(*cvc)->is_enabled(colvardeps::f_cvc_active)) continue;
        if (!(*cvc)->has_gpu_implementation()) {
          (*cvc)->calc_force_invgrads();
        }
      }
      if (cvmodule->debug()) {
        cvmodule->log("Done calculating total force of colvar \""+(*cvi)->name+"\".\n");
      }
    }
  }
  if (g.graph_exec_initialized) {
    error_code |= checkGPUError(cudaStreamSynchronize(stream));
    if (error_code != COLVARS_OK) return error_code;
  }
  for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
    // Calculate CVC total force
    if (!(*cvi)->is_enabled(colvardeps::f_cv_total_force_calc)) continue;
    const bool do_total_force =
      use_current_step ?
       (*cvi)->is_enabled(colvardeps::f_cv_total_force_current_step) :
      !(*cvi)->is_enabled(colvardeps::f_cv_total_force_current_step);
    const auto all_cvcs = (*cvi)->get_cvcs();
    if (do_total_force) {
      for (auto cvc = all_cvcs.begin(); cvc != all_cvcs.end(); ++cvc) {
        if (!(*cvc)->is_enabled(colvardeps::f_cvc_active)) continue;
        if ((*cvc)->has_gpu_implementation()) {
          error_code |= (*cvc)->calc_force_invgrads_after_gpu();
        }
      }
    }
  }
#if defined (COLVARS_NVTX_PROFILING)
  cvc_calc_total_force_prof.stop();
#endif // defined (COLVARS_NVTX_PROFILING)
  return error_code;
}

int colvarmodule_gpu_calc::atom_group_read_data_gpu(
  const std::vector<colvar*>& colvars,
  colvarmodule_gpu_calc::compute_gpu_graph_t& g,
  colvarmodule* cvmodule) {
  colvarproxy* p = cvmodule->proxy;
  cudaStream_t stream = p->get_default_stream();
  int error_code = COLVARS_OK;
  if (!g.graph_exec_initialized) {
    // Construct the graph if not initialized
    std::vector<cvm::atom_group *> all_unique_atom_groups;
    // Iterate over all active colvar {...}
    for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
      const auto all_cvcs = (*cvi)->get_cvcs();
      // Iterate over all colvarcomp objects
      for (auto cvc = all_cvcs.begin(); cvc != all_cvcs.end(); ++cvc) {
        if (!(*cvc)->is_enabled(colvardeps::f_cvc_active)) continue;
        if ((*cvc)->is_enabled(colvardeps::f_cvc_explicit_atom_groups)) {
          // Iterate over all atom groups
          // TODO: For the time being, a parent CVC own the atom groups from
          // its children CVCs, but this may be changed by
          // https://github.com/Colvars/colvars/pull/700
          std::vector<cvm::atom_group *>& ags = (*cvc)->atom_groups;
          for (auto ag = ags.begin(); ag != ags.end(); ++ag) {
            if ((*ag)->size() > 0 /*&&
                (*ag)->is_enabled(colvardeps::f_ag_active)*/) {
              all_unique_atom_groups.push_back(*ag);
            }
          }
        }
      }
    }
    // Find all unique atom groups
    std::sort(all_unique_atom_groups.begin(),
              all_unique_atom_groups.end());
    auto last = std::unique(all_unique_atom_groups.begin(),
                            all_unique_atom_groups.end());
    all_unique_atom_groups.erase(last, all_unique_atom_groups.end());
    // Build the CUDA graph for reading data
    using node_map_t = std::unordered_map<std::string, cudaGraphNode_t>;
    for (auto ag = all_unique_atom_groups.begin(); ag != all_unique_atom_groups.end(); ++ag) {
      node_map_t ag_node_map;
      cudaGraph_t ag_graph;
      error_code |= checkGPUError(cudaGraphCreate(&ag_graph, 0));
      if (error_code != COLVARS_OK) return error_code;
      error_code |= (*ag)->get_gpu_atom_group()->add_reset_atoms_data_nodes(*ag, ag_graph, ag_node_map);
      if (error_code != COLVARS_OK) return error_code;
      error_code |= (*ag)->get_gpu_atom_group()->add_read_positions_nodes(*ag, ag_graph, ag_node_map);
      if (error_code != COLVARS_OK) return error_code;
      error_code |= (*ag)->get_gpu_atom_group()->add_calc_required_properties_nodes(*ag, ag_graph, ag_node_map);
      if (error_code != COLVARS_OK) return error_code;
      // Check all parents of this ag to see if any of the parents require CPU buffers
      std::vector<colvardeps*> ag_parents = (*ag)->get_parents();
      bool require_cpu_buffers = false;
      for (auto it = ag_parents.begin(); it != ag_parents.end(); ++it) {
        // TODO: This is slow since the tree of colvardeps only works
        // partially as an AST, and has no way to know the actual type of it.
        try {
          colvar::cvc* cvc = dynamic_cast<colvar::cvc*>(*it);
          if (cvc->is_enabled(colvardeps::f_cvc_require_cpu_buffers)) {
            require_cpu_buffers = true;
            break;
          }
        } catch (const std::bad_cast& ex) {
          // Ignore
        }
      }
      if (require_cpu_buffers) {
        error_code |= (*ag)->get_gpu_atom_group()->add_update_cpu_buffers_nodes(*ag, ag_graph, ag_node_map);
        if (error_code != COLVARS_OK) return error_code;
      }
      cudaGraphNode_t child_graph_node;
      error_code |= checkGPUError(cudaGraphAddChildGraphNode(
        &child_graph_node, g.graph, NULL, 0, ag_graph));
      if (error_code != COLVARS_OK) return error_code;
      // According to the documentation of cudaGraphAddChildGraphNode
      // and https://github.com/NVIDIA/cudnn-frontend/blob/9793df569ce413f4b1844a9176f7ae24dd981603/samples/cpp/misc/cudagraphs.cpp#L188-L195,
      // the child graph has been cloned into the parent graph,
      // so I destroy the child graph here
      error_code |= checkGPUError(cudaGraphDestroy(ag_graph));
      if (error_code != COLVARS_OK) return error_code;
      g.nodes.push_back({*ag, child_graph_node, require_cpu_buffers});
    }
    cudaGraphInstantiateParams params{0};
    params.flags = cudaGraphInstantiateFlagUpload;
    params.uploadStream = stream;
    error_code |= checkGPUError(cudaGraphInstantiateWithParams(&g.graph_exec, g.graph, &params));
    if (params.result_out != cudaGraphInstantiateSuccess) {
      error_code |= cvmodule->error("Failed to instantiate CUDA graph!", COLVARS_ERROR);
    }
    if (error_code != COLVARS_OK) return error_code;
    g.graph_exec_initialized = true;
    if (cvmodule->debug()) {
      // Debug graph
      const std::string filename = cvmodule->output_prefix() + "_read_data.dot";
      cvmodule->log("Writing read data graph to " + filename);
      g.dump_graph(filename.c_str());
    }
  }
  // We need to reset the atom group data if it requires CPU buffers
  for (auto it = g.nodes.begin(); it != g.nodes.end(); ++it) {
    cvm::atom_group* atoms = dynamic_cast<cvm::atom_group*>(it->colvar_node);
    if (it->require_cpu_buffers) {
      atoms->reset_atoms_data();
    }
  }
#if defined (COLVARS_NVTX_PROFILING)
  ag_read_data_prof.start();
#endif // defined (COLVARS_NVTX_PROFILING)
  // Launch the graph
  error_code |= checkGPUError(cudaGraphLaunch(g.graph_exec, stream));
  if (error_code != COLVARS_OK) return error_code;
  error_code |= checkGPUError(cudaStreamSynchronize(stream));
  if (error_code != COLVARS_OK) return error_code;
#if defined (COLVARS_NVTX_PROFILING)
  ag_read_data_prof.stop();
#endif // defined (COLVARS_NVTX_PROFILING)
  // Synchronize the GPU data to CPU if required
  for (auto it = g.nodes.begin(); it != g.nodes.end(); ++it) {
    cvm::atom_group* atoms = dynamic_cast<cvm::atom_group*>(it->colvar_node);
    if (it->require_cpu_buffers) {
      error_code |= atoms->get_gpu_atom_group()->after_read_data_sync(
        atoms, true, stream);
      if (error_code != COLVARS_OK) return error_code;
    }
  }
  return error_code;
}

int colvarmodule_gpu_calc::cvc_calc_value(
  const std::vector<colvar*>& colvars,
  colvarmodule_gpu_calc::compute_gpu_graph_t& g,
  colvarmodule* cvmodule) {
  int error_code = COLVARS_OK;
  colvarproxy* p = cvmodule->proxy;
  cudaStream_t stream = p->get_default_stream();
#if defined (COLVARS_NVTX_PROFILING)
  cvc_calc_value_prof.start();
#endif // defined (COLVARS_NVTX_PROFILING)
  if (!g.graph_exec_initialized) {
    using node_map_t = std::unordered_map<std::string, cudaGraphNode_t>;
    for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
      const auto all_cvcs = (*cvi)->get_cvcs();
      for (auto cvc = all_cvcs.begin(); cvc != all_cvcs.end(); ++cvc) {
        node_map_t cvc_node_map;
        cudaGraph_t cvc_graph;
        if (!(*cvc)->is_enabled(colvardeps::f_cvc_active)) continue;
        if ((*cvc)->has_gpu_implementation()) {
          error_code |= checkGPUError(cudaGraphCreate(&cvc_graph, 0));
          error_code |= (*cvc)->add_calc_value_node(cvc_graph, cvc_node_map);
          cudaGraphNode_t child_graph_node;
          error_code |= checkGPUError(cudaGraphAddChildGraphNode(
            &child_graph_node, g.graph, NULL, 0, cvc_graph));
          error_code |= checkGPUError(cudaGraphDestroy(cvc_graph));
          if (error_code != COLVARS_OK) return error_code;
          // Let's assume that the outermost scope (aka colvars {...} blocks)
          // always require CPU buffers
          g.nodes.push_back({cvc->get(), child_graph_node, true});
        }
      }
    }
    if (!g.nodes.empty()) {
      cudaGraphInstantiateParams params{0};
      params.flags = cudaGraphInstantiateFlagUpload;
      params.uploadStream = stream;
      error_code |= checkGPUError(cudaGraphInstantiateWithParams(&g.graph_exec, g.graph, &params));
      if (params.result_out != cudaGraphInstantiateSuccess) {
        error_code |= cvmodule->error("Failed to instantiate CUDA graph!", COLVARS_ERROR);
      }
      g.graph_exec_initialized = true;
    }
  }
  if (cvmodule->debug()) {
    if (g.graph_exec_initialized) {
      // Debug graph
      const std::string filename = cvmodule->output_prefix() + "_calc_value.dot";
      cvmodule->log("Writing calc_value graph to " + filename);
      g.dump_graph(filename.c_str());
    }
  }
  if (g.graph_exec_initialized) {
    // Launch the graph
    error_code |= checkGPUError(cudaGraphLaunch(g.graph_exec, stream));
    if (error_code != COLVARS_OK) return error_code;
  }
  // Calculate the remaining CPU-only CVCs
  for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
    const auto all_cvcs = (*cvi)->get_cvcs();
    // Iterate over all colvarcomp objects that use only CPU
    for (auto cvc = all_cvcs.begin(); cvc != all_cvcs.end(); ++cvc) {
      if (!(*cvc)->is_enabled(colvardeps::f_cvc_active)) continue;
      if (!(*cvc)->has_gpu_implementation()) {
        (*cvc)->calc_value();
      }
      if (cvmodule->debug()) {
        cvmodule->log(
          "Colvar component "+(*cvc)->name+
          " within colvar \""+(*cvi)->name+"\" has value "+
          cvmodule->to_str((*cvc)->value(),
          cvmodule->cv_width, cvmodule->cv_prec)+".\n");
      }
    }
  }
  if (g.graph_exec_initialized) {
    error_code |= checkGPUError(cudaStreamSynchronize(stream));
    if (error_code != COLVARS_OK) return error_code;
  }
  for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
    const auto all_cvcs = (*cvi)->get_cvcs();
    for (auto cvc = all_cvcs.begin(); cvc != all_cvcs.end(); ++cvc) {
      if (!(*cvc)->is_enabled(colvardeps::f_cvc_active)) continue;
      if ((*cvc)->has_gpu_implementation()) {
        error_code |= (*cvc)->calc_value_after_gpu();
      }
    }
  }
#if defined (COLVARS_NVTX_PROFILING)
  cvc_calc_value_prof.stop();
#endif // defined (COLVARS_NVTX_PROFILING)
  return error_code;
}

int colvarmodule_gpu_calc::cvc_calc_gradients(
  const std::vector<colvar*>& colvars,
  colvarmodule_gpu_calc::compute_gpu_graph_t& g,
  colvarmodule* cvmodule) {
  int error_code = COLVARS_OK;
  colvarproxy* p = cvmodule->proxy;
  cudaStream_t stream = p->get_default_stream();
#if defined (COLVARS_NVTX_PROFILING)
  cvc_calc_gradients_prof.start();
#endif // defined (COLVARS_NVTX_PROFILING)
  if (!g.graph_exec_initialized) {
    using node_map_t = std::unordered_map<std::string, cudaGraphNode_t>;
    for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
      const auto all_cvcs = (*cvi)->get_cvcs();
      for (auto cvc = all_cvcs.begin(); cvc != all_cvcs.end(); ++cvc) {
        node_map_t cvc_node_map;
        cudaGraph_t cvc_graph;
        if (!(*cvc)->is_enabled(colvardeps::f_cvc_active)) continue;
        if ((*cvc)->has_gpu_implementation()) {
          error_code |= checkGPUError(cudaGraphCreate(&cvc_graph, 0));
          error_code |= (*cvc)->add_calc_gradients_node(cvc_graph, cvc_node_map);
          cudaGraphNode_t child_graph_node;
          error_code |= checkGPUError(cudaGraphAddChildGraphNode(
            &child_graph_node, g.graph, NULL, 0, cvc_graph));
          error_code |= checkGPUError(cudaGraphDestroy(cvc_graph));
          if (error_code != COLVARS_OK) return error_code;
          // Let's assume that the outermost scope (aka colvars {...} blocks)
          // always require CPU buffers
          g.nodes.push_back({cvc->get(), child_graph_node, true});
        }
      }
    }
    if (!g.nodes.empty()) {
      cudaGraphInstantiateParams params{0};
      params.flags = cudaGraphInstantiateFlagUpload;
      params.uploadStream = stream;
      error_code |= checkGPUError(cudaGraphInstantiateWithParams(&g.graph_exec, g.graph, &params));
      if (params.result_out != cudaGraphInstantiateSuccess) {
        error_code |= cvmodule->error("Failed to instantiate CUDA graph!", COLVARS_ERROR);
      }
      g.graph_exec_initialized = true;
    }
  }
  if (cvmodule->debug()) {
    if (g.graph_exec_initialized) {
      // Debug graph
      const std::string filename = cvmodule->output_prefix() + "_calc_gradients.dot";
      cvmodule->log("Writing calc_gradients graph to " + filename);
      g.dump_graph(filename.c_str());
    }
  }
  if (g.graph_exec_initialized) {
    // Launch the graph
    error_code |= checkGPUError(cudaGraphLaunch(g.graph_exec, stream));
    if (error_code != COLVARS_OK) return error_code;
  }
  for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
    const auto all_cvcs = (*cvi)->get_cvcs();
    // Iterate over all colvarcomp objects
    for (auto cvc = all_cvcs.begin(); cvc != all_cvcs.end(); ++cvc) {
      if (!(*cvc)->is_enabled(colvardeps::f_cvc_active)) continue;
      if ((*cvc)->is_enabled(colvardeps::f_cvc_gradient)) {
        if (!(*cvc)->has_gpu_implementation()) {
          (*cvc)->calc_gradients();
        }
      }
    }
  }
  if (g.graph_exec_initialized) {
    error_code |= checkGPUError(cudaStreamSynchronize(stream));
    if (error_code != COLVARS_OK) return error_code;
  }
#if defined (COLVARS_NVTX_PROFILING)
  cvc_calc_gradients_prof.stop();
#endif // defined (COLVARS_NVTX_PROFILING)
  return error_code;
}

int colvarmodule_gpu_calc::atom_group_calc_fit_gradients(
  const std::vector<colvar*>& colvars,
  colvarmodule_gpu_calc::compute_gpu_graph_t& g,
  colvarmodule* cvmodule) {
  int error_code = COLVARS_OK;
  colvarproxy* p = cvmodule->proxy;
  cudaStream_t stream = p->get_default_stream();
  if (!g.graph_exec_initialized) {
    // Find all atom groups requiring fit gradients
    std::vector<cvm::atom_group *> all_unique_atom_groups;
    // Iterate over all active colvar {...}
    for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
      const auto all_cvcs = (*cvi)->get_cvcs();
      // Iterate over all colvarcomp objects
      for (auto cvc = all_cvcs.begin(); cvc != all_cvcs.end(); ++cvc) {
        if (!(*cvc)->is_enabled(colvardeps::f_cvc_active)) continue;
        if (!(*cvc)->is_enabled(colvardeps::f_cvc_gradient)) continue;
        if ((*cvc)->is_enabled(colvardeps::f_cvc_explicit_atom_groups)) {
          // Iterate over all atom groups
          // TODO: For the time being, a parent CVC own the atom groups from
          // its children CVCs, but this may be changed by
          // https://github.com/Colvars/colvars/pull/700
          std::vector<cvm::atom_group *>& ags = (*cvc)->atom_groups;
          for (auto ag = ags.begin(); ag != ags.end(); ++ag) {
            if ((*ag)->size() > 0 &&
                // (*ag)->is_enabled(colvardeps::f_ag_active) &&
                (*ag)->is_enabled(colvardeps::f_ag_fit_gradients)) {
              all_unique_atom_groups.push_back(*ag);
            }
          }
        }
      }
    }
    std::sort(all_unique_atom_groups.begin(),
              all_unique_atom_groups.end());
    auto last = std::unique(all_unique_atom_groups.begin(),
                            all_unique_atom_groups.end());
    all_unique_atom_groups.erase(last, all_unique_atom_groups.end());
    using node_map_t = std::unordered_map<std::string, cudaGraphNode_t>;
    for (auto ag = all_unique_atom_groups.begin();
         ag != all_unique_atom_groups.end(); ++ag) {
      node_map_t ag_node_map;
      cudaGraph_t ag_graph;
      error_code |= checkGPUError(cudaGraphCreate(&ag_graph, 0));
      if (error_code != COLVARS_OK) return error_code;
      // Check all parents of this ag to see if any of the parents require CPU buffers
      std::vector<colvardeps*> ag_parents = (*ag)->get_parents();
      bool require_cpu_buffers = false;
      for (auto it = ag_parents.begin(); it != ag_parents.end(); ++it) {
        // TODO: This is slow since the tree of colvardeps only works
        // partially as an AST, and has no way to know the actual type of it.
        try {
          colvar::cvc* cvc = dynamic_cast<colvar::cvc*>(*it);
          if (cvc->is_enabled(colvardeps::f_cvc_require_cpu_buffers)) {
            require_cpu_buffers = true;
            break;
          }
        } catch (const std::bad_cast& ex) {
          // Ignore
        }
      }
      error_code |= (*ag)->get_gpu_atom_group()->add_calc_fit_gradients_nodes(
        *ag, ag_graph, ag_node_map, require_cpu_buffers);
      if (error_code != COLVARS_OK) return error_code;
      cudaGraphNode_t child_graph_node;
      error_code |= checkGPUError(cudaGraphAddChildGraphNode(
        &child_graph_node, g.graph, NULL, 0, ag_graph));
      if (error_code != COLVARS_OK) return error_code;
      error_code |= checkGPUError(cudaGraphDestroy(ag_graph));
      if (error_code != COLVARS_OK) return error_code;
      g.nodes.push_back({*ag, child_graph_node, require_cpu_buffers});
    }
    cudaGraphInstantiateParams params{0};
    params.flags = cudaGraphInstantiateFlagUpload;
    params.uploadStream = stream;
    error_code |= checkGPUError(cudaGraphInstantiateWithParams(&g.graph_exec, g.graph, &params));
    if (params.result_out != cudaGraphInstantiateSuccess) {
      error_code |= cvmodule->error("Failed to instantiate CUDA graph!", COLVARS_ERROR);
    }
    if (error_code != COLVARS_OK) return error_code;
    g.graph_exec_initialized = true;
    if (cvmodule->debug()) {
      // Debug graph
      const std::string filename = cvmodule->output_prefix() + "_fit_gradients.dot";
      cvmodule->log("Writing fit gradients graph to " + filename);
      g.dump_graph(filename.c_str());
    }
  }
#if defined (COLVARS_NVTX_PROFILING)
  ag_calc_fit_gradients_prof.start();
#endif // defined (COLVARS_NVTX_PROFILING)
  // Launch the graph only when there are actually nodes
  if (!g.nodes.empty()) {
    error_code |= checkGPUError(cudaGraphLaunch(g.graph_exec, stream));
    if (error_code != COLVARS_OK) return error_code;
    error_code |= checkGPUError(cudaStreamSynchronize(stream));
    if (error_code != COLVARS_OK) return error_code;
  }
#if defined (COLVARS_NVTX_PROFILING)
  ag_calc_fit_gradients_prof.stop();
#endif // defined (COLVARS_NVTX_PROFILING)
  return error_code;
}

int colvarmodule_gpu_calc::cvc_debug_gradients(const std::vector<colvar*>& colvars, colvarmodule* cvmodule) {
  int error_code = COLVARS_OK;
  for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
    const auto all_cvcs = (*cvi)->get_cvcs();
    for (auto cvc = all_cvcs.begin(); cvc != all_cvcs.end(); ++cvc) {
      if (!(*cvc)->is_enabled(colvardeps::f_cvc_active)) continue;
      if (!(*cvc)->is_enabled(colvardeps::f_cvc_gradient)) continue;
      if (!(*cvc)->is_enabled(colvardeps::f_cvc_debug_gradient)) continue;
      if ((*cvc)->has_gpu_implementation()) {
        (*cvc)->debug_gradients_gpu(calc_value_compute, calc_gradients_compute);
      } else {
        (*cvc)->debug_gradients();
      }
    }
  }
  return error_code;
}

int colvarmodule_gpu_calc::cvc_calc_Jacobian_derivative(
  const std::vector<colvar*>& colvars,
  colvarmodule_gpu_calc::compute_gpu_graph_t& g,
  colvarmodule* cvmodule) {
  int error_code = COLVARS_OK;
  colvarproxy* p = cvmodule->proxy;
  cudaStream_t stream = p->get_default_stream();
#if defined (COLVARS_NVTX_PROFILING)
  cvc_calc_Jacobian_derivative_prof.start();
#endif // defined (COLVARS_NVTX_PROFILING)
  if (!g.graph_exec_initialized) {
    using node_map_t = std::unordered_map<std::string, cudaGraphNode_t>;
    for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
      const auto all_cvcs = (*cvi)->get_cvcs();
      const bool calc_jacobian = (*cvi)->is_enabled(colvardeps::f_cv_Jacobian);
      for (auto cvc = all_cvcs.begin(); cvc != all_cvcs.end(); ++cvc) {
        node_map_t cvc_node_map;
        cudaGraph_t cvc_graph;
        if (!(*cvc)->is_enabled(colvardeps::f_cvc_active)) continue;
        if (calc_jacobian) {
          if ((*cvc)->has_gpu_implementation()) {
            error_code |= checkGPUError(cudaGraphCreate(&cvc_graph, 0));
            error_code |= (*cvc)->add_calc_Jacobian_derivative_node(cvc_graph, cvc_node_map);
            cudaGraphNode_t child_graph_node;
            error_code |= checkGPUError(cudaGraphAddChildGraphNode(
              &child_graph_node, g.graph, NULL, 0, cvc_graph));
            error_code |= checkGPUError(cudaGraphDestroy(cvc_graph));
            if (error_code != COLVARS_OK) return error_code;
            // Let's assume that the outermost scope (aka colvars {...} blocks)
            // always require CPU buffers
            g.nodes.push_back({cvc->get(), child_graph_node, true});
          }
        }
      }
    }
    if (!g.nodes.empty()) {
      cudaGraphInstantiateParams params{0};
      params.flags = cudaGraphInstantiateFlagUpload;
      params.uploadStream = stream;
      error_code |= checkGPUError(cudaGraphInstantiateWithParams(&g.graph_exec, g.graph, &params));
      if (params.result_out != cudaGraphInstantiateSuccess) {
        error_code |= cvmodule->error("Failed to instantiate CUDA graph!", COLVARS_ERROR);
      }
      g.graph_exec_initialized = true;
    }
  }
  if (cvmodule->debug()) {
    if (g.graph_exec_initialized) {
      // Debug graph
      const std::string filename = cvmodule->output_prefix() + "_calc_Jacobian_derivative.dot";
      cvmodule->log("Writing calc_Jacobian_derivative graph to " + filename);
      g.dump_graph(filename.c_str());
    }
  }
  if (g.graph_exec_initialized) {
    // Launch the graph
    error_code |= checkGPUError(cudaGraphLaunch(g.graph_exec, stream));
    if (error_code != COLVARS_OK) return error_code;
  }
  for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
    const auto all_cvcs = (*cvi)->get_cvcs();
    const bool calc_jacobian = (*cvi)->is_enabled(colvardeps::f_cv_Jacobian);
    for (auto cvc = all_cvcs.begin(); cvc != all_cvcs.end(); ++cvc) {
      if (!(*cvc)->is_enabled(colvardeps::f_cvc_active)) continue;
      if (calc_jacobian) {
        if (!(*cvc)->has_gpu_implementation()) {
          (*cvc)->calc_Jacobian_derivative();
        }
      }
    }
  }
  if (g.graph_exec_initialized) {
    error_code |= checkGPUError(cudaStreamSynchronize(stream));
    if (error_code != COLVARS_OK) return error_code;
  }
  for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
    const auto all_cvcs = (*cvi)->get_cvcs();
    const bool calc_jacobian = (*cvi)->is_enabled(colvardeps::f_cv_Jacobian);
    for (auto cvc = all_cvcs.begin(); cvc != all_cvcs.end(); ++cvc) {
      if (!(*cvc)->is_enabled(colvardeps::f_cvc_active)) continue;
      if (calc_jacobian) {
        if ((*cvc)->has_gpu_implementation()) {
          error_code |= (*cvc)->calc_Jacobian_derivative_after_gpu();
        }
      }
    }
  }
#if defined (COLVARS_NVTX_PROFILING)
  cvc_calc_Jacobian_derivative_prof.stop();
#endif // defined (COLVARS_NVTX_PROFILING)
  return error_code;
}

int colvarmodule_gpu_calc::cv_collect_cvc_data(const std::vector<colvar*>& colvars, colvarmodule* cvmodule) {
  int error_code = COLVARS_OK;
#if defined (COLVARS_NVTX_PROFILING)
  cv_collect_cvc_data_prof.start();
#endif // defined (COLVARS_NVTX_PROFILING)
  for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
    // If the total forces are not available, it will be reset in
    // collect_cvc_total_forces anyway (called by collect_cvc_data)
    error_code |= (*cvi)->collect_cvc_data();
    if (cvmodule->get_error()) {
      return COLVARS_ERROR;
    }
  }
#if defined (COLVARS_NVTX_PROFILING)
  cv_collect_cvc_data_prof.stop();
#endif // defined (COLVARS_NVTX_PROFILING)
  return error_code;
}

colvarmodule_gpu_calc::compute_gpu_graph_t::compute_gpu_graph_t():
graph_exec_initialized(false), graph(NULL), graph_exec(NULL) {
}

int colvarmodule_gpu_calc::compute_gpu_graph_t::init() {
  int error_code = COLVARS_OK;
  graph_exec_initialized = false;
  if (graph_exec) {
    error_code |= checkGPUError(cudaGraphExecDestroy(graph_exec));
    graph_exec = NULL;
  }
  if (graph) {
    error_code |= checkGPUError(cudaGraphDestroy(graph));
    graph = NULL;
  }
  nodes.clear();
  error_code |= checkGPUError(cudaGraphCreate(&graph, 0));
  return error_code;
}

colvarmodule_gpu_calc::compute_gpu_graph_t::~compute_gpu_graph_t() {
  if (graph_exec) {
    checkGPUError(cudaGraphExecDestroy(graph_exec));
    graph_exec = NULL;
  }
  if (graph) {
    checkGPUError(cudaGraphDestroy(graph));
    graph = NULL;
  }
  nodes.clear();
  graph_exec_initialized = false;
}

void colvarmodule_gpu_calc::compute_gpu_graph_t::dump_graph(const char* filename) {
  cudaGraphDebugDotFlags dotFlags = cudaGraphDebugDotFlagsVerbose;
  checkGPUError(cudaGraphDebugDotPrint(graph, filename, dotFlags));
}

colvarmodule_gpu_calc::colvarmodule_gpu_calc() {
#if defined (COLVARS_NVTX_PROFILING)
  ag_read_data_prof.set_name_color("ag_read_data", 0xc36eff);
  cvc_calc_value_prof.set_name_color("cvc_calc_value", 0x99dfff);
  cvc_calc_gradients_prof.set_name_color("cvc_calc_gradients", 0x0a78ff);
  ag_calc_fit_gradients_prof.set_name_color("ag_calc_fit_gradients", 0xff6c6e);
  cvc_calc_Jacobian_derivative_prof.set_name_color("cvc_calc_Jacobian_derivative", 0x91ff5e);
  cvc_calc_total_force_prof.set_name_color("cvc_calc_total_force", 0xfff242);
  cv_collect_cvc_data_prof.set_name_color("cv_collect_cvc_data", 0xff742f);
  apply_forces_prof.set_name_color("apply_forces", 0xff89a3);
#endif
}

int colvarmodule_gpu_calc::init() {
  int error_code = COLVARS_OK;
  error_code |= read_data_compute.init();
  error_code |= calc_value_compute.init();
  error_code |= calc_gradients_compute.init();
  error_code |= calc_Jacobian_derivative_compute.init();
  error_code |= calc_total_force_compute.init();
  error_code |= calc_fit_gradients_compute.init();
  error_code |= apply_forces_compute.init();
  forced_atom_groups.clear();
  return error_code;
}

int colvarmodule_gpu_calc::calc_cvs(const std::vector<colvar*>& colvars, colvarmodule* cvmodule) {
  int error_code = COLVARS_OK;
#define checkColvarsError(stmt) do { \
  error_code |= stmt; \
  if (error_code != COLVARS_OK) return error_code;\
} while (0);
  colvarproxy* p = cvmodule->proxy;
  // Update flags
  checkColvarsError(cv_update_flags(colvars));
  // Calculate total force
  if (cvmodule->step_relative() > 0) {
    checkColvarsError(cvc_calc_total_force(colvars, calc_total_force_compute, cvmodule, false));
  }
  // Read data to atom groups
  checkColvarsError(atom_group_read_data_gpu(
    colvars, read_data_compute, cvmodule));
  // Wait for extra information (for example, lattice) from the MD engine
  // before the CPU loop
  checkColvarsError(p->wait_for_extra_info_ready());
  // Calculate CVC values
  checkColvarsError(cvc_calc_value(colvars, calc_value_compute, cvmodule));
  // Calculate CVC gradients
  checkColvarsError(cvc_calc_gradients(colvars, calc_gradients_compute,  cvmodule));
  // Calculate fit gradients for atom groups
  checkColvarsError(atom_group_calc_fit_gradients(
    colvars, calc_fit_gradients_compute, cvmodule));
  // Debug gradients
  checkColvarsError(cvc_debug_gradients(colvars, cvmodule));
  // Calculate the Jacobian terms
  checkColvarsError(cvc_calc_Jacobian_derivative(colvars, calc_Jacobian_derivative_compute, cvmodule));
  // Calculate total force
  checkColvarsError(cvc_calc_total_force(colvars, calc_total_force_compute, cvmodule, true));
  // Collect CVC data
  checkColvarsError(cv_collect_cvc_data(colvars, cvmodule));
#undef checkColvarsError
  return error_code;
}

int colvarmodule_gpu_calc::apply_forces(const std::vector<colvar*>& colvars, colvarmodule* cvmodule) {
  int error_code = COLVARS_OK;
#define checkColvarsError(stmt) do { \
  error_code |= stmt; \
  if (error_code != COLVARS_OK) return error_code;\
} while (0);
  colvarproxy* p = cvmodule->proxy;
  cudaStream_t stream = p->get_default_stream();
  if (!apply_forces_compute.graph_exec_initialized) {
    forced_atom_groups.clear();
    // Find all unique atom groups requiring forces
    std::vector<cvm::atom_group*> all_unique_atom_groups;
    for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
      const auto all_cvcs = (*cvi)->get_cvcs();
      // Iterate over all colvarcomp objects
      if (!(*cvi)->is_enabled(colvardeps::f_cv_apply_force)) continue;
      for (auto cvc = all_cvcs.begin(); cvc != all_cvcs.end(); ++cvc) {
        if (!(*cvc)->is_enabled(colvardeps::f_cvc_active)) continue;
        // Iterate over all atom groups
        // TODO: For the time being, a parent CVC own the atom groups from
        // its children CVCs, but this may be changed by
        // https://github.com/Colvars/colvars/pull/700
        std::vector<cvm::atom_group *>& ags = (*cvc)->atom_groups;
        for (auto ag = ags.begin(); ag != ags.end(); ++ag) {
          if ((*ag)->size() > 0) {
            all_unique_atom_groups.push_back(*ag);
          }
        }
      }
    }
    std::sort(all_unique_atom_groups.begin(),
            all_unique_atom_groups.end());
    auto last = std::unique(all_unique_atom_groups.begin(),
                            all_unique_atom_groups.end());
    all_unique_atom_groups.erase(last, all_unique_atom_groups.end());
    // Save the atom groups
    forced_atom_groups.insert(forced_atom_groups.end(),
                              all_unique_atom_groups.begin(),
                              all_unique_atom_groups.end());
    /*
     * NOTE: There are two CPU routines for applying forces on
     * an atom group, namely (i) apply_colvar_force for a scalar force,
     * and (ii) apply_force for a vector force on the COM. However, at
     * this moment, I have no idea which one is called for a specific
     * atom group, so I need to call communicate_forces() once.
     */
    // Clear the h_sum_applied_colvar_force, use_group_force and use_apply_colvar_force
    for (auto it = all_unique_atom_groups.begin(); it != all_unique_atom_groups.end(); ++it) {
      cvm::atom_group* ag = *it;
      checkColvarsError(ag->get_gpu_atom_group()->begin_apply_force_gpu());
    }
    for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
      if ((*cvi)->is_enabled(colvardeps::f_cv_apply_force)) {
        (*cvi)->communicate_forces();
        if (cvmodule->get_error()) {
          return COLVARS_ERROR;
        }
      }
    }
    /*
     * Now I expect either use_apply_colvar_force or use_group_force
     * is set in any of the atom groups, so build the CUDA graph.
     */
    using node_map_t = std::unordered_map<std::string, cudaGraphNode_t>;
    // Add apply force CUDA graph nodes
    for (auto ag = all_unique_atom_groups.begin(); ag != all_unique_atom_groups.end(); ++ag) {
      node_map_t ag_node_map;
      cudaGraph_t ag_graph;
      checkColvarsError(checkGPUError(cudaGraphCreate(&ag_graph, 0)));
      checkColvarsError((*ag)->get_gpu_atom_group()->add_apply_force_nodes(*ag, ag_graph, ag_node_map));
      cudaGraphNode_t child_graph_node;
      checkColvarsError(checkGPUError(cudaGraphAddChildGraphNode(
        &child_graph_node, apply_forces_compute.graph, NULL, 0, ag_graph)));
      checkColvarsError(checkGPUError(cudaGraphDestroy(ag_graph)));
    }
    cudaGraphInstantiateParams params{0};
    params.flags = cudaGraphInstantiateFlagUpload;
    params.uploadStream = stream;
    checkColvarsError(checkGPUError(cudaGraphInstantiateWithParams(&apply_forces_compute.graph_exec, apply_forces_compute.graph, &params)));
    if (params.result_out != cudaGraphInstantiateSuccess) {
      error_code |= cvmodule->error("Failed to instantiate CUDA graph!", COLVARS_ERROR);
    }
    apply_forces_compute.graph_exec_initialized = true;
    // Debug graph
    if (cvmodule->debug()) {
      const std::string filename = cvmodule->output_prefix() + "_apply_forces.dot";
      cvmodule->log("Writing apply forces graph to " + filename);
      apply_forces_compute.dump_graph(filename.c_str());
    }
  } else {
    for (auto it = forced_atom_groups.begin(); it != forced_atom_groups.end(); ++it) {
      cvm::atom_group* ag = dynamic_cast<cvm::atom_group*>((*it));
      checkColvarsError(ag->get_gpu_atom_group()->begin_apply_force_gpu());
    }
    for (auto cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
      if ((*cvi)->is_enabled(colvardeps::f_cv_apply_force)) {
        (*cvi)->communicate_forces();
        if (cvmodule->get_error()) {
          return COLVARS_ERROR;
        }
      }
    }
  }
#if defined (COLVARS_NVTX_PROFILING)
  apply_forces_prof.start();
#endif
  checkColvarsError(checkGPUError(cudaGraphLaunch(apply_forces_compute.graph_exec, stream)));
  // NOTE: The synchronization here should be unnecessary because we use the proxy stream,
  // which will be synchronized by the MD engine anyway.
  // checkColvarsError(checkGPUError(cudaStreamSynchronize(stream)));
#if defined (COLVARS_NVTX_PROFILING)
  apply_forces_prof.stop();
#endif
#undef checkColvarsError
  return error_code;
}

#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)
