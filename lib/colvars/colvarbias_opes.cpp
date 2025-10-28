// This code is mainly adapted from the PLUMED opes module, which uses the
// LGPLv3 license as shown below:
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2020-2021 of Michele Invernizzi.

   This file is part of the OPES plumed module.

   The OPES plumed module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The OPES plumed module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#if defined(NAMD_VERSION)
#include "CkLoopAPI.h"
#include "CkLambda.h"
#define OPES_THREADING
#endif

#include "colvarbias_opes.h"
#include "colvarbias.h"
#include "colvardeps.h"
#include "colvarproxy.h"
#include "colvars_memstream.h"
#include "colvargrid.h"

#include <exception>
#include <iomanip>
#include <ios>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <limits>
#include <sstream>


colvarbias_opes::colvarbias_opes(char const *key):
  colvarbias(key), m_kbt(0), m_barrier(0), m_biasfactor(0),
  m_bias_prefactor(0), m_temperature(0),
  m_pace(0), m_adaptive_sigma_stride(0),
  m_adaptive_counter(0), m_counter(1),
  m_compression_threshold(0), m_compression_threshold2(0),
  m_adaptive_sigma(false), m_fixed_sigma(false),
  m_no_zed(false), m_nlist(false), m_recursive_merge(true),
  m_nlist_param(2, 0), m_epsilon(0), m_sum_weights(0),
  m_sum_weights2(0), m_cutoff(0), m_cutoff2(0),
  m_zed(1), m_old_kdenorm(0), m_kdenorm(0),
  m_val_at_cutoff(0), m_nlist_center(0), m_nlist_index(0),
  m_nlist_steps(0), m_nlist_update(false),
  m_nlist_pace_reset(false), m_nker(0), m_calc_work(false),
  m_work(0), comm(single_replica), m_num_walkers(1),
  m_num_threads(1), m_nlker(0), m_traj_output_frequency(0),
  m_traj_line(traj_line{0}), m_is_first_step(true),
  m_pmf_grid_on(false), m_reweight_grid(nullptr),
  m_pmf_grid(nullptr), m_pmf_hist_freq(0), m_pmf_shared(true),
  m_explore(false), m_inf_biasfactor(false)
{
#ifdef OPES_THREADING
  provide(f_cvb_smp, cvm::proxy->get_smp_mode() == colvarproxy_smp::smp_mode_t::inner_loop);
  if (is_available(f_cvb_smp)){
    enable(f_cvb_smp); // Enabled by default
  }
#endif
}

int colvarbias_opes::init(const std::string& conf) {
  int error_code = colvarbias::init(conf);
  enable(f_cvb_scalar_variables);
  get_keyval_feature(this, conf, "applyBias", f_cvb_apply_force, true);
  m_temperature = cvm::proxy->target_temperature();
  m_kbt = m_temperature * cvm::proxy->boltzmann();
  get_keyval(conf, "newHillFrequency", m_pace);
  if (m_pace % time_step_factor != 0) {
    error_code |= cvm::error("newHillFrequency (currently " + cvm::to_str(m_pace) +
                                 ") must be a multiple of timeStepFactor (" +
                                 cvm::to_str(time_step_factor) + ").\n",
                             COLVARS_INPUT_ERROR);
  }
  get_keyval(conf, "barrier", m_barrier);
  get_keyval(conf, "explore", m_explore, false);
  if (m_barrier < 0) {
    return cvm::error("the barrier should be greater than zero", COLVARS_INPUT_ERROR);
  }
  std::string biasfactor_str;
  get_keyval(conf, "biasfactor", biasfactor_str);
  if ((cvm::proxy->target_temperature() == 0.0) && cvm::proxy->simulation_running()) {
    cvm::log("WARNING: OPES should not be run without a thermostat or at 0 Kelvin!\n");
  }
  m_biasfactor = m_barrier / m_kbt;
  m_inf_biasfactor = biasfactor_str == "inf" || biasfactor_str == "INF";
  if (m_inf_biasfactor) {
    m_biasfactor = std::numeric_limits<cvm::real>::infinity();
    m_bias_prefactor = 1;
    if (m_explore) {
      return cvm::error("biasfactor cannot be infinity in the explore mode.");
    }
  } else {
    if (biasfactor_str.size() > 0) {
      try {
        m_biasfactor = std::stod(biasfactor_str);
      } catch (const std::exception& e) {
        return cvm::error(e.what(), COLVARS_INPUT_ERROR);
      }
    }
    if (m_biasfactor <= 1.0) {
      return cvm::error("biasfactor must be greater than one (use \"inf\" for uniform target)");
    }
    m_bias_prefactor = 1 - 1.0 / m_biasfactor;
  }
  if (m_explore) {
    m_bias_prefactor = m_biasfactor - 1;
  }
  get_keyval(conf, "adaptiveSigma", m_adaptive_sigma, false);
  m_sigma0.resize(num_variables());
  get_keyval(conf, "gaussianSigma", m_sigma0, std::vector<cvm::real>(num_variables()));
  m_av_cv.assign(num_variables(), 0);
  m_av_M2.assign(num_variables(), 0);
  if (m_adaptive_sigma) {
    get_keyval(conf, "adaptiveSigmaStride", m_adaptive_sigma_stride, 0);
    if (m_inf_biasfactor) {
      return cvm::error("cannot use infinite biasfactor with adaptive sigma",
                        COLVARS_INPUT_ERROR);
    }
    if (m_adaptive_sigma_stride == 0) {
      m_adaptive_sigma_stride = m_pace * 10;
    }
    if (m_adaptive_sigma_stride < m_pace) {
      return cvm::error("It is better to choose an adaptiveSigmaStride >= newHillFrequency.\n", COLVARS_INPUT_ERROR);
    }
  } else {
    if (m_sigma0.size() != num_variables()) {
      return cvm::error("number of sigma parameters does not match the number of variables",
                        COLVARS_INPUT_ERROR);
    }
    if (m_explore) {
      for (size_t i = 0; i < num_variables(); ++i) {
        m_sigma0[i] *= std::sqrt(m_biasfactor);
      }
    }
  }
  get_keyval(conf, "gaussianSigmaMin", m_sigma_min);
  if ((m_sigma_min.size() != 0) && (m_sigma_min.size() != num_variables())) {
    return cvm::error("incorrect number of parameters of gaussianSigmaMin");
  }
  if (m_sigma_min.size() > 0 && !m_adaptive_sigma) {
    for (size_t i = 0; i < num_variables(); ++i) {
      if (m_sigma_min[i] > m_sigma0[i]) {
        return cvm::error("gaussianSigmaMin of variable " + cvm::to_str(i) + " should be smaller than sigma");
      }
    }
  }
  get_keyval(conf, "epsilon", m_epsilon, std::exp(-m_barrier/m_bias_prefactor/m_kbt));
  if (m_epsilon <= 0) {
    return cvm::error("you must choose a value of epsilon greater than zero");
  }
  m_sum_weights = std::pow(m_epsilon, m_bias_prefactor);
  m_sum_weights2 = m_sum_weights * m_sum_weights;
  if (m_explore) {
    get_keyval(conf, "kernelCutoff", m_cutoff, std::sqrt(2.0*m_barrier/m_kbt));
  } else {
    get_keyval(conf, "kernelCutoff", m_cutoff, std::sqrt(2.0*m_barrier/m_bias_prefactor/m_kbt));
  }
  if (m_cutoff <= 0) {
    return cvm::error("you must choose a value of kernelCutoff greater than zero");
  }
  m_cutoff2 = m_cutoff * m_cutoff;
  m_val_at_cutoff = std::exp(-0.5 * m_cutoff2);
  get_keyval(conf, "compressionThreshold", m_compression_threshold, 1);
  if (m_compression_threshold != 0) {
    if (m_compression_threshold < 0 || m_compression_threshold > m_cutoff) {
      return cvm::error("compressionThreshold cannot be smaller than 0 or larger than kernelCutoff", COLVARS_INPUT_ERROR);
    }
  }
  m_compression_threshold2 = m_compression_threshold * m_compression_threshold;
  get_keyval(conf, "neighborList", m_nlist, false);
  if (m_nlist) {
    get_keyval(conf, "neighborListNewHillReset", m_nlist_pace_reset, false);
    std::vector<cvm::real> nlist_param;
    get_keyval(conf, "neighborListParameters", nlist_param, std::vector<cvm::real>());
    if (nlist_param.empty()) {
      m_nlist_param[0] = 3.0; //*cutoff2_ -> max distance of neighbors
      m_nlist_param[1] = 0.5; //*nlist_dev2_[i] -> condition for rebuilding
    } else {
      if (nlist_param.size() != 2) {
        return cvm::error("two cutoff parameters are needed for the neighbor list", COLVARS_INPUT_ERROR);
      }
      if (nlist_param[0] <= 1.0) {
        return cvm::error("the first of neighborListParam must be greater than 1.0. The smaller the first, the smaller should be the second as well", COLVARS_INPUT_ERROR);
      }
      const cvm::real min_PARAM_1 = (1.-1./std::sqrt(nlist_param[0]))+0.16;
      if (nlist_param[1] <= 0) {
        return cvm::error("the second of neighborListParam must be greater than 0", COLVARS_INPUT_ERROR);
      }
      if (nlist_param[1] > min_PARAM_1) {
        return cvm::error("the second of neighborListParam must be smaller to avoid systematic errors. Largest suggested value is: 1.16-1/sqrt(param_0) = " + cvm::to_str(min_PARAM_1), COLVARS_INPUT_ERROR);
      }
      m_nlist_param = nlist_param;
    }
    m_nlist_center.resize(num_variables());
    m_nlist_dev2.resize(num_variables(), 0);
    m_nlist_steps = 0;
    m_nlist_update = true;
  }
  get_keyval(conf, "noZed", m_no_zed, false);
  if (m_no_zed) {
    m_sum_weights = 1;
    m_sum_weights2 = 1;
  }
  get_keyval(conf, "fixedGaussianSigma", m_fixed_sigma, false);
  get_keyval(conf, "recursiveMerge", m_recursive_merge, true);
  get_keyval(conf, "calcWork", m_calc_work, false);
  bool b_replicas = false;
  get_keyval(conf, "multipleReplicas", b_replicas, false);

#ifdef OPES_THREADING
  get_keyval_feature(this, conf, "smp", f_cvb_smp, is_enabled(f_cvb_smp));
  if (is_enabled(f_cvb_smp)) {
    m_num_threads = cvm::proxy->smp_num_threads();
  } else {
    m_num_threads = 1;
  }
#else
  if (m_num_threads > 1) {
    return cvm::error("Multi-threading requested in OPES, which is not supported by this build.\n");
  }
#endif
  bool serial = false;
  get_keyval(conf, "serial", serial, false);
  if (serial) m_num_threads = 1;
  comm = b_replicas ? multiple_replicas : single_replica;
  if (comm == multiple_replicas) {
    colvarproxy *proxy = cvm::main()->proxy;
    get_keyval(conf, "replicaID", replica_id, replica_id);
    get_keyval(conf, "sharedFreq", shared_freq, output_freq);
    if (shared_freq % time_step_factor != 0) {
      error_code |= cvm::error("sharedFreq (currently " + cvm::to_str(shared_freq) +
                                   ") must be a multiple of timeStepFactor (" +
                                   cvm::to_str(time_step_factor) + ").\n",
                               COLVARS_INPUT_ERROR);
    }
    if (!replica_id.size()) {
      if (proxy->check_replicas_enabled() == COLVARS_OK) {
        // Obtain replicaID from the communicator
        replica_id = cvm::to_str(proxy->replica_index());
        cvm::log("Setting replicaID from communication layer: replicaID = "+
                 replica_id+".\n");
      } else {
        error_code |=
            cvm::error("Error: using more than one replica, but replicaID could not be obtained.\n",
                       COLVARS_INPUT_ERROR);
      }
    }
    m_num_walkers = proxy->num_replicas();
  }
  get_keyval(conf, "pmf", m_pmf_grid_on, false);
  if (m_pmf_grid_on) {
    std::vector<std::string> pmf_cv_name;
    get_keyval(conf, "pmfColvars", pmf_cv_name);
    for (auto it = pmf_cv_name.begin(); it != pmf_cv_name.end(); ++it) {
      bool found = false;
      for (size_t i = 0; i < num_variables(); ++i) {
        if (variables(i)->name == (*it)) {
          if (variables(i)->enable(f_cv_grid) != COLVARS_OK) {
            return cvm::error("CV " + (*it) + " does not support grid\n");
          }
          m_pmf_cvs.push_back(variables(i));
          found = true;
          break;
        }
      }
      if (!found) {
        return cvm::error("CV " + (*it) + " not found\n");
      }
    }
    key_lookup(conf, "grid", &grid_conf);
    m_reweight_grid.reset(new colvar_grid_scalar(m_pmf_cvs, nullptr, false, grid_conf));
    m_pmf_grid.reset(new colvar_grid_scalar(m_pmf_cvs, m_reweight_grid));
    get_keyval(conf, "pmfHistoryFrequency", m_pmf_hist_freq, output_freq);
    if ((m_pmf_hist_freq % output_freq) != 0) {
      error_code |= cvm::error("Error: pmfHistoryFrequency must be a multiple of outputFreq.\n",
                               COLVARS_INPUT_ERROR);
    }
    if (comm == multiple_replicas) {
      get_keyval(conf, "pmfShared", m_pmf_shared, true);
      if (m_pmf_shared) {
        m_global_reweight_grid.reset(new colvar_grid_scalar(m_pmf_cvs, m_reweight_grid));
        m_global_pmf_grid.reset(new colvar_grid_scalar(m_pmf_cvs, m_reweight_grid));
      }
    }
  }
  m_kdenorm = m_explore? m_counter : m_sum_weights;
  m_old_kdenorm = m_kdenorm;
  m_traj_line.rct = m_kbt * cvm::logn(m_sum_weights / m_counter);
  m_traj_line.zed = m_zed;
  m_traj_line.neff = (1 + m_sum_weights) * (1 + m_sum_weights) / (1 + m_sum_weights2);
  m_traj_line.nker = m_kernels.size();
  get_keyval(conf, "printTrajectoryFrequency", m_traj_output_frequency, cvm::cv_traj_freq);
  if (m_traj_output_frequency % time_step_factor != 0) {
    error_code |= cvm::error(
        "printTrajectoryFrequency (currently " + cvm::to_str(m_traj_output_frequency) +
            ") must be a multiple of timeStepFactor (" + cvm::to_str(time_step_factor) + ").\n",
        COLVARS_INPUT_ERROR);
  }
  m_cv.resize(num_variables(), 0);
  showInfo();
  return error_code;
}

void colvarbias_opes::showInfo() const {
  // Print information about this bias
  auto printInfo = [&](const std::string& info, const std::string& val){
    cvm::log(this->name + ": " + info + val + "\n");
  };
  printInfo("temperature = ", cvm::to_str(m_kbt / cvm::proxy->boltzmann()));
  printInfo("beta = ", cvm::to_str(1.0 / m_kbt));
  printInfo("depositing new kernels with newHillFrequency = ", cvm::to_str(m_pace));
  printInfo("expected barrier is ", cvm::to_str(m_barrier));
  printInfo("using target distribution with biasfactor (gamma) = ", m_inf_biasfactor ? "inf" : cvm::to_str(m_biasfactor));
  if (m_inf_biasfactor) {
    cvm::log("  (thus a uniform flat target distribution, no well-tempering)\n");
    cvm::log(this->name + ": " + "the equivalent bias temperature = inf\n");
  } else {
    cvm::log(this->name + ": " + "the equivalent bias temperature = " + cvm::to_str(cvm::proxy->target_temperature() * (m_biasfactor - 1)));
  }
  if (m_adaptive_sigma) {
    printInfo("adaptive sigma will be used, with adaptiveSigmaStride = ", cvm::to_str(m_adaptive_sigma_stride));
    size_t x = std::ceil(m_adaptive_sigma_stride / m_pace);
    printInfo("  thus the first x kernel depositions will be skipped, x = adaptiveSigmaStride/newHillFrequency = ", cvm::to_str(x));
  } else {
    std::string sigmas;
    for (size_t i = 0; i < num_variables(); ++i) {
      sigmas += " " + cvm::to_str(m_sigma0[i]);
    }
    cvm::log(this->name + ": kernels have initial gaussianSigma = " + sigmas + "\n");
  }
  if (m_fixed_sigma) {
    cvm::log(this->name + " fixedGaussianSigma: gaussianSigma will not decrease as the simulation proceeds\n");
  }
  printInfo("kernels are truncated with kernelCutoff = ", cvm::to_str(m_cutoff));
  if (m_cutoff < 3.5) {
    cvm::log(this->name + " +++ WARNING +++ probably kernels are truncated too much\n");
  }
  printInfo("the value at cutoff is = ", cvm::to_str(m_val_at_cutoff));
  printInfo("regularization epsilon = ", cvm::to_str(m_epsilon));
  if (m_val_at_cutoff > m_epsilon*(1+1e-6)) {
    cvm::log(this->name + " +++ WARNING +++ the kernelCutoff might be too small for the given epsilon\n");
  }
  printInfo("kernels will be compressed when closer than compression_threshold = ", cvm::to_str(m_compression_threshold));
  if (m_compression_threshold2 == 0) {
    cvm::log(this->name + " +++ WARNING +++ kernels will never merge, expect slowdowns\n");
  }
  if (!m_recursive_merge) {
    cvm::log(this->name + " -- RECURSIVE_MERGE_OFF: only one merge for each new kernel will be attempted. This is faster only if total number of kernels does not grow too much\n");
  }
  if (m_nlist) {
    cvm::log(this->name + " neighborList: using neighbor list for kernels, with parameters: " + cvm::to_str(m_nlist_param[0]) + " " + cvm::to_str(m_nlist_param[1]) + "\n");
    if (m_nlist_pace_reset) {
      cvm::log(this->name + " neighborListNewHillReset: forcing the neighbor list to update every time when depositing a new hill\n");
    }
  }
  if (m_no_zed) {
    printInfo("noZed: using fixed normalization factor = ", cvm::to_str(m_zed));
  }
  if (comm == multiple_replicas && m_num_walkers > 1) {
    cvm::log(this->name + " if multiple replicas are present, they will share the same bias\n");
  }
  if (m_num_threads > 1) {
    printInfo("using multiple threads per simulation: ", cvm::to_str(m_num_threads));
  }
  cvm::main()->cite_feature("OPES");
  if (m_adaptive_sigma || m_explore) {
    cvm::main()->cite_feature("OPES explore or adaptive kernels");
  }
}

cvm::real colvarbias_opes::evaluateKernel(
  const colvarbias_opes::kernel& G,
  const std::vector<cvm::real>& x) const {
  cvm::real norm2 = 0;
  for (size_t i = 0; i < num_variables(); ++i) {
    const cvm::real dist2_i = variables(i)->dist2(G.m_center[i], x[i]) / (G.m_sigma[i] * G.m_sigma[i]);
    norm2 += dist2_i;
    if (norm2 >= m_cutoff2) {
      return 0;
    }
  }
  return G.m_height * (std::exp(-0.5 * norm2) - m_val_at_cutoff);
}

cvm::real colvarbias_opes::evaluateKernel(
  const colvarbias_opes::kernel& G,
  const std::vector<cvm::real>& x,
  std::vector<cvm::real>& accumulated_derivative,
  std::vector<cvm::real>& dist) const {
  cvm::real norm2 = 0;
  for (size_t i = 0; i < num_variables(); ++i) {
    dist[i] = 0.5 * variables(i)->dist2_lgrad(x[i], G.m_center[i]) / G.m_sigma[i];
    norm2 += dist[i] * dist[i];
    if (norm2 >= m_cutoff2) {
      return 0;
    }
  }
  const cvm::real val = G.m_height * (std::exp(-0.5 * norm2) - m_val_at_cutoff);
  // The derivative of norm2 with respect to x
  for (size_t i = 0; i < num_variables(); ++i) {
    accumulated_derivative[i] -= val * dist[i] / G.m_sigma[i];
  }
  return val;
}

cvm::real colvarbias_opes::getProbAndDerivatives(
  const std::vector<cvm::real>& cv, std::vector<cvm::real>& der_prob) const {
  double prob = 0.0;
  std::vector<cvm::real> dist(num_variables(), 0);
  if (!m_nlist) {
    if (m_num_threads == 1 || m_kernels.size() < 2 * m_num_threads) {
      for (size_t k = 0; k < m_kernels.size(); ++k) {
        prob += evaluateKernel(m_kernels[k], cv, der_prob, dist);
      }
    } else {
#if defined(_OPENMP)
      #pragma omp parallel num_threads(m_num_threads)
      {
        std::vector<cvm::real> omp_deriv(der_prob.size(), 0);
        std::vector<cvm::real> tmp_dist(num_variables());
        #pragma omp for reduction(+:prob) nowait
        for (int k = 0; k < static_cast<int>(m_kernels.size()); ++k) {
          prob += evaluateKernel(m_kernels[k], cv, omp_deriv, tmp_dist);
        }
        #pragma omp critical
        for (int i = 0; i < static_cast<int>(num_variables()); ++i) {
          der_prob[i]+=omp_deriv[i];
        }
        #pragma omp single
        for (int i = 0; i < static_cast<int>(num_variables()); ++i) {
          dist[i] = tmp_dist[i];
        }
      }
#elif CMK_SMP && USE_CKLOOP
      // TODO: Test this once fine-grained parallelization is enabled
      std::vector<std::vector<cvm::real>> derivs(m_num_threads, std::vector<cvm::real>(num_variables(), 0));
      std::vector<std::vector<cvm::real>> dists(m_num_threads, std::vector<cvm::real>(num_variables(), 0));
      auto worker = [&](int start, int end, void* result){
        const int tid = cvm::proxy->smp_thread_id();
        double tmp_prob = 0;
        for (int i = start; i <= end; ++i) {
          tmp_prob += evaluateKernel(m_kernels[i], cv, derivs[tid], dists[tid]);
        }
        *(double *)result = tmp_prob;
      };
      const size_t lowerRange = 0;
      const size_t upperRange = m_kernels.size() - 1;
      CkLoop_Parallelize(
        m_num_threads, lowerRange, upperRange,
        worker, &prob, CKLOOP_DOUBLE_SUM, NULL);
      for (size_t i = 0; i < num_variables(); ++i) {
        for (size_t j = 0; j < m_num_threads; ++j) {
          if (j == 0) dist[i] = dists[j][i];
          der_prob[i] += derivs[j][i];
        }
      }
#else
      cvm::error("multiple threads required in OPES, but this binary is not linked with a supported threading library.\n");
#endif
    }
  } else {
    if (m_num_threads == 1 || m_nlist_index.size() < 2 * m_num_threads) {
      for (size_t nk = 0; nk < m_nlist_index.size(); ++nk) {
        const size_t k = m_nlist_index[nk];
        prob += evaluateKernel(m_kernels[k], cv, der_prob, dist);
      }
    } else {
#if defined(_OPENMP)
      #pragma omp parallel num_threads(m_num_threads)
      {
        std::vector<cvm::real> omp_deriv(der_prob.size(), 0);
        std::vector<cvm::real> tmp_dist(num_variables());
        #pragma omp for reduction(+:prob) nowait
        for (int nk = 0; nk < static_cast<int>(m_nlist_index.size()); ++nk) {
          const size_t k = m_nlist_index[nk];
          prob += evaluateKernel(m_kernels[k], cv, omp_deriv, tmp_dist);
        }
        #pragma omp critical
        for (int i = 0; i < static_cast<int>(num_variables()); ++i) {
          der_prob[i]+=omp_deriv[i];
        }
        #pragma omp single
        for (int i = 0; i < static_cast<int>(num_variables()); ++i) {
          dist[i] = tmp_dist[i];
        }
      }
#elif CMK_SMP && USE_CKLOOP
      // TODO: Test this once fine-grained parallelization is enabled
      std::vector<std::vector<cvm::real>> derivs(m_num_threads, std::vector<cvm::real>(num_variables(), 0));
      std::vector<std::vector<cvm::real>> dists(m_num_threads, std::vector<cvm::real>(num_variables(), 0));
      auto worker = [&](int start, int end, void* result){
        const int tid = cvm::proxy->smp_thread_id();
        double tmp_prob = 0;
        for (int i = start; i <= end; ++i) {
          const size_t k = m_nlist_index[i];
          tmp_prob += evaluateKernel(m_kernels[k], cv, derivs[tid], dists[tid]);
        }
        *(double *)result = tmp_prob;
      };
      const size_t lowerRange = 0;
      const size_t upperRange = m_nlist_index.size() - 1;
      CkLoop_Parallelize(
        m_num_threads, lowerRange, upperRange,
        worker, &prob, CKLOOP_DOUBLE_SUM, NULL);
      for (size_t i = 0; i < num_variables(); ++i) {
        for (size_t j = 0; j < m_num_threads; ++j) {
          if (j == 0) dist[i] = dists[j][i];
          der_prob[i] += derivs[j][i];
        }
      }
#else
      cvm::error("multiple threads required in OPES, but this binary is not linked with a supported threading library.\n");
#endif
    }
  }
  prob /= m_kdenorm;
  for (size_t i = 0; i < num_variables(); ++i) {
    der_prob[i] /= m_kdenorm;
  }
  return prob;
}

int colvarbias_opes::calculate_opes() {
  if (m_nlist) {
    ++m_nlist_steps;
    const bool exchange_step =
      (comm == multiple_replicas) &&
      cvm::step_absolute() % shared_freq == 0;
    if (exchange_step) {
      m_nlist_update = true;
    } else {
      for (size_t i = 0; i < num_variables(); ++i) {
        const cvm::real diff_i2 = variables(i)->dist2(m_cv[i], m_nlist_center[i]);
        if (diff_i2 > m_nlist_param[1] * m_nlist_dev2[i]) {
          m_nlist_update = true;
          break;
        }
      }
    }
    if (m_nlist_update) {
      updateNlist(m_cv);
    }
  }
  std::vector<cvm::real> der_prob(num_variables(), 0);
  const cvm::real prob = getProbAndDerivatives(m_cv, der_prob);
  const cvm::real bias = m_kbt * m_bias_prefactor * cvm::logn(prob / m_zed + m_epsilon);
  bias_energy = bias;
  if (is_enabled(f_cvb_apply_force)) {
    for (size_t i = 0; i < num_variables(); ++i) {
      colvar_forces[i] = -m_kbt * m_bias_prefactor / (prob / m_zed + m_epsilon) * der_prob[i] / m_zed;
    }
  }
  return COLVARS_OK;
}

int colvarbias_opes::update_opes() {
  if (m_adaptive_sigma) {
    m_adaptive_counter++;
    cvm::step_number tau = m_adaptive_sigma_stride;
    if (m_adaptive_counter < m_adaptive_sigma_stride) tau = m_adaptive_counter;
    for (size_t i = 0; i < num_variables(); ++i) {
      // Welford's online algorithm for standard deviation
      const cvm::real diff_i = 0.5 * variables(i)->dist2_lgrad(m_cv[i], m_av_cv[i]);
      m_av_cv[i] += diff_i / tau;
      m_av_M2[i] += diff_i * 0.5 * variables(i)->dist2_lgrad(m_cv[i], m_av_cv[i]);
    }
    if (m_adaptive_counter < m_adaptive_sigma_stride && m_counter == 1) {
      return COLVARS_OK;;
    }
  }
  if (cvm::step_absolute() % m_pace == 0) {
    m_old_kdenorm = m_kdenorm;
    m_delta_kernels.clear();
    const size_t old_nker = m_kernels.size();
    // TODO: how could I account for extra biases in Colvars?
    const cvm::real log_weight = bias_energy / m_kbt;
    cvm::real height = cvm::exp(log_weight);
    cvm::real sum_heights = height;
    cvm::real sum_heights2 = height * height;
    if (m_num_walkers > 1) {
      std::vector<cvm::real> replica_sum_heights(cvm::proxy->num_replicas() - 1, 0);
      // Send all sum_heights to PE 0
      if (cvm::proxy->replica_index() == 0) {
        for (int p = 1; p < cvm::proxy->num_replicas(); ++p) {
          if (cvm::proxy->replica_comm_recv((char*)&(replica_sum_heights[p - 1]), sizeof(cvm::real), p) != sizeof(cvm::real)) {
            return cvm::error("Error: receiving sum of weights from replica " + cvm::to_str(p));
          }
        }
      } else {
        if (cvm::proxy->replica_comm_send((char*)&sum_heights, sizeof(cvm::real), 0) != sizeof(cvm::real)) {
          return cvm::error("Error: sending sum of weights to replica 0.");
        }
      }
      cvm::proxy->replica_comm_barrier();
      // PE 0 sum all sum_heights and broadcast
      if (cvm::proxy->replica_index() == 0) {
        for (auto it = replica_sum_heights.begin(); it != replica_sum_heights.end(); ++it) {
          sum_heights += (*it);
        }
        for (int p = 1; p < cvm::proxy->num_replicas(); ++p) {
          if (cvm::proxy->replica_comm_send((char*)&sum_heights, sizeof(cvm::real), p) != sizeof(cvm::real)) {
            return cvm::error("Error: sending sum of weights to replica " + cvm::to_str(p));
          }
        }
      } else {
        if (cvm::proxy->replica_comm_recv((char*)&sum_heights, sizeof(cvm::real), 0) != sizeof(cvm::real)) {
          return cvm::error("Error: receiving sum of weights from replica 0.");
        }
      }
      cvm::proxy->replica_comm_barrier();
      // Send all sum_heights2 to PE 0
      std::vector<cvm::real> replica_sum_heights2(cvm::proxy->num_replicas() - 1, 0);
      if (cvm::proxy->replica_index() == 0) {
        for (int p = 1; p < cvm::proxy->num_replicas(); ++p) {
          if (cvm::proxy->replica_comm_recv((char*)&(replica_sum_heights2[p - 1]), sizeof(cvm::real), p) != sizeof(cvm::real)) {
            return cvm::error("Error: getting sum of weights2 from replica " + cvm::to_str(p));
          }
        }
      } else {
        if (cvm::proxy->replica_comm_send((char*)&sum_heights2, sizeof(cvm::real), 0) != sizeof(cvm::real)) {
          return cvm::error("Error: sending sum of weights2 from replica.");
        }
      }
      cvm::proxy->replica_comm_barrier();
      // PE 0 sum all sum_heights2 and broadcast
      if (cvm::proxy->replica_index() == 0) {
        for (auto it = replica_sum_heights2.begin(); it != replica_sum_heights2.end(); ++it) {
          sum_heights2 += (*it);
        }
        for (int p = 1; p < cvm::proxy->num_replicas(); ++p) {
          if (cvm::proxy->replica_comm_send((char*)&sum_heights2, sizeof(cvm::real), p) != sizeof(cvm::real)) {
            return cvm::error("Error: sending sum of weights2 to replica " + cvm::to_str(p));
          }
        }
      } else {
        if (cvm::proxy->replica_comm_recv((char*)&sum_heights2, sizeof(cvm::real), 0) != sizeof(cvm::real)) {
          return cvm::error("Error: receiving sum of weights2 from replica.");
        }
      }
      cvm::proxy->replica_comm_barrier();
    }
    m_counter += m_num_walkers;
    m_sum_weights += sum_heights;
    m_sum_weights2 += sum_heights2;
    m_neff = (1 + m_sum_weights) * (1 + m_sum_weights) / (1 + m_sum_weights2);
    m_rct = m_kbt * cvm::logn(m_sum_weights / m_counter);
    m_traj_line.neff = m_neff;
    m_traj_line.rct = m_rct;
    if (m_explore) {
      m_kdenorm = m_counter;
      height = 1.0;
    } else {
      m_kdenorm = m_sum_weights;
    }
    std::vector<cvm::real> sigma = m_sigma0;
    if (m_adaptive_sigma) {
      const cvm::real factor = m_explore ? 1.0 : m_biasfactor;
      if (m_counter == 1 + m_num_walkers) {
        for (size_t i = 0; i < num_variables(); ++i) {
          m_av_M2[i] *= m_biasfactor;
        }
        for (size_t i = 0; i < num_variables(); ++i) {
          m_sigma0[i] = std::sqrt(m_av_M2[i] / m_adaptive_counter / factor);
        }
        if (m_sigma_min.size() == 0) {
          for (size_t i = 0; i < num_variables(); ++i) {
            if (m_sigma0[i] < 1e-6) {
              cvm::error("Adaptive sigma is suspiciously small for CV " + cvm::to_str(i) + "\nManually provide sigma or set a safe sigma_min to avoid possible issues\n");
              return COLVARS_ERROR;
            }
          }
        } else {
          for (size_t i = 0; i < num_variables(); ++i) {
            m_sigma0[i] = std::max(m_sigma0[i], m_sigma_min[i]);
          }
        }
      }
      for (size_t i = 0; i < num_variables(); ++i) {
        sigma[i] = std::sqrt(m_av_M2[i] / m_adaptive_counter / factor);
      }
      if (m_sigma_min.size() == 0) {
        bool sigma_less_than_threshold = false;
        for (size_t i = 0; i < num_variables(); ++i) {
          if (sigma[i] < 1e-6) {
            cvm::log("The adaptive sigma is suspiciously small, you should set a safe sigma_min. 1e-6 will be used here\n");
            sigma[i] = 1e-6;
            sigma_less_than_threshold = true;
          }
        }
        if (sigma_less_than_threshold) {
          m_sigma_min.assign(num_variables(), 1e-6);
        }
      } else {
        for (size_t i = 0; i < num_variables(); ++i) {
          sigma[i] = std::max(sigma[i], m_sigma_min[i]);
        }
      }
    }
    if (!m_fixed_sigma) {
      const cvm::real size = m_explore ? m_counter : m_neff;
      const size_t ncv = num_variables();
      const cvm::real s_rescaling = std::pow(size * (ncv + 2.0) / 4, -1.0 / (4.0 + ncv));
      for (size_t i = 0; i < num_variables(); ++i) {
        sigma[i] *= s_rescaling;
      }
      if (m_sigma_min.size() > 0) {
        for (size_t i = 0; i < num_variables(); ++i) {
          sigma[i] = std::max(sigma[i], m_sigma_min[i]);
        }
      }
    }
    // the height should be divided by sqrt(2*pi)*sigma0_,
    // but this overall factor would be canceled when dividing by Zed
    // thus we skip it altogether, but keep any other sigma rescaling
    for (size_t i = 0; i < num_variables(); ++i) {
      height *= (m_sigma0[i] / sigma[i]);
    }
    if (m_num_walkers == 1) {
      addKernel(height, m_cv, sigma, log_weight);
    } else {
      std::vector<cvm::real> all_height(m_num_walkers, 0.0);
      std::vector<cvm::real> all_center(m_num_walkers * num_variables(), 0.0);
      std::vector<cvm::real> all_sigma(m_num_walkers * num_variables(), 0.0);
      std::vector<cvm::real> all_logweight(m_num_walkers, 0.0);
      const int my_replica = cvm::proxy->replica_index();

      // Allgather of heights
      if (my_replica == 0) {
        all_height[0] = height;
        for (int p = 1; p < cvm::proxy->num_replicas(); ++p) {
          if (cvm::proxy->replica_comm_recv((char*)&(all_height[p]), sizeof(decltype(all_height)::value_type), p) != sizeof(decltype(all_height)::value_type)) {
            return cvm::error("Error: on receiving height on replica 0 from replica " + cvm::to_str(p));
          }
        }
      } else {
        if (cvm::proxy->replica_comm_send((char*)&height, sizeof(decltype(height)), 0) != sizeof(cvm::real)) {
          return cvm::error("Error: on sending height to replica 0 from replica " + cvm::to_str(my_replica));
        }
      }
      cvm::proxy->replica_comm_barrier();
      // Broadcast heights
      if (my_replica == 0) {
        const int send_size = sizeof(decltype(all_height)::value_type) * all_height.size();
        for (int p = 1; p < cvm::proxy->num_replicas(); ++p) {
          if (cvm::proxy->replica_comm_send((char*)all_height.data(), send_size, p) != send_size) {
            return cvm::error("Error: on sending heights from replica 0 to replica " + cvm::to_str(p));
          }
        }
      } else {
        const int recv_size = sizeof(decltype(all_height)::value_type) * all_height.size();
        if (cvm::proxy->replica_comm_recv((char*)all_height.data(), recv_size, 0) != recv_size) {
          return cvm::error("Error: on receiving heights from replica 0 to replica " + cvm::to_str(my_replica));
        }
      }
      cvm::proxy->replica_comm_barrier();

      // Allgather of centers
      if (my_replica == 0) {
        std::copy(m_cv.begin(), m_cv.end(), all_center.begin());
        const int recv_size = sizeof(decltype(m_cv)::value_type) * m_cv.size();
        for (int p = 1; p < cvm::proxy->num_replicas(); ++p) {
          cvm::real* recv_start_ptr = &(all_center[p * m_cv.size()]);
          if (cvm::proxy->replica_comm_recv((char*)recv_start_ptr, recv_size, p) != recv_size) {
            return cvm::error("Error on receiving centers from replica 0 to replica " + cvm::to_str(p));
          }
        }
      } else {
        const int send_size = sizeof(decltype(m_cv)::value_type) * m_cv.size();
        if (cvm::proxy->replica_comm_send((char*)m_cv.data(), send_size, 0) != send_size) {
          return cvm::error("Error on sending centers to replica 0 from replica " + cvm::to_str(my_replica));
        }
      }
      cvm::proxy->replica_comm_barrier();
      // Broadcast centers
      if (my_replica == 0) {
        const int send_size = sizeof(decltype(all_center)::value_type) * all_center.size();
        for (int p = 1; p < cvm::proxy->num_replicas(); ++p) {
          if (cvm::proxy->replica_comm_send((char*)all_center.data(), send_size, p) != send_size) {
            return cvm::error("Error on sending centers from replica 0 to replica " + cvm::to_str(p));
          }
        }
      } else {
        const int recv_size = sizeof(decltype(all_center)::value_type) * all_center.size();
        if (cvm::proxy->replica_comm_recv((char*)all_center.data(), recv_size, 0) != recv_size) {
          return cvm::error("Error on receiving centers from replica 0 to replica " + cvm::to_str(my_replica));
        }
      }
      cvm::proxy->replica_comm_barrier();

      // Allgather of sigmas
      if (my_replica == 0) {
        std::copy(sigma.begin(), sigma.end(), all_sigma.begin());
        const int recv_size = sizeof(decltype(sigma)::value_type) * sigma.size();
        for (int p = 1; p < cvm::proxy->num_replicas(); ++p) {
          cvm::real* recv_start_ptr = &(all_sigma[p * m_cv.size()]);
          if (cvm::proxy->replica_comm_recv((char*)recv_start_ptr, recv_size, p) != recv_size) {
            return cvm::error("Error on receiving sigmas from replica 0 to replica " + cvm::to_str(p));
          }
        }
      } else {
        const int send_size = sizeof(decltype(sigma)::value_type) * sigma.size();
        if (cvm::proxy->replica_comm_send((char*)sigma.data(), send_size, 0) != send_size) {
          return cvm::error("Error on sending sigmas to replica 0 from replica " + cvm::to_str(my_replica));
        }
      }
      cvm::proxy->replica_comm_barrier();
      // Broadcast sigmas
      if (my_replica == 0) {
        const int send_size = sizeof(decltype(all_sigma)::value_type) * all_sigma.size();
        for (int p = 1; p < cvm::proxy->num_replicas(); ++p) {
          if (cvm::proxy->replica_comm_send((char*)all_sigma.data(), send_size, p) != send_size) {
            return cvm::error("Error on sending sigmas from replica 0 to replica " + cvm::to_str(p));
          }
        }
      } else {
        const int recv_size = sizeof(decltype(all_sigma)::value_type) * all_sigma.size();
        if (cvm::proxy->replica_comm_recv((char*)all_sigma.data(), recv_size, 0) != recv_size) {
          return cvm::error("Error on receiving sigmas from replica 0 to replica " + cvm::to_str(my_replica));
        }
      }
      cvm::proxy->replica_comm_barrier();

      // Allgather of logweights
      if (my_replica == 0) {
        all_logweight[0] = log_weight;
        for (int p = 1; p < cvm::proxy->num_replicas(); ++p) {
          if (cvm::proxy->replica_comm_recv((char*)&(all_logweight[p]), sizeof(decltype(all_logweight)::value_type), p) != sizeof(decltype(all_logweight)::value_type)) {
            return cvm::error("Error on receiving log_weight on replica 0 from replica " + cvm::to_str(p));
          }
        }
      } else {
        if (cvm::proxy->replica_comm_send((char*)&log_weight, sizeof(decltype(log_weight)), 0) != sizeof(cvm::real)) {
          return cvm::error("Error on sending log_weight to replica 0 from replica " + cvm::to_str(my_replica));
        }
      }
      cvm::proxy->replica_comm_barrier();
      // Broadcast log_weight
      if (my_replica == 0) {
        const int send_size = sizeof(decltype(all_logweight)::value_type) * all_logweight.size();
        for (int p = 1; p < cvm::proxy->num_replicas(); ++p) {
          if (cvm::proxy->replica_comm_send((char*)all_logweight.data(), send_size, p) != send_size) {
            return cvm::error("Error on sending log_weight from replica 0 to replica " + cvm::to_str(p));
          }
        }
      } else {
        const int recv_size = sizeof(decltype(all_logweight)::value_type) * all_logweight.size();
        if (cvm::proxy->replica_comm_recv((char*)all_logweight.data(), recv_size, 0) != recv_size) {
          return cvm::error("Error on receiving log_weight from replica 0 to replica " + cvm::to_str(my_replica));
        }
      }
      cvm::proxy->replica_comm_barrier();

      if (m_nlist) {
        std::vector<int> all_nlist_size(m_num_walkers);
        const int my_replica = cvm::proxy->replica_index();
        // Get the size of the neighbor list of each replica
        if (my_replica == 0) {
          all_nlist_size[0] = m_nlist_index.size();
          for (int p = 1; p < cvm::proxy->num_replicas(); ++p) {
            if (cvm::proxy->replica_comm_recv((char*)&(all_nlist_size[p]), sizeof(int), p) != sizeof(int)) {
              return cvm::error("Error on receiving neighbor list size from replica " + cvm::to_str(p));
            }
          }
        } else {
          const int nlist_size = m_nlist_index.size();
          if (cvm::proxy->replica_comm_send((char*)&nlist_size, sizeof(int), 0) != sizeof(int)) {
            return cvm::error("Error on sending neighbor list size from replica " + cvm::to_str(my_replica));
          }
        }
        cvm::proxy->replica_comm_barrier();
        // Broadcast the neighbor list sizes to all replicas
        if (my_replica == 0) {
          const int send_size = sizeof(int) * all_nlist_size.size();
          for (int p = 1; p < cvm::proxy->num_replicas(); ++p) {
            if (cvm::proxy->replica_comm_send((char*)all_nlist_size.data(), send_size, p) != send_size) {
              return cvm::error("Error on sending neighbor list sizes from replica 0 to replica " + cvm::to_str(p));
            }
          }
        } else {
          const int recv_size = sizeof(int) * all_nlist_size.size();
          if (cvm::proxy->replica_comm_recv((char*)all_nlist_size.data(), recv_size, 0) != recv_size) {
            return cvm::error("Error on receiving neighbor list sizes to replica " + cvm::to_str(my_replica));
          }
        }
        cvm::proxy->replica_comm_barrier();
        // Gather all neighbor lists from replicas
        const int tot_size = std::accumulate(all_nlist_size.begin(), all_nlist_size.end(), 0);
        if (tot_size > 0) {
          // Allgatherv all neighbor lists from replicas
          std::vector<size_t> all_nlist_index(tot_size);
          if (my_replica == 0) {
            std::vector<int> recv_start(m_num_walkers);
            // Accumulative sum
            recv_start[0] = 0;
            std::partial_sum(all_nlist_size.begin(), all_nlist_size.end() - 1, recv_start.begin() + 1);
            std::copy(m_nlist_index.begin(), m_nlist_index.end(), all_nlist_index.begin());
            for (int p = 1; p < cvm::proxy->num_replicas(); ++p) {
              size_t* recv_start_ptr = &(all_nlist_index[recv_start[p]]);
              const int recv_size = all_nlist_size[p] * sizeof(decltype(all_nlist_index)::value_type);
              if (cvm::proxy->replica_comm_recv((char*)recv_start_ptr, recv_size, p) != recv_size) {
                return cvm::error("Error on receiving neighbor list from replica " + cvm::to_str(p));
              }
            }
          } else {
            const int send_size = sizeof(decltype(m_nlist_index)::value_type) * m_nlist_index.size();
            if (cvm::proxy->replica_comm_send((char*)m_nlist_index.data(), send_size, 0) != send_size) {
              return cvm::error("Error on sending neighbor list from replica " + cvm::to_str(my_replica));
            }
          }
          cvm::proxy->replica_comm_barrier();
          // Broadcast the neighbor list
          if (my_replica == 0) {
            const int send_size = sizeof(decltype(all_nlist_index)::value_type) * tot_size;
            for (int p = 1; p < cvm::proxy->num_replicas(); ++p) {
              if (cvm::proxy->replica_comm_send((char*)all_nlist_index.data(), send_size, p) != send_size) {
                return cvm::error("Error on sending total neighbor list to replica " + cvm::to_str(p));
              }
            }
          } else {
            const int recv_size = sizeof(decltype(all_nlist_index)::value_type) * tot_size;
            if (cvm::proxy->replica_comm_recv((char*)all_nlist_index.data(), recv_size, 0) != recv_size) {
              return cvm::error("Error on receiving total neighbor list on replica " + cvm::to_str(my_replica));
            }
          }
          cvm::proxy->replica_comm_barrier();
          // Deduplicate and sort the merged neighbor list
          std::unordered_set<size_t> all_nlist_index_set;
          for (auto it = all_nlist_index.cbegin(); it != all_nlist_index.cend(); ++it) {
            all_nlist_index_set.insert(*it);
          }
          m_nlist_index.assign(all_nlist_index_set.begin(), all_nlist_index_set.end());
          std::sort(m_nlist_index.begin(), m_nlist_index.end());
        }
      }
      for (size_t w = 0; w < m_num_walkers; ++w) {
        std::vector<cvm::real> center_w(
          all_center.begin() + num_variables() * w,
          all_center.begin() + num_variables() * (w + 1));
        std::vector<cvm::real> sigma_w(
          all_sigma.begin() + num_variables() * w,
          all_sigma.begin() + num_variables() * (w + 1));
        addKernel(all_height[w], center_w, sigma_w, all_logweight[w]);
      }
    }
    m_nker = m_kernels.size();
    m_traj_line.nker = m_nker;
    if (m_nlist) {
      m_nlker = m_nlist_index.size();
      m_traj_line.nlker = m_nlker;
      if (m_nlist_pace_reset) {
        m_nlist_update = true;
      }
    }
    if (!m_no_zed) {
      cvm::real sum_uprob = 0;
      const size_t ks = m_kernels.size();
      const size_t ds = m_delta_kernels.size();
      const int num_parallel = 1; // Always 1
      const bool few_kernels = (ks * ks < (3 * ks * ds + 2 * ds * ds * num_parallel + 100));
      if (few_kernels) {
        if (m_num_threads == 1) {
          for (size_t k = 0; k < m_kernels.size(); ++k) {
            for (size_t kk = 0; kk < m_kernels.size(); ++kk) {
              sum_uprob += evaluateKernel(m_kernels[kk], m_kernels[k].m_center);
            }
          }
        } else {
#if defined(_OPENMP)
          #pragma omp parallel num_threads(m_num_threads)
          {
            #pragma omp for reduction(+:sum_uprob) nowait
            for (int k = 0; k < static_cast<int>(m_kernels.size()); ++k) {
              for (int kk = 0; kk < static_cast<int>(m_kernels.size()); ++kk) {
                sum_uprob += evaluateKernel(m_kernels[kk], m_kernels[k].m_center);
              }
            }
          }
#elif CMK_SMP && USE_CKLOOP
          // TODO: Does this work??
          auto worker = [&](int start, int end, void* result) {
            double tmp_prob = 0;
            for (int i = start; i <= end; ++i) {
              for (size_t kk = 0; kk < m_kernels.size(); ++kk) {
                tmp_prob += evaluateKernel(m_kernels[kk], m_kernels[i].m_center);
              }
            }
            *(double *)result = tmp_prob;
          };
          const size_t lowerRange = 0;
          const size_t upperRange = m_kernels.size() - 1;
          CkLoop_Parallelize(
            m_num_threads, lowerRange, upperRange,
            worker, &sum_uprob, CKLOOP_DOUBLE_SUM, NULL);
#else
          cvm::error("OPES cannot run because this binary is not linked with a supported threading library.\n");
#endif
        }
        if (num_parallel > 1) {
          return cvm::error("Unimplemented feature: OPES in parallel running.\n");
        }
      } else {
        cvm::real delta_sum_uprob = 0;
        if (!m_nlist) {
          if (m_num_threads == 1) {
            for (size_t i = 0; i < m_kernels.size(); ++i) {
              for (size_t d = 0; d < m_delta_kernels.size(); ++d) {
                const int sign = m_delta_kernels[d].m_height < 0 ? -1 : 1;
                delta_sum_uprob += evaluateKernel(m_delta_kernels[d], m_kernels[i].m_center) + sign * evaluateKernel(m_kernels[i], m_delta_kernels[d].m_center);
              }
            }
          } else {
#if defined(_OPENMP)
            #pragma omp parallel num_threads(m_num_threads)
            {
              #pragma omp for reduction(+:delta_sum_uprob) nowait
              for (int i = 0; i < static_cast<int>(m_kernels.size()); ++i) {
                for (int d = 0; d < static_cast<int>(m_delta_kernels.size()); ++d) {
                  const int sign = m_delta_kernels[d].m_height < 0 ? -1 : 1;
                  delta_sum_uprob += evaluateKernel(m_delta_kernels[d], m_kernels[i].m_center) + sign * evaluateKernel(m_kernels[i], m_delta_kernels[d].m_center);
                }
              }
            }
#elif CMK_SMP && USE_CKLOOP
            auto worker = [&](int start, int end, void* result) {
              double tmp_prob = 0;
              for (int i = start; i <= end; ++i) {
                for (size_t d = 0; d < m_delta_kernels.size(); ++d) {
                  const int sign = m_delta_kernels[d].m_height < 0 ? -1 : 1;
                  tmp_prob += evaluateKernel(m_delta_kernels[d], m_kernels[i].m_center) + sign * evaluateKernel(m_kernels[i], m_delta_kernels[d].m_center);
                }
              }
              *(double *)result = tmp_prob;
            };
            const size_t lowerRange = 0;
            const size_t upperRange = m_kernels.size() - 1;
            CkLoop_Parallelize(
              m_num_threads, lowerRange, upperRange,
              worker, &delta_sum_uprob, CKLOOP_DOUBLE_SUM, NULL);
#else
            cvm::error("OPES cannot run because this binary is not linked with a supported threading library.\n");
#endif
          }
        } else {
          if (m_num_threads == 1) {
            for (size_t i = 0; i < m_nlist_index.size(); ++i) {
              const size_t k = m_nlist_index[i];
              for (size_t d = 0; d < m_delta_kernels.size(); ++d) {
                const double sign = m_delta_kernels[d].m_height < 0 ? -1 : 1;
                delta_sum_uprob += evaluateKernel(m_delta_kernels[d], m_kernels[k].m_center) + sign * evaluateKernel(m_kernels[k], m_delta_kernels[d].m_center);
              }
            }
          } else {
#if defined(_OPENMP)
            #pragma omp parallel num_threads(m_num_threads)
            {
              #pragma omp for reduction(+:delta_sum_uprob) nowait
              for (int i = 0; i < static_cast<int>(m_nlist_index.size()); ++i) {
                const size_t k = m_nlist_index[i];
                for (int d = 0; d < static_cast<int>(m_delta_kernels.size()); ++d) {
                  const double sign = m_delta_kernels[d].m_height < 0 ? -1 : 1;
                  delta_sum_uprob += evaluateKernel(m_delta_kernels[d], m_kernels[k].m_center) + sign * evaluateKernel(m_kernels[k], m_delta_kernels[d].m_center);
                }
              }
            }
#elif CMK_SMP && USE_CKLOOP
            auto worker = [&](int start, int end, void* result) {
              double tmp_prob = 0;
              for (int i = start; i <= end; ++i) {
                const size_t k = m_nlist_index[i];
                for (size_t d = 0; d < m_delta_kernels.size(); ++d) {
                  const double sign = m_delta_kernels[d].m_height < 0 ? -1 : 1;
                  tmp_prob += evaluateKernel(m_delta_kernels[d], m_kernels[k].m_center) + sign * evaluateKernel(m_kernels[k], m_delta_kernels[d].m_center);
                }
              }
              *(double *)result = tmp_prob;
            };
            const size_t lowerRange = 0;
            const size_t upperRange = m_nlist_index.size() - 1;
            CkLoop_Parallelize(
              m_num_threads, lowerRange, upperRange,
              worker, &delta_sum_uprob, CKLOOP_DOUBLE_SUM, NULL);
#else
            cvm::error("OPES cannot run because this binary is not linked with a supported threading library.\n");
#endif
          }
        }
        if (num_parallel > 1) {
          return cvm::error("Unimplemented feature: OPES in parallel running.\n");
        }
        if (m_num_threads == 1) {
          for (size_t d = 0; d < m_delta_kernels.size(); ++d) {
            for (size_t dd = 0; dd < m_delta_kernels.size(); ++dd) {
              const int sign = m_delta_kernels[d].m_height < 0 ? -1 : 1;
              delta_sum_uprob -= sign *evaluateKernel(m_delta_kernels[dd], m_delta_kernels[d].m_center);
            }
          }
        } else {
#if defined(_OPENMP)
          #pragma omp parallel num_threads(m_num_threads)
          {
            #pragma omp for reduction(+:delta_sum_uprob)
            for (int d = 0; d < static_cast<int>(m_delta_kernels.size()); ++d) {
              for (int dd = 0; dd < static_cast<int>(m_delta_kernels.size()); ++dd) {
                const int sign = m_delta_kernels[d].m_height < 0 ? -1 : 1;
                delta_sum_uprob -= sign * evaluateKernel(m_delta_kernels[dd], m_delta_kernels[d].m_center);
              }
            }
          }
#elif CMK_SMP && USE_CKLOOP
          auto worker = [&](int start, int end, void* result) {
            double tmp_prob = 0;
            for (int d = start; d <= end; ++d) {
              for (size_t dd = 0; dd < m_delta_kernels.size(); ++dd) {
                const int sign = m_delta_kernels[d].m_height < 0 ? -1 : 1;
                tmp_prob += sign * evaluateKernel(m_delta_kernels[dd], m_delta_kernels[d].m_center);
              }
            }
            *(double *)result = tmp_prob;
          };
          const size_t lowerRange = 0;
          const size_t upperRange = m_delta_kernels.size() - 1;
          double tmp = 0;
          CkLoop_Parallelize(
            m_num_threads, lowerRange, upperRange,
            worker, &tmp, CKLOOP_DOUBLE_SUM, NULL);
          delta_sum_uprob -= tmp;
#else
          cvm::error("OPES cannot run because this binary is not linked with a supported threading library.\n");
#endif
        }
        sum_uprob = m_zed * m_old_kdenorm * old_nker + delta_sum_uprob;
      }
      m_zed = sum_uprob / m_kdenorm / m_kernels.size();
      m_traj_line.zed = m_zed;
    }
    if (m_calc_work) {
      std::vector<cvm::real> dummy(num_variables());
      const cvm::real prob = getProbAndDerivatives(m_cv, dummy);
      const cvm::real new_bias = m_kbt * m_bias_prefactor * cvm::logn(prob / m_zed + m_epsilon);
      m_work += new_bias - bias_energy;
      m_traj_line.work = m_work;
    }
  }
  return COLVARS_OK;
}

void colvarbias_opes::save_state() {
  if (cvm::step_absolute() % cvm::restart_out_freq == 0) {
    m_saved_zed = m_zed;
    m_saved_sum_weights = m_sum_weights;
    m_saved_sum_weights2 = m_sum_weights2;
    m_saved_kernels = m_kernels;
  }
}

int colvarbias_opes::update() {
  int error_code = COLVARS_OK;
  for (size_t i = 0; i < num_variables(); ++i) {
    m_cv[i] = variables(i)->value();
  }
  error_code |= calculate_opes();
  // NOTE: I don't think that calling dumpStateToFile() after update in
  //       the PLUMED implementation is correct for step 0, so I save the
  //       data after calculate() that does not modify the internal state
  //       of the bias.
  save_state();
  if (error_code != COLVARS_OK) return error_code;
  if (m_is_first_step) {
    // NOTE: Colvars does not allow chainned biases, so I have to implement
    //       the PRINT here. Even if OPESmetad::update() is skipped we should
    //       still call Print::update()
    writeTrajBuffer();
    if (m_pmf_grid_on) error_code |= collectSampleToPMFGrid();
    m_is_first_step = false;
    return COLVARS_OK;
  }
  error_code |= update_opes();
  if (error_code != COLVARS_OK) return error_code;
  writeTrajBuffer(); // Print::update()
  if (m_pmf_grid_on) error_code |= collectSampleToPMFGrid();
  return error_code;
}

int colvarbias_opes::collectSampleToPMFGrid() {
  if (m_reweight_grid) {
    // Get the bin index
    std::vector<int> bin(m_pmf_cvs.size(), 0);
    for (size_t i = 0; i < m_pmf_cvs.size(); ++i) {
      bin[i] = m_reweight_grid->current_bin_scalar(i);
    }
    const cvm::real reweighting_factor = cvm::exp(bias_energy / m_kbt);
    if (m_reweight_grid->index_ok(bin)) {
      m_reweight_grid->acc_value(bin, reweighting_factor);
    }
  }
  return COLVARS_OK;
}

template <typename OST> OST& colvarbias_opes::write_state_data_template_(OST &os) const {
  std::ios_base::fmtflags f;
  const bool formatted = !std::is_same<OST, cvm::memory_stream>::value;
  if (formatted) {
    f = os.flags();
    os.setf(std::ios::scientific, std::ios::floatfield);
  }
  write_state_data_key(os, "opes_metad_" + this->name);
  auto printFieldReal = [&](const std::string& s, cvm::real x){
    write_state_data_key(os, s, false);
    if (formatted)
      os << std::setprecision(cvm::en_prec) << std::setw(cvm::en_width);
    os << x;
    if (formatted)
      os << "\n";
  };
  auto printFieldULL = [&](const std::string& s, unsigned long long x){
    write_state_data_key(os, s, false);
    if (formatted)
      os << std::setprecision(cvm::en_prec) << std::setw(cvm::en_width);
    os << x;
    if (formatted)
      os << "\n";
  };
  auto printFieldString = [&](const std::string& s, const std::string& x){
    write_state_data_key(os, s, false);
    if (formatted)
      os << std::setprecision(cvm::en_prec) << std::setw(cvm::en_width);
    os << x;
    if (formatted)
      os << "\n";
  };
  std::ostringstream oss;
  if (m_inf_biasfactor) {
    oss << "inf";
  } else {
    oss << m_biasfactor;
  }
  printFieldString("biasfactor", oss.str());
  printFieldReal("epsilon", m_epsilon);
  printFieldReal("kernel_cutoff", cvm::sqrt(m_cutoff2));
  printFieldReal("compression_threshold", m_compression_threshold);
  printFieldReal("zed", m_saved_zed);
  printFieldReal("sum_weights", m_saved_sum_weights);
  printFieldReal("sum_weights2", m_saved_sum_weights2);
  printFieldULL("counter", m_counter);
  if (m_adaptive_sigma) {
    printFieldULL("adaptive_counter", m_adaptive_counter);
    for (size_t i = 0; i < num_variables(); ++i) {
      printFieldReal("sigma0_" + variables(i)->name, m_sigma0[i]);
      printFieldReal("av_cv_" + variables(i)->name, m_av_cv[i]);
      printFieldReal("av_M2_" + variables(i)->name, m_av_M2[i]);
    }
  }
  printFieldULL("num_hills", m_saved_kernels.size());
  write_state_data_key(os, "hills", false);
  if (formatted) os << "{\n";
  for (size_t k = 0; k < m_saved_kernels.size(); ++k) {
    if (formatted) os << "{ ";
    os << k;
    if (formatted) os << " ";
    for (size_t i = 0; i < num_variables(); ++i) {
      os << m_saved_kernels[k].m_center[i];
      if (formatted) os << " ";
    }
    for (size_t i = 0; i < num_variables(); ++i) {
      os << m_saved_kernels[k].m_sigma[i];
      if (formatted) os << " ";
    }
    os << m_saved_kernels[k].m_height;
    if (formatted) os << " }\n";
  }
  if (formatted) os << "}\n";
  if (formatted) os.setf(f);
  if (m_pmf_grid_on) {
    write_state_data_key(os, "probability_grid");
    m_reweight_grid->write_raw(os, 8);
  }
  return os;
}

std::ostream& colvarbias_opes::write_state_data(std::ostream &os) {
  try {
    auto& s = write_state_data_template_<std::ostream>(os);
    return s;
  } catch (const std::exception& e) {
    cvm::error(e.what());
  }
  return os;
}

cvm::memory_stream& colvarbias_opes::write_state_data(cvm::memory_stream& os) {
  try {
    auto& s = write_state_data_template_<cvm::memory_stream>(os);
    return s;
  } catch (const std::exception& e) {
    cvm::error(e.what());
  }
  return os;
}

template <typename IST> IST& colvarbias_opes::read_state_data_template_(IST &is) {
  bool const formatted = !std::is_same<IST, cvm::memory_stream>::value;
  std::string tmp_name;
  is >> tmp_name;
  if (tmp_name.rfind("opes_metad_", 0) != 0) {
    throw std::runtime_error("Unknown action name: " + tmp_name + "\n");
  }
  auto readFieldString = [&](const std::string& s, std::string& x){
    std::string field;
    is >> field;
    if (field.compare(s) == 0) {
      is >> x;
    } else {
      throw std::runtime_error("Expect field \"" + s + "\" , but got \"" + field + "\"\n");
    }
  };
  auto readFieldReal = [&](const std::string& s, cvm::real& x){
    std::string field;
    is >> field;
    if (field.compare(s) == 0) {
      is >> x;
    } else {
      throw std::runtime_error("Expect field \"" + s + "\" , but got \"" + field + "\"\n");
    }
  };
  auto readFieldULL = [&](const std::string& s, unsigned long long& x){
    std::string field;
    is >> field;
    if (field.compare(s) == 0) {
      is >> x;
    } else {
      throw std::runtime_error("Expect field \"" + s + "\" , but got \"" + field + "\"\n");
    }
  };
  std::string old_biasfactor_str;
  cvm::real old_biasfactor;
  readFieldString("biasfactor", old_biasfactor_str);
  if (old_biasfactor_str == "inf" || old_biasfactor_str == "-inf" || old_biasfactor_str == "+inf" ||
      old_biasfactor_str == "INF" || old_biasfactor_str == "-INF" || old_biasfactor_str == "+INF") {
    old_biasfactor = std::numeric_limits<cvm::real>::infinity();
    m_inf_biasfactor = true;
  } else {
    old_biasfactor = std::stod(old_biasfactor_str);
    m_inf_biasfactor = false;
  }
  if (std::abs(old_biasfactor - m_biasfactor) > 1e-6 * m_biasfactor) {
    cvm::log("WARNING: previous bias factor was " + cvm::to_str(old_biasfactor) +
             " while now it is " + cvm::to_str(m_biasfactor) +
             " (the new one is used).\n");
  }
  cvm::real old_epsilon;
  readFieldReal("epsilon", old_epsilon);
  if (std::abs(old_epsilon - m_epsilon) > 1e-6 * m_epsilon) {
    cvm::log("WARNING: previous epsilon was " + cvm::to_str(old_epsilon) +
             " while now it is " + cvm::to_str(m_epsilon) +
             " (the new one is used).\n");
  }
  cvm::real old_cutoff;
  readFieldReal("kernel_cutoff", old_cutoff);
  if (std::abs(old_cutoff - m_cutoff) > 1e-6 * m_cutoff) {
    cvm::log("WARNING: previous cutoff was " + cvm::to_str(old_cutoff) +
             " while now it is " + cvm::to_str(m_cutoff) +
             " (the new one is used).\n");
  }
  m_cutoff2 = m_cutoff * m_cutoff;
  cvm::real old_compression_threshold;
  readFieldReal("compression_threshold", old_compression_threshold);
  if (std::abs(old_compression_threshold - m_compression_threshold) > 1e-6 * m_compression_threshold) {
    cvm::log("WARNING: previous cutoff was " + cvm::to_str(old_compression_threshold) +
             " while now it is " + cvm::to_str(m_compression_threshold) +
             " (the new one is used).\n");
  }
  m_compression_threshold2 = m_compression_threshold * m_compression_threshold;
  readFieldReal("zed", m_zed);
  readFieldReal("sum_weights", m_sum_weights);
  readFieldReal("sum_weights2", m_sum_weights2);
  unsigned long long tmp_counter = 1;
  readFieldULL("counter", tmp_counter);
  m_counter = tmp_counter;
  if (m_adaptive_sigma) {
    readFieldULL("adaptive_counter", tmp_counter);
    m_adaptive_counter = tmp_counter;
    for (size_t i = 0; i < num_variables(); ++i) {
      readFieldReal("sigma0_" + variables(i)->name, m_sigma0[i]);
      readFieldReal("av_cv_" + variables(i)->name, m_av_cv[i]);
      readFieldReal("av_M2_" + variables(i)->name, m_av_M2[i]);
    }
  }
  unsigned long long kernel_size = 0;
  readFieldULL("num_hills", kernel_size);
  if (kernel_size > 0) m_kernels.resize(kernel_size);
  read_state_data_key(is, "hills");
  auto consume = [&](const std::string& expected_token){
    if (formatted) {
      std::string field;
      is >> field;
      if (field.compare(expected_token) != 0) {
        throw std::runtime_error("Expect " + expected_token + " but got " + field + "\n");
      }
    }
  };
  consume("{");
  for (size_t k = 0; k < m_kernels.size(); ++k) {
    consume("{");
    unsigned long long tmp_k = 0;
    is >> tmp_k;
    if (formatted && k != tmp_k) {
      throw std::runtime_error("Corrupt hill data\n");
    }
    kernel current_kernel;
    current_kernel.m_center.resize(num_variables());
    current_kernel.m_sigma.resize(num_variables());
    for (size_t i = 0; i < num_variables(); ++i) {
      is >> current_kernel.m_center[i];
    }
    for (size_t i = 0; i < num_variables(); ++i) {
      is >> current_kernel.m_sigma[i];
    }
    is >> current_kernel.m_height;
    m_kernels[k] = current_kernel;
    consume("}");
  }
  consume("}");
  if (m_pmf_grid_on) {
    read_state_data_key(is, "probability_grid");
    m_reweight_grid->read_raw(is);
  }
  m_kdenorm = m_explore ? m_counter : m_sum_weights;
  m_traj_line.rct = m_kbt * cvm::logn(m_sum_weights / m_counter);
  m_traj_line.zed = m_zed;
  m_traj_line.neff = (1 + m_sum_weights) * (1 + m_sum_weights) / (1 + m_sum_weights2);
  m_traj_line.nker = m_kernels.size();
  showInfo();
  return is;
}

std::istream& colvarbias_opes::read_state_data(std::istream &is) {
  try {
    auto& ret = read_state_data_template_<std::istream>(is);
    return ret;
  } catch (const std::exception& e) {
    cvm::error(e.what());
  }
  return is;
}

cvm::memory_stream& colvarbias_opes::read_state_data(cvm::memory_stream &is) {
  try {
    auto& ret = read_state_data_template_<cvm::memory_stream>(is);
    return ret;
  } catch (const std::exception& e) {
    cvm::error(e.what());
  }
  return is;
}

void colvarbias_opes::addKernel(const double height, const std::vector<cvm::real>& center, const std::vector<cvm::real>& sigma, const double logweight) {
  addKernel(height,center,sigma);
  const std::ios_base::fmtflags f = m_kernels_output.flags();
  m_kernels_output << std::right;
  // simulation time in ps
  m_kernels_output << std::setw(24) << (cvm::step_absolute() * cvm::dt()) * 1e-3;
  for (size_t i = 0; i < num_variables(); ++i) {
    m_kernels_output << " " << std::setw(24) << std::setprecision(16) <<  center[i];
  }
  for (size_t i = 0; i < num_variables(); ++i) {
    m_kernels_output << " " << std::setw(24) << std::setprecision(16) << sigma[i];
  }
  m_kernels_output << " " << std::setw(24) << std::setprecision(16) << height;
  m_kernels_output << " " << std::setw(24) << std::setprecision(16) << logweight;
  m_kernels_output << std::endl;
  m_kernels_output.flags(f);
}

void colvarbias_opes::addKernel(const double height, const std::vector<cvm::real>& center, const std::vector<cvm::real>& sigma) {
  bool no_match = true;
  if (m_compression_threshold2 != 0) {
    size_t taker_k = getMergeableKernel(center, m_kernels.size());
    if (taker_k < m_kernels.size()) {
      no_match = false;
      m_delta_kernels.emplace_back(-1 * m_kernels[taker_k].m_height, m_kernels[taker_k].m_center, m_kernels[taker_k].m_sigma);
      mergeKernels(m_kernels[taker_k], kernel(height, center, sigma));
      m_delta_kernels.push_back(m_kernels[taker_k]);
      if (m_recursive_merge) {
        size_t giver_k = taker_k;
        taker_k = getMergeableKernel(m_kernels[giver_k].m_center, giver_k);
        while (taker_k < m_kernels.size()) {
          m_delta_kernels.pop_back();
          m_delta_kernels.emplace_back(-1 * m_kernels[taker_k].m_height, m_kernels[taker_k].m_center, m_kernels[taker_k].m_sigma);
          if (taker_k > giver_k) std::swap(taker_k, giver_k);
          mergeKernels(m_kernels[taker_k], m_kernels[giver_k]);
          m_delta_kernels.push_back(m_kernels[taker_k]);
          m_kernels.erase(m_kernels.begin() + giver_k);
          if (m_nlist) {
            size_t giver_nk = 0;
            bool found_giver = false;
            for (size_t nk = 0; nk < m_nlist_index.size(); ++nk) {
              if (found_giver) m_nlist_index[nk]--;
              if (m_nlist_index[nk] == giver_k) {
                giver_nk = nk;
                found_giver = true;
              }
            }
            if (found_giver == false) {
              cvm::error("problem with merging and nlist\n");
            }
            m_nlist_index.erase(m_nlist_index.begin() + giver_nk);
          }
          giver_k = taker_k;
          taker_k = getMergeableKernel(m_kernels[giver_k].m_center, giver_k);
        }
      }
    }
  }
  if (no_match) {
    m_kernels.emplace_back(height, center, sigma);
    m_delta_kernels.emplace_back(height, center, sigma);
    if (m_nlist) m_nlist_index.push_back(m_kernels.size() - 1);
  }
}

void colvarbias_opes::mergeKernels(kernel& k1, const kernel& k2) const {
  const double h = k1.m_height + k2.m_height;
  for (size_t i = 0; i < k1.m_center.size(); ++i) {
    const bool isPeriodic_i = variables(i)->is_enabled(f_cv_periodic);
    if (isPeriodic_i) {
      k1.m_center[i] = k2.m_center[i] + 0.5 * variables(i)->dist2_lgrad(k1.m_center[i], k2.m_center[i]).real_value;
    }
    const cvm::real c_i = (k1.m_height * k1.m_center[i] +
                           k2.m_height * k2.m_center[i]) / h;
    const cvm::real ss_k1_part = k1.m_height * (k1.m_sigma[i] * k1.m_sigma[i] + k1.m_center[i] * k1.m_center[i]);
    const cvm::real ss_k2_part = k2.m_height * (k2.m_sigma[i] * k2.m_sigma[i] + k2.m_center[i] * k2.m_center[i]);
    const cvm::real ss_i = (ss_k1_part + ss_k2_part) / h - c_i * c_i;
    if (isPeriodic_i) {
      colvarvalue tmp(c_i);
      variables(i)->wrap(tmp);
      k1.m_center[i] = tmp.real_value;
    } else {
      k1.m_center[i] = c_i;
    }
    k1.m_sigma[i] = cvm::sqrt(ss_i);
  }
  k1.m_height = h;
}

size_t colvarbias_opes::getMergeableKernel(const std::vector<cvm::real>& giver_center, const size_t giver_k) const {
  size_t min_k = m_kernels.size();
  cvm::real min_norm2 = m_compression_threshold2;
  const int num_parallel = 1;
  if (!m_nlist) {
    if (m_num_threads == 1) {
      for (size_t k = 0; k < m_kernels.size(); ++k) {
        if (k == giver_k) continue;
        double norm2 = 0;
        for (size_t i = 0; i < num_variables(); ++i) {
          norm2 += variables(i)->dist2(giver_center[i], m_kernels[k].m_center[i]) / (m_kernels[k].m_sigma[i] * m_kernels[k].m_sigma[i]);
          if (norm2 >= min_norm2) break;
        }
        if (norm2 < min_norm2) {
          min_norm2 = norm2;
          min_k = k;
        }
      }
    } else {
#if defined(_OPENMP)
      #pragma omp parallel num_threads(m_num_threads)
      {
        int min_k_omp = min_k;
        cvm::real min_norm2_omp = m_compression_threshold2;
        #pragma omp for nowait
        for (int k = 0; k < static_cast<int>(m_kernels.size()); ++k) {
          if (k == static_cast<int>(giver_k)) continue;
          double norm2 = 0;
          for (int i = 0; i < static_cast<int>(num_variables()); ++i) {
            norm2 += variables(i)->dist2( giver_center[i], m_kernels[k].m_center[i]) / (m_kernels[k].m_sigma[i] * m_kernels[k].m_sigma[i]);
            if (norm2 >= min_norm2_omp) break;
          }
          if (norm2 < min_norm2_omp) {
            min_norm2_omp = norm2;
            min_k_omp = k;
          }
        }
        #pragma omp critical
        {
          if (min_norm2_omp < min_norm2) {
            min_norm2 = min_norm2_omp;
            min_k = min_k_omp;
          }
        }
      }
#elif CMK_SMP && USE_CKLOOP
      // NOTE: No existing reduction type for finding the minimum, so I have
      //       to use such a workaround.
      std::vector<size_t> min_k_smp(m_num_threads, min_k);
      std::vector<cvm::real> min_norm2_smp(m_num_threads, m_compression_threshold2);
      auto worker = [&](int start, int end, void* unused) {
        const int tid = cvm::proxy->smp_thread_id();
        for (int k = start; k <= end; ++k) {
          if (k == giver_k) continue;
          double norm2 = 0;
          for (size_t i = 0; i < num_variables(); ++i) {
            norm2 += variables(i)->dist2( giver_center[i], m_kernels[k].m_center[i]) / (m_kernels[k].m_sigma[i] * m_kernels[k].m_sigma[i]);
            if (norm2 >= min_norm2_smp[tid]) break;
          }
          if (norm2 < min_norm2_smp[tid]) {
            min_norm2_smp[tid] = norm2;
            min_k_smp[tid] = k;
          }
        }
      };
      const size_t lowerRange = 0;
      const size_t upperRange = m_kernels.size() - 1;
      CkLoop_Parallelize(
        m_num_threads, lowerRange, upperRange,
        worker, NULL, CKLOOP_NONE, NULL);
      const auto it_min = std::min_element(min_norm2_smp.begin(), min_norm2_smp.end());
      min_norm2 = *it_min;
      min_k = min_k_smp[std::distance(min_norm2_smp.begin(), it_min)];
#else
      cvm::error("OPES cannot run because this binary is not linked with a supported threading library.\n");
#endif
    }
  } else {
    if (m_num_threads == 1) {
      // size_t min_k_omp = min_k;
      // cvm::real min_norm2_omp = m_compression_threshold2;
      for (size_t nk = 0; nk < m_nlist_index.size(); ++nk) {
        const size_t k = m_nlist_index[nk];
        if (k == giver_k) continue;
        double norm2 = 0;
        for (size_t i = 0; i < num_variables(); ++i) {
          norm2 += variables(i)->dist2(giver_center[i], m_kernels[k].m_center[i]) / (m_kernels[k].m_sigma[i] * m_kernels[k].m_sigma[i]);
          if (norm2 >= min_norm2) break;
        }
        if (norm2 < min_norm2) {
          min_norm2 = norm2;
          min_k = k;
        }
      }
    } else {
#if defined(_OPENMP)
      #pragma omp parallel num_threads(m_num_threads)
      {
        size_t min_k_omp = min_k;
        cvm::real min_norm2_omp = m_compression_threshold2;
        #pragma omp for nowait
        for (int nk = 0; nk < static_cast<int>(m_nlist_index.size()); ++nk) {
          const size_t k = m_nlist_index[nk];
          if (k == giver_k) continue;
          double norm2 = 0;
          for (int i = 0; i < static_cast<int>(num_variables()); ++i) {
            norm2 += variables(i)->dist2(giver_center[i], m_kernels[k].m_center[i]) / (m_kernels[k].m_sigma[i] * m_kernels[k].m_sigma[i]);
            if (norm2 >= min_norm2_omp) break;
          }
          if (norm2 < min_norm2_omp) {
            min_norm2_omp = norm2;
            min_k_omp = k;
          }
        }
        #pragma omp critical
        {
          if (min_norm2_omp < min_norm2) {
            min_norm2 = min_norm2_omp;
            min_k = min_k_omp;
          }
        }
      }
#elif CMK_SMP && USE_CKLOOP
      // NOTE: No existing reduction type for finding the minimum, so I have
      //       to use such a workaround.
      std::vector<size_t> min_k_smp(m_num_threads, min_k);
      std::vector<cvm::real> min_norm2_smp(m_num_threads, m_compression_threshold2);
      auto worker = [&](int start, int end, void* unused) {
        const int tid = cvm::proxy->smp_thread_id();
        for (int nk = start; nk <= end; ++nk) {
          const size_t k = m_nlist_index[nk];
          if (k == giver_k) continue;
          double norm2 = 0;
          for (size_t i = 0; i < num_variables(); ++i) {
            norm2 += variables(i)->dist2( giver_center[i], m_kernels[k].m_center[i]) / (m_kernels[k].m_sigma[i] * m_kernels[k].m_sigma[i]);
            if (norm2 >= min_norm2_smp[tid]) break;
          }
          if (norm2 < min_norm2_smp[tid]) {
            min_norm2_smp[tid] = norm2;
            min_k_smp[tid] = k;
          }
        }
      };
      const size_t lowerRange = 0;
      const size_t upperRange = m_nlist_index.size() - 1;
      CkLoop_Parallelize(
        m_num_threads, lowerRange, upperRange,
        worker, NULL, CKLOOP_NONE, NULL);
      const auto it_min = std::min_element(min_norm2_smp.begin(), min_norm2_smp.end());
      min_norm2 = *it_min;
      min_k = min_k_smp[std::distance(min_norm2_smp.begin(), it_min)];
#else
      cvm::error("OPES cannot run because this binary is not linked with a supported threading library.\n");
#endif
    }
  }
  if (num_parallel > 1) {
    cvm::error("The Colvars OPES implementation does not support running OPES in parallel across nodes.\n");
  }
  return min_k;
}

std::string const colvarbias_opes::traj_file_name(const std::string& suffix) const {
  return std::string(cvm::output_prefix()+
                     ".colvars."+this->name+
                     ( (comm != single_replica) ?
                       ("."+replica_id) :
                       ("") )+
                     suffix);
}

int colvarbias_opes::write_output_files() {
  int error_code = COLVARS_OK;
  thread_local static bool firsttime = true;
  // Write the kernels
  const std::string kernels_filename = traj_file_name(".kernels.dat");
  std::ostream& os_kernels = cvm::proxy->output_stream(kernels_filename, "kernels file");
  const std::ios_base::fmtflags format_kernels = os_kernels.flags();
  if (firsttime) {
    os_kernels << "#! FIELDS time ";
    for (size_t i = 0; i < num_variables(); ++i) {
      os_kernels << variables(i)->name + " ";
    }
    for (size_t i = 0; i < num_variables(); ++i) {
      os_kernels << "sigma_" + variables(i)->name + " ";
    }
    os_kernels << "height logweight\n";
    // Make sure the action name compatible with the script in https://github.com/invemichele/opes/blob/master/postprocessing/State_from_Kernels.py
    if (m_explore) os_kernels << "#! SET action OPES_METAD_EXPLORE_kernels\n";
    else os_kernels << "#! SET action OPES_METAD_kernels\n";
    if (m_inf_biasfactor) {
      os_kernels << "#! SET biasfactor " << "inf" << "\n";
    } else {
      os_kernels << "#! SET biasfactor " << m_biasfactor << "\n";
    }
    os_kernels << "#! SET epsilon " << m_epsilon << "\n";
    os_kernels << "#! SET kernel_cutoff " << m_cutoff << "\n";
    os_kernels << "#! SET compression_threshold " << m_compression_threshold << "\n";
    for (size_t i = 0; i < num_variables(); ++i) {
      if (variables(i)->is_enabled(f_cv_periodic)) {
        if (variables(i)->is_enabled(f_cv_lower_boundary)) {
          os_kernels << "#! SET min_" + variables(i)->name + " " << variables(i)->lower_boundary.real_value << "\n";
        }
        if (variables(i)->is_enabled(f_cv_upper_boundary)) {
          os_kernels << "#! SET max_" + variables(i)->name + " " << variables(i)->upper_boundary.real_value << "\n";
        }
      }
    }
  }
  os_kernels << m_kernels_output.str();
  os_kernels.setf(format_kernels);
  error_code |= cvm::proxy->flush_output_stream(kernels_filename);
  m_kernels_output.str("");
  m_kernels_output.clear();

  // Write the trajectory
  const std::string traj_filename = traj_file_name(".misc.traj");
  std::ostream& os_traj = cvm::proxy->output_stream(traj_filename, "trajectory of various OPES properties");
  const std::ios_base::fmtflags format_traj = os_traj.flags();
  if (firsttime) {
    os_traj << "#! FIELDS time ";
    for (size_t i = 0; i < num_variables(); ++i) {
      os_traj << variables(i)->name + " ";
    }
    os_traj << this->name + ".bias ";
    os_traj << this->name + ".rct ";
    if (!m_no_zed) os_traj << this->name + ".zed ";
    os_traj << this->name + ".neff ";
    if (m_calc_work) if (!m_no_zed) os_traj << this->name + ".work ";
    os_traj << this->name + ".nker ";
    if (m_nlist) os_traj << this->name + ".nlker ";
    if (m_nlist) os_traj << this->name + ".nlsteps ";
    os_traj << "\n";
    for (size_t i = 0; i < num_variables(); ++i) {
      if (variables(i)->is_enabled(f_cv_lower_boundary)) {
        os_traj << "#! SET min_" + variables(i)->name + " " << variables(i)->lower_boundary.real_value << "\n";
      }
      if (variables(i)->is_enabled(f_cv_upper_boundary)) {
        os_traj << "#! SET max_" + variables(i)->name + " " << variables(i)->upper_boundary.real_value << "\n";
      }
    }
  }
  os_traj << m_traj_oss.str();
  os_traj.setf(format_traj);
  error_code |= cvm::proxy->flush_output_stream(traj_filename);
  m_traj_oss.str("");
  m_traj_oss.clear();
  if (firsttime) firsttime = false;
  if (m_pmf_grid_on) {
    error_code |= computePMF();
    const std::string pmf_filename = traj_file_name(".pmf");
    error_code |= writePMF(m_pmf_grid, pmf_filename, false);
    if (comm == multiple_replicas && m_pmf_shared) {
      if (cvm::proxy->replica_index() == 0) {
        const std::string global_pmf_filename = traj_file_name(".global.pmf");
        error_code |= writePMF(m_global_pmf_grid, global_pmf_filename, false);
      }
    }
    if (m_pmf_hist_freq > 0 && cvm::step_absolute() % m_pmf_hist_freq == 0) {
      const std::string pmf_hist_filename = traj_file_name(".hist.pmf");
      error_code |= writePMF(m_pmf_grid, pmf_hist_filename, true);
      if (comm == multiple_replicas && m_pmf_shared) {
        if (cvm::proxy->replica_index() == 0) {
          const std::string global_hist_pmf_filename = traj_file_name(".global.hist.pmf");
          error_code |= writePMF(m_global_pmf_grid, global_hist_pmf_filename, true);
        }
      }
    }
  }
  // To prevent the case that one replica exits earlier and then destroys all streams
  if (comm == multiple_replicas) cvm::proxy->replica_comm_barrier();
  return error_code;
}

void hist_to_pmf(const cvm::real kbt, const colvar_grid_scalar *hist, std::unique_ptr<colvar_grid_scalar>& pmf) {
  // Get the sum of probabilities of all grids
  cvm::real norm_factor = 0;
  cvm::real max_prob = 0;
  auto& prob_data = hist->data;
  for (auto it = prob_data.begin(); it != prob_data.end(); ++it) {
    norm_factor += (*it);
    if ((*it) > max_prob) max_prob = (*it);
  }
  if (norm_factor > 0) {
    const cvm::real min_pmf = (max_prob > 0) ? -1.0 * kbt * cvm::logn(max_prob / norm_factor) : 0;
    auto& pmf_data = pmf->data;
    for (size_t i = 0; i < pmf_data.size(); ++i) {
      if (prob_data[i] > 0) {
        pmf_data[i] = -1.0 * kbt * cvm::logn(prob_data[i] / norm_factor) - min_pmf;
      }
    }
    auto max_pmf = *std::max_element(pmf_data.begin(), pmf_data.end());
    for (size_t i = 0; i < pmf_data.size(); ++i) {
      if (!(prob_data[i] > 0)) {
        pmf_data[i] = max_pmf;
      }
    }
  }
}

int colvarbias_opes::computePMF() {
  // Multiple replica: collect all samples from other replicas
  if (comm == multiple_replicas && m_pmf_shared) {
    const size_t samples_n = m_reweight_grid->raw_data_num();
    const int msg_size = samples_n * sizeof(cvm::real);
    std::vector<cvm::real> buffer;
    if (cvm::main()->proxy->replica_index() == 0) {
      buffer.resize(samples_n * (cvm::proxy->num_replicas() - 1));
      for (int p = 1; p < cvm::proxy->num_replicas(); p++) {
        const size_t start_pos = (p - 1) * samples_n;
        if (cvm::proxy->replica_comm_recv((char*)&(buffer[start_pos]), msg_size, p) != msg_size) {
          return cvm::error("Error getting shared OPES reweighting histogram from replica " + cvm::to_str(p));
        }
      }
    } else {
      if (cvm::proxy->replica_comm_send((char*)(&(m_reweight_grid->data[0])), msg_size, 0) != msg_size) {
        return cvm::error("Error sending shared OPES reweighting histogram from replica " + cvm::to_str(cvm::main()->proxy->replica_index()));
      }
    }
    cvm::proxy->replica_comm_barrier();
    // Broadcast m_reweight_grid to all replicas
    auto& global_data = m_global_reweight_grid->data;
    if (cvm::main()->proxy->replica_index() == 0) {
      global_data = m_reweight_grid->data;
      // Sum the samples on PE 0
      for (int p = 1; p < cvm::proxy->num_replicas(); p++) {
        const size_t start_pos = (p - 1) * samples_n;
        for (size_t i = 0 ; i < samples_n; ++i) {
          global_data[i] += buffer[start_pos+i];
        }
      }
    }
  }
  // Get the sum of probabilities of all grids
  hist_to_pmf(m_kbt, m_reweight_grid.get(), m_pmf_grid);
  if (comm == multiple_replicas && m_pmf_shared) {
    if (cvm::main()->proxy->replica_index() == 0) {
      hist_to_pmf(m_kbt, m_global_reweight_grid.get(), m_global_pmf_grid);
    }
  }
  if (comm == multiple_replicas) {
    cvm::proxy->replica_comm_barrier();
  }
  return COLVARS_OK;
}

int colvarbias_opes::writePMF(const std::unique_ptr<colvar_grid_scalar>& pmf_grid, const std::string &filename, bool keep_open) {
  std::ostream& os = cvm::proxy->output_stream(filename, "output stream of " + filename);
  if (!os) {
    return COLVARS_FILE_ERROR;
  }
  pmf_grid->write_multicol(os);
  if (!keep_open) {
    cvm::proxy->close_output_stream(filename);
  } else {
    cvm::proxy->flush_output_stream(filename);
  }
  return COLVARS_OK;
}

void colvarbias_opes::writeTrajBuffer() {
  if (m_traj_output_frequency > 0 && cvm::step_absolute() % m_traj_output_frequency == 0) {
    m_traj_oss << std::right;
    m_traj_oss << std::scientific << " " << std::setw(cvm::cv_width) << std::setprecision(cvm::cv_prec) << (cvm::step_absolute() * cvm::dt()) * 1e-3;
    for (size_t i = 0; i < num_variables(); ++i) {
      m_traj_oss << std::scientific << " " << std::setw(cvm::cv_width) << std::setprecision(cvm::cv_prec) << variables(i)->value().real_value;
    }
    m_traj_oss << std::scientific << " " << std::setw(cvm::cv_width) << std::setprecision(cvm::cv_prec) << bias_energy;
    m_traj_oss << std::scientific << " " << std::setw(cvm::cv_width) << std::setprecision(cvm::cv_prec) << m_traj_line.rct;
    if (!m_no_zed) m_traj_oss << std::scientific << " " << std::setw(cvm::cv_width) << std::setprecision(cvm::cv_prec) << m_traj_line.zed;
    m_traj_oss << std::scientific << " " << std::setw(cvm::cv_width) << std::setprecision(cvm::cv_prec) << m_traj_line.neff;
    if (m_calc_work) m_traj_oss << std::scientific << " " << std::setw(cvm::cv_width) << std::setprecision(cvm::cv_prec) << m_traj_line.work;
    m_traj_oss << " " << m_traj_line.nker;
    if (m_nlist) m_traj_oss << " " << m_traj_line.nlker;
    if (m_nlist) m_traj_oss << " " << m_traj_line.nlsteps;
    m_traj_oss << "\n";
  }
}

void colvarbias_opes::updateNlist(const std::vector<cvm::real>& center) {
  if (m_kernels.empty()) return;
  m_nlist_center = center;
  m_nlist_index.clear();
  if (m_num_threads == 1 || m_kernels.size() < 2 * m_num_threads) {
    for (size_t k = 0; k < m_kernels.size(); ++k) {
      cvm::real norm2_k = 0;
      for (size_t i = 0; i < num_variables(); ++i) {
        norm2_k += variables(i)->dist2(m_nlist_center[i], m_kernels[k].m_center[i]) / (m_kernels[k].m_sigma[i] * m_kernels[k].m_sigma[i]);
      }
      if (norm2_k <= m_nlist_param[0] * m_cutoff2) {
        m_nlist_index.push_back(k);
      }
    }
  } else {
#if defined (_OPENMP)
    #pragma omp parallel num_threads(m_num_threads)
    {
      std::vector<size_t> private_nlist_index;
      #pragma omp for nowait
      for (int k = 0; k < static_cast<int>(m_kernels.size()); ++k) {
        cvm::real norm2_k = 0;
        for (int i = 0; i < static_cast<int>(num_variables()); ++i) {
          norm2_k += variables(i)->dist2(m_nlist_center[i], m_kernels[k].m_center[i]) / (m_kernels[k].m_sigma[i] * m_kernels[k].m_sigma[i]);
        }
        if (norm2_k <= m_nlist_param[0] * m_cutoff2) {
          private_nlist_index.push_back(static_cast<size_t>(k));
        }
      }
      #pragma omp critical
      m_nlist_index.insert(m_nlist_index.end(), private_nlist_index.begin(), private_nlist_index.end());
    }
#elif CMK_SMP && USE_CKLOOP
    std::vector<std::vector<size_t>> private_nlist_index(m_num_threads);
    auto worker = [&](int start, int end, void* unused){
      const int tid = cvm::proxy->smp_thread_id();
      for (int k = start; k <= end; ++k) {
        cvm::real norm2_k = 0;
        for (size_t i = 0; i < num_variables(); ++i) {
          norm2_k += variables(i)->dist2(m_nlist_center[i], m_kernels[k].m_center[i]) / (m_kernels[k].m_sigma[i] * m_kernels[k].m_sigma[i]);
        }
        if (norm2_k <= m_nlist_param[0] * m_cutoff2) {
          private_nlist_index[tid].push_back(k);
        }
      }
    };
    const size_t lowerRange = 0;
    const size_t upperRange = m_kernels.size() - 1;
    CkLoop_Parallelize(
      m_num_threads, lowerRange, upperRange,
      worker, NULL, CKLOOP_NONE, NULL);
    for (size_t i = 0; i < m_num_threads; ++i) {
      m_nlist_index.insert(m_nlist_index.end(), private_nlist_index[i].begin(), private_nlist_index[i].end());
    }
#else
    cvm::error("OPES cannot run because this binary is not linked with a supported threading library.\n");
#endif
    if (m_recursive_merge) {
      std::sort(m_nlist_index.begin(), m_nlist_index.end());
    }
  }
  std::vector<cvm::real> dev2(num_variables(), 0);
  for (size_t k = 0; k < m_nlist_index.size(); ++k) {
    for (size_t i = 0; i < num_variables(); ++i) {
      dev2[i] += variables(i)->dist2(m_nlist_center[i], m_kernels[m_nlist_index[k]].m_center[i]);
    }
  }
  for (size_t i = 0; i < num_variables(); ++i) {
    if (m_nlist_index.empty()) {
      m_nlist_dev2[i] = m_kernels.back().m_sigma[i] * m_kernels.back().m_sigma[i];
    } else {
      m_nlist_dev2[i] = dev2[i] / m_nlist_index.size();
    }
  }
  m_traj_line.nlker = m_nlist_index.size();
  m_traj_line.nlsteps = m_nlist_steps;
  m_nlist_steps = 0;
  m_nlist_update = false;
}
