// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvar_gpu_support.h"
#include "colvarmodule.h"
#include "colvartypes.h"
#include "colvaratoms.h"
#include "colvar_rotation_derivative.h"
#include "colvarproxy.h"
#include "cuda/colvartypes_kernel.h"

#ifdef COLVARS_LAMMPS
// Use open-source Jacobi implementation
#include "math_eigen_impl.h"
#else
// Fall back to NR routine
#include "nr_jacobi.h"
#endif

bool      colvarmodule::rotation::monitor_crossings = false;
cvm::real colvarmodule::rotation::crossing_threshold = 1.0E-02;


std::string cvm::rvector::to_simple_string() const
{
  std::ostringstream os;
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(cvm::main()->cv_prec);
  os << x << " " << y << " " << z;
  return os.str();
}


int cvm::rvector::from_simple_string(std::string const &s)
{
  std::stringstream stream(s);
  if ( !(stream >> x) ||
       !(stream >> y) ||
       !(stream >> z) ) {
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
}


std::ostream & operator << (std::ostream &os, colvarmodule::rvector const &v)
{
  std::streamsize const w = os.width();
  std::streamsize const p = os.precision();

  os.width(2);
  os << "( ";
  os.width(w); os.precision(p);
  os << v.x << " , ";
  os.width(w); os.precision(p);
  os << v.y << " , ";
  os.width(w); os.precision(p);
  os << v.z << " )";
  return os;
}


std::istream & operator >> (std::istream &is, colvarmodule::rvector &v)
{
  std::streampos const start_pos = is.tellg();
  char sep;
  if ( !(is >> sep) || !(sep == '(') ||
       !(is >> v.x) || !(is >> sep)  || !(sep == ',') ||
       !(is >> v.y) || !(is >> sep)  || !(sep == ',') ||
       !(is >> v.z) || !(is >> sep)  || !(sep == ')') ) {
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
    return is;
  }
  return is;
}

std::string cvm::quaternion::to_simple_string() const
{
  std::ostringstream os;
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(cvm::main()->cv_prec);
  os << q0 << " " << q1 << " " << q2 << " " << q3;
  return os.str();
}

int cvm::quaternion::from_simple_string(std::string const &s)
{
  std::stringstream stream(s);
  if ( !(stream >> q0) ||
       !(stream >> q1) ||
       !(stream >> q2) ||
       !(stream >> q3) ) {
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
}

std::ostream & operator << (std::ostream &os, colvarmodule::quaternion const &q)
{
  std::streamsize const w = os.width();
  std::streamsize const p = os.precision();

  os.width(2);
  os << "( ";
  os.width(w); os.precision(p);
  os << q.q0 << " , ";
  os.width(w); os.precision(p);
  os << q.q1 << " , ";
  os.width(w); os.precision(p);
  os << q.q2 << " , ";
  os.width(w); os.precision(p);
  os << q.q3 << " )";
  return os;
}


std::istream & operator >> (std::istream &is, colvarmodule::quaternion &q)
{
  std::streampos const start_pos = is.tellg();
  char sep;
  if ( !(is >> sep)  || !(sep == '(') ||
       !(is >> q.q0) || !(is >> sep)  || !(sep == ',') ||
       !(is >> q.q1) || !(is >> sep)  || !(sep == ',') ||
       !(is >> q.q2) || !(is >> sep)  || !(sep == ',') ||
       !(is >> q.q3) || !(is >> sep)  || !(sep == ')') ) {
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
  }
  return is;
}


#ifdef COLVARS_LAMMPS
namespace {
  inline void *new_Jacobi_solver(int size) {
    return reinterpret_cast<void *>(new MathEigen::Jacobi<cvm::real,
                                    cvm::vector1d<cvm::real> &,
                                    cvm::matrix2d<cvm::real> &>(4));
  }
}
#endif


int colvarmodule::rotation::init()
{
  b_debug_gradients = false;
  // lambda = 0.0;
  cvm::main()->cite_feature("Optimal rotation via flexible fitting");
  return COLVARS_OK;
}


colvarmodule::rotation::rotation()
{
  init();
#ifdef COLVARS_LAMMPS
  jacobi = new_Jacobi_solver(4);
#else
  jacobi = NULL;
#endif
}


colvarmodule::rotation::rotation(cvm::quaternion const &qi)
  : q(qi)
{
  init();
#ifdef COLVARS_LAMMPS
  jacobi = new_Jacobi_solver(4);
#else
  jacobi = NULL;
#endif
}


colvarmodule::rotation::rotation(cvm::real angle, cvm::rvector const &axis)
{
  init();
  cvm::rvector const axis_n = axis.unit();
  cvm::real const sina = cvm::sin(angle/2.0);
  q = cvm::quaternion(cvm::cos(angle/2.0),
                      sina * axis_n.x, sina * axis_n.y, sina * axis_n.z);
#ifdef COLVARS_LAMMPS
  jacobi = new_Jacobi_solver(4);
#else
  jacobi = NULL;
#endif
}


colvarmodule::rotation::~rotation()
{
#ifdef COLVARS_LAMMPS
  delete reinterpret_cast<
    MathEigen::Jacobi<cvm::real,
                      cvm::vector1d<cvm::real> &,
                      cvm::matrix2d<cvm::real> &> *>(jacobi);
#endif
}


void colvarmodule::rotation::build_correlation_matrix(
                                        std::vector<cvm::atom_pos> const &pos1,
                                        std::vector<cvm::atom_pos> const &pos2)
{
  // build the correlation matrix
  size_t i;
  for (i = 0; i < pos1.size(); i++) {
    C.xx += pos1[i].x * pos2[i].x;
    C.xy += pos1[i].x * pos2[i].y;
    C.xz += pos1[i].x * pos2[i].z;
    C.yx += pos1[i].y * pos2[i].x;
    C.yy += pos1[i].y * pos2[i].y;
    C.yz += pos1[i].y * pos2[i].z;
    C.zx += pos1[i].z * pos2[i].x;
    C.zy += pos1[i].z * pos2[i].y;
    C.zz += pos1[i].z * pos2[i].z;
  }
}


void colvarmodule::rotation::compute_overlap_matrix()
{
  // build the "overlap" matrix, whose eigenvectors are stationary
  // points of the RMSD in the space of rotations
  S[0][0] =    C.xx + C.yy + C.zz;
  S[1][0] =    C.yz - C.zy;
  S[0][1] = S[1][0];
  S[2][0] =  - C.xz + C.zx ;
  S[0][2] = S[2][0];
  S[3][0] =    C.xy - C.yx;
  S[0][3] = S[3][0];
  S[1][1] =    C.xx - C.yy - C.zz;
  S[2][1] =    C.xy + C.yx;
  S[1][2] = S[2][1];
  S[3][1] =    C.xz + C.zx;
  S[1][3] = S[3][1];
  S[2][2] = - C.xx + C.yy - C.zz;
  S[3][2] =   C.yz + C.zy;
  S[2][3] = S[3][2];
  S[3][3] = - C.xx - C.yy + C.zz;
}


#ifndef COLVARS_LAMMPS
namespace NR {

int diagonalize_matrix(cvm::real m[4][4],
                        cvm::real eigval[4],
                        cvm::real eigvec[4][4])
{
  std::memset(eigval, 0, sizeof(cvm::real) * 4);
  std::memset(eigvec, 0, sizeof(cvm::real) * 4 * 4);

  // diagonalize
  int jac_nrot = 0;
  if (NR_Jacobi::jacobi(m, eigval, eigvec, &jac_nrot) !=
      COLVARS_OK) {
    return COLVARS_ERROR;
  }
  NR_Jacobi::eigsrt(eigval, eigvec);
  // jacobi saves eigenvectors by columns
  NR_Jacobi::transpose(eigvec);

  // normalize eigenvectors
  for (size_t ie = 0; ie < 4; ie++) {
    cvm::real norm2 = 0.0;
    size_t i;
    for (i = 0; i < 4; i++) {
      norm2 += eigvec[ie][i] * eigvec[ie][i];
    }
    cvm::real const norm = cvm::sqrt(norm2);
    for (i = 0; i < 4; i++) {
      eigvec[ie][i] /= norm;
    }
  }
  return COLVARS_OK;
}

}
#endif


void colvarmodule::rotation::debug_gradients(
  cvm::rotation &rot,
  const cvm::ag_vector_real_t &pos1,
  const cvm::ag_vector_real_t &pos2,
  const size_t num_atoms_pos1,
  const size_t num_atoms_pos2) {
  // eigenvalues and eigenvectors
  cvm::real const L0 = rot.S_eigval[0];
  cvm::real const L1 = rot.S_eigval[1];
  cvm::real const L2 = rot.S_eigval[2];
  cvm::real const L3 = rot.S_eigval[3];
  cvm::quaternion const Q0(rot.S_eigvec[0]);
  cvm::quaternion const Q1(rot.S_eigvec[1]);
  cvm::quaternion const Q2(rot.S_eigvec[2]);
  cvm::quaternion const Q3(rot.S_eigvec[3]);

  cvm::log_static("L0 = "+cvm::to_str(L0, cvm::main()->cv_width, cvm::main()->cv_prec)+
            ", Q0 = "+cvm::to_str(Q0, cvm::main()->cv_width, cvm::main()->cv_prec)+
            ", Q0*Q0 = "+cvm::to_str(Q0.inner(Q0), cvm::main()->cv_width, cvm::main()->cv_prec)+
            "\n");
  cvm::log_static("L1 = "+cvm::to_str(L1, cvm::main()->cv_width, cvm::main()->cv_prec)+
            ", Q1 = "+cvm::to_str(Q1, cvm::main()->cv_width, cvm::main()->cv_prec)+
            ", Q0*Q1 = "+cvm::to_str(Q0.inner(Q1), cvm::main()->cv_width, cvm::main()->cv_prec)+
            "\n");
  cvm::log_static("L2 = "+cvm::to_str(L2, cvm::main()->cv_width, cvm::main()->cv_prec)+
            ", Q2 = "+cvm::to_str(Q2, cvm::main()->cv_width, cvm::main()->cv_prec)+
            ", Q0*Q2 = "+cvm::to_str(Q0.inner(Q2), cvm::main()->cv_width, cvm::main()->cv_prec)+
            "\n");
  cvm::log_static("L3 = "+cvm::to_str(L3, cvm::main()->cv_width, cvm::main()->cv_prec)+
            ", Q3 = "+cvm::to_str(Q3, cvm::main()->cv_width, cvm::main()->cv_prec)+
            ", Q0*Q3 = "+cvm::to_str(Q0.inner(Q3), cvm::main()->cv_width, cvm::main()->cv_prec)+
            "\n");
  rotation_derivative deriv(rot, pos1, pos2, num_atoms_pos1, num_atoms_pos2);
  cvm::rvector dl0_2;
  std::array<cvm::rvector, 4> dq0_2;
  std::array<std::array<cvm::rvector, 4>, 4> ds_2;
#ifdef COLVARS_LAMMPS
    MathEigen::Jacobi<cvm::real,
                      cvm::real[4],
                      cvm::real[4][4]> *ecalc =
        reinterpret_cast<MathEigen::Jacobi<cvm::real,
                                           cvm::real[4],
                                           cvm::real[4][4]> *>(rot.jacobi);
#endif
  deriv.prepare_derivative(rotation_derivative_dldq::use_dl | rotation_derivative_dldq::use_dq);
  cvm::real S_new[4][4];
  cvm::real S_new_eigval[4];
  cvm::real S_new_eigvec[4][4];
  // const size_t num_atoms = soa ? deriv.m_num_atoms_pos2 : pos2.size();
  for (size_t ia = 0; ia < num_atoms_pos2; ++ia) {
    deriv.calc_derivative_wrt_group2<true, true, true>(ia, &dl0_2, &dq0_2, &ds_2);
    // make an infitesimal move along each cartesian coordinate of
    // this atom, and solve again the eigenvector problem
    for (size_t comp = 0; comp < 3; comp++) {
      std::memcpy(S_new, rot.S_backup, sizeof(cvm::real) * 4 * 4);
      std::memset(S_new_eigval, 0, sizeof(cvm::real) * 4);
      std::memset(S_new_eigvec, 0, sizeof(cvm::real) * 4 * 4);
      for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
          S_new[i][j] +=
            colvarmodule::debug_gradients_step_size * ds_2[i][j][comp];
        }
      }
#ifdef COLVARS_LAMMPS
      ecalc->Diagonalize(S_new, S_new_eigval, S_new_eigvec);
#else
      NR::diagonalize_matrix(S_new, S_new_eigval, S_new_eigvec);
#endif
      cvm::real const &L0_new = S_new_eigval[0];
      cvm::quaternion const Q0_new(S_new_eigvec[0]);

      cvm::real const DL0 = (dl0_2[comp]) * colvarmodule::debug_gradients_step_size;
      cvm::quaternion const DQ0(dq0_2[0][comp] * colvarmodule::debug_gradients_step_size,
                                dq0_2[1][comp] * colvarmodule::debug_gradients_step_size,
                                dq0_2[2][comp] * colvarmodule::debug_gradients_step_size,
                                dq0_2[3][comp] * colvarmodule::debug_gradients_step_size);

      cvm::log_static(  "|(l_0+dl_0) - l_0^new|/l_0 = "+
                cvm::to_str(cvm::fabs(L0+DL0 - L0_new)/L0, cvm::main()->cv_width, cvm::main()->cv_prec)+
                ", |(q_0+dq_0) - q_0^new| = "+
                cvm::to_str((Q0+DQ0 - Q0_new).norm(), cvm::main()->cv_width, cvm::main()->cv_prec)+
                "\n");
    }
  }
}

// Calculate the rotation, plus its derivatives

void colvarmodule::rotation::calc_optimal_rotation(
                                        std::vector<cvm::atom_pos> const &pos1,
                                        std::vector<cvm::atom_pos> const &pos2)
{
  C.reset();
  build_correlation_matrix(pos1, pos2);

  calc_optimal_rotation_impl();

  if (b_debug_gradients) {
    // debug_gradients<cvm::atom_pos, cvm::atom_pos, false>(*this, pos1, pos2, pos1.size(), pos2.size());
    const cvm::ag_vector_real_t pos1_soa = cvm::atom_group::pos_aos_to_soa(pos1);
    const cvm::ag_vector_real_t pos2_soa = cvm::atom_group::pos_aos_to_soa(pos2);
    debug_gradients(*this, pos1_soa, pos2_soa, pos1.size(), pos2.size());
  }
}

void colvarmodule::rotation::calc_optimal_rotation_soa(
  cvm::ag_vector_real_t const &pos1,
  cvm::ag_vector_real_t const &pos2,
  const size_t num_atoms_pos1,
  const size_t num_atoms_pos2) {
  C.reset();
  // build_correlation_matrix SOA
  // TODO: I expect std::span from C++20 to simplify this
  const cvm::real* pos1_x = pos1.data();
  const cvm::real* pos1_y = pos1_x + num_atoms_pos1;
  const cvm::real* pos1_z = pos1_y + num_atoms_pos1;
  const cvm::real* pos2_x = pos2.data();
  const cvm::real* pos2_y = pos2_x + num_atoms_pos2;
  const cvm::real* pos2_z = pos2_y + num_atoms_pos2;
  for (size_t i = 0; i < num_atoms_pos1; ++i) {
    C.xx += pos1_x[i] * pos2_x[i];
    C.xy += pos1_x[i] * pos2_y[i];
    C.xz += pos1_x[i] * pos2_z[i];
    C.yx += pos1_y[i] * pos2_x[i];
    C.yy += pos1_y[i] * pos2_y[i];
    C.yz += pos1_y[i] * pos2_z[i];
    C.zx += pos1_z[i] * pos2_x[i];
    C.zy += pos1_z[i] * pos2_y[i];
    C.zz += pos1_z[i] * pos2_z[i];
  }
  calc_optimal_rotation_impl();
  if (b_debug_gradients) debug_gradients(*this, pos1, pos2, num_atoms_pos1, num_atoms_pos2);
}

// Calculate the optimal rotation between two groups, and implement it
// as a quaternion.  Uses the method documented in: Coutsias EA,
// Seok C, Dill KA.  Using quaternions to calculate RMSD.  J Comput
// Chem. 25(15):1849-57 (2004) DOI: 10.1002/jcc.20110 PubMed: 15376254
void colvarmodule::rotation::calc_optimal_rotation_impl() {
  compute_overlap_matrix();

  // S_backup = S;
  std::memcpy(&S_backup[0][0], &S, 4*4*sizeof(cvm::real));

  if (b_debug_gradients) {
    cvm::matrix2d<cvm::real> S_backup_out(4, 4);
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        S_backup_out[i][j] = S_backup[i][j];
      }
    }
    cvm::log_static("S     = "+cvm::to_str(S_backup_out, cvm::main()->cv_width, cvm::main()->cv_prec)+"\n");
  }


#ifdef COLVARS_LAMMPS
  MathEigen::Jacobi<cvm::real,
                    cvm::real[4],
                    cvm::real[4][4]> *ecalc =
    reinterpret_cast<MathEigen::Jacobi<cvm::real,
                                       cvm::real[4],
                                       cvm::real[4][4]> *>(jacobi);

  int ierror = ecalc->Diagonalize(S, S_eigval, S_eigvec);
#else
  int ierror = NR::diagonalize_matrix(S, S_eigval, S_eigvec);
#endif
  if (ierror) {
    cvm::log_static("Failed to diagonalize the following overlapping matrix:\n");
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        cvm::log_static(cvm::to_str(S[i][j]) + " ");
      }
      cvm::log_static("\n");
    }
    cvm::log_static("The corresponding correlation matrix is:\n");
    cvm::log_static(" " + cvm::to_str(C.xx) + " " + cvm::to_str(C.xy) + " " + cvm::to_str(C.xz));
    cvm::log_static(" " + cvm::to_str(C.yx) + " " + cvm::to_str(C.yy) + " " + cvm::to_str(C.yz));
    cvm::log_static(" " + cvm::to_str(C.zx) + " " + cvm::to_str(C.zy) + " " + cvm::to_str(C.zz) + "\n");
    cvm::error_static("Too many iterations in jacobi diagonalization.\n"
               "This is usually the result of an ill-defined set of atoms for "
               "rotational alignment (RMSD, rotateReference, etc).\n");
  }
  q = cvm::quaternion{S_eigvec[0][0], S_eigvec[0][1], S_eigvec[0][2], S_eigvec[0][3]};

  if (cvm::rotation::monitor_crossings) {
    if (q_old.norm2() > 0.0) {
      q.match(q_old);
      if (q_old.inner(q) < (1.0 - crossing_threshold)) {
        cvm::log_static("Warning: one molecular orientation has changed by more than "+
                 cvm::to_str(crossing_threshold)+": discontinuous rotation ?\n");
      }
    }
    q_old = q;
  }
}

#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
namespace colvars_gpu {

rotation_gpu::rotation_gpu():
  d_S(nullptr), d_S_eigval(nullptr),
  d_S_eigvec(nullptr), tbcount(nullptr),
  d_q(nullptr), d_q_old(nullptr),
  discontinuous_rotation(nullptr),
  max_iteration_reached(nullptr), b_initialized(false),
  h_C(nullptr), h_S(nullptr), h_S_eigval(nullptr),
  h_S_eigvec(nullptr)
{}

rotation_gpu::~rotation_gpu() {
  colvarproxy* p = cvm::main()->proxy;
  p->deallocate_device(&d_S);
  p->deallocate_device(&d_S_eigval);
  p->deallocate_device(&d_S_eigvec);
  p->deallocate_device(&tbcount);
  p->deallocate_device(&d_q);
  p->deallocate_device(&d_q_old);
  p->deallocate_host(&discontinuous_rotation);
  p->deallocate_host(&max_iteration_reached);
  p->deallocate_host(&h_C);
  p->deallocate_host(&h_S);
  p->deallocate_host(&h_S_eigval);
  p->deallocate_host(&h_S_eigvec);
  b_initialized = false;
}

int rotation_gpu::init(/*const cudaStream_t& stream_in*/) {
  int error_code = COLVARS_OK;
  // stream = stream_in;
  colvarproxy* p = cvm::main()->proxy;
  if (!b_initialized) {
    error_code |= p->allocate_device(&d_S, 4 * 4);
    error_code |= p->allocate_device(&d_S_eigval, 4);
    error_code |= p->allocate_device(&d_S_eigvec, 4 * 4);
    error_code |= p->allocate_device(&tbcount, 1);
    error_code |= p->allocate_device(&d_q, 1);
    error_code |= p->allocate_device(&d_q_old, 1);
    error_code |= p->allocate_host(&discontinuous_rotation, 1);
    error_code |= p->allocate_host(&max_iteration_reached, 1);
    error_code |= p->allocate_host(&h_C, 1);
    error_code |= p->allocate_host(&h_S, 4 * 4);
    error_code |= p->allocate_host(&h_S_eigval, 4);
    error_code |= p->allocate_host(&h_S_eigvec, 4 * 4);
    error_code |= p->clear_device_array(tbcount, 1);
    max_iteration_reached[0] = 0;
    discontinuous_rotation[0] = 0;
    if (colvarmodule::rotation::monitor_crossings) {
      error_code |= p->clear_device_array(&d_q_old, 1);
    }
    cvm::main()->cite_feature("Optimal rotation via flexible fitting");
    b_initialized = true;
  }
  return error_code;
}

int rotation_gpu::add_optimal_rotation_nodes(
  cvm::real* const d_pos1,
  cvm::real* const d_pos2,
  const size_t num_atoms_pos1,
  const size_t num_atoms_pos2,
  // cudaGraphNode_t& node,
  cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map) {
  int error_code = COLVARS_OK;
  // Add memset nodes
  cudaGraphNode_t d_SSetNode;
  error_code |= colvars_gpu::add_clear_array_node(
    d_S, 4*4, d_SSetNode, graph, {});
  nodes_map["rotation_gpu_S_set"] = d_SSetNode;
  const cvm::real* d_pos1_x = d_pos1;
  const cvm::real* d_pos1_y = d_pos1_x + num_atoms_pos1;
  const cvm::real* d_pos1_z = d_pos1_y + num_atoms_pos1;
  const cvm::real* d_pos2_x = d_pos2;
  const cvm::real* d_pos2_y = d_pos2_x + num_atoms_pos2;
  const cvm::real* d_pos2_z = d_pos2_y + num_atoms_pos2;
  // Kernel node for building S
  cudaGraphNode_t build_S_node;
  std::vector<cudaGraphNode_t> dependencies;
  // The coordinates are not always moved to origin, so these dependencies are conditional
  colvars_gpu::prepare_dependencies(
    {{"read_positions_main", true},
     {"read_positions_fitting", true},
     {"move_to_origin_main", true},
     {"move_to_origin_fitting", true}}, dependencies, nodes_map,
     "build_overlapping_matrix");
  dependencies.push_back(d_SSetNode);
  error_code |= colvars_gpu::build_overlapping_matrix(
    d_pos1_x, d_pos1_y, d_pos1_z,
    d_pos2_x, d_pos2_y, d_pos2_z,
    d_S, d_S_eigvec, h_C, tbcount, num_atoms_pos1,
    build_S_node, graph, dependencies);
  nodes_map["build_overlapping_matrix"] = build_S_node;
  // Copy the overlapping matrix to host for CPU buffers
  cudaGraphNode_t copy_DtoH_S_node;
  error_code |= colvars_gpu::add_copy_node(
    d_S, h_S, 4*4, cudaMemcpyDeviceToHost, copy_DtoH_S_node,
    graph, {build_S_node});
  nodes_map["copy_DtoH_S"] = copy_DtoH_S_node;
  // Kernel node for eigensystem
  cudaGraphNode_t Jacobi4x4Node;
  error_code |= colvars_gpu::jacobi_4x4(
    d_S_eigvec, d_S_eigval,
    max_iteration_reached,
    d_q, colvarmodule::rotation::monitor_crossings,
    colvarmodule::rotation::crossing_threshold,
    d_q_old, discontinuous_rotation,
    Jacobi4x4Node, graph, {build_S_node});
  nodes_map["calc_optimal_rotation"] = Jacobi4x4Node;
  cudaGraphNode_t copy_DtoH_S_eigvec_node;
  error_code |= colvars_gpu::add_copy_node(
    d_S_eigvec, h_S_eigvec, 4*4,
    cudaMemcpyDeviceToHost, copy_DtoH_S_eigvec_node,
    graph, {Jacobi4x4Node});
  nodes_map["copy_DtoH_S_eigvec"] = copy_DtoH_S_eigvec_node;
  cudaGraphNode_t copy_DtoH_S_eigval_node;
  error_code |= colvars_gpu::add_copy_node(
    d_S_eigval, h_S_eigval, 4,
    cudaMemcpyDeviceToHost, copy_DtoH_S_eigval_node,
    graph, {Jacobi4x4Node});
  nodes_map["copy_DtoH_S_eigval"] = copy_DtoH_S_eigval_node;
  return error_code;
}

void rotation_gpu::after_sync_check() const {
  if (max_iteration_reached[0]) {
    cvm::error_static("Too many iterations in jacobi diagonalization.\n"
               "This is usually the result of an ill-defined set of atoms for "
               "rotational alignment (RMSD, rotateReference, etc).\n");
    max_iteration_reached[0] = 0;
  }
  if (colvarmodule::rotation::monitor_crossings) {
    if (discontinuous_rotation[0]) {
      cvm::main()->log("Warning: one molecular orientation has changed by more than "+
                cvm::to_str(colvarmodule::rotation::crossing_threshold)+
                ": discontinuous rotation ?\n");
    }
    discontinuous_rotation[0] = 0;
  }
}

void rotation_gpu::to_cpu(cvm::rotation& rot) const {
  std::memcpy(rot.get_S(), h_S, 4*4*sizeof(cvm::real));
  std::memcpy(rot.get_S_backup(), h_S, 4*4*sizeof(cvm::real));
  std::memcpy(rot.get_eigenvalues(), h_S_eigval, 4*sizeof(cvm::real));
  std::memcpy(rot.get_eigenvectors(), h_S_eigvec, 4*4*sizeof(cvm::real));
  *(rot.get_C()) = *h_C;
  rot.q = cvm::quaternion{h_S_eigvec[0], h_S_eigvec[1], h_S_eigvec[2], h_S_eigvec[3]};
}
}
#endif // defined(COLVARS_CUDA) || defined(COLVARS_HIP)
