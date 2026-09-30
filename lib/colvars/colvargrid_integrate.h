#ifndef COLVARGRID_INTEGRATE_H
#define COLVARGRID_INTEGRATE_H

#include "colvargrid.h"

/// Integrate (1D, 2D or 3D) gradients

class colvargrid_integrate : public colvar_grid_scalar
{
  public:

  colvargrid_integrate();

  virtual ~colvargrid_integrate()
  {}

  /// Constructor from a vector of colvars + gradient grid
  colvargrid_integrate(std::vector<colvar *> &colvars,
                      std::shared_ptr<colvar_grid_gradient> gradients);

  /// Constructor from a gradient grid (for processing grid files without a Colvars config)
  colvargrid_integrate(std::shared_ptr<colvar_grid_gradient> gradients);

  /// \brief Calculate potential from divergence (in 2D); return number of steps
  int integrate(const int itmax, const cvm::real & tol, cvm::real & err, bool verbose = true);

  /// \brief Update matrix containing divergence and boundary conditions
  /// based on new gradient point value, in neighboring bins
  void update_div_neighbors(const std::vector<int> &ix);

  /// \brief Update matrix containing divergence and boundary conditions
  /// called by update_div_neighbors and by colvarbias_abf::adiabatic_reweighting_update_gradient_pmf
  void update_div_local(const std::vector<int> &ix);

  /// \brief Set matrix containing divergence and boundary conditions
  /// based on complete gradient grid
  void set_div();

  /// \brief Add constant to potential so that its minimum value is zero
  /// Useful e.g. for output
  inline void set_zero_minimum() {
    add_constant(-1.0 * minimum_value());
  }

  /// \brief Flag requesting the use of a smoothed version of the gradient (default: false)
  bool b_smoothed;


  protected:

  // Reference to gradient grid
  std::shared_ptr<colvar_grid_gradient> gradients;

  /// Array holding divergence + boundary terms (modified Neumann) if not periodic
  std::vector<cvm::real> divergence;

//   std::vector<cvm::real> inv_lap_diag; // Inverse of the diagonal of the Laplacian; for conditioning

  /// Obtain the gradient vector at given location ix, if available
  /// or zero if it is on the edge of the gradient grid
  /// ix gets wrapped in PBC
  void get_grad(cvm::real * g, std::vector<int> &ix);

  /// \brief Solve linear system based on CG, valid for symmetric matrices only
  void nr_linbcg_sym(const std::vector<cvm::real> &b, std::vector<cvm::real> &x,
                     const cvm::real &tol, const int itmax, int &iter, cvm::real &err);

  /// l2 norm of a vector
  cvm::real l2norm(const std::vector<cvm::real> &x);

  /// Multiplication by sparse matrix representing Lagrangian (or its transpose)
  void atimes(const std::vector<cvm::real> &x, std::vector<cvm::real> &r);

//   /// Inversion of preconditioner matrix
//   void asolve(const std::vector<cvm::real> &b, std::vector<cvm::real> &x);
};

#endif
