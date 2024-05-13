//*****************************************************************
// Iterative template routine -- CG
//
// CG solves the symmetric positive definite linear
// system Ax=b using the Conjugate Gradient method.
//
// CG follows the algorithm described on p. 15 in the
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
//*****************************************************************

  /**
   *  @class  CG
   *  @brief  Base class for solving the linear system Ax=b using the Conjugate Gradient method
   */

template < class Matrix, class Vector, class DataVector, class Preconditioner, class Real >
int
CG(const Matrix &A, Vector &x, const DataVector &b, const Preconditioner &M, int &max_iter, Real &tol) {
  Real resid;
  DenseVector<Real> p, z, q;
  Real alpha, beta, rho, rho_1(0);
  DenseVector<Real> tmp;
  tmp.reset(b.size());

  p.reset(b.size());
  z.reset(b.size());
  q.reset(b.size());

  Real normb = b.norm();
  DenseVector<Real> r;
  tmp = A*x;
  r = b - tmp;
  // Implicit assumption that only diagonal matrices are being used for preconditioning
  Preconditioner Minv = M.inv();

  if (normb == 0.0)
    normb = 1;

  if ((resid = r.norm() / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (int i = 0; i < max_iter; i++) {
    z = Minv*r;
    rho = r.dot(z);

    if (i == 0)
      p = z;
    else {
      beta = rho / rho_1;
      tmp = p*beta;
      p = z + tmp;
    }

    q = A*p;
    alpha = rho / p.dot(q);

    x += p*alpha;
    r -= q*alpha;

    if ((resid = r.norm() / normb) <= tol)
    {
      tol = resid;
      max_iter = i+1;
      return 0;
    }
    rho_1 = rho;
  }
  tol = resid;
  return 1;
}
