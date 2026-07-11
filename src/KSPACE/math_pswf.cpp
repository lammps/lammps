#include "math_pswf.h"

#include "math_const.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

using LAMMPS_NS::MathConst::MY_PI;

static constexpr int MAX_CHEB_ORDER = 40;

// start of legendre functions
static inline void legepol(double x, int n, double &pol, double &der)
{
  double pkm1 = 1.0;
  double pk = x;
  double pkp1;

  if (n == 0) {
    pol = 1.0;
    der = 0.0;
    return;
  }

  if (n == 1) {
    pol = x;
    der = 1.0;
    return;
  }

  pk = 1.0;
  pkp1 = x;

  for (int k = 1; k < n; ++k) {
    pkm1 = pk;
    pk = pkp1;
    pkp1 = ((2 * k + 1) * x * pk - k * pkm1) / (k + 1);
  }

  pol = pkp1;
  der = n * (x * pkp1 - pk) / (x * x - 1);
}

static inline void legetayl(double pol, double der, double x, double h, int n, int k, double &sum,
                            double &sumder)
{
  double done = 1.0;
  double q0 = pol;
  double q1 = der * h;
  double q2 = (2 * x * der - n * (n + done) * pol) / (1 - x * x);
  q2 = q2 * h * h / 2;

  sum = q0 + q1 + q2;
  sumder = q1 / h + q2 * 2 / h;

  if (k <= 2) return;

  double qi = q1;
  double qip1 = q2;

  for (int i = 1; i <= k - 2; ++i) {
    double d = 2 * x * (i + 1) * (i + 1) / h * qip1 - (n * (n + done) - i * (i + 1)) * qi;
    d = d / (i + 1) / (i + 2) * h * h / (1 - x * x);
    double qip2 = d;

    sum += qip2;
    sumder += d * (i + 2) / h;

    qi = qip1;
    qip1 = qip2;
  }
}

/* Constructs Gaussian quadrature of order n.
   itype=1 => both roots (ts) and weights (whts) are computed.
   itype=0 => only roots (ts) are computed. */
static inline void legerts(int itype, int n, double *ts, double *whts)
{
  int k = 30;
  double d = 1.0;
  double d2 = d + 1.0e-24;
  if (d2 != d) { k = 54; }

  int half = n / 2;
  int ifodd = n - 2 * half;
  double pi_val = atan(1.0) * 4.0;
  double h = pi_val / (2.0 * n);

  /* Initial approximations (for i >= n/2+1) */
  int ii = 0;
  for (int i = 1; i <= n; i++) {
    if (i < (n / 2 + 1)) { continue; }
    ii++;
    double t = (2.0 * i - 1.0) * h;
    ts[ii - 1] = -cos(t);
  }

  /* Start from center: find roots one by one via Newton updates */
  double pol = 1.0, der = 0.0;
  double x0 = 0.0;
  legepol(x0, n, pol, der);
  double x1 = ts[0];

  int n2 = (n + 1) / 2;
  double pol3 = pol, der3 = der;

  for (int kk = 1; kk <= n2; kk++) {
    if ((ifodd == 1) && (kk == 1)) {
      ts[kk - 1] = x0;
      if (itype > 0) { whts[kk - 1] = der; }
      x0 = x1;
      x1 = ts[kk];
      pol3 = pol;
      der3 = der;
      continue;
    }

    /* Newton iteration */
    int ifstop = 0;
    for (int i = 1; i <= 10; i++) {
      double hh = x1 - x0;

      legetayl(pol3, der3, x0, hh, n, k, pol, der);
      x1 = x1 - pol / der;

      if (fabs(pol) < 1.0e-12) { ifstop++; }
      if (ifstop == 3) { break; }
    }

    ts[kk - 1] = x1;
    if (itype > 0) { whts[kk - 1] = der; }

    x0 = x1;
    x1 = ts[kk];
    pol3 = pol;
    der3 = der;
  }

  /* Mirror roots around 0: fill second half of ts[] */
  for (int i = n2; i >= 1; i--) { ts[i - 1 + half] = ts[i - 1]; }
  for (int i = 1; i <= half; i++) { ts[i - 1] = -ts[n - i]; }
  if (itype <= 0) { return; }

  /* Mirror weights similarly */
  for (int i = n2; i >= 1; i--) { whts[i - 1 + half] = whts[i - 1]; }
  for (int i = 1; i <= half; i++) { whts[i - 1] = whts[n - i]; }

  /* Compute final weights = 2 / (1 - ts[i]^2) / (der[i]^2) */
  for (int i = 0; i < n; i++) {
    double tmp = 1.0 - ts[i] * ts[i];
    whts[i] = 2.0 / tmp / (whts[i] * whts[i]);
  }
}

static inline void legepols(double x, int n, double *pols)
{
  double pkm1 = 1.0;
  double pk = x;

  if (n == 0) {
    pols[0] = 1.0;
    return;
  }

  if (n == 1) {
    pols[0] = 1.0;
    pols[1] = x;
    return;
  }

  pols[0] = 1.0;
  pols[1] = x;

  for (int k = 1; k < n; ++k) {
    double pkp1 = ((2 * k + 1) * x * pk - k * pkm1) / (k + 1);
    pols[k + 1] = pkp1;
    pkm1 = pk;
    pk = pkp1;
  }
}

// TODO: legepols() is not tested yet.
// only itype !=2 is tested.
static inline void legeexps(int itype, int n, double *x, std::vector<std::vector<double>> &u,
                            std::vector<std::vector<double>> &v, double *whts)
{
  int itype_rts = (itype > 0) ? 1 : 0;

  // Call legerts to construct the nodes and weights of the n-point Gaussian quadrature
  legerts(itype_rts, n, x, whts);

  // If itype is not 2, return early
  if (itype != 2) return;

  // Construct the matrix of values of the Legendre polynomials at these nodes
  for (int i = 0; i < n; ++i) {
    std::vector<double> pols(n);
    legepols(x[i], n - 1, pols.data());
    for (int j = 0; j < n; ++j) { u[j][i] = pols[j]; }
  }

  // Transpose u to get v
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) { v[i][j] = u[j][i]; }
  }

  // Construct the inverse u, converting the values of a function at Gaussian nodes into the
  // coefficients of a Legendre expansion of that function
  for (int i = 0; i < n; ++i) {
    double d = 1.0 * (2 * (i + 1) - 1) / 2;
    for (int j = 0; j < n; ++j) { u[i][j] = v[j][i] * whts[j] * d; }
  }
}

static inline void legeexev(double x, double &val, const double *pexp, int n)
{
  double pjm2 = 1.0;
  double pjm1 = x;

  val = pexp[0] * pjm2 + pexp[1] * pjm1;

  for (int j = 2; j <= n; ++j) {
    double pj = ((2 * j - 1) * x * pjm1 - (j - 1) * pjm2) / j;
    val += pexp[j] * pj;
    pjm2 = pjm1;
    pjm1 = pj;
  }
}

static inline void legeFDER(double x, double &val, double &der, const double *pexp, int n)
{
  double pjm2 = 1.0;
  double pjm1 = x;
  double derjm2 = 0.0;
  double derjm1 = 1.0;

  val = pexp[0] * pjm2 + pexp[1] * pjm1;
  der = pexp[1];

  for (int j = 2; j <= n; ++j) {
    double pj = ((2 * j - 1) * x * pjm1 - (j - 1) * pjm2) / j;
    val += pexp[j] * pj;

    double derj = (2 * j - 1) * (pjm1 + x * derjm1) - (j - 1) * derjm2;
    derj /= j;
    der += pexp[j] * derj;

    pjm2 = pjm1;
    pjm1 = pj;
    derjm2 = derjm1;
    derjm1 = derj;
  }
}

// Function to compute the Legendre polynomial p_n(x) and its derivative p'_n(x)
static inline void legendre(int n, double x, double &pn, double &pn_prime)
{
  if (n == 0) {
    pn = 1.0;
    pn_prime = 0.0;
    return;
  }

  if (n == 1) {
    pn = x;
    pn_prime = 1.0;
    return;
  }

  double pn_minus1 = 1.0;    // P_0(x)
  double pn_minus2 = 0.0;    // P_-1(x)
  pn = x;                    // P_1(x)

  for (int k = 2; k <= n; ++k) {
    pn_minus2 = pn_minus1;
    pn_minus1 = pn;
    pn = ((2.0 * k - 1.0) * x * pn_minus1 - (k - 1.0) * pn_minus2) / k;
  }

  pn_prime = n * (x * pn - pn_minus1) / (x * x - 1.0);
}

// Function to compute the nodes and weights of the n-point Gaussian quadrature
static inline void gaussian_quadrature(int n, double *nodes, double *weights)
{
  static std::vector<std::vector<double>> precomp_nodes(1000);
  static std::vector<std::vector<double>> precomp_weights(1000);

  {    // Precompute
    if (precomp_nodes[n].size() == 0 && precomp_weights[n].size() == 0) {
#pragma omp critical(GAUSS_QUAD)
      if (precomp_nodes[n].size() == 0 && precomp_weights[n].size() == 0) {
        std::vector<double> nodes1(n);
        std::vector<double> weights1(n);

        if (true) {    // usr original code
          legerts(1, n, nodes1.data(), weights1.data());
        } else {    // self implementation test code
          const double tolerance = 1e-16;
          const int max_iterations = 100;
          for (int i = 0; i < (n + 1) / 2; ++i) {
            double x = std::cos(MY_PI * (i + 0.75) / (n + 0.5));
            double pn, pn_prime;

            for (int iter = 0; iter < max_iterations; ++iter) {
              legendre(n, x, pn, pn_prime);
              double delta_x = -pn / pn_prime;
              x += delta_x;

              if (std::abs(delta_x) < tolerance) { break; }
            }

            nodes1[i] = x;
            nodes1[n - 1 - i] = -x;
            weights1[i] = 2.0 / ((1.0 - x * x) * pn_prime * pn_prime);
            weights1[n - 1 - i] = weights1[i];
          }
        }

        precomp_nodes[n].swap(nodes1);
        precomp_weights[n].swap(weights1);
      }
    }
  }

  memcpy(nodes, precomp_nodes[n].data(), n * sizeof(double));
  memcpy(weights, precomp_weights[n].data(), n * sizeof(double));
}
// end of legendre functions

// start of monomial utils

static inline void cheb_nodes_1d(int order, std::vector<double> &nodes, double a, double b)
{
  if ((int) nodes.size() != order) nodes.resize(order);
  for (int i = 0; i < order; i++) {
    nodes[i] = -cos((i + (double) 0.5) * MY_PI / order) * (double) 0.5 + (double) 0.5;
    nodes[i] = nodes[i] * (b - a) + a;
  }
}

static inline void monomial_interp_1d(int nnodes, std::vector<double> &fn_v,
                                      std::vector<double> &coeff, double a = 0, double b = 1)
{
  std::vector<double> x;
  cheb_nodes_1d(nnodes, x, a, b);

  auto multiply_x = [](const std::vector<double> &p, double x) {
    std::vector<double> r(p.size() + 1, double(0));
    for (std::size_t i = 0; i < p.size(); ++i) {
      r[i] += -x * p[i];
      r[i + 1] += p[i];
    }
    return r;
  };

  int dof = fn_v.size() / nnodes;

  std::vector<double> newton_coeffs = fn_v;

  if ((int)coeff.size() != dof * nnodes) coeff.resize((std::size_t)dof * nnodes);
  coeff.assign((std::size_t)dof * nnodes, 0.0);

  // Turner, L. R.. “Inverse of the Vandermonde matrix with applications.” (1966)
  for (int idof = 0; idof < dof; ++idof) {
    int offset = idof * nnodes;
    for (int j = 1; j < nnodes; ++j)
      for (int i = nnodes - 1; i >= j; --i)
        newton_coeffs[offset + i] =
            (newton_coeffs[offset + i] - newton_coeffs[offset + i - 1]) / (x[i] - x[i - j]);
    coeff[offset] = newton_coeffs[offset];
  }
  std::vector<double> basis{1.0};
  for (int j = 1; j < nnodes; ++j) {
    basis = multiply_x(basis, x[j - 1]);
    for (int idof = 0; idof < dof; idof++) {
      int offset = idof * nnodes;
      for (std::size_t m = 0; m < basis.size(); ++m)
        coeff[offset + m] += newton_coeffs[offset + j] * basis[m];
    }
  }
}

static inline void cheb_basis_1d(int order, const std::vector<double> &x, std::vector<double> &y,
                                 double a = 0, double b = 1)
{
  int n = x.size();
  if ((int) y.size() != order * n) y.resize((std::size_t)order * n);

  if (order > 0) {
    for (int i = 0; i < n; i++) y[i] = 1.0;
  }
  if (order > 1) {
    for (int i = 0; i < n; i++) y[i + n] = x[i] * 2 / (b - a) - 2 * a / (b - a) - 1;
  }
  for (int i = 2; i < order; i++) {
    for (int j = 0; j < n; j++) {
      y[i * n + j] = 2 * y[n + j] * y[i * n - 1 * n + j] - y[i * n - 2 * n + j];
    }
  }
}

static void cheb_interp_1d(int order, std::vector<double> &fn_v, std::vector<double> &coeff)
{
  static std::vector<std::vector<double>> precomp(1000);
  {    // Precompute
    if (precomp[order].size() == 0) {
#pragma omp critical(CHEB_BASIS_APPROX)
      if (precomp[order].size() == 0) {
        std::vector<double> x, p;
        cheb_nodes_1d(order, x, -1, 1);
        cheb_basis_1d(order, x, p, -1, 1);
        precomp[order].swap(p);
      }
    }
  }
  std::vector<double> &p = precomp[order];

  int dof = fn_v.size() / order;
  if ((int) coeff.size() != dof * order) coeff.resize((std::size_t)dof * order);

  const double inv_order = 1.0 / static_cast<double>(order);
  const double two_inv_order = 2.0 * inv_order;
  for (int id = 0; id < dof; ++id) {
    const int base = id * order;
    for (int k = 0; k < order; ++k) {
      double sum = 0.0;
      const double *pk = &p[(std::size_t)k * order];
      for (int j = 0; j < order; ++j) sum += fn_v[base + j] * pk[j];
      coeff[base + k] = (k == 0) ? (sum * inv_order) : (sum * two_inv_order);
    }
  }
}
// end of monomial utils

static inline void prosinin(double c, const double *ts, const double *whts, const double *fs,
                            double x, int n, double &rint, double &derrint)
{
  rint = 0.0;
  derrint = 0.0;

  for (int i = 0; i < n; ++i) {
    double diff = x - ts[i];
    double sin_term = sin(c * diff);
    double cos_term = cos(c * diff);

    rint += whts[i] * fs[i] * sin_term / diff;

    derrint += whts[i] * fs[i] / (diff * diff) * (c * diff * cos_term - sin_term);
  }
}

static inline void prolcoef(double rlam, int k, double c, double &alpha0, double &beta0,
                            double &gamma0, double &alpha, double &beta, double &gamma)
{
  double d = k * (k - 1);
  d = d / (2 * k + 1) / (2 * k - 1);
  double uk = d;

  d = (k + 1) * (k + 1);
  d = d / (2 * k + 3);
  double d2 = k * k;
  d2 = d2 / (2 * k - 1);
  double vk = (d + d2) / (2 * k + 1);

  d = (k + 1) * (k + 2);
  d = d / (2 * k + 1) / (2 * k + 3);
  double wk = d;

  alpha = -c * c * uk;
  beta = rlam - k * (k + 1) - c * c * vk;
  gamma = -c * c * wk;

  alpha0 = uk;
  beta0 = vk;
  gamma0 = wk;
}

static inline void prolmatr(double *as, double *bs, double *cs, int n, double c, double rlam,
                            int ifsymm, int ifodd)
{
  double done = 1.0;
  double half = done / 2.0;
  int k = 0;

  if (ifodd > 0) {
    for (int k0 = 1; k0 <= n + 2; k0 += 2) {
      k++;
      double alpha0, beta0, gamma0, alpha, beta, gamma;
      prolcoef(rlam, k0, c, alpha0, beta0, gamma0, alpha, beta, gamma);

      as[k - 1] = alpha;
      bs[k - 1] = beta;
      cs[k - 1] = gamma;

      if (ifsymm != 0) {
        if (k0 > 1) { as[k - 1] = as[k - 1] / std::sqrt(k0 - 2 + half) * std::sqrt(k0 + half); }
        cs[k - 1] = cs[k - 1] * std::sqrt(k0 + half) / std::sqrt(k0 + half + 2);
      }
    }
  } else {
    for (int k0 = 0; k0 <= n + 2; k0 += 2) {
      k++;
      double alpha0, beta0, gamma0, alpha, beta, gamma;
      prolcoef(rlam, k0, c, alpha0, beta0, gamma0, alpha, beta, gamma);

      as[k - 1] = alpha;
      bs[k - 1] = beta;
      cs[k - 1] = gamma;

      if (ifsymm != 0) {
        if (k0 != 0) { as[k - 1] = as[k - 1] / std::sqrt(k0 - 2 + half) * std::sqrt(k0 + half); }
        cs[k - 1] = cs[k - 1] * std::sqrt(k0 + half) / std::sqrt(k0 + half + 2);
      }
    }
  }
}

static inline void prolql1(int n, double *d, double *e, int &ierr)
{
  ierr = 0;
  if (n == 1) return;

  for (int i = 1; i < n; ++i) { e[i - 1] = e[i]; }
  e[n - 1] = 0.0;

  for (int l = 0; l < n; ++l) {
    int j = 0;
    while (true) {
      int m;
      for (m = l; m < n - 1; ++m) {
        double tst1 = std::abs(d[m]) + std::abs(d[m + 1]);
        double tst2 = tst1 + std::abs(e[m]);
        if (tst2 == tst1) break;
      }

      if (m == l) break;
      if (j == 30) {
        ierr = l + 1;
        return;
      }
      ++j;

      double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
      double r = std::sqrt(g * g + 1.0);
      g = d[m] - d[l] + e[l] / (g + std::copysign(r, g));
      double s = 1.0;
      double c = 1.0;
      double p = 0.0;

      for (int i = m - 1; i >= l; --i) {
        double f = s * e[i];
        double b = c * e[i];
        r = std::sqrt(f * f + g * g);
        e[i + 1] = r;
        if (r == 0.0) {
          d[i + 1] -= p;
          e[m] = 0.0;
          break;
        }
        s = f / r;
        c = g / r;
        g = d[i + 1] - p;
        r = (d[i] - g) * s + 2.0 * c * b;
        p = s * r;
        d[i + 1] = g + p;
        g = c * r - b;
      }

      if (r == 0.0) break;
      d[l] -= p;
      e[l] = g;
      e[m] = 0.0;
    }

    if (l == 0) continue;
    for (int i = l; i > 0; --i) {
      if (d[i] >= d[i - 1]) break;
      std::swap(d[i], d[i - 1]);
    }
  }
}

static inline void prolfact(double *a, double *b, double *c, int n, double *u, double *v, double *w)
{
  // Eliminate down
  for (int i = 0; i < n - 1; ++i) {
    double d = c[i + 1] / a[i];
    a[i + 1] -= b[i] * d;
    u[i] = d;
  }

  // Eliminate up
  for (int i = n - 1; i > 0; --i) {
    double d = b[i - 1] / a[i];
    v[i] = d;
  }

  // Scale the diagonal
  double done = 1.0;
  for (int i = 0; i < n; ++i) { w[i] = done / a[i]; }
}

static inline void prolsolv(const double *u, const double *v, const double *w, int n, double *rhs)
{
  // Eliminate down
  for (int i = 0; i < n - 1; ++i) { rhs[i + 1] -= u[i] * rhs[i]; }

  // Eliminate up
  for (int i = n - 1; i > 0; --i) { rhs[i - 1] -= rhs[i] * v[i]; }

  // Scale
  for (int i = 0; i < n; ++i) { rhs[i] *= w[i]; }
}

static inline void prolfun0(int &ier, int n, double c, double *as, double *bs, double *cs,
                            double *xk, double *u, double *v, double *w, double eps, int &nterms,
                            double &rkhi)
{
  ier = 0;
  double delta = 1.0e-8;
  int ifsymm = 1;
  int numit = 4;
  double rlam = 0;
  int ifodd = -1;

  prolmatr(as, bs, cs, n, c, rlam, ifsymm, ifodd);

  prolql1(n / 2, bs, as, ier);
  if (ier != 0) {
    ier = 2048;
    return;
  }

  rkhi = -bs[n / 2 - 1];
  rlam = -bs[n / 2 - 1] + delta;

  std::fill(xk, xk + n, 1.0);

  prolmatr(as, bs, cs, n, c, rlam, ifsymm, ifodd);

  prolfact(bs, cs, as, n / 2, u, v, w);

  for (int ijk = 0; ijk < numit; ++ijk) {
    prolsolv(u, v, w, n / 2, xk);

    double d = 0;
    for (int j = 0; j < n / 2; ++j) { d += xk[j] * xk[j]; }

    d = std::sqrt(d);
    for (int j = 0; j < n / 2; ++j) { xk[j] /= d; }

    double err = 0;
    for (int j = 0; j < n / 2; ++j) {
      err += (as[j] - xk[j]) * (as[j] - xk[j]);
      as[j] = xk[j];
    }
    err = std::sqrt(err);
  }

  double half = 0.5;
  for (int i = 0; i < n / 2; ++i) {
    if (std::abs(xk[i]) > eps) nterms = i + 1;
    xk[i] *= std::sqrt(i * 2 + half);
    cs[i] = xk[i];
  }

  int j = 0;
  for (int i = 0; i <= nterms; ++i) {
    xk[j++] = cs[i];
    xk[j++] = 0;
  }

  nterms *= 2;
}

static inline void prolps0i(int &ier, double c, double *w, int lenw, int &nterms, int &ltot,
                            double &rkhi)
{
  static const std::array<int, 20> ns = {48,  64,  80,  92,  106, 120, 130, 144, 156, 168,
                                         178, 190, 202, 214, 224, 236, 248, 258, 268, 280};

  double eps = 1.0e-16;
  int n = static_cast<int>(c * 3);
  n = n / 2;

  int i = static_cast<int>(c / 10);
  if (i <= 19) n = ns[i];

  ier = 0;
  int ixk = 1;
  int lxk = n + 2;

  int ias = ixk + lxk;
  int las = n + 2;

  int ibs = ias + las;
  int lbs = n + 2;

  int ics = ibs + lbs;
  int lcs = n + 2;

  int iu = ics + lcs;
  int lu = n + 2;

  int iv = iu + lu;
  int lv = n + 2;

  int iw = iv + lv;
  int lw = n + 2;

  ltot = iw + lw;

  if (ltot >= lenw) {
    ier = 512;
    return;
  }

  prolfun0(ier, n, c, w + ias - 1, w + ibs - 1, w + ics - 1, w + ixk - 1, w + iu - 1, w + iv - 1,
           w + iw - 1, eps, nterms, rkhi);

  if (ier != 0) return;
}

static inline void prol0ini(int &ier, double c, double *w, double &rlam20, double &rkhi, int lenw,
                            int &keep, int &ltot)
{
  ier = 0;
  double thresh = 45;
  int iw = 11;
  w[0] = iw + 0.1;
  w[8] = thresh;

  // Create the data to be used in the evaluation of the function ψ^c_0(x) for x ∈ [-1,1]
  int nterms = 0;
  prolps0i(ier, c, w + iw - 1, lenw, nterms, ltot, rkhi);

  if (ier != 0) return;

  // If c > thresh, do not prepare data for the evaluation of ψ^c_0 outside the interval [-1,1]
  if (c >= thresh) {
    w[7] = c;
    w[4] = nterms + 0.1;
    keep = nterms + 3;
    return;
  }

  // Create the data to be used in the evaluation of the function ψ^c_0(x) for x outside the
  // interval [-1,1]
  int ngauss = nterms * 2;
  int lw = nterms + 2;
  int its = iw + lw;
  int lts = ngauss + 2;
  int iwhts = its + lts;
  int lwhts = ngauss + 2;
  int ifs = iwhts + lwhts;
  int lfs = ngauss + 2;

  keep = ifs + lfs;
  if (keep > ltot) ltot = keep;
  if (keep >= lenw) {
    ier = 1024;
    return;
  }

  w[1] = its + 0.1;
  w[2] = iwhts + 0.1;
  w[3] = ifs + 0.1;

  int itype = 1;
  std::vector<std::vector<double>> u;
  std::vector<std::vector<double>> v;
  legeexps(itype, ngauss, w + its - 1, u, v, w + iwhts - 1);

  // Evaluate the prolate function at the Gaussian nodes
  for (int i = 0; i < ngauss; ++i) {
    legeexev(w[its + i - 1], w[ifs + i - 1], w + iw - 1, nterms - 1);
  }

  // Calculate the eigenvalue corresponding to ψ^c_0
  double rlam = 0;
  double x0 = 0;
  double f0;
  legeexev(x0, f0, w + iw - 1, nterms - 1);
  double der;
  prosinin(c, w + its - 1, w + iwhts - 1, w + ifs - 1, x0, ngauss, rlam, der);

  rlam = rlam / f0;
  rlam20 = rlam;

  w[4] = nterms + 0.1;
  w[5] = ngauss + 0.1;
  w[6] = rlam;
  w[7] = c;
}

static inline void prol0eva(double x, const double *w, double &psi0, double &derpsi0)
{
  int iw = static_cast<int>(w[0]);
  int its = static_cast<int>(w[1]);
  int iwhts = static_cast<int>(w[2]);
  int ifs = static_cast<int>(w[3]);

  int nterms = static_cast<int>(w[4]);
  int ngauss = static_cast<int>(w[5]);
  double rlam = w[6];
  double c = w[7];
  double thresh = w[8];

  if (std::abs(x) > 1) {
    if (c >= thresh - 1.0e-10) {
      psi0 = 0;
      derpsi0 = 0;
      return;
    }

    prosinin(c, &w[its - 1], &w[iwhts - 1], &w[ifs - 1], x, ngauss, psi0, derpsi0);
    psi0 /= rlam;
    derpsi0 /= rlam;
    return;
  }

  legeFDER(x, psi0, derpsi0, &w[iw - 1], nterms - 2);
  // to match chebfun psi0, needs a factor of sqrt(2)
  psi0 = sqrt(2.0) * psi0;
  derpsi0 = sqrt(2.0) * derpsi0;
}

static inline void prol0int0r(const double *w, double r, double &val)
{
  static int npts = 200;
  static int itype = 1;
  double derpsi0;
  static std::vector<double> xs(npts, 0), ws(npts, 0), fvals(npts, 0);
  static int need_init = 1;
  std::vector<std::vector<double>> u;
  std::vector<std::vector<double>> v;

  // since xs, ws, fval of size 200 are static
  // only need to get nodes and weights once
  if (need_init) {
#pragma omp critical(PROL0INT0R)
    if (need_init) {
      legeexps(itype, npts, xs.data(), u, v, ws.data());
      need_init = 0;
    }
  }

  // Scale the nodes and weights to [0, r]
  double xs_r;
  for (int i = 0; i < npts; ++i) {
    xs_r = (xs[i] + 1) * r / 2;
    prol0eva(xs_r, w, fvals[i], derpsi0);
  }

  val = 0;
  for (int i = 0; i < npts; ++i) { val += ws[i] * r / 2 * fvals[i]; }
}

namespace LAMMPS_NS::MathPSWF {

struct Prolate0Fun {
  Prolate0Fun() = default;

  inline Prolate0Fun(double c_, int lenw_) : c(c_), lenw(lenw_)
  {
    int ier;
    workarray.resize(lenw);
    prol0ini(ier, c, workarray.data(), rlam20, rkhi, lenw, keep, ltot);
    //if (ier) error->all(FLERR,"Unable to init Prolate0Fun");
  }

  // evaluate prolate0 function val and derivative
  [[nodiscard]] inline std::pair<double, double> eval_val_derivative(double x) const
  {
    double psi0, derpsi0;
    prol0eva(x, workarray.data(), psi0, derpsi0);
    return {psi0, derpsi0};
  }

  // evaluate prolate0 function value
  [[nodiscard]] inline double eval_val(double x) const
  {
    auto [val, dum] = eval_val_derivative(x);
    return val;
  }

  // evaluate prolate0 function derivative
  [[nodiscard]] inline double eval_derivative(double x) const
  {
    auto [dum, der] = eval_val_derivative(x);
    return der;
  }

  // int_0^r prolate0(x) dx
  [[nodiscard]] inline double int_eval(double r) const
  {
    double val;
    prol0int0r(workarray.data(), r, val);
    return val;
  }

  double c;
  int lenw, keep, ltot;
  std::vector<double> workarray;
  double rlam20, rkhi;
};

void prolc180(double eps, double &c)
{
  static const std::array<double, 180> cs = {
      0.43368E-16, 0.10048E+01, 0.17298E+01, 0.22271E+01, 0.26382E+01, 0.30035E+01, 0.33409E+01,
      0.36598E+01, 0.39658E+01, 0.42621E+01, 0.45513E+01, 0.48347E+01, 0.51136E+01, 0.53887E+01,
      0.56606E+01, 0.59299E+01, 0.61968E+01, 0.64616E+01, 0.67247E+01, 0.69862E+01, 0.72462E+01,
      0.75049E+01, 0.77625E+01, 0.80189E+01, 0.82744E+01, 0.85289E+01, 0.87826E+01, 0.90355E+01,
      0.92877E+01, 0.95392E+01, 0.97900E+01, 0.10040E+02, 0.10290E+02, 0.10539E+02, 0.10788E+02,
      0.11036E+02, 0.11284E+02, 0.11531E+02, 0.11778E+02, 0.12024E+02, 0.12270E+02, 0.12516E+02,
      0.12762E+02, 0.13007E+02, 0.13251E+02, 0.13496E+02, 0.13740E+02, 0.13984E+02, 0.14228E+02,
      0.14471E+02, 0.14714E+02, 0.14957E+02, 0.15200E+02, 0.15443E+02, 0.15685E+02, 0.15927E+02,
      0.16169E+02, 0.16411E+02, 0.16652E+02, 0.16894E+02, 0.17135E+02, 0.17376E+02, 0.17617E+02,
      0.17858E+02, 0.18098E+02, 0.18339E+02, 0.18579E+02, 0.18819E+02, 0.19059E+02, 0.19299E+02,
      0.19539E+02, 0.19778E+02, 0.20018E+02, 0.20257E+02, 0.20496E+02, 0.20736E+02, 0.20975E+02,
      0.21214E+02, 0.21452E+02, 0.21691E+02, 0.21930E+02, 0.22168E+02, 0.22407E+02, 0.22645E+02,
      0.22884E+02, 0.23122E+02, 0.23360E+02, 0.23598E+02, 0.23836E+02, 0.24074E+02, 0.24311E+02,
      0.24549E+02, 0.24787E+02, 0.25024E+02, 0.25262E+02, 0.25499E+02, 0.25737E+02, 0.25974E+02,
      0.26211E+02, 0.26448E+02, 0.26685E+02, 0.26922E+02, 0.27159E+02, 0.27396E+02, 0.27633E+02,
      0.27870E+02, 0.28106E+02, 0.28343E+02, 0.28580E+02, 0.28816E+02, 0.29053E+02, 0.29289E+02,
      0.29526E+02, 0.29762E+02, 0.29998E+02, 0.30234E+02, 0.30471E+02, 0.30707E+02, 0.30943E+02,
      0.31179E+02, 0.31415E+02, 0.31651E+02, 0.31887E+02, 0.32123E+02, 0.32358E+02, 0.32594E+02,
      0.32830E+02, 0.33066E+02, 0.33301E+02, 0.33537E+02, 0.33773E+02, 0.34008E+02, 0.34244E+02,
      0.34479E+02, 0.34714E+02, 0.34950E+02, 0.35185E+02, 0.35421E+02, 0.35656E+02, 0.35891E+02,
      0.36126E+02, 0.36362E+02, 0.36597E+02, 0.36832E+02, 0.37067E+02, 0.37302E+02, 0.37537E+02,
      0.37772E+02, 0.38007E+02, 0.38242E+02, 0.38477E+02, 0.38712E+02, 0.38947E+02, 0.39181E+02,
      0.39416E+02, 0.39651E+02, 0.39886E+02, 0.40120E+02, 0.40355E+02, 0.40590E+02, 0.40824E+02,
      0.41059E+02, 0.41294E+02, 0.41528E+02, 0.41763E+02, 0.41997E+02, 0.42232E+02, 0.42466E+02,
      0.42700E+02, 0.42935E+02, 0.43169E+02, 0.43404E+02, 0.43638E+02, 0.43872E+02, 0.44107E+02,
      0.44341E+02, 0.44575E+02, 0.44809E+02, 0.45044E+02, 0.45278E+02};

  if (eps < 1.0e-18) eps = 1e-18;
  double d = -log10(eps);
  int i = static_cast<int>(d * 10 + 0.1);
  c = cs[i - 1];
}

/*
evaluate prolate0c derivative at x, i.e., \psi_0^c(x)
*/
double prolate0_eval_derivative(double c, double x)
{
  static std::unordered_map<double, Prolate0Fun> prolate0_funcs_cache;
  if (prolate0_funcs_cache.find(c) == prolate0_funcs_cache.end()) {
#pragma omp critical(PROLATE0_EVAL)
    if (prolate0_funcs_cache.find(c) == prolate0_funcs_cache.end()) {
      prolate0_funcs_cache.emplace(c, Prolate0Fun(c, 10000));
    }
  }
  return prolate0_funcs_cache[c].eval_derivative(x);
}

/*
evaluate prolate0c at x, i.e., \psi_0^c(x)
*/
double prolate0_eval(double c, double x)
{
  static std::unordered_map<double, Prolate0Fun> prolate0_funcs_cache;
  if (prolate0_funcs_cache.find(c) == prolate0_funcs_cache.end()) {
#pragma omp critical(PROLATE0_EVAL)
    if (prolate0_funcs_cache.find(c) == prolate0_funcs_cache.end()) {
      prolate0_funcs_cache.emplace(c, Prolate0Fun(c, 10000));
    }
  }
  return prolate0_funcs_cache[c].eval_val(x);
}

/*
evaluate prolate0c function integral of \int_0^r \psi_0^c(x) dx
*/
double prolate0_int_eval(double c, double r)
{
  static std::unordered_map<double, Prolate0Fun> prolate0_funcs_cache;
  if (prolate0_funcs_cache.find(c) == prolate0_funcs_cache.end()) {
#pragma omp critical(PROLATE0_INT_EVAL)
    if (prolate0_funcs_cache.find(c) == prolate0_funcs_cache.end()) {
      prolate0_funcs_cache.emplace(c, Prolate0Fun(c, 10000));
    }
  }
  return prolate0_funcs_cache[c].int_eval(r);
}
// end of prolate functions

// start of approximation functions
void force_poly(double tol_coeff, const double &c, std::vector<double> &coeffs)
{
  double c0 = prolate0_int_eval(c, 1.0);

  int order = MAX_CHEB_ORDER;
  std::vector<double> nodes;
  cheb_nodes_1d(order, nodes, -1.0, 1.0);

  auto f = [](double c0, double c, double x) {
    x = (x + 1.0) / 2.0;
    double val = prolate0_int_eval(c, x) / c0 - x / c0 * prolate0_eval(c, x);
    val = 1 - val;
    return val;
  };

  int dof = 1;
  std::vector<double> fn_v((std::size_t)dof * order);
  for (int idof = 0; idof < dof; idof++) {
    for (int i = 0; i < order; i++) { fn_v[idof * order + i] = f(c0, c, nodes[i]); }
  }

  std::vector<double> cheb_coeff;
  cheb_interp_1d(order, fn_v, cheb_coeff);
  int max_order = -1;
  // filter chebyshev coefficients
  double max_coeffs = 0.0;
  for (int idof = 0; idof < dof; idof++) {
    max_coeffs = 0.0;
    for (int i = 0; i < order; i++) {
      max_coeffs = std::max<double>(max_coeffs, std::abs(cheb_coeff[idof * order + i]));
    }
    for (int i = 0; i < order; i++) {
      if (std::abs(cheb_coeff[idof * order + i]) > tol_coeff * max_coeffs) {
        max_order = std::max<int>(max_order, i + 1);
      }
    }
  }

  coeffs.resize((std::size_t)dof * max_order, 0.0);

  std::vector<double> coeffs_tmp(max_order);
  int nnodes = (int) max_order;    // * 1.75;
  cheb_nodes_1d(nnodes, nodes, -1.0, 1.0);
  fn_v.resize((std::size_t)dof * nnodes);

  for (int i = 0; i < nnodes; i++) { fn_v[i] = f(c0, c, nodes[i]); }

  monomial_interp_1d(nnodes, fn_v, coeffs_tmp, -1.0, 1.0);
  coeffs.resize(max_order, 0.0);

  for (int i = 0; i < max_order; i++) { coeffs[i] = coeffs_tmp[i]; }
}

void energy_poly(double tol_coeff, const double &c, std::vector<double> &coeffs)
{
  double c0 = prolate0_int_eval(c, 1.0);

  int order = MAX_CHEB_ORDER;
  std::vector<double> nodes;
  cheb_nodes_1d(order, nodes, -1.0, 1.0);

  auto f = [](double c0, double c, double x) {
    x = (x + 1.0) / 2.0;
    double val = prolate0_int_eval(c, x) / c0;
    val = 1 - val;
    return val;
  };

  int dof = 1;
  std::vector<double> fn_v((std::size_t)dof * order);
  for (int idof = 0; idof < dof; idof++) {
    for (int i = 0; i < order; i++) { fn_v[idof * order + i] = f(c0, c, nodes[i]); }
  }

  std::vector<double> cheb_coeff;
  cheb_interp_1d(order, fn_v, cheb_coeff);
  int max_order = -1;
  // filter chebyshev coefficients
  double max_coeffs = 0.0;
  for (int idof = 0; idof < dof; idof++) {
    max_coeffs = 0.0;
    for (int i = 0; i < order; i++) {
      max_coeffs = std::max<double>(max_coeffs, std::abs(cheb_coeff[idof * order + i]));
    }
    for (int i = 0; i < order; i++) {
      if (std::abs(cheb_coeff[idof * order + i]) > tol_coeff * max_coeffs) {
        max_order = std::max<int>(max_order, i + 1);
      }
    }
  }

  coeffs.resize((std::size_t)dof * max_order, 0.0);

  std::vector<double> coeffs_tmp((std::size_t)dof * max_order);
  int nnodes = (int) max_order;    // * 1.75;
  cheb_nodes_1d(nnodes, nodes, -1.0, 1.0);
  fn_v.resize((std::size_t)dof * nnodes);
  for (int i = 0; i < nnodes; i++) fn_v[i] = f(c0, c, nodes[i]);

  monomial_interp_1d(nnodes, fn_v, coeffs_tmp, -1.0, 1.0);

  for (int i = 0; i < max_order; i++) coeffs[i] = coeffs_tmp[i];
}

void fourier_poly(double tol_coeff, const double &c, double &lambda,
                  std::vector<double> &coeffs)
{
  double c0 = prolate0_int_eval(c, 1.0);

  int quad_npts = 200;
  std::vector<double> xs(quad_npts, 0), ws(quad_npts, 0);
  gaussian_quadrature(quad_npts, xs.data(), ws.data());
  lambda = 0.0;
  for (int i = 0; i < quad_npts; i++) {
    lambda += ws[i] * prolate0_eval(c, xs[i]) * std::cos(c * xs[i] * 0.5);
  }
  lambda /= prolate0_eval(c, 0.5);

  int order = MAX_CHEB_ORDER;
  std::vector<double> nodes;
  cheb_nodes_1d(order, nodes, -1.0, 1.0);

  auto f = [](double lambda, double c0, double c, double x) {
    x = (x + 1.0) / 2.0;
    double val = lambda * prolate0_eval(c, x) / c0;
    return val;
  };

  int dof = 1;
  std::vector<double> fn_v((std::size_t)dof * order);
  for (int idof = 0; idof < dof; idof++) {
    for (int i = 0; i < order; i++) { fn_v[idof * order + i] = f(lambda, c0, c, nodes[i]); }
  }

  std::vector<double> cheb_coeff;
  cheb_interp_1d(order, fn_v, cheb_coeff);
  int max_order = -1;
  // filter chebyshev coefficients
  double max_coeffs = 0.0;
  for (int idof = 0; idof < dof; idof++) {
    max_coeffs = 0.0;
    for (int i = 0; i < order; i++) {
      max_coeffs = std::max<double>(max_coeffs, std::abs(cheb_coeff[idof * order + i]));
    }

    for (int i = 0; i < order; i++) {
      if (std::abs(cheb_coeff[idof * order + i]) > tol_coeff * max_coeffs) {
        max_order = std::max<int>(max_order, i + 1);
      }
    }
  }

  coeffs.resize((std::size_t)dof * max_order, 0.0);

  std::vector<double> coeffs_tmp((std::size_t)dof * max_order);
  int nnodes = (int) max_order;    // * 1.75;
  cheb_nodes_1d(nnodes, nodes, -1.0, 1.0);
  fn_v.resize((std::size_t)dof * nnodes);
  for (int i = 0; i < nnodes; i++) fn_v[i] = f(lambda, c0, c, nodes[i]);

  monomial_interp_1d(nnodes, fn_v, coeffs_tmp, -1.0, 1.0);

  for (int i = 0; i < max_order; i++) coeffs[i] = coeffs_tmp[i];
}

void spread_fourier_poly(double tol_coeff, const double &c, double &lambda,
                         std::vector<double> &coeffs)
{
  int quad_npts = 200;
  std::vector<double> xs(quad_npts, 0), ws(quad_npts, 0);
  gaussian_quadrature(quad_npts, xs.data(), ws.data());
  lambda = 0.0;
  for (int i = 0; i < quad_npts; i++) {
    lambda += ws[i] * prolate0_eval(c, xs[i]) * std::cos(c * xs[i] * 0.5);
  }
  lambda /= prolate0_eval(c, 0.5);

  int order = MAX_CHEB_ORDER;
  std::vector<double> nodes;
  cheb_nodes_1d(order, nodes, -1.0, 1.0);

  auto f = [](double lambda, double c, double x) {
    x = (x + 1.0) / 2.0;
    double val = lambda * prolate0_eval(c, x);
    return val;
  };

  int dof = 1;
  std::vector<double> fn_v((std::size_t)dof * order);
  for (int idof = 0; idof < dof; idof++) {
    for (int i = 0; i < order; i++) { fn_v[idof * order + i] = f(lambda, c, nodes[i]); }
  }

  std::vector<double> cheb_coeff;
  cheb_interp_1d(order, fn_v, cheb_coeff);
  int max_order = -1;
  // filter chebyshev coefficients
  double max_coeffs = 0.0;
  for (int idof = 0; idof < dof; idof++) {
    max_coeffs = 0.0;
    for (int i = 0; i < order; i++) {
      max_coeffs = std::max<double>(max_coeffs, std::abs(cheb_coeff[idof * order + i]));
    }

    for (int i = 0; i < order; i++) {
      if (std::abs(cheb_coeff[idof * order + i]) > tol_coeff * max_coeffs) {
        max_order = std::max<int>(max_order, i + 1);
      }
    }
  }

  coeffs.resize((std::size_t)dof * max_order, 0.0);
  std::vector<double> coeffs_tmp((std::size_t)dof * max_order);
  int nnodes = (int) max_order;    // * 1.75;
  //monomial_nodes_1d(nnodes, nodes, 0, 1);
  cheb_nodes_1d(nnodes, nodes, -1.0, 1.0);
  fn_v.resize((std::size_t)dof * nnodes);
  for (int i = 0; i < nnodes; i++) { fn_v[i] = f(lambda, c, nodes[i]); }

  //monomial_interp_1d(max_order, nnodes, fn_v, coeffs_tmp);
  monomial_interp_1d(nnodes, fn_v, coeffs_tmp, -1.0, 1.0);

  //for (int i = 0; i < max_order; i++) { coeffs[max_order - i - 1] = coeffs_tmp[i]; }
  for (int i = 0; i < max_order; i++) { coeffs[i] = coeffs_tmp[i]; }
}

void spread_real_poly(int P, double tol_coeff, const double &c,
                      std::vector<double> &coeffs)
{
  int order = MAX_CHEB_ORDER;
  std::vector<double> nodes;
  cheb_nodes_1d(order, nodes, -0.5, 0.5);

  auto f = [](int P, int iP, double c, double x) {
    double arg = x - P / 2.0 + iP + 0.5;
    arg /= P / 2.0;
    double val = prolate0_eval(c, arg);
    return val;
  };
  int dof = P;
  std::vector<double> fn_v((std::size_t)dof * order);
  for (int idof = 0; idof < dof; idof++) {
    for (int i = 0; i < order; i++) { fn_v[idof * order + i] = f(P, idof, c, nodes[i]); }
  }

  std::vector<double> cheb_coeff;
  cheb_interp_1d(order, fn_v, cheb_coeff);
  int max_order = -1;
  // filter chebyshev coefficients
  double max_coeffs = 0.0;
  for (int idof = 0; idof < dof; idof++) {
    max_coeffs = 0.0;
    for (int i = 0; i < order; i++) {
      max_coeffs = std::max<double>(max_coeffs, std::abs(cheb_coeff[idof * order + i]));
    }
    for (int i = 0; i < order; i++) {
      if (std::abs(cheb_coeff[idof * order + i]) > tol_coeff * max_coeffs) {
        max_order = std::max<int>(max_order, i + 1);
      }
    }
  }

  coeffs.resize((std::size_t)dof * max_order, 0.0);

  int nnodes = (int) max_order;    // * 1.75;
  //monomial_nodes_1d(nnodes, nodes, -0.5, 0.5);
  cheb_nodes_1d(nnodes, nodes, -0.5, 0.5);    // use Chebyshev nodes for better accuracy
  fn_v.resize((std::size_t)dof * nnodes);          // sample Chebyshev points > order (polynomial interpolation)
  for (int iP = 0; iP < P; iP++) {    // loop over P
    for (int i = 0; i < nnodes; i++) {    // loop over nodes
      fn_v[iP * nnodes + i] =
          f(P, iP, c, nodes[i]);    // evaluate prolate0 function at each LS point
    }
  }

  std::vector<double> coeffs_tmp((std::size_t)P * order);    // coefficients for each P, order polynomial
  //monomial_interp_1d(max_order, nnodes, fn_v, coeffs_tmp, -0.5,
  //0.5);    // interpolate each polynomial for each P
  monomial_interp_1d(nnodes, fn_v, coeffs_tmp, -0.5,
                     0.5);    // interpolate each polynomial for each P

  // copy coeffs_tmp to coeffs
  for (int i = 0; i < P; i++) {
    for (int j = 0; j < max_order; j++) {
      //coeffs[j * P + i] = coeffs_tmp[i * max_order + max_order - j - 1];
      coeffs[j * P + i] = coeffs_tmp[i * max_order + j];
    }
  }
}
}    // namespace MathPSWF
