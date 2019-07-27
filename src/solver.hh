#ifndef SOLVER_H
#define SOLVER_H

#include <Eigen/Dense>
using namespace Eigen;

typedef const Eigen::Matrix<double, Eigen::Dynamic, 1> c_vector_t;
typedef const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> c_matrix_t;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_t;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_t;


class qpsolver {
  private:
    c_matrix_t& Qmat; // QP symmetric matrix
    c_vector_t& qvec; // qp vector
    c_matrix_t& Cmat; // equality constraint matrix
    c_vector_t& cvec; // equality constraint vector
    c_matrix_t& Dmat; // inequality constraint matrix
    c_vector_t& dvec; // inequality constraint vector
    c_vector_t& w0;   // starting value
    unsigned int maxiter;
    double tol;

  public:
    qpsolver(c_matrix_t& Qmat, c_vector_t& qvec, c_matrix_t& Cmat,
             c_vector_t& cvec, c_matrix_t& Dmat, c_vector_t& dvec,
             c_vector_t& w0, unsigned int maxiter = 500, double tol = .5e-5):
             Qmat(Qmat), qvec(qvec), Cmat(Cmat), cvec(cvec), Dmat(Dmat),
             dvec(dvec), w0(w0), maxiter(maxiter), tol(tol) {}
    vector_t solve(void) const;
};

vector_t qpsolver::solve(void) const {
  unsigned int n = Cmat.cols();
  vector_t w_best(n), w_prev(n), w_k(n);
  vector_t chi_next(n), xi_next(n);
  vector_t chi = vector_t::Zero(n);
  vector_t chi_prev = vector_t::Zero(n);
  vector_t xi = vector_t::Zero(n);
  vector_t xi_prev = vector_t::Zero(n);
  LLT<MatrixXd> lltOfQ(Qmat);
  matrix_t B(Cmat.rows() + Dmat.rows(), Cmat.cols());
  B << Cmat, Dmat;
  double LC = (B * lltOfQ.solve(B.transpose())).norm(), fac;
  w_prev = w0;
  for (unsigned int i = 1; i < maxiter; ++i) {
    fac = (i - 1.)/(i + 2.);
    w_best = -lltOfQ.solve(qvec + Cmat.transpose() * xi + Dmat.transpose() * chi);
    w_k = w_best + fac * (w_best - w_prev);
    xi_next = xi + fac * (xi - xi_prev) + (Cmat * w_k - cvec) / LC;
    chi_next = (chi + fac * (chi - chi_prev) + (Dmat * w_k - dvec) / LC).array().max(0);
    chi_prev = chi;
    xi_prev = xi;
    chi = chi_next;
    xi = xi_next;
    if(((w_best - w_prev).array().abs() <= tol * (w_best.array().abs() + w_prev.array().abs())).all())
      break;
    w_prev = w_best;
  }
  return w_best;
}
#endif
