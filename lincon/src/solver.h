#ifndef SOLVER_H
#define SOLVER_H

#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

typedef const Eigen::Matrix<double, Eigen::Dynamic, 1> c_vector_t;
typedef const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> c_matrix_t;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_t;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_t;

// save the state of the dual variables
class dual_state {
  private:
    c_vector_t& dual_mu_prev;
    c_vector_t& dual_mu;
    c_vector_t& dual_lambda_prev;
    c_vector_t& dual_lambda;

  public:
    dual_state(c_vector_t& dual_mu_prev, c_vector_t& dual_mu,
               c_vector_t& dual_lambda_prev, c_vector_t& dual_lambda):
              dual_mu(dual_mu), dual_lambda(dual_lambda){}
    c_vector_t& get_dual_mu_prev(void) const {return dual_mu_prev;}
    c_vector_t& get_dual_mu(void) const {return dual_mu;}
    c_vector_t& get_dual_lambda_prev(void) const {return dual_lambda_prev;}
    c_vector_t& get_dual_lambda(void) const {return dual_lambda;}
};


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
    unsigned int niter; // number of iterations taken
    bool convergence;     // did it converged?
    void set_niter(unsigned int value) { niter = value; }
    void set_convergence(bool value) { convergece = value; }

  public:
    qpsolver(c_matrix_t& Qmat, c_vector_t& qvec, c_matrix_t& Cmat,
             c_vector_t& cvec, c_matrix_t& Dmat, c_vector_t& dvec,
             c_vector_t& w0, unsigned int maxiter = 500, double tol = .5e-5):
            Qmat(Qmat), qvec(qvec), Cmat(Cmat), cvec(cvec), Dmat(Dmat),
            dvec(dvec), w0(w0), maxiter(maxiter), tol(tol), niter(0),
            convergence(false) {}
    vector_t solve(void);
    c_matrix_t& get_Qmat(void) const { return Qmat;}
    c_vector_t& get_qvec(void) const { return qvec;}
    c_matrix_t& get_Cmat(void) const { return Cmat;}
    c_vector_t& get_cvec(void) const { return cvec;}
    c_matrix_t& get_Dmat(void) const { return Dmat;}
    c_vector_t& get_dvec(void) const { return dvec;}
    c_vector_t& get_w0(void) const { return w0;}
    unsigned int get_maxiter(void) const { return maxiter; }
    double get_tol(void) const { return tol; }
    unsigned int get_niter(void) const { return niter; }
    bool get_convergence(void) const { return convergence; }
    void set_maxiter(unsigned int value) { maxiter = value; }
    void set_tol(double value) { tol = value; }
};

vector_t qpsolver::solve(void) const {
  unsigned int n = Cmat.cols();
  unsigned int xi_len = Cmat.rows();
  unsigned int chi_len = Dmat.rows();
  vector_t w_best(n), w_prev(n), w_k(n);
  vector_t chi_next(chi_len), xi_next(xi_len);
  vector_t chi = vector_t::Zero(chi_len);
  vector_t chi_prev = vector_t::Zero(chi_len);
  vector_t xi = vector_t::Zero(xi_len);
  vector_t xi_prev = vector_t::Zero(xi_len);
  LLT<matrix_t> lltOfQ(Qmat);
  matrix_t B(Cmat.rows() + Dmat.rows(), Cmat.cols());
  B << Cmat, Dmat;
  double LC = (B * lltOfQ.solve(B.transpose())).norm(), fac;
  w_prev = w0;
  for (unsigned int i = 0; i < maxiter; ++i) {
    fac = (i - 1.)/(i + 2.);
    w_best = -lltOfQ.solve(qvec + Cmat.transpose() * xi + Dmat.transpose() * chi);
    w_k = w_best + fac * (w_best - w_prev);
    xi_next = xi + fac * (xi - xi_prev) + (Cmat * w_k - cvec) / LC;
    chi_next = (chi + fac * (chi - chi_prev) + (Dmat * w_k - dvec) / LC).array().max(0);
    chi_prev = chi;
    xi_prev = xi;
    chi = chi_next;
    xi = xi_next;
    if(((w_best - w_prev).array().abs() <= tol *
        (w_best.array().abs() + w_prev.array().abs())).all()) {
      set_convergence(true);
      set_niter(i);
      break;
    }
    w_prev = w_best;
  }
  return w_best;
}
#endif
