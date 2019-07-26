using namespace Eigen;
using namespace std;

typedef Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> c_vector_t;
typedef Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> c_matrix_t;
typedef Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>> vector_t;
typedef Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> matrix_t;

class qpsolver {
  private:
    c_matrix_t& Qmat; // QP symmetric matrix
    c_matrix_t& Cmat; // equality constraint matrix
    c_matrix_t& Dmat; // inequality constraint matrix
    c_vector_t& qvec; // qp vector
    c_vector_t& cvec; // equality constraint vector
    c_vector_t& dvec; // inequality constraint vector
    c_vector_t& w0; // starting value
    unsigned int maxiter = 500;
    double tol = .5e-5
}

Eigen::VectorXd qpsolver::solve(void) {
  unsigned int n = this->Cmat.cols();
  vector_t w_best(n), w_prev(n), w_k(n);
  vector_t chi_next(n), xi_next(n);
  vector_t chi(n) = vector_t::Zero(n);
  vector_t chi_prev(n) = vector_t::Zero(n);
  vector_t xi(n) = vector_t::Zero(n);
  vector_t xi_prev(n) = vector_t::Zero(n);
  LLT<MatrixXd> lltOfQ(this->Q);
  matrix_t B(this->Cmat.rows() + this->Dmat.rows(), this->Cmat.cols());
  B << Cmat, Dmat;
  double LC = (B * lltOfQ.solve(B.transpose())).norm(), fac;
  w_prev = this->w0;
  unsigned int i = 1;
  while (true) {
    fac = (i - 1.)/(i + 2.);
    w_best = -lltOfQ.solve(q + Cmat.transpose() * xi + Dmat.transpose() * chi);
    w_k = w_best +  fac * (w_best - w_prev);
    xi_next = xi + fac * (xi - xi_prev) + (Cmat * w_k - cvec) / LC;
    chi_next = (chi + fac * (chi - chi_prev) + (Dmat * w_k - dvec) / LC).array().max(0);
    chi_prev = chi;
    xi_prev = xi;
    chi = chi_next;
    xi = xi_next;
    if(((w_best - w_prev).array().abs() <= tol * (w_best.array().abs() + w_prev.array().abs())).all())
      break;
    w_prev = w_best;
    ++i;
  }
  return w_best;
}
