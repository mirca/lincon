#include "solver.hh"
#include <iostream>

using namespace Eigen;
using namespace std;

int main() {
  Eigen::MatrixXd Qmat(2, 2), Dmat(2, 2);
  Eigen::MatrixXd Cmat(2, 2);
  Eigen::VectorXd qvec(2), cvec(2), dvec(2), w0(2);
  Qmat << 2.,.5,  .5,1;
  w0 << .5, .5;
  qvec << 1.,1.;
  Cmat << 1.,1,  0., 0.;
  Dmat << -1.,0.,  0.,-1.;
  cvec << 1.,0.;
  dvec << 0.,0.;
  qpsolver* my_solver = new qpsolver(Qmat, qvec, Cmat, cvec, Dmat, dvec, w0,
                                     100, .5e-7);
  std::cout << my_solver->solve() << std::endl;
}
