#include "solver.h"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

py::tuple solve(c_matrix_t Qmat, c_vector_t qvec, c_matrix_t Cmat,
               c_vector_t cvec, c_matrix_t Dmat, c_vector_t dvec,
               c_vector_t w0, unsigned int maxiter, double tol) {
  qpsolver* qpsol = new qpsolver(Qmat, qvec, Cmat, cvec, Dmat, dvec, w0,
                                     maxiter, tol);
  py::tuple args = py::make_tuple(qpsol->solve(),
                                  qpsol->get_niter(),
                                  qpsol->get_convergence());
  return args;
}

PYBIND11_MODULE(qp, m) {
  m.doc() = "a collection of constrained quadratic program solvers";
  m.def("solve", &solve,
        py::arg("Qmat"), py::arg("qvec"), py::arg("Cmat"), py::arg("cvec"),
        py::arg("Dmat"), py::arg("dvec"), py::arg("w0"),
        py::arg("maxiter") = 100, py::arg("tol") = .5E-4,
        R"pbdoc(
          Solves quadratic program with linear equality and inequality constraints.
         )pbdoc"
        );
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
