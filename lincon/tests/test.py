#import quadprog
import numpy as np
from lincon import qp

Qmat = np.array([[2., .5], [.5, 1]], dtype=np.float64)
w0 = .5 * np.ones(2, dtype=np.float64)
qvec = np.ones(2, dtype=np.float64)
Cmat = np.atleast_2d(np.ones(2, dtype=np.float64))
Dmat = np.array([[-1, 0], [0, -1]], dtype=np.float64)
cvec = np.ones(1, dtype=np.float64)
dvec = np.zeros(2, dtype=np.float64)

Amat = np.vstack((Cmat, -Dmat)).T
bvec = np.concatenate((cvec, dvec))
#print(quadprog.solve_qp(Qmat, -qvec, Amat, bvec, meq = 1)[0])

def test_simple_qp():
    x = qp.solve(Qmat=Qmat, qvec=qvec, Cmat=Cmat, cvec=cvec, Dmat=Dmat,
                 dvec=dvec, w0=w0, maxiter=500, tol=.1e-7)
    np.assert_almost_equal(x[0], [0.25, 0.75])
