import quadprog
import numpy as np

Q = np.array([[2., .5], [.5, 1]], dtype=np.float64)
w0 = .5 * np.ones(2)
qvec = np.ones(2)
Cmat = np.atleast_2d(np.ones(2))
Dmat = np.array([[1, 0], [0, 1]], dtype=np.float64)
cvec = np.ones(1)
dvec = np.zeros(2)

Amat = np.vstack((Cmat, Dmat)).T
bvec = np.concatenate((cvec, dvec))
print(quadprog.solve_qp(Q, -qvec, Amat, bvec, meq = 1)[0])
