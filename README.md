# lincon
a quadratic programming solver with support for linear equality and inequality constraints, i.e.

```
minimize (1/2)x^TQx + q^Tx
subject to Cx = c, Dx <= d
```

**convergence rate: O (1/k^2)**

**references**

* [Giselsson (2013)](https://www.sciencedirect.com/science/article/pii/S0005109813000101)
* [Nesterov (1983)](http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=dan&paperid=46009&option_lang=eng)

**documentation: [https://mirca.github.io/lincon](https://mirca.github.io/lincon)**

### Instalation

#### C++

just copy the file ``lincon/src/solver.h`` =)

#### Python

```
pip install lincon
```

### Basic usage

#### C++

```{c++}
#include "solver.h"
#include <iostream>

using namespace Eigen;
using namespace std;

int main() {
  Eigen::MatrixXd Qmat(2, 2), Dmat(2, 2);
  Eigen::MatrixXd Cmat(1, 2);
  Eigen::VectorXd qvec(2), cvec(1), dvec(2), w0(2);
  Qmat << 2.,.5,  .5,1.;
  w0 << .5, .5;
  qvec << 1.,1.;
  Cmat << 1.,1.;
  Dmat << -1.,0.,  0.,-1.;
  cvec << 1.;
  dvec << 0.,0.;
  qpsolver* my_solver = new qpsolver(Qmat, qvec, Cmat, cvec, Dmat, dvec, w0);
  std::cout << my_solver->solve() << std::endl;
}
```

#### Python

```{python}
import numpy as np
from lincon import qp

Qmat = np.array([[2., .5], [.5, 1]], dtype=np.float64)
w0 = .5 * np.ones(2, dtype=np.float64)
qvec = np.ones(2, dtype=np.float64)
Cmat = np.atleast_2d(np.ones(2, dtype=np.float64))
Dmat = np.array([[-1, 0], [0, -1]], dtype=np.float64)
cvec = np.ones(1, dtype=np.float64)
dvec = np.zeros(2, dtype=np.float64)

x = qp.solve(Qmat=Qmat, qvec=qvec, Cmat=Cmat, cvec=cvec,
             Dmat=Dmat, dvec=dvec, w0=w0)
```

### License

Copyright 2019 Ze Vinicius

This project is licensed under the terms of the MIT License.
