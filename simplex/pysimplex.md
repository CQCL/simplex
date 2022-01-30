# pysimplex

This is a fast Clifford circuit simulator based on [BH21][1].

Known issues:

* Not yet optimized for sparsity.
* Measurements in the Y basis, while correct, are less efficient than they
  should be because the method described in the paper was not quite working.

Work is in progress on both these issues and the package will be updated in due
course.

Wheels are currently available for Linux, MacOS and Windows, for Python versions
3.8, 3.9 and 3.10.

Example usage:

```python
from pysimplex import Simplex
S = Simplex(2) # a 2-qubit system
S.H(0)
S.CX(0, 1)
b0 = S.MeasZ(0, coin=0) # 'coin' is an optional argument that fixes the result
assert b0 == 0 # without 'coin' it would be 0 or 1 with equal probability
b1 = S.MeasZ(1)
assert b0 == b1
assert not S.is_deterministic() # 'coin' is irrelevant here
```

The available operations are: `X`, `Y`, `Z`, `H`, `S`, `Sdg`, `CX`, `CZ`,
`MeasX`, `MeasY` and `MeasZ`.

[1]: https://arxiv.org/abs/2109.08629
