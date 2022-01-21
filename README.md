# Simplex

This software implements the methods described in
https://arxiv.org/abs/2109.08629 for simulating Clifford circuits.

## Reference implementation

A reference implementation in Python is in the `reference` directory. The
implementation is in `simplex.py`. It requires `numpy` and `bidict`:

```shell
pip install numpy bidict
```

To run the unit tests:

```shell
cd reference
./test-simplex
```

There is an option to perform some additional internal state validation during
testing:

```shell
cd reference
./test-simplex --validate
```

However, this requires that `sagemath` be installed on your system and that sage
scripts be executable with `/usr/bin/python`. (It uses sage for a matrix rank
computation.) It also slows down the tests considerably.

### API

The main class is the `Simplex` class, defined and documented in `simplex.py`.
Example usage:

```python
from simplex import Simplex
S = Simplex(2) # 2 qubits, initially in the all-zero state
S.H(0) # apply a Hadamard on qubit 0
S.CX(0,1) # apply a CNOT between qubits 0 and 1
b0 = S.MeasZ(0) # measure qubit 0 in the Z basis, return 0 or 1
b1 = S.MeasZ(1) # measure qubit 1 in the Z basis, return 0 or 1
assert b0 == b1 # it's a Bell state
assert not S.is_deterministic() # at least one non-deterministic measurement has been performed
```

By default, non-deterministic measurements use a PRNG to decide the result.
However, for reproduciblility and testing you may specify "coins" by passing an
extra argument to the measure  operation. For example:

```python
b0 = S.MeasZ(0, coin=0)
assert b0 == 0
```

Note that `is_deterministic()` will still return `False` regardless of whether a
coin was actually specified.

## C++ implementation

The C++ implementation is in the `simplex` directory. It is built using `cmake`. To build the shared library and the unit tests:

```shell
cd simplex
mkdir build
cd build
cmake ..
cmake --build .
```

Then to run the tests:

```shell
./test/simplex-test
```

### API

The C++ API is similar to the API of the Python reference implementation.

```cpp
#include <simplex.hpp>

Simplex S(2);
S.H(0);
S.CX(0, 1);
```

Measurements return integer values 0 or 1. There is also the same optional
argument to the measure operations for specifying coins:

```cpp
int b0 = S.MeasY(0, 1); // b0 will be 1
int b1 = S.MeasZ(1);
```

The internal state can be viewed by means of the `operator<<` method:

```cpp
#include <iostream>
// ...
std::cout << S;
```

This will print out the matrices _A_ and _Q_, the vector _b_, etc. Note that
_Q_ is stored as a 0,1-matrix representing the off-diagonal elements only; the
diagonal elements are represented by two 0,1-vectors _R0_ and _R1_ representing
the low and high bits respectively of the diagonal elements.