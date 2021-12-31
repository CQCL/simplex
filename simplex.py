import numpy as np
from bidict import bidict
from random import randrange
import subprocess


class QFE:
    """Quadratic form expansion as in equation (8)"""

    def __init__(self, n):
        """n: number of qubits"""
        self.n = n  # n = number of rows of A and b
        self.r = 0  # n >= r = number of rows/cols of Q, number of cols of A
        self.g = 0  # exponent of tau
        self.Q = np.zeros(
            (n + 1, n + 1), dtype=int
        )  # rxr matrix (mod 4), extra row and col are for temporary use
        self.A = np.zeros(
            (n, n + 1), dtype=int
        )  # nxr matrix (mod 2), extra col is for temporary use
        self.b = np.zeros(n, dtype=int)  # nx1 column vector (mod 2)
        self.p = bidict()  # principal index map for A (col < r <--> row < n)
        self.deterministic = True

    def copy(self):
        other = QFE(self.n)
        other.r = self.r
        other.g = self.g
        other.Q = np.copy(self.Q)
        other.A = np.copy(self.A)
        other.b = np.copy(self.b)
        other.p = bidict(self.p)
        other.deterministic = self.deterministic
        return other

    def show(self):
        print("n =", self.n)
        print("r =", self.r)
        print("g =", self.g)
        print("Q:")
        print(self.Q[: self.r, : self.r])
        print("A:")
        print(self.A[:, : self.r])
        print("b:", self.b)

    def rankA(self):
        Astring = str(self.A[:, : self.r].tolist())
        out = subprocess.run(
            ["/usr/bin/python", "rank.sage", Astring], capture_output=True
        ).stdout
        return int(out.strip())

    def validate(self):
        assert self.r <= self.n
        assert all(self.Q[i, i] in range(4) for i in range(self.r))
        assert all(
            self.Q[i, j] in range(2)
            for i in range(self.r)
            for j in range(self.r)
            if i != j
        )
        assert all(
            self.A[i, j] in range(2) for i in range(self.n) for j in range(self.r)
        )
        assert all(self.b[i] in range(2) for i in range(self.n))
        assert self.rankA() == self.r
        assert all(
            (k in range(self.r) and v in range(self.n)) for k, v in self.p.items()
        )

    def toss_coin(self, coin=None):
        self.deterministic = False
        if coin is None:
            return randrange(2)
        else:
            assert coin in [0, 1]
            return coin

    def ReduceGramRowCol(self, c):
        assert c < self.r
        for k in range(self.r):
            if k == c:
                self.Q[k, k] %= 4
            else:
                self.Q[c, k] %= 2
                self.Q[k, c] %= 2

    def ReindexSubtColumn(self, k, c):
        assert k < self.r
        assert c < self.r
        if k == c:
            return
        self.A[:, k] ^= self.A[:, c]
        self.Q[: self.r, k] -= self.Q[: self.r, c]
        self.Q[k, : self.r] -= self.Q[c, : self.r]
        self.ReduceGramRowCol(k)

    def ReindexSwapColumns(self, k, c):
        assert k < self.r
        assert c < self.r
        if k == c:
            return
        self.A[:, [k, c]] = self.A[:, [c, k]]
        self.Q[:, [k, c]] = self.Q[:, [c, k]]
        self.Q[[k, c], :] = self.Q[[c, k], :]
        # Swap p[k] and p[c]:
        pk, pc = self.p.get(k), self.p.get(c)
        if pk is not None:
            del self.p[k]
        if pc is not None:
            del self.p[c]
        if pk is not None:
            self.p[c] = pk
        if pc is not None:
            self.p[k] = pc

    def MakePrincipal(self, c, j):
        assert c < self.r
        assert j < self.n
        if self.A[j, c] == 1:
            for k in range(self.r):
                if k != c and self.A[j, k] == 1:
                    self.ReindexSubtColumn(k, c)
            self.p[c] = j

    def ReselectPrincipalRow(self, j, c):
        # 0 <= c < r, 0 <= j < n or j is None
        assert c < self.r
        assert j is None or j < self.n
        n0 = None  # infinity
        for j1 in range(self.n):
            if j1 != j and self.A[j1, c] == 1:
                n1 = sum(self.A[j1, : self.r])
                if (n0 is None) or (n1 < n0):
                    j0 = j1
                    n0 = n1
        if n0 is not None:
            self.MakePrincipal(c, j0)

    def decrement_r(self):
        if self.r - 1 in self.p:
            del self.p[self.r - 1]
        self.r -= 1

    def FixFinalBit(self, z):
        """Reduces r by 1."""
        assert z in [0, 1]
        r = self.r
        assert r > 0
        a = np.copy(self.A[:, r - 1])
        q = np.copy(self.Q[: r - 1, r - 1])
        u = self.Q[r - 1, r - 1]
        self.decrement_r()
        for i in range(self.r):
            self.Q[i, i] += 2 * z * q[i]
            self.Q[i, i] %= 4
        self.b ^= z * a
        self.g += 2 * z * u

    def ZeroColumnElim(self, c):
        """Reduces r by 1 or 2."""
        r = self.r
        assert r > 0
        assert c < r
        assert all(self.A[:, c] == 0)
        self.ReindexSwapColumns(c, r - 1)
        q = np.copy(self.Q[: r - 1, [r - 1]])
        u = self.Q[r - 1, r - 1]
        self.decrement_r()
        if u % 2 == 1:
            self.Q[: self.r, : self.r] += (u - 2) * q @ q.transpose()
            self.Q[: self.r, : self.r] %= 4
            self.g -= u - 2
        else:
            l = 0
            while l < self.r:
                if q[l, 0] == 1:
                    break
                l += 1
            else:
                return
            for k in range(self.r):
                if k != l and q[k, 0] != 0:
                    self.ReindexSubtColumn(k, l)
            self.ReindexSwapColumns(self.r - 1, l)
            self.FixFinalBit(u // 2)

    def SimulateX(self, j):
        self.b[j] ^= 1

    def SimulateZ(self, j):
        self.g += 4 * self.b[j]
        self.g %= 8
        for i in range(self.r):
            self.Q[i, i] += 2 * self.A[j, i]
            self.Q[i, i] %= 4

    def SimulateY(self, j):
        self.g += 2
        self.g %= 8
        self.SimulateZ(j)
        self.SimulateX(j)

    def SimulateH(self, j):
        c = self.p.inverse.get(j)
        if c is not None:
            self.ReselectPrincipalRow(j, c)
            if j != self.p[c]:
                c = None
        a = np.copy(self.A[j, : self.r])
        self.A[j, : self.r] = 0
        self.A[:, self.r] = 0
        self.A[j, self.r] = 1
        self.p.inverse[j] = self.r
        self.Q[self.r, : self.r] = a
        self.Q[: self.r, self.r] = a
        self.Q[self.r, self.r] = 2 * self.b[j]
        self.b[j] = 0
        self.r += 1
        if c is not None:
            self.ZeroColumnElim(c)

    def SimulateS(self, j):
        a = np.copy(self.A[[j], : self.r])
        self.Q[: self.r, : self.r] += (1 - 2 * self.b[j]) * a.transpose() @ a
        for k in range(self.r):
            if self.A[j, k] == 1:
                self.ReduceGramRowCol(k)
        self.g += 2 * self.b[j]

    def SimulateSdg(self, j):
        a = np.copy(self.A[[j], : self.r])
        self.Q[: self.r, : self.r] -= (1 - 2 * self.b[j]) * a.transpose() @ a
        for k in range(self.r):
            if self.A[j, k] == 1:
                self.ReduceGramRowCol(k)
        self.g -= 2 * self.b[j]

    def SimulateCZ(self, j, k):
        assert j != k
        a_j = np.copy(self.A[[j], : self.r])
        a_k = np.copy(self.A[[k], : self.r])
        self.Q[: self.r, : self.r] += a_j.transpose() @ a_k + a_k.transpose() @ a_j
        for h in range(self.r):
            self.Q[h, h] += 2 * (self.b[k] * a_j[0, h] + self.b[j] * a_k[0, h])
        for h in range(self.r):
            if self.A[j, h] == 1 or self.A[k, h] == 1:
                self.ReduceGramRowCol(h)
        self.g += 4 * self.b[j] * self.b[k]

    def SimulateCX(self, h, j):
        assert h != j
        for k in range(self.r):
            if self.A[h, k] == 1:
                self.A[j, k] ^= 1
        self.b[j] ^= self.b[h]
        c = self.p.inverse.get(j)
        if c is not None:
            self.ReselectPrincipalRow(None, c)

    def SimulateMeasZ(self, j, coin=None):
        if all(self.A[j, : self.r] == 0):
            return self.b[j]
        else:
            beta = self.toss_coin(coin)
            k, n = None, None
            for k0 in range(self.r):
                if self.A[j, k0] == 1:
                    n0 = sum(self.A[:, k0])
                    if (n is None) or (n0 < n):
                        k = k0
                        n = n0
            self.ReindexSwapColumns(k, self.r - 1)
            self.MakePrincipal(self.r - 1, j)
            self.FixFinalBit(beta ^ self.b[j])
            return beta

    def SimulateMeasX(self, j, coin=None):
        c = self.p.inverse.get(j)
        if c is not None:
            self.ReselectPrincipalRow(j, c)
            if j != self.p[c]:
                c = None
        if (c is None) or any(self.Q[c, k] != 0 for k in range(self.r) if k != c):
            beta = self.toss_coin(coin)
        else:
            if self.Q[c, c] == 0:
                return 0
            elif self.Q[c, c] == 2:
                return 1
            else:
                beta = self.toss_coin(coin)
                self.Q[c, c] = 2 * beta
                return beta
        a = np.copy(self.A[j, : self.r])
        self.A[j, : self.r] = 0
        self.A[:, self.r] = 0
        self.A[j, self.r] = 1
        self.p.inverse[j] = self.r
        self.Q[self.r, : self.r + 1] = 0
        self.Q[: self.r + 1, self.r] = 0
        for h in range(self.r):
            self.Q[h, h] += 2 * beta * a[h]
        self.Q[self.r, self.r] = 2 * beta
        for k in range(self.r):
            if a[k] == 1:
                self.Q[k, k] %= 4
        self.r += 1
        self.b[j] = 0
        if c is not None:
            self.ZeroColumnElim(c)
        return beta

    def SimulateMeasY(self, j, coin=None):
        c = self.p.inverse.get(j)
        if c is not None:
            self.ReselectPrincipalRow(j, c)
            if j != self.p[c]:
                c = None
        if (c is None) or any(self.Q[c, k] != 0 for k in range(self.r) if k != c):
            beta = self.toss_coin(coin)
        else:
            if self.Q[c, c] == 1:
                return 0
            elif self.Q[c, c] == 3:
                return 1
            else:
                beta = self.toss_coin(coin)
                self.Q[c, c] = 2 * beta + 1
                return beta
        a = np.copy(self.A[[j], : self.r])
        self.A[j, : self.r] = 0
        self.A[:, self.r] = 0
        self.A[j, self.r] = 1
        self.p.inverse[j] = self.r
        self.Q[self.r, : self.r + 1] = 0
        self.Q[: self.r + 1, self.r] = 0
        self.Q[: self.r, : self.r] += (2 * self.b[j] + 2 * beta + 1) * a.transpose() @ a
        self.Q[self.r, self.r] = 2 * beta + 1
        for k in range(self.r):
            if a[0, k] == 1:
                self.ReduceGramRowCol(k)
        self.r += 1
        self.b[j] = 0
        if c is not None:
            self.ZeroColumnElim(c)
        return beta


class Simplex:
    """Clifford circuit simulator"""

    def __init__(self, n):
        """Initialize a simulator with `n` qubits."""
        self.n = n
        self.E = QFE(n)

    def copy(self):
        other = Simplex(self.n)
        other.E = self.E.copy()
        return other

    def X(self, j):
        """Apply an X gate to qubit `j`."""
        self.E.SimulateX(j)
        return self

    def Y(self, j):
        """Apply a Y gate to qubit `j`."""
        self.E.SimulateY(j)
        return self

    def Z(self, j):
        """Apply a Z gate to qubit `j`."""
        self.E.SimulateZ(j)
        return self

    def H(self, j):
        """Apply an H gate to qubit `j`."""
        self.E.SimulateH(j)
        return self

    def S(self, j):
        """Apply an S gate to qubit `j`."""
        self.E.SimulateS(j)
        return self

    def Sdg(self, j):
        """Apply an inverse S gate to qubit `j`."""
        self.E.SimulateSdg(j)
        return self

    def CX(self, j, k):
        """Apply a CX gate to qubits `j` and `k`."""
        self.E.SimulateCX(j, k)
        return self

    def CZ(self, j, k):
        """Apply a CZ gate to qubits `j` and `k`."""
        self.E.SimulateCZ(j, k)
        return self

    def MeasX(self, j, coin=None):
        """Measure qubit `j` in the X basis.

        If the outcome is non-deterministic, then the outcome is the value of `coin` if
        specified (must be 0 or 1), otherwise random."""
        return self.E.SimulateMeasX(j, coin=coin)

    def MeasY(self, j, coin=None):
        """Measure qubit `j` in the Y basis.

        If the outcome is non-deterministic, then the outcome is the value of `coin` if
        specified (must be 0 or 1), otherwise random."""
        return self.E.SimulateMeasY(j, coin=coin)

    def MeasZ(self, j, coin=None):
        """Measure qubit `j` in the Z basis.

        If the outcome is non-deterministic, then the outcome is the value of `coin` if
        specified (must be 0 or 1), otherwise random."""
        return self.E.SimulateMeasZ(j, coin=coin)

    def is_deterministic(self):
        """Report whether all measurements are deterministic."""
        return self.E.deterministic

    def show(self):
        self.E.show()

    def validate(self):
        self.E.validate()
