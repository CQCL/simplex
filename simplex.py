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
        self.Q = np.zeros(
            (n + 1, n + 1), dtype=int
        )  # rxr off-diagonal symmetric matrix (mod 2), extra row and col are for temporary use
        self.R = np.zeros(n + 1, dtype=int)  # diagonal of Q (mod 4)
        self.A = np.zeros(
            (n, n + 1), dtype=int
        )  # nxr matrix (mod 2), extra col is for temporary use
        self.b = np.zeros(n, dtype=int)  # nx1 column vector (mod 2)
        self.p = bidict()  # principal index map for A (col < r <--> row < n)
        self.deterministic = True

    def copy(self):
        other = QFE(self.n)
        other.r = self.r
        other.Q = np.copy(self.Q)
        other.R = np.copy(self.R)
        other.A = np.copy(self.A)
        other.b = np.copy(self.b)
        other.p = bidict(self.p)
        other.deterministic = self.deterministic
        return other

    def show(self):
        print("n =", self.n)
        print("r =", self.r)
        print("Q:")
        print(self.Q[: self.r, : self.r])
        print("R:")
        print(self.R[: self.r])
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
        assert all(self.R[i] in range(4) for i in range(self.r))
        assert all(
            self.Q[i, j] in range(2) for i in range(self.r) for j in range(self.r)
        )
        assert all(
            self.Q[i, j] == self.Q[j, i] for i in range(self.r) for j in range(self.r)
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

    def ReindexSubtColumn(self, k, c):
        assert k < self.r
        assert c < self.r
        if k == c:
            return
        self.A[:, k] ^= self.A[:, c]
        self.R[k] += 2 * self.Q[c, k]  # Q is symmetric
        self.R[k] %= 4
        self.Q[: self.r, k] ^= self.Q[: self.r, c]
        self.Q[k, : self.r] ^= self.Q[c, : self.r]

    def ReindexSwapColumns(self, k, c):
        assert k < self.r
        assert c < self.r
        if k == c:
            return
        self.A[:, [k, c]] = self.A[:, [c, k]]
        self.R[[k, c]] = self.R[[c, k]]
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

    def principate(self, j):
        c = self.p.inverse.get(j)
        if c is not None:
            self.ReselectPrincipalRow(j, c)
            if j != self.p[c]:
                c = None
        return c

    def decrement_r(self):
        if self.r - 1 in self.p:
            del self.p[self.r - 1]
        self.r -= 1

    def FixFinalBit(self, z):
        """Reduces r by 1."""
        assert z in [0, 1]
        r = self.r
        assert r > 0
        a = self.A[:, r - 1]
        u = self.R[r - 1]
        self.decrement_r()
        for i in range(self.r):
            self.R[i] += 2 * z * self.Q[i, self.r]
            self.R[i] %= 4
        self.b ^= z * a

    def ZeroColumnElim(self, c):
        """Reduces r by 1 or 2."""
        r = self.r
        assert r > 0
        assert c < r
        assert all(self.A[:, c] == 0)
        self.ReindexSwapColumns(c, r - 1)
        H = [h for h in range(r - 1) if self.Q[h, r - 1] == 1]
        u = self.R[r - 1]
        self.decrement_r()
        if u % 2 == 1:
            for h1 in H:
                for h2 in H:
                    self.Q[h1, h2] ^= 1
            for h in H:
                self.R[h] += u + 2
                self.R[h] %= 4
        else:
            if len(H) == 0:
                return
            else:
                l = H[0]
            for h in H[1:]:
                self.ReindexSubtColumn(h, l)
            self.ReindexSwapColumns(self.r - 1, l)
            self.FixFinalBit(u // 2)

    def SimulateX(self, j):
        self.b[j] ^= 1

    def SimulateZ(self, j):
        for i in range(self.r):
            self.R[i] += 2 * self.A[j, i]
            self.R[i] %= 4

    def SimulateY(self, j):
        self.SimulateZ(j)
        self.SimulateX(j)

    def SimulateH(self, j):
        c = self.principate(j)
        self.Q[self.r, : self.r] = self.A[j, : self.r]
        self.Q[: self.r, self.r] = self.A[j, : self.r]
        self.R[self.r] = 2 * self.b[j]
        self.A[j, : self.r] = 0
        self.A[:, self.r] = 0
        self.A[j, self.r] = 1
        self.p.inverse[j] = self.r
        self.b[j] = 0
        self.r += 1
        if c is not None:
            self.ZeroColumnElim(c)

    def SimulateS(self, j):
        H = [h for h in range(self.r) if self.A[j, h] == 1]
        for h1 in H:
            for h2 in H:
                self.Q[h1, h2] ^= 1
        v = 1 - 2 * self.b[j]
        for h in H:
            self.R[h] += v
            self.R[h] %= 4

    def SimulateSdg(self, j):
        H = [h for h in range(self.r) if self.A[j, h] == 1]
        v = 1 - 2 * self.b[j]
        for h1 in H:
            for h2 in H:
                self.Q[h1, h2] ^= 1
        for h in H:
            self.R[h] -= v
            self.R[h] %= 4

    def SimulateCZ(self, j, k):
        assert j != k
        H_j = set(h for h in range(self.r) if self.A[j, h] == 1)
        H_k = set(h for h in range(self.r) if self.A[k, h] == 1)
        for h1 in H_j:
            for h2 in H_k:
                self.Q[h1, h2] ^= 1
                self.Q[h2, h1] ^= 1
        for h in H_j & H_k:
            self.R[h] += 2
            self.R[h] %= 4
        v_k = 2 * self.b[k]
        for h in H_j:
            self.R[h] += v_k
            self.R[h] %= 4
        v_j = 2 * self.b[j]
        for h in H_k:
            self.R[h] += v_j
            self.R[h] %= 4

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
        c = self.principate(j)
        if (c is None) or any(self.Q[c, k] == 1 for k in range(self.r) if k != c):
            beta = self.toss_coin(coin)
        else:
            if self.R[c] == 0:
                return 0
            elif self.R[c] == 2:
                return 1
            else:
                beta = self.toss_coin(coin)
                self.R[c] = 2 * beta
                return beta
        self.Q[self.r, : self.r] = 0
        self.Q[: self.r, self.r] = 0
        for h in range(self.r):
            self.R[h] += 2 * beta * self.A[j, h]
            self.R[h] %= 4
        self.R[self.r] = 2 * beta
        self.A[j, : self.r] = 0
        self.A[:, self.r] = 0
        self.A[j, self.r] = 1
        self.p.inverse[j] = self.r
        self.r += 1
        self.b[j] = 0
        if c is not None:
            self.ZeroColumnElim(c)
        return beta

    def SimulateMeasY(self, j, coin=None):
        c = self.principate(j)
        if (c is None) or any(self.Q[c, k] == 1 for k in range(self.r) if k != c):
            beta = self.toss_coin(coin)
        else:
            if self.R[c] == 1:
                return 0
            elif self.R[c] == 3:
                return 1
            else:
                beta = self.toss_coin(coin)
                self.R[c] = 2 * beta + 1
                return beta
        H = [h for h in range(self.r) if self.A[j, h] == 1]
        self.Q[self.r, : self.r] = 0
        self.Q[: self.r, self.r] = 0
        v = 2 * self.b[j] + 2 * beta - 1
        for h1 in H:
            for h2 in H:
                self.Q[h1, h2] ^= 1
        for h in H:
            self.R[h] += v
            self.R[h] %= 4
        self.R[self.r] = 2 * beta + 1
        self.A[j, : self.r] = 0
        self.A[:, self.r] = 0
        self.A[j, self.r] = 1
        self.p.inverse[j] = self.r
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
