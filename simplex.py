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
        self.R0 = np.zeros(n + 1, dtype=int)  # diagonal of Q (mod 4), lsb (0 or 1)
        self.R1 = np.zeros(n + 1, dtype=int)  # diagonal of Q (mod 4), msb (0 or 1)
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
        other.R0 = np.copy(self.R0)
        other.R1 = np.copy(self.R1)
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
        print("R0:")
        print(self.R0[: self.r])
        print("R1:")
        print(self.R1[: self.r])
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
        assert all(self.R0[i] in range(2) for i in range(self.r))
        assert all(self.R1[i] in range(2) for i in range(self.r))
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
        self.R1[k] ^= self.Q[c, k]  # NB Q is symmetric
        self.Q[: self.r, k] ^= self.Q[: self.r, c]
        self.Q[k, : self.r] ^= self.Q[c, : self.r]

    def ReindexSwapColumns(self, k, c):
        assert k < self.r
        assert c < self.r
        if k == c:
            return
        self.A[:, [k, c]] = self.A[:, [c, k]]
        self.R0[[k, c]] = self.R0[[c, k]]
        self.R1[[k, c]] = self.R1[[c, k]]
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
        assert self.r > 0
        self.b ^= z & self.A[:, self.r - 1]
        self.decrement_r()
        self.R1[: self.r] ^= z & self.Q[: self.r, self.r]

    def ZeroColumnElim(self, c):
        """Reduces r by 1 or 2."""
        assert c < self.r
        assert all(self.A[:, c] == 0)
        self.ReindexSwapColumns(c, self.r - 1)
        H = [h for h in range(self.r - 1) if self.Q[h, self.r - 1] == 1]
        u0, u1 = self.R0[self.r - 1], self.R1[self.r - 1]
        self.decrement_r()
        if u0 == 1:
            self.Q[np.ix_(H, H)] ^= 1
            self.R0[H] ^= 1
            self.R1[H] ^= self.R0[H] ^ u1
        else:
            if len(H) == 0:
                return
            else:
                l = H[0]
                for h in H[1:]:
                    self.ReindexSubtColumn(h, l)
                self.ReindexSwapColumns(self.r - 1, l)
                self.FixFinalBit(u1)

    def SimulateX(self, j):
        self.b[j] ^= 1

    def SimulateZ(self, j):
        self.R1[: self.r] ^= self.A[j, : self.r]

    def SimulateY(self, j):
        self.SimulateZ(j)
        self.SimulateX(j)

    def SimulateH(self, j):
        c = self.principate(j)
        self.Q[self.r, : self.r] = self.A[j, : self.r]
        self.Q[: self.r, self.r] = self.A[j, : self.r]
        self.R0[self.r] = 0
        self.R1[self.r] = self.b[j]
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
        self.Q[np.ix_(H, H)] ^= 1
        self.R1[H] ^= self.R0[H] ^ self.b[j]
        self.R0[H] ^= 1

    def SimulateSdg(self, j):
        H = [h for h in range(self.r) if self.A[j, h] == 1]
        self.Q[np.ix_(H, H)] ^= 1
        self.R0[H] ^= 1
        self.R1[H] ^= self.R0[H] ^ self.b[j]

    def SimulateCZ(self, j, k):
        assert j != k
        H_j = [h for h in range(self.r) if self.A[j, h] == 1]
        H_k = [h for h in range(self.r) if self.A[k, h] == 1]
        self.Q[np.ix_(H_j, H_k)] ^= 1
        self.Q[np.ix_(H_k, H_j)] ^= 1
        self.R1[: self.r] ^= self.A[j, : self.r] & self.A[k, : self.r]
        self.R1[H_j] ^= self.b[k]
        self.R1[H_k] ^= self.b[j]

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
            if self.R0[c] == 0:
                return self.R1[c]
            else:
                beta = self.toss_coin(coin)
                self.R0[c] = 0
                self.R1[c] = beta
                return beta
        self.Q[self.r, : self.r] = 0
        self.Q[: self.r, self.r] = 0
        self.R1[: self.r] ^= beta & self.A[j, : self.r]
        self.R1[self.r] = beta
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
            if self.R0[c] == 1:
                return self.R1[c]
            else:
                beta = self.toss_coin(coin)
                self.R0[c] = 1
                self.R1[c] = beta
                return beta
        H = [h for h in range(self.r) if self.A[j, h] == 1]
        self.Q[np.ix_(H, H)] ^= 1
        self.R0[H] ^= 1
        self.R1[H] ^= self.R0[H] ^ self.b[j] ^ beta
        self.Q[self.r, : self.r] = 0
        self.Q[: self.r, self.r] = 0
        self.R0[self.r] = 1
        self.R1[self.r] = beta
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
