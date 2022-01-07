#include "simplex.hpp"
#include "bimap.hpp"

#include <algorithm>
#include <iostream>
#include <list>
#include <memory>
#include <optional>
#include <random>
#include <vector>

/* Implementation */

class A_matrix {
public:
  A_matrix(unsigned n) : n(n), r(0), data(n, std::vector<int>(n+1, 0)) {}

  // Reference to one row
  std::vector<int>& row(unsigned j) {
    return data[j];
  }

  // Const reference to one row
  const std::vector<int>& constrow(unsigned j) const {
    return data[j];
  }

  // XOR column k into column h
  void add_col(unsigned h, unsigned k) {
    for (unsigned j = 0; j < n; j++) {
      data[j][h] ^= data[j][k];
    }
  }

  // XOR row k into row j
  void add_row(unsigned j, unsigned k) {
    for (unsigned h = 0; h < r; h++) {
      data[j][h] ^= data[k][h];
    }
  }

  // Number of elements in row j containing 1
  unsigned row_weight(unsigned j) {
    unsigned c = 0;
    for (unsigned h = 0; h < r; h++) {
      if (data[j][h]) c++;
    }
    return c;
  }

  // Number of elements in column h containing 1
  unsigned col_weight(unsigned h) {
    unsigned c = 0;
    for (unsigned j = 0; j < n; j++) {
      if (data[j][h]) c++;
    }
    return c;
  }

  // Swap (distinct) columns
  void swap_cols(unsigned h1, unsigned h2) {
    for (unsigned j = 0; j < n; j++) {
      std::vector<int>& A_j = row(j);
      std::iter_swap(A_j.begin() + h1, A_j.begin() + h2);
    }
  }

  // Set a row to zero
  void zero_row(unsigned j) {
    std::vector<int>& A_j = row(j);
    for (unsigned h = 0; h < r; h++) {
      A_j[h] = 0;
    }
  }

  // Write e_j to column r
  void append_basis_col(unsigned j) {
    for (unsigned k = 0; k < n; k++) {
      data[k][r] = 0;
    }
    data[j][r] = 1;
  }

  // List of column indices h s.t. A[j,h] = 1
  std::list<unsigned> cols_where_one(unsigned j) const {
    const std::vector<int>& A_j = constrow(j);
    std::list<unsigned> H;
    for (unsigned h = 0; h < r; h++) {
      if (A_j[h]) {
        H.push_back(h);
      }
    }
    return H;
  }

  void increment_r() { r++; }
  void decrement_r() { r--; }

  friend std::ostream& operator<<(std::ostream& os, const A_matrix& A);

private:
  unsigned n;
  unsigned r;
  std::vector<std::vector<int>> data;
};

std::ostream& operator<<(std::ostream& os, const A_matrix& A) {
  for (unsigned j = 0; j < A.n; j++) {
    for (unsigned h = 0; h < A.r; h++) {
      std::cout << A.data[j][h];
    }
    std::cout << std::endl;
  }
  return os;
}

// A symmetric 0,1 off-diagonal matrix
class Q_matrix {
public:
  Q_matrix(unsigned n) : n(n), r(0), data(n+1, std::vector<int>(n+1, 0)) {}

  int entry(unsigned h1, unsigned h2) { return data[h1][h2]; }

  // XOR row/column k into row/column h
  void add_rowcol(unsigned h, unsigned k) {
    for (unsigned j = 0; j < r; j++) {
      data[j][h] ^= data[j][k];
    }
    for (unsigned j = 0; j < r; j++) {
      data[h][j] ^= data[k][j];
    }
  }

  // Swap (distinct) rows and columns
  void swap_rowcol(unsigned h1, unsigned h2) {
    for (unsigned j = 0; j < r; j++) {
      int t = data[j][h1];
      data[j][h1] = data[j][h2];
      data[j][h2] = t;
    }
    for (unsigned j = 0; j < r; j++) {
      int t = data[h1][j];
      data[h1][j] = data[h2][j];
      data[h2][j] = t;
    }
  }

  // List of h s.t. Q[h][r-1] = 1
  std::list<unsigned> rows_with_terminal_1() const {
    std::list<unsigned> H;
    unsigned r1 = r - 1;
    for (unsigned h = 0; h < r1; h++) {
      if (data[h][r1]) H.push_back(h);
    }
    return H;
  }

  // Flip the (h1, h2) entry for all h1, h2 in H
  void flip_submatrix(const std::list<unsigned>& H) {
    for (auto h1 : H) {
      for (auto h2 : H) {
        data[h1][h2] ^= 1;
      }
    }
  }

  // Flip the (h1, h2) entry for all h1 in H1, h2 in H2 (preserving symmetry)
  void flip_submatrix(
    const std::list<unsigned>& H1, const std::list<unsigned>& H2)
  {
    for (auto h1 : H1) {
      for (auto h2 : H2) {
        data[h1][h2] ^= 1;
        data[h2][h1] ^= 1;
      }
    }
  }

  // Append the given row/column, without incrementing r
  void append_rowcol(const std::vector<int>& a) {
    for (unsigned h = 0; h < r; h++) {
      data[h][r] = data[r][h] = a[h];
    }
  }

  // Whether a given row/column is all-zero
  bool rowcol_is_zero(unsigned h) {
    for (unsigned j = 0; j < r; j++) {
      if (j != h && data[h][j]) {
        return false;
      }
    }
    return true;
  }

  // Append a zero row/column, without incrementing r
  void append_zero_rowcol() {
    for (unsigned h = 0; h < r; h++) {
      data[h][r] = data[r][h] = 0;
    }
  }

  void increment_r() { r++; }
  void decrement_r() { r--; }

  friend std::ostream& operator<<(std::ostream& os, const Q_matrix& Q);

private:
  unsigned n;
  unsigned r;
  std::vector<std::vector<int>> data;
};

std::ostream& operator<<(std::ostream& os, const Q_matrix& Q) {
  for (unsigned j = 0; j < Q.r; j++) {
    for (unsigned h = 0; h < Q.r; h++) {
      std::cout << Q.data[j][h];
    }
    std::cout << std::endl;
  }
  return os;
}

class RBG {
public:
  RBG(int seed = 0) : gen(seed), distrib(0, 1) {}
  int get() { return distrib(gen); }
private:
  std::mt19937 gen;
  std::uniform_int_distribution<> distrib;
};

struct Simplex::impl {
  impl(unsigned n, int seed = 0)
    : n(n), r(0), A(n), b(n, 0), Q(n), R0(n+1, 0), R1(n+1, 0), p(),
    deterministic(true), rbg(seed) {}

  /* Data */

  unsigned n;
  unsigned r;
  A_matrix A;
  std::vector<int> b;
  Q_matrix Q;
  std::vector<int> R0;
  std::vector<int> R1;
  Bimap p;
  bool deterministic;
  RBG rbg;

  /* Methods */

  void ReindexSubtColumn(unsigned k, unsigned c) {
    if (k == c) return;
    A.add_col(k, c);
    R1[k] ^= Q.entry(c, k);
    Q.add_rowcol(k, c);
  }

  void MakePrincipal(unsigned c, unsigned j) {
    std::vector<int>& A_j = A.row(j);
    if (A_j[c]) {
      for (unsigned k = 0; k < r; k++) {
        if ((k != c) && A_j[k]) {
          ReindexSubtColumn(k, c);
        }
      }
      p.make_match(c, j);
    }
  }

  void ReselectPrincipalRow(
    unsigned c, std::optional<unsigned> j = std::nullopt)
  {
    std::optional<unsigned> n0;
    unsigned j0;
    for (unsigned j1 = 0; j1 < n; j1++) {
      if (!j || j1 != *j) {
        std::vector<int>& A_j1 = A.row(j1);
        if (A_j1[c]) {
          unsigned n1 = A.row_weight(j1);
          if (!n0 || n1 < *n0) {
            j0 = j1;
            n0 = n1;
          }
        }
      }
    }
    if (n0) {
      MakePrincipal(c, j0);
    }
  }

  std::optional<unsigned> principate(unsigned j) {
    std::optional<unsigned> c = p.inv_at(j);
    if (c) {
      ReselectPrincipalRow(*c, j);
      if (j != *p.fwd_at(*c)) {
        c = std::nullopt;
      }
    }
    return c;
  }

  void ReindexSwapColumns(unsigned k, unsigned c) {
    if (k == c) return;
    A.swap_cols(k, c);
    std::iter_swap(R0.begin() + k, R0.begin() + c);
    std::iter_swap(R1.begin() + k, R1.begin() + c);
    Q.swap_rowcol(k, c);
    p.swap_fwd(k, c);
  }

  void increment_r() {
    A.increment_r();
    Q.increment_r();
    r++;
  }

  void decrement_r() {
    p.fwd_erase(r - 1);
    A.decrement_r();
    Q.decrement_r();
    r--;
  }

  void FixFinalBit(int z) {
    for (unsigned j = 0; j < n; j++) {
      b[j] ^= z & A.row(j)[r -1];
    }
    decrement_r();
    for (unsigned h = 0; h < r; h++) {
      R1[h] ^= z & Q.entry(h, r);
    }
  }

  void new_principal_column(
    unsigned j, std::optional<unsigned> c = std::nullopt)
  {
    A.zero_row(j);
    A.append_basis_col(j);
    p.make_match(r, j);
    b[j] = 0;
    increment_r();
    if (c) {
      ZeroColumnElim(*c);
    }
  }

  void ZeroColumnElim(unsigned c) {
    ReindexSwapColumns(c, r - 1);
    std::list<unsigned> H = Q.rows_with_terminal_1();
    int u0 = R0[r - 1];
    int u1 = R1[r - 1];
    decrement_r();
    if (u0) {
      Q.flip_submatrix(H);
      for (auto h : H) {
        R0[h] ^= 1;
        R1[h] ^= R0[h] ^ u1;
      }
    } else if (!H.empty()) {
      unsigned l = H.front();
      H.pop_front();
      for (auto h : H) {
        ReindexSubtColumn(h, l);
      }
      ReindexSwapColumns(r - 1, l);
      FixFinalBit(u1);
    }
  }

  void SimulateX(unsigned j) { b[j] ^= 1; }

  void SimulateY(unsigned j) { SimulateZ(j); SimulateX(j); }

  void SimulateZ(unsigned j) {
    const std::vector<int>& A_j = A.row(j);
    for (unsigned h = 0; h < r; h++) {
      R1[h] ^= A_j[h];
    }
  }

  void SimulateH(unsigned j) {
    std::optional<unsigned> c = principate(j);
    std::vector<int>& A_j = A.row(j);
    Q.append_rowcol(A_j);
    R0[r] = 0;
    R1[r] = b[j];
    new_principal_column(j, c);
  }

  void SimulateS(unsigned j) {
    const std::list<unsigned> H = A.cols_where_one(j);
    Q.flip_submatrix(H);
    int z = b[j];
    for (auto h : H) {
      R1[h] ^= R0[h] ^ z;
      R0[h] ^= 1;
    }
  }

  void SimulateSdg(unsigned j) {
    const std::list<unsigned> H = A.cols_where_one(j);
    Q.flip_submatrix(H);
    const int z = b[j];
    for (auto h : H) {
      R0[h] ^= 1;
      R1[h] ^= R0[h] ^ z;
    }
  }
  void SimulateCX(unsigned j, unsigned k) {
    A.add_row(k, j);
    b[k] ^= b[j];
    std::optional<unsigned> c = p.inv_at(k);
    if (c) {
      ReselectPrincipalRow(*c);
    }
  }

  void SimulateCZ(unsigned j, unsigned k) {
    const std::list<unsigned> H_j = A.cols_where_one(j);
    const std::list<unsigned> H_k = A.cols_where_one(k);
    Q.flip_submatrix(H_j, H_k);
    const std::vector<int>& A_j = A.constrow(j);
    const std::vector<int>& A_k = A.constrow(k);
    for (unsigned h = 0; h < r; h++) {
      R1[h] ^= A_j[h] & A_k[h];
    }
    const int z_j = b[j];
    const int z_k = b[k];
    for (auto h : H_j) {
      R1[h] ^= z_k;
    }
    for (auto h : H_k) {
      R1[h] ^= z_j;
    }
  }

  int toss_coin(std::optional<int> coin) {
    deterministic = false;
    if (coin) {
      return *coin;
    } else {
      return rbg.get();
    }
  }

  int SimulateMeasX(unsigned j, std::optional<int> coin) {
    int beta;
    std::optional<unsigned> c = principate(j);
    if (c && Q.rowcol_is_zero(*c)) {
      if (R0[*c] == 0) {
        return R1[*c];
      } else {
        beta = toss_coin(coin);
        R0[*c] = 0;
        R1[*c] = beta;
        return beta;
      }
    } else {
      beta = toss_coin(coin);
    }
    Q.append_zero_rowcol();
    const std::vector<int>& A_j = A.row(j);
    for (unsigned h = 0; h < r; h++) {
      R1[h] ^= beta & A_j[h];
    }
    R1[r] = beta;
    new_principal_column(j, c);
    return beta;
  }

  int SimulateMeasY(unsigned j, std::optional<int> coin) {
    int beta;
    std::optional<unsigned> c = principate(j);
    if (c && Q.rowcol_is_zero(*c)) {
      if (R0[*c] == 1) {
        return R1[*c];
      } else {
        beta = toss_coin(coin);
        R0[*c] = 1;
        R1[*c] = beta;
        return beta;
      }
    } else {
      beta = toss_coin(coin);
    }
    const std::list<unsigned> H = A.cols_where_one(j);
    Q.flip_submatrix(H);
    const int z = b[j] ^ beta;
    for (auto h : H) {
      R0[h] ^= 1;
      R1[h] ^= R0[h] ^ z;
    }
    Q.append_zero_rowcol();
    R0[r] = 1;
    R1[r] = beta;
    new_principal_column(j, c);
    return beta;
  }

  int SimulateMeasZ(unsigned j, std::optional<int> coin) {
    if (A.row_weight(j) == 0) {
      return b[j];
    } else {
      int beta = toss_coin(coin);
      const std::list<unsigned> H = A.cols_where_one(j);
      unsigned k;
      unsigned m = n + 1;
      for (auto h : H) {
        unsigned c = A.col_weight(h);
        if (c < m) {
          k = h;
          m = c;
        }
      }
      ReindexSwapColumns(k, r - 1);
      MakePrincipal(r - 1, j);
      FixFinalBit(beta ^ b[j]);
      return beta;
    }
  }

  bool is_deterministic() const { return deterministic; }
};

/* Public interface */

Simplex::Simplex(unsigned n, int seed)
  : pImpl(std::make_unique<impl>(n, seed)) {}

Simplex::~Simplex() = default;
Simplex::Simplex(const Simplex& other)
  : pImpl(std::make_unique<impl>(*other.pImpl)) {}
Simplex::Simplex(Simplex&& other) = default;
Simplex& Simplex::operator=(const Simplex& other) {
  return *this = Simplex(other);
}
Simplex& Simplex::operator=(Simplex&& other) = default;

unsigned Simplex::n() const { return pImpl->n; }
void Simplex::X(unsigned j) { pImpl->SimulateX(j); }
void Simplex::Y(unsigned j) { pImpl->SimulateY(j); }
void Simplex::Z(unsigned j) { pImpl->SimulateZ(j); }
void Simplex::H(unsigned j) { pImpl->SimulateH(j); }
void Simplex::S(unsigned j) { pImpl->SimulateS(j); }
void Simplex::Sdg(unsigned j) { pImpl->SimulateSdg(j); }
void Simplex::CX(unsigned j, unsigned k) { pImpl->SimulateCX(j, k); }
void Simplex::CZ(unsigned j, unsigned k) { pImpl->SimulateCZ(j, k); }
int Simplex::MeasX(unsigned j, std::optional<int> coin) {
  return pImpl->SimulateMeasX(j, coin);
}
int Simplex::MeasY(unsigned j, std::optional<int> coin) {
  return pImpl->SimulateMeasY(j, coin);
}
int Simplex::MeasZ(unsigned j, std::optional<int> coin) {
  return pImpl->SimulateMeasZ(j, coin);
}
bool Simplex::is_deterministic() const { return pImpl->is_deterministic(); }

std::ostream& operator<<(std::ostream& os, const Simplex& S) {
  os << "n: " << S.n() << std::endl;
  os << "A:" << std::endl << S.pImpl->A;
  os << "b:" << std::endl;
  for (unsigned j = 0; j < S.n(); j++) {
    os << S.pImpl->b[j];
  }
  os << std::endl;
  os << "Q:" << std::endl << S.pImpl->Q;
  os << "R0:" << std::endl;
  for (unsigned h = 0; h < S.pImpl->r; h++) {
    os << S.pImpl->R0[h];
  }
  os << std::endl;
  os << "R1:" << std::endl;
  for (unsigned h = 0; h < S.pImpl->r; h++) {
    os << S.pImpl->R1[h];
  }
  os << std::endl;
  os << "p:" << std::endl << S.pImpl->p;
  return os;
}
