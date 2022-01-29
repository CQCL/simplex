#include "Q_matrix.hpp"
#include <iostream>
#include <list>
#include <memory>
#include <vector>

struct Q_matrix::impl {
  impl(unsigned n) : n(n), r(0), data(n+1, std::vector<int>(n+1, 0)) {}

  /* Data */

  unsigned n;
  unsigned r;
  std::vector<std::vector<int>> data;

  /* Methods */

  int entry(unsigned h1, unsigned h2) const { return data[h1][h2]; }

  void add_rowcol(unsigned h, unsigned k) {
    for (unsigned j = 0; j < r; j++) {
      data[j][h] ^= data[j][k];
    }
    for (unsigned j = 0; j < r; j++) {
      data[h][j] ^= data[k][j];
    }
  }

  void swap_rowcol(unsigned h) {
    const unsigned r1 = r - 1;
    for (unsigned j = 0; j < r; j++) {
      int t = data[j][h];
      data[j][h] = data[j][r1];
      data[j][r1] = t;
    }
    for (unsigned j = 0; j < r; j++) {
      int t = data[h][j];
      data[h][j] = data[r1][j];
      data[r1][j] = t;
    }
  }

  std::list<unsigned> rows_with_terminal_1() const {
    std::list<unsigned> H;
    unsigned r1 = r - 1;
    for (unsigned h = 0; h < r1; h++) {
      if (data[h][r1]) H.push_back(h);
    }
    return H;
  }

  void flip_submatrix(const std::list<unsigned>& H) {
    for (unsigned h1 : H) {
      for (unsigned h2 : H) {
        data[h1][h2] ^= 1;
      }
    }
  }

  void flip_submatrix(
    const std::list<unsigned>& H1, const std::list<unsigned>& H2)
  {
    for (unsigned h1 : H1) {
      for (unsigned h2 : H2) {
        data[h1][h2] ^= 1;
        data[h2][h1] ^= 1;
      }
    }
  }

  void append_rowcol(const std::list<unsigned>& H) {
    for (unsigned h = 0; h < r; h++) {
      data[h][r] = data[r][h] = 0;
    }
    for (unsigned h : H) {
      data[h][r] = data[r][h] = 1;
    }
    r++;
  }

  bool rowcol_is_zero(unsigned h) const {
    for (unsigned j = 0; j < r; j++) {
      if (j != h && data[h][j]) {
        return false;
      }
    }
    return true;
  }

  void drop_final_rowcol() { r--; }
};

Q_matrix::Q_matrix(unsigned n) : pImpl(std::make_unique<impl>(n)) {}
Q_matrix::~Q_matrix() = default;
Q_matrix::Q_matrix(const Q_matrix& other)
  : pImpl(std::make_unique<impl>(*other.pImpl)) {}
Q_matrix::Q_matrix(Q_matrix&& other) = default;
Q_matrix& Q_matrix::operator=(const Q_matrix& other) {
  return *this = Q_matrix(other);
}
Q_matrix& Q_matrix::operator=(Q_matrix&& other) = default;

int Q_matrix::entry(unsigned h1, unsigned h2) const {
  return pImpl->entry(h1, h2);
}
void Q_matrix::add_rowcol(unsigned h, unsigned k) { pImpl->add_rowcol(h, k); }
void Q_matrix::swap_rowcol(unsigned h) { pImpl->swap_rowcol(h); }
std::list<unsigned> Q_matrix::rows_with_terminal_1() const {
  return pImpl->rows_with_terminal_1();
}
void Q_matrix::flip_submatrix(const std::list<unsigned>& H) {
  pImpl->flip_submatrix(H);
}
void Q_matrix::flip_submatrix(
  const std::list<unsigned>& H1, const std::list<unsigned>& H2) {
  pImpl->flip_submatrix(H1, H2);
}
void Q_matrix::append_rowcol(const std::list<unsigned>& H) {
  pImpl->append_rowcol(H);
}
bool Q_matrix::rowcol_is_zero(unsigned h) const {
  return pImpl->rowcol_is_zero(h);
}
void Q_matrix::drop_final_rowcol() {pImpl->drop_final_rowcol(); }
unsigned Q_matrix::r() const { return pImpl->r; }

std::ostream& operator<<(std::ostream& os, const Q_matrix& Q) {
  for (unsigned j = 0; j < Q.pImpl->r; j++) {
    os << "[ ";
    for (unsigned h = 0; h < Q.pImpl->r; h++) {
      os << Q.pImpl->data[j][h] << " ";
    }
    os << "]" << std::endl;
  }
  return os;
}
