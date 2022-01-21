#include "Q_matrix.hpp"
#include <iostream>
#include <list>
#include <vector>

Q_matrix::Q_matrix(unsigned n) : n(n), r(0), data(n+1, std::vector<int>(n+1, 0))
{}

int Q_matrix::entry(unsigned h1, unsigned h2) const { return data[h1][h2]; }

// XOR row/column k into row/column h
void Q_matrix::add_rowcol(unsigned h, unsigned k) {
  for (unsigned j = 0; j < r; j++) {
    data[j][h] ^= data[j][k];
  }
  for (unsigned j = 0; j < r; j++) {
    data[h][j] ^= data[k][j];
  }
}

// Swap row/column with row/column r-1
void Q_matrix::swap_rowcol(unsigned h) {
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

// List of h s.t. Q[h][r-1] = 1
std::list<unsigned> Q_matrix::rows_with_terminal_1() const {
  std::list<unsigned> H;
  unsigned r1 = r - 1;
  for (unsigned h = 0; h < r1; h++) {
    if (data[h][r1]) H.push_back(h);
  }
  return H;
}

// Flip the (h1, h2) entry for all h1, h2 in H
void Q_matrix::flip_submatrix(const std::list<unsigned>& H) {
  for (unsigned h1 : H) {
    for (unsigned h2 : H) {
      data[h1][h2] ^= 1;
    }
  }
}

// Flip the (h1, h2) entry for all h1 in H1, h2 in H2 (preserving symmetry)
void Q_matrix::flip_submatrix(
  const std::list<unsigned>& H1, const std::list<unsigned>& H2)
{
  for (unsigned h1 : H1) {
    for (unsigned h2 : H2) {
      data[h1][h2] ^= 1;
      data[h2][h1] ^= 1;
    }
  }
}

// Append new row/column (given as list of indices where 1)
void Q_matrix::append_rowcol(const std::list<unsigned>& H) {
  for (unsigned h = 0; h < r; h++) {
    data[h][r] = data[r][h] = 0;
  }
  for (unsigned h : H) {
    data[h][r] = data[r][h] = 1;
  }
  r++;
}

// Whether a given row/column is all-zero
bool Q_matrix::rowcol_is_zero(unsigned h) const {
  for (unsigned j = 0; j < r; j++) {
    if (j != h && data[h][j]) {
      return false;
    }
  }
  return true;
}

void Q_matrix::drop_final_rowcol() { r--; }

std::ostream& operator<<(std::ostream& os, const Q_matrix& Q) {
  for (unsigned j = 0; j < Q.r; j++) {
    for (unsigned h = 0; h < Q.r; h++) {
      std::cout << Q.data[j][h];
    }
    std::cout << std::endl;
  }
  return os;
}
