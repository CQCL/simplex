#pragma once

#include <iostream>
#include <list>
#include <vector>

class A_matrix {
public:
  A_matrix(unsigned n);

  int entry(unsigned j, unsigned h) const;

  // XOR column k into column h
  void add_col(unsigned h, unsigned k);

  // XOR row k into row j
  void add_row(unsigned j, unsigned k);

  // Number of elements in row j containing 1
  unsigned row_weight(unsigned j) const;

  // Number of elements in column h containing 1
  unsigned col_weight(unsigned h) const;

  // Swap column h with column r-1
  void swap_col(unsigned h);

  // Zero row j and append new column e_j
  void zero_append_basis_col(unsigned j);

  // List of column indices h s.t. A[j,h] = 1
  std::list<unsigned> cols_where_one(unsigned j) const;

  // List of column indices h s.t. A[j,h] = A[j,k] = 1
  std::list<unsigned> cols_where_one(unsigned j, unsigned k) const;

  void drop_final_col();

  friend std::ostream& operator<<(std::ostream& os, const A_matrix& A);

private:
  unsigned n;
  unsigned r;
  std::vector<std::vector<int>> data;
};
