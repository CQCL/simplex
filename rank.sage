from sage.all import *
import sys

# Accepts a string representation of a 0,1-matrix as list of lists.
# Prints the rank of the matrix over GF(2).
# E.g.:
#   /usr/bin/python rank.sage "[[0,1],[1,1]]"
# prints "2".

A = matrix(GF(2), eval(sys.argv[1]))

print(A.rank())
