#include <iostream>
#include <cassert>

#include "matrix.h"

int main()
{
  MatrixQ A(2, 3);

  A(0, 0) = 1;
  A(0, 1) = 0;
  A(0, 2) = 1;
  A(1, 0) = 0;
  A(1, 1) = 1;
  A(1, 2) = 1;

  std::cout << A;

  MatrixQ B(3, 2);

  B(0, 0) = 1;
  B(0, 1) = 0;
  B(1, 0) = 0;
  B(1, 1) = 1;
  B(2, 0) = 1;
  B(2, 1) = 1;

  std::cout << B;

  MatrixQ C = A * B;

  MatrixQ C_ref(2, 2);

  C_ref(0, 0) = 2;
  C_ref(0, 1) = 1;
  C_ref(1, 0) = 1;
  C_ref(1, 1) = 2;

  std::cout << C;

  MatrixQ Id = MatrixQ::identity(2);

  std::cout << Id;

  assert(Id * C * Id == C_ref);

  MatrixQ D(3, 3);
  D(1, 1, 2, 2) = MatrixQ::identity(2);
  std::cout << "D: " << D;

  MatrixQ E = D(0, 0, 2, 2);
  std::cout << "E: " << E;

  MatrixQ F(3, 3);
  F(1, 1, 2, 2) = D(0, 0, 2, 2);
  std::cout << "F: " << F;

  MatrixQ G = MatrixQ::identity(6);
  G(0, 0, 2, 2) = G(0, 2, 2, 2);
  std::cout << "G: " << G << std::flush;
}
