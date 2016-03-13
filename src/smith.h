#pragma once

#include "matrix.h"
#include "types.h"

template <typename T>
void smith_reduce_p(const mod_t p, Matrix<T>& f, MatrixRefList<T>& to_X,
                    MatrixRefList<T>& from_X, MatrixRefList<T>& to_Y,
                    MatrixRefList<T>& from_Y);

#include "smith_impl.h"
