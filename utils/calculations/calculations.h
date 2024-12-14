#ifndef CALCULATIONS_H
#define CALCULATIONS_H

#include "inmost.h"
using namespace INMOST;


template<typename MatrixType>
void gradient(Cell c, BlockEntry& UVW, Matrix<MatrixType>& G);

template<typename MatrixType>
void gradient_to_strain(const Matrix<MatrixType>& G, Matrix<MatrixType>& strain);

template<typename MatrixType>
void stress_to_traction(const Matrix<MatrixType>& stress, const rMatrix& n, Matrix<MatrixType>& trac)


#endif //CALCULATIONS_H
