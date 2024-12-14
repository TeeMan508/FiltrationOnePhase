#include "calculations.h"

using namespace INMOST;

template<typename MatrixType>
void gradient(Cell c, BlockEntry& UVW, Matrix<MatrixType>& G)
{
    const MatrixUnit<double> I(3);
    ElementArray<Node> nodes = c.getNodes();
    rMatrix A(9, 9), v(3, 1), xc(3, 1), xf(3, 1), xn(3, 1);
    Matrix<MatrixType> b(9, 1), U(3, 1), dU(3, 1);
    c.Centroid(xc.data());
    U.Zero();
    for (unsigned l = 0; l < nodes.size(); ++l)
        U += UVW.Access<MatrixType>(nodes[l]);
    U /= (double)nodes.size();
    A.Zero();
    b.Zero();
    for (unsigned l = 0; l < nodes.size(); ++l)
    {
        nodes[l].Centroid(xn.data());
        v = xn - xc;
        dU = UVW.Access<MatrixType>(nodes[l]) - U;
        A += I.Kronecker(v * v.Transpose());
        b += dU.Kronecker(v);
    }
    G = A.CholeskySolve(b);
}

template<typename MatrixType>
void gradient_to_strain(const Matrix<MatrixType>& G, Matrix<MatrixType>& strain)
{
    strain[0] = G[0];
    strain[1] = G[4];
    strain[2] = G[8];
    strain[3] = G[5] + G[7];
    strain[4] = G[2] + G[6];
    strain[5] = G[1] + G[3];
}

template<typename MatrixType>
void stress_to_traction(const Matrix<MatrixType>& stress, const rMatrix& n, Matrix<MatrixType>& trac)
{
    trac[0] = n[0] * stress[0] + n[1] * stress[5] + n[2] * stress[4];
    trac[1] = n[0] * stress[5] + n[1] * stress[1] + n[2] * stress[3];
    trac[2] = n[0] * stress[4] + n[1] * stress[3] + n[2] * stress[2];

}

