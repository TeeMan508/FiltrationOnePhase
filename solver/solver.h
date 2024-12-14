#ifndef SOLVER_H
#define SOLVER_H

#include "inmost.h"
#include "../utils/calculations/calculations.h"

using namespace INMOST;


class CustomSolver {
    Solver solver = Solver(Solver::INNER_ILU2);
    Residual residual;
    Sparse::Vector update;

public:
    CustomSolver(Residual& resid, Sparse::Vector& upd);
    void setup();
    bool solve();
};


#endif //SOLVER_H
