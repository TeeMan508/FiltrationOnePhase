#include "solver.h"


CustomSolver::CustomSolver(Residual &resid, Sparse::Vector &upd) {
    this->residual = resid;
    this->update = upd;
}




void CustomSolver::setup() {
    this->solver.SetParameterReal("absolute_tolerance",1.0e-17);
    this->solver.SetParameterReal("relative_tolerance",1.0e-11);
    this->solver.SetParameterEnum("verbosity", 1);
}

bool CustomSolver::solve() {
    bool success = false;
    const bool solve_output = true;

    if (this->solver.Solve(this->residual.GetResidual(), this->update))
    {
        success = true;
        if (solve_output)
            std::cout << "Solve system iterations " << this->solver.Iterations()
        << " reason: " << this->solver.GetReason() << std::endl;
    }
    else std::cout << "System solution failed" << std::endl;
    return success;
}
