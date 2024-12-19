#include "solver.h"


CustomSolver::CustomSolver(Residual &resid, Sparse::Vector &upd) {
    // std::cout << "11111" << std::endl;

    // this->solver = Solver(Solver::INNER_ILU2);
    // std::cout << solver.GetParameter("absolute_tolerance");

    this->residual = resid;
    this->solver.SetMatrix(residual.GetJacobian());
    this->update = upd;
}





void CustomSolver::setup() {
    this->solver.SetParameterReal("absolute_tolerance",1.0e-17);
    this->solver.SetParameterReal("relative_tolerance",1.0e-11);
    this->solver.SetParameterEnum("verbosity", 1);
    // std::cout << "setuped" << std::endl;

}

bool CustomSolver::solve() {
    bool success = false;
    const bool solve_output = true;
    std::cout << "11111" << std::endl;
    Solver s = this->solver;
    std::cout << "11111" << std::endl;
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
