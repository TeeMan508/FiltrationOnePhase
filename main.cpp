#include <iostream>
#include "utils/params/params.h"
#include "model/one_phase_model.h"

using namespace INMOST;


int main(int argc, char** argv) {
    if (argc < 1) {
        std::cout << "Usage: " << argv[0] << " mesh_file " << std::endl;
        return 1;
    }
    Mesh *m = new Mesh;
    try {
        m->Load(argv[1]);
    }
    catch (...) {
        std::cout << "Cannot read mesh file: " << argv[1] << std::endl;
        return 1;
    }
    // rMatrix buf(3,1);
    // TagRealArray tag_bc = m.GetTag("BOUNDARY_CONDITION");
    // for (Mesh::iteratorFace face = m.BeginFace(); face != m.EndFace(); ++face) {
    //     if (!face->Boundary()) continue;
    //     buf = rMatrix::FromVector(tag_bc[*face].data(), tag_bc[*face].size());
    //     std::cout << buf[0] << " " << buf[1] << " " << buf[2] << std::endl;
    // }

    OnePhaseModel model(m);
    model.setup();
    model.assemble_system();

    Residual* residual = model.get_residual() ;
    Sparse::Vector* update = model.get_update();

    Solver solver(Solver::INNER_ILU2);
    solver.SetParameterEnum("verbosity", 1);
    solver.SetMatrix(residual->GetJacobian());

    if (solver.Solve(residual->GetResidual(), *update)) {
        model.update_solution();
    }
    else std::cout << "Solve failed!" << std::endl;

    model.save_result("../result.vtk");

    return 0;


}

