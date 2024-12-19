#include <iostream>
#include "utils/params/params.h"
#include "model/darcy_one_phase.h"

using namespace INMOST;


int main(int argc, char** argv) {
    Mesh::Initialize(&argc, &argv);
    Solver::Initialize(&argc, &argv);

    Parameters params;
    params.load("../data/example.txt");

    if (argc < 1) {
        std::cout << "Usage: " << argv[0] << " mesh_file " << std::endl;
        return 1;
    };

    Mesh m;
    try { m.Load(argv[1]); }
    catch (...) { std::cout << "Cannot read mesh file: " << argv[1] << std::endl; return 1;}

    OnePhaseModel model;
    model.init_mesh(m);
    model.init_parameters(params);
    model.setup_model();

    model.solve_elasticity_model();
    model.update_stress_and_strain_tags();

    model.save_result("out.vtu");
    Solver::Finalize();
    Mesh::Finalize();
    return 0;
}
// Вопросы
// по iterable objects
// про передачу ссылки на mesh
// про инициализацию Solver для обертки
// про sparse при создании тега
// про parameters
//
