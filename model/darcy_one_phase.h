#ifndef DARCYONEPHASE_H
#define DARCYONEPHASE_H

#include "inmost.h"
#include "../utils/params/params.h"
#include "../utils/calculations/calculations.h"

using namespace INMOST;


class OnePhaseModel : AbstractSubModel {
    Mesh mesh;
    Parameters parameters;
    static void calculate_stiffnes_tensor(const double &E, const double &nu, rMatrix &C);

public:
    OnePhaseModel() = default;

    void init_parameters(const Parameters &params);
    void init_mesh(const Mesh &mesh);
    void setup_tags();
    void setup_boundary_conditions();
    void setup_initial_conditions();
    void extract_boundary_conditions(
        Face face, rMatrix& alpha, rMatrix& beta, rMatrix& gamma);
    void assemble_system();
};




#endif //DARCYONEPHASE_H
