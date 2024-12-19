#ifndef DARCYONEPHASE_H
#define DARCYONEPHASE_H

#include "inmost.h"
#include "../utils/params/params.h"

using namespace INMOST;


class OnePhaseModel {
    Mesh mesh;
    Parameters parameters;
    static void calculate_stiffnes_tensor(const double &E, const double &nu, rMatrix &C);
    void assemble_internal_part(BlockEntry& uvw, Residual& residual);
    void assemble_boundary_part(BlockEntry& uvw, Residual& residual);
    void setup_tags();
    void setup_boundary_conditions();
    void setup_initial_conditions();

public:
    OnePhaseModel() = default;
    Tag get_tag(const std::string& name);
    void init_parameters(const Parameters &params);
    void init_mesh(const Mesh &mesh);
    void setup_model();
    void extract_boundary_conditions(Face face, rMatrix& alpha, rMatrix& beta, rMatrix& gamma);
    void assemble_system(BlockEntry& uvw, Residual& residual);
    void update_stress_and_strain_tags();
    bool solve_elasticity_model();
    void save_result(const std::string &filename);

};




#endif //DARCYONEPHASE_H
