#ifndef ONE_PHASE_MODEL_H
#define ONE_PHASE_MODEL_H

#include "inmost.h"

using namespace INMOST;

// Model
class OnePhaseModel {
    Mesh* mesh;
    Automatizator aut;
    Residual residual;
    Sparse::Vector update;
    BlockEntry entry_pressure;
    Tag tag_pressure;
    MarkerType boundary_marker;

    double lambda(const Cell &c, const rMatrix& normal);
    double gamma_b(const Cell &c, const Face &f);
    double boundary_conductivity(const Cell &c, const Face &f);
    double internal_conductivity(const Face &f);

public:
    explicit OnePhaseModel(Mesh* init_mesh);

    void setup();
    void fill_boundary_conditions(const std::string& tag_bc_name);
    void fill_permeability(const std::string& tag_K_name);
    void assemble_system();
    void update_solution();

    Tag get_tag(const std::string& name);
    Residual* get_residual();
    Sparse::Vector* get_update();
    void save_result(const std::string &filename);

};




#endif //ONE_PHASE_MODEL_H
