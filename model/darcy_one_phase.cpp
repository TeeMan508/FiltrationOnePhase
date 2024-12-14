#include "darcy_one_phase.h"


void OnePhaseModel::init_parameters(const Parameters &params) {
    this->parameters = params;
}

void OnePhaseModel::init_mesh(const Mesh &mesh) {
    this->mesh = mesh;
}

void OnePhaseModel::calculate_stiffnes_tensor(const double &E, const double &nu, rMatrix &C) {
    double lambda = E * nu / (1 + nu) / (1 - 2*nu), mu = E / 2 / (1 + nu);
    C.Zero();

    C(0, 0) = C(1, 1) = C(2, 2) = lambda + 2 * mu;
    C(0, 1) = C(1, 0) = lambda;
    C(0, 2) = C(2, 0) = lambda;
    C(1, 2) = C(2, 1) = lambda;
    C(3, 3) = C(4, 4) = C(5, 5) = mu;

}

void OnePhaseModel::setup_tags() {
    TagRealArray tag_displacement = this->mesh.CreateTag(
        "Displacement", DATA_REAL, NODE, NONE, 3);

    // for (auto &node : this->mesh.BeginNode()) {
    // TODO: Question about iteration
    //

    TagRealArray tag_boundary_conditions = this->mesh.CreateTag(
        "Boundary Conditions", DATA_REAL, FACE, FACE, 3);

    TagRealArray tag_stiffness_tensor = this->mesh.CreateTag(
        "Stiffness Tensor", DATA_REAL, CELL, NONE, 36);
}

void OnePhaseModel::setup_boundary_conditions() {
    rMatrix normal(3,1);
    TagRealArray tag_boundary_conditions = this->mesh.GetTag("Boundary Conditions");

    for (Mesh::iteratorFace face = this->mesh.BeginFace(); face != this->mesh.EndFace(); ++face) {
        face->UnitNormal(normal.data());

        if (fabs(normal[0] + 1.) < 1e-5) // Left Boundary
            tag_boundary_conditions(*face, 7, 1) = rMatrix::Make(
                7, 1, 1., 1., 0., 0., 0., 0., 0.);

        if (fabs(normal[0] - 1.) < 1e-5) // Right Boundary
            tag_boundary_conditions(*face, 7, 1) = rMatrix::Make(
                7, 1, 1., 1., 0., 0., 0., 0., 0.);
    }

}

void OnePhaseModel::setup_initial_conditions() {
    TagRealArray tag_displacement = this->mesh.GetTag("Displacement");
    TagRealArray tag_stiffness_tensor = this->mesh.GetTag("Stiffness Tensor");

    for (Mesh::iteratorNode node = this->mesh.BeginNode(); node != this->mesh.EndNode(); ++node)
        tag_displacement(*node, 3, 1).Zero();

    // Lame Constants
    constexpr double E = 25., nu = 0.3;
    rMatrix C(6, 6);
    calculate_stiffnes_tensor(E, nu, C);

    for (Mesh::iteratorCell cell = this->mesh.BeginCell(); cell != this->mesh.EndCell(); ++cell)
        tag_stiffness_tensor(*cell, 6, 6) = C;
}

void OnePhaseModel::extract_boundary_conditions(Face face, rMatrix &alpha, rMatrix &beta, rMatrix &gamma) {
    TagRealArray tag_boundary_conditions = this->mesh.GetTag("Boundary Conditions");
    const MatrixUnit<double> I(3);
    double alpha_perp = 0, beta_perp = 1;
    double alpha_parallel = 0, beta_parallel = 1;
    rMatrix normal(3,1);
    face->UnitNormal(normal.data());
    gamma.Zero();

    if (tag_boundary_conditions.isValid() && face.HaveData(tag_boundary_conditions))
    {
        Storage::real_array boundary_conditions = tag_boundary_conditions[face];
        alpha_perp = boundary_conditions[0];
        alpha_parallel = boundary_conditions[1];
        beta_perp = boundary_conditions[2];
        beta_parallel = boundary_conditions[3];
        gamma[0] = boundary_conditions[4];
        gamma[1] = boundary_conditions[5];
        gamma[2] = boundary_conditions[6];
    }
    alpha = (alpha_perp - alpha_parallel) * normal * normal.Transpose() + alpha_parallel * I;
    beta = (beta_perp - beta_parallel) * normal * normal.Transpose() + beta_parallel * I;
}

void OnePhaseModel::assemble_system() {

}


