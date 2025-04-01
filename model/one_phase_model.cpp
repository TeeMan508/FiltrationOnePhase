#include "one_phase_model.h"


OnePhaseModel::OnePhaseModel(Mesh* init_mesh) {
    mesh = init_mesh;
    tag_pressure = mesh->CreateTag("pressure", DATA_REAL, CELL, NONE, 1);
    boundary_marker = mesh->CreateMarker();

    entry_pressure = BlockEntry();
    entry_pressure.AddTag(tag_pressure);
    entry_pressure.SetElementType(CELL);

    aut = Automatizator("");
    aut.RegisterEntry(entry_pressure);
    aut.EnumerateEntries();

    residual = Residual("", aut.GetFirstIndex(), aut.GetLastIndex());
    update = Sparse::Vector("", aut.GetFirstIndex(), aut.GetLastIndex());

}

Tag OnePhaseModel::get_tag(const std::string &name) {
    return mesh->GetTag(name);
}


void OnePhaseModel::setup() {
    mesh->MarkBoundaryFaces(boundary_marker);

    mesh->CreateTag("permeability", DATA_REAL, CELL, NONE, 9);
    mesh->CreateTag("boundary_conditions", DATA_REAL, FACE, FACE, 3);
    this->fill_boundary_conditions("BOUNDARY_CONDITION");
    this->fill_permeability("PERM");

}


void OnePhaseModel::fill_boundary_conditions(const std::string& tag_bc_name) {
    if (mesh->HaveTag(tag_bc_name)) {
        std::cout << "HAVE!!!";
    }

    const TagRealArray tag_bc_prepared = mesh->GetTag(tag_bc_name);
    const TagRealArray tag_bc = mesh->GetTag("boundary_conditions");

    rMatrix buf(3, 1);

    for (Mesh::iteratorFace face = mesh->BeginFace(); face != mesh->EndFace(); ++face) {
        if (!face->Boundary()) continue;
        buf = rMatrix::FromVector(tag_bc_prepared[*face].data(), tag_bc_prepared[*face].size());
        std::cout << buf[0] << " " << buf[1] << " " << buf[2] << std::endl;
        tag_bc(*face, 3, 1) = buf;

    }
    // for (Mesh::iteratorFace face = mesh->BeginFace(); face != mesh->EndFace(); ++face) {
        // if (face->getFaces(boundary_marker).empty()) continue;
        // std::cout << tag_bc(*face, 3, 1)[0] << " " << tag_bc(*face, 3, 1)[1] << " " << tag_bc(*face, 3, 1)[2] << std::endl;
    // }
}


void OnePhaseModel::fill_permeability(const std::string &tag_K_name) {
    const TagRealArray tag_perm_prepared = mesh->GetTag(tag_K_name);
    const TagRealArray tag_perm = mesh->GetTag("permeability");

    rMatrix buf;

    for (Mesh::iteratorCell c = mesh->BeginCell(); c!= mesh->EndCell(); ++c) {
        buf = rMatrix::FromTensor(tag_perm_prepared[*c].data(), tag_perm_prepared[*c].size());
        tag_perm(*c, 3, 3) = buf;

    }
}

double OnePhaseModel::lambda(const Cell &c, const rMatrix& normal) {
    const TagRealArray tag_permeability = mesh->GetTag("permeability");

    rMatrix const K = tag_permeability(c, 3, 3);
    return normal.DotProduct(K * normal);
}

double OnePhaseModel::gamma_b(const Cell &c, const Face &f) {
    const TagRealArray tag_boundary_cond = mesh->GetTag("boundary_conditions");

    rMatrix normal(3,1), x_face(3,1), x_cell(3,1);
    f.UnitNormal(normal.data());

    f.Centroid(x_face.data());
    c.Centroid(x_cell.data());

    const rMatrix delta = (x_face - x_cell) * 1.;
    const double r = normal.DotProduct(delta);

    Storage::real_array bc = tag_boundary_cond[f];
    const double alpha = bc[0];
    const double beta = bc[1];
    const double gamma_bc = bc[2];

    const double lamda_i = lambda(c, normal);

    return lamda_i * gamma_bc / (alpha * r + beta * lamda_i);

}

double OnePhaseModel::boundary_conductivity(const Cell &c, const Face &f) {
    const TagRealArray tag_boundary_cond = mesh->GetTag("boundary_conditions");

    rMatrix normal(3,1), x_face(3,1), x_cell(3,1);
    f.UnitNormal(normal.data());
    f.Centroid(x_face.data());
    c.Centroid(x_cell.data());

    const rMatrix delta = (x_face - x_cell) * 1.;
    const double r = normal.DotProduct(delta);

    Storage::real_array bc = tag_boundary_cond[f];
    const double alpha = bc[0];
    const double beta = bc[1];

    const double lamda_i = lambda(c, normal);

    return alpha*lamda_i / (alpha*r + beta*lamda_i);
}

double OnePhaseModel::internal_conductivity(const Face &f) {
    TagRealArray tag_permeability = mesh->GetTag("permeability");
    TagRealArray tag_boundary_cond = mesh->GetTag("boundary_conditions");

    rMatrix normal(3,1), normal_2(3, 1), x_face(3,1), x_cell_1(3,1), x_cell_2(3,1);
    const Cell c1 = f->BackCell();
    const Cell c2 = f->FrontCell();

    f.UnitNormal(normal.data());
    f.Centroid(x_face.data());
    c1.Centroid(x_cell_1.data());
    c2.Centroid(x_cell_2.data());

    const rMatrix delta_1 = (x_face - x_cell_1) * 1.;
    const rMatrix delta_2 = (x_face - x_cell_2) * -1.;
    const double r1 = normal.DotProduct(delta_1);
    const double r2 = normal.DotProduct(delta_2);

    const double lambda_1 = lambda(c1, normal);
    const double lambda_2 = lambda(c2, normal * -1);

    return lambda_1 * lambda_2 / (lambda_1*r2 + lambda_2*r1);

}

void OnePhaseModel::assemble_system() {
    TagRealArray tag_permeability = mesh->GetTag("permeability");
    TagRealArray tag_boundary_cond = mesh->GetTag("boundary_conditions");

    for (Mesh::iteratorCell c = mesh->BeginCell(); c!= mesh->EndCell(); ++c) {

        ElementArray<Face> cell_faces = c->getFaces();
        for (ElementArray<Face>::iterator f = cell_faces.begin(); f != cell_faces.end(); ++f) {

            if (!f->getFaces(boundary_marker).empty()) {
                double T_b_ = boundary_conductivity(c->self(), f->self());
                const double gamma_b_ = gamma_b(c->self(), f->self());
                rMatrix gamma_b_rmatrix(1,1);
                gamma_b_rmatrix[0, 0] = gamma_b_;

                residual[entry_pressure.Index(c->self())] += f->Area()*(T_b_ * entry_pressure.Unknown(c->self()) - gamma_b_rmatrix);

            }
            else {
                ElementArray<Cell> cells = f->getCells();
                Cell neighbor_cell = c->Neighbour(f->self());
                const double T_ij = internal_conductivity(f->self());

                residual[entry_pressure.Index(c->self())] += f->Area() * T_ij * (entry_pressure.Unknown(c->self()) - entry_pressure.Unknown(neighbor_cell->self()));
            }

        }

    }

}

void OnePhaseModel::save_result(const std::string &filename) {
    try {
        mesh->Save(filename);
        std::cout << "Saving to " << filename << std::endl;
    } catch (...) {
        std::cout << "Wrong output path, saving to ../out.vtk" << std::endl;
        mesh->Save("../out.vtk");
    }
}


void OnePhaseModel::update_solution() {
    for (Mesh::iteratorCell c = mesh->BeginCell(); c != mesh->EndCell(); ++c) {
        entry_pressure.Value(c->self()) -= update[entry_pressure.Index(c->self())];
    }
}

Residual *OnePhaseModel::get_residual() {
    return &residual;
}

Sparse::Vector *OnePhaseModel::get_update() {
    return &update;
}




