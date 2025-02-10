#include "inmost.h"
#include "iostream"

using namespace INMOST;


void fill_tag_K(const TagRealArray &K, Mesh &mesh) {
    rMatrix k_buf(3,3);
    for (Mesh::iteratorCell c = mesh.BeginCell(); c!= mesh.EndCell(); ++c) {
        K(*c, 3, 3) = rMatrix::Make(3,3, 2., 0., 0., 0., 2., 0., 0., 0., 2.);
    }
}

void fill_tag_K_prepared(const TagRealArray &K, Mesh &mesh) {
    TagRealArray tag_perm = mesh.GetTag("PERM");
    for (Mesh::iteratorCell c = mesh.BeginCell(); c!= mesh.EndCell(); ++c) {
        K(*c, 3, 3) = rMatrix::FromTensor(tag_perm[*c].data(), tag_perm[*c].size());
    }
}

void fill_tag_bc(TagRealArray &bc, Mesh &mesh, MarkerType &boundary_marker) {
    rMatrix normal(3,1);
    for (Mesh::iteratorFace face = mesh.BeginFace(); face != mesh.EndFace(); ++face) {
        if (face->getFaces(boundary_marker).empty()) continue;
        face->UnitNormal(normal.data());

        if (fabs(normal[0] + 1.) < 1e-5) {
            // Left Boundary
            bc(*face, 3, 1) = rMatrix::Make(
                3, 1, 1., 0., 1.); // TODO: REFACTOR BC
            continue;
        }
        if (fabs(normal[0] - 1.) < 1e-5) {
            // Right Boundary
            bc(*face, 3, 1) = rMatrix::Make(
                3, 1, 1., 0., 0);
            continue;
        }
        // Other boundaries
        bc(*face, 3, 1) = rMatrix::Make(
                3, 1, 0., 1., 0.);
    }
}

void fill_tag_bc_prepared(TagRealArray &bc, Mesh &mesh, MarkerType &boundary_marker) {
    TagRealArray tag_bc_prepared = mesh.GetTag("BOUNDARY_CONDITION");

    for (Mesh::iteratorFace face = mesh.BeginFace(); face != mesh.EndFace(); ++face) {

        if (face->getFaces(boundary_marker).empty()) continue;
        bc(*face, 3, 1) = rMatrix::FromVector(tag_bc_prepared[*face].data(), tag_bc_prepared[*face].size());
    }
}


void fill_initial_tag_p(TagReal &p, Mesh &mesh) {
    // for (Mesh::iteratorCell c = mesh.BeginCell(); c!= mesh.EndCell(); ++c) {
    //     p[*c] = 0.;
    // }
}

double lambda(const TagRealArray &tag_K, const Cell &c, const Face &f) {
    rMatrix normal(3,1);
    f.OrientedUnitNormal(c, normal.data());
    rMatrix const K = tag_K(c, 3, 3);
    // std::cout << K(0, 0) << std::endl;
    // std::cout << K(1, 1) << std::endl;
    // std::cout << K(2, 2) << std::endl;
    return normal.DotProduct(K * normal);
}

// rMatrix gamma(const TagRealArray &tag_K, const Cell &c, const Face &f, double neightbor_multiplier = 1.) {
//     rMatrix normal(3,1), x_face(3,1), x_cell(3,1);;
//     f.UnitNormal(normal.data());
//     f.Centroid(x_face.data());
//     c.Centroid(x_cell.data());
//     rMatrix const K = tag_K(c, 3, 3);
//     rMatrix delta = (x_face - x_cell) * neightbor_multiplier;
//     double r f= fabs(normal.DotProduct(delta));
//     return K * normal - lambda(tag_K, c, f) / r * delta;
// }

double internal_conductivity(const TagRealArray &tag_K, const Face &f) {
    rMatrix normal(3,1), normal_2(3, 1), x_face(3,1), x_cell_1(3,1), x_cell_2(3,1);
    // f.UnitNormal(normal_1.data());
    Cell c1 = f->BackCell();
    Cell c2 = f->FrontCell();

    f.UnitNormal(normal.data()); // todo: почему так
    // f.UnitNormal(normal.data());
    f.Centroid(x_face.data());
    c1.Centroid(x_cell_1.data());
    c2.Centroid(x_cell_2.data());

    rMatrix delta_1 = (x_face - x_cell_1) * 1.;
    rMatrix delta_2 = (x_face - x_cell_2) * -1.;

    // double r1 = fabs(normal.DotProduct(delta_1));
    // double r2 = fabs(normal.DotProduct(delta_2));

    double r1 = normal.DotProduct(delta_1);
    double r2 = normal.DotProduct(delta_2);
    // std::cout << r1 << " " << r2 << std::endl;

    double lambda_1 = lambda(tag_K, c1, f);
    double lambda_2 = lambda(tag_K, c2, f);
    // std::cout << lambda_1 << " " << lambda_2 << std::endl;
    // rMatrix gamma_1 = gamma(tag_K, c1, f, 1.);
    // rMatrix gamma_2 = gamma(tag_K, c2, f, -1.);
    // std::cout<< lambda_1 * lambda_2 / (lambda_1*r2 + lambda_2*r1) << std::endl;
    // return lambda_1 / 2 / r1;

    return lambda_1 * lambda_2 / (lambda_1*r2 + lambda_2*r1);
}

double gamma_b(const TagRealArray &tag_K, const TagRealArray &tag_bc, const Cell &c, const Face &f) {
    rMatrix normal(3,1), x_face(3,1), x_cell(3,1);
    f.UnitNormal(normal.data());
    // normal /= sqrt(normal.DotProduct(normal));
    f.Centroid(x_face.data());
    c.Centroid(x_cell.data());

    rMatrix delta = (x_face - x_cell) * 1.;
    double r = normal.DotProduct(delta);
    // double r = normal.DotProduct(delta);

    Storage::real_array bc = tag_bc[f];
    double alpha = bc[0];
    double beta = bc[1];
    double gamma_bc = bc[2];
    // std::cout << alpha << " " << beta << " " << gamma_bc << std::endl;
    double lamda_i = lambda(tag_K, c, f);

    return lamda_i * gamma_bc / (alpha * r + beta * lamda_i);

}

double boundary_conductivity(const TagRealArray &tag_K, const TagRealArray &tag_bc, const Cell &c, const Face &f) {
    rMatrix normal(3,1), x_face(3,1), x_cell(3,1);
    f.UnitNormal(normal.data());
    f.Centroid(x_face.data());
    c.Centroid(x_cell.data());

    rMatrix delta = (x_face - x_cell) * 1.;
    double r = normal.DotProduct(delta);

    // double r = normal.DotProduct(delta);
    // std::cout << r << std::endl; 0.015625
    Storage::real_array bc = tag_bc[f];
    double alpha = bc[0];
    double beta = bc[1];
    double gamma_bc = bc[2];

    double lamda_i = lambda(tag_K, c, f);

    return alpha*lamda_i / (alpha*r + beta*lamda_i);
}


int main(int argc, char ** argv) {
    if (argc < 1) {
        std::cout << "Usage: " << argv[0] << " mesh_file " << std::endl;
        return 1;
    }
    Mesh *m = new Mesh;
    // m->SetFileOption("VTK_OUTPUT_FACES","1");
    try {
        m->Load(argv[1]);
    }
    catch (...) {
        std::cout << "Cannot read mesh file: " << argv[1] << std::endl;
        return 1;
    }

    MarkerType boundary_marker = m->CreateMarker();
    m->MarkBoundaryFaces(boundary_marker);

    TagRealArray tag_permeability = m->CreateTag("perm", DATA_REAL, CELL, NONE, 9);
    TagRealArray tag_boundary_cond = m->CreateTag("bc", DATA_REAL, FACE, FACE, 3);
    TagReal pressure = m->CreateTag("pressure", DATA_REAL, CELL, NONE, 1);

    if (m->HaveTag("PERM")) fill_tag_K_prepared(tag_permeability, *m);
    else fill_tag_K(tag_permeability, *m);

    if (m->HaveTag("BOUNDARY_CONDITION")) fill_tag_bc_prepared(tag_boundary_cond, *m, boundary_marker);
    else fill_tag_bc(tag_boundary_cond, *m, boundary_marker);

    // fill_initial_tag_p(pressure, *m);

    BlockEntry p;
    p.AddTag(pressure);
    p.SetElementType(CELL);

    Automatizator aut("OnePhase");
    aut.RegisterEntry(p);
    aut.EnumerateEntries();
    Sparse::Vector Update("", aut.GetFirstIndex(), aut.GetLastIndex());
    Residual Resid("", aut.GetFirstIndex(), aut.GetLastIndex());


    for (Mesh::iteratorCell c = m->BeginCell(); c!= m->EndCell(); ++c) {
        ElementArray<Face> cell_faces = c->getFaces();
        for (ElementArray<Face>::iterator f = cell_faces.begin(); f != cell_faces.end(); ++f) {
            if (!f->getFaces(boundary_marker).empty()) {
                double T_b_ = boundary_conductivity(tag_permeability, tag_boundary_cond, c->self(), f->self());
                double gamma_b_ = gamma_b(tag_permeability, tag_boundary_cond, c->self(), f->self());
                rMatrix gamma_b_rmatrix(1,1);
                gamma_b_rmatrix[0, 0] = gamma_b_;

                Resid[p.Index(c->self())] += f->Area()*(T_b_ * p.Unknown(c->self()) - gamma_b_rmatrix);
                // f->RemMarker(boundary_marker);
            }
            else {
                ElementArray<Cell> cells = f->getCells();
                // todo: c->Neighbour(face)
                Cell neighbor_cell;
                for (ElementArray<Cell>::iterator c2 = cells.begin(); c2 != cells.end(); ++c2) {
                    if (c2->self() != c->self()) neighbor_cell = c2->self();
                }
                double T_ij = internal_conductivity(tag_permeability, f->self());

                Resid[p.Index(c->self())] += f->Area() * T_ij * (p.Unknown(c->self()) - p.Unknown(neighbor_cell->self()));
            }

        }

    }

    Solver solver(Solver::INNER_ILU2);
    solver.SetParameterEnum("verbosity", 1);
    solver.SetMatrix(Resid.GetJacobian());

    if (solver.Solve(Resid.GetResidual(), Update))
    {
        for (Mesh::iteratorCell c = m->BeginCell(); c != m->EndCell(); ++c)
            p.Value(c->self()) -= Update[p.Index(c->self())];

        std::cout << "Saved ../out.vtk" << std::endl;
    }
    else std::cout << "Solve failed!" << std::endl;

    TagReal tag_reference = m->GetTag("REFERENCE_SOLUTION");
    TagReal tag_error = m->CreateTag("error", DATA_REAL, CELL, NONE, 1);
    // TagReal tag_bc_2 = m->CreateTag("bc_2", DATA_REAL, CELL, CELL, 1);

    // for (Mesh::iteratorCell c = m->BeginCell(); c != m->EndCell(); ++c) {
    //     tag_error[*c] = fabs(pressure[*c] - tag_reference[*c]);
    //     if (!c->getFaces(boundary_marker).empty()) {
    //         ElementArray<Face> faces = c->getFaces();
    //         Face f = faces.begin()->self();
    //         tag_error[*c] = rMatrix::FromVector(tag_boundary_cond[f], 3);
    //     }
    // }


    m->Save("../out.vtk");
    return 0;
}