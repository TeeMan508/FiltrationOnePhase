#include "inmost.h"
#include "iostream"

using namespace INMOST;


void fill_tag_K(const TagRealArray &K, Mesh &mesh) {
    rMatrix k_buf(3,3);
    for (Mesh::iteratorCell c = mesh.BeginCell(); c!= mesh.EndCell(); ++c) {
        K(*c, 3, 3) = rMatrix::Make(3,3, 1., 0., 0., 0., 1., 0., 0., 0., 1.);
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

void fill_initial_tag_p(TagReal &p, Mesh &mesh) {
    for (Mesh::iteratorCell c = mesh.BeginCell(); c!= mesh.EndCell(); ++c) {
        p[*c] = 0.;
    }
}

double lambda(const TagRealArray &tag_K, const Cell &c, const Face &f) {
    rMatrix normal(3,1);
    f.UnitNormal(normal.data());
    rMatrix const K = tag_K(c, 3, 3);

    return normal.DotProduct(K* normal);
}

rMatrix gamma(const TagRealArray &tag_K, const Cell &c, const Face &f, double neightbor_multiplier = 1.) {
    rMatrix normal(3,1), x_face(3,1), x_cell(3,1);;
    f.UnitNormal(normal.data());
    f.Centroid(x_face.data());
    c.Centroid(x_cell.data());
    rMatrix const K = tag_K(c, 3, 3);
    rMatrix delta = (x_face - x_cell) * neightbor_multiplier;
    double r = fabs(normal.DotProduct(delta));
    return K * normal - lambda(tag_K, c, f) / r * delta;
}

double internal_conductivity(const TagRealArray &tag_K, const Cell &c1, const Cell &c2, const Face &f) {
    rMatrix normal(3,1), x_face(3,1), x_cell_1(3,1), x_cell_2(3,1);;
    f.UnitNormal(normal.data());
    f.Centroid(x_face.data());
    c1.Centroid(x_cell_1.data());
    c2.Centroid(x_cell_2.data());

    rMatrix delta_1 = (x_face - x_cell_1) * 1.;
    rMatrix delta_2 = (x_face - x_cell_1) * -1.;

    double r1 = fabs(normal.DotProduct(delta_1));
    double r2 = fabs(normal.DotProduct(delta_2));

    double lambda_1 = lambda(tag_K, c1, f);
    double lambda_2 = lambda(tag_K, c2, f);
    // rMatrix gamma_1 = gamma(tag_K, c1, f, 1.);
    // rMatrix gamma_2 = gamma(tag_K, c2, f, -1.);

    return lambda_1 * lambda_2 / (lambda_1*r2 + lambda_2*r1);
}

double gamma_b(const TagRealArray &tag_K, const TagRealArray &tag_bc, const Cell &c, const Face &f) {
    rMatrix normal(3,1), x_face(3,1), x_cell(3,1);
    f.UnitNormal(normal.data());
    f.Centroid(x_face.data());
    c.Centroid(x_cell.data());

    rMatrix delta = (x_face - x_cell) * 1.;
    double r = fabs(normal.DotProduct(delta));

    Storage::real_array bc = tag_bc[f];
    double alpha = bc[0];
    double beta = bc[1];
    double gamma_bc = bc[2];

    double lamda_i = lambda(tag_K, c, f);

    return r * gamma_bc / (alpha * r + beta * lamda_i);

}

double boundary_conductivity(const TagRealArray &tag_K, const TagRealArray &tag_bc, const Cell &c, const Face &f) {
    rMatrix normal(3,1), x_face(3,1), x_cell(3,1);
    f.UnitNormal(normal.data());
    f.Centroid(x_face.data());
    c.Centroid(x_cell.data());

    rMatrix delta = (x_face - x_cell) * 1.;
    double r = fabs(normal.DotProduct(delta));

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
    Mesh m;
    try {
        m.Load(argv[1]);
    }
    catch (...) {
        std::cout << "Cannot read mesh file: " << argv[1] << std::endl;
        return 1;
    }

    MarkerType boundary_marker = m.CreateMarker();
    m.MarkBoundaryFaces(boundary_marker);

    TagRealArray tag_permeability = m.CreateTag("perm", DATA_REAL, CELL, NONE, 9);
    TagRealArray tag_boundary_cond = m.CreateTag("bc", DATA_REAL, FACE, FACE, 3);
    TagReal pressure = m.CreateTag("pressure", DATA_REAL, CELL, NONE, 1);

    fill_tag_K(tag_permeability, m);

    fill_tag_bc(tag_boundary_cond, m, boundary_marker);

    fill_initial_tag_p(pressure, m);

    BlockEntry p;
    p.AddTag(pressure);
    p.SetElementType(CELL);

    Automatizator aut("OnePhase");
    aut.RegisterEntry(p);
    aut.EnumerateEntries();
    Sparse::Vector Update("", aut.GetFirstIndex(), aut.GetLastIndex());
    Residual Resid("", aut.GetFirstIndex(), aut.GetLastIndex());


    for (Mesh::iteratorCell c = m.BeginCell(); c!= m.EndCell(); ++c) {
        ElementArray<Face> cell_faces = c->getFaces();
        for (ElementArray<Face>::iterator f = cell_faces.begin(); f != cell_faces.end(); ++f) {
            if (!f->getFaces(boundary_marker).empty()) {
                double T_b_ = boundary_conductivity(tag_permeability, tag_boundary_cond, c->self(), f->self());
                double gamma_b_ = gamma_b(tag_permeability, tag_boundary_cond, c->self(), f->self());
                rMatrix gamma_b_rmatrix(1,1);
                gamma_b_rmatrix[0, 0] = gamma_b_;

                Resid[p.Index(c->self())] += f->Area()*(T_b_ * p.Unknown(c->self()) - gamma_b_rmatrix);
                f->RemMarker(boundary_marker);
            }
            else {
                ElementArray<Cell> cells = f->getCells();
                Cell neighbor_cell;
                for (ElementArray<Cell>::iterator c2 = cells.begin(); c2 != cells.end(); ++c2) {
                    if (c2->self() != c->self()) neighbor_cell = c2->self();
                }
                double T_ij = internal_conductivity(tag_permeability, c->self(), neighbor_cell, f->self());

                Resid[p.Index(c->self())] += f->Area() * T_ij * (p.Unknown(c->self()) - p.Unknown(neighbor_cell->self()));
            }

        }

    }

    Solver solver(Solver::INNER_ILU2);
    solver.SetParameterEnum("verbosity", 1);
    solver.SetMatrix(Resid.GetJacobian());


    if (solver.Solve(Resid.GetResidual(), Update))
    {
        for (Mesh::iteratorCell c = m.BeginCell(); c != m.EndCell(); ++c)
            p.Value(c->self()) -= Update[p.Index(c->self())];
    }
    else std::cout << "Solve failed!" << std::endl;

    m.Save("../out.vtk");
    return 0;
}