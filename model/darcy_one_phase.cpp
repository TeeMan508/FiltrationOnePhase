#include "darcy_one_phase.h"

#include <utility>
#include "../utils/calculations/calculations.h"
// #include "../solver/solver.h"


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
        "Boundary Conditions", DATA_REAL, FACE, FACE, 7);

    TagRealArray tag_stiffness_tensor = this->mesh.CreateTag(
        "Stiffness Tensor", DATA_REAL, CELL, NONE, 36);

	TagReal tag_dual_volume = this->mesh.CreateTag(
		"Dual Volume", DATA_REAL, NODE, NODE, 1);

	TagRealArray tag_stress = this->mesh.CreateTag(
		"Stress", DATA_REAL, CELL, CELL, 6);

	TagRealArray tag_strain = this->mesh.CreateTag(
		"Strain", DATA_REAL, CELL, CELL, 6);
	// this->save_result("out.vtu");
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
                7, 1, 1., 0., 0., 1., -5., 1., -1.);
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

void OnePhaseModel::assemble_internal_part(BlockEntry &uvw, Residual &residual) {
	TagRealArray tag_displacement = this->mesh.GetTag("Displacement");
	TagRealArray tag_stiffness_tensor = this->mesh.GetTag("Stiffness Tensor");
	TagReal tag_dual_volume = this->mesh.GetTag("Dual Volume");
	TagRealArray tag_boundary_conditions = this->mesh.GetTag("Boundary Conditions");

	MarkerType belong_to_cell_marker = this->mesh.CreateMarker();

	ElementArray<Face> faces, neightbor_faces;
	ElementArray<Edge> edges, neightbor_edges;
	ElementArray<Node> nodes;

	double A, volume, distance;

	rMatrix
	stiffness_tensor(6, 6),
	x_cell_center(3, 1), x_edge_center(3, 1), x_face_center(3, 1), xt(3, 1), x_begin_node(3, 1), x_end_node(3, 1),
	l(3, 1), normal(3, 1);

	vMatrix gradient_(9, 1), gradient_corrected(9, 1),
	traction(3, 1), strain(6, 1), stress(6, 1),
	dU(3, 1), gradient_l(3, 1);

	for (Mesh::iteratorCell c = this->mesh.BeginCell(); c != this->mesh.EndCell(); ++c)
	{
		c->Centroid(x_cell_center.data());
		edges = c->getEdges();
		faces = c->getFaces();
		faces.SetMarker(belong_to_cell_marker); // поставили пометку граням ячейки что они ей принадлежат
		gradient(c->self(), uvw, gradient_); //расчет градиента
		stiffness_tensor = rMatrix::FromTensor(tag_stiffness_tensor[*c].data(),
			tag_stiffness_tensor[*c].size(), 6); // вытакиваем 6х6 матрицу
		for (ElementArray<Edge>::iterator e = edges.begin(); e != edges.end(); ++e)
		{
			// выделяем только 2 соседние которые образуют треугольники
			Node begin_node = e->getBeg(), end_node = e->getEnd();
			//edge-vector between nodes
			begin_node.Centroid(x_begin_node.data());
			end_node.Centroid(x_end_node.data());
			l = x_end_node - x_begin_node;
			e->Centroid(x_edge_center.data());
			distance = l.DotProduct(l);
			//faces of the edge sharing the cell
			neightbor_faces = e->getFaces(belong_to_cell_marker); // грани с пометкой берем (те что принадлежат ячейке)
			// только грани которые принадлежат ячейке
			for (ElementArray<Face>::iterator f = neightbor_faces.begin(); f != neightbor_faces.end(); ++f)
			{

				f->Barycenter(x_face_center.data());
				normal = (x_cell_center - x_edge_center).CrossProduct(x_face_center - x_edge_center); //считаем нормаль
				A = normal.FrobeniusNorm();
				if (not A) continue;//internal flux

				normal /= A;
				if (l.DotProduct(normal) < 0.0) normal *= -1.0; // ориентация
				A *= 0.5; // площадь
				xt = (x_cell_center + x_edge_center + x_face_center) / 3.0;
				volume = A * normal.DotProduct(xt) / 3.0; // объем дуальной ячейки
				//добавим ГУ
				gradient_l = gradient_.Repack(3, 3) * l;  // 9х1 в 3х3 (G*l)
				dU = uvw.Unknown(end_node) - uvw.Unknown(begin_node) - gradient_l; //разность неизвестных 2-х соседних узлов
				gradient_corrected = gradient_ + dU.Kronecker(l / distance); // к G добавляем поправку
				gradient_to_strain(gradient_corrected, strain);
				stress = stiffness_tensor * strain;
				stress_to_traction(stress, normal, traction);

				tag_dual_volume[end_node] -= volume;
				tag_dual_volume[begin_node] += volume;

				residual[uvw.Index(end_node)] += A * traction; // невязки второго элемента
				residual[uvw.Index(begin_node)] -= A * traction;
			}
		}
		faces.RemMarker(belong_to_cell_marker); // стираем метки когда закончили с ячейкой
	}
	mesh.ReleaseMarker(belong_to_cell_marker, FACE);
}

void OnePhaseModel::assemble_boundary_part(BlockEntry &uvw, Residual &residual) {
	TagRealArray tag_displacement = this->mesh.GetTag("Displacement");
	TagRealArray tag_stiffness_tensor = this->mesh.GetTag("Stiffness Tensor");
	TagReal tag_dual_volume = this->mesh.GetTag("Dual Volume");
	TagRealArray tag_boundary_conditions = this->mesh.GetTag("Boundary Conditions");

	MarkerType belong_to_cell_marker = this->mesh.CreateMarker();

	ElementArray<Face> faces, neightbor_faces;
	ElementArray<Edge> edges, neightbor_edges;
	ElementArray<Node> nodes;

	double A, volume, distance;
	rMatrix
	stiffness_tensor(6, 6),
	x_cell_center(3, 1), x_edge_center(3, 1), x_face_center(3, 1), xt(3, 1), x_begin_node(3, 1), x_end_node(3, 1),
	l(3, 1), normal(3, 1);

	vMatrix gradient_(9, 1), gradient_corrected(9, 1),
	traction(3, 1), strain(6, 1), stress(6, 1),
	dU(3, 1), gradient_l(3, 1);

	MarkerType boundary_marker;
	mesh.MarkBoundaryFaces(boundary_marker); // граничные грани метим
	rMatrix xn(3, 1), alpha(3, 3), beta(3, 3), gamma(3, 1), alphasum(3, 3), ialphasum(3, 3), ibeta(3, 3);
	vMatrix bccond(3, 1);

	for (Mesh::iteratorNode v = mesh.BeginNode(); v != mesh.EndNode(); ++v)
	{
		faces = v->getFaces(boundary_marker); // берем граничные с маркером
		if (!faces.empty())
		{
			v->Centroid(xn.data());
			edges = v->getEdges();
			edges.SetMarker(belong_to_cell_marker); // помечаем ребра
			alphasum.Zero();
			bccond.Zero();

			for (ElementArray<Face>::iterator f = faces.begin(); f != faces.end(); ++f)
			{
				// итерации по граням
				Cell c = f->BackCell(); // достаем ячейку
				f->Centroid(x_face_center.data()); // центр ребра
				f->UnitNormal(normal.data()); // вытаскиваем единичная нормаль к грани
				extract_boundary_conditions(f->self(), alpha, beta, gamma); // вытаскиваем ГУ
				gradient(c, uvw, gradient_); //градиент
				stiffness_tensor = rMatrix::FromTensor(tag_stiffness_tensor[c].data(),
					tag_stiffness_tensor[c].size(), 6); // тензор жесткости

				gradient_to_strain(gradient_, strain); // тенз деформации
				stress = stiffness_tensor * strain;
				stress_to_traction(stress, normal, traction); // сила
				if (beta.FrobeniusNorm()) // если бета ненулевая
				{
					ibeta = beta.PseudoInvert();
					traction += ibeta * (alpha * uvw.Unknown(v->self()) - gamma - beta * traction); //модифицируем силу с учетом ГУ
				}
				// выделим ребра которые и к узлу и к грани одновременно относятся
				neightbor_edges = f->getEdges(belong_to_cell_marker);
				for (ElementArray<Edge>::iterator e = neightbor_edges.begin(); e != neightbor_edges.end(); ++e)
				{
					e->Centroid(x_edge_center.data()); // центр ребра
					A = (xn - x_edge_center).CrossProduct(x_face_center - x_edge_center).FrobeniusNorm() * 0.5; // площадь треугольничка
					xt = (xn + x_edge_center + x_face_center) / 3.0;
					volume = A * normal.DotProduct(xt) / 3.0;
					tag_dual_volume[*v] += volume;
					residual[uvw.Index(v->self())] -= A * traction;
					if (alpha.FrobeniusNorm())
					{
						bccond += (alpha * uvw.Unknown(v->self()) - beta * traction - gamma) * A; // сумма граничных условий - невязка
						alphasum += alpha * A; // суммируем
					}
				}
			}

			if (alphasum.FrobeniusNorm())
			{
				ialphasum = alphasum.PseudoInvert(); // если условия на перемещение были, делаем псевдообращение
				residual[uvw.Index(v->self())] += ialphasum * (bccond - alphasum * residual[uvw.Index(v->self())]); // учитываем ГУ
			}
			edges.RemMarker(belong_to_cell_marker);

		}
	}
	mesh.ReleaseMarker(boundary_marker, FACE);
	mesh.ReleaseMarker(belong_to_cell_marker, NODE);
}

void OnePhaseModel::assemble_system(BlockEntry& uvw, Residual& residual) {
	this->assemble_internal_part(uvw, residual);
	this->assemble_boundary_part(uvw, residual);

}

void OnePhaseModel::setup_model() {
	this->setup_tags();
	this->setup_boundary_conditions();
	this->setup_initial_conditions();
}

void OnePhaseModel::update_stress_and_strain_tags() {
	TagRealArray tag_displacement = this->mesh.GetTag("Displacement");
	TagRealArray tag_stiffness_tensor = this->mesh.GetTag("Stiffness Tensor");
	TagRealArray tag_stress = this->mesh.GetTag("Stress");
	TagRealArray tag_strain = this->mesh.GetTag("Strain");
	BlockEntry uvw;
	uvw.AddTag(tag_displacement);
	uvw.SetElementType(NODE);
	vMatrix gradient_(9, 1);
	rMatrix stiffness_tensor(6, 6), stress(6, 1), strain(6, 1);

	for (Mesh::iteratorCell c = mesh.BeginCell(); c != mesh.EndCell(); ++c)
	{
		gradient(c->self(), uvw, gradient_);
		stiffness_tensor = rMatrix::FromTensor(tag_stiffness_tensor[*c].data(), tag_stiffness_tensor[*c].size(), 6);
		strain[0] = gradient_[0].GetValue();
		strain[1] = gradient_[4].GetValue();
		strain[2] = gradient_[8].GetValue();
		strain[3] = gradient_[5].GetValue() + gradient_[7].GetValue();
		strain[4] = gradient_[2].GetValue() + gradient_[6].GetValue();
		strain[5] = gradient_[1].GetValue() + gradient_[3].GetValue();

		stress = stiffness_tensor * strain;

		std::swap(strain[3], strain[5]);
		std::swap(strain[4], strain[5]);
		std::swap(stress[3], stress[5]);
		std::swap(stress[4], stress[5]);

		tag_stress(*c, 6, 1) = stress;
		tag_strain(*c, 6, 1) = strain;
	}

}

bool OnePhaseModel::solve_elasticity_model() {
	bool success = false;
	TagRealArray tag_displacement = this->mesh.GetTag("Displacement");

	BlockEntry uvw; // неизвестные с производными
	uvw.AddTag(tag_displacement);
	uvw.SetElementType(NODE);

	Automatizator aut("Elasticity");
	aut.RegisterEntry(uvw);
	aut.EnumerateEntries(); // нумеруем
	Sparse::Vector Update("", aut.GetFirstIndex(), aut.GetLastIndex());
	Residual Resid("", aut.GetFirstIndex(), aut.GetLastIndex());

	Solver solver(Solver::INNER_ILU2);
	solver.SetParameterEnum("verbosity", 1);
	solver.SetMatrix(Resid.GetJacobian());
	// CustomSolver solver(Resid, Update);

	// solver.setup();

	this->assemble_system(uvw, Resid);
	solver.SetMatrix(Resid.GetJacobian());

	if (solver.Solve(Resid.GetResidual(), Update))
	{

		for (Mesh::iteratorNode v = this->mesh.BeginNode(); v != this->mesh.EndNode(); ++v)
			uvw.Value(v->self()) -= Update[uvw.Index(v->self())];
		success = true;
	}
	else std::cout << "Solve failed!" << std::endl;
	return success;
}

Tag OnePhaseModel::get_tag(const std::string& name) {
	return this->mesh.GetTag(name);
}

void OnePhaseModel::save_result(const std::string &filename) {

	this->mesh.Save(filename);
}







