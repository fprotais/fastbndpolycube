#include <string>
#include <iostream>
#include <fstream>
#include <ultimaille/all.h>
#include <OpenNL_psm/OpenNL_psm.h>


using namespace UM;
#define FOR(i, n) for(int i = 0; i < n; i++)


void soft_deformation(Triangles& m, FacetAttribute<int>& flag) {

	auto begin = std::chrono::steady_clock::now();
	auto context = nlNewContext();

	nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
	nlSolverParameteri(NL_NB_VARIABLES, NLint(3 * m.nverts()));
	//nlEnable(NL_VERBOSE);
	int N = m.nverts();
	nlBegin(NL_SYSTEM);

	nlBegin(NL_MATRIX);
	FOR(f, m.nfacets()) {
		int dim = flag[f] / 2;
		vec3 n = m.util.normal(f);
		FOR(d, 3) {
			FOR(fv, 3) {
				nlRowScaling(std::sqrt(n.norm()));
				nlBegin(NL_ROW);
				nlCoefficient(d * N + m.vert(f, fv), 1);
				nlCoefficient(d * N + m.vert(f, (fv + 1) % 3), -1);
				if (d != dim) nlRightHandSide(m.points[m.vert(f, fv)][d] - m.points[m.vert(f, (fv + 1) % 3)][d]);
				nlEnd(NL_ROW);
			}
		}
	}
	nlEnd(NL_MATRIX);
	nlEnd(NL_SYSTEM);

	auto buidlend = std::chrono::steady_clock::now();

	nlSolve();
	FOR(v, m.nverts()) FOR(d, 3) m.points[v][d] = nlGetVariable(d * N + v);
	nlDeleteContext(context);
	auto solveend = std::chrono::steady_clock::now();

	double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds> (solveend - begin).count()) / 1000;
	double buildtime = double(std::chrono::duration_cast<std::chrono::milliseconds> (buidlend - begin).count()) / 1000;
	double solvetime = double(std::chrono::duration_cast<std::chrono::milliseconds> (solveend - buidlend).count()) / 1000;
	std::cerr << "Running time: " << buildtime << " + " << solvetime << " = " << totaltime << "sec." << std::endl;
}


void hard_deformation(Triangles& m, FacetAttribute<int>& flag) {
	int N = m.nverts();

	auto begin = std::chrono::steady_clock::now();
	auto context = nlNewContext();

	DisjointSet ds(3 * m.nverts());
	FOR(f, m.nfacets()) {
		int dim = flag[f] / 2;
		FOR(fv, 3) ds.merge(dim * N + m.vert(f, fv), dim * N + m.vert(f, (fv + 1) % 3));
	}

	std::vector<int> idmap;
	int nb_variables = ds.get_sets_id(idmap);

	nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
	nlSolverParameteri(NL_NB_VARIABLES, NLint(nb_variables));
	nlSolverParameteri(NL_PRECONDITIONER, NL_PRECOND_JACOBI);
	//      nlSolverParameteri(NL_SYMMETRIC, NL_TRUE);
	//nlSolverParameteri(NL_MAX_ITERATIONS, NLint(100));
	//nlSolverParameterd(NL_THRESHOLD, 1e-5);
	//nlEnable(NL_VERBOSE);
	
	nlBegin(NL_SYSTEM);

	nlBegin(NL_MATRIX);
	FOR(f, m.nfacets()) {
		int dim = flag[f] / 2;
		double parea = std::sqrt(m.util.unsigned_area(f));
		FOR(d, 3) {
			if (dim == d) continue;
			FOR(fv, 3) {
				nlRowScaling(parea);
				nlBegin(NL_ROW);
				nlCoefficient(idmap[d * N + m.vert(f, fv)], 1);
				nlCoefficient(idmap[d * N + m.vert(f, (fv + 1) % 3)], -1);
				nlRightHandSide(m.points[m.vert(f, fv)][d] - m.points[m.vert(f, (fv + 1) % 3)][d]);
				nlEnd(NL_ROW);
			}
		}
	}
	nlEnd(NL_MATRIX);
	nlEnd(NL_SYSTEM);

	auto buidlend = std::chrono::steady_clock::now();

	nlSolve();

	FOR(v, m.nverts()) FOR(d, 3) m.points[v][d] = nlGetVariable(idmap[d * N + v]);
	nlDeleteContext(context);
	auto solveend = std::chrono::steady_clock::now();

	double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds> (solveend - begin).count()) / 1000;
	double buildtime = double(std::chrono::duration_cast<std::chrono::milliseconds> (buidlend - begin).count()) / 1000;
	double solvetime = double(std::chrono::duration_cast<std::chrono::milliseconds> (solveend - buidlend).count()) / 1000;
	std::cerr << "Running time: " << buildtime << " + " << solvetime << " = " << totaltime << "sec." << std::endl;

}



void hilbert_sort_mesh(Triangles& m) {
	HilbertSort sort(*m.points.data);
	std::vector<int> newid(m.nverts());
	std::iota(newid.begin(), newid.end(), 0);
	sort.apply(newid);
	std::vector<vec3> pts = *m.points.data;
	FOR(v, m.nverts()) m.points[newid[v]] = pts[v];
	FOR(f, m.nfacets()) FOR(fv, 3) m.vert(f, fv) = newid[m.vert(f, fv)];
}


void naive_flag_mesh(const Triangles& m, FacetAttribute<int>& flags) {
	const std::array<vec3, 6> AXES = { {{1.,0.,0.},{-1.,0.,0.},{0.,1.,0.},{0.,-1.,0.},{0.,0.,1.},{0.,0.,-1.}} };
	FOR(f, m.nfacets()) {
		vec3 n = m.util.normal(f);
		flags[f] = 0;
		double best = AXES[0] * n;
		FOR(i, 5) if (n * AXES[i + 1] > best) {
			best = n * AXES[i + 1];
			flags[f] = i + 1;
		}
	}
}

void fast() {
	Triangles m;
	read_by_extension("C:/fprotais/marchinghex/meshes/big.mesh", m);
	FacetAttribute<int> flag(m);
	naive_flag_mesh(m, flag);
	write_by_extension("flagging.geogram", m, { {},{{"flag", flag.ptr}},{} });
	write_by_extension("../bigmesh.obj", m, { });

	std::ofstream ofs("../bigmesh.flag");
	FOR(f, m.nfacets()) {
		ofs << flag[f] << std::endl;
	}
	ofs.close();
}

int main(int argc, char** argv) {

	if (argc != 4) {
		std::cout << "Usage is: " << argv[0] << " triangularmesh.ext flagfile deformedoutput.ext" << std::endl;
		return 1;
	}
	std::string inputfile = argv[1];
	std::string flagfile = argv[2];
	std::string outputfile = argv[3];

	
	Triangles m;
	read_by_extension(inputfile, m);
	m.delete_isolated_vertices();

	//hilbert_sort_mesh(m);

	FacetAttribute<int> flag(m, 0);
	std::ifstream ifs(flagfile);
	if (!ifs.is_open()) {
		std::cerr << "Failed opening of flags at : " << flagfile << std::endl;
		abort();
	}
	FOR(f, m.nfacets()) {
		if (ifs.eof()) break;
		ifs >> flag[f];
	}
	ifs.close();
	write_by_extension("flagging.geogram", m, { {},{{"flag", flag.ptr}},{} });

	hard_deformation(m, flag);

	write_by_extension(outputfile, m);

	return 0;
}

