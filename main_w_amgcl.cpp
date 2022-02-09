#include <string>
#include <iostream>
#include <fstream>
#include <ultimaille/all.h>
#include <OpenNL_psm/OpenNL_psm.h>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/preconditioner/dummy.hpp>
#include <amgcl/preconditioner/cpr.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/coarsening/ruge_stuben.hpp>
#include <amgcl/coarsening/aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/ILUP.hpp>
#include <amgcl/relaxation/Chebyshev.hpp>
#include <amgcl/relaxation/damped_jacobi.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/io/mm.hpp>


#include <algorithm>
#include <execution>

using namespace UM;
#define FOR(i, n) for(int i = 0; i < n; i++)

#pragma warning( disable : 4244)


struct triplet {
	int i, j;
	double v;
};


void triplet2CRS(
	std::vector<triplet>& trip, int nb_variables,
	std::vector<int>& ptr,
	std::vector<int>& col,
	std::vector<double>& val
) {
	ptr.resize(nb_variables + 1); ptr[0] = 0;
	col.reserve(nb_variables);
	val.reserve(nb_variables);

	std::sort(std::execution::par_unseq, trip.begin(), trip.end(), [](const triplet& a, const triplet& b) -> bool {
		if (a.i < b.i) return true;
		if (a.i == b.i) return a.j < b.j;
		return false;
		});

	int actual_i = 0;
	int actual_j = -1;
	for (auto t : trip) {
		if (t.i == actual_i && t.j == actual_j) {
			val.back() += t.v;
		} 
		else {
			if (t.i > actual_i) {
				um_assert(t.i < nb_variables);
				for (int i = actual_i; i < t.i; i++) ptr[i + 1] = col.size();
				actual_i = t.i;
			}
			um_assert(t.j < nb_variables);
			col.push_back(t.j);
			val.push_back(t.v);
			actual_j = t.j;
		}
	}
	ptr.back() = col.size();
}

void fast_hard_deformation(Triangles& m, FacetAttribute<int>& flag) {
	auto begin = std::chrono::steady_clock::now();

	int N = m.nverts();
	DisjointSet ds(3 * m.nverts());
	FOR(f, m.nfacets()) {
		int dim = flag[f] / 2;
		FOR(fv, 3) ds.merge(dim * N + m.vert(f, fv), dim * N + m.vert(f, (fv + 1) % 3));
	}

	std::vector<int> idmap;
	int nb_variables = ds.get_sets_id(idmap);

	std::vector<triplet> trip; trip.reserve(3 * 3 * 4 * m.nfacets());
	std::vector<double> rhs(nb_variables);

	FOR(f, m.nfacets()) {
		int dim = flag[f] / 2;
		//double parea = 1; std::sqrt(m.util.unsigned_area(f));
		FOR(d, 3) {
			if (dim == d) continue;
			FOR(fv, 3) {
				int i = idmap[d * N + m.vert(f, fv)];
				int j = idmap[d * N + m.vert(f, (fv + 1) % 3)];
				double edge_in_d = m.points[m.vert(f, fv)][d] - m.points[m.vert(f, (fv + 1) % 3)][d];
				trip.push_back({ i, i, 1 });
				trip.push_back({ i, j, -1 });
				trip.push_back({ j, j, 1 });
				trip.push_back({ j, i, -1 });
				rhs[i] += 1 * edge_in_d;
				rhs[j] -= 1 * edge_in_d;
			}
		}
	}

	std::vector<int>    ptr, col;
	std::vector<double> val;
	triplet2CRS(trip, nb_variables, ptr, col, val);


	//amgcl::io::mm_write("matrix.mm",std::tie(nb_variables, ptr, col, val));
	auto buidlend = std::chrono::steady_clock::now();

	typedef amgcl::backend::builtin<double> Backend;
	typedef amgcl::make_solver<
		// Use AMG as preconditioner:
		amgcl::preconditioner::dummy<Backend>,
		//amgcl::amg<
		//Backend,
		//amgcl::coarsening::smoothed_aggregation,
		////amgcl::coarsening::ruge_stuben,
		////amgcl::coarsening::aggregation,
		////amgcl::relaxation::damped_jacobi
		////amgcl::relaxation::ilup
		////amgcl::relaxation::chebyshev
		//amgcl::relaxation::spai0
		//>,
		// And BiCGStab as iterative solver:
		//amgcl::solver::bicgstab<Backend>
		amgcl::solver::cg<Backend>
	> Solver;


	Solver solve(std::tie(nb_variables, ptr, col, val));
	std::vector<double> x(nb_variables, 0.0);
	int    iters;
	double error;
	

	std::tie(iters, error) = solve(rhs, x);
	FOR(v, m.nverts()) FOR(d, 3) m.points[v][d] = x[idmap[d * N + v]];

	auto solveend = std::chrono::steady_clock::now();

	double totaltime = double(std::chrono::duration_cast<std::chrono::milliseconds> (solveend - begin).count()) / 1000;
	double buildtime = double(std::chrono::duration_cast<std::chrono::milliseconds> (buidlend - begin).count()) / 1000;
	double solvetime = double(std::chrono::duration_cast<std::chrono::milliseconds> (solveend - buidlend).count()) / 1000;
	std::cerr << "Running time: " << buildtime << " + " << solvetime << " = " << totaltime << "sec." << std::endl;
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

	fast_hard_deformation(m, flag);

	write_by_extension(outputfile, m);

	return 0;
}

