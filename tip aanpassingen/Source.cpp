#include "Source.h"
#include <omp.h>
#include <chrono>
#include <random>

using namespace std;

void Model::init() {
	for (int h = 501; h < n + 501; h++) {
		radii.push_back(h * (radius_gal / n));
	}

	r_x.clear();
	r_y.clear();
    vector<double> phi;
    random_device phidev;
    std::mt19937 phirng(phidev());
    int phi_decimals = 10000;
    uniform_int_distribution<std::mt19937::result_type> phidist10000(0, phi_decimals);


	for (int i = 0; i < n; i++) {
		int r_phi_int = phidist10000(phirng);
		double phi_i = 2 * PI * (double) r_phi_int / phi_decimals;
		phi.push_back(phi_i);

        r_x.push_back(cos(phi[i]) * radii[i]);
        r_y.push_back(sin(phi[i]) * radii[i]);
	}

	rd_x = vector <double>(n);
	rd_y = vector <double>(n);
	rdd_x = vector <double>(n);
	rdd_y = vector <double>(n);
	m.clear();
	
	typedef mt19937 G;
	typedef gamma_distribution<> D;
	G g;
	double k = .7;
	double theta = 0.47;
	D d(k, theta);

	for (int i = 0; i < n; i++) {
		m.push_back(d(g) + 0.05);
	}
	

	Cell quad_tree = Cell(-radius_gal * 8, radius_gal * 8, double(16 * radius_gal), nullptr /*cp */,&cells);

	for (int i = 0; i < n; i++) {
		quad_tree.insert(r_x[i], r_y[i], m[i]);
	}

	quad_tree.compute_mass();

	for (int i = 0; i < n; i++) {
		double r_x_i_0 = r_x[i];
		double r_y_i_0 = r_y[i];
		double l_i_factor = sqrt(Gvar * m_sol * 1e3 / radii[i]) * 0;
		double rd_x_i_0 = -sin(phi[i]) * l_i_factor;
		double rd_y_i_0 = cos(phi[i]) * l_i_factor;
		double mdt_x_i_0, mdt_y_i_0;
		tie(mdt_x_i_0, mdt_y_i_0) = quad_tree.get_mass_distance_term(r_x_i_0, r_y_i_0);
		double rdd_x_i_0 = -Gvar * m_sol * mdt_x_i_0;
		double rdd_y_i_0 = -Gvar * m_sol * mdt_y_i_0;
		double r_x_i_1 = r_x_i_0 + rd_x_i_0 * dt + 0.5 * rdd_x_i_0 * pow(dt, 2);
		double r_y_i_1 = r_y_i_0 + rd_y_i_0 * dt + 0.5 * rdd_y_i_0 * pow(dt, 2);
		r_x[i] = r_x_i_1;
		r_y[i] = r_y_i_1;
		rd_x[i] = rd_x_i_0;
		rd_y[i] = rd_y_i_0;
		rdd_x[i] = rdd_x_i_0;
		rdd_y[i] = rdd_y_i_0;
	}


	for (Cell* c : cells) {
		free(c);
	}

	cells.clear();

	alpha = 1.5;
	f = 0.000001;
	double msum = 0;
	for (double i : m) {
		msum += i;
	}
	M_d = m_sol * msum;
	M_h = M_d / f - M_d;
	R = radius_gal;
	r_0 = radius_gal;
	rho_0 = (2 - alpha) * M_h / (2 * PI * pow(r_0, alpha) * pow(R, (2 - alpha)));

	Gvar = 6.67408e-11;
	omega_0 = sqrt(2 * PI * Gvar * rho_0 / (2 - alpha));

}

double Model::M(double rr) {
	const double PI = 3.141592653589793238463;
	double r = 2 * PI * rho_0 * pow(r_0, alpha) * pow(rr, (2.0 - alpha)) / (2.0 - alpha);
	return r;
}


#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
void Model::model_loop() {
	

	Cell quad_tree = Cell(-radius_gal * 5, radius_gal * 5, double(10 * radius_gal), nullptr /*cp */,&cells);

	for (int i = 0; i < n; i++) {
		quad_tree.insert(r_x[i], r_y[i], m[i]);
		
	}


	quad_tree.compute_mass();

	vector<tuple<double, double>> mdi_j(n, make_tuple(0.0,0.0));
	vector<double> rr(n,0.0);
	vector<double> rd_x_i_j(n, 0.0);
	vector<double> rd_y_i_j(n, 0.0);
	vector<double> r_x_i_jp1(n, 0.0);
	vector<double> r_y_i_jp1(n, 0.0);
	vector<double> rdd_x_i_j(n, 0.0);
	vector<double> rdd_y_i_j(n, 0.0);
	
	/*
	fill(mdi_j.begin(), mdi_j.end(), make_tuple(0.0,0.0));
	fill(rr.begin(), rr.end(), 0.0);
	fill(rd_x_i_j.begin(), rd_x_i_j.end(), 0.0);
	fill(rd_y_i_j.begin(), rd_y_i_j.end(), 0.0);
	fill(r_x_i_jp1.begin(), r_x_i_jp1.end(), 0.0);
	fill(r_y_i_jp1.begin(), r_y_i_jp1.end(), 0.0);
	fill(rdd_x_i_j.begin(), rdd_x_i_j.end(), 0.0);
	fill(rdd_y_i_j.begin(), rdd_y_i_j.end(), 0.0);
	*/
	


	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		
		rr[i] = sqrt(pow(r_x[i], 2) + pow(r_y[i], 2));

		//calculate rdd_i_j
		mdi_j[i] = quad_tree.get_mass_distance_term(r_x[i], r_y[i]);

		rdd_x_i_j[i] = -Gvar * m_sol * get<0>(mdi_j[i]);
		rdd_y_i_j[i] = -Gvar * m_sol * get<1>(mdi_j[i]);

		//calculate rd_i_j
		rd_x_i_j[i] = rd_x[i] + 0.5 * (rdd_x[i] + rdd_x_i_j[i]) * dt;
		rd_y_i_j[i] = rd_y[i] + 0.5 * (rdd_y[i] + rdd_y_i_j[i]) * dt;

		//calculate r_i_(j + 1)
		r_x_i_jp1[i] = r_x[i] + rd_x_i_j[i] * dt + 0.5 * rdd_x_i_j[i] * pow(dt,2);
		r_y_i_jp1[i] = r_y[i] + rd_y_i_j[i] * dt + 0.5 * rdd_y_i_j[i] * pow(dt, 2);
		r_x[i] = r_x_i_jp1[i];
		r_y[i] = r_y_i_jp1[i];
		rd_x[i] = rd_x_i_j[i];
		rd_y[i] = rd_y_i_j[i];
		rdd_x[i] = rdd_x_i_j[i];
		rdd_y[i] = rdd_y_i_j[i];
	}
	for (Cell* c : cells) {
		free(c);
	}
	cells.clear();

}
#pragma clang diagnostic pop


int main() {
    auto start = chrono::high_resolution_clock::now();
    Model model = Model();

	model.init();
	
	ofstream outFile;
	outFile.open("./huts.csv");
	int Iterations = 500;
	


	cout << "Starting now..." << endl;

    cout << model.n << endl;

	for (int h = 0; h < Iterations; h++) {

		cout << h << "\n";
		
		for (unsigned int j = 0; j < model.r_x.size(); j++) {
			outFile << model.r_x[j] << ";" << model.r_y[j] << ";";
		}
		outFile << endl;

		model.model_loop();
		
	}
    auto stop = chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << duration.count() / 1000 << std::endl;
	cout << endl << "Done";
	
	outFile.close();

}







