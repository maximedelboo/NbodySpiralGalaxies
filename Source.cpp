#include "Source.h"
#include <omp.h>
using namespace std;

void Model::init() {
	for (int h = 0; h < 50; h++) {
		radii.push_back(max_radius / n_radius * double(h + 1.0));
		circumferences.push_back(2 * PI * radii.back());
		total_circumference += circumferences.back();
	}

	r_x.clear();
	r_y.clear();

	for (int i = 0; i < 50; i++) {
		double radius = radii[i];
		double j_tot = floor(rough_n / pow(n_radius, 2) * ((double) rand() / RAND_MAX + 0.5));
		for (int j = 0; j < j_tot; j++) {
			double phi = ((double) rand() / RAND_MAX) * 2 * PI;
			r_x.push_back(cos(phi)* radius);
			r_y.push_back(sin(phi) * radius);
		}
		
	}
	
	n = r_x.size();
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
	

	Cell quad_tree = Cell(-max_radius * 5, max_radius * 5, double(10 * max_radius), nullptr /*cp */,&cells);

	for (int i = 0; i < n; i++) {
		quad_tree.insert(r_x[i], r_y[i], m[i]);
	}

	quad_tree.compute_mass();

	for (int i = 0; i < n; i++) {
		double r_x_i_0 = r_x[i];
		double r_y_i_0 = r_y[i];
		double l_i_factor = ((2 * ((double)rand() / RAND_MAX) + 1) * 1e5) / sqrt(pow(r_x_i_0, 2) + pow(r_y_i_0, 2)) * 0.1;
		double rd_x_i_0 = r_y_i_0 * l_i_factor;
		double rd_y_i_0 = -r_x_i_0 * l_i_factor;
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
	m_sol = 1.98847e30;
	f = 0.000001;
	double msum = 0;
	for (double i : m) {
		msum += i;
	}
	M_d = m_sol * msum;
	M_h = M_d / f - M_d;
	R = max_radius;
	r_0 = max_radius;
	rho_0 = (2 - alpha) * M_h / (2 * PI * pow(r_0, alpha) * pow(R, (2 - alpha)));

	Gvar = 6.67408e-11;
	omega_0 = sqrt(2 * PI * Gvar * rho_0 / (2 - alpha));

}

double Model::M(double rr) {
	const double PI = 3.141592653589793238463;
	double r = 2 * PI * rho_0 * pow(r_0, alpha) * pow(rr, (2.0 - alpha)) / (2.0 - alpha);
	return r;
}


void Model::model_loop() {
	

	Cell quad_tree = Cell(-max_radius * 5, max_radius * 5, double(10 * max_radius), nullptr /*cp */,&cells);

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
		rdd_x_i_j[i] = -Gvar * m_sol * get<0>(mdi_j[i]) - (r_x[i] / rr[i]) * Gvar * M(rr[i]) / pow(rr[i], 2);
		rdd_y_i_j[i] = -Gvar * m_sol * get<1>(mdi_j[i]) - (r_y[i] / rr[i]) * Gvar * M(rr[i]) / pow(rr[i], 2);

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


int main() {
	Model model = Model();

	model.init();
	
	ofstream outFile;
	outFile.open("./huts.csv");
	int Iterations = 1000;
	


	cout << "Starting now..." << endl;

	for (int h = 0; h < Iterations; h++) {

		cout << h << "\n";
		
		for (unsigned int j = 0; j < model.r_x.size(); j++) {
			outFile << model.r_x[j] << ";" << model.r_y[j] << ";";
		}
		outFile << endl;

		model.model_loop();
		
	}
	cout << endl << "Done";
	
	outFile.close();

}







