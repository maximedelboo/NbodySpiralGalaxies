#ifndef SOURCE_H
#define SOURCE_H
#include <vector>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <vector>
#include "quadtree.h"

using namespace std;
class Model {
public:
	vector<Cell*> cells;
	const double PI = 3.141592653589793238463;
	int n = 10000;
	double radius_gal = 1e18;
	double m_center_mass = 700;
	double pert_distance = 1e18;
	vector<double> radii;
	vector<double> r_x;
	vector<double> r_y;
	vector<double> rd_x;
	vector<double> rd_y;
	vector<double> rdd_x;
	vector<double> rdd_y;
	vector<double> m;
	

	vector<tuple<double, double>> mdi_j;
	vector<double> rr;
	vector<double> rd_x_i_j;
	vector<double> rd_y_i_j;
	vector<double> r_x_i_jp1;
	vector<double> r_y_i_jp1;
	vector<double> rdd_x_i_j;
	vector<double> rdd_y_i_j;

	double Gvar = 6.67408e-11;
	double m_sol = 1.98847e30;
	double dt = 35e12;
	double f = 0.000001;
	double alpha = 1.5;
	double M_d = 0;
	double M_h = 0;
	double R = 0;
	double r_0 = 0;
	double rho_0 = 0;
	double omega_0 = 0;
	double epsilon = radius_gal/10;
    double cte = m_center_mass * m_sol / 20 / 100000000000;
    double omega = sqrt(Gvar*m_center_mass*m_sol/ pow((radius_gal / 2),3));
    double t = 0;
    double T = 0;
    vector<double> AX;
    vector<double> AY;


	void init();
	double M(double rr);
	void model_loop();
};
#endif