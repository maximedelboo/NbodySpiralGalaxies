#ifndef QUADTREE_H
#define QUADTREE_H
#include <vector>
#include <functional>
#include <numeric>
#include <math.h>
#include <random>
#include <tuple>

#define THETA 0.5

using namespace std;
class Cell {
public:
	vector<Cell*>* cells;
	double x;
	double y;
	double size;
	Cell* parent;
	vector<Cell*> children = { NULL,NULL,NULL,NULL};
	double point_x = 0;
	double point_y = 0;
	double point_m = 0;
	double cm_x = 0;
	double cm_y = 0;
	double total_m = 0;

	Cell(double top_left_x, double top_left_y, double s, Cell* p,vector<Cell*>* c);
	~Cell();

	void quadruple();

	void find_child_cells(double new_point_x, double new_point_y, double new_point_m);

	void insert(double new_point_x, double new_point_y, double new_point_m);

	double compute_mass();

	tuple<double,double> get_mass_distance_term(double probe_x, double probe_y);

	void draw();

};
#endif