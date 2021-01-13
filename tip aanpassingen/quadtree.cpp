#include "quadtree.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
using namespace std;
Cell::Cell(double top_left_x, double top_left_y, double s, Cell* p,vector<Cell*>* c){
	cells = c;
	x = top_left_x;
	y = top_left_y;
	size = s;
	parent = p;
}

Cell::~Cell() {
}

void Cell::quadruple() {
	double new_size = size / 2;
	children = { NULL,NULL,NULL,NULL };
	double x_new_tl = x;
	double y_new_tl = y;
	double x_new_tr = x_new_tl + new_size;
	double y_new_tr = y_new_tl;
	double x_new_bl = x_new_tl;
	double y_new_bl = y_new_tl - new_size;
	double x_new_br = x_new_tl + new_size;
	double y_new_br = y_new_tl - new_size;
	Cell* cell_tl = new Cell(x_new_tl, y_new_tl, new_size, this,cells);
	Cell* cell_tr = new Cell(x_new_tr, y_new_tr, new_size, this,cells);
	Cell* cell_bl = new Cell(x_new_bl, y_new_bl, new_size, this,cells);
	Cell* cell_br = new Cell(x_new_br, y_new_br, new_size, this,cells);
	cells->push_back(cell_tl);
	cells->push_back(cell_tr);
	cells->push_back(cell_bl);
	cells->push_back(cell_br);

	children = { cell_tl, cell_tr, cell_bl, cell_br };
}

void Cell::find_child_cells(double new_point_x, double new_point_y, double new_point_m) {
	for (Cell* child : children) {
		if (point_x and child->x <= point_x and point_x <= child->x + child->size and child->y - child->size <= point_y and point_y <= child->y) {
			child->insert(point_x,point_y,point_m);
			point_x = NULL;
			point_y = NULL;
			point_m = NULL;
		}
		else if (new_point_x and child->x <= new_point_x and new_point_x <= child->x + child->size and child->y - child->size <= new_point_y and new_point_y <= child->y) {
			child->insert(new_point_x, new_point_y, new_point_m);
			new_point_x = NULL;
			
		}

	}
}

void Cell::insert(double new_point_x, double new_point_y, double new_point_m) {
	if(point_x){
		quadruple();
		find_child_cells(new_point_x, new_point_y, new_point_m);
		
	}
	else if (children[0]) {
		find_child_cells(new_point_x, new_point_y, new_point_m);
		
	}
	else {
		point_x = new_point_x;
		point_y = new_point_y;
		point_m = new_point_m;
	}
}

double Cell::compute_mass() {
	if(x) {
		cm_x = point_x;
		cm_y = point_y;
		total_m = point_m;
		
	}
	else if (children[0]) {
		
		double c1_cm_x = children[0]->point_x;
		double c1_cm_y = children[0]->point_y;
		double c1_m = children[0]->compute_mass();
		double c2_cm_x = children[1]->point_x;
		double c2_cm_y = children[1]->point_y;
		double c2_m = children[1]->compute_mass();
		double c3_cm_x = children[2]->point_x;
		double c3_cm_y = children[2]->point_y;
		double c3_m = children[2]->compute_mass();
		double c4_cm_x = children[3]->point_x;
		double c4_cm_y = children[3]->point_y;
		double c4_m = children[3]->compute_mass();

		total_m = (c1_m + c2_m + c3_m + c4_m);

		cm_x = (c1_cm_x * c1_m + c2_cm_x * c2_m + c3_cm_x * c3_m + c4_cm_x * c4_m) / total_m;
		cm_y = (c1_cm_y * c1_m + c2_cm_y * c2_m + c3_cm_y * c3_m + c4_cm_y * c4_m) / total_m;
	}
	else {
		return 0;
	}

	return total_m;
}

tuple<double, double> Cell::get_mass_distance_term(double probe_x, double probe_y) {
	if (point_x) {
		double d = sqrt(pow(point_x - probe_x, 2) + pow(point_y - probe_y, 2));
		if(d == 0){
			return make_tuple(0, 0);
		}
		double md_x = point_m * (point_x - probe_x) / pow(d, 3);
		double md_y = point_m * (point_y - probe_y) / pow(d, 3);

		return make_tuple(md_x, md_y);

	}
	else if (children[0]) {

		double d = sqrt(pow(cm_x - probe_x,2) + pow(cm_y - probe_y,2));
		if (d == 0) {
			return make_tuple(0, 0);
		}
		if (size / d < THETA) {
			double md_x = total_m * (cm_x - probe_x) / pow(d, 3);
			double md_y = total_m * (cm_y - probe_y) / pow(d, 3);
			return make_tuple(md_x, md_y);
		}
		else {
			double md_x = 0;
			double md_y = 0;
			for (Cell* child : children) {
				double p_md_x, p_md_y;
				tie(p_md_x,p_md_y) = child->get_mass_distance_term(probe_x, probe_y);
				md_x += p_md_x;
				md_y += p_md_y;
			}
			return make_tuple(md_x, md_y);
		}
	}
	else {
		return make_tuple(0, 0);
	}
}

void Cell::draw() {

}