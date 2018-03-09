#ifndef _PROTOS_H
#define _PROTOS_H

#include "defines.h"

int read_patterns(FILE*, double**, int**, int*, int*, int);
double sqdist(double*, double*, int);
void compute_mean(treenode*, double*, int);
void compute_radius(treenode*, double*, int);
void dump_cluster(treenode*, double*, int);
void three_means(treenode*, double*, int, double*, int);
treenode* tree_start(double*, int, int);
void tree_build(int, treenode*, double*, int, double*, int);
double search_tree(double*, int, treenode*, double*, int, double*, int, int*,
		   double*, double, int*);
int predict(double*, int, double*, int, int*, int, int*, treenode*,
	    double*, int, int, int*, int*);
void usage(char*, char*);

#endif
