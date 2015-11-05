//
//	main.cpp
//	Main N-Body file. Important parameters may be edited here. 
//
//	Created by Noah Muldavin
//	Reed College
//	2013
/******************************************************************************/


#include "dataio.hpp"
#include "Vector.hpp"
#include "maketree.hpp"


#ifndef helperfunctions_headers
#define helperfunctions_headers


double randomreal(double min, double max);
double findtopbound(double (*f)(double x), double xmin, double xmax);
double rejectionsample(double (*f)(double x), double xmin, double xmax, double topbound);
double gaussianrandom(double sigma);
Vector randomsphere(double radius);
Vector cyltocart(Vector cyl);
Vector carttocyl(Vector cart);
Vector spheretocart(Vector sphere);
Vector carttosphere(Vector cart);
double epicyclicfrequency(bodynode* bodies, int N, int i, cellnode* root, double eps, double theta, double dr);


#endif