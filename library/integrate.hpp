//
//  integrate.hpp
//  Header containing integration routines for N-Body simulation
//
//  Created by Noah Muldavin
//  Reed College
//  2013
/******************************************************************************/


#include "treeforce.hpp"

#ifndef integrate_hpp
#define integrate_hpp

double tpmin(bodynode* bodies, int N);
double dtideal(Vector force, double eps, double alpha, double dtmax);
int roundtoint(double roundme);
void setmidpos(bodynode* bodies, int N, double tpm);
void leapfrogvar(bodynode* body, double eps,
                    cellnode* root, double theta, double alpha, double dtmax);
void advance(bodynode* bodies, int N, cellnode* root,
                double eps, double theta, double alpha, double dtmax, double& t);

#endif


/******************************************************************************/