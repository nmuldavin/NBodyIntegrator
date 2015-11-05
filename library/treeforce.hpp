//
//  treeforce.cpp
//  header for force calculation given a tree datastructure
//
//  Created by Noah Muldavin
//  Reed College
//  Spring 2013
/******************************************************************************/


#include "maketree.hpp"

#ifndef treeforce_hpp
#define treeforce_hpp

Vector gravity(double mass, Vector rij, double eps);
Vector treeforce(bodynode* body, treenode* root, double eps, double theta);

#endif

/******************************************************************************/