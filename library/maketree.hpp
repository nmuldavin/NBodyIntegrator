//
//  maketree.hpp
//  header containing routines for creation and organization of the tree
//
//  Created by Noah Muldavin for Physics 367: Scientific Computing
//  Reed College
//  Spring 2013
/******************************************************************************/


#include "treenode.hpp"
#include <vector>
#include <ctime>

#ifndef maketree_hpp
#define maketree_hpp

cellnode* makecell();
void collectcells(treenode* node);
void makeboundingvol(bodynode* bodies, int N, cellnode* root);
int childindex(bodynode body, cellnode* cell);
Vector offsetdir(int index);
void insertbody(bodynode body, cellnode cell);
void setcoms(cellnode* root);
void maketree(bodynode* bodies, int N, cellnode* root);

#endif


/******************************************************************************/