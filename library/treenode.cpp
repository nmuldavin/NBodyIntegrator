//
//  treenode.cpp
//  Constructors for treenode class and subclasses
//
//  Created by Noah Muldavin
//  Reed College
//  2013
/******************************************************************************/


#include "treenode.hpp"


//  Default constructor for treenode class

treenode::treenode() : pos(3)
{
    isbody = false;
    d = 0.0;
    mass = 0.0;
}


//  Default constructor for bodynode subclass

bodynode::bodynode() : treenode(), vel(3)
{
    isbody = true;
    DM = false;
    t = 0.0;
    dt = 0.0;
    tp = 0.0;
    mass = 1.0;
    id = 0;
    d = 0.0;
}


//  Default constructor for cellnode subclass

cellnode::cellnode() : treenode()
{
    isbody = false;
    for (int i = 0; i < 8; i++)
    {
        child[i] = 0;
    }
}


//

void cellnode::clear()
{
    for (int i = 0; i < 8; i++)
    {
        child[i] = 0;
    }
}


/******************************************************************************/