//
//  treenode.hpp
//  defines the treenode class as well as bodynode and cellnode subclasses
//
//  Created by Noah Muldavin
//  Reed College
//  Spring 2013
/******************************************************************************/


#include <iostream>
#include "Vector.hpp"


#ifndef Classes_treenode_h
#define Classes_treenode_h


//  defines a general treenode class

class treenode
{
public:
    treenode();         //  Constructor
    bool isbody;        //  "true" if node is a bodynode
    double d;           //  "size" of treenode
    double mass;        //  Total mass of particles
    Vector pos;         //  Position of mass or center of mass
};


//  Defines bodynode subclass. This is where the information about
//  particles will be stored. The bodynode class inherits all of the
//  elements of the treenode class plus a few extra.

class bodynode : public treenode
{
public:
    bodynode();         //  Constructor
    double t;           //  current time
    double dt;          //  timestep
    double tp;          //  prediction time
    Vector vel;         //  Velocity of particle
    unsigned int id;    //  Particle id
    bool DM;
};


//  Defines a cellnode subclass. Each cellnode has links to eight
//  child cellnodes.

class cellnode : public treenode
{
public:
    cellnode();         //  Constructor
    void clear();
    treenode* child[8]; //  Pointer to child nodes. 
};


#endif


/******************************************************************************/
