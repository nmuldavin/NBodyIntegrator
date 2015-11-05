//
//  treeforce.cpp
//  methods for force calculation given a tree datastructure
//
//  Created by Noah Muldavin
//  Reed College
//  2013
/******************************************************************************/


#include "treeforce.hpp"


//  gravity returns the force on an object i from another object j given its
//  mass and a vector pointing from i to j. Includes the softening length eps

Vector gravity(double mass, Vector rij, double eps)
{
    Vector force(3);                        //  force variable
    double d2;                              //  distance squared
    d2 = rij*rij + eps*eps;                 //  setting distance squared
    force = -rij*mass/pow(d2, 1.5);         //  setting force
    return force;                           //  output
}


//  treeforce returns the force on a body due to the collection of bodies
//  within a tree located at root.

Vector treeforce(bodynode* body, treenode* root, double eps, double theta)
{
    Vector rij(3);                          //  separation vector
    Vector force(3);                        //  force vector
    double r;                               //  separation magnitude
    force = 0.0;                            //  zeros force vector
    rij = body->pos - root->pos;            //  sets separation vector
    if (!root->isbody)                      //  if the node is a cell
    {
        r = rij.norm();                     //  sets separation magnitude
        if ((2.0*root->d)/r < theta)        //  if acceptable
        {
            force += gravity(root->mass, rij, eps);  //  add force from COM
        }
        else                                //  if not acceptible
        {
            for (int i = 0; i < 8; i++)     //  add force from each child node
            {
                if (((cellnode*) root)->child[i] != 0)  //  if it exists
                {
                    force += treeforce(body, ((cellnode*) root)->child[i],
                                                                    eps, theta);
                }
            }
        }
    }
    else                                    //  if the node is a bodynode
    {
        if (body->id != ((bodynode*) root)->id) //  add force from body
        {
            force += gravity(root->mass, rij, eps);
        }
    }
    return force;
}


/******************************************************************************/