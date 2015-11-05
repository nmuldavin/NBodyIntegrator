//
//  maketree.cpp
//  defines routines for the creation and organizatoin of the tree
//
//  Created by Noah Muldavin
//  Reed College
//  2013
/******************************************************************************/


#include "maketree.hpp"


//  global variables:

std::vector<treenode*> freecells;           //  Stores a list of pointers to
                                            //  all free cells

Vector offset(3);                           //  a unit vector which containing
                                            //  the direction from the position
                                            //  of a cell node to one of its
                                            //  child nodes


//  makecell returns a pointer to a free, allocated cellnode
//  replaces memory allocation with "new" command such that old cells are
//  reused

cellnode* makecell()
{
    cellnode* cell;
    if (freecells.empty())                  //  if there are no free cells,
    {                                       //  allocates a new one
        cell = new cellnode;
    }
    else                                    //  otherwise takes the last free
    {                                       //  cell
        cell = (cellnode*) freecells.back();
        freecells.pop_back();
    }
    cell->clear();                          //  clearing mass and linkage data
    return cell;
}


//  collectcells scans the tree starting at a root node, collecting all
//  cellnodes into the freecells list

void collectcells(treenode* node)
{
    cellnode* p;
    if (!node->isbody)                      //  if the node is a cell node
    {
        freecells.push_back(node);          //  adds it to the freecells list
        p = (cellnode*) node;
        for (int i = 0; i < 8; i++)         //  then repeats scan for each
        {                                   //  unempty child node
           if(p->child[i] != 0)
           {
                collectcells(p->child[i]);
           }
        }
    }
}


//  makeboundingvol sets the box size of the rootnode box so that it
//  contains all particles.

void makeboundingvol(bodynode* bodies, int N, cellnode* root)
{
    double d = 0.0;
    double dmax = 0.0;
    for (int i = 0; i < N; i++)             //  for all particles
    {
        for(int j = 0; j < 3; j++)          //  if coordinate is greater 
        {                                   //  than max
            d = fabs(bodies[i].pos[j] - root->pos[j]);
            if (d > dmax)
            {
                dmax = d;                   //  resets max
            }
        }
    }
    root->d = 1.001*dmax;                      //  updates root size
}


//  childindex finds the index of the child node containing a given bodynode.

int childindex(bodynode* body, cellnode* cell)
{
    int index = 0;
    offset = -1.0;
    for (int i = 0; i < 3; i++)
    {;
        if (cell->pos[i] <= body->pos[i])
        {
            index += pow(2, 2-i);
            offset[i] = 1.0;
        }
    }
    return index;
}


//  insertbody places a body into its appropriate location in a cell by
//  scanning down the tree until a free childnode is reached. 

void insertbody(bodynode* body, cellnode* cell)
{
    int index;                              
    int oldbodyindex;
    index = childindex(body, cell);         //  calculates index of child
                                            //  which contains the body
    if(cell->child[index] == 0)             //  if child cell is free
    {
        cell->child[index] = (treenode*) body;  //  place body
    }
    else if (cell->child[index]->isbody)    //  if child is a body
    {
        cellnode* newcell; 
        newcell = makecell();               //  creates a new cell
        newcell->d = cell->d/2.0;           //  new cell is half as big
        newcell->pos = cell->pos + offset*(newcell->d);     //  newcell
                                            //  position offset by d/2 in
                                            //  the appropriate direction
        oldbodyindex = childindex((bodynode*) cell->child[index], newcell); //
                                            //  finds the index of the old body
                                            //  in the new cell
        newcell->child[oldbodyindex] = cell->child[index];  //  inserts old body
                                            //  in new cell
        cell->child[index] = (treenode*) newcell;   //  links newcell to
                                            //  the original cellnode
        insertbody(body, newcell);
    }
    else                                    //  if child is a cellnode
    {                                       //  inserts body in that cell
        insertbody(body, (cellnode*) cell->child[index]);   
    }
}


//  setcoms sets the mass and center of mass of each cellnode recursively, by
//  first setting the mass and center of mass of each child node, then summing
//  appropriately.

void setcoms(cellnode* cell)
{
    cell->pos = 0.0;                        //  zeros cell position
    cell->mass = 0.0;                       //  zeros cell mass
    for(int i = 0; i < 8; i++)              //  for each child node
    {
        if (cell->child[i] != 0)            //  if child exists
        {
            if (!cell->child[i]->isbody)    //  if it is a cellnode
            {
                setcoms((cellnode*) cell->child[i]);    //  setcoms for child
            }
            cell->mass += cell->child[i]->mass; //  updates mass
            cell->pos += (cell->child[i]->pos)*cell->child[i]->mass;    //
                                            //  updates position
        }
    }
    cell->pos /= cell->mass;                //  divides by total mass
}


//  maketree is the main tree construction function which includes all above
//  routines. It takes as aguments a table of bodynodes, the number of bodynodes
//  and a pointer to a rootnode which will contain the cells.

void maketree(bodynode* bodies, int N, cellnode* root)
{
    clock_t buildtime = clock();            //  notes the time
    for(int i = 0; i < 8; i++)              //  collects all allocated cells
    {                                       //  starting one level below the
        if (root->child[i] != 0)            //  root
        {
            collectcells(root->child[i]);
        }
    }
    root->clear();                          //  clears root linkage
    root->pos = 0.0;                        //  zeros rootnode position
    makeboundingvol(bodies, N, root);
    for (int i = 0; i < N; i++)             //  inserts each body in root by
    {
        insertbody(&bodies[i], root);
    }
    setcoms(root);                          //  sets COMs for all cellnodes
    buildtime = clock() - buildtime;        //  notes the total time and prints
}


/******************************************************************************/