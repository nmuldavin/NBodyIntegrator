//
//  dataio.cpp
//  Header defining routines for communicating with the system
//
//  Created by Noah Muldavin for Physics 367: Scientific Computing
//  Reed College
//  Spring 2013
/******************************************************************************/


#include "integrate.hpp"


#ifndef dataio_hpp
#define dataio_Vector_hpp


void writeinitfile(int N, int NDM, bodynode* bodies, std::string fileloc);
bodynode* readinitfile(int& N, int& NDM, std::string fileloc);
void initialize (bodynode* bodies, int N, cellnode* root,
					double dtmax, double eps, double alpha, double theta);
void writeifritfile(int N, int NDM, int step, bodynode* bodies,
                        double boxsize, bool plotDM, std::string datadir);
std::string inttostring (int x);
void writefullsave(int N, int NDM, bodynode* bodies, std::string fileloc);
bodynode* readfullsave(int& N, int& NDM, std::string fileloc);


#endif


/******************************************************************************/