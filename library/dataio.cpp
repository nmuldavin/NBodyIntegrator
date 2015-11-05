//
//  dataio.cpp
//  Defines routines for saving data + communicating with the system
//
//  Created by Noah Muldavin for Physics 367: Scientific Computing
//  Reed College
//  Spring 2013
/******************************************************************************/


#include "dataio.hpp"
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <cassert>


//  writes a startfile at "filename" which includes the number of particles N,
//  the number which are dark matter particles, and for each particle the
//  position, velocity, and boolean value indicating whether it is a dark matter
//  particle or not. 

void writeinitfile(int N, int NDM, bodynode* bodies, std::string filename)
{
    std::ofstream output(filename.c_str()); //  opens file using "output" stream
    assert(output.is_open());               //  checks that it's open
    
    output << N << " " << NDM << "\n";      //  writes N and NDM on first line
    
    output.setf(std::ios::scientific);      //  sets scientific notation output
    output.setf(std::ios::showpos);         //  show +/- sign
    output.precision(7);                    //  set precision to 7 decimals
    
    for (int i = 0; i < N; i++)             //  for each particle
    {
        output << bodies[i].mass << " ";
        output << bodies[i].pos << " " << bodies[i].vel
            << " " << bodies[i].DM << "\n"; //  write position, velocity, and 
    }                                       //  boolean DMs
    
    output.close();                         //  closes output stream
    std::cout << "Wrote file " << filename << "\n";
}


//  reads a startfile located at "filename" in the form outputted by
//  "writeinitfile" above. Returns an allocated array of bodynodes, and assigns
//  N and NDM appropriately.

bodynode* readinitfile(int& N, int& NDM, std::string filename)
{
    std::cout << "Reading from file " << filename << "\n";
    std::ifstream input(filename.c_str());  //  opens file in "input" stream
    assert(input.is_open());                //  checks that it's open
    input >> N >> NDM;                      //  assigns N and NDM
    bodynode* bodies = new bodynode[N];     //  allocates necessary memory
    for (int i = 0; i < N; i++)             //  for each particle
    {
        input >> bodies[i].mass;
        for (int j = 0; j < 3; j++)
        {
            input >> bodies[i].pos[j];      //  reads position
        }
        for (int j = 0; j < 3; j++)
        {
            input >> bodies[i].vel[j];      //  reads velocity
        }
        input >> bodies[i].DM;              //  reads DM
        bodies[i].id = i;                   //  gives particle an ID
    }
    input.close();                          //  closes input stream
    return bodies;
}


//  initializes the timesteps for a group of particles. Designed to be used
//  after "readinitfile" to initialize an integration.

void initialize(bodynode* bodies, int N, cellnode* root,
					double dtmax, double eps, double alpha, double theta)
{
	std::cout << "Initializing " << N << " particles.\n";
	maketree(bodies, N, root);              //  makes tree for all bodies
	Vector force(3);                        //  force variable
	double dtprov;                          //  provisional timestep variable
	for (int i = 0; i < N; i++)             //  for each particle
	{
		force = treeforce(&bodies[i], root, eps, theta);    //  finds force
		bodies[i].dt = dtideal(force, eps, alpha, 10e10);   //  sets ideal dt
		dtprov = dtmax;                     //  sets provisional timetsep to max
		while (bodies[i].dt < dtprov)       //  while ideal < provisional
		{
			dtprov *= 0.5;                  //  divides by two
		}
		bodies[i].dt = dtprov;              //  sets timestep to final value
		bodies[i].tp = 0.5*dtprov;          //  prediction time is half timestep
	}
}


//  takes an integer input at returns a string containing that number, as well as
//  '0's as fillers if there are less than 4 digits. This will be used to write
//  file names in "writeifritfile()"

std::string inttostring (int x)
{
        std::stringstream xstream; 
            //allows string "xstream" to be manipulated as an output
        xstream << std::setw(4) << std::setfill( '0' ) << x;
            //sets character # to 4, fills in 0s, imports x into xstream variable
        return xstream.str();
}


//  writes a file to be used by the ifrit visualization program

void writeifritfile(int N, int NDM, int step, bodynode* bodies,
                        double boxsize, bool plotDM, std::string datadir)
{
    std::string sstep = inttostring (step); //  makes string containing file #
    std::string filetype = ".txt";          //  file will be .dat
    sstep+=filetype;                        //  adds file extension
    sstep = datadir + sstep;
    std::ofstream output(sstep.c_str());    //  opens file
    assert(output.is_open());               //  makes sure it's open
    
    if(plotDM)                              //  if plotting DM
    {
        output << N << "\n";                //  write total # of particles
    }
    else
    {
        output << (N - NDM) << "\n";        //  otherwise write number of non-dark
    }
    
    
    output.setf(std::ios::showpoint);       //  show decimal
    output.precision(2);                    //  two decimals
                                            //  writing box size
    output << -boxsize << " " << -boxsize << " " << -boxsize << " "
            << boxsize << " " << boxsize << " " << boxsize << "\n";
    
    output.setf(std::ios::scientific);      //  sets scientific notation output
    output.setf(std::ios::showpos);         //  show +/- sign
    output.precision(7);                    //  set precision to 7 decimals
    
    if(plotDM)                              //  if plotting dark matter
    {
        for (int i = 0; i < N; i++)         //  writes position and DM boolean
        {
            output << bodies[i].pos << " " << bodies[i].DM << "\n";
        }
    }
    else                                    //  if not
    {
        for (int i = 0; i < N; i++)         //  writes position of all non dark
        {
            if(!bodies[i].DM)
            {
                output << bodies[i].pos << "\n";
            }
        }
    }
    
    
    output.close();                         //  closes file
    
    std::cout << "Wrote ifrit file " << sstep << "\n";
}


//  writes a full savefile containing all data for the array of N particles

void writefullsave(int N, int NDM, bodynode* bodies, std::string filename)
{
    std::ofstream output(filename.c_str()); //  opens file
    assert(output.is_open());               //  makes sure it's open
    
    output << N << " " << NDM << "\n";      //  writes N and NDM
    
    output.setf(std::ios::scientific);      //  sets scientific notation output
    output.setf(std::ios::showpos);         //  show +/- sign
    output.precision(7);                    //  set precision to 7 decimals
    
    for (int i = 0; i < N; i++)             //  for each particle, writes all
    {                                       //  data
        output << bodies[i].mass << " ";
        output << bodies[i].pos << " " << bodies[i].vel << " "
                << bodies[i].t << " " << bodies[i].tp << " "
                << bodies[i].dt << " " << bodies[i].DM << "\n";
    }
    
    output.close();                         //  closes file
    std::cout << "Wrote savefile " << filename << "\n";
}


//  reads a full savefile in the form outputted above. Returns a pointer to
//  an allocated array of bodynodes, and updates N and NDM

bodynode* readfullsave(int& N, int& NDM, std::string filename)
{
    std::cout << "Reading from savefile " << filename << "\n";
    std::ifstream input(filename.c_str());  //  opens file
    assert(input.is_open());                //  makes sure it's open
    input >> N >> NDM;                      //  read N and NDM
    bodynode* bodies = new bodynode[N];
    for (int i = 0; i < N; i++)             //  for each particle
    {
        input >> bodies[i].mass;
        for (int j = 0; j < 3; j++)         //  read position
        {
            input >> bodies[i].pos[j];
        }
        for (int j = 0; j < 3; j++)         //  velocity
        {
            input >> bodies[i].vel[j];
        }
        input >> bodies[i].t;               //  time
        input >> bodies[i].tp;              //  prediction time
        input >> bodies[i].dt;              //  timestep
        input >> bodies[i].DM;              //  Dark Matter boolean
    }
    input.close();                          //  closes file
    return bodies;
}


/******************************************************************************/