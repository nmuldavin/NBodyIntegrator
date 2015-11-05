//
//	main.cpp
//	Main N-Body file. Important parameters may be edited here. 
//
//	Created by Noah Muldavin
//	Reed College
//	2013
/******************************************************************************/


#include "library/dataio.hpp"


//////////////////////////////////////////////	IMPORTANT PARAMETERS
std::string startfile = "data/initfile.txt";//	location of startfile
std::string datadir = "data/";				//	directory where data will go
double dtmax = 0.0787;						//	maximum timestep
double tfinal = 1000000;						//	final time in units of dtmax
double alpha = 0.1;							//	timestep tolerance param.
double theta = 1.0;							//	MAC theta
double renderbox = 10.0;					//	size of rendering box
bool plotDM = true;						//	plot dark matter particles
//////////////////////////////////////////////



// other variables - DONT EDIT 				
int N;										//	number of particles
int NDM;									//	number of dark particles
double eps;									//	softening length
double t = 0.0;								//	current time
int writecount = 0;							//	count of IFrIT out files
int iter = 1;								//	iteration count


//	data locations - DONT EDIT				
bodynode* bodies;							//	array of bodynodes
cellnode* root = new cellnode;				//	root cell



int main()
{
	bodies = readinitfile(N, NDM, startfile);	//	reads from startfile
	std::cout << "N = " << N << "\n";			//	reports N
	std::cout << "NDM = " << NDM << "\n";		//	reports NDM
	
	eps = 0.98*pow((double)N, -0.26);			//	setting softening length
	std::cout << "Softening Length epsilon = " << eps << "\n";
	
	initialize(bodies, N, root, dtmax, eps, alpha, theta);	//	initializing
	
	
	while (t/dtmax < tfinal)					//	while time less than final
	{											//	reports iteration
		std::cout << "Iteration " << iter << ", t = " << t/dtmax << ": ";
		
		advance(bodies, N, root, eps, theta, alpha, dtmax, t);	//	advance
		
		
		if (t/dtmax >= (double)writecount)		//	if time advanced past dtmax
		{										//	write IFrIT file
			writeifritfile(N, NDM, writecount, bodies,
								renderbox, plotDM, datadir);
			writecount += 1;					//	increase writecount
		}
		iter += 1;								//	next iteration
	}
												//	write fullsave file
	writefullsave(N, NDM, bodies, datadir+"finalsave.txt");	
}



/******************************************************************************/
