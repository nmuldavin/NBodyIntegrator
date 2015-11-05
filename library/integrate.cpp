//
//  integrate.cpp
//  Contains integration routines for N-Body simulation
//
//  Created by Noah Muldavin
//  Reed College
//  2013
/******************************************************************************/


#include "integrate.hpp"


//	calculates the ideal timestep given the force from the previous integration,
//	the softening length and the tolerance parameter

double dtideal(Vector force, double eps, double alpha, double dtmax)
{
    double tpi;
    tpi = sqrt(alpha*eps/force.norm());		//	ideal timestep
	if(tpi > dtmax)
	{
		tpi = dtmax;						//	caps at maximum
	}
    return tpi;
}


//	finds the minimum prediction time in a list of N particles

double tpmin(bodynode* bodies, int N)
{
    double min = bodies[0].tp;
    for (int i = 1; i < N; i++)				//	for all particles
    {
        if (bodies[i].tp < min)				//	if tp < current min
        {
            min = bodies[i].tp;				//	resets min
        }
    }
    return min;
}


//	Advances all particles to the time at which the velocity of the next group
//	of particles will be advanced.

void setmidpos(bodynode* bodies, int N, double tpm)
{
    for (int i = 0; i < N; i++)				//	moves all particles to tpmin
    {
        bodies[i].pos += bodies[i].vel*(tpm - bodies[i].t);
        bodies[i].t = tpm;					
    }
}


//	Takes a double precision floating point number and rounds it to an integer

int roundtoint(double roundme)				
{
    return floor(roundme+0.5);				//	returns integer
}


//	Advances a bodynode in time given the force calculated from the group of
//	particles stored in a tree with a given root.

void leapfrogvar(bodynode* body, double eps, cellnode* root,
					double theta, double alpha, double dtmax)
{
    Vector force(3);
    double dtp;
    force = treeforce(body, root, eps, theta);	//	finds force with treeforce()
    body->vel += force*body->dt;				//	updates velocity
    body->pos += body->vel*body->dt/2.0;		//	updates position
    body->t = body->tp + body->dt/2.0;			//	updates time
    dtp = dtideal(force, eps, alpha, dtmax);	//	finds ideal timestep
	if (dtp > dtmax)
	{
		dtp = dtmax;
	}
    if (dtp < body->dt)							//	if ideal < previous
    {
        body->dt /= 2.0;						//	divides previous by two
    }
    else if ((dtp > 2.0*body->dt) &&			//	if ideal > 2*previous
             (roundtoint(body->t/body->dt) % 2 == 0))	//	and t/dt is even
	{
        body->dt *= 2.0;					//	new is double old
    }
    body->tp = body->t + body->dt/2.0;		//	setting new prediction time
}


//	advances the entire collection of N bodies, using tree construction
//	and integration routines

void advance(bodynode* bodies, int N, cellnode* root, double eps,
				double theta, double alpha, double dtmax, double& t)
{							
	clock_t forcetime;						//	to note total time
    int count = 0;
	forcetime = clock();					//	notes time
    t = tpmin(bodies, N);					//	minimum prediction time
    setmidpos(bodies, N, t);				//	advances all masses to tpmin
	maketree(bodies, N, root);				//	builds the tree
    for (int i = 0; i < N; i++)
    {
        if (fabs(t - bodies[i].tp) < 10.0e-10)	// if within rounding error-
											//	bound of tpmin
        {
            leapfrogvar(&bodies[i], eps, root, theta, alpha, dtmax); //	advances
            count += 1;						//	adds to count
        }
    }
	forcetime = clock() - forcetime;		//	notes total time and prints: 
    std::cout << "Advanced " << count << " particles in "
				<< (double) forcetime / (double) CLOCKS_PER_SEC << " seconds\n";
}


/******************************************************************************/