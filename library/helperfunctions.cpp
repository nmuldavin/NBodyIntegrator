//
//	main.cpp
//	Main N-Body file. Important parameters may be edited here. 
//
//	Created by Noah Muldavin
//	Reed College
//	2013
/******************************************************************************/


#include "helperfunctions.hpp"
#include <stdlib.h>
#include <cmath>

const double PI = 3.141592653589;

//  generates a random number between min and max with a uniform distribution

double randomreal(double min, double max)
{
    double random = ((double) rand()) / (double) RAND_MAX;  // returns a
                                                //  random number in [0,1] 
    double diff = max - min;                    //  expands to [0, max-min]
    double r = random * diff;
        
    return min + r;                             //  adjusts to [min, max]
        
}


//  generates a gaussian random number with a box-muller transform

double gaussianrandom(double sigma)
{
    double x1 = ((double) rand()) / (double) RAND_MAX;
    double x2 = ((double) rand()) / (double) RAND_MAX;
    
    double out  = sqrt(-2.0 * sigma * log(x1)) * cos (2.0*PI*x2);
    
    return out;
}


//  rejectionsample will generate a random number according to a specified
//  probability distribution. The probability distribution must be a
//  funcion of one variable
double rejectionsample(double (*f)(double x), double xmin,
                        double xmax, double topbound)
{
    while(true)
    {
        double xr = randomreal(xmin, xmax);
        double yr = randomreal(0.0, topbound);
        if (yr <= f(xr))                        //  if y random # is greater 
        {                                       //  than height of the function
            return xr;                          //  at that point, throw it out
            break;                              //  and repeat, otherwise
        }                                       //  return x random #
    }
}


//  will find the maximum value of a one dimensional function between
//  xmin and xmax
double findtopbound(double (*f)(double x), double xmin,
                    double xmax)
{
    int numsteps = 1000000;
    double ymax = f(xmin);
    double xr, y;
    for(int i = 0; i < numsteps; i++)
    {
        xr = xmin + (xmax - xmin)*(double)i / (double)numsteps;
        y = f(xr);
  
        if (y > ymax)
        {
            ymax = y;                       //  if value exceeds current max,
        }                                   //  reset
    }
    return ymax;
}


//  randomizes position on a sphere of radius r. Returns a vector in
//  cartesian coordinates
Vector randomsphere(double radius)
{
	double inc = PI*(3.0 - sqrt(5.0));
	double off = 2.0/(double)(RAND_MAX);
	double y, r, phi;
    double pos = (double)rand();
    y = pos*off - 1.0 + 0.5*off;
    phi = pos*inc;
    r = sqrt(1.0 - y*y);
    Vector out(3);
    out[0] = radius*r*cos(phi);
    out[1] = radius*y;
    out[2] = radius*r*sin(phi);
    return out;
}


//  transforms cylindrical vector to cartesian vector
Vector cyltocart(Vector cyl)
{
    Vector cart(3);
    cart[0] = cyl[0] * cos(cyl[1]);
    cart[1] = cyl[0] * sin(cyl[1]);
    cart[2] = cyl[2];
    return cart;
}


//  transforms cartesian vector to cylindrical vector
Vector carttocyl(Vector cart)
{
    Vector cyl(3);
    cyl[0] = sqrt(cart[0]*cart[0]+cart[1]*cart[1]);
    cyl[1] = atan(cart[1]/cart[0]);
    cyl[2] = cart[2];
    return cyl;
}


//  transforms spherical vector to cartesian vector
Vector spheretocart(Vector sphere)
{
    Vector cart(3);
    cart[0] = sphere[0]*sin(sphere[1])*cos(sphere[2]);
    cart[1] = sphere[0]*sin(sphere[1])*sin(sphere[2]);
    cart[2] = sphere[0]*cos(sphere[1]);
    return cart;
}


//  transforms cartesian vector to spherical vector
Vector carttosphere(Vector cart)
{
    Vector sphere(3);
    sphere[0] = sqrt(cart[0]*cart[0] + cart[1]*cart[1] + cart[2]*cart[2]);
    sphere[1] = acos(cart[2]/sphere[0]);
    sphere[2] = atan(cart[1]/cart[0]);
    return sphere;
}


//  finds the epicylic frequency at the location of a bodynode with index i.
//  since this uses the treecode force calculation, a tree must be constructed
//  at "root" before use. 
double epicyclicfrequency(bodynode* bodies, int N, int i, cellnode* root,
                          double eps, double theta, double dr)
{
    
    double r = sqrt(bodies[i].pos[0]*bodies[i].pos[0]
                        + bodies[i].pos[1]*bodies[i].pos[1]);
    double frad, frad2;
    double k2;
    bodynode testbody;
    Vector force(3);
    Vector force2(3);
    force = treeforce(&bodies[i], root, eps, theta);    //  acceleration
    frad = (force[0]*bodies[i].pos[0] + force[1]*bodies[i].pos[1])/r;//
                                            //  radial force
    
    testbody.id = bodies[i].id;             //  setting testbody id
    testbody.pos = bodies[i].pos*(1.0 + dr);//  position is slightly out
    force2 = treeforce(&testbody, root, eps, theta);    //  new force
    frad2 = (force2[0]*bodies[i].pos[0] + force2[1]*bodies[i].pos[1])/r;
    k2 = (3.0 / r)*(frad) + (frad2 - frad)/dr;  // epicyclic frequency squared
    return sqrt(fabs(k2));
}
