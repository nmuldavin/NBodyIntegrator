//
//	main.cpp
//	Main N-Body file. Important parameters may be edited here. 
//
//	Created by Noah Muldavin
//	Reed College
//	2013
/******************************************************************************/

#include "library/helperfunctions.hpp"
#include <cmath>
#include <cassert>
#include <fstream>
#include <vector>
#include <algorithm>

const double PI = 3.141592653589;           //  Pi, dont change Pi



//////////////////////////////////////////////	DISK PARAMETERS
int Ndisk = 40000;                         //  Number of particles
const double Md = 1.0;                      //  Mass of disk
const double h = 1.0;                       //  radial scale length
const double z0 = 0.2;                      //  vertical scale length
const double diskcut = 7.5;                 //  radial cutoff
const double thickcut = 3.0*z0;             //  vertical cutoff
                                            //  normalization constant (set
                                            //  later)
double Q = 1.2;                             //  Toomre Q parameter
double rref = 2.5*h;                        //  reference radius at which to
                                            //  normalize radial dispersion
double surfacedist (double r)               //  surface probability distribution
{                                           //  to be sampled
    return exp(-r / h)*r;                   //  
}                                           //
                                            //
double Sigma(double r)                      //  normalized surface density 
{                                           //  to use when calculating Q
    return Md / (2.0*PI*h*h)*exp(-r / h);//
}                                           //
                                            //
double diskthick (double z)                 //  vertical mass distribution
{                                           //
    return 1.0/pow(cosh(z / z0), 2.0);      //
}                                           //
                                            //
//////////////////////////////////////////////


//////////////////////////////////////////////	HALO PARAMETERS
int Nhalo = 120000;                         //  Number of particles
const double Mh = 5.0;                      //  Halo Mass
const double rc = 10.0;                     //  scale length
const double g = 2.0;                       //  another scale length
const double halocut = 15.0;                //  cutoff radius
double vr2;                                 //  radial velocity dispersion
                                            //
double hernquisthalo (double r)             //  halo probability distribution
{                                           //  as a function of radius
    return exp(-pow(r/rc, 2.0))/(r*r + g*g)*r*r;    //
}                                           //  
                                            //
double halospeeddist(double v)              //  distribution of absolute speed
{                                           //
    return v*v*exp(-v*v/2.0/vr2);           //
}                                           //
                                            //
//////////////////////////////////////////////


//////////////////////////////////////////////	BULGE PARAMETERS
int Nbulge = 20000;                        //  Number of particles
const double Mb = 0.6;                      //  Bulge Mass
const double a = 0.4;                       //  scale length
const double bulgecut = 5.0;                //  cutoff radius
                                            //
double bulgedist(double r)                  //  bulge probability distribution
{                                           //  as a function of radius in 
    return r / (a*a*pow(1.0 + r/a, 3.0));   //  spherical coordinates
}                                           //
                                            //
//////////////////////////////////////////////


int Ntot = Ndisk + Nhalo + Nbulge;          //  total number of masses

int massbins = 1000;                        //  number of bins to create for
                                            //  numerical integration of M(r)
double dr;                                  //  step in integration of M(r)
double r;                                   //  holds particle radius
int startbin;                               //  step to start integration

double eps = 0.98*pow((double)Ntot, -0.26); //  softening length
double theta = 0.5;                         //  Barnes-Hut theta parameter

double min;                                 //  min of probability range
double max;                                 //  max of probability range
double topbound;                            //  top bound of distributionss

cellnode* root = new cellnode;
bodynode* bodies = new bodynode[Ntot];      //  allocating space for all bodies






main()
{
      
/******************************************************************************/
/*  BUILDING DISK POSITION DISTRIBUTION                                       */
/******************************************************************************/

    std::cout << "Building disk\n";         //  printing disk parameters
    std::cout << "  Mdisk = " << Md << "\n";
    std::cout << "  Ndisk = " << Ndisk << "\n";
    std::cout << "  Radial Scale Length h = " << h << "\n";
    std::cout << "  Vertical Scale Length z0 = " << z0 << "\n";
    std::cout << "  Radial Cutoff = " << diskcut << "\n";
    std::cout << "  Vertical Cutoff = " << thickcut << "\n";
    
    std::vector<Vector> disk(Ndisk, Vector(3)); //  array of disk positions
    
    min = 0.0, max = diskcut;               //  setting max and min radii
    topbound = findtopbound(surfacedist, min, max);   //  finding topbound
    
    for (int i = 0; i < Ndisk; i++)         //  for each particle in disk
    {
        r = rejectionsample(surfacedist, min, max, topbound); //  choose
                                            //  random radius
        disk[i][0] = r;                     //  assigning as such
        disk[i][1] = randomreal(0.0, 2.0*PI); //  randomize azimuthal angle
    }                                       
    
    min = -thickcut;                        //  setting min and max
    max = thickcut;                         //  heights
    topbound = findtopbound(diskthick, min, max);   //  finding topbound
    for (int i = 0; i < Ndisk; i++)         //  for each disk particle
    {
        disk[i][2] = rejectionsample(diskthick, min, max, topbound);    //
                                            //  choosing height randomly
        bodies[i].mass = Md / ((double)Ndisk);  //  assigning masses such that
                                            //  total mass is correct
        bodies[i].pos = cyltocart(disk[i]); //  transforming to cartesian coords
        bodies[i].id = i;                   //  assigning particle id
        bodies[i].DM = false;               //  these are not dark
    }
    

/******************************************************************************/
/*  BUILDING HALO POSITION DISTRIBUTION                                       */
/******************************************************************************/
    
    std::cout << "Building Halo\n";         //  printing parameters
    std::cout << "  Mhalo = " << Mh << "\n";
    std::cout << "  Nhalo = " << Nhalo << "\n";
    std::cout << "  Scale Length rc = " << rc << "\n";
    std::cout << "  Scale Length gamma = " << g << "\n";
    std::cout << "  Cutoff Radius = " << halocut << "\n";

    std::vector<Vector> halo(Nhalo, Vector(3)); //  array of halo positions
    
    min = 0.0, max = halocut;               //  max and min for distribution
    topbound = findtopbound(hernquisthalo, min, max);   //  finding topbound
    for (int i = 0; i < Nhalo; i++)         //  for each bulge particle
    {
        r = rejectionsample(hernquisthalo, min, max, topbound);  //  select r from
                                            //  distribution
        bodies[Ndisk + i].pos = randomsphere(r);    //  randomize position
                                            //  on sphere
        bodies[Ndisk + i].mass = Mh / ((double)Nhalo);  //  normalizing mass
        bodies[Ndisk + i].id = Ndisk + i;   //  setting appropriate ID
        bodies[Ndisk + i].DM = true;        //  these are dark
        halo[i] = carttosphere(bodies[Ndisk + i].pos);  //  saving copy in
                                            //  spherical coords for later
    }
    

/******************************************************************************/
/*  BUILDING BULGE POSITION DISTRIBUTION                                      */
/******************************************************************************/

    std::cout << "Building Bulge\n";
    std::cout << "  Mbulge = " << Mb << "\n";
    std::cout << "  Nbulge = " << Nbulge << "\n";
    std::cout << "  Scale Length a = " << a << "\n";
    std::cout << "  Cutoff Radius = " << bulgecut << "\n";
    
    std::vector<Vector> bulge(Nbulge, Vector(3));   //  array of bulge positions
    
    min = 0.0; max = bulgecut;              //  distribution max and min
    topbound = findtopbound(bulgedist, min, max);   //  finding topbound
    for (int i = 0; i < Nbulge; i++)        //  for each particle
    {
        r = rejectionsample(bulgedist, min, max, topbound);  //  select r from
                                            //  distribution
        bodies[Ndisk + Nhalo + i].pos = randomsphere(r);    //  randomize sphere
                                            //  position
        bodies[Ndisk + Nhalo + i].mass = Mb / (double)Nbulge;   //  setting mass
        bodies[Ndisk + Nhalo + i].id = Ndisk + Nhalo + i;   //  setting IDs
        bulge[i] = carttosphere(bodies[Ndisk + Nhalo + i].pos); //  saving copy
                                            //  in spherical coordinates
    }

    
/******************************************************************************/
/*  Approximating Cumulative Mass Distribution M(r)                           */
/******************************************************************************/
    
    dr = halocut / ((double) massbins);     //  setting separation between mass
                                            //  bins
    
    std::cout << "Approximating Cumulative Mass Distribution M(r)\n";
    std::cout << "  Number of bins = " << massbins << "\n";
    std::cout << "  dr = " << dr << "\n";
    
    std::vector <Vector> Massatr(massbins, Vector(2));  //  Array to hold
                                            //  radius and value of cumulative
                                            //  mass distribution
    
    for (int i = 0; i < massbins; i++)      //  for each mass bin
    {
        Massatr[i][1] = ((double)(i+1))*dr; //  setting radius
        Massatr[i][0] = 0.0;                //  clearing total mass
        
        for (int j = 0; j < Ndisk; j++)     //  for each disk mass
        {
            r = sqrt(disk[j][0]*disk[j][0] + disk[j][2]*disk[j][2]);    // radius
                                            //  in spherical coordinates
            if((r < (double)(i+1)*dr))      //  if radius less than bin radius
            {
                Massatr[i][0] += Md / (double)Ndisk;    //  add mass
            }
        }
        
        for (int j = 0; j < Nhalo; j++)     //  for each halo mass
        {
            r = halo[j][0];                 //  radius
            if((r < (double)(i+1)*dr))      // if radius less than bin radius
            {
                Massatr[i][0] += Mh / (double)Nhalo;    //  add mass
            }
        }
        
        for (int j = 0; j < Nbulge; j++)    // for each bulge mass
        {
            r = bulge[j][0];                //  radius
            if((r < (double)(i+1)*dr))      //  if radius less than bin radius
            {
                Massatr[i][0] += Mb / (double)Nbulge;   //  add mass
            }
        }
    }
    
    
/******************************************************************************/
/*  Setting Halo Velocities                                                   */
/******************************************************************************/
    
    std::cout << "Setting Halo Velocities\n";
    
    double v;                               //  variable to hold speed
    double vesc;                            //  variable to hold escape speed
    
    std::vector<Vector> halovels(Nhalo, Vector(3)); //  Vector to hold halo
                                            //  velocities
    
    for (int i = 0; i < Nhalo; i++)         //  for each halo mass
    {
        r = halo[i][0];                     //  radius
        
        startbin = floor(r/dr);             //  starting index is floor of r/dr
        
        vesc = sqrt(2.0*Massatr[startbin][0]/r);    //  escape velocity
        
        vr2 = 0.0;                          //  clearing radial dispersion
        for (int j = startbin; j < massbins; j++)   //  for each mass bin
        {
            vr2 += hernquisthalo(Massatr[j][1])*dr*Massatr[j][0];   //  add 
        }                                   //  contribution
        vr2 /= (hernquisthalo(r)/(r*r));    //  dividing by halo density at r
        
        
        min = 0.0;                          //  distribution min
        max = 0.95*vesc;                    //  distribution max is 0.95 of
                                            //  escape velocity
        topbound = vr2 / 2.71828;           //  topbound is vr^2 / e
        
        v = rejectionsample(halospeeddist, min, max, topbound); //  selecting
                                            //  absolute speed
        halovels[i] = randomsphere(v);      //  randomizing cartesian velocities
                                            //  on a sphere of radius v
        bodies[Ndisk + i].vel = halovels[i];//  assigning velocities as such
        
    }
    
    
/******************************************************************************/
/*  Setting Bulge Velocities                                                  */
/******************************************************************************/

    
    std::cout << "Setting Bulge Velocities\n";
    
    std::vector<Vector> bulgevels(Nbulge, Vector(3));   // Vector to hold bulge
                                             // velocities
    
    for (int i = 0; i < Nbulge; i++)        //  for each particle
    {
        r = bulge[i][0];                    //  radius
        
        startbin = floor(r / dr);           //  starting index
        
        vesc = sqrt(2.0*Massatr[startbin][0]/r);    //  escape velocity
        
        vr2 = 0.0;                          //  clearing radial dispersion
        
        for (int j = startbin; j < massbins; j++)   //  for each mass bin
        {
            vr2 += bulgedist(Massatr[j][1])*dr*Massatr[j][0];   // add
                                            //  contribution
        }
        vr2 /= bulgedist(r)/(r*r);          //  dividing by halo density at r
        
        
        min = 0.0;                          //  distribution min
        max = 0.95*vesc;                    //  max is 0.95 of escape velocity
        topbound = vr2 / 2.71828;           //  topbound is vr^2 /e
        
        v = rejectionsample(bulgedist, min, max, topbound); //  selecting absolute
                                            //  absolute speed
        bulgevels[i] = randomsphere(v);     //  randomizing constrained cartesian
                                            //  velocities
        bodies[Ndisk + Nhalo + i].vel = bulgevels[i];   //  assining data as such
    }

    
/******************************************************************************/
/*  Setting Disk Velocities                                                   */
/******************************************************************************/

    std::cout << "Setting Disk Velocities\n";
    std::cout << "  Q = " << Q << "\n";
    std::cout << "  Reference radius = " << rref << "\n";
    std::vector<Vector> diskvels(Ndisk, Vector(3)); //  Vector to hold disk
                                            //  velocities
    
    double vz2;                             //  vertical dispersion
    double vc;                              //  circular speed
    double ar;                              //  radial acceleration
    double k;                               //  epicyclic frequency
    double sigmaz2;                         //  azimuthal velocity dispersion
    double sigmaavg;                        //  average dispersion
    double Omega;                           //  angular frequency
    double A;                               //  radial dispersion constant
    int count;                              //  count variable for averaging
    Vector acc(3);                          //  vector to store acceleration
    double as = 0.25*h;
    
    
    maketree(bodies, Ntot, root);           //  making tree
    
    // NORMALIZING RADIAL DISPERSION
    std::cout << "  Normalizing Radial Distribution\n";
    dr = diskcut / 1000.0;                  //  width of annulus in which
                                            //  to average dispersion
    sigmaavg = 0.0;                         //  zeroing average
    for (int i = 0; i < Ndisk; i++)         //  for each disk particle
    {
        r = disk[i][0];                     //  radius
        if (fabs(r - rref) < dr)            //  if radius in annulus
        {                                   //  calculate epicylclic frequency
            k = epicyclicfrequency(bodies, Ntot, i, root, eps, theta, 0.05*dr);
            sigmaavg += 3.36*Sigma(r)/k;    //  calculate dispersion and add to
                                            //  average
            count += 1;                     //  up count
        }  
    }
    
    sigmaavg /= (double)count;              //  divide total by count
    sigmaavg *= Q;                          //  adjust by Q
    
    A = sigmaavg*sigmaavg / Sigma(rref);    //  setting norm constant
    
    
    //  ASSIGNING VELOCITIES
    std::cout << "  Setting particle velocities\n";
    for (int i = 0; i < Ndisk; i++)         //  for every particle
    {
        r = disk[i][0];                     //  radius
        vz2 = PI*z0*Sigma(sqrt(r*r + 2.0*as*as));    //  vertical dispersion
        diskvels[i][2] = gaussianrandom(sqrt(vz2)); //  randomizing vertical
                                            //  with this dispersion
        
        vr2 = A*Sigma(sqrt(r*r + 2.0*as*as));   //  assigning radial dispersion
        diskvels[i][0] = gaussianrandom(sqrt(vr2)); //  randomizing radial
                                            //  dispersion
        
        acc = treeforce(&bodies[i], root, eps, theta);  //  acceleration
        ar = (acc[0]*bodies[i].pos[0] + acc[1]*bodies[i].pos[1])/r; //
                                            //  radial acceleration
        Omega = sqrt(fabs(ar)/r);           //  angular frequency
        k = epicyclicfrequency(bodies, Ntot, i, root, eps, theta, dr);  //
                                            //  epicyclic frequency
        vc = Omega*r;                       //  circular speed
        v = sqrt(fabs(vc*vc + vr2*(1.0 - (k*k)/(4.0*Omega*Omega) - 2.0*r/h)));  //
                                            //  azimuthal streaming velocity
        sigmaz2 = vr2*k*k/(4.0*Omega*Omega);//  azimuthal dispersion
        v += gaussianrandom(sqrt(sigmaz2)); //  adding random azimuthal component
        diskvels[i][1] = v;                 //  assigning azimuthal velocity 
        
                                            //  transforming to cartesian coords
        bodies[i].vel[0] = diskvels[i][0]*cos(disk[i][1])
                            - diskvels[i][1]*sin(disk[i][1]);
        bodies[i].vel[1] = diskvels[i][0]*sin(disk[i][1])
                            + diskvels[i][1]*cos(disk[i][1]);
        bodies[i].vel[2] = diskvels[i][2];
    }
    
    
/******************************************************************************/
/*  Reporting Average Disk Speed                                              */
/******************************************************************************/
    
    v = 0.0;
    
    for (int i = 0; i < Ndisk; i++)
    {
        v += bodies[i].vel.norm();
    }
    
    v /= (double)Ndisk;
    
    std::cout << "Average Disk Particle Speed: " << v << "\n";
    std::cout << "Disk Size: " << diskcut << "\n";
    std::cout << "Disk Crossing Time: " << diskcut / v << "\n";
    
    
    writeinitfile(Ntot, Nhalo, bodies, "data/initfile.txt");
    
    
 
}