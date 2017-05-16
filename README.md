# NBodyIntegrator

This code and accompanying documentation was the culmination of a couple years of independent research on N-Body simulations in galactic astrophyics that I undertook as part of my education at Reed College from 2011 to 2013. It was designed to be useful to students at Reed who want either a pre-built code to use to study some problem, or a comprehensive guide to building one's own N-Body simulation code (which can be surprisingly hard to find). In this light, I made an effort to document everything extremely well and to include detailed conceptual guides (see the attached PDFs.) 

The code itself, which uses an adaptive-timestep algorithm and a Barnes-Hut Treecode force calculation routine, is designed to be used as easily as possible, with the user only needing to modify a few input parameters. 
For a PDF guide to Barnes-Hut Treecodes in general and how to edit this code, see the [treecode guide](treecodeguide.pdf).

Additionally included is code used to construct the initial conditions of a disk galaxy using the Hernquist Method. The methods used there are described [here](GalacticInitialConditions.pdf).

## Use

To use you will need the [GNU compiler](https://gcc.gnu.org/). On Macs, this comes bundled with XCode developer tools.

To run the main integrator, first edit the config options in [main.cpp](main.cpp) then compile and run with

```
make -f nbd.make
./nbody
```

The initial conditions routine is found in [makegalaxy.cpp](makegalaxy.cpp). To use, compile and run with

```
make -f init.make
./init
```

Enjoy!
