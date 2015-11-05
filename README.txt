This code and accompanying documentation was the culmination of a couple years of independent research on N-Body simulations in galactic astrophyics that I undertook as part of my education at Reed College from 2011 to 2013. It was designed to be useful to students at Reed who want either a pre-built code to use to study some problem, or a comprehensive guide to building one's own N-Body simulation code (which can be surprisingly hard to find). In this light, I made an effort to document everything extremely well and to include detailed conceptual guides (see the attached PDFs.) The code itself, which uses an adaptive-timestep algorithm and a Barnes-Hut Treecode force calculation routine, is designed to be used as easily as possible, with the user only needing to modify a few input parameters. 
For a PDF guide to Barnes-Hut Treecodes in general and how to edit this code, see treecodeguide.pdf

In addition, I included code that I used to construct the initial conditions of a disk galaxy using the hernquist method. For a PDF guide to the Hernquist Method for constructing galactic initial conditions, see GalacticInitialConditions.pdf

RUNNING THE CODE

The initial conditions routine is found in "makegalaxy.cpp". It may be compiled with "make -f init.make" and run with ./init

The N-Body integrator is contained in "main.cpp". It can be compiled with "make -f nbd.make" and run with ./nbody

Enjoy!

 - Noah Muldavin cerca 2015
