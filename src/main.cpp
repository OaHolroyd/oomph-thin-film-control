//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
//LIC//
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC//
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC//
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC//
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC//
//LIC//====================================================================
//A demo driver that solves the classic fluid flow problem of flow
//of a fluid film along an inclined plane. Stability analysis performed by
//Yih (1963), Benjamin (1957), ... and Blyth & Pozrikidis (2004).

//This is an example of the subtleties involved in even a seemingly simple
//free surface problem.

//Standard C++ library includes
#include <iostream>
#include <cmath>

//Finite-Element library routines
#include "generic.h"

// Project-specific includes
#include "SpineInclinedPlaneProblem.h"

using namespace std;

using namespace oomph;


//start of main
int main(int argc, char **argv) {
  using namespace Global_Physical_Variables;

  //Set the constitutive law for the mesh deformation
  Constitutive_law_pt = new GeneralisedHookean(&Nu);

#ifdef CR_ELEMENT
#define FLUID_ELEMENT QCrouzeixRaviartElement<2>
#else
#define FLUID_ELEMENT QTaylorHoodElement<2>
#endif

  //Initialise physical parameters
  //Scale Reynolds number to be independent of alpha.
  Re = 4.0 / sin(Alpha);

  //Set the direction of gravity
  G[0] = sin(Alpha);
  G[1] = -cos(Alpha);

  //The wall normal to the inlet is in the negative x direction
  Wall_normal.resize(2);
  Wall_normal[0] = -1.0;
  Wall_normal[1] = 0.0;

  // Run once using Spine formulation...
  {
    // Create the problem
    SpineInclinedPlaneProblem<SpineElement<FLUID_ELEMENT >, BDF<2> > problem(30, 4, Length);

    // Solve the steady problem
    problem.solve_steady();

    // Prepare the problem for timestepping
    // (assume that it's been at the flat-film solution for all previous time)
    double dt = 0.1;
    problem.assign_initial_values_impulsive(dt);

    //Timestep it
    problem.timestep(dt, 50, 10, 10);
  }
}
