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

#ifdef CR_ELEMENT
#define FLUID_ELEMENT QCrouzeixRaviartElement<2>
#else
#define FLUID_ELEMENT QTaylorHoodElement<2>
#endif

  // Set the direction of gravity
  // TODO: should be in the Global_Physical_Variables namespace
  G[0] = 2.0;
  G[1] = -2.0 / tan(Alpha);

  // Create the control problem
  unsigned nx = 100;
  unsigned ny = 6;
  int n_control = 100;
  int m_control = 7;
  int p_control = 1;
  SpineControlledFilmProblem<SpineElement<FLUID_ELEMENT >, BDF<2> > problem(nx, ny, n_control, m_control, p_control);

  // Step up to the start of the controls
  problem.initial_condition(1, 0.01);
  double tburn = 200.0;
  double dtburn = 0.1;
  problem.assign_initial_values_impulsive(dtburn);
  problem.timestep(dtburn, static_cast<int>(tburn / dtburn), 10, 0);

  // Step with controls turned on
  double tcontrol = 100.0;
  double dt = 0.1;
  problem.timestep(dt, static_cast<int>(tcontrol / dt), 10, 1);
}
