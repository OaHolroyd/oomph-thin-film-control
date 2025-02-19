// A driver that solves the control problem of stabilising the flat
// interface of a thin film falling down an inclined plane.

// Standard C++ library includes
#include <iostream>

// Finite-Element library routines
#include "generic.h"

// Project-specific includes
#include "SpineProblem.h"

using namespace std;

using namespace oomph;

// start of main
int main(int argc, char **argv) {
  CommandLineArgs::setup(argc, argv);


  using namespace Global_Physical_Variables;
  // Create the control problem
  unsigned nx = 5;
  unsigned ny = 4;
  unsigned nz = 3;
  int nx_control = 90;
  int ny_control = 30;
  int m_control = 7;
  int p_control = 1;
  SpineControlledFilmProblem<SpineElement<QTaylorHoodElement<3>>, BDF<2>>
      problem(nx, ny, nz, nx_control, ny_control, m_control, p_control);

  // Step up to the start of the controls
  problem.initial_condition(1, 1, 0.01);
  double tburn = 1.0;
  double dtburn = 0.1;
  problem.assign_initial_values_impulsive(dtburn);
  problem.timestep(dtburn, static_cast<int>(tburn / dtburn), 10, 0);

  // // Step with controls turned on
  // double tcontrol = 100.0;
  // double dt = 0.1;
  // problem.timestep(dt, static_cast<int>(tcontrol / dt), 10, 1);
}
