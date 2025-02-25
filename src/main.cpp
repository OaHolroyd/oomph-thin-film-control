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
  using namespace Global_Physical_Variables;

#ifdef CR_ELEMENT
#define FLUID_ELEMENT QCrouzeixRaviartElement<2>
#else
#define FLUID_ELEMENT QTaylorHoodElement<2>
#endif

  // Create the control problem
  unsigned nx = 100;
  unsigned ny = 6;
  int n_control = 100;
  int m_control = 7;
  int p_control = 1;
  SpineControlledFilmProblem<SpineElement<FLUID_ELEMENT>, BDF<2>> problem(
      nx, ny, n_control, m_control, p_control);

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
