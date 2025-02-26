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

#ifdef OOMPH_HAS_MPI
  // Initialise MPI
  MPI_Helpers::init(argc, argv);
#endif

  CommandLineArgs::setup(argc, argv);

  using namespace Global_Physical_Variables;

  // Create the control problem
  unsigned nx = 32;
  unsigned ny = 32;
  unsigned nz = 4;
  int nx_control = 32;
  int ny_control = 32;
  int m_control = 7;
  int p_control = 1;
  SpineControlledFilmProblem<SpineElement<QTaylorHoodElement<3>>, BDF<2>>
      problem(nx, ny, nz, nx_control, ny_control, m_control, p_control);

#ifdef OOMPH_HAS_MUMPS
  // Use mumps if available
  problem.linear_solver_pt() = new MumpsSolver;
#endif

  // Step up to the start of the controls
  problem.initial_condition(1, 1, 0.01, 0.8);
  double tburn = 200.0;
  double dtburn = 0.5;
  problem.assign_initial_values_impulsive(
      dtburn); // TODO: this might be mucking up the initial condition
  problem.timestep(dtburn, static_cast<int>(tburn / dtburn), 1, 0);

  // // Step with controls turned on
  // double tcontrol = 100.0;
  // double dt = 0.1;
  // problem.timestep(dt, static_cast<int>(tcontrol / dt), 10, 1);

  // Finalise MPI
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif
}
