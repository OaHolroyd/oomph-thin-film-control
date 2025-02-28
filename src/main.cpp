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

  using namespace Global_Physical_Variables;

  // allow overwriting of the default values from the command line
  CommandLineArgs::setup(argc, argv);
  CommandLineArgs::specify_command_line_flag("--lx", &Lx, "streamwise length");
  CommandLineArgs::specify_command_line_flag("--ly", &Ly, "spanwise length");

  CommandLineArgs::specify_command_line_flag("--re", &Re, "Reynolds number");
  CommandLineArgs::specify_command_line_flag("--ca", &Ca, "capillary number");
  CommandLineArgs::specify_command_line_flag("--theta", &Theta,
                                             "incline angle");

  CommandLineArgs::specify_command_line_flag("--tburn", &tburn, "burn in time");
  CommandLineArgs::specify_command_line_flag("--dtburn", &dtburn,
                                             "time step during burn");
  CommandLineArgs::specify_command_line_flag("--tcontrol", &tcontrol,
                                             "control time");
  CommandLineArgs::specify_command_line_flag("--dtcontrol", &dtcontrol,
                                             "time step during control");

  CommandLineArgs::specify_command_line_flag("--nx", &nx,
                                             "streamwise discretisation");
  CommandLineArgs::specify_command_line_flag("--ny", &ny,
                                             "spanwise discretisation");
  CommandLineArgs::specify_command_line_flag("--nz", &nz,
                                             "vertical discretisation");
  CommandLineArgs::specify_command_line_flag(
      "--nx_control", &nx_control, "streamwise control discretisation");
  CommandLineArgs::specify_command_line_flag("--ny_control", &ny_control,
                                             "spanwise control discretisation");

#ifdef OOMPH_HAS_MPI
  // only output if this is rank 0
  if (MPI_Helpers::communicator_pt()->my_rank() == 0) {
    CommandLineArgs::doc_specified_flags();
  } else {
    fprintf(stderr, "Rank %d: not outputting command line arguments\n",
            MPI_Helpers::communicator_pt()->my_rank());
  }
#else
  CommandLineArgs::doc_specified_flags();
#endif

  CommandLineArgs::parse_and_assign();

#ifdef OOMPH_HAS_MPI
  // only output if this is rank 0
  if (MPI_Helpers::communicator_pt()->my_rank() == 0) {
#endif

    // output the values
    fprintf(stderr, "Lx = %g\n", Lx);
    fprintf(stderr, "Ly = %g\n", Ly);
    fprintf(stderr, "Re = %g\n", Re);
    fprintf(stderr, "Ca = %g\n", Ca);
    fprintf(stderr, "Theta = %g\n", Theta);
    fprintf(stderr, "tburn = %g\n", tburn);
    fprintf(stderr, "dtburn = %g\n", dtburn);
    fprintf(stderr, "tcontrol = %g\n", tcontrol);
    fprintf(stderr, "dtcontrol = %g\n", dtcontrol);
    fprintf(stderr, "nx = %d\n", nx);
    fprintf(stderr, "ny = %d\n", ny);
    fprintf(stderr, "nz = %d\n", nz);
    fprintf(stderr, "nx_control = %d\n", nx_control);
    fprintf(stderr, "ny_control = %d\n", ny_control);

#ifdef OOMPH_HAS_MPI
  }
#endif

#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

  return 0;

  // Create the control problem
  int m_control = 7;
  int p_control = 1;
  SpineControlledFilmProblem<SpineElement<QTaylorHoodElement<3>>, BDF<2>>
      problem(nx, ny, nz, nx_control, ny_control, m_control, p_control);

#ifdef OOMPH_HAS_MPI
  problem.distribute();
#endif

  // Step up to the start of the controls
  problem.initial_condition(1, 1, 0.01, 0.8);
  problem.assign_initial_values_impulsive(
      dtburn); // TODO: mucks up the initial condition
  problem.timestep(dtburn, static_cast<int>(tburn / dtburn), 1, 0);

  // Step with controls turned on
  // if (tcontrol > 0.0) {
  //   problem.timestep(dtcontrol, static_cast<int>(tcontrol / dtcontrol), 1,
  //   1);
  // }

  // Finalise MPI
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif
}
