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

  using namespace Global_Variables;

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
      "--nxcontrol", &nx_control, "streamwise control discretisation");
  CommandLineArgs::specify_command_line_flag("--nycontrol", &ny_control,
                                             "spanwise control discretisation");
  CommandLineArgs::specify_command_line_flag("--mcontrol", &m_control,
                                             "number of actuators");
  CommandLineArgs::specify_command_line_flag("--pcontrol", &p_control,
                                             "number of observers");

  CommandLineArgs::specify_command_line_flag("--nomumps",
                                             "Don't use the MUMPS solver");

  CommandLineArgs::parse_and_assign();

  if (CommandLineArgs::command_line_flag_has_been_set("--no_mumps")) {
    use_mumps = false;
  } else {
    use_mumps = true;
  }

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
    fprintf(stderr, "m_control = %d\n", m_control);
    fprintf(stderr, "p_control = %d\n", p_control);
    fprintf(stderr, "MUMPS: %d\n", use_mumps);

#ifdef OOMPH_HAS_MPI
  }
#endif

  // Create the control problem
  SpineControlledFilmProblem<SpineElement<QTaylorHoodElement<3>>, BDF<2>>
      problem(nx, ny, nz, nx_control, ny_control, m_control, p_control,
              use_mumps);

  // Initial condition
  problem.initial_condition(1, 1, 0.01, 0.8);
  problem.assign_initial_values_impulsive(
      dtburn); // TODO: mucks up the initial condition

  // Step up to the start of the controls
  if (tburn > 0.0) {
    problem.timestep(dtburn, static_cast<int>(tburn / dtburn), 1, UNCONTROLLED);
  }

  // Step with controls turned on
  if (tcontrol > 0.0) {
    problem.timestep(dtcontrol, static_cast<int>(tcontrol / dtcontrol), 1,
                     PAIR);
  }

  // Finalise MPI
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif
}
