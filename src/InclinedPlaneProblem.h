//
// Created by Oscar Holroyd on 09/01/2025.
//

#ifndef INCLINEDPLANEPROBLEM_H
#define INCLINEDPLANEPROBLEM_H

//Standard C++ library includes
#include <iostream>
#include <filesystem>
#include <fstream>
#include <string>
#include <cmath>

//Finite-Element library routines
#include "generic.h"
#include "navier_stokes.h"
#include "fluid_interface.h"

#include "GlobalPhysicalVariables.h"
#include "control.h"

using namespace std;

using namespace oomph;


//=====================================================================
/// Generic problem class that will form the base class for both
/// spine and elastic mesh-updates of the problem.
/// Templated by the bulk element and interface element types
//====================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
class InclinedPlaneProblem : public Problem {
protected:
  /// Bulk fluid mesh
  Mesh *Bulk_mesh_pt;

  /// Mesh for the free surface elements
  Mesh *Surface_mesh_pt;

  /// Prefix for output files
  std::string Output_prefix;

public:
  /// Generic Constructor (empty)
  InclinedPlaneProblem(const unsigned &nx, const unsigned &ny, const double &length) : Output_prefix("Unset") {
    this->Shut_up_in_newton_solve = true; // don't print loads of solver details
  }

  /// Set the output prefix
  void set_output_prefix(const std::string &prefix) {
    // Set up the output directory (always delete the old one)
    std::filesystem::path out_dir = std::filesystem::path("output");
    if (std::filesystem::exists(out_dir)) {
      // delete the directory
      std::filesystem::remove_all(out_dir);
    }
    std::filesystem::create_directory(out_dir);

    // Open an output file
    this->Output_prefix = out_dir.c_str() + std::string("/") + prefix;
  }

  /// Solve the steady problem
  void solve_steady();

  /// Take n_tsteps timesteps of size dt
  void timestep(
    const double &dt, const unsigned &n_tsteps, int full_out_step = 2, int interface_out_step = 1,
    int control_strategy = 0
  );

  // /// Actions before the timestep
  // /// (update the time-dependent boundary conditions)
  // void actions_before_implicit_timestep() {
  //   //Read out the current time
  //   double time = this->time_pt()->time();
  //   //Now add a temporary sinusoidal suction and blowing to the base
  //   //Amplitude of the perturbation
  //   double epsilon = 0.01;
  //   //Loop over the nodes on the base
  //   unsigned n_node = this->Bulk_mesh_pt->nboundary_node(0);
  //   for (unsigned n = 0; n < n_node; n++) {
  //     Node *nod_pt = this->Bulk_mesh_pt->boundary_node_pt(0, n);
  //     double arg = Global_Physical_Variables::K * nod_pt->x(0);
  //     double value = sin(arg) * epsilon * time * exp(-time);
  //     nod_pt->set_value(1, value);
  //   }
  // } //end_of_actions_before_implicit_timestep

  //Make the free surface elements on the top surface
  void make_free_surface_elements() {
    //Create the (empty) meshes
    Surface_mesh_pt = new Mesh;

    //The free surface is on the boundary 2
    unsigned b = 2;
    unsigned n_boundary_element = Bulk_mesh_pt->nboundary_element(b);
    //Loop over the elements and create the appropriate interface elements
    for (unsigned e = 0; e < n_boundary_element; e++) {
      INTERFACE_ELEMENT *surface_element_pt = new INTERFACE_ELEMENT(
        Bulk_mesh_pt->boundary_element_pt(b, e),
        Bulk_mesh_pt->face_index_at_boundary(b, e)
      );
      //Add elements to the mesh
      Surface_mesh_pt->add_element_pt(surface_element_pt);
      //Assign the capillary number to the free surface
      surface_element_pt->ca_pt() = &Global_Physical_Variables::Ca;
    }
  } //end of make_free_surface_elements

  /// Complete the build of the problem setting all standard
  /// parameters and boundary conditions
  void complete_build() {
    using namespace Global_Physical_Variables;

    //Complete the build of the fluid elements by passing physical parameters
    //Find the number of bulk elements
    unsigned n_element = Bulk_mesh_pt->nelement();
    //Loop over all the fluid elements
    for (unsigned e = 0; e < n_element; e++) {
      //Cast to a fluid element
      ELEMENT *temp_pt = dynamic_cast<ELEMENT *>(Bulk_mesh_pt->element_pt(e));

      //Set the Reynolds number
      temp_pt->re_pt() = &Re;
      //The Strouhal number is 1, so ReSt = Re
      temp_pt->re_st_pt() = &Re;
      //Set the Reynolds number / Froude number
      temp_pt->re_invfr_pt() = &ReInvFr;
      //Set the direction of gravity
      temp_pt->g_pt() = &G;
    }

    //------------Set the boundary conditions for this problem----------
    {
      //Loop over the bottom of the mesh (the wall of the channel)
      unsigned n_node = Bulk_mesh_pt->nboundary_node(0);
      for (unsigned j = 0; j < n_node; j++) {
        //Pin the u- and v- velocities
        Bulk_mesh_pt->boundary_node_pt(0, j)->pin(0);
        Bulk_mesh_pt->boundary_node_pt(0, j)->pin(1);
      }
    }

    //Attach the boundary conditions to the mesh
    std::cout << assign_eqn_numbers() << " in the main problem" << std::endl;
  } //end of complete_build

  /// Generic desructor to clean up the memory allocated in the problem
  ~InclinedPlaneProblem() {
    //Clear node storage and then delete mesh
    this->Surface_mesh_pt->flush_node_storage();
    delete this->Surface_mesh_pt;
    //Delete the bulk mesh (no need to clear node storage)
    delete this->Bulk_mesh_pt;
    //Delete the time stepper
    delete this->time_stepper_pt();
  }
};


//-------------------------------------------------------------------------
/// Solve the steady problem
//-------------------------------------------------------------------------
template<class ELEMENT, class INTERFACE_ELEMENT>
void InclinedPlaneProblem<ELEMENT, INTERFACE_ELEMENT>::solve_steady() {
  //Load the namespace
  using namespace Global_Physical_Variables;

  //Initially set all nodes to the Nusselt flat-film solution
  {
    unsigned n_node = Bulk_mesh_pt->nnode();
    for (unsigned n = 0; n < n_node; n++) {
      double x = Bulk_mesh_pt->node_pt(n)->x(0);
      double y = Bulk_mesh_pt->node_pt(n)->x(1);

      // perturb the y values
      Bulk_mesh_pt->node_pt(n)->x(1) *= 1.0 + 0.01 * sin(2.0 * M_PI * (x + 10.0) / Length);

      // set the velocity
      Bulk_mesh_pt->node_pt(n)->set_value(0, y * (2.0 - y));
      Bulk_mesh_pt->node_pt(n)->set_value(1, 0.0);
    }
  }
} //end of solve_steady


//----------------------------------------------------------------------
/// Perform n_tsteps timesteps of size dt
//----------------------------------------------------------------------
template<class ELEMENT, class INTERFACE_ELEMENT>
void InclinedPlaneProblem<ELEMENT, INTERFACE_ELEMENT>::timestep(
  const double &dt, const unsigned &n_tsteps, int full_out_step, int interface_out_step, int control_strategy
) {
  // Need to use the Global variables here
  using namespace Global_Physical_Variables;

  if (1) {
    {
      std::ofstream file;
      std::ostringstream filename;
      filename << Output_prefix << "_step_" << 0 << ".dat";
      file.open(filename.str().c_str());
      Bulk_mesh_pt->output(file, 5);
      file.close();
    } {
      std::ofstream file;
      std::ostringstream filename;
      filename << Output_prefix << "_interface_" << 0 << ".dat";
      file.open(filename.str().c_str());
      Surface_mesh_pt->output(file, 5);
      file.close();
    }
  }

  // if required, set up control variables
  int control_n = 100;
  int control_m = 7;
  double *actuators;
  double *h;
  double *q;
  if (control_strategy > 0) {
    control_set(LQR, WR, control_m, 1, 0.1, 1.0, 0.5, 0.0, Length, control_n, Re, Ca, Alpha);
    actuators = new double[control_m];
    h = new double[control_n];
    q = new double[control_n];
  }

  //Loop over the desired number of timesteps
  for (unsigned t = 1; t <= n_tsteps; t++) {
    //Increase the counter
    cout << std::endl;
    cout << "--------------TIMESTEP (" << control_strategy << ") " << t << " ------------------" << std::endl;


    // Use the control scheme to get the basal forcing
    if (control_strategy > 0) {
      /* set the values of h and q */
      unsigned int j = 0; // keep track of where we are along the surface mesh
      for (int i = 0; i < control_n; i++) {
        // find the coordinate of the ith measurement point
        double DX = Length / control_n;
        double xi = (DX * (static_cast<double>(i) + 0.5));

        // loop over the free surface elements to find the one containing the point
        //   assuming the surface elements are ordered in the x direction, we can start the search from the same place as
        //   we found the last point
        FaceElement *element;
        Node *n0, *n1;
        for (; j < Surface_mesh_pt->nelement(); j++) {
          element = dynamic_cast<FaceElement *>(Surface_mesh_pt->element_pt(j));

          // assert that the element has 3 nodes TODO: why?
          assert(element->nnode() == 3);

          // get the x coordinates of the nodes
          double x0 = element->node_pt(0)->x(0);
          double x1 = element->node_pt(1)->x(0);
          double x2 = element->node_pt(2)->x(0);

          // find the enclosing nodes
          if (xi >= x0 && xi <= x2) {
            if (xi < x1) {
              n0 = element->node_pt(0);
              n1 = element->node_pt(1);
            } else {
              n0 = element->node_pt(1);
              n1 = element->node_pt(2);
            }
            break;
          }
        }

        // linearly interpolate to find h(xi)
        double h0 = n0->x(1);
        double h1 = n1->x(1);
        h[i] = h0 + (h1 - h0) * (xi - n0->x(0)) / (n1->x(0) - n0->x(0));

        // use the leading order approximation to get q(xi)
        q[i] = 2.0 / 3.0 * h[i]; // TODO: do the integration properly
      }

      /* compute the actuator strengths */
      control_step(dt, h, q);

      /* set basal velocity from actuator strengths */
      unsigned n_node = this->Bulk_mesh_pt->nboundary_node(0);
      for (unsigned n = 0; n < n_node; n++) {
        Node *node = this->Bulk_mesh_pt->boundary_node_pt(0, n);
        node->set_value(1, control(node->x(0)));
      }
    }

    //Take a timestep of size dt
    unsteady_newton_solve(dt);

    //Uncomment to get full solution output
    //Change this number to get output every n steps
    if (t % full_out_step == 0) {
      std::ofstream file;
      std::ostringstream filename;
      filename << Output_prefix << "_step_" << t << ".dat";
      file.open(filename.str().c_str());
      Bulk_mesh_pt->output(file, 5);
      file.close();
    }

    //Always output the interface
    if (t % interface_out_step == 0) {
      std::ofstream file;
      std::ostringstream filename;
      filename << Output_prefix << "_interface_" << t << ".dat";
      file.open(filename.str().c_str());
      Surface_mesh_pt->output(file, 5);
      file.close();
    }
  }

  // if required, tear down control variables
  if (control_strategy > 0) {
    delete[] actuators;
    delete[] h;
    delete[] q;
  }
} //end of timestep


#endif //INCLINEDPLANEPROBLEM_H
