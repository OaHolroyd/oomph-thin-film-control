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

  /// Mesh for the traction elements that are added at inlet and outlet
  Mesh *Traction_mesh_pt;

  /// Mesh for the free surface elements
  Mesh *Surface_mesh_pt;

  /// Mesh for the point elements at each end of the free surface
  Mesh *Point_mesh_pt;

  /// Prefix for output files
  std::string Output_prefix;

public:
  /// Generic Constructor (empty)
  InclinedPlaneProblem(const unsigned &nx, const unsigned &ny, const double &length) : Output_prefix("Unset") {
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
  void timestep(const double &dt, const unsigned &n_tsteps, int full_out_step = 2, int interface_out_step = 1);

  /// Actions before the timestep
  /// (update the time-dependent boundary conditions)
  void actions_before_implicit_timestep() {
    //Read out the current time
    double time = this->time_pt()->time();
    //Now add a temporary sinusoidal suction and blowing to the base
    //Amplitude of the perturbation
    double epsilon = 0.01;
    //Loop over the nodes on the base
    unsigned n_node = this->Bulk_mesh_pt->nboundary_node(0);
    for (unsigned n = 0; n < n_node; n++) {
      Node *nod_pt = this->Bulk_mesh_pt->boundary_node_pt(0, n);
      double arg = Global_Physical_Variables::K * nod_pt->x(0);
      double value = sin(arg) * epsilon * time * exp(-time);
      nod_pt->set_value(1, value);
    }
  } //end_of_actions_before_implicit_timestep

  /// Function to add the traction boundary elements to boundaries
  /// 3(inlet) and 1(outlet) of the mesh
  void make_traction_elements() {
    //Create a new (empty mesh)
    Traction_mesh_pt = new Mesh;
    //Inlet boundary conditions (boundary 3)
    {
      unsigned b = 3;
      //Find the number of elements adjacent to mesh boundary
      unsigned n_boundary_element = Bulk_mesh_pt->nboundary_element(b);
      //Loop over these elements and create the traction elements
      for (unsigned e = 0; e < n_boundary_element; e++) {
        NavierStokesTractionElement<ELEMENT> *surface_element_pt = new NavierStokesTractionElement<ELEMENT>(
          Bulk_mesh_pt->boundary_element_pt(b, e),
          Bulk_mesh_pt->face_index_at_boundary(b, e)
        );
        //Add the elements to the mesh
        Traction_mesh_pt->add_element_pt(surface_element_pt);
        //Set the traction function
        surface_element_pt->traction_fct_pt() = &Global_Physical_Variables::hydrostatic_pressure_inlet;
      }
    }

    //Outlet boundary conditions (boundary 1)
    {
      unsigned b = 1;
      //Find the number of elements adjacent to mesh boundary
      unsigned n_boundary_element = Bulk_mesh_pt->nboundary_element(b);
      //Loop over these elements and create the traction elements
      for (unsigned e = 0; e < n_boundary_element; e++) {
        NavierStokesTractionElement<ELEMENT> *surface_element_pt = new NavierStokesTractionElement<ELEMENT>(
          Bulk_mesh_pt->boundary_element_pt(b, e),
          Bulk_mesh_pt->face_index_at_boundary(b, e)
        );
        //Add the elements to the mesh
        Traction_mesh_pt->add_element_pt(surface_element_pt);
        //Set the traction function
        surface_element_pt->traction_fct_pt() = &Global_Physical_Variables::hydrostatic_pressure_outlet;
      }
    }
  } //end of make_traction_elements

  //Make the free surface elements on the top surface
  void make_free_surface_elements() {
    //Create the (empty) meshes
    Surface_mesh_pt = new Mesh;
    Point_mesh_pt = new Mesh;

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

      //Make a point element from left-hand side of the
      //first surface element (note that this relies on knowledge of
      //the element order within the mesh)
      if (e == 0) {
        FluidInterfaceBoundingElement *point_element_pt = surface_element_pt->make_bounding_element(-1);
        //Add element to the point mesh
        Point_mesh_pt->add_element_pt(point_element_pt);
        //Set the capillary number
        point_element_pt->ca_pt() = &Global_Physical_Variables::Ca;
        //Set the wall normal
        point_element_pt->wall_unit_normal_fct_pt() = &Global_Physical_Variables::wall_unit_normal_inlet_fct;
        //Set the contact angle (using the strong version of the constraint)
        point_element_pt->set_contact_angle(&Global_Physical_Variables::Inlet_Angle);
      }

      //Make another point element from the right-hand side of the
      //last surface element (note that this relies on knowledge of
      //the element order within the mesh)
      if (e == n_boundary_element - 1) {
        FluidInterfaceBoundingElement *point_element_pt = surface_element_pt->make_bounding_element(1);
        //Add element to the mesh
        Point_mesh_pt->add_element_pt(point_element_pt);
        //Set the capillary number
        point_element_pt->ca_pt() = &Global_Physical_Variables::Ca;
        // Set the function that specifies the wall normal
        point_element_pt->wall_unit_normal_fct_pt() = &Global_Physical_Variables::wall_unit_normal_outlet_fct;
      }
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
      //Determine whether we are solving an elastic problem or not
      bool elastic = false;
      if (dynamic_cast<SolidNode *>(Bulk_mesh_pt->node_pt(0))) { elastic = true; }

      //Loop over the bottom of the mesh (the wall of the channel)
      unsigned n_node = Bulk_mesh_pt->nboundary_node(0);
      for (unsigned j = 0; j < n_node; j++) {
        //Pin the u- and v- velocities
        Bulk_mesh_pt->boundary_node_pt(0, j)->pin(0);
        Bulk_mesh_pt->boundary_node_pt(0, j)->pin(1);

        //If we are formulating the elastic problem pin both positions
        //of nodes
        if (elastic) {
          static_cast<SolidNode *>(Bulk_mesh_pt->boundary_node_pt(0, j))->pin_position(0);
          static_cast<SolidNode *>(Bulk_mesh_pt->boundary_node_pt(0, j))->pin_position(1);
        }
      }

      //Loop over the inlet and set the Dirichlet condition
      //of no vertical velocity
      n_node = Bulk_mesh_pt->nboundary_node(3);
      for (unsigned j = 0; j < n_node; j++) {
        Bulk_mesh_pt->boundary_node_pt(3, j)->pin(1);

        //If elastic pin horizontal position of nodes
        if (elastic) {
          static_cast<SolidNode *>(Bulk_mesh_pt->boundary_node_pt(3, j))->pin_position(0);
        }
      }

      //Loop over the outlet and set the Dirichlet condition
      //of no vertical velocity
      n_node = Bulk_mesh_pt->nboundary_node(1);
      for (unsigned j = 0; j < n_node; j++) {
        Bulk_mesh_pt->boundary_node_pt(1, j)->pin(1);

        //If elastic pin horizontal position
        if (elastic) {
          static_cast<SolidNode *>(Bulk_mesh_pt->boundary_node_pt(1, j))->pin_position(0);
        }
      }
    }

    //Attach the boundary conditions to the mesh
    std::cout << assign_eqn_numbers() << " in the main problem" << std::endl;
  } //end of complete_build

  /// Generic desructor to clean up the memory allocated in the problem
  ~InclinedPlaneProblem() {
    //Clear node storage before the mesh can be deleted,
    //to avoid double deletion
    this->Point_mesh_pt->flush_node_storage();
    //Then delete the mesh
    delete this->Point_mesh_pt;
    //Clear node storage and then delete mesh
    this->Surface_mesh_pt->flush_node_storage();
    delete this->Surface_mesh_pt;
    //Clear node storage and then delete mesh
    this->Traction_mesh_pt->flush_node_storage();
    delete this->Traction_mesh_pt;
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
      double y = Bulk_mesh_pt->node_pt(n)->x(1);
      //Top row
      Bulk_mesh_pt->node_pt(n)->set_value(0, 0.5 * ReInvFr * sin(Alpha) * (2.0 * y - y * y));
    }
  }

  //Do one steady solve
  steady_newton_solve();

  //Output the full flow field
  std::string filename = Output_prefix;;
  filename.append("_output.dat");
  ofstream file(filename.c_str());
  Bulk_mesh_pt->output(file, 5);
  file.close();
} //end of solve_steady


//----------------------------------------------------------------------
/// Perform n_tsteps timesteps of size dt
//----------------------------------------------------------------------
template<class ELEMENT, class INTERFACE_ELEMENT>
void InclinedPlaneProblem<ELEMENT, INTERFACE_ELEMENT>::timestep(
  const double &dt, const unsigned &n_tsteps, int full_out_step, int interface_out_step
) {
  // Need to use the Global variables here
  using namespace Global_Physical_Variables;

  std::string filename = Output_prefix;
  filename.append("_time_trace.dat");
  ofstream trace(filename.c_str());

  // Initial output of the time and the value of the vertical position at the
  // left and right-hand end of the free surface
  trace << time_pt()->time() << " "
      << Bulk_mesh_pt->boundary_node_pt(2, 0)->value(1)
      << " "
      << Bulk_mesh_pt->
      boundary_node_pt(2, Bulk_mesh_pt->nboundary_node(2) - 1)->x(1)
      << " "
      << std::endl;

  //Loop over the desired number of timesteps
  for (unsigned t = 1; t <= n_tsteps; t++) {
    //Increase the counter
    cout << std::endl;
    cout << "--------------TIMESTEP " << t << " ------------------" << std::endl;

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

    //Output the time and value of the vertical position of the free surface
    //at the left- and right-hand ends
    trace << time_pt()->time() << " "
        << Bulk_mesh_pt->boundary_node_pt(2, 0)->x(1) << " "
        <<
        Bulk_mesh_pt->
        boundary_node_pt(2, Bulk_mesh_pt->nboundary_node(2) - 1)->x(1) << " "
        << std::endl;
  }
} //end of timestep


#endif //INCLINEDPLANEPROBLEM_H
