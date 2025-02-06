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
#include "progress_bar.h"

using namespace std;

using namespace oomph;


//=====================================================================
/// Generic problem class that will form the base class for both
/// spine and elastic mesh-updates of the problem.
/// Templated by the bulk element and interface element types
//====================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
class ControlledFilmProblem : public Problem {
protected:
  /// Bulk fluid mesh
  Mesh *Bulk_mesh_pt;

  /// Mesh for the free surface elements
  Mesh *Surface_mesh_pt;

  /// Prefix for output files
  std::string Output_prefix;

  /// fill the h, q, f arrays with the current values
  void set_hqf(int use_control = 0);

  /// output surface information
  void output_surface();

  /// output mesh information
  void output_mesh();

public:
  /// Storage of the interface/flux/forcing at regular intervals
  double *h, *q, *f;

  /// Information for the control system
  int n_control; // dimension of the control system
  int m_control; // number of actuators
  int p_control; // number of observers

  /// record time, number of time steps, and number of output steps
  double time;
  int step;
  int out_step;


  /**
   * Constructor for the controlled film problem
   *
   * @param n_control the dimension of the control system
   * @param m_control the number of actuators
   * @param p_control the number of observers
   */
  ControlledFilmProblem(const int &n_control, const int &m_control, const int &p_control = 0) {
    // Set up the output directory (always delete the old one)
    std::filesystem::path out_dir = std::filesystem::path("output");
    if (std::filesystem::exists(out_dir)) {
      // delete the directory
      std::filesystem::remove_all(out_dir);
    }
    std::filesystem::create_directory(out_dir);

    // Ensure we use the directory as a prefix
    this->Output_prefix = out_dir.c_str() + std::string("/");

    // don't print loads of internal solver details
    this->Shut_up_in_newton_solve = true;
    this->linear_solver_pt()->disable_doc_time();

    // set control details
    this->n_control = n_control;
    this->m_control = m_control;
    this->p_control = p_control;

    this->h = new double[n_control];
    this->q = new double[n_control];
    this->f = new double[n_control];

    // mesh details
  }

  /**
   * Set the initial condition for the problem
   *
   * @param K the wavenumber of the initial condition
   * @param delta the amplitude of the initial condition
   */
  void initial_condition(int K = 1, double delta = 0.01);

  /**
   * Take a timestep of size dt for n_tsteps
   *
   * @param dt the size of the timestep
   * @param nsteps the number of timesteps to take
   * @param out_step the number of timesteps between outputs
   * @param control_strategy the control strategy to use (0 for none)
   */
  void timestep(
    const double &dt, const unsigned &nsteps, int out_step = 1, int control_strategy = 0
  );

  //Make the free surface elements on the top surface
  void make_free_surface_elements() {
    // Create the (empty) meshes
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

    // Complete the build of the fluid elements by passing physical parameters
    // Find the number of bulk elements
    unsigned n_element = Bulk_mesh_pt->nelement();
    // Loop over all the fluid elements
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
  ~ControlledFilmProblem() {
    //Clear node storage and then delete mesh
    this->Surface_mesh_pt->flush_node_storage();
    delete this->Surface_mesh_pt;
    //Delete the bulk mesh (no need to clear node storage)
    delete this->Bulk_mesh_pt;
    //Delete the time stepper
    delete this->time_stepper_pt();

    delete[] this->h;
    delete[] this->q;
    delete[] this->f;
  }
};


template<class ELEMENT, class INTERFACE_ELEMENT>
void ControlledFilmProblem<ELEMENT, INTERFACE_ELEMENT>::initial_condition(int K, double delta) {
  //Load the namespace
  using namespace Global_Physical_Variables;

  // start at t = 0
  time = 0.0;
  step = 0;
  out_step = 0;

  //Initially set all nodes to the Nusselt flat-film solution
  {
    unsigned n_node = Bulk_mesh_pt->nnode();
    for (unsigned n = 0; n < n_node; n++) {
      double x = Bulk_mesh_pt->node_pt(n)->x(0);
      double y = Bulk_mesh_pt->node_pt(n)->x(1);

      // perturb the y values using the wavenumber and amplitude specified
      Bulk_mesh_pt->node_pt(n)->x(1) *= 1.0 + delta * sin(K * 2.0 * M_PI * (x + 10.0) / Length);

      // set the velocity to the Nusselt flat-film solution
      Bulk_mesh_pt->node_pt(n)->set_value(0, y * (2.0 - y));
      Bulk_mesh_pt->node_pt(n)->set_value(1, 0.0);
    }
  }
} //end of initial_condition


template<class ELEMENT, class INTERFACE_ELEMENT>
void ControlledFilmProblem<ELEMENT, INTERFACE_ELEMENT>::set_hqf(int use_control) {
  using namespace Global_Physical_Variables;

  unsigned int j = 0; // keep track of where we are along the surface mesh
  for (int i = 0; i < n_control; i++) {
    // find the coordinate of the ith measurement point
    double DX = Length / n_control;
    double xi = (DX * (static_cast<double>(i) + 0.5));

    // loop over the free surface elements to find the one containing the point
    //   assuming the surface elements are ordered in the x direction, we can start the search from the same place as
    //   we found the last point
    FaceElement *element;
    Node *n0 = nullptr, *n1 = nullptr;
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

    // ensure that we found a pair of nodes
    assert(n0 != nullptr && n1 != nullptr);

    // linearly interpolate to find h(xi)
    double h0 = n0->x(1);
    double h1 = n1->x(1);
    h[i] = h0 + (h1 - h0) * (xi - n0->x(0)) / (n1->x(0) - n0->x(0));

    // use the leading order approximation to get q(xi)
    q[i] = 2.0 / 3.0 * h[i]; // TODO: do the integration properly

    // use control to set f
    f[i] = (use_control > 0) ? control(xi) : 0.0;
  }
}


template<class ELEMENT, class INTERFACE_ELEMENT>
void ControlledFilmProblem<ELEMENT, INTERFACE_ELEMENT>::output_surface() {
  using namespace Global_Physical_Variables;

  // open the file
  std::ofstream file;
  std::ostringstream filename;
  filename << Output_prefix << "1d_" << out_step << ".dat";
  file.open(filename.str().c_str());

  // write the time
  file << "# " << step << " " << time << std::endl;

  // write the data
  for (unsigned i = 0; i < n_control; i++) {
    double DX = Length / n_control;
    double xi = (DX * (static_cast<double>(i) + 0.5));
    file << xi << " " << h[i] << " " << q[i] << " " << f[i] << std::endl;
  }

  // close the file
  file.close();
}


template<class ELEMENT, class INTERFACE_ELEMENT>
void ControlledFilmProblem<ELEMENT, INTERFACE_ELEMENT>::output_mesh() {
  using namespace Global_Physical_Variables;

  // open the file
  std::ofstream file;
  std::ostringstream filename;
  filename << Output_prefix << "2d_" << out_step << ".dat";
  file.open(filename.str().c_str());

  // write the time
  file << "# " << step << " " << time << std::endl;

  // write all of the mesh data
  Bulk_mesh_pt->output(file, 5);

  // close the file
  file.close();
}


template<class ELEMENT, class INTERFACE_ELEMENT>
void prog_bar_print(void *problem) {
  auto p = (ControlledFilmProblem<ELEMENT, INTERFACE_ELEMENT> *)(problem);

  // find the maximum interfacial deviation
  double max_dh = 0.0;
  for (int i = 0; i < p->n_control; i++) {
    max_dh = std::max(max_dh, std::abs(p->h[i] - 1.0));
  }

  fprintf(stderr, " [t = %3.2f, dhmax = %5g]", p->time, max_dh);
}


template<class ELEMENT, class INTERFACE_ELEMENT>
void ControlledFilmProblem<ELEMENT, INTERFACE_ELEMENT>::timestep(
  const double &dt, const unsigned &nsteps, int out_step, int control_strategy
) {
  // Need to use the Global variables here
  using namespace Global_Physical_Variables;

  // output the initial condition
  set_hqf(control_strategy);
  this->output_surface();
  this->out_step++;

  // if required, set up control variables
  if (control_strategy > 0) {
    control_set(LQR, WR, m_control, p_control, 0.1, 1.0, 0.5, 0.0, Length, n_control, Re, Ca, Alpha);
  }

  //Loop over the desired number of timesteps
  ProgressBar pbar = ProgressBar(nsteps, 50, &prog_bar_print<ELEMENT, INTERFACE_ELEMENT>);
  pbar.start();
  pbar.update(this->step);
  for (unsigned t = 0; t < nsteps; t++) {
    /* Use the control scheme to get the basal forcing */
    // NOTE h and q must be set to the current values
    if (control_strategy > 0) {
      /* compute the actuator strengths */
      control_step(dt, h, q);

      /* set basal velocity from actuator strengths */
      unsigned n_node = this->Bulk_mesh_pt->nboundary_node(0);
      for (unsigned n = 0; n < n_node; n++) {
        Node *node = this->Bulk_mesh_pt->boundary_node_pt(0, n);
        node->set_value(1, control(node->x(0)));
      }
    }

    /* take a timestep of size dt */
    unsteady_newton_solve(dt);
    this->time += dt;
    this->step++;
    set_hqf(control_strategy); // update the h, q, f arrays
    pbar.update(this->step, this);

    // output interface information if required
    if (step % out_step == 0) {
      this->output_surface();
      this->out_step++;
    }
  }

  pbar.end(this);
} //end of timestep


#endif //INCLINEDPLANEPROBLEM_H
