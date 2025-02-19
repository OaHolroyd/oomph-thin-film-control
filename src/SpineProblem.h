//
// Created by Oscar Holroyd on 09/01/2025.
//

#ifndef SPINEINCLINEDPLANEPROBLEM_H
#define SPINEINCLINEDPLANEPROBLEM_H

#include "navier_stokes.h"
#include "fluid_interface.h"

#include "single_layer_cubic_spine_mesh.h"
#include "Problem.h"

// ======================================================================
//   Create a spine mesh for the problem
// ======================================================================
template<class ELEMENT>
class SpineFilmMesh : public SingleLayerCubicSpineMesh<ELEMENT> {
public:
  /**
   * Constructor for the spine inclined plane mesh
   *
   * @param nx Number of elements in the x direction
   * @param ny Number of elements in the y direction
   * @param nz Number of elements in the z direction
   * @param lx Length of the domain in the x direction
   * @param ly Length of the domain in the y direction
   * @param lz Length of the domain in the z direction
   * @param time_stepper_pt Pointer to the time stepper
   */
  SpineFilmMesh(
    const unsigned &nx, const unsigned &ny, const unsigned &nz,
    const double &lx, const double &ly, const double &lz,
    TimeStepper *time_stepper_pt
  ) : SingleLayerCubicSpineMesh<ELEMENT>(nx, ny, nz, lx, ly, lz, true, true, time_stepper_pt) {
  } //end of constructor
};


// ============================================================================
//  Specific class for controlled film problem using spines
// ============================================================================
template<class ELEMENT, class TIMESTEPPER>
class SpineControlledFilmProblem :
    public ControlledFilmProblem<ELEMENT, SpineSurfaceFluidInterfaceElement<ELEMENT> > {
public:
  /**
   * Constructor for the spine controlled film problem
   *
   * @param nx Number of elements in the x direction
   * @param ny Number of elements in the y direction
   * @param nz Number of elements in the z direction
   * @param nx_control Number of elements in the x direction of the control system
   * @param ny_control Number of elements in the y direction of the control system
   * @param m_control the number of actuators
   * @param p_control the number of observers
   */
  SpineControlledFilmProblem(
    const unsigned &nx, const unsigned &ny, const unsigned &nz,
    const int &nx_control, const int &ny_control, const int &m_control, const int &p_control
  ): ControlledFilmProblem<ELEMENT, SpineSurfaceFluidInterfaceElement<ELEMENT> >(nx_control, ny_control, m_control, p_control) {
    using namespace Global_Physical_Variables; // to access the length

    // create our one and only timestepper, with adaptive timestepping
    this->add_time_stepper_pt(new TIMESTEPPER);

    // create the bulk mesh
    this->Bulk_mesh_pt = new SpineFilmMesh<ELEMENT>(nx, ny, nz, Lx, Ly, 1.0, this->time_stepper_pt());

    // create the free surface elements
    this->make_free_surface_elements();

    // add all sub meshes to the problem
    this->add_sub_mesh(this->Bulk_mesh_pt);
    this->add_sub_mesh(this->Surface_mesh_pt);

    // create the global mesh
    this->build_global_mesh();

    // complete the build of the problem
    this->complete_build();
  }

  /// Spine heights/lengths are unknowns in the problem so their
  /// values get corrected during each Newton step. However,
  /// changing their value does not automatically change the
  /// nodal positions, so we need to update all of them
  void actions_before_newton_convergence_check() {
    this->Bulk_mesh_pt->node_update();
  }
};

#endif //SPINEINCLINEDPLANEPROBLEM_H
