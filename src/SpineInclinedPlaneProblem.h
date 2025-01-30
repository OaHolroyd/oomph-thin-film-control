//
// Created by Oscar Holroyd on 09/01/2025.
//

#ifndef SPINEINCLINEDPLANEPROBLEM_H
#define SPINEINCLINEDPLANEPROBLEM_H

#include "meshes/single_layer_spine_mesh.h"

#include "InclinedPlaneProblem.h"

//======================================================================
/// Create a spine mesh for the problem
//======================================================================
template<class ELEMENT>
class SpineInclinedPlaneMesh : public SingleLayerSpineMesh<ELEMENT> {
public:
  SpineInclinedPlaneMesh(
    const unsigned &nx, const unsigned &ny, const double &lx, const double &ly, TimeStepper *time_stepper_pt
  ) : SingleLayerSpineMesh<ELEMENT>(nx, ny, lx, ly, true, time_stepper_pt) {
  } //end of constructor

  /// General node update function implements pure virtual function
  /// defined in SpineMesh base class and performs specific node update
  /// actions:  along vertical spines
  virtual void spine_node_update(SpineNode *spine_node_pt) {
    //Get fraction along the spine
    double W = spine_node_pt->fraction();
    //Get spine height
    double H = spine_node_pt->h();
    //Set the value of y
    spine_node_pt->x(1) = W * H;
  }
};


//============================================================================
//Specific class for inclined plane problem using spines
//============================================================================
template<class ELEMENT, class TIMESTEPPER>
class SpineInclinedPlaneProblem :
    public InclinedPlaneProblem<ELEMENT, SpineLineFluidInterfaceElement<ELEMENT> > {
public:
  //Constructor
  SpineInclinedPlaneProblem(const unsigned &nx, const unsigned &ny, const double &length): InclinedPlaneProblem<
    ELEMENT, SpineLineFluidInterfaceElement<ELEMENT> >(nx, ny, length) {
    //Set the name
    this->set_output_prefix("spine");

    //Create our one and only timestepper, with adaptive timestepping
    this->add_time_stepper_pt(new TIMESTEPPER);

    //Create the bulk mesh
    this->Bulk_mesh_pt = new SpineInclinedPlaneMesh<ELEMENT>(nx, ny, length, 1.0, this->time_stepper_pt());

    //Create the free surface elements
    this->make_free_surface_elements();

    //Add all sub meshes to the problem
    this->add_sub_mesh(this->Bulk_mesh_pt);
    this->add_sub_mesh(this->Surface_mesh_pt);
    //Create the global mesh
    this->build_global_mesh();

    //Complete the build of the problem
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
