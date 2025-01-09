//
// Created by Oscar Holroyd on 09/01/2025.
//

#ifndef ELASTICINCLINEDPLANEPROBLEM_H
#define ELASTICINCLINEDPLANEPROBLEM_H

#include "meshes/simple_rectangular_quadmesh.h"

#include "InclinedPlaneProblem.h"

//======================================================================
/// Create an Elastic mesh for the problem
//======================================================================
template<class ELEMENT>
class ElasticInclinedPlaneMesh :
    public SimpleRectangularQuadMesh<ELEMENT>,
    public SolidMesh {
  //Public functions
public:
  ElasticInclinedPlaneMesh(
    const unsigned &nx, const unsigned &ny, const double &lx, const double &ly, TimeStepper *time_stepper_pt
  ) : SimpleRectangularQuadMesh<ELEMENT>(nx, ny, lx, ly, time_stepper_pt), SolidMesh() {
    // Make the current configuration the undeformed one
    set_lagrangian_nodal_coordinates();
  }
};


//============================================================================
// Specific class for inclined plane problem using pseudo-elastic formulation
//============================================================================
template<class ELEMENT, class TIMESTEPPER>
class ElasticInclinedPlaneProblem : public InclinedPlaneProblem<ELEMENT, ElasticLineFluidInterfaceElement<ELEMENT> > {
public:
  // Constructor
  ElasticInclinedPlaneProblem(
    const unsigned &nx, const unsigned &ny, const double &length
  ) : InclinedPlaneProblem<ELEMENT, ElasticLineFluidInterfaceElement<ELEMENT> >(nx, ny, length) {
    // Set the name
    this->set_output_prefix("elastic");

    // Create our one and only timestepper, with adaptive timestepping
    this->add_time_stepper_pt(new TIMESTEPPER);

    // Create the bulk mesh
    this->Bulk_mesh_pt = new ElasticInclinedPlaneMesh<ELEMENT>(nx, ny, length, 1.0, this->time_stepper_pt());

    // Set the consititutive law for the elements
    unsigned n_element = this->Bulk_mesh_pt->nelement();
    // Loop over all the fluid elements
    for (unsigned e = 0; e < n_element; e++) {
      // Cast to a fluid element
      ELEMENT *temp_pt = dynamic_cast<ELEMENT *>(this->Bulk_mesh_pt->element_pt(e));
      // Set the constitutive law
      temp_pt->constitutive_law_pt() = Global_Physical_Variables::Constitutive_law_pt;
    }

    // Create the traction elements
    this->make_traction_elements();
    // Create the free surface element
    this->make_free_surface_elements();

    // Add all sub meshes to the problem
    this->add_sub_mesh(this->Bulk_mesh_pt);
    this->add_sub_mesh(this->Traction_mesh_pt);
    this->add_sub_mesh(this->Surface_mesh_pt);
    this->add_sub_mesh(this->Point_mesh_pt);
    // Create the global mesh
    this->build_global_mesh();

    // Complete the rest of the build
    this->complete_build();
  } // end of constructor

  /// Update Lagrangian positions after each timestep
  /// (updated-lagrangian approach)
  void actions_after_implicit_timestep() {
    // Now loop over all the nodes and reset their Lagrangian coordinates
    unsigned n_node = this->Bulk_mesh_pt->nnode();
    for (unsigned n = 0; n < n_node; n++) {
      // Cast node to an elastic node
      SolidNode *temp_pt = static_cast<SolidNode *>(this->Bulk_mesh_pt->node_pt(n));
      for (unsigned j = 0; j < 2; j++) {
        temp_pt->xi(j) = temp_pt->x(j);
      }
    }
  } // end of actions_after_implicit_timestep
};

#endif //ELASTICINCLINEDPLANEPROBLEM_H
