//
// Created by Oscar Holroyd on 09/01/2025.
//

#ifndef SPINEINCLINEDPLANEPROBLEM_H
#define SPINEINCLINEDPLANEPROBLEM_H

#include "meshes/rectangular_quadmesh.h"

#include "InclinedPlaneProblem.h"

//======================================================================
/// Create a spine mesh for the problem
//======================================================================
template<class ELEMENT>
class SpineInclinedPlaneMesh : public RectangularQuadMesh<ELEMENT>, public SpineMesh {
public:
  SpineInclinedPlaneMesh(
    const unsigned &nx, const unsigned &ny, const double &lx, const double &ly, TimeStepper *time_stepper_pt
  ) : RectangularQuadMesh<ELEMENT>(nx, ny, lx, ly, true, time_stepper_pt), SpineMesh() {
    //Find the number of linear points in the element
    unsigned n_p = dynamic_cast<ELEMENT *>(finite_element_pt(0))->nnode_1d();
    //Reserve storage for the number of spines
    Spine_pt.reserve((n_p - 1) * nx + 1); // TODO: can we remove 1 because th s;ine at the start and end are the same?

    //Create single pointer to a spine
    Spine *new_spine_pt = 0;

    //Now loop over the elements horizontally
    for (unsigned long j = 0; j < nx; j++) {
      //In most elements, we don't assign a spine to the last column,
      //beacuse that will be done by the next element
      unsigned n_pmax = n_p - 1;
      //In the last element, however, we must assign the final spine
      if (j == nx - 1) { n_pmax = n_p; }

      //Loop over all nodes horizontally
      for (unsigned l2 = 0; l2 < n_pmax; l2++) {
        //Create a new spine with unit height and add to the mesh
        new_spine_pt = new Spine(1.0);
        Spine_pt.push_back(new_spine_pt);

        // Get the node
        SpineNode *nod_pt = element_node_pt(j, l2);
        //Set the pointer to spine
        nod_pt->spine_pt() = new_spine_pt;
        //Set the fraction
        nod_pt->fraction() = 0.0;
        // Pointer to the mesh that implements the update fct
        nod_pt->spine_mesh_pt() = this;

        //Loop vertically along the spine
        //Loop over the elements
        for (unsigned long i = 0; i < ny; i++) {
          //Loop over the vertical nodes, apart from the first
          for (unsigned l1 = 1; l1 < n_p; l1++) {
            // Get the node
            SpineNode *nod_pt = element_node_pt(i * nx + j, l1 * n_p + l2);
            //Set the pointer to the spine
            nod_pt->spine_pt() = new_spine_pt;
            //Set the fraction
            nod_pt->fraction() = (double(i) + double(l1) / double(n_p - 1)) / double(ny);
            // Pointer to the mesh that implements the update fct
            nod_pt->spine_mesh_pt() = this;
          }
        }
      }
    } //End of horizontal loop over elements
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
