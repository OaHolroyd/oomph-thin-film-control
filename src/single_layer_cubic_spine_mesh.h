//
// Created by Oscar Holroyd on 17/02/2025.
//

#ifndef SINGLE_LAYER_CUBIC_SPINE_MESH_H
#define SINGLE_LAYER_CUBIC_SPINE_MESH_H

// oomph-lib includes
#include "cubic_brick_mesh.h"
#include "generic/spines.h"

namespace oomph {
//======================================================================
/// Spine mesh class derived from standard cubic 3D mesh.
/// The mesh contains a layer of spinified fluid elements (of type ELEMENT;
/// e.g  SpineElement<QCrouzeixRaviartElement<3>)
/// for 3D problems, in which the interface's vertical position can vary
//======================================================================
template <class ELEMENT>
class SingleLayerCubicSpineMesh : public CubicBrickMesh<ELEMENT>,
                                  public SpineMesh {
public:
  /// Constructor: Pass number of elements in x-direction, number of
  /// elements in y-direction, number of elements in z-direction,
  /// lengths in x- and y- directions, height of layer, and pointer
  /// to timestepper (defaults to Steady timestepper)
  SingleLayerCubicSpineMesh(
      const unsigned &nx, const unsigned &ny, const unsigned &nz,
      const double &lx, const double &ly, const double &h,
      TimeStepper *time_stepper_pt = &Mesh::Default_TimeStepper);

  /// Constructor: Pass number of elements in x-direction, number of
  /// elements in y-direction, number of elements in z-direction,
  /// lengths in x- and y- directions, height of layer, whether the
  /// mesh is periodic in x, whether the mesh is periodic in y, and pointer
  /// to timestepper (defaults to Steady timestepper)
  SingleLayerCubicSpineMesh(
      const unsigned &nx, const unsigned &ny, const unsigned &nz,
      const double &lx, const double &ly, const double &h,
      const bool &periodic_in_x, const bool &periodic_in_y,
      TimeStepper *time_stepper_pt = &Mesh::Default_TimeStepper);

  /// General node update function implements pure virtual function
  /// defined in SpineMesh base class and performs specific node update
  /// actions:  along vertical spines
  virtual void spine_node_update(SpineNode *spine_node_pt) {
    // Get fraction along the spine
    double W = spine_node_pt->fraction();
    // Get spine height
    double H = spine_node_pt->h();
    // Set the value of z
    spine_node_pt->x(2) = this->Zmin + W * H;
  }

protected:
  /// Helper function to build a single spine for the nth node of element e
  virtual void build_spine(unsigned e, unsigned n);

  /// Helper function to match a spine on a periodic boundary with it's matching
  /// counterpart. The mode is a character that can be 'x', 'y', or 'c' to
  /// denote x-periodicity, y-periodicity, or both (corner), respectively.
  virtual void match_spine_periodic(unsigned e, unsigned n, char mode);

  /// Helper function to actually build the single-layer spine mesh
  /// (called from various constructors)
  virtual void build_single_layer_mesh(TimeStepper *time_stepper_pt);
};

//===========================================================================
/// Constructor for spine 3D mesh: Pass number of elements in x-direction,
/// number of elements in y-direction, number elements in z-direction,
/// length, width and height of layer,
/// and pointer to timestepper (defaults to Static timestepper).
//===========================================================================
template <class ELEMENT>
SingleLayerCubicSpineMesh<ELEMENT>::SingleLayerCubicSpineMesh(
    const unsigned &nx, const unsigned &ny, const unsigned &nz,
    const double &lx, const double &ly, const double &h,
    TimeStepper *time_stepper_pt)
    : CubicBrickMesh<ELEMENT>(nx, ny, nz, lx, ly, h, true, time_stepper_pt) {
  // Mesh can only be built with 3D Qelements.
  MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

  // Mesh can only be built with spine elements
  MeshChecker::assert_geometric_element<SpineFiniteElement, ELEMENT>(3);

  // Now build the single layer mesh on top of the existing mesh
  build_single_layer_mesh(time_stepper_pt);
}

//===========================================================================
/// Constructor for spine 3D mesh: Pass number of elements in x-direction,
/// number of elements in y-direction, number elements in z-direction,
/// length, width and height of layer, whether the mesh is periodic in x,
/// whether the mesh is periodic in y,  and pointer to timestepper (defaults
/// to Static timestepper).
//===========================================================================
template <class ELEMENT>
SingleLayerCubicSpineMesh<ELEMENT>::SingleLayerCubicSpineMesh(
    const unsigned &nx, const unsigned &ny, const unsigned &nz,
    const double &lx, const double &ly, const double &h,
    const bool &periodic_in_x, const bool &periodic_in_y,
    TimeStepper *time_stepper_pt)
    : CubicBrickMesh<ELEMENT>(nx, ny, nz, lx, ly, h, periodic_in_x,
                              periodic_in_y, true, time_stepper_pt) {
  // Mesh can only be built with 3D Qelements.
  MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

  // Mesh can only be built with spine elements
  MeshChecker::assert_geometric_element<SpineFiniteElement, ELEMENT>(3);

  // Now build the single layer mesh on top of the existing mesh
  build_single_layer_mesh(time_stepper_pt);
}

template <class ELEMENT>
void SingleLayerCubicSpineMesh<ELEMENT>::build_spine(const unsigned e,
                                                     const unsigned n) {
  // Read out the number of elements in the x-direction
  unsigned n_x = this->Nx;
  // Read out the number of elements in the y-direction
  unsigned n_y = this->Ny;
  // Read out the number of elements in the z-direction
  unsigned n_z = this->Nz;
  // Read out number of linear points in the element
  unsigned n_p = dynamic_cast<ELEMENT *>(finite_element_pt(0))->nnode_1d();

  // Check if this spine is actually on a periodic boundary
  if (this->Xperiodic && this->Yperiodic && (e == (n_x * n_y - 1)) &&
      (n == (n_p * n_p - 1))) {
    match_spine_periodic(e, n, 'c');
    return;
  }

  if (this->Xperiodic && ((e % n_x) == (n_x - 1)) && ((n % n_p) == (n_p - 1))) {
    match_spine_periodic(e, n, 'x');
    return;
  }

  if (this->Yperiodic && (e >= ((n_y - 1) * n_x)) && (n >= ((n_p - 1) * n_p))) {
    match_spine_periodic(e, n, 'y');
    return;
  }

  // Assign the new spine with length h
  Spine *new_spine_pt = new Spine(1.0);
  Spine_pt.push_back(new_spine_pt);

  // Get pointer to node
  {
    SpineNode *nod_pt = element_node_pt(e, n);
    // Set the pointer to the spine
    nod_pt->spine_pt() = new_spine_pt;
    // Set the fraction
    nod_pt->fraction() = 0.0;
    // Pointer to the mesh that implements the update fct
    nod_pt->spine_mesh_pt() = this;
  }

  // Loop vertically along the spine
  // Loop over the elements
  for (unsigned long k = 0; k < n_z; k++) {
    // Loop over the vertical nodes, apart from the first
    for (unsigned l3 = 1; l3 < n_p; l3++) {
      // Get pointer to node
      SpineNode *nod_pt =
          element_node_pt(e + k * n_x * n_y, l3 * n_p * n_p + n);
      // Set the pointer to the spine
      nod_pt->spine_pt() = new_spine_pt;
      // Set the fraction
      nod_pt->fraction() =
          (double(k) + double(l3) / double(n_p - 1)) / double(n_z);
      // Pointer to the mesh that implements the update fct
      nod_pt->spine_mesh_pt() = this;
    }
  }
}

template <class ELEMENT>
void SingleLayerCubicSpineMesh<ELEMENT>::match_spine_periodic(const unsigned e,
                                                              const unsigned n,
                                                              const char mode) {
  // Read out the number of elements in the x-direction
  unsigned n_x = this->Nx;
  // Read out the number of elements in the y-direction
  unsigned n_y = this->Ny;
  // Read out the number of elements in the z-direction
  unsigned n_z = this->Nz;
  // Read out number of linear points in the element
  unsigned n_p = dynamic_cast<ELEMENT *>(finite_element_pt(0))->nnode_1d();

  // find the corresponding element/node/spine on the matching boundary
  unsigned e_periodic, n_periodic;
  Spine *spine_pt_periodic;
  if (mode == 'x') {
    e_periodic = e - (this->Nx - 1);
    n_periodic = n - (this->Np - 1);
  } else if (mode == 'y') {
    e_periodic = e % n_x;
    n_periodic = n % n_p;
  } else if (mode == 'c') {
    e_periodic = 0;
    n_periodic = 0;
  } else {
    throw OomphLibError("Invalid mode in match_spine_periodic",
                        OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }
  spine_pt_periodic = element_node_pt(e_periodic, n_periodic)->spine_pt();

  // Get pointer to node
  {
    // Get pointer to node and periodic counterpart
    SpineNode *nod_pt = element_node_pt(e, n);
    SpineNode *nod_pt_periodic = element_node_pt(e_periodic, n_periodic);
    // Set the pointer to the spine
    nod_pt->spine_pt() = spine_pt_periodic;
    // Set the fraction
    nod_pt->fraction() = nod_pt_periodic->fraction();
    // Pointer to the mesh that implements the update fct
    nod_pt->spine_mesh_pt() = nod_pt_periodic->spine_mesh_pt();
  }

  // Loop vertically along the spine
  // Loop over the elements
  for (unsigned long k = 0; k < n_z; k++) {
    // Loop over the vertical nodes, apart from the first
    for (unsigned l3 = 1; l3 < n_p; l3++) {
      // Get pointer to node and periodic counterpart
      SpineNode *nod_pt =
          element_node_pt(e + k * n_x * n_y, l3 * n_p * n_p + n);
      SpineNode *nod_pt_periodic = element_node_pt(e_periodic + k * n_x * n_y,
                                                   l3 * n_p * n_p + n_periodic);
      // Set the pointer to the spine
      nod_pt->spine_pt() = spine_pt_periodic;
      // Set the fraction
      nod_pt->fraction() = nod_pt_periodic->fraction();
      // Pointer to the mesh that implements the update fct
      nod_pt->spine_mesh_pt() = nod_pt_periodic->spine_mesh_pt();
    }
  }
}

//===========================================================================
/// Helper function that actually builds the single-layer spine mesh
/// based on the parameters set in the various constructors
//===========================================================================
template <class ELEMENT>
void SingleLayerCubicSpineMesh<ELEMENT>::build_single_layer_mesh(
    TimeStepper *time_stepper_pt) {
  // Read out the number of elements in the x-direction
  unsigned n_x = this->Nx;
  // Read out the number of elements in the y-direction
  unsigned n_y = this->Ny;
  // Read out number of linear points in the element
  unsigned n_p = dynamic_cast<ELEMENT *>(finite_element_pt(0))->nnode_1d();

  // Allocate store for the spines: (different in the case of periodic meshes
  // !!)
  if (this->Xperiodic && this->Yperiodic) {
    Spine_pt.reserve((n_p - 1) * n_x * (n_p - 1) * n_y);
  } else if (this->Xperiodic) {
    Spine_pt.reserve((n_p - 1) * n_x * ((n_p - 1) * n_y + 1));
  } else if (this->Yperiodic) {
    Spine_pt.reserve(((n_p - 1) * n_x + 1) * (n_p - 1) * n_y);
  } else {
    Spine_pt.reserve(((n_p - 1) * n_x + 1) * ((n_p - 1) * n_y + 1));
  }

  // Now we loop over all the elements and attach the spines

  // FIRST ELEMENT: Element 0
  // Loop over the nodes on the base of the element
  for (unsigned ly = 0; ly < n_p; ly++) {
    for (unsigned lx = 0; lx < n_p; lx++) {
      this->build_spine(0, lx + ly * n_p);
    }
  }

  // LOOP OVER OTHER ELEMENTS IN THE FIRST ROW
  //-----------------------------------------

  // The procedure is the same but we have to identify the
  // before defined spines for not defining them two times
  for (unsigned col = 1; col < n_x; col++) {
    for (unsigned ly = 0; ly < n_p; ly++) {
      // First we copy the last row of nodes into the
      // first row of the new element (and extend to the third dimension)
      for (unsigned lx = 1; lx < n_p; lx++) {
        this->build_spine(col, lx + ly * n_p);
      }
    }
  }

  // REST OF THE ELEMENTS
  // Now we loop over the rest of the elements.
  // We will separate the first of each row being al the rest equal
  for (unsigned long row = 1; row < n_y; row++) {
    // FIRST ELEMENT OF THE ROW

    // First line of nodes is copied from the element of the bottom
    for (unsigned ly = 1; ly < n_p; ly++) {
      for (unsigned lx = 0; lx < n_p; lx++) {
        this->build_spine(row * n_x, lx + ly * n_p);
      }
    }

    // REST OF THE ELEMENTS OF THE ROW
    for (unsigned col = 1; col < n_x; col++) {
      // First line of nodes is copied from the element of the bottom
      for (unsigned ly = 1; ly < n_p; ly++) {
        for (unsigned lx = 1; lx < n_p; lx++) {
          this->build_spine(col + row * n_x, lx + ly * n_p);
        }
      }
    }
  }
}

} // namespace oomph

#endif // SINGLE_LAYER_CUBIC_SPINE_MESH_H
