//
// Created by Oscar Holroyd on 07/02/2025.
//

#ifndef SINGLELAYERSPINEMESH3D_H
#define SINGLELAYERSPINEMESH3D_H

// oomph-lib includes
#include "generic/spines.h"

// Project-specific includes
#include "cubic_brick_mesh.h"

template <class ELEMENT>
class SingleLayerSpineMesh3D : public CubicBrickMesh<ELEMENT>,
                               public SpineMesh {
public:
  // basic constructor
  SingleLayerSpineMesh3D(
      const unsigned &nx, const unsigned &ny, const unsigned &nz,
      const double &lx, const double &ly, const double &h,
      TimeStepper *time_stepper_pt = &Mesh::Default_TimeStepper);

  // constructor with periodicity
  SingleLayerSpineMesh3D(
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
    // Set the value of y
    spine_node_pt->x(1) = this->Ymin + W * H;
  }

protected:
  /// Helper function to build a single spine for the (lx, ly)-node of element
  /// (ex, ey)
  virtual void build_spine(unsigned ex, unsigned ey, unsigned ly, unsigned lx);

  /// Helper function to actually build the single-layer spine mesh
  /// (called from various constructors)
  virtual void build_single_layer_mesh(TimeStepper *time_stepper_pt);
};

// IMPLEMENTATION

template <class ELEMENT>
SingleLayerSpineMesh3D<ELEMENT>::SingleLayerSpineMesh3D(
    const unsigned &nx, const unsigned &ny, const unsigned &nz,
    const double &lx, const double &ly, const double &h,
    TimeStepper *time_stepper_pt)
    : CubicBrickMesh<ELEMENT>(nx, ny, nz, lx, ly, h, time_stepper_pt) {

  // Mesh can only be built with 3D Qelements.
  MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

  // Mesh can only be built with spine elements
  MeshChecker::assert_geometric_element<SpineFiniteElement, ELEMENT>(3);

  // Build the single layer mesh
  build_single_layer_mesh(time_stepper_pt);
}

template <class ELEMENT>
SingleLayerSpineMesh3D<ELEMENT>::SingleLayerSpineMesh3D(
    const unsigned &nx, const unsigned &ny, const unsigned &nz,
    const double &lx, const double &ly, const double &h,
    const bool &periodic_in_x, const bool &periodic_in_y,
    TimeStepper *time_stepper_pt)
    : CubicBrickMesh<ELEMENT>(nx, ny, nz, lx, ly, h, periodic_in_x,
                              periodic_in_y, time_stepper_pt) {

  // Mesh can only be built with 3D Qelements.
  MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

  // Mesh can only be built with spine elements
  MeshChecker::assert_geometric_element<SpineFiniteElement, ELEMENT>(3);

  // Build the single layer mesh
  build_single_layer_mesh(time_stepper_pt);
}

template <class ELEMENT>
void SingleLayerSpineMesh3D<ELEMENT>::build_spine(const unsigned ex,
                                                  const unsigned ey,
                                                  const unsigned ly,
                                                  const unsigned lx) {

  // Read out the number of elements in the x-direction
  const unsigned nx = this->Nx;
  // Read out the number of elements in the y-direction
  const unsigned ny = this->Ny;
  // Read out the number of elements in the z-direction
  const unsigned nz = this->Nz;
  // Read out number of linear points in the element
  const unsigned np = dynamic_cast<ELEMENT *>(finite_element_pt(0))->nnode_1d();

  // Find the element number
  const unsigned e = ex + ey * nx;
  // Find the node number
  const unsigned n = lx + ly * np;

  // Check if this is a truely new spine or the image of a periodic one that
  // has already been created
  // NOTE: finding the periodic counterpart relies on creating the spines in
  //       a particular order (see build_single_layer_mesh).
  bool is_periodic = true;
  Spine *periodic_spine_pt = nullptr;
  if (this->Xperiodic && this->Yperiodic && (e == nx * ny - 1) &&
      (lx == np - 1) && (ly == np - 1)) {
    // corner case where the mesh is periodic in both x and y
    periodic_spine_pt = Spine_pt[0];
  } else if (this->Xperiodic && (ex == nx - 1) && (lx == np - 1)) {
    // periodic in x
    // have to treat the first row differently
    if (ey == 0) {
      periodic_spine_pt = Spine_pt[ly * np];
    } else {
      // number of spines in the first row of elements
      unsigned row0 = np * (np - 1) * nx;
      // number of spines in all subsequent rows of elements
      unsigned row = (np - 1) * (np - 1) * nx;

      periodic_spine_pt = Spine_pt[(ly - 1) * np + row0 + row * (ey - 1)];
    }
  } else if (this->Yperiodic && (ey == ny - 1) && (ly == np - 1)) {
    // periodic in y
    // have to treat the first column differently
    if (ex == 0) {
      periodic_spine_pt = Spine_pt[lx];
    } else {
      // number of spines in the first element
      unsigned el0 = np * np;
      // number of spines in all subsequent elements
      unsigned el = (np - 1) * np;

      periodic_spine_pt = Spine_pt[(lx - 1) + el0 + el * (ex - 1)];
    }
  } else {
    is_periodic = false;
  }

  // catch the corner case where the mesh is periodic in both x and y

  if (!is_periodic) {
    // NON-PERIODIC SPINE

    // Assign the new spine with length h
    // TODO: should the height be h or 1.0?
    Spine *new_spine_pt = new Spine(1.0);
    Spine_pt.push_back(new_spine_pt);

    // Add the node at z = Zmin (usually this would be the top node of the
    // element down but for the bottom element there is no previous element)
    {
      // Get pointer to node
      SpineNode *nod_pt = element_node_pt(e, n);
      // Set the pointer to the spine
      nod_pt->spine_pt() = new_spine_pt;
      // Set the fraction
      nod_pt->fraction() = 0.0;
      // Pointer to the mesh that implements the update fct
      nod_pt->spine_mesh_pt() = this;

      // Check that the node is not on a periodic boundary
      assert(!this->Xperiodic || (nod_pt->x(0) != this->Xmax));
      assert(!this->Yperiodic || (nod_pt->x(1) != this->Ymax));
    }

    // Loop vertically up the spine
    // Loop over the elements
    for (unsigned long k = 0; k < nz; k++) {
      // Loop over the vertical nodes (apart from the first, because it is the
      // same as the last from the element down)
      for (unsigned l3 = 1; l3 < np; l3++) {
        // Get pointer to node
        SpineNode *nod_pt = element_node_pt(e + k * nx * ny, l3 * np * np + n);
        // Set the pointer to the spine
        nod_pt->spine_pt() = new_spine_pt;
        // Set the fraction
        nod_pt->fraction() = (double(k + l3) / double(np - 1)) / double(nz);
        // Pointer to the mesh that implements the update fct
        nod_pt->spine_mesh_pt() = this;
      } // end l3-loop over vertical nodes
    }   // end k-loop over vertical elements
  } else {
    // PERIODIC SPINE
    // should not create a new spine since it is periodic

    // find the periodic element/node indices
    unsigned periodic_e = -1;
    unsigned periodic_n = -1;
    this->get_periodic_element_node(e, n, &periodic_e, &periodic_n);

    // Add the node at z = Zmin (usually this would be the top node of the
    // element down but for the bottom element there is no previous element)
    {
      // Get pointer to node
      SpineNode *nod_pt = element_node_pt(e, n);
      // Get the pointer to its periodic counterpart
      SpineNode *periodic_nod_pt = element_node_pt(periodic_e, periodic_n);
      // Set the pointer to the spine
      nod_pt->spine_pt() = periodic_spine_pt;
      // Set the fraction (fixed to the periodic counterpart)
      nod_pt->fraction() = periodic_nod_pt->fraction();
      // Pointer to the mesh that implements the update fct (fixed to the
      // periodic counterpart)
      nod_pt->spine_mesh_pt() = periodic_nod_pt->spine_mesh_pt();
    }

    // Loop vertically up the spine
    // Loop over the elements
    for (unsigned long k = 0; k < nz; k++) {
      // Loop over the vertical nodes (apart from the first, because it is the
      // same as the last from the element down)
      for (unsigned l3 = 1; l3 < np; l3++) {
        // Get pointer to node
        SpineNode *nod_pt = element_node_pt(e + k * nx * ny, l3 * np * np + n);
        // Get the pointer to its periodic counterpart
        SpineNode *periodic_nod_pt = element_node_pt(periodic_e + k * nx * ny,
                                                     l3 * np * np + periodic_n);
        // Set the pointer to the spine
        nod_pt->spine_pt() = periodic_spine_pt;
        // Set the fraction (fixed to the periodic counterpart)
        nod_pt->fraction() = periodic_nod_pt->fraction();
        // Pointer to the mesh that implements the update fct (fixed to the
        // periodic counterpart)
        nod_pt->spine_mesh_pt() = periodic_nod_pt->spine_mesh_pt();
      } // end l3-loop over vertical nodes
    }   // end k-loop over vertical elements
  }
}

template <class ELEMENT>
void SingleLayerSpineMesh3D<ELEMENT>::build_single_layer_mesh(
    TimeStepper *time_stepper_pt) {
  // mesh can only be built with 3D Qelements
  MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

  // build the underlying brick mesh
  CubicBrickMesh<ELEMENT>::build_mesh(time_stepper_pt);

  // Read out the number of elements in the x-direction
  unsigned n_x = this->Nx;
  // Read out the number of elements in the y-direction
  unsigned n_y = this->Ny;
  // Read out number of linear points in the element
  unsigned n_p = dynamic_cast<ELEMENT *>(finite_element_pt(0))->nnode_1d();

  // Allocate memory for the spines and fractions along spines
  //---------------------------------------------------------
  // Allocate store for the spines:
  if (this->Xperiodic && this->Yperiodic) {
    Spine_pt.reserve((n_p - 1) * n_x * n_y);
  } else if (this->Xperiodic) {
    Spine_pt.reserve((n_p - 1) * n_x * (n_y + 1));
  } else if (this->Yperiodic) {
    Spine_pt.reserve((n_p - 1) * (n_x + 1) * n_y);
  } else {
    Spine_pt.reserve((n_p - 1) * (n_x + 1) * (n_y + 1));
  }

  // Now we loop over all the elements and attach the spines
  // TODO: the indexing would be much simpler if we just did the lower-left
  //       (np-1) x (np-1) spines for each element and treated the right/upper
  //       boundary nodes as special cases rather than the other way around

  // FIRST ELEMENT: Element 0
  // Loop over the nodes on the base of the element
  for (unsigned l1 = 0; l1 < n_p; l1++) {   // y loop over the nodes
    for (unsigned l2 = 0; l2 < n_p; l2++) { // x loop over the nodes
      build_spine(0, 0, l1, l2);
    } // end l2-loop over x nodes
  }   // end l1-loop over y nodes
  // END OF FIRST ELEMENT

  // LOOP OVER OTHER ELEMENTS IN THE FIRST ROW
  //-----------------------------------------
  // The procedure is the same but we have to identify the
  // before defined spines for not defining them two times
  for (unsigned j = 1; j < n_x; j++) { // loop over elements in the first row
    for (unsigned l1 = 0; l1 < n_p; l1++) { // y loop over the nodes
      // First we copy the last row of nodes into the
      // first row of the new element (and extend to the third dimension)
      for (unsigned l2 = 1; l2 < n_p; l2++) { // x loop over the nodes
        build_spine(j, 0, l1, l2);
      } // end l2-loop over x nodes
    }   // end l1-loop over y nodes
  }     // end j-loop over elements in the first row
  // END OF FIRST ROW

  // REST OF THE ELEMENTS
  // Now we loop over the rest of the elements.
  // We will separate the first of each row being al the rest equal
  for (unsigned long i = 1; i < n_y; i++) {
    // FIRST ELEMENT OF THE ROW

    // First line of nodes is copied from the element of the bottom
    for (unsigned l1 = 1; l1 < n_p; l1++) {   // y loop over the nodes
      for (unsigned l2 = 0; l2 < n_p; l2++) { // x loop over the nodes
        build_spine(0, i, l1, l2);
      } // end l2-loop over x nodes
    }   // end l1-loop over y nodes
    // END OF FIRST ELEMENT OF THE ROW

    // REST OF THE ELEMENTS OF THE ROW
    for (unsigned j = 1; j < n_x; j++) {
      // First line of nodes is copied from the element of the bottom
      for (unsigned l1 = 1; l1 < n_p; l1++) {   // y loop over the nodes
        for (unsigned l2 = 1; l2 < n_p; l2++) { // x loop over the nodes
          build_spine(j, i, l1, l2);
        } // end l2-loop over x nodes
      }   // end l1-loop over y nodes
    }     // end j-loop over elements in the row
    // END OF REST OF THE ELEMENTS OF THE ROW
  } // end i-loop over rows
}

#endif // SINGLELAYERSPINEMESH3D_H
