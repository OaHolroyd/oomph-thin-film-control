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
  /// Helper function to build a single spine for the nth node of element e
  virtual void build_spine(unsigned e, unsigned n);

  /// Helper function to build a single spine for the nth node of element e
  /// in the case of periodicity in x. (e, n) must be on the right boundary (2)
  virtual void build_spine_periodic_x(unsigned e, unsigned n);

  /// Helper function to build a single spine for the nth node of element e
  /// in the case of periodicity in y. (e, n) must be on the rear boundary (3)
  virtual void build_spine_periodic_y(unsigned e, unsigned n);

  /// Helper function to build a single spine for the nth node of element e
  /// in the case of periodicity in x and y. (e, n) must be at the rear right
  virtual void build_spine_periodic_xy(unsigned e, unsigned n);

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
    : CubicBrickMesh<ELEMENT>(nx, ny, nz, lx, ly, h, false, time_stepper_pt) {

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
                              periodic_in_y, false, time_stepper_pt) {

  // Mesh can only be built with 3D Qelements.
  MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

  // Mesh can only be built with spine elements
  MeshChecker::assert_geometric_element<SpineFiniteElement, ELEMENT>(3);

  // Build the single layer mesh
  build_single_layer_mesh(time_stepper_pt);
}

template <class ELEMENT>
void SingleLayerSpineMesh3D<ELEMENT>::build_spine(const unsigned e,
                                                  const unsigned n) {
  // Read out the number of elements in the x-direction
  unsigned n_x = this->Nx;
  // Read out the number of elements in the y-direction
  unsigned n_y = this->Ny;
  // Read out the number of elements in the z-direction
  unsigned n_z = this->Nz;
  // Read out number of linear points in the element
  unsigned n_p = dynamic_cast<ELEMENT *>(finite_element_pt(0))->nnode_1d();

  // create a new spine
  Spine *new_spine_pt = new Spine(1.0);
  Spine_pt.push_back(new_spine_pt);

  // Add the node at z = Zmin (usually this would be the top node of the
  // element down but for the bottom element there is no previous
  // element)
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

  // Loop vertically up the elements
  for (unsigned long k = 0; k < n_z; k++) {
    // Loop over the vertical nodes (apart from the first, because it is
    // the same as the last from the element down)
    for (unsigned l3 = 1; l3 < n_p; l3++) {
      // Get pointer to node
      SpineNode *nod_pt =
          element_node_pt(e + k * n_x * n_y, l3 * n_p * n_p + n);
      // Set the pointer to the spine
      nod_pt->spine_pt() = new_spine_pt;
      // Set the fraction
      nod_pt->fraction() = (double(k + l3) / double(n_p - 1)) / double(n_z);
      // Pointer to the mesh that implements the update fct
      nod_pt->spine_mesh_pt() = this;
    } // end l3-loop over vertical nodes
  }   // end k-loop over vertical elements
}

template <class ELEMENT>
void SingleLayerSpineMesh3D<ELEMENT>::build_spine_periodic_x(const unsigned e,
                                                             const unsigned n) {
  // Read out the number of elements in the x-direction
  unsigned n_x = this->Nx;
  // Read out the number of elements in the y-direction
  unsigned n_y = this->Ny;
  // Read out the number of elements in the z-direction
  unsigned n_z = this->Nz;
  // Read out number of linear points in the element
  unsigned n_p = dynamic_cast<ELEMENT *>(finite_element_pt(0))->nnode_1d();

  // find the corresponding element on the left boundary
  unsigned e_left = e - (this->Nx - 1);
  // find the corresponding node on the left boundary
  unsigned n_left = n - (this->Np - 1);
  // find the corresponding spine on the left boundary
  Spine *spine_pt_left = element_node_pt(e_left, n_left)->spine_pt();

  // check that this matches the one that we expect
  Spine *spine_pt_test = Spine_pt[e_left * (this->Np - 1) * (this->Np - 1) +
                                  (n_left / this->Np) * (this->Np - 1)];
  assert(spine_pt_left == spine_pt_test);

  // Add the node at z = Zmin (usually this would be the top node of the
  // element down but for the bottom element there is no previous
  // element)
  {
    // Get pointer to node
    SpineNode *nod_pt = element_node_pt(e, n);
    // Get pointer to corresponding node on the left boundary
    SpineNode *nod_pt_left = element_node_pt(e_left, n_left);
    // Set the pointer to the spine
    nod_pt->spine_pt() = spine_pt_left;
    // Set the fraction
    nod_pt->fraction() = nod_pt_left->fraction();
    // Pointer to the mesh that implements the update fct
    nod_pt->spine_mesh_pt() = nod_pt_left->spine_mesh_pt();

    // Check that the node is on a periodic boundary
    assert(nod_pt->x(0) == this->Xmax);
    assert(nod_pt_left->x(0) == this->Xmin);
    assert(nod_pt->x(1) == nod_pt_left->x(1));
    assert(nod_pt->x(2) == nod_pt_left->x(2));
  }

  // Loop vertically up the elements
  for (unsigned long k = 0; k < n_z; k++) {
    // Loop over the vertical nodes (apart from the first, because it is
    // the same as the last from the element down)
    for (unsigned l3 = 1; l3 < n_p; l3++) {
      // Get pointer to node
      SpineNode *nod_pt =
          element_node_pt(e + k * n_x * n_y, l3 * n_p * n_p + n);
      // Get pointer to corresponding node on the left boundary
      SpineNode *nod_pt_left =
          element_node_pt(e_left + k * n_x * n_y, l3 * n_p * n_p + n_left);
      // Set the pointer to the spine
      nod_pt->spine_pt() = spine_pt_left;
      // Set the fraction
      nod_pt->fraction() = nod_pt_left->fraction();
      // Pointer to the mesh that implements the update fct
      nod_pt->spine_mesh_pt() = nod_pt_left->spine_mesh_pt();

      // Check that the node is on a periodic boundary
      assert(nod_pt->x(0) == this->Xmax);
      assert(nod_pt_left->x(0) == this->Xmin);
      assert(nod_pt->x(1) == nod_pt_left->x(1));
      assert(nod_pt->x(2) == nod_pt_left->x(2));
    } // end l3-loop over vertical nodes
  }   // end k-loop over vertical elements
}

template <class ELEMENT>
void SingleLayerSpineMesh3D<ELEMENT>::build_spine_periodic_y(const unsigned e,
                                                             const unsigned n) {
  // Read out the number of elements in the x-direction
  unsigned n_x = this->Nx;
  // Read out the number of elements in the y-direction
  unsigned n_y = this->Ny;
  // Read out the number of elements in the z-direction
  unsigned n_z = this->Nz;
  // Read out number of linear points in the element
  unsigned n_p = dynamic_cast<ELEMENT *>(finite_element_pt(0))->nnode_1d();

  // find the corresponding element on the front boundary
  unsigned e_front = e % n_x;
  // find the corresponding node on the front boundary
  unsigned n_front = n % n_p;
  // find the corresponding spine on the front boundary
  Spine *spine_pt_front = element_node_pt(e_front, n_front)->spine_pt();

  // check that this matches the one that we expect
  Spine *spine_pt_test = Spine_pt[e_front * (this->Np - 1) * (this->Np - 1) + n_front];
  assert(spine_pt_front == spine_pt_test);

  // Add the node at z = Zmin (usually this would be the top node of the
  // element down but for the bottom element there is no previous
  // element)
  {
    // Get pointer to node
    SpineNode *nod_pt = element_node_pt(e, n);
    // Get pointer to corresponding node on the front boundary
    SpineNode *nod_pt_front = element_node_pt(e_front, n_front);
    // Set the pointer to the spine
    nod_pt->spine_pt() = spine_pt_front;
    // Set the fraction
    nod_pt->fraction() = nod_pt_front->fraction();
    // Pointer to the mesh that implements the update fct
    nod_pt->spine_mesh_pt() = nod_pt_front->spine_mesh_pt();

    // Check that the node is on a periodic boundary
    assert(nod_pt->x(1) == this->Ymax);
    assert(nod_pt_front->x(1) == this->Ymin);
    assert(nod_pt->x(0) == nod_pt_front->x(0));
    assert(nod_pt->x(2) == nod_pt_front->x(2));
  }

  // Loop vertically up the elements
  for (unsigned long k = 0; k < n_z; k++) {
    // Loop over the vertical nodes (apart from the first, because it is
    // the same as the last from the element down)
    for (unsigned l3 = 1; l3 < n_p; l3++) {
      // Get pointer to node
      SpineNode *nod_pt =
          element_node_pt(e + k * n_x * n_y, l3 * n_p * n_p + n);
      // Get pointer to corresponding node on the front boundary
      SpineNode *nod_pt_front =
          element_node_pt(e_front + k * n_x * n_y, l3 * n_p * n_p + n_front);
      // Set the pointer to the spine
      nod_pt->spine_pt() = spine_pt_front;
      // Set the fraction
      nod_pt->fraction() = nod_pt_front->fraction();
      // Pointer to the mesh that implements the update fct
      nod_pt->spine_mesh_pt() = nod_pt_front->spine_mesh_pt();

      // Check that the node is on a periodic boundary
      assert(nod_pt->x(1) == this->Ymax);
      assert(nod_pt_front->x(1) == this->Ymin);
      assert(nod_pt->x(0) == nod_pt_front->x(0));
      assert(nod_pt->x(2) == nod_pt_front->x(2));
    } // end l3-loop over vertical nodes
  }   // end k-loop over vertical elements
}

template <class ELEMENT>
void SingleLayerSpineMesh3D<ELEMENT>::build_spine_periodic_xy(
    const unsigned e, const unsigned n) {
  // Read out the number of elements in the x-direction
  unsigned n_x = this->Nx;
  // Read out the number of elements in the y-direction
  unsigned n_y = this->Ny;
  // Read out the number of elements in the z-direction
  unsigned n_z = this->Nz;
  // Read out number of linear points in the element
  unsigned n_p = dynamic_cast<ELEMENT *>(finite_element_pt(0))->nnode_1d();

  // find the corresponding element at the front left corner
  unsigned e_fl = 0;
  // find the corresponding node at the front left corner
  unsigned n_fl = 0;
  // find the corresponding spine at the front left corner
  Spine *spine_pt_fl = element_node_pt(e_fl, n_fl)->spine_pt();

  // check that this matches the one that we expect
  Spine *spine_pt_test = Spine_pt[0];
  assert(spine_pt_fl == spine_pt_test);

  // Add the node at z = Zmin (usually this would be the top node of the
  // element down but for the bottom element there is no previous
  // element)
  {
    // Get pointer to node
    SpineNode *nod_pt = element_node_pt(e, n);
    // Get pointer to corresponding node at the front left corner
    SpineNode *nod_pt_fl = element_node_pt(e_fl, n_fl);
    // Set the pointer to the spine
    nod_pt->spine_pt() = spine_pt_fl;
    // Set the fraction
    nod_pt->fraction() = nod_pt_fl->fraction();
    // Pointer to the mesh that implements the update fct
    nod_pt->spine_mesh_pt() = nod_pt_fl->spine_mesh_pt();

    // Check that the node is on a doubley periodic boundary
    assert(nod_pt->x(0) == this->Xmax);
    assert(nod_pt->x(1) == this->Ymax);
    assert(nod_pt_fl->x(0) == this->Xmin);
    assert(nod_pt_fl->x(1) == this->Ymin);
    assert(nod_pt->x(2) == nod_pt_fl->x(2));
  }

  // Loop vertically up the elements
  for (unsigned long k = 0; k < n_z; k++) {
    // Loop over the vertical nodes (apart from the first, because it is
    // the same as the last from the element down)
    for (unsigned l3 = 1; l3 < n_p; l3++) {
      // Get pointer to node
      SpineNode *nod_pt =
          element_node_pt(e + k * n_x * n_y, l3 * n_p * n_p + n);
      // Get pointer to corresponding node at the front left corner
      SpineNode *nod_pt_fl =
          element_node_pt(e_fl + k * n_x * n_y, l3 * n_p * n_p + n_fl);
      // Set the pointer to the spine
      nod_pt->spine_pt() = spine_pt_fl;
      // Set the fraction
      nod_pt->fraction() = nod_pt_fl->fraction();
      // Pointer to the mesh that implements the update fct
      nod_pt->spine_mesh_pt() = nod_pt_fl->spine_mesh_pt();

      // Check that the node is on a periodic boundary
      assert(nod_pt->x(0) == this->Xmax);
      assert(nod_pt->x(1) == this->Ymax);
      assert(nod_pt_fl->x(0) == this->Xmin);
      assert(nod_pt_fl->x(1) == this->Ymin);
      assert(nod_pt->x(2) == nod_pt_fl->x(2));
    } // end l3-loop over vertical nodes
  }   // end k-loop over vertical elements
}

template <class ELEMENT>
void SingleLayerSpineMesh3D<ELEMENT>::build_single_layer_mesh(
    TimeStepper *time_stepper_pt) {
  // mesh can only be built with 3D Qelements
  MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

  // build the underlying brick mesh
  fprintf(stderr, "build_single_layer_mesh calling build mesh\n");
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
    Spine_pt.reserve((n_p - 1) * n_x * (n_p - 1) * n_y);
  } else if (this->Xperiodic) {
    Spine_pt.reserve((n_p - 1) * n_x * ((n_p - 1) * n_y + 1));
  } else if (this->Yperiodic) {
    Spine_pt.reserve(((n_p - 1) * n_x + 1) * (n_p - 1) * n_y);
  } else {
    Spine_pt.reserve(((n_p - 1) * n_x + 1) * ((n_p - 1) * n_y + 1));
  }

  // loop over all of the elements
  for (unsigned i = 0; i < n_y; i++) {
    for (unsigned j = 0; j < n_x; j++) {
      unsigned e = j + i * n_x; // linear element number

      // for each element we add the front left np-1 x np-1 spines (so that they
      // don't overlap with those of a previous element)
      for (unsigned l1 = 0; l1 < n_p - 1; l1++) {
        for (unsigned l2 = 0; l2 < n_p - 1; l2++) {
          unsigned n = l2 + l1 * n_p; // linear node number

          this->build_spine(e, n);
        }
      } // end l1-loop over y nodes
    }
  } // end i-loop over y nodes

  // the main loop over the elements has missed out the right/rear boundary
  // nodes of the last row/column of elements

  // right boundary (2)
  {
    unsigned j = n_x - 1;  // last column of elements
    unsigned l2 = n_p - 1; // last column of nodes inside each element

    // loop over the last column of elements
    for (unsigned i = 0; i < n_y; i++) {
      unsigned e = j + i * n_x; // linear element number

      // loop over the last column of nodes inside each element
      for (unsigned l1 = 0; l1 < n_p - 1; l1++) {
        unsigned n = l2 + l1 * n_p; // linear node number

        if (this->Xperiodic) {
          this->build_spine_periodic_x(e, n);
        } else {
          this->build_spine(e, n);
        }
      }
    }
  }

  // rear boundary (3)
  {
    unsigned i = n_y - 1;  // last row of elements
    unsigned l1 = n_p - 1; // last row of nodes inside each element

    // loop over the last row of elements
    for (unsigned j = 0; j < n_x; j++) {
      unsigned e = j + i * n_x; // linear element number

      // loop over the last column of nodes inside each element
      for (unsigned l2 = 0; l2 < n_p - 1; l2++) {
        unsigned n = l2 + l1 * n_p; // linear node number

        if (this->Yperiodic) {
          this->build_spine_periodic_y(e, n);
        } else {
          this->build_spine(e, n);
        }
      }
    }
  }

  // right/rear node
  {
    unsigned i = n_y - 1;  // last row of elements
    unsigned j = n_x - 1;  // last column of elements
    unsigned l1 = n_p - 1; // last row of nodes inside each element
    unsigned l2 = n_p - 1; // last column of nodes inside each element

    unsigned e = j + i * n_x;   // linear element number
    unsigned n = l2 + l1 * n_p; // linear node number

    if (this->Xperiodic && this->Yperiodic) {
      this->build_spine_periodic_xy(e, n);
    } else if (this->Xperiodic) {
      this->build_spine_periodic_x(e, n);
    } else if (this->Yperiodic) {
      this->build_spine_periodic_y(e, n);
    } else {
      this->build_spine(e, n);
    }
  }
}

#endif // SINGLELAYERSPINEMESH3D_H
